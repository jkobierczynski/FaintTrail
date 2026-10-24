// SPDX-License-Identifier: GPL-3.0-or-later
// Copyright 2026 Jurgen Kobierczynski
//
// Tiled search (see tiled.hpp). The per-pixel arithmetic is the untiled code
// itself (build_stage, the backend, refine); this file only decides which part
// of the image each step sees, and gathers the whole-frame quantities first.
//
// Per time order (the real one, then each null permutation):
//   A  residual noise per frame: residuals at the sample positions only, plus
//      (real order) the static sky of the whole field for the PSF width;
//   B  matched-filter normalisation per frame: tiles with the filter halo,
//      filtered but not normalised, samples collected over the core;
//   C  tiles with filter and stacking halos, fully prepared and stacked; the
//      cores of the detection maps are assembled into whole-field maps.
// Peaks are found on the whole-field maps, the threshold is fitted, and the
// candidates are refined on tiles with a refinement halo.
#include "fainttrail/tiled.hpp"

#include <algorithm>
#include <atomic>
#include <chrono>
#include <cmath>
#include <stdexcept>

#include "common.hpp"
#include "fainttrail/backend.hpp"
#include "fainttrail/fits.hpp"
#include "internal.hpp"
#include "psf.hpp"

namespace fainttrail {

// ---------------------------------------------------------------- sources
namespace {

class MemorySource final : public FrameSource {
public:
    explicit MemorySource(const std::vector<Frame>& f) : frames_(f) {
        if (f.empty()) throw std::runtime_error("no frames");
        for (const auto& fr : f) {
            if (fr.image.width != f[0].image.width || fr.image.height != f[0].image.height)
                throw std::runtime_error("all frames must have the same size (are they registered?)");
            FrameInfo in;
            in.mjd_mid = fr.mjd_mid;
            in.exptime_s = fr.exptime_s;
            in.wcs = fr.wcs;
            in.has_reg_offset = fr.has_reg_offset;
            in.reg_dx = fr.reg_dx;
            in.reg_dy = fr.reg_dy;
            in.source = fr.source;
            info_.push_back(in);
        }
    }
    int count() const override { return static_cast<int>(frames_.size()); }
    int width() const override { return frames_[0].image.width; }
    int height() const override { return frames_[0].image.height; }
    const FrameInfo& info(int i) const override { return info_[static_cast<std::size_t>(i)]; }
    void read(int i, int x0, int y0, int w, int h, Image& out) const override {
        const Image& src = frames_[static_cast<std::size_t>(i)].image;
        if (x0 < 0 || y0 < 0 || x0 + w > src.width || y0 + h > src.height) throw std::runtime_error("region outside the frame");
        out = Image(w, h);
        for (int y = 0; y < h; ++y)
            std::copy_n(&src.data[static_cast<std::size_t>(y0 + y) * src.width + x0], w,
                        &out.data[static_cast<std::size_t>(y) * w]);
    }
    void set_wcs(const Wcs& w) override {
        for (auto& in : info_) in.wcs = w;
    }

private:
    const std::vector<Frame>& frames_;
    std::vector<FrameInfo> info_;
};

class FitsSource final : public FrameSource {
public:
    explicit FitsSource(const std::vector<std::string>& paths) {
        if (paths.empty()) throw std::runtime_error("no input FITS files given");
        for (const auto& p : paths) {
            FitsLayout L = fits_layout(p);
            FitsImage header_only;
            header_only.header = L.header;
            const Frame f = frame_from_fits(std::move(header_only), p);
            FrameInfo in;
            in.mjd_mid = f.mjd_mid;
            in.exptime_s = f.exptime_s;
            in.wcs = f.wcs;
            in.has_reg_offset = f.has_reg_offset;
            in.reg_dx = f.reg_dx;
            in.reg_dy = f.reg_dy;
            in.source = p;
            if (!layouts_.empty() && (L.width != layouts_[0].width || L.height != layouts_[0].height))
                throw std::runtime_error("all frames must have the same size (are they registered?)");
            layouts_.push_back(std::move(L));
            info_.push_back(in);
        }
    }
    int count() const override { return static_cast<int>(layouts_.size()); }
    int width() const override { return layouts_[0].width; }
    int height() const override { return layouts_[0].height; }
    const FrameInfo& info(int i) const override { return info_[static_cast<std::size_t>(i)]; }
    void read(int i, int x0, int y0, int w, int h, Image& out) const override {
        read_fits_region(layouts_[static_cast<std::size_t>(i)], x0, y0, w, h, out);
    }
    void set_wcs(const Wcs& w) override {
        for (auto& in : info_) in.wcs = w;
    }

private:
    std::vector<FitsLayout> layouts_;
    std::vector<FrameInfo> info_;
};

struct Rect {
    int x0 = 0, y0 = 0, x1 = 0, y1 = 0;  // [x0,x1) x [y0,y1)
    int w() const { return x1 - x0; }
    int h() const { return y1 - y0; }
};

Rect grow_clip(const Rect& r, int by, int W, int H) {
    return {std::max(0, r.x0 - by), std::max(0, r.y0 - by), std::min(W, r.x1 + by), std::min(H, r.y1 + by)};
}

}  // namespace

std::unique_ptr<FrameSource> memory_frame_source(const std::vector<Frame>& frames) {
    return std::make_unique<MemorySource>(frames);
}

std::unique_ptr<FrameSource> fits_frame_source(const std::vector<std::string>& paths) {
    return std::make_unique<FitsSource>(paths);
}

TransferReport apply_wcs_from_file(const std::string& solved_fits, FrameSource& frames, double scale_hint) {
    FitsImage solved = read_fits(solved_fits);
    const bool distortion = solved.header.has("A_ORDER") || solved.header.has("PV1_1") || solved.header.has("PV2_1");
    const Wcs sw = wcs_from_header(solved.header);
    if (!sw.valid) throw std::runtime_error(solved_fits + " has no TAN plate solution (CTYPE RA---TAN with CD or CDELT keywords)");
    Image ref;
    frames.read(frames.count() / 2, 0, 0, frames.width(), frames.height(), ref);
    TransferReport rep;
    const Wcs w = transfer_wcs(sw, solved.image, ref, &rep, scale_hint);
    rep.source_has_distortion = distortion;
    frames.set_wcs(w);
    return rep;
}

// ---------------------------------------------------------------- planning
namespace {

constexpr double kMB = 1024.0 * 1024.0;

struct Geometry {
    int stack_halo = 0, filter_halo = 0, refine_halo = 0;
};

Geometry geometry(const Prepared& meta, const std::vector<std::pair<double, double>>& grid, double dv, double fwhm) {
    Geometry g;
    double tmax = 0.0;
    for (double t : meta.t_h) tmax = std::max(tmax, std::fabs(t));
    int s = 0;
    for (const auto& v : grid)
        for (double t : meta.t_h)
            s = std::max({s, std::abs(static_cast<int>(std::lround(v.first * t))),
                          std::abs(static_cast<int>(std::lround(v.second * t)))});
    g.stack_halo = s;
    const double psf_sigma = fwhm * detail::kFwhmToSigma;
    const int R = std::max(1, static_cast<int>(std::ceil(3.0 * psf_sigma)));
    const int dil = std::max(1, static_cast<int>(std::lround(0.5 * fwhm)));
    g.filter_halo = R + 2 + dil + 1;  // matched filter, cosmic-ray test (+-2), mask dilation
    double reg_range = 0.0;
    if (!meta.reg_dx.empty()) {
        const auto [ax, bx] = std::minmax_element(meta.reg_dx.begin(), meta.reg_dx.end());
        const auto [ay, by] = std::minmax_element(meta.reg_dy.begin(), meta.reg_dy.end());
        reg_range = std::max(*bx - *ax, *by - *ay);
    }
    // refinement moves position by up to 1.5 px and motion by up to 1.25 dv, then
    // samples the track (bilinear), photometry and edge boxes around it, and the
    // detector hypothesis, which follows the registration offsets
    const int track = s + static_cast<int>(std::ceil(1.25 * dv * tmax)) + 2;
    g.refine_halo = 2 * track + static_cast<int>(std::ceil(reg_range)) + R + 2 +
                    static_cast<int>(std::ceil(2.0 * fwhm)) + 1 + static_cast<int>(std::ceil(std::max(2.0, 2.0 * fwhm))) +
                    4 + g.filter_halo;
    return g;
}

double tiled_memory_mb(int W, int H, int n, int nt, int tile, const Geometry& g) {
    const double full = 16.0 * W * H;  // two detection maps, motion indices, static sky
    const double side_c = std::min<double>(W, tile + 2.0 * (g.stack_halo + g.filter_halo));
    const double side_c2 = std::min<double>(H, tile + 2.0 * (g.stack_halo + g.filter_halo));
    const double region = side_c * side_c2;
    // background-subtracted + residual + filtered frames of the region (the
    // residuals are released while the filtered frames are built), per-thread
    // working images, the region's detection maps
    const double work = (2.3 * n + 3.0 * nt + 4.0) * 4.0 * region;
    return (full + work) / kMB + 30.0;
}

}  // namespace

double estimate_untiled_memory_mb(int width, int height, int nframes) {
    // loaded frames are handed over; peak during a null run: 3 copies + working images
    const double set = 4.0 * width * height * nframes;
    return (3.3 * set + 24.0 * width * height) / kMB + 30.0;
}

namespace {

Prepared make_meta(const FrameSource& src) {
    const int n = src.count();
    if (n < 5) throw std::runtime_error("synthetic tracking needs at least 5 frames");
    const int w = src.width(), h = src.height();
    if (w < 16 || h < 16) throw std::runtime_error("frames are too small");
    Prepared p;
    p.width = w;
    p.height = h;
    p.nframes = n;
    double tmin = src.info(0).mjd_mid, tmax = tmin;
    for (int i = 0; i < n; ++i) {
        tmin = std::min(tmin, src.info(i).mjd_mid);
        tmax = std::max(tmax, src.info(i).mjd_mid);
    }
    p.t_ref_mjd = 0.5 * (tmin + tmax);
    p.span_h = (tmax - tmin) * 24.0;
    if (p.span_h <= 0.0) throw std::runtime_error("frames must span a non-zero time interval");
    for (int i = 0; i < n; ++i) {
        p.t_h.push_back((src.info(i).mjd_mid - p.t_ref_mjd) * 24.0);
        p.exptime_s.push_back(src.info(i).exptime_s);
    }
    for (int i = 0; i < n; ++i)
        if (src.info(i).wcs.valid) {
            p.wcs = src.info(i).wcs;
            break;
        }
    bool all = true;
    for (int i = 0; i < n; ++i) all = all && src.info(i).has_reg_offset;
    if (all)
        for (int i = 0; i < n; ++i) {
            p.reg_dx.push_back(src.info(i).reg_dx);
            p.reg_dy.push_back(src.info(i).reg_dy);
        }
    return p;
}

}  // namespace

TilePlan plan_tiles(const FrameSource& src, const SearchParams& params, const TileOptions& opt, double fwhm_px) {
    const Prepared meta = make_meta(src);
    double dv = 0.0;
    const auto grid = velocity_grid(meta, params, &dv);
    const Geometry g = geometry(meta, grid, dv, fwhm_px > 0.0 ? fwhm_px : 3.0);
    const int W = meta.width, H = meta.height, n = meta.nframes, nt = detail::thread_count(params.threads);
    TilePlan plan;
    plan.stack_halo = g.stack_halo;
    plan.filter_halo = g.filter_halo;
    if (opt.tile > 0) {
        plan.tile = opt.tile;
    } else {
        const double budget = opt.max_memory_mb > 0.0 ? opt.max_memory_mb : 2048.0;
        int t = (std::max(W, H) + 63) / 64 * 64;
        while (t > 64 && tiled_memory_mb(W, H, n, nt, t, g) > budget) t -= 64;
        plan.tile = std::max(64, t);
    }
    plan.tile = std::max(16, plan.tile);
    plan.tiles_x = (W + plan.tile - 1) / plan.tile;
    plan.tiles_y = (H + plan.tile - 1) / plan.tile;
    plan.memory_mb = tiled_memory_mb(W, H, n, nt, plan.tile, g);
    return plan;
}

// ---------------------------------------------------------------- search
SearchResult search_tiled(const FrameSource& src, const SearchParams& params_in, const TileOptions& opt,
                          Prepared* meta_out, const ProgressFn& progress) {
    using namespace detail;
    SearchParams params = params_in;
    params.keep_residuals = false;  // tiles never need them
    const auto t_start = std::chrono::steady_clock::now();
    Prepared meta = make_meta(src);
    const int n = meta.nframes, W = meta.width, H = meta.height;
    const int nt = thread_count(params.threads);
    Stopwatch psw(&meta.timings);
    auto report = [&](double f, const std::string& stage) {
        if (progress && !progress(std::clamp(f, 0.0, 1.0), stage)) throw Cancelled();
    };

    // ---- background level per frame: median of the same samples as prepare()
    report(0.0, "background");
    meta.background.assign(static_cast<std::size_t>(n), 0.0);
    {
        const SampleGrid sg = sample_grid(W, H, 200000, 0);
        std::atomic<bool> failed{false};
        std::string error;
#pragma omp parallel for num_threads(nt) schedule(dynamic)
        for (int i = 0; i < n; ++i) {
            try {
                std::vector<float> s;
                Image band;
                const int bh = 64;
                for (int y0 = 0; y0 < H; y0 += bh) {
                    const int h = std::min(bh, H - y0);
                    src.read(i, 0, y0, W, h, band);
                    for (int y = y0; y < y0 + h; ++y) {
                        if (y % sg.step) continue;
                        for (int x = 0; x < W; x += sg.step) {
                            const float v = band.at(x, y - y0);
                            if (std::isfinite(v)) s.push_back(v);
                        }
                    }
                }
                meta.background[static_cast<std::size_t>(i)] = median_inplace(s);
            } catch (const std::exception& e) {
#pragma omp critical
                error = e.what();
                failed = true;
            }
        }
        if (failed) throw std::runtime_error(error);
    }
    psw.lap("background");

    // background-subtracted frames of a rectangle, blank areas grown as in prepare()
    auto load_sub = [&](const Rect& r) {
        std::vector<Image> sub(static_cast<std::size_t>(n));
        const Rect m = grow_clip(r, kBlankGrow, W, H);
        std::atomic<bool> failed{false};
        std::string error;
#pragma omp parallel for num_threads(nt) schedule(dynamic)
        for (int i = 0; i < n; ++i) {
            try {
                Image raw;
                src.read(i, m.x0, m.y0, m.w(), m.h(), raw);
                subtract_and_grow(raw, meta.background[static_cast<std::size_t>(i)]);
                Image out(r.w(), r.h());
                for (int y = 0; y < r.h(); ++y)
                    std::copy_n(&raw.data[static_cast<std::size_t>(y + r.y0 - m.y0) * m.w() + (r.x0 - m.x0)], r.w(),
                                &out.data[static_cast<std::size_t>(y) * r.w()]);
                sub[static_cast<std::size_t>(i)] = std::move(out);
            } catch (const std::exception& e) {
#pragma omp critical
                error = e.what();
                failed = true;
            }
        }
        if (failed) throw std::runtime_error(error);
        return sub;
    };
    // a Prepared for a rectangle, with the given frame times
    auto region_prepared = [&](const Rect& r, const std::vector<double>& t) {
        Prepared q;
        q.width = r.w();
        q.height = r.h();
        q.nframes = n;
        q.t_ref_mjd = meta.t_ref_mjd;
        q.span_h = meta.span_h;
        q.t_h = t;
        q.background = meta.background;
        q.fwhm_px = meta.fwhm_px;
        q.exptime_s = meta.exptime_s;
        q.wcs = meta.wcs;
        if (q.wcs.valid) {
            q.wcs.crpix1 -= r.x0;
            q.wcs.crpix2 -= r.y0;
        }
        q.reg_dx = meta.reg_dx;  // detector = registered + reg; registered = region + origin
        q.reg_dy = meta.reg_dy;
        for (auto& v : q.reg_dx) v += r.x0;
        for (auto& v : q.reg_dy) v += r.y0;
        return q;
    };
    auto stage_region = [&](const Rect& region, const Rect& core) {
        StageRegion s;
        s.ox = region.x0;
        s.oy = region.y0;
        s.full_w = W;
        s.full_h = H;
        s.cx0 = core.x0 - region.x0;
        s.cy0 = core.y0 - region.y0;
        s.cx1 = core.x1 - region.x0;
        s.cy1 = core.y1 - region.y0;
        return s;
    };

    // ---- search setup
    SearchResult res;
    double dv = 0.0;
    const auto grid = velocity_grid(meta, params, &dv);
    res.n_velocities = grid.size();
    res.dv_pxh = dv;
    res.t_ref_mjd = meta.t_ref_mjd;
    const std::unique_ptr<SearchBackend> backend = make_backend(params.backend);
    res.backend = backend->name();
    res.backend_description = backend->description();
    const bool auto_thr = params.threshold_snr <= 0.0;
    const int runs = auto_thr ? std::max(0, params.null_runs) : 0;
    std::vector<std::vector<double>> orders{meta.t_h};
    for (auto& t : null_time_orders(meta.t_h, runs, params.seed)) orders.push_back(std::move(t));
    const int edge = std::max(params.edge_px, 1);
    const float floor = static_cast<float>(auto_thr ? params.min_snr_floor : params.threshold_snr);
    Stopwatch sw(&res.timings);

    // The tile size needs the PSF width; plan with a provisional width and redo
    // once it is known (only the halos depend on it).
    TilePlan plan = plan_tiles(src, params, opt, 3.0);
    auto tiles = [&]() {
        std::vector<Rect> v;
        for (int ty = 0; ty < H; ty += plan.tile)
            for (int tx = 0; tx < W; tx += plan.tile)
                v.push_back({tx, ty, std::min(W, tx + plan.tile), std::min(H, ty + plan.tile)});
        return v;
    };
    std::vector<Rect> cores = tiles();
    const double units = static_cast<double>(orders.size()) * 3.0 + 1.0;
    std::vector<Peak> peaks;
    std::vector<float> pooled;
    std::vector<int> cr_total(static_cast<std::size_t>(n), 0);
    std::vector<double> sigma_real, sg_real;

    for (std::size_t o = 0; o < orders.size(); ++o) {
        const std::vector<double>& t = orders[o];
        const bool real = o == 0;
        const std::string name = real ? "real" : "null " + std::to_string(o) + "/" + std::to_string(runs);
        const double base = static_cast<double>(o) * 3.0;

        // ---- A: residual noise per frame (and the static sky of the field for the PSF)
        std::vector<std::vector<float>> rs(static_cast<std::size_t>(n));
        Image static_full;
        if (real) static_full = Image(W, H, NAN);
        {
            std::vector<int> rank(static_cast<std::size_t>(n));
            std::vector<int> order(static_cast<std::size_t>(n));
            for (int i = 0; i < n; ++i) order[static_cast<std::size_t>(i)] = i;
            std::stable_sort(order.begin(), order.end(), [&](int a, int b) {
                return t[static_cast<std::size_t>(a)] < t[static_cast<std::size_t>(b)];
            });
            for (int r = 0; r < n; ++r) rank[static_cast<std::size_t>(order[static_cast<std::size_t>(r)])] = r;
            const int exclude = params.static_exclude >= 0 ? params.static_exclude
                                                           : std::max(1, static_cast<int>(std::lround(n / 8.0)));
            const SampleGrid grid2 = sample_grid(W, H, 200000, 2);
            for (std::size_t k = 0; k < cores.size(); ++k) {
                report((base + static_cast<double>(k) / cores.size()) / units, name + ": noise");
                const Rect& c = cores[k];
                const std::vector<Image> sub = load_sub(c);
                std::vector<std::vector<float>> local(static_cast<std::size_t>(nt) * n);
#pragma omp parallel num_threads(nt)
                {
                    std::vector<std::pair<float, int>> sorted;
                    std::vector<float> vals(static_cast<std::size_t>(n)), res_col(static_cast<std::size_t>(n));
                    auto* mine = &local[static_cast<std::size_t>(thread_id()) * n];
#pragma omp for schedule(static)
                    for (int y = 0; y < c.h(); ++y) {
                        const int gy = c.y0 + y;
                        const bool row_sampled = gy >= grid2.border && gy < H - grid2.border && (gy - grid2.border) % grid2.step == 0;
                        if (!real && !row_sampled) continue;
                        for (int x = 0; x < c.w(); ++x) {
                            const std::size_t idx = static_cast<std::size_t>(y) * c.w() + x;
                            const bool sampled = row_sampled && grid2.hit(c.x0 + x, gy);
                            if (sampled) {
                                const float s = static_sky_column(sub, idx, rank, exclude, sorted, vals, res_col.data());
                                if (real) static_full.at(c.x0 + x, gy) = s;
                                for (int i = 0; i < n; ++i)
                                    if (std::isfinite(res_col[static_cast<std::size_t>(i)]))
                                        mine[i].push_back(res_col[static_cast<std::size_t>(i)]);
                            } else if (real) {
                                static_full.at(c.x0 + x, gy) = static_sky_column(sub, idx, rank, exclude, sorted, vals, nullptr);
                            }
                        }
                    }
                }
                for (int th = 0; th < nt; ++th)
                    for (int i = 0; i < n; ++i) {
                        auto& from = local[static_cast<std::size_t>(th) * n + i];
                        rs[static_cast<std::size_t>(i)].insert(rs[static_cast<std::size_t>(i)].end(), from.begin(), from.end());
                    }
            }
        }
        std::vector<double> sigma(static_cast<std::size_t>(n));
        for (int i = 0; i < n; ++i) sigma[static_cast<std::size_t>(i)] = sigma_from_samples(rs[static_cast<std::size_t>(i)]);
        rs.clear();
        if (real) {
            meta.fwhm_px = estimate_psf_fwhm(static_full, sigma, params);
            static_full = Image();
            plan = plan_tiles(src, params, opt, meta.fwhm_px);
            cores = tiles();
            res.fwhm_px = meta.fwhm_px;
            res.slowest_reliable_pxh = 2.0 * meta.fwhm_px / meta.span_h;
        }
        release_free_memory();
        (real ? psw : sw).lap(real ? "noise (tiles, sampled)" : "null: noise (tiles, sampled)");

        // ---- B: matched-filter normalisation per frame
        std::vector<std::vector<float>> fs(static_cast<std::size_t>(n));
        {
            StageGlobals gl;
            gl.sigma = &sigma;
            for (std::size_t k = 0; k < cores.size(); ++k) {
                report((base + 1.0 + static_cast<double>(k) / cores.size()) / units, name + ": normalisation");
                const Rect r = grow_clip(cores[k], plan.filter_halo, W, H);
                std::vector<Image> sub = load_sub(r);
                Prepared q = region_prepared(r, t);
                const StageRegion sr = stage_region(r, cores[k]);
                StageSamples samples;
                build_stage(sub, q, params, {}, false, StageMode::FilterSamples, &sr, &gl, &samples);
                for (int i = 0; i < n; ++i) {
                    auto& from = samples.per_frame[static_cast<std::size_t>(i)];
                    fs[static_cast<std::size_t>(i)].insert(fs[static_cast<std::size_t>(i)].end(), from.begin(), from.end());
                }
            }
        }
        std::vector<double> filter_sg(static_cast<std::size_t>(n));
        for (int i = 0; i < n; ++i) filter_sg[static_cast<std::size_t>(i)] = robust_sigma(fs[static_cast<std::size_t>(i)]);
        fs.clear();
        release_free_memory();
        (real ? psw : sw).lap(real ? "normalisation (tiles)" : "null: normalisation (tiles)");
        if (real) {
            sigma_real = sigma;
            sg_real = filter_sg;
            meta.sigma = sigma;
        }

        // ---- C: prepare and stack every tile, assemble the detection maps
        StackResult full;
        full.max_snr = Image(W, H, 0.0f);
        full.best_v.assign(static_cast<std::size_t>(W) * H, -1);
        double t_prep = 0.0, t_stack = 0.0;
        {
            StageGlobals gl;
            gl.sigma = &sigma;
            gl.filter_sg = &filter_sg;
            for (std::size_t k = 0; k < cores.size(); ++k) {
                report((base + 2.0 + static_cast<double>(k) / cores.size()) / units, name + ": stack");
                const Rect& c = cores[k];
                const Rect r = grow_clip(c, plan.stack_halo + plan.filter_halo, W, H);
                const auto ta = std::chrono::steady_clock::now();
                Prepared q = region_prepared(r, t);
                {
                    std::vector<Image> sub = load_sub(r);
                    const StageRegion sr = stage_region(r, c);
                    build_stage(sub, q, params, {}, false, StageMode::Full, &sr, &gl, nullptr);
                }
                if (real)
                    for (int i = 0; i < n; ++i) cr_total[static_cast<std::size_t>(i)] += q.cr_pixels[static_cast<std::size_t>(i)];
                const auto tb = std::chrono::steady_clock::now();
                const PassOutput part = run_pass(*backend, q, grid, t, params, {}, 0.0, 1.0, name);
                for (int y = c.y0; y < c.y1; ++y) {
                    const std::size_t from = static_cast<std::size_t>(y - r.y0) * r.w() + (c.x0 - r.x0);
                    const std::size_t to = static_cast<std::size_t>(y) * W + c.x0;
                    std::copy_n(&part.max_snr.data[from], c.w(), &full.max_snr.data[to]);
                    std::copy_n(&part.best_v[from], c.w(), &full.best_v[to]);
                }
                t_prep += std::chrono::duration<double>(tb - ta).count();
                t_stack += std::chrono::duration<double>(std::chrono::steady_clock::now() - tb).count();
            }
        }
        sw.lap(real ? "tiles: prepare + stack (real)" : "tiles: prepare + stack (null)");
        res.timings.push_back({"  prepare (tiles)", t_prep});
        res.timings.push_back({"  stack (tiles)", t_stack});
        release_free_memory();

        if (real) {
            const int nms = std::max(2, static_cast<int>(std::lround(meta.fwhm_px)));
            peaks = find_peaks(full, floor, nms, edge, 100000);
            res.max_snr = std::move(full.max_snr);
        } else {
            const int nms = std::max(2, static_cast<int>(std::lround(meta.fwhm_px)));
            const std::vector<Peak> npk = find_peaks(full, 3.0f, nms, edge, 20000);
            // artefact test of the strongest null peaks, on tiles with the refinement halo
            std::vector<char> art(static_cast<std::size_t>(null_tested_count(params, npk.size())), 0);
            if (!art.empty()) {
                const Geometry g = geometry(meta, grid, dv, meta.fwhm_px);
                StageGlobals gl;
                gl.sigma = &sigma;
                gl.filter_sg = &filter_sg;
                for (const auto& c : cores) {
                    std::vector<int> mine;
                    for (std::size_t j = 0; j < art.size(); ++j)
                        if (npk[j].x >= c.x0 && npk[j].x < c.x1 && npk[j].y >= c.y0 && npk[j].y < c.y1)
                            mine.push_back(static_cast<int>(j));
                    if (mine.empty()) continue;
                    const Rect r = grow_clip(c, g.refine_halo, W, H);
                    Prepared q = region_prepared(r, t);
                    q.sub = load_sub(r);
                    const StageRegion sr = stage_region(r, c);
                    build_stage(q.sub, q, params, {}, false, StageMode::Full, &sr, &gl, nullptr);
#pragma omp parallel for num_threads(nt) schedule(dynamic)
                    for (int m = 0; m < static_cast<int>(mine.size()); ++m) {
                        const std::size_t j = static_cast<std::size_t>(mine[static_cast<std::size_t>(m)]);
                        Peak local = npk[j];
                        local.x -= r.x0;
                        local.y -= r.y0;
                        art[j] = is_artefact_peak(q, local, grid[static_cast<std::size_t>(npk[j].vidx)], params);
                    }
                }
            }
            pool_null_peaks(npk, art, pooled);
        }
        sw.lap(real ? "peaks (real)" : "peaks (null)");
    }
    meta.cr_pixels = cr_total;
    set_threshold(meta, params, grid.size(), runs, pooled, res);
    sw.lap("threshold fit");

    // ---- refine the candidates on tiles with a refinement halo
    report((units - 1.0) / units, "refine");
    const std::vector<Peak> keep = select_peaks(peaks, res.threshold_snr, params);
    std::vector<RefinedCandidate> rc(keep.size());
    {
        const Geometry g = geometry(meta, grid, dv, meta.fwhm_px);
        StageGlobals gl;
        gl.sigma = &sigma_real;
        gl.filter_sg = &sg_real;
        for (const auto& c : cores) {
            std::vector<int> mine;
            for (std::size_t k = 0; k < keep.size(); ++k)
                if (keep[k].x >= c.x0 && keep[k].x < c.x1 && keep[k].y >= c.y0 && keep[k].y < c.y1)
                    mine.push_back(static_cast<int>(k));
            if (mine.empty()) continue;
            const Rect r = grow_clip(c, g.refine_halo, W, H);
            Prepared q = region_prepared(r, meta.t_h);
            q.sub = load_sub(r);
            const StageRegion sr = stage_region(r, c);
            build_stage(q.sub, q, params, {}, false, StageMode::Full, &sr, &gl, nullptr);
#pragma omp parallel for num_threads(nt) schedule(dynamic)
            for (int m = 0; m < static_cast<int>(mine.size()); ++m) {
                const std::size_t k = static_cast<std::size_t>(mine[static_cast<std::size_t>(m)]);
                Peak local = keep[k];
                local.x -= r.x0;
                local.y -= r.y0;
                RefinedCandidate out = refine_with_contrib(q, local, grid[static_cast<std::size_t>(keep[k].vidx)], dv, params);
                out.c.x += r.x0;
                out.c.y += r.y0;
                rc[k] = std::move(out);
            }
        }
    }
    sw.lap("refine (tiles)");
    finish_candidates(meta, std::move(rc), res);
    sw.lap("ghosts + flags");
    report(1.0, "done");
    res.elapsed_s = std::chrono::duration<double>(std::chrono::steady_clock::now() - t_start).count();
    if (meta_out) *meta_out = std::move(meta);
    return res;
}

}  // namespace fainttrail
