// SPDX-License-Identifier: GPL-3.0-or-later
// Copyright 2026 Jurgen Kobierczynski
//
// Preprocessing: from registered frames to unit-variance, matched-filtered
// difference images in which only things that change or move remain.
#include <algorithm>
#include <chrono>
#include <atomic>
#include <cmath>
#include <numeric>
#include <stdexcept>
#include <string>

#include "common.hpp"
#include "internal.hpp"
#include "psf.hpp"
#include "fainttrail/registration.hpp"
#include "fainttrail/search.hpp"

namespace fainttrail {

using detail::kFwhmToSigma;

double estimate_fwhm(const Image& img, double noise_sigma) {
    if (img.empty() || noise_sigma <= 0.0) return 0.0;
    const int w = img.width, h = img.height;
    float vmax = -1e30f;
    for (float v : img.data)
        if (std::isfinite(v)) vmax = std::max(vmax, v);
    const double lo = 25.0 * noise_sigma, hi = 0.6 * vmax;
    const int R = 5, B = 9;
    std::vector<double> fwhms;
    for (int y = B; y < h - B && fwhms.size() < 400; ++y) {
        for (int x = B; x < w - B && fwhms.size() < 400; ++x) {
            const float p = img.at(x, y);
            if (!(p > lo && p < hi)) continue;
            bool is_max = true;
            for (int dy = -1; dy <= 1 && is_max; ++dy)
                for (int dx = -1; dx <= 1; ++dx)
                    if ((dx || dy) && img.at(x + dx, y + dy) >= p) { is_max = false; break; }
            if (!is_max) continue;
            // local background from an annulus, isolation test
            std::vector<float> ann;
            bool isolated = true;
            for (int dy = -B; dy <= B; ++dy)
                for (int dx = -B; dx <= B; ++dx) {
                    const double r = std::hypot(dx, dy);
                    const float v = img.at(x + dx, y + dy);
                    if (r > R + 1 && r <= B) ann.push_back(v);
                    // another star nearby: a brighter pixel, or a separate local maximum
                    if (r > 0 && v > p) isolated = false;
                    if (r > 2 && r <= B && v > 0.3f * p && x + dx > 0 && y + dy > 0 && x + dx < w - 1 &&
                        y + dy < h - 1) {
                        bool local = true;
                        for (int ey = -1; ey <= 1 && local; ++ey)
                            for (int ex = -1; ex <= 1; ++ex)
                                if ((ex || ey) && img.at(x + dx + ex, y + dy + ey) > v) { local = false; break; }
                        if (local) isolated = false;
                    }
                }
            if (!isolated) continue;
            const double bg = detail::median_inplace(ann);
            const double peak = p - bg;
            // hot pixels and other sharp defects have no PSF wings
            const double n4 = 0.25 * (img.at(x - 1, y) + img.at(x + 1, y) + img.at(x, y - 1) + img.at(x, y + 1)) - bg;
            if (n4 < 0.2 * peak) continue;
            // adaptive (Gaussian-weighted) second moments: for a Gaussian source the
            // weighted moment converges to sigma^2 / 2 when the weight matches the source
            double sw2 = 1.5 * 1.5;
            bool good = true;
            for (int it = 0; it < 12 && good; ++it) {
                double s = 0, sxx = 0;
                for (int dy = -B; dy <= B; ++dy)
                    for (int dx = -B; dx <= B; ++dx) {
                        const double r2 = dx * dx + dy * dy;
                        const double wgt = std::exp(-0.5 * r2 / sw2);
                        const double v = (img.at(x + dx, y + dy) - bg) * wgt;
                        s += v;
                        sxx += v * r2;
                    }
                if (s <= 0) { good = false; break; }
                const double m = 0.5 * sxx / s;  // per-axis weighted second moment
                const double next = std::clamp(2.0 * m, 0.25, 49.0);
                if (std::fabs(next - sw2) < 1e-4) { sw2 = next; break; }
                sw2 = next;
            }
            if (!good) continue;
            const double var = sw2 - 1.0 / 12.0;  // remove pixel-integration width
            if (var <= 0.1) continue;
            fwhms.push_back(2.354820045 * std::sqrt(var));
        }
    }
    if (fwhms.size() < 5) return 0.0;
    std::vector<float> f(fwhms.begin(), fwhms.end());
    // the truncated moment window underestimates wide PSFs slightly; that is acceptable
    return std::clamp(detail::median_inplace(f), 1.0, 20.0);
}

namespace {

// Separable Gaussian matched filter normalised to unit noise variance:
// out = (z * K) / sqrt(sum K^2) with K = k k^T, so sqrt(sum K^2) = sum k^2.
void matched_filter(const Image& in, Image& out, double sigma, Image& tmp) {
    const int w = in.width, h = in.height;
    const int R = std::max(1, static_cast<int>(std::ceil(3.0 * sigma)));
    std::vector<float> k(static_cast<std::size_t>(2 * R + 1));
    double k2 = 0.0;
    for (int i = -R; i <= R; ++i) {
        const double v = std::exp(-0.5 * i * i / (sigma * sigma));
        k[static_cast<std::size_t>(i + R)] = static_cast<float>(v);
        k2 += v * v;
    }
    const float norm = static_cast<float>(1.0 / k2);
    tmp = Image(w, h);
    out = Image(w, h);
    for (int y = 0; y < h; ++y) {
        const float* row = &in.data[static_cast<std::size_t>(y) * w];
        float* t = &tmp.data[static_cast<std::size_t>(y) * w];
        for (int x = 0; x < w; ++x) {
            const int a = std::max(0, x - R), b = std::min(w - 1, x + R);
            float acc = 0.0f;
            for (int j = a; j <= b; ++j) acc += row[j] * k[static_cast<std::size_t>(j - x + R)];
            t[x] = acc;
        }
    }
    for (int y = 0; y < h; ++y) {
        const int a = std::max(0, y - R), b = std::min(h - 1, y + R);
        float* o = &out.data[static_cast<std::size_t>(y) * w];
        for (int j = a; j <= b; ++j) {
            const float kj = k[static_cast<std::size_t>(j - y + R)] * norm;
            const float* t = &tmp.data[static_cast<std::size_t>(j) * w];
            for (int x = 0; x < w; ++x) o[x] += kj * t[x];
        }
    }
}

}  // namespace

namespace detail {

float static_sky_column(const std::vector<Image>& sub, std::size_t idx, const std::vector<int>& rank, int exclude,
                        std::vector<std::pair<float, int>>& sorted, std::vector<float>& vals, float* residual) {
    const int n = static_cast<int>(sub.size());
    sorted.clear();
    for (int i = 0; i < n; ++i) {
        const float v = sub[static_cast<std::size_t>(i)].data[idx];
        vals[static_cast<std::size_t>(i)] = v;
        if (std::isfinite(v)) sorted.emplace_back(v, rank[static_cast<std::size_t>(i)]);
    }
    // too few frames cover this pixel for a reliable static sky (field edges)
    if (sorted.size() < std::max<std::size_t>(5, static_cast<std::size_t>(0.6 * n))) {
        if (residual)
            for (int i = 0; i < n; ++i) residual[i] = NAN;
        return NAN;
    }
    std::sort(sorted.begin(), sorted.end());
    const std::size_t m_all = sorted.size();
    const double global = 0.5 * (sorted[(m_all - 1) / 2].first + sorted[m_all / 2].first);
    if (!residual) return static_cast<float>(global);
    for (int i = 0; i < n; ++i) {
        const float v = vals[static_cast<std::size_t>(i)];
        if (!std::isfinite(v)) {
            residual[i] = NAN;
            continue;
        }
        double st = global;
        if (exclude > 0) {
            const int ri = rank[static_cast<std::size_t>(i)];
            std::size_t m = 0;
            for (const auto& e : sorted)
                if (std::abs(e.second - ri) > exclude) ++m;
            if (m >= 5) {
                const std::size_t lo = (m - 1) / 2, hi = m / 2;
                std::size_t k = 0;
                double a = 0.0, b = 0.0;
                for (const auto& e : sorted) {
                    if (std::abs(e.second - ri) <= exclude) continue;
                    if (k == lo) a = e.first;
                    if (k == hi) {
                        b = e.first;
                        break;
                    }
                    ++k;
                }
                st = 0.5 * (a + b);
            }
        }
        residual[i] = static_cast<float>(v - st);
    }
    return static_cast<float>(global);
}

Image residual_patch(const Prepared& p, int frame, int x0, int y0, int w, int h) {
    Image out(w, h, NAN);
    if (static_cast<int>(p.residual.size()) == p.nframes) {  // kept: just copy
        const Image& r = p.residual[static_cast<std::size_t>(frame)];
        for (int y = 0; y < h; ++y)
            for (int x = 0; x < w; ++x)
                if (r.inside(x0 + x, y0 + y)) out.at(x, y) = r.at(x0 + x, y0 + y);
        return out;
    }
    const std::vector<Image>& subs = sub_frames(p);
    if (static_cast<int>(subs.size()) != p.nframes || static_cast<int>(p.time_rank.size()) != p.nframes)
        throw std::runtime_error("residuals are neither kept nor recomputable (no background-subtracted frames)");
    std::vector<std::pair<float, int>> sorted;
    std::vector<float> vals(static_cast<std::size_t>(p.nframes)), res(static_cast<std::size_t>(p.nframes));
    for (int y = 0; y < h; ++y)
        for (int x = 0; x < w; ++x) {
            if (x0 + x < 0 || y0 + y < 0 || x0 + x >= p.width || y0 + y >= p.height) continue;
            const std::size_t idx = static_cast<std::size_t>(y0 + y) * p.width + (x0 + x);
            static_sky_column(subs, idx, p.time_rank, p.static_exclude, sorted, vals, res.data());
            out.at(x, y) = res[static_cast<std::size_t>(frame)];
        }
    return out;
}

// Steps that depend on the time order of the frames: static sky, residuals,
// noise, mask and matched filter. Used for the real data and, with permuted
// times, for every null run, so that both see exactly the same estimator.
void subtract_and_grow(Image& im, double bg) {
    const int w = im.width, h = im.height;
    for (auto& v : im.data) v = std::isfinite(v) ? static_cast<float>(v - bg) : NAN;
    // Registered frames carry blank (NaN) borders where the field drifted out.
    // Interpolation leaves damaged pixels along those borders, so the blank area
    // is grown by a few pixels.
    const int grow = kBlankGrow;
    std::vector<unsigned char> bad(im.data.size(), 0);
    bool any = false;
    for (std::size_t k = 0; k < im.data.size(); ++k)
        if (!std::isfinite(im.data[k])) { bad[k] = 1; any = true; }
    if (!any) return;
    Image grown = im;
    for (int y = 0; y < h; ++y)
        for (int x = 0; x < w; ++x) {
            if (!bad[static_cast<std::size_t>(y) * w + x]) continue;
            for (int dy = -grow; dy <= grow; ++dy)
                for (int dx = -grow; dx <= grow; ++dx)
                    if (grown.inside(x + dx, y + dy)) grown.at(x + dx, y + dy) = NAN;
        }
    im = std::move(grown);
}

double estimate_psf_fwhm(const Image& static_sky, std::vector<double> sg, const SearchParams& params) {
    const int n = static_cast<int>(sg.size());
    std::nth_element(sg.begin(), sg.begin() + n / 2, sg.end());
    const double static_noise = sg[static_cast<std::size_t>(n / 2)] / std::sqrt(static_cast<double>(n));
    double f = params.fwhm_px > 0.0 ? params.fwhm_px : estimate_fwhm(static_sky, static_noise);
    if (!(f > 0.0)) f = 3.0;  // sensible default when no stars are found
    return f;
}

double sigma_from_samples(std::vector<float>& s) {
    const double sg = detail::robust_sigma(s);
    return sg > 0.0 ? sg : 1.0;
}

bool SampleGrid::hit(int gx, int gy) const {
    if (gx < border || gy < border || gx >= full_w - border || gy >= full_h - border) return false;
    return (gx - border) % step == 0 && (gy - border) % step == 0;
}

SampleGrid sample_grid(int full_w, int full_h, std::size_t max_samples, int border) {
    SampleGrid g;
    g.full_w = full_w;
    g.full_h = full_h;
    g.border = border;
    const std::size_t usable = static_cast<std::size_t>(std::max(0, full_w - 2 * border)) *
                               static_cast<std::size_t>(std::max(0, full_h - 2 * border));
    g.step = usable == 0 ? 1 : std::max(1, static_cast<int>(std::sqrt(static_cast<double>(usable) / max_samples)));
    return g;
}

void build_stage(const std::vector<Image>& sub, Prepared& p, const SearchParams& params,
                 const ProgressFn& progress, bool estimate_psf, StageMode mode, const StageRegion* region,
                 const StageGlobals* globals, StageSamples* samples) {
    const int n = p.nframes, w = p.width, h = p.height;
    // A region of a larger image (tiled search): per-frame quantities measured on the
    // whole image come from `globals`, samples for them are collected over the core
    // only, at the same positions sample_pixels() would use on the whole image.
    const bool regional = region != nullptr;
    const int ox = regional ? region->ox : 0, oy = regional ? region->oy : 0;
    const int fw = regional ? region->full_w : w, fh = regional ? region->full_h : h;
    const int cx0 = regional ? region->cx0 : 0, cy0 = regional ? region->cy0 : 0;
    const int cx1 = regional ? region->cx1 : w, cy1 = regional ? region->cy1 : h;
    auto collect = [&](const Image& img, const SampleGrid& grid, const Image* mask, std::vector<float>& out) {
        const int st = grid.step, b = grid.border;
        for (int y = cy0; y < cy1; ++y) {
            const int gy = oy + y;
            if (gy < b || gy >= grid.full_h - b || (gy - b) % st != 0) continue;
            int gx = std::max(b, ox + cx0);
            gx += ((b - gx) % st + st) % st;  // first sampled column at or after the core start
            for (; gx < std::min(ox + cx1, grid.full_w - b); gx += st) {
                const int x = gx - ox;
                if (mask && mask->at(x, y) == 0.0f) continue;
                const float v = img.at(x, y);
                if (std::isfinite(v)) out.push_back(v);
            }
        }
    };
    auto report = [&](double f, const char* stage) {
        if (progress && !progress(f, stage)) throw Cancelled();
    };
    const int nt = detail::thread_count(params.threads);
    Stopwatch sw(&p.timings);

    // static sky = per-pixel median over time. For the residual of frame i the
    // median excludes the frames closest in time to i, so an object that moves
    // more than about a PSF width in `exclude` frames does not subtract itself.
    report(0.1, "static sky");
    std::vector<int> rank(static_cast<std::size_t>(n));
    {
        std::vector<int> order(static_cast<std::size_t>(n));
        std::iota(order.begin(), order.end(), 0);
        std::stable_sort(order.begin(), order.end(), [&](int a, int b) {
            return p.t_h[static_cast<std::size_t>(a)] < p.t_h[static_cast<std::size_t>(b)];
        });
        for (int r = 0; r < n; ++r) rank[static_cast<std::size_t>(order[static_cast<std::size_t>(r)])] = r;
    }
    const int exclude = params.static_exclude >= 0 ? params.static_exclude
                                                   : std::max(1, static_cast<int>(std::lround(n / 8.0)));
    p.static_sky = Image(w, h);
    p.residual.assign(static_cast<std::size_t>(n), Image(w, h));
#pragma omp parallel num_threads(nt)
    {
        std::vector<std::pair<float, int>> sorted;
        std::vector<float> vals(static_cast<std::size_t>(n)), res(static_cast<std::size_t>(n));
        sorted.reserve(static_cast<std::size_t>(n));
#pragma omp for schedule(static)
        for (int y = 0; y < h; ++y) {
            for (int x = 0; x < w; ++x) {
                const std::size_t idx = static_cast<std::size_t>(y) * w + x;
                p.static_sky.data[idx] = static_sky_column(sub, idx, rank, exclude, sorted, vals, res.data());
                for (int i = 0; i < n; ++i) p.residual[static_cast<std::size_t>(i)].data[idx] = res[static_cast<std::size_t>(i)];
            }
        }
    }
    p.time_rank = rank;
    p.static_exclude = exclude;

    sw.lap("static sky + residuals");
    // per-frame residual noise
    report(0.3, "residuals");
    if (globals && globals->sigma) {
        p.sigma = *globals->sigma;
    } else {
        p.sigma.assign(static_cast<std::size_t>(n), 1.0);
#pragma omp parallel for num_threads(nt) schedule(dynamic)
        for (int i = 0; i < n; ++i) {
            auto s = detail::sample_pixels(p.residual[static_cast<std::size_t>(i)], 200000, nullptr, 2);
            p.sigma[static_cast<std::size_t>(i)] = sigma_from_samples(s);
        }
    }

    // PSF width
    if (estimate_psf) p.fwhm_px = estimate_psf_fwhm(p.static_sky, p.sigma, params);

    sw.lap("noise + PSF");
    // temporal scatter map -> weights and mask
    report(0.45, "mask");
    Image scatter(w, h);
#pragma omp parallel num_threads(nt)
    {
        std::vector<float> col;
        col.reserve(static_cast<std::size_t>(n));
#pragma omp for schedule(static)
        for (int y = 0; y < h; ++y)
            for (int x = 0; x < w; ++x) {
                col.clear();
                const std::size_t idx = static_cast<std::size_t>(y) * w + x;
                for (int i = 0; i < n; ++i) {
                    const float v = p.residual[static_cast<std::size_t>(i)].data[idx];
                    if (std::isfinite(v)) col.push_back(static_cast<float>(std::fabs(v) / p.sigma[static_cast<std::size_t>(i)]));
                }
                scatter.data[idx] = col.size() < static_cast<std::size_t>(n) / 2
                                        ? 1e9f
                                        : static_cast<float>(1.4826 * detail::median_inplace(col));
            }
    }
    p.mask = Image(w, h, 1.0f);
    const int dil = std::max(1, static_cast<int>(std::lround(0.5 * p.fwhm_px)));
    for (int y = 0; y < h; ++y)
        for (int x = 0; x < w; ++x) {
            if (scatter.at(x, y) <= params.mask_scatter) continue;
            for (int dy = -dil; dy <= dil; ++dy)
                for (int dx = -dil; dx <= dil; ++dx)
                    if (p.mask.inside(x + dx, y + dy)) p.mask.at(x + dx, y + dy) = 0.0f;
        }

    sw.lap("scatter mask");
    // normalise, weight, reject cosmic rays, clip, matched filter, renormalise
    report(0.55, "matched filter");
    const double psf_sigma = p.fwhm_px * kFwhmToSigma;
    const float clip = static_cast<float>(params.clip_sigma);
    const float wscat = static_cast<float>(params.weight_scatter);
    p.filtered.assign(static_cast<std::size_t>(n), Image());
    p.cr_pixels.assign(static_cast<std::size_t>(n), 0);
    if (mode == StageMode::FilterSamples) samples->per_frame.assign(static_cast<std::size_t>(n), {});
    std::atomic<int> done{0};
    std::atomic<bool> cancel{false};
#pragma omp parallel for num_threads(nt) schedule(dynamic)
    for (int i = 0; i < n; ++i) {
        if (cancel.load()) continue;
        const Image& r = p.residual[static_cast<std::size_t>(i)];
        Image z(w, h);
        const float inv = static_cast<float>(1.0 / p.sigma[static_cast<std::size_t>(i)]);
        for (std::size_t k = 0; k < z.data.size(); ++k) {
            const float v = r.data[k];
            if (p.mask.data[k] == 0.0f || !std::isfinite(v)) continue;
            // mild excess scatter (e.g. along the path of a bright mover) is not penalised
            const float e = scatter.data[k] > wscat ? scatter.data[k] : 1.0f;
            z.data[k] = v * inv / e;
        }
        // Cosmic rays and other single-frame defects are sharper than the PSF:
        // their direct neighbours are much fainter than a stellar profile allows.
        int cr = 0;
        {
            const float ratio = static_cast<float>(0.5 * std::exp(-0.5 / (psf_sigma * psf_sigma)));
            const float hit = static_cast<float>(params.cr_sigma);
            std::vector<std::size_t> kill;
            for (int y = 1; y < h - 1; ++y)
                for (int x = 1; x < w - 1; ++x) {
                    const std::size_t k = static_cast<std::size_t>(y) * w + x;
                    const float v = z.data[k];
                    if (v < hit) continue;
                    const float n4 = 0.25f * (z.data[k - 1] + z.data[k + 1] + z.data[k - w] + z.data[k + w]);
                    if (n4 >= ratio * v) continue;
                    if (x >= cx0 && x < cx1 && y >= cy0 && y < cy1) ++cr;
                    // remove the hit and any bright pixels touching it
                    for (int dy = -1; dy <= 1; ++dy)
                        for (int dx = -1; dx <= 1; ++dx) {
                            const std::size_t q = static_cast<std::size_t>(y + dy) * w + (x + dx);
                            if ((dx == 0 && dy == 0) || z.data[q] > 2.5f) kill.push_back(q);
                        }
                }
            for (std::size_t q : kill) z.data[q] = 0.0f;
        }
        p.cr_pixels[static_cast<std::size_t>(i)] = cr;
        for (auto& v : z.data) v = std::clamp(v, -clip, clip);
        Image g, tmp;
        matched_filter(z, g, psf_sigma, tmp);
        const int fborder = static_cast<int>(std::ceil(3 * psf_sigma)) + 1;
        if (mode == StageMode::FilterSamples) {
            collect(g, sample_grid(fw, fh, 200000, fborder), &p.mask, samples->per_frame[static_cast<std::size_t>(i)]);
            if (!params.keep_residuals) p.residual[static_cast<std::size_t>(i)] = Image();
            ++done;
            continue;
        }
        double sg = 0.0;
        if (globals && globals->filter_sg) {
            sg = (*globals->filter_sg)[static_cast<std::size_t>(i)];
        } else {
            auto s = detail::sample_pixels(g, 200000, &p.mask, fborder);
            sg = detail::robust_sigma(s);
        }
        if (sg > 0.0) {
            const float ig = static_cast<float>(1.0 / sg);
            for (auto& v : g.data) v *= ig;
        }
        p.filtered[static_cast<std::size_t>(i)] = std::move(g);
        // the residual of this frame is no longer needed (the scatter map is done)
        if (!params.keep_residuals) p.residual[static_cast<std::size_t>(i)] = Image();
        const int d = ++done;
        if (detail::thread_id() == 0 && progress && !progress(0.55 + 0.45 * d / n, "matched filter"))
            cancel = true;
    }
    if (cancel) throw Cancelled();
    if (!params.keep_residuals) {
        p.residual.clear();
        p.residual.shrink_to_fit();
    }
    release_free_memory();
    sw.lap("matched filter");
}

Prepared with_times(const Prepared& base, const std::vector<double>& t_h, const SearchParams& params) {
    if (base.sub.empty()) throw std::runtime_error("prepared data has no input frames for the null runs");
    Prepared q;
    q.width = base.width;
    q.height = base.height;
    q.nframes = base.nframes;
    q.t_ref_mjd = base.t_ref_mjd;
    q.span_h = base.span_h;
    q.t_h = t_h;
    q.background = base.background;
    q.fwhm_px = base.fwhm_px;
    q.wcs = base.wcs;
    q.exptime_s = base.exptime_s;
    q.reg_dx = base.reg_dx;
    q.reg_dy = base.reg_dy;
    SearchParams np = params;
    np.keep_residuals = false;  // the null runs only need the filtered frames
    build_stage(base.sub, q, np, {}, false);
    q.shared_sub = &sub_frames(base);
    return q;
}

}  // namespace detail

namespace {

// steal = true: the frames' pixel buffers are moved into the prepared data.
Prepared prepare_core(std::vector<Frame>& frames, bool steal, const SearchParams& params,
                      const ProgressFn& progress) {
    const int n = static_cast<int>(frames.size());
    if (n < 5) throw std::runtime_error("synthetic tracking needs at least 5 frames");
    const int w = frames[0].image.width, h = frames[0].image.height;
    if (w < 16 || h < 16) throw std::runtime_error("frames are too small");
    for (const auto& f : frames)
        if (f.image.width != w || f.image.height != h)
            throw std::runtime_error("all frames must have the same size (are they registered?)");
    if (progress && !progress(0.0, "background")) throw Cancelled();
    const int nt = detail::thread_count(params.threads);

    Prepared p;
    p.width = w;
    p.height = h;
    p.nframes = n;
    double tmin = frames[0].mjd_mid, tmax = tmin;
    for (const auto& f : frames) {
        tmin = std::min(tmin, f.mjd_mid);
        tmax = std::max(tmax, f.mjd_mid);
    }
    p.t_ref_mjd = 0.5 * (tmin + tmax);
    p.span_h = (tmax - tmin) * 24.0;
    if (p.span_h <= 0.0) throw std::runtime_error("frames must span a non-zero time interval");
    for (const auto& f : frames) {
        p.t_h.push_back((f.mjd_mid - p.t_ref_mjd) * 24.0);
        p.exptime_s.push_back(f.exptime_s);
    }
    for (const auto& f : frames)
        if (f.wcs.valid) { p.wcs = f.wcs; break; }
    {
        bool all = true;
        for (const auto& f : frames) all = all && f.has_reg_offset;
        if (all)
            for (const auto& f : frames) {
                p.reg_dx.push_back(f.reg_dx);
                p.reg_dy.push_back(f.reg_dy);
            }
    }

    // per-frame background level
    detail::Stopwatch sw(&p.timings);
    p.background.resize(static_cast<std::size_t>(n));
    p.sub.resize(static_cast<std::size_t>(n));
#pragma omp parallel for num_threads(nt) schedule(dynamic)
    for (int i = 0; i < n; ++i) {
        auto s = detail::sample_pixels(frames[static_cast<std::size_t>(i)].image, 200000);
        const double bg = detail::median_inplace(s);
        p.background[static_cast<std::size_t>(i)] = bg;
        Image im = steal ? std::move(frames[static_cast<std::size_t>(i)].image) : frames[static_cast<std::size_t>(i)].image;
        detail::subtract_and_grow(im, bg);
        p.sub[static_cast<std::size_t>(i)] = std::move(im);
    }
    sw.lap("background + borders");
    detail::build_stage(p.sub, p, params, progress, true);
    if (progress && !progress(1.0, "prepared")) throw Cancelled();
    return p;
}

// Registration check and, when needed, alignment (SearchParams::align), then the
// preprocessing proper.
Prepared prepare_frames(std::vector<Frame>& frames, bool steal, const SearchParams& params,
                        const ProgressFn& progress) {
    const std::string& mode = params.align;
    if (mode != "auto" && mode != "always" && mode != "check" && mode != "off")
        throw std::runtime_error("align must be auto, always, check or off (not '" + mode + "')");
    if (frames.size() < 5) throw std::runtime_error("synthetic tracking needs at least 5 frames");
    using clk = std::chrono::steady_clock;
    DriftCheck chk;
    AlignmentReport rep;
    std::vector<StageTime> times;
    std::vector<Frame> local;
    std::vector<Frame>* src = &frames;
    bool st = steal;
    bool same = true;
    for (const auto& f : frames)
        same = same && f.image.width == frames[0].image.width && f.image.height == frames[0].image.height;
    if (mode != "off") {
        const auto t0 = clk::now();
        if (same) {
            std::vector<const Image*> ims;
            std::vector<double> t;
            for (const auto& f : frames) {
                ims.push_back(&f.image);
                t.push_back((f.mjd_mid - frames[0].mjd_mid) * 24.0);
            }
            chk = check_registration(ims, t);
        }
        times.push_back({"registration check", std::chrono::duration<double>(clk::now() - t0).count()});
        const bool problem = !same || !chk.registered;
        if (mode == "check" && problem)
            throw std::runtime_error(same ? chk.message : "the frames differ in size: they are not registered on a common grid");
        if (mode == "always" || (mode == "auto" && problem)) {
            if (progress && !progress(0.0, "registering frames")) throw Cancelled();
            const auto t1 = clk::now();
            if (!steal) {  // the caller's frames stay untouched
                local = frames;
                src = &local;
                st = true;
            }
            rep = register_frames(*src, params.threads);
            times.push_back({"register frames", std::chrono::duration<double>(clk::now() - t1).count()});
        }
    } else if (!same) {
        throw std::runtime_error("all frames must have the same size (are they registered?)");
    }
    Prepared p = prepare_core(*src, st, params, progress);
    p.registration = chk;
    p.alignment = std::move(rep);
    p.timings.insert(p.timings.begin(), times.begin(), times.end());
    return p;
}

}  // namespace

Prepared prepare(const std::vector<Frame>& frames, const SearchParams& params, const ProgressFn& progress) {
    // not modified: steal = false only reads the frames
    return prepare_frames(const_cast<std::vector<Frame>&>(frames), false, params, progress);
}

Prepared prepare(std::vector<Frame>&& frames, const SearchParams& params, const ProgressFn& progress) {
    return prepare_frames(frames, true, params, progress);
}

}  // namespace fainttrail
