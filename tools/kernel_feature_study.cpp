// SPDX-License-Identifier: GPL-3.0-or-later
// Copyright 2026 Jurgen Kobierczynski
//
// Which features should the stacking kernel have? (docs/KERNEL_FEATURES.md)
//
//   fainttrail_kernel_feature_study shifts      whole-pixel vs sub-pixel shifts, grid spacing
//   fainttrail_kernel_feature_study curvature   S/N lost by straight-line stacking of curved tracks
//   fainttrail_kernel_feature_study weights DIR... [--wcs-from F]
//                                               per-frame efficiency in real sequences
//   fainttrail_kernel_feature_study weights-sim  end-to-end gain of per-frame weights
//
// All comparisons that involve a threshold run the full null calibration with the
// same method, so detections are compared at equal false-alarm rates.
#include <algorithm>
#include <atomic>
#include <chrono>
#include <cmath>
#include <cstdio>
#include <cstdlib>
#include <filesystem>
#include <limits>
#include <numeric>
#include <string>
#include <vector>

#include "../src/common.hpp"
#include "../src/internal.hpp"
#include "../src/psf.hpp"
#include "fainttrail/fainttrail.hpp"

using namespace fainttrail;
namespace fs = std::filesystem;
using clk = std::chrono::steady_clock;

namespace {

double secs(clk::time_point t) { return std::chrono::duration<double>(clk::now() - t).count(); }

// ---------------------------------------------------------------- sub-pixel kernel
// Experimental backend: shifts v*t are not rounded; every frame is sampled with
// bilinear interpolation. A frame covers a pixel when all four samples are inside.
// Straightforward (not cache-blocked) implementation, for accuracy, not speed.
class SubpixelBackend final : public SearchBackend {
public:
    std::string name() const override { return "cpu-subpixel"; }
    std::string description() const override { return "experimental: bilinear sub-pixel shifts"; }
    StackResult run(const StackPass& pass, const ProgressFn& = {}) override {
        const auto& F = *pass.frames;
        const auto& t = *pass.t_h;
        const auto& V = *pass.velocities;
        const int w = F[0].width, h = F[0].height, n = static_cast<int>(F.size()), nv = static_cast<int>(V.size());
        const int nt = detail::thread_count(pass.threads);
        std::vector<Image> lmax(static_cast<std::size_t>(nt));
        std::vector<std::vector<int>> larg(static_cast<std::size_t>(nt));
#pragma omp parallel num_threads(nt)
        {
            const int tid = detail::thread_id();
            Image& M = lmax[static_cast<std::size_t>(tid)];
            auto& A = larg[static_cast<std::size_t>(tid)];
            M = Image(w, h, -std::numeric_limits<float>::infinity());
            A.assign(static_cast<std::size_t>(w) * h, -1);
            std::vector<float> acc(static_cast<std::size_t>(w) * h);
            std::vector<int> cov(static_cast<std::size_t>(w) * h);
#pragma omp for schedule(dynamic, 2)
            for (int vi = 0; vi < nv; ++vi) {
                std::fill(acc.begin(), acc.end(), 0.0f);
                std::fill(cov.begin(), cov.end(), 0);
                for (int i = 0; i < n; ++i) {
                    const double sx = V[static_cast<std::size_t>(vi)].first * t[static_cast<std::size_t>(i)];
                    const double sy = V[static_cast<std::size_t>(vi)].second * t[static_cast<std::size_t>(i)];
                    const int ix = static_cast<int>(std::floor(sx)), iy = static_cast<int>(std::floor(sy));
                    const float fx = static_cast<float>(sx - ix), fy = static_cast<float>(sy - iy);
                    const float w00 = (1 - fx) * (1 - fy), w01 = fx * (1 - fy), w10 = (1 - fx) * fy, w11 = fx * fy;
                    const int x0 = std::max(0, -ix), x1 = std::min(w, w - ix - 1);
                    const int y0 = std::max(0, -iy), y1 = std::min(h, h - iy - 1);
                    const float* g = F[static_cast<std::size_t>(i)].data.data();
                    for (int y = y0; y < y1; ++y) {
                        const float* r0 = g + static_cast<std::ptrdiff_t>(y + iy) * w + ix;
                        const float* r1 = r0 + w;
                        float* o = &acc[static_cast<std::size_t>(y) * w];
                        int* c = &cov[static_cast<std::size_t>(y) * w];
                        for (int x = x0; x < x1; ++x) {
                            o[x] += w00 * r0[x] + w01 * r0[x + 1] + w10 * r1[x] + w11 * r1[x + 1];
                            ++c[x];
                        }
                    }
                }
                for (std::size_t q = 0; q < acc.size(); ++q) {
                    if (cov[q] < pass.min_count) continue;
                    const float s = acc[q] / std::sqrt(static_cast<float>(cov[q]));
                    if (s > M.data[q]) {
                        M.data[q] = s;
                        A[q] = vi;
                    }
                }
            }
        }
        StackResult out;
        out.max_snr = Image(w, h, 0.0f);
        out.best_v.assign(static_cast<std::size_t>(w) * h, -1);
        for (int k = 0; k < nt; ++k)
            for (std::size_t q = 0; q < out.max_snr.data.size(); ++q) {
                const int a = larg[static_cast<std::size_t>(k)][q];
                if (a < 0) continue;
                const float v = lmax[static_cast<std::size_t>(k)].data[q];
                if (out.best_v[q] < 0 || v > out.max_snr.data[q] || (v == out.max_snr.data[q] && a < out.best_v[q])) {
                    out.max_snr.data[q] = v;
                    out.best_v[q] = a;
                }
            }
        return out;
    }
};

// the strongest detection-map value within r pixels of (x, y)
float map_peak(const Image& M, double x, double y, int r = 2) {
    float best = 0.0f;
    const int cx = static_cast<int>(std::lround(x)), cy = static_cast<int>(std::lround(y));
    for (int dy = -r; dy <= r; ++dy)
        for (int dx = -r; dx <= r; ++dx)
            if (M.inside(cx + dx, cy + dy)) best = std::max(best, M.at(cx + dx, cy + dy));
    return best;
}

SimParams base_sim(std::uint64_t seed, int size, double fwhm) {
    SimParams s;
    s.width = size;
    s.height = size;
    s.nframes = 60;
    s.exptime_s = 45;
    s.cadence_s = 60;
    s.fwhm_px = fwhm;
    s.nstars = static_cast<int>(800.0 * size * size / (512.0 * 512.0));
    s.hot_pixels = static_cast<int>(50.0 * size * size / (512.0 * 512.0));
    s.seed = seed;
    return s;
}

// ---------------------------------------------------------------- 1. shifts
int study_shifts() {
    register_backend("cpu-subpixel", [] { return std::make_unique<SubpixelBackend>(); });
    struct Variant {
        const char* label;
        const char* backend;
        double grid_tol;
    };
    const Variant variants[] = {{"whole-pixel, grid 0.5 (current)", "cpu", 0.5},
                                {"whole-pixel, grid 0.25", "cpu", 0.25},
                                {"sub-pixel,   grid 0.5", "cpu-subpixel", 0.5},
                                {"sub-pixel,   grid 1.0", "cpu-subpixel", 1.0}};
    std::printf("%-6s %-32s %7s %9s %10s %10s %9s %10s\n", "FWHM", "kernel", "motions", "threshold", "S/N ratio",
                "S/N / thr", "found", "s / pass");
    for (double fwhm : {1.6, 3.0}) {
        std::vector<SimResult> sims;
        for (int k = 0; k < 2; ++k) {
            SimParams s = base_sim(300 + static_cast<std::uint64_t>(k) + static_cast<std::uint64_t>(fwhm * 10), 384, fwhm);
            s.nmovers = 40;
            s.rate_min_arcsec_min = 0.15;
            s.rate_max_arcsec_min = 0.7;
            // movers near the detection limit: about S/N 6-12 in the stack
            s.mover_mag_bright = fwhm < 2 ? 20.9 : 20.6;
            s.mover_mag_faint = fwhm < 2 ? 21.7 : 21.4;
            sims.push_back(simulate(s));
        }
        for (const auto& v : variants) {
            SearchParams p;
            p.vmin_pxh = 5;
            p.vmax_pxh = 30;
            p.grid_tol_px = v.grid_tol;
            p.backend = v.backend;
            double ratio_sum = 0, rel_sum = 0, thr_sum = 0, t_pass = 0;
            int nobj = 0, found = 0;
            std::size_t nv = 0;
            for (const auto& sim : sims) {
                const Prepared prep = prepare(sim.frames, p);
                const auto t0 = clk::now();
                const SearchResult r = search(prep, p);
                t_pass += secs(t0) / (1 + p.null_runs) / sims.size();
                nv = r.n_velocities;
                thr_sum += r.threshold_snr / sims.size();
                std::vector<SimMover> movers;
                for (const auto& tr : sim.truth) movers.push_back(tr.mover);
                const auto m = match(r.candidates, movers, 2.0, std::max(3.0, 2 * r.dv_pxh), true);
                found += static_cast<int>(m.size());
                for (const auto& tr : sim.truth) {
                    const float s = map_peak(r.max_snr, tr.mover.x_ref, tr.mover.y_ref);
                    ratio_sum += s / tr.expected_snr_stack;
                    rel_sum += s / r.threshold_snr;
                    ++nobj;
                }
            }
            std::printf("%-6.1f %-32s %7zu %9.2f %10.3f %10.3f %5d/%-3d %10.1f\n", fwhm, v.label, nv, thr_sum,
                        ratio_sum / nobj, rel_sum / nobj, found, nobj, t_pass);
            std::fflush(stdout);
        }
    }
    std::printf("\nS/N ratio: detection-map S/N at the true position / ideal stacked S/N (simulator)\n");
    std::printf("S/N / thr: the same S/N relative to the calibrated threshold (higher = easier to detect)\n");
    std::printf("found: injected objects recovered as passing candidates\n");
    return 0;
}

// ---------------------------------------------------------------- 2. curvature
void add_gauss(Image& img, double x0, double y0, double flux, double sigma) {
    const int R = static_cast<int>(std::ceil(4 * sigma));
    const double norm = flux / (2 * 3.14159265358979 * sigma * sigma);
    for (int y = static_cast<int>(y0) - R; y <= static_cast<int>(y0) + R; ++y)
        for (int x = static_cast<int>(x0) - R; x <= static_cast<int>(x0) + R; ++x)
            if (img.inside(x, y)) img.at(x, y) += static_cast<float>(norm * std::exp(-((x - x0) * (x - x0) + (y - y0) * (y - y0)) / (2 * sigma * sigma)));
}

int study_curvature() {
    // geometry: how far a track departs from a straight line over a sequence of T hours
    std::printf("Departure of a track from a straight line over T hours (aT^2/12, arcsec)\n");
    std::printf("diurnal parallax at latitude 51 deg, object at distance D (AU):\n");
    std::printf("%8s", "D \\ T");
    for (double T : {1.0, 2.0, 4.0, 6.0}) std::printf(" %7.0f h", T);
    std::printf("\n");
    for (double D : {0.01, 0.02, 0.05, 0.1, 0.2, 0.5, 1.0, 2.0}) {
        std::printf("%8.2f", D);
        const double P = 8.794 * std::cos(51.0 * 3.14159265 / 180) / D;  // parallax amplitude, arcsec
        const double w = 2 * 3.14159265 / 23.934;                          // rad/h
        for (double T : {1.0, 2.0, 4.0, 6.0}) std::printf(" %9.2f", P * w * w * T * T / 12);
        std::printf("\n");
    }
    std::printf("\nS/N kept by straight-line stacking of a curved track (FWHM 2.4 px, 60 frames over 59 min)\n");
    std::printf("%12s %12s %12s %12s\n", "departure", "departure", "S/N straight", "kept");
    std::printf("%12s %12s %12s %12s\n", "(px)", "(FWHM)", "/ curved", "");
    SimParams s = base_sim(77, 256, 2.4);
    s.nmovers = 0;
    const SimResult base = simulate(s);
    SearchParams p;
    p.vmin_pxh = 5;
    p.vmax_pxh = 30;
    const double sigma = 2.4 / 2.3548;
    std::vector<double> t;
    for (const auto& f : base.frames) t.push_back((f.mjd_mid - base.t_ref_mjd) * 24.0);
    const double T = (t.back() - t.front());
    for (double dev : {0.0, 0.5, 1.0, 1.5, 2.0, 3.0, 5.0}) {
        // a = 12 dev / T^2 along y; motion 15 px/h along x; bright enough to measure well
        const double a = 12.0 * dev / (T * T);
        std::vector<Frame> frames = base.frames;
        const double x0 = 128, y0 = 128, vx = 15.0, vy = 3.0;
        for (std::size_t i = 0; i < frames.size(); ++i)
            add_gauss(frames[i].image, x0 + vx * t[i], y0 + vy * t[i] + 0.5 * a * (t[i] * t[i] - T * T / 12.0), 2500.0,
                      sigma);
        const Prepared prep = prepare(frames, p);
        const auto grid = velocity_grid(prep, p);
        auto be = make_backend("cpu");
        const auto pass = detail::run_pass(*be, prep, grid, prep.t_h, p, {}, 0, 1, "s");
        const float lin = map_peak(pass.max_snr, x0, y0, 3);
        // stack along the true curved track (bilinear), the ideal of a shift-table kernel
        double sum = 0;
        int cnt = 0;
        for (std::size_t i = 0; i < frames.size(); ++i) {
            const double ti = prep.t_h[i];
            const float v = detail::bilinear(prep.filtered[i], x0 + vx * ti, y0 + vy * ti + 0.5 * a * (ti * ti - T * T / 12.0));
            if (std::isfinite(v)) {
                sum += v;
                ++cnt;
            }
        }
        const double curved = sum / std::sqrt(static_cast<double>(cnt));
        std::printf("%12.1f %12.2f %6.1f /%5.1f %11.0f%%\n", dev, dev / 2.4, lin, curved, 100.0 * lin / curved);
        std::fflush(stdout);
    }
    return 0;
}

// ---------------------------------------------------------------- 3. weights
// Relative point-source efficiency of every frame: the matched-filter response of
// moderately bright stars in the background-subtracted frame, divided by the frame
// noise, normalised per star by its median over the frames.
std::vector<double> frame_efficiency(const Prepared& p, int* nstars_used = nullptr) {
    const auto stars = detect_stars(p.static_sky, 150, 10.0);
    const double s = p.fwhm_px * detail::kFwhmToSigma;
    const int R = static_cast<int>(std::ceil(3 * s));
    std::vector<std::vector<double>> resp;  // per star, per frame
    const auto& subs = detail::sub_frames(p);
    for (std::size_t k = 15; k < stars.size(); ++k) {  // skip the brightest (saturation)
        const int cx = static_cast<int>(std::lround(stars[k].x)), cy = static_cast<int>(std::lround(stars[k].y));
        std::vector<double> r(static_cast<std::size_t>(p.nframes));
        bool ok = true;
        for (int i = 0; i < p.nframes && ok; ++i) {
            double acc = 0;
            for (int dy = -R; dy <= R && ok; ++dy)
                for (int dx = -R; dx <= R; ++dx) {
                    const Image& im = subs[static_cast<std::size_t>(i)];
                    if (!im.inside(cx + dx, cy + dy) || !std::isfinite(im.at(cx + dx, cy + dy))) {
                        ok = false;
                        break;
                    }
                    acc += im.at(cx + dx, cy + dy) * std::exp(-(dx * dx + dy * dy) / (2 * s * s));
                }
            r[static_cast<std::size_t>(i)] = acc / p.sigma[static_cast<std::size_t>(i)];
        }
        if (!ok) continue;
        std::vector<double> sorted = r;
        std::nth_element(sorted.begin(), sorted.begin() + sorted.size() / 2, sorted.end());
        const double med = sorted[sorted.size() / 2];
        if (med <= 0) continue;
        for (auto& v : r) v /= med;
        resp.push_back(r);
    }
    if (nstars_used) *nstars_used = static_cast<int>(resp.size());
    std::vector<double> eff(static_cast<std::size_t>(p.nframes), 1.0);
    if (resp.size() < 5) return eff;
    for (int i = 0; i < p.nframes; ++i) {
        std::vector<double> col;
        for (const auto& r : resp) col.push_back(r[static_cast<std::size_t>(i)]);
        std::nth_element(col.begin(), col.begin() + col.size() / 2, col.end());
        eff[static_cast<std::size_t>(i)] = std::max(0.0, col[col.size() / 2]);
    }
    return eff;
}

double weight_gain(const std::vector<double>& s) {
    double sum = 0, sum2 = 0;
    for (double v : s) {
        sum += v;
        sum2 += v * v;
    }
    return std::sqrt(s.size() * sum2) / sum;
}

int study_weights_real(const std::vector<std::string>& dirs, const std::string& wcs_from) {
    std::printf("%-28s %6s %6s %8s %8s %8s %10s\n", "data", "frames", "stars", "eff min", "eff max", "eff CV",
                "S/N gain");
    for (const auto& d : dirs) {
        std::vector<std::string> files;
        for (const auto& e : fs::directory_iterator(d))
            if (e.path().extension() == ".fits" || e.path().extension() == ".fit") files.push_back(e.path().string());
        std::sort(files.begin(), files.end());
        auto frames = load_frames(files);
        SearchParams p;
        const Prepared prep = prepare(std::move(frames), p);
        int ns = 0;
        const auto eff = frame_efficiency(prep, &ns);
        double mean = std::accumulate(eff.begin(), eff.end(), 0.0) / eff.size(), var = 0;
        for (double v : eff) var += (v - mean) * (v - mean);
        const double cv = std::sqrt(var / eff.size()) / mean;
        std::printf("%-28s %6d %6d %8.2f %8.2f %8.3f %9.1f%%\n", fs::path(d).filename().string().c_str(), prep.nframes,
                    ns, *std::min_element(eff.begin(), eff.end()), *std::max_element(eff.begin(), eff.end()), cv,
                    100 * (weight_gain(eff) - 1));
        std::fflush(stdout);
    }
    (void)wcs_from;
    std::printf("\nefficiency: relative point-source S/N of each frame (seeing, transparency, noise)\n");
    std::printf("S/N gain: optimal weighting (w ~ efficiency) against equal weights, sqrt(n sum s^2) / sum s\n");
    return 0;
}

// Weighted stacking emulated with the unchanged kernel: frame i is scaled by
// w_i sqrt(n) / sqrt(sum w^2), which gives sum(w g) / sqrt(sum w^2) wherever all
// frames cover a pixel.
void apply_weights(Prepared& p, const std::vector<double>& w) {
    double s2 = 0;
    for (double v : w) s2 += v * v;
    for (int i = 0; i < p.nframes; ++i) {
        const float k = static_cast<float>(w[static_cast<std::size_t>(i)] * std::sqrt(static_cast<double>(p.nframes) / s2));
        for (auto& v : p.filtered[static_cast<std::size_t>(i)].data) v *= k;
    }
}

int study_weights_sim() {
    std::printf("%-26s %10s %10s %10s %8s\n", "sequence", "weights", "threshold", "S/N / thr", "detected");
    for (double jit : {0.1, 0.25, 0.4}) {
        SimParams s = base_sim(900 + static_cast<std::uint64_t>(jit * 100), 384, 2.5);
        s.fwhm_jitter = jit;
        s.transparency_jitter = jit;
        s.nmovers = 40;
        s.rate_min_arcsec_min = 0.15;
        s.rate_max_arcsec_min = 0.7;
        s.mover_mag_bright = 20.6;
        s.mover_mag_faint = 21.4;
        const SimResult sim = simulate(s);
        SearchParams p;
        p.vmin_pxh = 5;
        p.vmax_pxh = 30;
        const Prepared base = prepare(sim.frames, p);
        const auto eff = frame_efficiency(base);
        for (int weighted = 0; weighted < 2; ++weighted) {
            std::vector<double> w = weighted ? eff : std::vector<double>(eff.size(), 1.0);
            Prepared prep = base;
            apply_weights(prep, w);
            double dv = 0;
            const auto grid = velocity_grid(prep, p, &dv);
            auto be = make_backend("cpu");
            const int nms = std::max(2, static_cast<int>(std::lround(prep.fwhm_px)));
            const auto real = detail::run_pass(*be, prep, grid, prep.t_h, p, {}, 0, 1, "real");
            std::vector<float> pooled;
            const auto orders = detail::null_time_orders(prep.t_h, p.null_runs, p.seed);
            for (const auto& t : orders) {
                Prepared q = detail::with_times(base, t, p);
                apply_weights(q, w);  // weights belong to the frame, not to its time
                const auto null = detail::run_pass(*be, q, grid, t, p, {}, 0, 1, "null");
                const auto npk = detail::find_peaks(null, 3.0f, nms, std::max(p.edge_px, 1), 20000);
                std::vector<char> art(static_cast<std::size_t>(detail::null_tested_count(p, npk.size())), 0);
                for (std::size_t j = 0; j < art.size(); ++j)
                    art[j] = detail::is_artefact_peak(q, npk[j], grid[static_cast<std::size_t>(npk[j].vidx)], p);
                detail::pool_null_peaks(npk, art, pooled);
            }
            SearchResult r;
            detail::set_threshold(prep, p, grid.size(), p.null_runs, pooled, r);
            double rel = 0;
            int det = 0;
            for (const auto& tr : sim.truth) {
                const float v = map_peak(real.max_snr, tr.mover.x_ref, tr.mover.y_ref);
                rel += v / r.threshold_snr;
                det += v >= r.threshold_snr;
            }
            char label[64];
            std::snprintf(label, sizeof label, "seeing+transp. jitter %.0f%%", 100 * jit);
            std::printf("%-26s %10s %10.2f %10.3f %5d/%zu\n", label, weighted ? "by eff." : "equal", r.threshold_snr,
                        rel / sim.truth.size(), det, sim.truth.size());
            std::fflush(stdout);
        }
        std::printf("%-26s predicted gain from the efficiencies: %.1f%%\n", "", 100 * (weight_gain(eff) - 1));
    }
    std::printf("\ndetected: injected objects whose detection-map S/N reaches the calibrated threshold\n");
    return 0;
}

}  // namespace

int main(int argc, char** argv) {
    tune_allocator_for_images();
    const std::string mode = argc > 1 ? argv[1] : "";
    if (mode == "shifts") return study_shifts();
    if (mode == "curvature") return study_curvature();
    if (mode == "weights-sim") return study_weights_sim();
    if (mode == "weights") {
        std::vector<std::string> dirs;
        std::string wcs;
        for (int i = 2; i < argc; ++i) {
            if (std::string(argv[i]) == "--wcs-from" && i + 1 < argc) wcs = argv[++i];
            else dirs.push_back(argv[i]);
        }
        return study_weights_real(dirs, wcs);
    }
    std::fprintf(stderr, "usage: %s shifts | curvature | weights DIR... | weights-sim\n", argv[0]);
    return 1;
}
