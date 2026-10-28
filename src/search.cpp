// SPDX-License-Identifier: GPL-3.0-or-later
// Copyright 2026 Jurgen Kobierczynski
//
// CPU reference implementation of the synthetic-tracking search.
//
// For every trial velocity v the matched-filtered S/N frames g_i are shifted by
// round(v * t_i) and summed; S_v(p) = sum_i g_i(p + v t_i) / sqrt(n(p)) is unit-
// variance noise wherever no object moves with velocity v. The per-pixel maximum
// over all velocities is kept, local maxima above the threshold become
// candidates, and every candidate is refined and tested individually.
//
// The threshold is calibrated empirically: the same search is repeated with the
// frame times randomly permuted, which destroys any coherent motion but keeps
// every other property of the data (noise, artefacts, star residuals). The tail
// of the resulting noise-peak distribution sets the threshold for the requested
// number of false detections per search.
#include <algorithm>
#include <atomic>
#include <chrono>
#include <cmath>
#include <functional>
#include <limits>
#include <numeric>
#include <stdexcept>

#include "common.hpp"
#include "internal.hpp"
#include "psf.hpp"
#include "random.hpp"
#include "fainttrail/backend.hpp"
#include "fainttrail/search.hpp"

namespace fainttrail {

using detail::kFwhmToSigma;
using detail::kPi;

std::vector<std::pair<double, double>> velocity_grid(const Prepared& p, const SearchParams& params,
                                                     double* dv_out) {
    if (p.span_h <= 0.0) throw std::runtime_error("invalid time span");
    if (params.vmax_pxh <= 0.0 || params.vmax_pxh < params.vmin_pxh)
        throw std::runtime_error("invalid velocity range");
    const double dv = 4.0 * std::max(0.05, params.grid_tol_px) / p.span_h;
    if (dv_out) *dv_out = dv;
    const int K = static_cast<int>(std::ceil(params.vmax_pxh / dv));
    std::vector<std::pair<double, double>> grid;
    const double vlo = std::max(0.0, params.vmin_pxh - 0.5 * dv);
    const double vhi = params.vmax_pxh + 0.5 * dv;
    for (int j = -K; j <= K; ++j)
        for (int i = -K; i <= K; ++i) {
            const double vx = i * dv, vy = j * dv;
            const double v = std::hypot(vx, vy);
            if (v >= vlo && v <= vhi) grid.emplace_back(vx, vy);
        }
    if (grid.size() > 5000000) throw std::runtime_error("velocity grid too large; narrow the range");
    return grid;
}

namespace detail {


// One stacking pass through the selected backend; the backend's own progress
// (0..1) is mapped onto [f0, f1] of the caller's progress.
PassOutput run_pass(SearchBackend& backend, const Prepared& p, const std::vector<std::pair<double, double>>& grid,
                    const std::vector<double>& t, const SearchParams& params, const ProgressFn& progress, double f0,
                    double f1, const std::string& stage) {
    StackPass pass;
    pass.frames = &p.filtered;
    pass.t_h = &t;
    pass.velocities = &grid;
    pass.min_count = std::max(3, static_cast<int>(std::ceil(params.min_coverage * p.nframes)));
    pass.threads = params.threads;
    ProgressFn cb;
    if (progress) cb = [&](double f, const std::string&) { return progress(f0 + (f1 - f0) * f, stage); };
    return backend.run(pass, cb);
}

std::vector<Peak> find_peaks(const PassOutput& pass, float floor, int radius, int edge,
                             std::size_t max_peaks) {
    const Image& M = pass.max_snr;
    const int w = M.width, h = M.height;
    std::vector<Peak> peaks;
    for (int y = edge; y < h - edge; ++y) {
        for (int x = edge; x < w - edge; ++x) {
            const float v = M.at(x, y);
            if (v < floor) continue;
            bool is_max = true;
            for (int dy = -radius; dy <= radius && is_max; ++dy) {
                const int yy = y + dy;
                if (yy < 0 || yy >= h) continue;
                for (int dx = -radius; dx <= radius; ++dx) {
                    const int xx = x + dx;
                    if ((!dx && !dy) || xx < 0 || xx >= w) continue;
                    const float u = M.at(xx, yy);
                    // strict for later pixels, non-strict for earlier ones: one peak per plateau
                    if (u > v || (u == v && (dy < 0 || (dy == 0 && dx < 0)))) {
                        is_max = false;
                        break;
                    }
                }
            }
            if (is_max) peaks.push_back({x, y, v, pass.best_v[static_cast<std::size_t>(y) * w + x]});
        }
    }
    std::sort(peaks.begin(), peaks.end(), [](const Peak& a, const Peak& b) { return a.snr > b.snr; });
    if (peaks.size() > max_peaks) peaks.resize(max_peaks);
    return peaks;
}

// Upper-tail probability of a standard normal.
double gauss_q(double z) { return 0.5 * std::erfc(z / std::sqrt(2.0)); }


// S/N along a sub-pixel track. which: 0 all, 1 first half (t<0), 2 second half.
TrackEval track_snr(const Prepared& p, double x, double y, double vx, double vy, int which,
                    std::vector<double>* per_frame) {
    double sum = 0.0;
    int count = 0;
    if (per_frame) per_frame->assign(static_cast<std::size_t>(p.nframes), NAN);
    for (int i = 0; i < p.nframes; ++i) {
        const double t = p.t_h[static_cast<std::size_t>(i)];
        if (which == 1 && t >= 0.0) continue;
        if (which == 2 && t < 0.0) continue;
        const float v = detail::bilinear(p.filtered[static_cast<std::size_t>(i)], x + vx * t, y + vy * t);
        if (!std::isfinite(v)) continue;
        sum += v;
        ++count;
        if (per_frame) (*per_frame)[static_cast<std::size_t>(i)] = v;
    }
    TrackEval e;
    e.count = count;
    e.snr = count >= 3 ? sum / std::sqrt(static_cast<double>(count)) : 0.0;
    return e;
}

// Field edge: in registered frames the blank border moves with the drift; tracks
// that touch it in more than a few frames pick up edge residuals.
bool edge_flag(const Prepared& p, double bx, double by, double bvx, double bvy) {
    const int r = static_cast<int>(std::ceil(2.0 * p.fwhm_px)) + 1;
    int near_blank = 0, total = 0;
    for (int i = 0; i < p.nframes; ++i) {
        const double t = p.t_h[static_cast<std::size_t>(i)];
        const int xi = static_cast<int>(std::lround(bx + bvx * t)), yi = static_cast<int>(std::lround(by + bvy * t));
        // a residual pixel is blank exactly where the background-subtracted frame or
        // the static sky is undefined (see build_stage)
        const Image& sub = detail::sub_frames(p)[static_cast<std::size_t>(i)];
        ++total;
        bool blank = false;
        for (int dy = -r; dy <= r && !blank; ++dy)
            for (int dx = -r; dx <= r; ++dx) {
                if (!sub.inside(xi + dx, yi + dy) || !std::isfinite(sub.at(xi + dx, yi + dy)) ||
                    !std::isfinite(p.static_sky.at(xi + dx, yi + dy))) {
                    blank = true;
                    break;
                }
            }
        near_blank += blank;
    }
    return total > 0 && near_blank > std::max(1, total / 10);
}

// Detector-fixed: with registration offsets known, compare the linear sky track
// with the hypothesis "fixed on the sensor" (warm pixel or defect carried along by
// the registration, which follows the mount's drift and periodic error). If a
// fixed detector position explains the signal about as well, reject.
bool detector_flag(const Prepared& p, double bx, double by, double bvx, double bvy, double track_snr_all) {
    if (static_cast<int>(p.reg_dx.size()) != p.nframes) return false;
    double mx = 0, my = 0;
    for (int i = 0; i < p.nframes; ++i) {
        const double t = p.t_h[static_cast<std::size_t>(i)];
        mx += bx + bvx * t + p.reg_dx[static_cast<std::size_t>(i)];
        my += by + bvy * t + p.reg_dy[static_cast<std::size_t>(i)];
    }
    mx /= p.nframes;
    my /= p.nframes;
    auto fixed_snr = [&](double dxp, double dyp) {
        double sum = 0.0;
        int count = 0;
        for (int i = 0; i < p.nframes; ++i) {
            const float v = detail::bilinear(p.filtered[static_cast<std::size_t>(i)], dxp - p.reg_dx[static_cast<std::size_t>(i)],
                                             dyp - p.reg_dy[static_cast<std::size_t>(i)]);
            if (!std::isfinite(v)) continue;
            sum += v;
            ++count;
        }
        return count >= 3 ? sum / std::sqrt(static_cast<double>(count)) : 0.0;
    };
    double best_fixed = 0.0;
    const double span = std::max(2.0, 2.0 * p.fwhm_px);
    for (double dy = -span; dy <= span; dy += 0.5)
        for (double dx = -span; dx <= span; dx += 0.5) best_fixed = std::max(best_fixed, fixed_snr(mx + dx, my + dy));
    return best_fixed >= 0.8 * track_snr_all;
}

// The artefact tests (single, edge, detector) at a peak's grid position and
// motion, without refinement: used to keep artefacts out of the null tail.
bool is_artefact_peak(const Prepared& p, const Peak& pk, const std::pair<double, double>& v, const SearchParams& params) {
    std::vector<double> per_frame;
    const TrackEval all = track_snr(p, pk.x, pk.y, v.first, v.second, 0, &per_frame);
    double pos_sum = 0.0, maxc = 0.0;
    for (double x : per_frame) {
        if (!std::isfinite(x)) continue;
        pos_sum += std::max(0.0, x);
        maxc = std::max(maxc, x);
    }
    const double frac = pos_sum > 0.0 ? maxc / pos_sum : 1.0;
    if (p.nframes >= 8 && frac > params.max_frame_fraction) return true;
    if (edge_flag(p, pk.x, pk.y, v.first, v.second)) return true;
    return detector_flag(p, pk.x, pk.y, v.first, v.second, all.snr);
}

Candidate refine(const Prepared& p, const Peak& pk, const std::pair<double, double>& v0, double dv,
                 const SearchParams& params) {
    Candidate c;
    c.snr_coarse = pk.snr;
    double bx = pk.x, by = pk.y, bvx = v0.first, bvy = v0.second;
    double best = track_snr(p, bx, by, bvx, bvy, 0).snr;

    // two-level local search: (velocity step, velocity half-range, position step, position half-range)
    const double levels[2][4] = {{dv / 4, dv, 0.5, 1.0}, {dv / 8, dv / 4, 0.25, 0.5}};
    for (const auto& L : levels) {
        const double cx = bx, cy = by, cvx = bvx, cvy = bvy;
        const int nvs = static_cast<int>(std::lround(L[1] / L[0]));
        const int nps = static_cast<int>(std::lround(L[3] / L[2]));
        for (int a = -nvs; a <= nvs; ++a)
            for (int b = -nvs; b <= nvs; ++b)
                for (int i = -nps; i <= nps; ++i)
                    for (int j = -nps; j <= nps; ++j) {
                        const double vx = cvx + a * L[0], vy = cvy + b * L[0];
                        const double x = cx + i * L[2], y = cy + j * L[2];
                        const double s = track_snr(p, x, y, vx, vy, 0).snr;
                        if (s > best) {
                            best = s;
                            bx = x; by = y; bvx = vx; bvy = vy;
                        }
                    }
    }
    // sub-pixel position from a parabola through three samples per axis
    {
        const double h = 0.25;
        const double sx0 = track_snr(p, bx - h, by, bvx, bvy, 0).snr;
        const double sx2 = track_snr(p, bx + h, by, bvx, bvy, 0).snr;
        const double sy0 = track_snr(p, bx, by - h, bvx, bvy, 0).snr;
        const double sy2 = track_snr(p, bx, by + h, bvx, bvy, 0).snr;
        const double denx = sx0 - 2 * best + sx2, deny = sy0 - 2 * best + sy2;
        if (denx < 0) bx += std::clamp(0.5 * h * (sx0 - sx2) / denx, -h, h);
        if (deny < 0) by += std::clamp(0.5 * h * (sy0 - sy2) / deny, -h, h);
    }

    // sub-step motion from a parabola as well
    {
        const double h = dv / 8;
        const double s0 = track_snr(p, bx, by, bvx - h, bvy, 0).snr, s2 = track_snr(p, bx, by, bvx + h, bvy, 0).snr;
        const double t0 = track_snr(p, bx, by, bvx, bvy - h, 0).snr, t2 = track_snr(p, bx, by, bvx, bvy + h, 0).snr;
        const double sc = track_snr(p, bx, by, bvx, bvy, 0).snr;
        const double denx = s0 - 2 * sc + s2, deny = t0 - 2 * sc + t2;
        if (denx < 0) bvx += std::clamp(0.5 * h * (s0 - s2) / denx, -h, h);
        if (deny < 0) bvy += std::clamp(0.5 * h * (t0 - t2) / deny, -h, h);
    }

    std::vector<double> per_frame;
    const TrackEval all = track_snr(p, bx, by, bvx, bvy, 0, &per_frame);
    c.x = bx;
    c.y = by;
    c.vx_pxh = bvx;
    c.vy_pxh = bvy;
    c.snr = all.snr;
    {
        // centroid precision of a matched-filter fit: sqrt(2) * sigma_psf / (S/N) per axis
        const double s_psf = p.fwhm_px * kFwhmToSigma;
        // The 1.25 factor is empirical: on simulated sequences the ideal matched-filter
        // precision underestimates the actual scatter by about 25 % (integer-shift
        // grid, static-sky noise, PSF variations). See docs/VALIDATION.md.
        constexpr double kEmpirical = 1.25;
        c.pos_err_px = c.snr > 0.0 ? kEmpirical * std::sqrt(2.0) * s_psf / c.snr : 0.0;
        double t2 = 0.0;
        for (double t : p.t_h) t2 += t * t;
        t2 /= std::max(1, p.nframes);
        c.vel_err_pxh = t2 > 0.0 ? c.pos_err_px / std::sqrt(t2) : 0.0;
    }
    c.snr_half1 = track_snr(p, bx, by, bvx, bvy, 1).snr;
    c.snr_half2 = track_snr(p, bx, by, bvx, bvy, 2).snr;

    double pos_sum = 0.0, maxc = 0.0;
    int masked = 0, inside = 0;
    for (int i = 0; i < p.nframes; ++i) {
        const double v = per_frame[static_cast<std::size_t>(i)];
        if (!std::isfinite(v)) continue;
        ++inside;
        pos_sum += std::max(0.0, v);
        maxc = std::max(maxc, v);
        const double t = p.t_h[static_cast<std::size_t>(i)];
        const int xi = static_cast<int>(std::lround(bx + bvx * t)), yi = static_cast<int>(std::lround(by + bvy * t));
        if (p.mask.inside(xi, yi) && p.mask.at(xi, yi) == 0.0f) ++masked;
    }
    c.max_frame_fraction = pos_sum > 0.0 ? maxc / pos_sum : 1.0;

    // PSF photometry on the stacked residuals
    const double s = p.fwhm_px * kFwhmToSigma;
    const Image cut = stack_cutout(p, bx, by, bvx, bvy, static_cast<int>(std::ceil(3.0 * s)), false, 0);
    double num = 0.0, den = 0.0;
    const int R = cut.width / 2;
    for (int y = 0; y < cut.height; ++y)
        for (int x = 0; x < cut.width; ++x) {
            const float v = cut.at(x, y);
            if (!std::isfinite(v)) continue;
            const double k = std::exp(-0.5 * ((x - R) * (x - R) + (y - R) * (y - R)) / (s * s));
            num += v * k;
            den += k * k;
        }
    c.flux = den > 0.0 ? num / den * 2.0 * kPi * s * s : 0.0;
    c.mag_inst = c.flux > 0.0 ? 25.0 - 2.5 * std::log10(c.flux) : 99.0;

    // flags
    std::vector<std::string> flags;
    const double expect_half = c.snr_coarse / std::sqrt(2.0);
    const double half_min = std::max(2.0, expect_half - params.half_tolerance);
    if (c.snr_half1 < half_min || c.snr_half2 < half_min) flags.push_back("half");
    if (p.nframes >= 8 && c.max_frame_fraction > params.max_frame_fraction) flags.push_back("single");
    if (inside > 0 && masked > inside / 2) flags.push_back("masked");
    if (edge_flag(p, bx, by, bvx, bvy)) flags.push_back("edge");
    if (detector_flag(p, bx, by, bvx, bvy, all.snr)) flags.push_back("detector");
    if (std::hypot(bvx, bvy) > params.vmax_pxh - 0.5 * dv) flags.push_back("vedge");
    if (std::hypot(bvx, bvy) < 2.0 * p.fwhm_px / p.span_h) flags.push_back("slow");
    for (std::size_t i = 0; i < flags.size(); ++i) c.flags += (i ? "," : "") + flags[i];
    return c;
}

bool is_rejecting(const std::string& flag) {
    return flag == "half" || flag == "single" || flag == "ghost" || flag == "edge" || flag == "detector";
}

void finalise_flags(Candidate& c) {
    c.ok = true;
    std::size_t start = 0;
    while (start <= c.flags.size() && !c.flags.empty()) {
        const std::size_t end = c.flags.find(',', start);
        const std::string f = c.flags.substr(start, end == std::string::npos ? std::string::npos : end - start);
        if (is_rejecting(f)) c.ok = false;
        if (end == std::string::npos) break;
        start = end + 1;
    }
}

void add_flag(Candidate& c, const std::string& f) {
    if (c.flags.find(f) != std::string::npos) return;
    c.flags += c.flags.empty() ? f : "," + f;
}

void sky_coordinates(const Prepared& p, Candidate& c) {
    if (!p.wcs.valid) return;
    double ra, dec, ra2, dec2;
    p.wcs.pix2sky(c.x, c.y, ra, dec);
    const double dt = 0.1;  // hours
    p.wcs.pix2sky(c.x + c.vx_pxh * dt, c.y + c.vy_pxh * dt, ra2, dec2);
    c.ra_deg = ra;
    c.dec_deg = dec;
    const double d2r = kPi / 180.0;
    double dra = ra2 - ra;
    if (dra > 180.0) dra -= 360.0;
    if (dra < -180.0) dra += 360.0;
    const double east = dra * std::cos(dec * d2r) * 3600.0;  // arcsec
    const double north = (dec2 - dec) * 3600.0;
    c.rate_arcsec_min = std::hypot(east, north) / (dt * 60.0);
    double pa = std::atan2(east, north) / d2r;
    if (pa < 0) pa += 360.0;
    c.pa_deg = pa;
}

// Fit log N(>s) = a - b s to the pooled tail of the null peaks.
// Tail of the peaks of a smooth Gaussian random field: by the Euler-characteristic
// heuristic the expected number of local maxima above u falls like
// u^(D-1) exp(-u^2 / 2) for a unit-variance field in D dimensions. The null
// detection map is the maximum over a family of correlated Gaussian fields
// (position and motion), with a variance that is not exactly 1 after the robust
// normalisation, so the power k and the variance are fitted:
//   log N(>u) = a + k log u - inv_s2 u^2 / 2
// by least squares on the top `max_per_run` peaks per run (all ranks from 1).
bool fit_null_tail_gauss(std::vector<float> pooled, double runs, double& a, double& k, double& inv_s2,
                         double max_per_run) {
    std::sort(pooled.begin(), pooled.end(), std::greater<float>());
    const std::size_t m = std::min(pooled.size(), static_cast<std::size_t>(max_per_run * runs));
    if (m < 8) return false;
    // normal equations for [1, log u, -u^2/2]
    double S[3][3] = {{0}}, R[3] = {0};
    for (std::size_t j = 0; j < m; ++j) {
        const double u = pooled[j];
        if (!(u > 0.0)) return false;
        const double f[3] = {1.0, std::log(u), -0.5 * u * u};
        const double y = std::log(static_cast<double>(j + 1) / runs);
        for (int r = 0; r < 3; ++r) {
            R[r] += f[r] * y;
            for (int c = 0; c < 3; ++c) S[r][c] += f[r] * f[c];
        }
    }
    // solve the 3x3 system (Gaussian elimination with partial pivoting)
    double A[3][4];
    for (int r = 0; r < 3; ++r) {
        for (int c = 0; c < 3; ++c) A[r][c] = S[r][c];
        A[r][3] = R[r];
    }
    for (int c = 0; c < 3; ++c) {
        int piv = c;
        for (int r = c + 1; r < 3; ++r)
            if (std::fabs(A[r][c]) > std::fabs(A[piv][c])) piv = r;
        if (std::fabs(A[piv][c]) < 1e-300) return false;
        for (int q = 0; q < 4; ++q) std::swap(A[c][q], A[piv][q]);
        for (int r = 0; r < 3; ++r) {
            if (r == c) continue;
            const double f = A[r][c] / A[c][c];
            for (int q = c; q < 4; ++q) A[r][q] -= f * A[c][q];
        }
    }
    a = A[0][3] / A[0][0];
    k = A[1][3] / A[1][1];
    inv_s2 = A[2][3] / A[2][2];
    return std::isfinite(a) && std::isfinite(k) && inv_s2 > 0.0 && std::isfinite(inv_s2);
}

double gauss_tail_log_rate(double a, double k, double inv_s2, double u) {
    return a + k * std::log(u) - 0.5 * inv_s2 * u * u;
}

// Level where the fitted tail reaches `rate` peaks per run, on its falling branch.
double gauss_tail_level(double a, double k, double inv_s2, double rate) {
    const double lt = std::log(rate);
    double lo = std::max(1.0, k > 0.0 ? std::sqrt(k / inv_s2) : 1.0), hi = 100.0;
    if (gauss_tail_log_rate(a, k, inv_s2, lo) <= lt) return lo;
    for (int it = 0; it < 200; ++it) {
        const double mid = 0.5 * (lo + hi);
        if (gauss_tail_log_rate(a, k, inv_s2, mid) > lt) lo = mid;
        else hi = mid;
    }
    return hi;
}

bool fit_null_tail(std::vector<float> pooled, double runs, double& a, double& b) {
    std::sort(pooled.begin(), pooled.end(), std::greater<float>());
    std::vector<double> xs, ys;
    for (std::size_t j = 0; j < pooled.size(); ++j) {
        const double n_per_run = static_cast<double>(j + 1) / runs;
        if (n_per_run < 0.67) continue;
        if (n_per_run > 30.0) break;
        xs.push_back(pooled[j]);
        ys.push_back(std::log(n_per_run));
    }
    if (xs.size() < 8) return false;
    const double mx = std::accumulate(xs.begin(), xs.end(), 0.0) / xs.size();
    const double my = std::accumulate(ys.begin(), ys.end(), 0.0) / ys.size();
    double sxy = 0, sxx = 0;
    for (std::size_t i = 0; i < xs.size(); ++i) {
        sxy += (xs[i] - mx) * (ys[i] - my);
        sxx += (xs[i] - mx) * (xs[i] - mx);
    }
    if (sxx <= 0.0) return false;
    const double slope = sxy / sxx;
    if (slope >= 0.0) return false;
    b = -slope;
    a = my + b * mx;
    return true;
}

}  // namespace detail

namespace detail {

std::vector<std::vector<double>> null_time_orders(const std::vector<double>& t_h, int runs, std::uint64_t seed) {
    detail::Rng rng(seed);
    std::vector<std::vector<double>> out;
    for (int k = 0; k < runs; ++k) {
        std::vector<double> t = t_h;
        for (std::size_t i = t.size() - 1; i > 0; --i) {
            const std::size_t j = static_cast<std::size_t>(rng.uniform() * (i + 1));
            std::swap(t[i], t[j]);
        }
        out.push_back(std::move(t));
    }
    return out;
}

void set_threshold(const Prepared& p, const SearchParams& params, std::size_t n_velocities, int runs,
                   const std::vector<float>& pooled, SearchResult& res) {
    if (params.threshold_snr > 0.0) {
        res.threshold_snr = params.threshold_snr;
        res.threshold_method = "fixed";
    } else if (runs > 0) {
        res.null_peaks.assign(pooled.begin(), pooled.end());
        const double target = std::max(1e-6, params.target_false_per_run);
        double a = 0, b = 0, k = 0, inv_s2 = 0;
        if (params.null_model == "gauss" && fit_null_tail_gauss(pooled, runs, a, k, inv_s2, 3000.0)) {
            res.null_a = a;
            res.null_k = k;
            res.null_inv_s2 = inv_s2;
            res.threshold_snr = gauss_tail_level(a, k, inv_s2, target);
            res.threshold_method = "null-fit gauss";
        } else if (fit_null_tail(pooled, runs, a, b)) {
            res.null_a = a;
            res.null_b = b;
            res.threshold_snr = (a - std::log(target)) / b;
            res.threshold_method = "null-fit linear";
        } else {
            const float mx = pooled.empty() ? 0.0f : *std::max_element(pooled.begin(), pooled.end());
            res.threshold_snr = mx + 0.5;
            res.threshold_method = "null-quantile";
        }
        res.threshold_snr = std::max(res.threshold_snr, params.min_snr_floor);
    } else {
        // Gaussian approximation: independent trials ~ resolution elements x velocities.
        const double resel = 4.0 * kPi * std::pow(p.fwhm_px * kFwhmToSigma, 2);
        const double trials = std::max(1.0, p.width * p.height / resel * n_velocities * 0.5);
        double lo = 3.0, hi = 12.0;
        for (int it = 0; it < 60; ++it) {
            const double mid = 0.5 * (lo + hi);
            if (trials * gauss_q(mid) > params.target_false_per_run) lo = mid; else hi = mid;
        }
        res.threshold_snr = std::max(hi, params.min_snr_floor);
        res.threshold_method = "gaussian";
        res.null_b = 0.0;
    }
}

std::vector<Peak> select_peaks(const std::vector<Peak>& peaks, double threshold, const SearchParams& params) {
    std::vector<Peak> keep;
    for (const auto& pk : peaks)
        if (pk.snr >= threshold && pk.vidx >= 0) keep.push_back(pk);
    if (static_cast<int>(keep.size()) > params.max_candidates) keep.resize(static_cast<std::size_t>(params.max_candidates));
    return keep;
}

// Duplicates and ghosts: a weaker candidate whose track coincides with a stronger
// one is either the same object (most frames) or a sidelobe of it (a large share
// of its signal comes from frames where the tracks overlap). Then final flags,
// false-alarm rates, sky coordinates and ids. `meta` needs only the frame times,
// FWHM and WCS; `contrib` holds each candidate's per-frame track values.
void finish_candidates(const Prepared& meta, std::vector<RefinedCandidate> rc, SearchResult& res) {
    res.has_sky = meta.wcs.valid;
    std::sort(rc.begin(), rc.end(), [](const RefinedCandidate& a, const RefinedCandidate& b) { return a.c.snr > b.c.snr; });
    const Prepared& p = meta;
    const double r_overlap = 1.5 * p.fwhm_px;
    std::vector<bool> drop(rc.size(), false);
    for (std::size_t j = 0; j < rc.size(); ++j) {
        double total = 0.0;
        for (double v : rc[j].contrib)
            if (std::isfinite(v) && v > 0.0) total += v;
        std::vector<bool> shared(static_cast<std::size_t>(p.nframes), false);
        int overlap_max = 0;
        for (std::size_t i = 0; i < j && !drop[j]; ++i) {
            if (drop[i]) continue;
            int overlap = 0;
            for (int f = 0; f < p.nframes; ++f) {
                const double t = p.t_h[static_cast<std::size_t>(f)];
                const double dx = (rc[i].c.x + rc[i].c.vx_pxh * t) - (rc[j].c.x + rc[j].c.vx_pxh * t);
                const double dy = (rc[i].c.y + rc[i].c.vy_pxh * t) - (rc[j].c.y + rc[j].c.vy_pxh * t);
                if (dx * dx + dy * dy < r_overlap * r_overlap) {
                    ++overlap;
                    shared[static_cast<std::size_t>(f)] = true;
                }
            }
            overlap_max = std::max(overlap_max, overlap);
            if (overlap >= p.nframes / 2) drop[j] = true;
        }
        if (drop[j] || overlap_max == 0) continue;
        double from_shared = 0.0;
        for (int f = 0; f < p.nframes; ++f) {
            const double v = rc[j].contrib[static_cast<std::size_t>(f)];
            if (shared[static_cast<std::size_t>(f)] && std::isfinite(v) && v > 0.0) from_shared += v;
        }
        if (total > 0.0 && from_shared / total > 0.3) add_flag(rc[j].c, "ghost");
    }
    for (std::size_t j = 0; j < rc.size(); ++j) {
        if (drop[j]) continue;
        Candidate c = rc[j].c;
        finalise_flags(c);
        c.false_alarms_per_run = expected_null_peaks(res, c.snr_coarse);
        sky_coordinates(p, c);
        res.candidates.push_back(c);
    }
    std::stable_sort(res.candidates.begin(), res.candidates.end(), [](const Candidate& a, const Candidate& b) {
        if (a.ok != b.ok) return a.ok;
        return a.snr > b.snr;
    });
    for (std::size_t i = 0; i < res.candidates.size(); ++i) res.candidates[i].id = static_cast<int>(i + 1);
}

bool is_artefact(const Candidate& c) {
    std::size_t s = 0;
    while (s <= c.flags.size() && !c.flags.empty()) {
        const std::size_t e = c.flags.find(',', s);
        const std::string f = c.flags.substr(s, e == std::string::npos ? std::string::npos : e - s);
        if (f == "edge" || f == "detector" || f == "single") return true;
        if (e == std::string::npos) break;
        s = e + 1;
    }
    return false;
}

void pool_null_peaks(const std::vector<Peak>& peaks, const std::vector<char>& artefact, std::vector<float>& pooled) {
    for (std::size_t k = 0; k < peaks.size(); ++k)
        if (k >= artefact.size() || !artefact[k]) pooled.push_back(peaks[k].snr);
}

int null_tested_count(const SearchParams& params, std::size_t npeaks) {
    if (params.null_model != "gauss" || params.null_tested <= 0) return 0;
    return static_cast<int>(std::min<std::size_t>(npeaks, static_cast<std::size_t>(params.null_tested)));
}

RefinedCandidate refine_with_contrib(const Prepared& p, const Peak& pk, const std::pair<double, double>& v0, double dv,
                                     const SearchParams& params) {
    RefinedCandidate r;
    r.c = refine(p, pk, v0, dv, params);
    track_snr(p, r.c.x, r.c.y, r.c.vx_pxh, r.c.vy_pxh, 0, &r.contrib);
    return r;
}

}  // namespace detail

double expected_null_peaks(const SearchResult& res, double snr) {
    if (res.null_inv_s2 > 0.0 && snr > 0.0) return std::exp(detail::gauss_tail_log_rate(res.null_a, res.null_k, res.null_inv_s2, snr));
    if (res.null_b > 0.0) return std::exp(res.null_a - res.null_b * snr);
    return 0.0;
}

SearchResult search(const Prepared& p, const SearchParams& params, const ProgressFn& progress) {
    using namespace detail;
    const auto t0 = std::chrono::steady_clock::now();
    SearchResult res;
    res.fwhm_px = p.fwhm_px;
    res.t_ref_mjd = p.t_ref_mjd;
    res.slowest_reliable_pxh = 2.0 * p.fwhm_px / p.span_h;
    double dv = 0.0;
    const auto grid = velocity_grid(p, params, &dv);
    res.n_velocities = grid.size();
    const std::unique_ptr<SearchBackend> backend = make_backend(params.backend);
    res.backend = backend->name();
    res.backend_description = backend->description();
    res.dv_pxh = dv;
    const int nt = detail::thread_count(params.threads);
    detail::Stopwatch sw(&res.timings);
    sw.lap("setup");

    const bool auto_thr = params.threshold_snr <= 0.0;
    const int runs = auto_thr ? std::max(0, params.null_runs) : 0;
    const double share = 0.9 / (1 + runs);
    const int nms = std::max(2, static_cast<int>(std::lround(p.fwhm_px)));
    const int edge = std::max(params.edge_px, 1);

    // the real search
    PassOutput real = run_pass(*backend, p, grid, p.t_h, params, progress, 0.0, share, "search");
    sw.lap("stack (real)");
    for (const auto& st : real.timings) res.timings.push_back({"  " + st.name, st.seconds});
    const float floor = static_cast<float>(auto_thr ? params.min_snr_floor : params.threshold_snr);
    std::vector<Peak> peaks = find_peaks(real, floor, nms, edge, 100000);
    res.max_snr = std::move(real.max_snr);
    sw.lap("peaks (real)");

    // null runs: the whole time-dependent pipeline is rerun with permuted times,
    // so the null sees exactly the same noise structure as the real search
    std::vector<float> pooled;
    const auto orders = null_time_orders(p.t_h, runs, params.seed);
    for (int k = 0; k < runs; ++k) {
        const std::vector<double>& t = orders[static_cast<std::size_t>(k)];
        const std::string stage = "null " + std::to_string(k + 1) + "/" + std::to_string(runs);
        if (progress && !progress(share * (k + 1), stage)) throw Cancelled();
        detail::release_free_memory();
        const Prepared q = detail::with_times(p, t, params);
        sw.lap("null rebuild");
        PassOutput null = run_pass(*backend, q, grid, t, params, progress, share * (k + 1), share * (k + 2), stage);
        sw.lap("stack (null)");
        for (const auto& st : null.timings) res.timings.push_back({"  " + st.name, st.seconds});
        const std::vector<Peak> npk = find_peaks(null, 3.0f, nms, edge, 20000);
        // artefact test of the strongest null peaks (the real search rejects these)
        std::vector<char> art(static_cast<std::size_t>(null_tested_count(params, npk.size())), 0);
#pragma omp parallel for num_threads(nt) schedule(dynamic)
        for (int j = 0; j < static_cast<int>(art.size()); ++j) {
            const Peak& pk = npk[static_cast<std::size_t>(j)];
            art[static_cast<std::size_t>(j)] = is_artefact_peak(q, pk, grid[static_cast<std::size_t>(pk.vidx)], params);
        }
        pool_null_peaks(npk, art, pooled);
        sw.lap("peaks (null)");
        // the null run's frames are released when q goes out of scope
    }
    detail::release_free_memory();
    set_threshold(p, params, grid.size(), runs, pooled, res);
    sw.lap("threshold fit");

    // refine in parallel
    const std::vector<Peak> keep = select_peaks(peaks, res.threshold_snr, params);
    std::vector<RefinedCandidate> rc(keep.size());
    std::atomic<int> done{0};
    std::atomic<bool> cancel{false};
#pragma omp parallel for num_threads(nt) schedule(dynamic)
    for (int k = 0; k < static_cast<int>(keep.size()); ++k) {
        if (cancel) continue;
        const Peak& pk = keep[static_cast<std::size_t>(k)];
        rc[static_cast<std::size_t>(k)] = refine_with_contrib(p, pk, grid[static_cast<std::size_t>(pk.vidx)], dv, params);
        const int d = ++done;
        if (detail::thread_id() == 0 && progress && !progress(0.9 + 0.1 * d / keep.size(), "refine"))
            cancel = true;
    }
    if (cancel) throw Cancelled();
    sw.lap("refine");

    finish_candidates(p, std::move(rc), res);
    sw.lap("ghosts + flags");
    if (progress) progress(1.0, "done");
    res.elapsed_s = std::chrono::duration<double>(std::chrono::steady_clock::now() - t0).count();
    return res;
}

SearchResult search(const std::vector<Frame>& frames, const SearchParams& params, const ProgressFn& progress) {
    ProgressFn prep_progress, search_progress;
    if (progress) {
        prep_progress = [&](double f, const std::string& s) { return progress(0.1 * f, s); };
        search_progress = [&](double f, const std::string& s) { return progress(0.1 + 0.9 * f, s); };
    }
    const Prepared p = prepare(frames, params, prep_progress);
    return search(p, params, search_progress);
}

static void require_residuals(const Prepared& p) {
    if (static_cast<int>(p.residual.size()) != p.nframes)
        throw std::runtime_error("residual images were not kept; prepare with SearchParams::keep_residuals = true");
}

Image stack_cutout(const Prepared& p, double x, double y, double vx, double vy, int half_size,
                   bool snr, int which) {
    const int size = 2 * half_size + 1;
    Image out(size, size, 0.0f);
    std::vector<int> count(static_cast<std::size_t>(size) * size, 0);
    for (int i = 0; i < p.nframes; ++i) {
        const double t = p.t_h[static_cast<std::size_t>(i)];
        if (which == 1 && t >= 0.0) continue;
        if (which == 2 && t < 0.0) continue;
        const double cx = x + vx * t, cy = y + vy * t;
        // residuals: the small area the bilinear samples touch, copied or recomputed
        const int px0 = static_cast<int>(std::floor(cx - half_size)) - 1, py0 = static_cast<int>(std::floor(cy - half_size)) - 1;
        Image patch;
        if (!snr) patch = detail::residual_patch(p, i, px0, py0, size + 3, size + 3);
        for (int v = 0; v < size; ++v)
            for (int u = 0; u < size; ++u) {
                const float s = snr ? detail::bilinear(p.filtered[static_cast<std::size_t>(i)], cx + u - half_size, cy + v - half_size)
                                    : detail::bilinear_patch(patch, px0, py0, p.width, p.height, cx + u - half_size,
                                                             cy + v - half_size);
                if (!std::isfinite(s)) continue;
                out.data[static_cast<std::size_t>(v) * size + u] += s;
                ++count[static_cast<std::size_t>(v) * size + u];
            }
    }
    for (std::size_t k = 0; k < out.data.size(); ++k) {
        const int c = count[k];
        if (c == 0) out.data[k] = NAN;
        else out.data[k] = snr ? out.data[k] / std::sqrt(static_cast<float>(c)) : out.data[k] / c;
    }
    return out;
}

Image stack_full(const Prepared& p, double vx, double vy) {
    require_residuals(p);
    Image out(p.width, p.height, 0.0f);
    std::vector<int> count(out.data.size(), 0);
    for (int i = 0; i < p.nframes; ++i) {
        const double t = p.t_h[static_cast<std::size_t>(i)];
        const Image& src = p.residual[static_cast<std::size_t>(i)];
        for (int y = 0; y < p.height; ++y)
            for (int x = 0; x < p.width; ++x) {
                const float s = detail::bilinear(src, x + vx * t, y + vy * t);
                if (!std::isfinite(s)) continue;
                out.data[static_cast<std::size_t>(y) * p.width + x] += s;
                ++count[static_cast<std::size_t>(y) * p.width + x];
            }
    }
    for (std::size_t k = 0; k < out.data.size(); ++k) out.data[k] = count[k] ? out.data[k] / count[k] : NAN;
    return out;
}

std::vector<TrackPoint> track_points(const Prepared& p, const Candidate& c) {
    std::vector<TrackPoint> pts;
    for (int i = 0; i < p.nframes; ++i) {
        TrackPoint tp;
        tp.frame = i;
        const double t = p.t_h[static_cast<std::size_t>(i)];
        tp.mjd = p.t_ref_mjd + t / 24.0;
        tp.x = c.x + c.vx_pxh * t;
        tp.y = c.y + c.vy_pxh * t;
        if (p.wcs.valid) p.wcs.pix2sky(tp.x, tp.y, tp.ra_deg, tp.dec_deg);
        pts.push_back(tp);
    }
    return pts;
}

}  // namespace fainttrail
