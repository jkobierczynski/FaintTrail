// SPDX-License-Identifier: GPL-3.0-or-later
// Copyright 2026 Jurgen Kobierczynski
#include "fainttrail/validate.hpp"

#include <algorithm>
#include <cmath>
#include <limits>

#include "psf.hpp"
#include "random.hpp"

namespace fainttrail {

using detail::kFwhmToSigma;
using detail::kPi;

void inject(std::vector<Frame>& frames, const std::vector<Injection>& objects, double t_ref_mjd,
            double fwhm_px) {
    const double s = fwhm_px * kFwhmToSigma;
    for (auto& f : frames) {
        const double t_mid = (f.mjd_mid - t_ref_mjd) * 24.0;
        const double half = 0.5 * f.exptime_s / 3600.0;
        for (const auto& o : objects) {
            detail::add_trailed(f.image, o.x_ref + o.vx_pxh * (t_mid - half), o.y_ref + o.vy_pxh * (t_mid - half),
                                o.x_ref + o.vx_pxh * (t_mid + half), o.y_ref + o.vy_pxh * (t_mid + half),
                                o.flux, s);
        }
    }
}

double flux_for_snr(double stack_snr, double sigma, double fwhm_px, int nframes) {
    const double s = fwhm_px * kFwhmToSigma;
    return stack_snr / std::sqrt(static_cast<double>(std::max(1, nframes))) * sigma * std::sqrt(4.0 * kPi) * s;
}

std::vector<MatchResult> match(const std::vector<Candidate>& cands, const std::vector<SimMover>& truth,
                               double max_dpos, double max_dv, bool ok_only) {
    std::vector<MatchResult> out;
    std::vector<bool> used(cands.size(), false);
    for (std::size_t t = 0; t < truth.size(); ++t) {
        double best = std::numeric_limits<double>::infinity();
        int bi = -1;
        double bdv = 0;
        for (std::size_t c = 0; c < cands.size(); ++c) {
            if (used[c] || (ok_only && !cands[c].ok)) continue;
            const double dp = std::hypot(cands[c].x - truth[t].x_ref, cands[c].y - truth[t].y_ref);
            const double dv = std::hypot(cands[c].vx_pxh - truth[t].vx_pxh, cands[c].vy_pxh - truth[t].vy_pxh);
            if (dp <= max_dpos && dv <= max_dv && dp < best) {
                best = dp;
                bi = static_cast<int>(c);
                bdv = dv;
            }
        }
        if (bi >= 0) {
            used[static_cast<std::size_t>(bi)] = true;
            out.push_back({static_cast<int>(t), bi, best, bdv});
        }
    }
    return out;
}

CompletenessReport completeness(const std::vector<Frame>& frames, const SearchParams& params,
                                const CompletenessParams& cp, const ProgressFn& progress) {
    CompletenessReport rep;
    const int levels = static_cast<int>(cp.snr_levels.size());
    const double steps = 1.0 + levels;
    auto sub = [&](int k) -> ProgressFn {
        if (!progress) return {};
        return [&, k](double f, const std::string& s) { return progress((k + f) / steps, s); };
    };

    // baseline: calibrate the threshold on the untouched data
    const Prepared base_prep = prepare(frames, params);
    const SearchResult base = search(base_prep, params, sub(0));
    rep.threshold_snr = base.threshold_snr;
    rep.threshold_method = base.threshold_method;
    for (const auto& c : base.candidates)
        if (c.ok) ++rep.baseline_ok;

    std::vector<double> sg(base_prep.sigma);
    std::nth_element(sg.begin(), sg.begin() + static_cast<std::ptrdiff_t>(sg.size() / 2), sg.end());
    const double sigma = sg[sg.size() / 2];
    const double fwhm = base_prep.fwhm_px;
    const int n = base_prep.nframes;
    const int w = base_prep.width, h = base_prep.height;
    const double margin = std::max(12.0, 3.0 * fwhm);

    SearchParams fixed = params;
    fixed.threshold_snr = rep.threshold_snr;
    fixed.fwhm_px = fwhm;
    fixed.null_runs = 0;

    detail::Rng rng(cp.seed);
    for (int li = 0; li < levels; ++li) {
        const double level = cp.snr_levels[static_cast<std::size_t>(li)];
        std::vector<Injection> inj;
        std::vector<SimMover> truth;
        for (int k = 0; k < cp.per_level; ++k) {
            Injection o;
            bool placed = false;
            for (int attempt = 0; attempt < 500 && !placed; ++attempt) {
                o.x_ref = rng.uniform(margin, w - margin);
                o.y_ref = rng.uniform(margin, h - margin);
                placed = true;
                for (const auto& q : inj)
                    if (std::hypot(q.x_ref - o.x_ref, q.y_ref - o.y_ref) < cp.min_separation) placed = false;
                // stay clear of objects already in the data
                for (const auto& c : base.candidates)
                    if (c.ok && std::hypot(c.x - o.x_ref, c.y - o.y_ref) < cp.min_separation) placed = false;
            }
            if (!placed) continue;
            const double vlo = std::max(params.vmin_pxh, 0.0) * 1.05, vhi = params.vmax_pxh * 0.95;
            const double v = rng.uniform(vlo, std::max(vlo, vhi));
            const double ang = rng.uniform(0.0, 2.0 * kPi);
            o.vx_pxh = v * std::cos(ang);
            o.vy_pxh = v * std::sin(ang);
            o.target_snr = level;
            o.flux = flux_for_snr(level, sigma, fwhm, n);
            inj.push_back(o);
            SimMover m;
            m.x_ref = o.x_ref; m.y_ref = o.y_ref; m.vx_pxh = o.vx_pxh; m.vy_pxh = o.vy_pxh;
            truth.push_back(m);
        }
        std::vector<Frame> work = frames;
        inject(work, inj, base_prep.t_ref_mjd, fwhm);
        const SearchResult r = search(work, fixed, sub(1 + li));
        const double max_dv = std::max(3.0, 2.0 * r.dv_pxh);
        const auto ok_m = match(r.candidates, truth, 2.0, max_dv, true);
        const auto any_m = match(r.candidates, truth, 2.0, max_dv, false);
        CompletenessRow row;
        row.target_snr = level;
        row.injected = static_cast<int>(inj.size());
        row.recovered = static_cast<int>(ok_m.size());
        row.detected = static_cast<int>(any_m.size());
        for (const auto& m : ok_m) {
            row.mean_snr += r.candidates[static_cast<std::size_t>(m.candidate_index)].snr;
            row.mean_dpos += m.dpos_px;
            row.mean_dv += m.dv_pxh;
        }
        if (!ok_m.empty()) {
            row.mean_snr /= ok_m.size();
            row.mean_dpos /= ok_m.size();
            row.mean_dv /= ok_m.size();
        }
        rep.rows.push_back(row);
    }

    // S/N at 50 % completeness by linear interpolation
    for (std::size_t i = 1; i < rep.rows.size(); ++i) {
        const auto& a = rep.rows[i - 1];
        const auto& b = rep.rows[i];
        if (a.injected == 0 || b.injected == 0) continue;
        const double fa = static_cast<double>(a.recovered) / a.injected;
        const double fb = static_cast<double>(b.recovered) / b.injected;
        if (fa < 0.5 && fb >= 0.5) {
            rep.snr50 = a.target_snr + (0.5 - fa) / (fb - fa) * (b.target_snr - a.target_snr);
            break;
        }
    }
    if (rep.snr50 == 0.0 && !rep.rows.empty() && rep.rows.front().injected > 0 &&
        static_cast<double>(rep.rows.front().recovered) / rep.rows.front().injected >= 0.5)
        rep.snr50 = rep.rows.front().target_snr;
    if (progress) progress(1.0, "done");
    return rep;
}

}  // namespace fainttrail
