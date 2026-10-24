// SPDX-License-Identifier: GPL-3.0-or-later
// Copyright 2026 Jurgen Kobierczynski
//
// Validation tools: injection-recovery (completeness) on any data set, real or
// simulated, and matching of candidates against known truth.
#pragma once

#include <cstdint>
#include <string>
#include <vector>

#include "fainttrail/search.hpp"
#include "fainttrail/simulate.hpp"

namespace fainttrail {

// A fake object to add to frames.
struct Injection {
    double x_ref = 0.0, y_ref = 0.0;   // 0-based at t_ref of the sequence
    double vx_pxh = 0.0, vy_pxh = 0.0;
    double flux = 0.0;                 // data units per frame
    double target_snr = 0.0;           // stack S/N it was scaled for (informational)
};

// Add PSF-shaped, trailed sources to the frames (in place).
// t_ref_mjd: reference epoch of x_ref/y_ref.
void inject(std::vector<Frame>& frames, const std::vector<Injection>& objects, double t_ref_mjd,
            double fwhm_px);

// Flux (data units per frame) that gives a stacked matched-filter S/N of
// `stack_snr` over `nframes` frames with per-pixel noise `sigma`.
double flux_for_snr(double stack_snr, double sigma, double fwhm_px, int nframes);

struct MatchResult {
    int truth_index = -1;
    int candidate_index = -1;
    double dpos_px = 0.0;     // position error at t_ref
    double dv_pxh = 0.0;      // velocity error
};

// Match candidates to known movers: within max_dpos pixels at t_ref and
// max_dv pixels/hour in velocity. Each truth object matches at most one candidate.
std::vector<MatchResult> match(const std::vector<Candidate>& cands, const std::vector<SimMover>& truth,
                               double max_dpos = 2.0, double max_dv = 3.0, bool ok_only = true);

struct CompletenessParams {
    std::vector<double> snr_levels{4.0, 5.0, 6.0, 7.0, 8.0, 10.0, 12.0, 15.0};
    int per_level = 10;          // injected objects per S/N level (one search per level)
    double min_separation = 40.0;
    std::uint64_t seed = 99;
};

struct CompletenessRow {
    double target_snr = 0.0;
    int injected = 0;
    int recovered = 0;           // found and passed every test
    int detected = 0;            // found at all (also flagged ones)
    double mean_snr = 0.0;       // mean refined S/N of the recovered ones
    double mean_dpos = 0.0;
    double mean_dv = 0.0;
};

struct CompletenessReport {
    std::vector<CompletenessRow> rows;
    double threshold_snr = 0.0;  // fixed for all injection runs
    std::string threshold_method;
    int baseline_ok = 0;         // candidates that passed in the data without injections
    double snr50 = 0.0;          // interpolated S/N at 50 % recovery (0 if not reached)
};

// Injection-recovery test. The threshold is calibrated once on the data
// without injections, then kept fixed so that every level is comparable.
// Objects already present in the data are counted in baseline_ok and do not
// affect the recovery numbers (matching only considers injected tracks).
CompletenessReport completeness(const std::vector<Frame>& frames, const SearchParams& params,
                                const CompletenessParams& cp, const ProgressFn& progress = {});

}  // namespace fainttrail
