// SPDX-License-Identifier: GPL-3.0-or-later
// Copyright 2026 Jurgen Kobierczynski
//
// Real-data regression cases.
//
// A case is a small text file that names a directory of registered frames,
// the search settings, and what the search must find there:
//
//   * known objects (asteroids identified independently, e.g. by MPChecker),
//     each with a position at a given epoch, a motion and tolerances. These are
//     hard requirements: a missing or rejected known object fails the case.
//   * limits on everything else (number of unexplained passing candidates,
//     threshold range).
//   * optionally a baseline: the complete candidate list of an earlier run.
//     Differences against it are reported as "changes"; they fail the case only
//     in strict mode. A baseline catches silent drift that the truth checks do
//     not (a new spurious detection that happens to stay under the limits, a
//     shifted threshold, a different backend that sums differently).
//
// Case file format (INI style, '#' starts a comment):
//
//   name        = leo-triplet-cut
//   description = (21996) 1999 XP31 in a crop of the Leo Triplet sequence
//   data        = data/leo-triplet-cut        # relative to the case file;
//                                             # ${VAR} expands environment variables
//   baseline    = leo-triplet-cut.baseline.csv  # optional
//
//   [search]      vmin, vmax, threshold, null_runs, false_per_run, fwhm, grid_tol,
//                 static_exclude, mask_scatter, weight_scatter, seed, wcs_from, wcs_scale
//   [expect]      max_unexpected (default 0), threshold_min, threshold_max
//   [object NAME] mjd, x, y (FITS 1-based pixels at mjd), vx, vy (px/h),
//                 pos_tol (px, default 1.5), vel_tol (px/h, default 3), min_snr,
//                 ra, dec (deg) and sky_tol (arcsec, default 5), source (free text)
//
// If the data directory does not exist the case is skipped, not failed, so the
// large full-frame cases can live next to the small fixtures and run only where
// the data are available.
#pragma once

#include <string>
#include <vector>

#include "fainttrail/search.hpp"

namespace fainttrail {

struct RegressionObject {
    std::string name;
    std::string source;      // where the truth comes from
    double mjd = 0.0;
    double x = 0.0, y = 0.0; // FITS 1-based at mjd
    double vx = 0.0, vy = 0.0;
    double pos_tol = 1.5, vel_tol = 3.0;
    double min_snr = 0.0;
    bool has_sky = false;
    double ra = 0.0, dec = 0.0, sky_tol_arcsec = 5.0;
};

struct RegressionCase {
    std::string name, description;
    std::string case_path;       // the .case file
    std::string data_dir;        // resolved
    std::string baseline_path;   // resolved, empty if none
    std::string wcs_from;        // resolved, empty if none
    double wcs_scale = 0.0;
    SearchParams params;
    int max_unexpected = 0;
    double threshold_min = 0.0, threshold_max = 0.0;  // 0 = unchecked
    std::vector<RegressionObject> objects;
};

RegressionCase load_regression_case(const std::string& path);

enum class RegressionStatus { Pass, Changed, Fail, Skipped };
const char* status_name(RegressionStatus s);

struct RegressionCheck {
    bool ok = true;
    bool hard = true;   // hard checks fail the case; soft ones only mark it changed
    std::string text;
};

struct RegressionOutcome {
    RegressionStatus status = RegressionStatus::Pass;
    std::string skip_reason;
    std::vector<RegressionCheck> checks;
    SearchResult result;
    std::size_t nframes = 0;
    int width = 0, height = 0;
    double load_s = 0.0, prepare_s = 0.0, search_s = 0.0;
};

struct RegressionOptions {
    std::string backend;  // empty: as in the case (normally auto)
    int threads = 0;
    bool strict = false;  // baseline changes fail the case
};

RegressionOutcome run_regression_case(const RegressionCase& c, const RegressionOptions& opt = {},
                                      const ProgressFn& progress = {});

// Baseline files: the candidate list of one run, with the settings that matter
// for comparing it. Written by `fainttrail regress --write-baseline`.
void write_regression_baseline(const std::string& path, const RegressionCase& c, const RegressionOutcome& o);

}  // namespace fainttrail
