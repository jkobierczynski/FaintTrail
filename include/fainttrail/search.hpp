// SPDX-License-Identifier: GPL-3.0-or-later
// Copyright 2026 Jurgen Kobierczynski
//
// Synthetic tracking: blind shift-and-stack search over a grid of linear
// motions, with an empirical false-alarm calibration and per-candidate
// consistency tests.
#pragma once

#include <cstdint>
#include <string>
#include <vector>

#include "fainttrail/registration.hpp"
#include "fainttrail/types.hpp"

namespace fainttrail {

struct SearchParams {
    // Motion range in pixels per hour. Objects slower than vmin are not searched
    // (they self-subtract in the static-sky model anyway).
    double vmin_pxh = 4.0;
    double vmax_pxh = 40.0;
    // Maximum positional drift (px) a wrong grid velocity may cause over half
    // the sequence; sets the grid spacing dv = 4 * grid_tol_px / span.
    double grid_tol_px = 0.5;
    // PSF FWHM in pixels; 0 = estimate from the stars.
    double fwhm_px = 0.0;
    // Residuals are clipped at this many sigma before stacking, which keeps
    // cosmic rays and other single-frame defects from dominating a stack.
    double clip_sigma = 8.0;
    // Pixels at least this many sigma above the residual noise that are much
    // sharper than the PSF are treated as cosmic rays and zeroed in that frame.
    double cr_sigma = 5.0;
    // Frames on each side of frame i left out of its static-sky median, so that
    // moving objects do not subtract themselves; -1 = automatic (N/8).
    int static_exclude = -1;
    // Pixels whose temporal scatter exceeds this many times the frame noise
    // (stars, bad columns) are masked.
    double mask_scatter = 2.5;
    // Pixels with temporal scatter above this (but below mask_scatter) are
    // down-weighted by their scatter.
    double weight_scatter = 1.5;
    // A stack pixel needs at least this fraction of frames overlapping.
    double min_coverage = 0.7;
    // Detection threshold on the coarse-grid S/N; 0 = calibrate automatically.
    double threshold_snr = 0.0;
    // Lower bound for an automatic threshold.
    double min_snr_floor = 5.0;
    // Automatic calibration: number of time-shuffled null searches and the
    // accepted number of noise detections per search.
    int null_runs = 3;
    double target_false_per_run = 0.02;
    // Tail model fitted to the null peaks (docs/CALIBRATION.md):
    //   "gauss"  log N(>s) = a + k log s - s^2 / (2 sigma^2), the tail of the peaks of
    //            a smooth Gaussian random field, fitted to the top 3000 peaks per run
    //   "linear" log N(>s) = a - b s fitted to the peaks between 0.67 and 30 per run
    //            (FaintTrail <= 0.7; conservative by a data-dependent amount)
    std::string null_model = "gauss";
    // gauss model: the strongest null peaks of each run, as many as the fit uses, are
    // checked with the artefact tests (edge, detector, single) at their grid position
    // and motion, and artefacts are left out of the fit, since the same tests reject
    // them in the real search. 0 = test none.
    int null_tested = 3000;
    // Consistency tests.
    double max_frame_fraction = 0.25; // max share of the signal from one frame
    double half_tolerance = 3.0;      // allowed shortfall (sigma) of each half
    int max_candidates = 300;
    int edge_px = 4;
    int threads = 0;                  // 0 = all available (CPU backend)
    // Stacking backend: "auto" (best available), "cpu", ... see backend.hpp.
    std::string backend = "auto";
    // Keep the residual images (frame - background - static sky) in Prepared after
    // preparation. Needed only for residual displays, stack_full() and
    // stack_cutout(snr = false); costs one more copy of the frames in memory.
    bool keep_residuals = false;
    // Registration of the input frames (registration.hpp):
    //   "auto"   check; align the frames when the stars move (> 1.5 px) or the sizes differ
    //   "always" align the frames on their middle frame
    //   "check"  check; refuse unregistered frames
    //   "off"    trust the frames as they are
    std::string align = "auto";
    std::uint64_t seed = 12345;
};

// Preprocessed, noise-normalised data; reusable for several searches,
// cutouts and the GUI.
// Wall-clock time of one pipeline stage (for profiling; see `--profile`).
struct StageTime {
    std::string name;
    double seconds = 0.0;
    double rss_mb = 0.0;  // resident memory of the process at the end of the stage (0 = not measured)
};

struct Prepared {
    int width = 0, height = 0, nframes = 0;
    double t_ref_mjd = 0.0;             // reference epoch (middle of the sequence)
    double span_h = 0.0;                // last minus first mid-exposure time, hours
    std::vector<double> t_h;            // frame times relative to t_ref, hours
    std::vector<double> sigma;          // per-frame residual noise, data units
    std::vector<double> background;     // per-frame background level, data units
    std::vector<Image> sub;             // frame - background (data units), kept for null runs
    std::vector<Image> residual;        // frame - background - static sky (data units); empty unless keep_residuals
    std::vector<Image> filtered;        // matched-filtered S/N images (unit variance)
    Image static_sky;                   // median of background-subtracted frames
    Image mask;                         // 1 = usable, 0 = masked
    std::vector<int> cr_pixels;         // cosmic-ray hits removed per frame
    double fwhm_px = 0.0;
    Wcs wcs;
    std::vector<double> exptime_s;
    std::vector<double> reg_dx, reg_dy; // registration offsets, empty if unknown
    std::vector<StageTime> timings;     // preprocessing stages, in order
    const std::vector<Image>* shared_sub = nullptr;  // null runs: the real data's `sub` (not owned)
    std::vector<int> time_rank;         // rank of each frame in time (for recomputing residuals)
    int static_exclude = 0;             // frames left out on each side by the leave-out median
    DriftCheck registration;            // registration check of the input frames
    AlignmentReport alignment;          // what register_frames() did, if it ran
};

struct Candidate {
    int id = 0;
    double x = 0.0, y = 0.0;            // 0-based position at t_ref
    double vx_pxh = 0.0, vy_pxh = 0.0;  // motion
    double pos_err_px = 0.0;            // 1-sigma position uncertainty per axis at t_ref
    double vel_err_pxh = 0.0;           // 1-sigma motion uncertainty per axis
    double snr = 0.0;                   // refined S/N (sub-pixel track)
    double snr_coarse = 0.0;            // grid S/N, the statistic the threshold applies to
    double snr_half1 = 0.0, snr_half2 = 0.0;  // first / second half of the frames
    double max_frame_fraction = 0.0;    // largest single-frame share of the signal
    double false_alarms_per_run = 0.0;  // expected noise detections at this S/N
    double flux = 0.0;                  // PSF flux in data units per frame
    double mag_inst = 0.0;              // -2.5 log10(flux) + 25
    double ra_deg = 0.0, dec_deg = 0.0; // at t_ref, if a WCS is present
    double mag = 99.0, mag_err = 0.0;   // calibrated magnitude in SearchResult::mag_band (photometry.hpp)
    double rate_arcsec_min = 0.0;
    double pa_deg = 0.0;                // position angle of motion, east of north
    bool ok = false;                    // passed every test
    std::string flags;                  // reasons for rejection, comma separated
};

struct TrackPoint {
    int frame = 0;
    double mjd = 0.0;
    double x = 0.0, y = 0.0;            // 0-based
    double ra_deg = 0.0, dec_deg = 0.0;
};

struct SearchResult {
    std::vector<Candidate> candidates; // sorted by S/N, ok ones first
    double threshold_snr = 0.0;
    std::string threshold_method;       // "fixed", "null-fit gauss", "null-fit linear", "null-quantile", "gaussian"
    // Null-tail model: expected noise peaks per run above s is
    //   exp(null_a + null_k log s - null_inv_s2 s^2 / 2)   (gauss model)
    //   exp(null_a - null_b s)                            (linear model)
    // see expected_null_peaks().
    double null_a = 0.0, null_b = 0.0;
    double null_k = 0.0, null_inv_s2 = 0.0;
    std::vector<double> null_peaks;     // pooled peak S/N of the null runs
    std::size_t n_velocities = 0;
    double dv_pxh = 0.0;
    double fwhm_px = 0.0;
    double t_ref_mjd = 0.0;
    bool has_sky = false;               // candidates carry RA/Dec (a plate solution was available)
    char mag_band = 0;                  // 'V', 'R' or 'G' once calibrated (photometry.hpp), else 0
    double mag_zero_point = 0.0, mag_zp_error = 0.0;
    int mag_stars = 0;
    // Motions slower than this move less than two PSF widths during the
    // sequence and are partly absorbed by the static-sky model.
    double slowest_reliable_pxh = 0.0;
    std::string backend;                // backend that ran the stacking passes
    std::string backend_description;
    double elapsed_s = 0.0;
    std::vector<StageTime> timings;     // search stages, in order (null runs include their rebuild)
    Image max_snr;                      // max over all velocities, per pixel
};

// Preprocess frames: background removal, static-sky subtraction, noise
// normalisation, masking and matched filtering. Frames must be registered.
Prepared prepare(const std::vector<Frame>& frames, const SearchParams& params,
                 const ProgressFn& progress = {});

// Expected number of noise peaks per search above `snr`, from the fitted null
// tail (0 if the threshold was not calibrated with null runs).
double expected_null_peaks(const SearchResult& res, double snr);

// For applications (not libraries embedding FaintTrail): call once at start-up.
// With glibc, frame-sized buffers are then always mapped separately, so memory
// freed between stages goes back to the system at once instead of staying with
// the process (about 30 % lower peak memory on large sequences). No-op elsewhere.
void tune_allocator_for_images();

// Same, but takes over the frames' pixel buffers instead of copying them, which
// saves one copy of the sequence in memory. The frames' images are left empty.
Prepared prepare(std::vector<Frame>&& frames, const SearchParams& params,
                 const ProgressFn& progress = {});

// Full blind search, including the null calibration if threshold_snr == 0.
SearchResult search(const Prepared& prep, const SearchParams& params,
                    const ProgressFn& progress = {});

// Convenience: prepare + search.
SearchResult search(const std::vector<Frame>& frames, const SearchParams& params,
                    const ProgressFn& progress = {});

// Velocity grid used by the search.
std::vector<std::pair<double, double>> velocity_grid(const Prepared& prep, const SearchParams& params,
                                                     double* dv_out = nullptr);

// Shift-and-stack along a known motion. `which`: 0 = all frames, 1 = first
// half, 2 = second half. With `snr` true the matched-filtered S/N frames are
// stacked (result in sigma), otherwise the residuals (mean, data units).
Image stack_cutout(const Prepared& prep, double x, double y, double vx_pxh, double vy_pxh,
                   int half_size, bool snr = false, int which = 0);

// Full-frame shift-and-stack of the residuals for a known motion (ephemeris).
Image stack_full(const Prepared& prep, double vx_pxh, double vy_pxh);

// Per-frame positions of a candidate.
std::vector<TrackPoint> track_points(const Prepared& prep, const Candidate& c);

// Estimate the stellar FWHM (px) from an image with stars; returns 0 on failure.
double estimate_fwhm(const Image& img, double noise_sigma);

}  // namespace fainttrail
