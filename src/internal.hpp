// SPDX-License-Identifier: GPL-3.0-or-later
// Copyright 2026 Jurgen Kobierczynski
//
// Internal entry points shared between implementation files.
#pragma once

#include <utility>
#include <vector>

#include <cstdint>
#include <string>

#include "fainttrail/backend.hpp"
#include "fainttrail/search.hpp"

namespace fainttrail::detail {

// Positions sample_pixels() uses on a full_w x full_h image (regular grid with
// the given border), so that samples gathered tile by tile are the same set.
struct SampleGrid {
    int full_w = 0, full_h = 0, border = 0, step = 1;
    bool hit(int gx, int gy) const;
};
SampleGrid sample_grid(int full_w, int full_h, std::size_t max_samples, int border);

// Residual noise of a frame from its residual samples (1 if undefined).
double sigma_from_samples(std::vector<float>& samples);

// What build_stage computes.
enum class StageMode {
    Full,           // everything up to the normalised matched-filtered frames
    FilterSamples,  // stop before normalising; collect matched-filter samples (for the normalisation)
};

// Where the prepared region lies in the full image (tiled search). The core is the
// part whose results are used; the rest is halo.
struct StageRegion {
    int ox = 0, oy = 0;          // region origin in the full image
    int full_w = 0, full_h = 0;  // full image size
    int cx0 = 0, cy0 = 0, cx1 = 0, cy1 = 0;  // core, region coordinates, [cx0,cx1) x [cy0,cy1)
};

// Per-frame quantities measured on the whole image.
struct StageGlobals {
    const std::vector<double>* sigma = nullptr;      // residual noise (sigma_from_samples)
    const std::vector<double>* filter_sg = nullptr;  // robust sigma of the unnormalised matched filter
};

struct StageSamples {
    std::vector<std::vector<float>> per_frame;
};

// Time-order dependent preprocessing (static sky ... matched filter). Without a
// region, the whole image is processed and every global quantity estimated from it.
void build_stage(const std::vector<Image>& sub, Prepared& p, const SearchParams& params,
                 const ProgressFn& progress, bool estimate_psf, StageMode mode = StageMode::Full,
                 const StageRegion* region = nullptr, const StageGlobals* globals = nullptr,
                 StageSamples* samples = nullptr);

// Background-subtract a raw frame (NaN stays NaN) and grow blank areas by 3 pixels,
// as prepare() does. For a region, read 3 extra pixels on each side and crop them.
void subtract_and_grow(Image& im, double background);
constexpr int kBlankGrow = 3;

// PSF width from the static sky (as in prepare()).
double estimate_psf_fwhm(const Image& static_sky, std::vector<double> sigma, const SearchParams& params);

// Static sky of one pixel (median over the frames) and every frame's residual
// against the leave-out median that excludes the frames within `exclude` ranks
// of it in time. Returns the static sky (NaN where too few frames cover the pixel).
// With residual == nullptr only the static sky is computed.
float static_sky_column(const std::vector<Image>& sub, std::size_t idx, const std::vector<int>& rank, int exclude,
                        std::vector<std::pair<float, int>>& sorted, std::vector<float>& vals, float* residual);

// Residuals of one frame in a rectangle (0-based origin x0, y0), NaN outside the
// image. Copied when Prepared keeps residuals, otherwise recomputed exactly.
Image residual_patch(const Prepared& p, int frame, int x0, int y0, int w, int h);

// ---------------------------------------------------------------- search steps
// (shared by search() and the tiled search)

struct Peak {
    int x = 0, y = 0;
    float snr = 0.0f;
    int vidx = -1;
};
using PassOutput = StackResult;

struct TrackEval {
    double snr = 0.0;
    int count = 0;
};

struct RefinedCandidate {
    Candidate c;
    std::vector<double> contrib;  // per-frame track values (for the ghost test)
};

PassOutput run_pass(SearchBackend& backend, const Prepared& p, const std::vector<std::pair<double, double>>& grid,
                    const std::vector<double>& t, const SearchParams& params, const ProgressFn& progress, double f0,
                    double f1, const std::string& stage);
std::vector<Peak> find_peaks(const PassOutput& pass, float floor, int radius, int edge, std::size_t max_peaks);
TrackEval track_snr(const Prepared& p, double x, double y, double vx, double vy, int which,
                    std::vector<double>* per_frame = nullptr);
Candidate refine(const Prepared& p, const Peak& pk, const std::pair<double, double>& v0, double dv,
                 const SearchParams& params);
RefinedCandidate refine_with_contrib(const Prepared& p, const Peak& pk, const std::pair<double, double>& v0, double dv,
                                     const SearchParams& params);
std::vector<std::vector<double>> null_time_orders(const std::vector<double>& t_h, int runs, std::uint64_t seed);
void set_threshold(const Prepared& p, const SearchParams& params, std::size_t n_velocities, int runs,
                   const std::vector<float>& pooled, SearchResult& res);
std::vector<Peak> select_peaks(const std::vector<Peak>& peaks, double threshold, const SearchParams& params);
void finish_candidates(const Prepared& meta, std::vector<RefinedCandidate> rc, SearchResult& res);
void sky_coordinates(const Prepared& p, Candidate& c);

// Re-run the time-dependent preprocessing with different frame times
// (used for the time-permuted null searches).
Prepared with_times(const Prepared& base, const std::vector<double>& t_h, const SearchParams& params);

}  // namespace fainttrail::detail
