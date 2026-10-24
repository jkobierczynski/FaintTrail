// SPDX-License-Identifier: GPL-3.0-or-later
// Copyright 2026 Jurgen Kobierczynski
//
// Internal entry points shared between implementation files.
#pragma once

#include <utility>
#include <vector>

#include "fainttrail/search.hpp"

namespace fainttrail::detail {

// Time-order dependent preprocessing (static sky ... matched filter).
void build_stage(const std::vector<Image>& sub, Prepared& p, const SearchParams& params,
                 const ProgressFn& progress, bool estimate_psf);

// Static sky of one pixel (median over the frames) and every frame's residual
// against the leave-out median that excludes the frames within `exclude` ranks
// of it in time. Returns the static sky (NaN where too few frames cover the pixel).
float static_sky_column(const std::vector<Image>& sub, std::size_t idx, const std::vector<int>& rank, int exclude,
                        std::vector<std::pair<float, int>>& sorted, std::vector<float>& vals, float* residual);

// Residuals of one frame in a rectangle (0-based origin x0, y0), NaN outside the
// image. Copied when Prepared keeps residuals, otherwise recomputed exactly.
Image residual_patch(const Prepared& p, int frame, int x0, int y0, int w, int h);

// Re-run the time-dependent preprocessing with different frame times
// (used for the time-permuted null searches).
Prepared with_times(const Prepared& base, const std::vector<double>& t_h, const SearchParams& params);

}  // namespace fainttrail::detail
