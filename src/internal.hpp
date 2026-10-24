// SPDX-License-Identifier: GPL-3.0-or-later
// Copyright 2026 Jurgen Kobierczynski
//
// Internal entry points shared between implementation files.
#pragma once

#include <vector>

#include "fainttrail/search.hpp"

namespace fainttrail::detail {

// Time-order dependent preprocessing (static sky ... matched filter).
void build_stage(const std::vector<Image>& sub, Prepared& p, const SearchParams& params,
                 const ProgressFn& progress, bool estimate_psf);

// Re-run the time-dependent preprocessing with different frame times
// (used for the time-permuted null searches).
Prepared with_times(const Prepared& base, const std::vector<double>& t_h, const SearchParams& params);

}  // namespace fainttrail::detail
