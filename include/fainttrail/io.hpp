// SPDX-License-Identifier: GPL-3.0-or-later
// Copyright 2026 Jurgen Kobierczynski
//
// CSV input/output. All pixel coordinates in files use the FITS convention
// (1-based, centre of the first pixel is 1,1).
#pragma once

#include <string>
#include <vector>

#include "fainttrail/search.hpp"
#include "fainttrail/simulate.hpp"
#include "fainttrail/validate.hpp"

namespace fainttrail {

void write_candidates_csv(const std::string& path, const SearchResult& res);
std::vector<Candidate> read_candidates_csv(const std::string& path);

// One line per candidate per frame: the tracklet a host turns into ADES records.
void write_tracklets_csv(const std::string& path, const Prepared& prep, const SearchResult& res,
                         bool ok_only = true);

void write_truth_csv(const std::string& path, const SimResult& sim);
std::vector<SimTruth> read_truth_csv(const std::string& path);

void write_completeness_csv(const std::string& path, const CompletenessReport& rep);

// Human-readable summary of a search.
std::string summary_text(const SearchResult& res, bool ok_only = false);

}  // namespace fainttrail
