// SPDX-License-Identifier: GPL-3.0-or-later
// Copyright 2026 Jurgen Kobierczynski
#pragma once

#include <string>
#include <vector>

namespace fainttrail::detail {

// Run `exe` (a full path) with `args`, no shell; stdout and stderr go to
// log_path. Returns the exit code, -1 if it could not be started, -2 on timeout.
int run_program(const std::string& exe, const std::vector<std::string>& args, int timeout_s, const std::string& log_path);

// First of `names` found in extra_dirs, then on the PATH; "" if none.
std::string find_program(const std::vector<std::string>& names, const std::vector<std::string>& extra_dirs = {});

}  // namespace fainttrail::detail
