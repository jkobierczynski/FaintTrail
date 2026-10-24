// SPDX-License-Identifier: GPL-3.0-or-later
// Copyright 2026 Jurgen Kobierczynski
//
// CPU reference backend (OpenMP + auto-vectorised row loops).
#pragma once

#include "fainttrail/backend.hpp"

namespace fainttrail::detail {

class CpuBackend final : public SearchBackend {
public:
    std::string name() const override { return "cpu"; }
    std::string description() const override;
    StackResult run(const StackPass& pass, const ProgressFn& progress = {}) override;
};

}  // namespace fainttrail::detail
