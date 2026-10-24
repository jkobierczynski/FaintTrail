// SPDX-License-Identifier: GPL-3.0-or-later
// Copyright 2026 Jurgen Kobierczynski
//
// CPU reference backend (OpenMP, cache-blocked, auto-vectorised row loops).
#pragma once

#include "fainttrail/backend.hpp"

namespace fainttrail::detail {

// Block sizes of the cache-blocked kernel. 0 = choose automatically.
struct CpuBlocking {
    int band_rows = 0;   // output rows per band
    int vel_block = 0;   // trial motions accumulated together per frame
};

// Block sizes the kernel will use for an image of width w and height h on nt threads.
CpuBlocking choose_blocking(int w, int h, int nv, int nt, const CpuBlocking& requested = {});

class CpuBackend final : public SearchBackend {
public:
    explicit CpuBackend(CpuBlocking blocking = {}) : blocking_(blocking) {}
    std::string name() const override { return "cpu"; }
    std::string description() const override;
    StackResult run(const StackPass& pass, const ProgressFn& progress = {}) override;

private:
    CpuBlocking blocking_;
};

}  // namespace fainttrail::detail
