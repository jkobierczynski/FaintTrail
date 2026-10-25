// SPDX-License-Identifier: GPL-3.0-or-later
// Copyright 2026 Jurgen Kobierczynski
//
// Search backends: the interchangeable implementation of the shift-and-stack
// loop, which is where almost all of the search time goes.
//
// Every backend must satisfy the same contract, so that a faster one (CUDA,
// Vulkan compute, ...) can be checked pixel by pixel against the CPU reference.
//
// ---------------------------------------------------------------- contract
//
// Input (StackPass):
//   frames     n images g_i of w x h finite floats (unit-variance S/N images)
//   t_h        n frame times in hours relative to the reference epoch
//   velocities m trial motions (vx, vy) in pixels per hour
//   shifts     optional: an explicit table of integer shifts (dx, dy) for every
//              trial motion and frame (curved, ephemeris-guided or parallax-
//              corrected tracks); replaces the straight-line rule below
//   min_count  minimum number of frames that must cover a pixel
//
// Without a shift table, the integer shift of frame i for velocity v = (vx, vy) is
//   dx_i = round(vx * t_i),  dy_i = round(vy * t_i)
// where round() rounds half away from zero (std::lround), evaluated in double
// precision; linear_shift_table() builds exactly these. With a table, trial motion
// j uses its shifts (dx_ji, dy_ji) and "velocity index" below means j.
// Frame i covers output pixel (x, y) when (x + dx_i, y + dy_i) lies
// inside the frame. With C(x,y) the set of covering frames and c = |C|:
//   S_v(x, y) = ( sum_{i in C, ascending i} g_i(x + dx_i, y + dy_i) ) / sqrt(c)
// accumulated in single precision in ascending frame order, and undefined when
// c < min_count.
//
// Output (StackResult):
//   max_snr(x, y) = max over v of S_v(x, y)     (0 where S_v is undefined for all v)
//   best_v(x, y)  = index of the maximising v    (-1 where undefined)
//   Ties are broken towards the lowest velocity index.
//
// Tolerance: other backends may sum in a different order or use fused
// multiply-adds. They must match the reference within |d max_snr| <= 1e-4 *
// max(1, |max_snr|); best_v may differ only where the two velocities' S/N values
// are within that tolerance of each other. check_backend() tests exactly this.
#pragma once

#include <cmath>
#include <functional>
#include <memory>
#include <string>
#include <utility>
#include <vector>

#include "fainttrail/search.hpp"
#include "fainttrail/types.hpp"

namespace fainttrail {

// Integer shifts for every trial motion and frame: motion j, frame i at j * frames + i.
struct ShiftTable {
    int motions = 0, frames = 0;
    std::vector<int> dx, dy;
    void resize(int m, int n) {
        motions = m;
        frames = n;
        dx.assign(static_cast<std::size_t>(m) * n, 0);
        dy.assign(static_cast<std::size_t>(m) * n, 0);
    }
    int& x(int j, int i) { return dx[static_cast<std::size_t>(j) * frames + i]; }
    int& y(int j, int i) { return dy[static_cast<std::size_t>(j) * frames + i]; }
    int x(int j, int i) const { return dx[static_cast<std::size_t>(j) * frames + i]; }
    int y(int j, int i) const { return dy[static_cast<std::size_t>(j) * frames + i]; }
};

struct StackPass {
    const std::vector<Image>* frames = nullptr;
    const std::vector<double>* t_h = nullptr;          // not needed with a shift table
    const std::vector<std::pair<double, double>>* velocities = nullptr;  // optional with a shift table
    const ShiftTable* shifts = nullptr;                // optional explicit shifts
    int min_count = 3;
    int threads = 0;  // CPU backends: 0 = all cores
};

// The straight-line shifts round(v * t) of the contract, as a table.
ShiftTable linear_shift_table(const std::vector<double>& t_h,
                              const std::vector<std::pair<double, double>>& velocities);

// Number of trial motions of a pass, and the shift of frame i for motion j
// (from the table if there is one, otherwise round(v * t)). For backends.
int pass_motions(const StackPass& pass);
inline void pass_shift(const StackPass& pass, int j, int i, int& dx, int& dy) {
    if (pass.shifts) {
        dx = pass.shifts->x(j, i);
        dy = pass.shifts->y(j, i);
    } else {
        const auto& v = (*pass.velocities)[static_cast<std::size_t>(j)];
        const double t = (*pass.t_h)[static_cast<std::size_t>(i)];
        dx = static_cast<int>(std::lround(v.first * t));
        dy = static_cast<int>(std::lround(v.second * t));
    }
}

struct StackResult {
    Image max_snr;
    std::vector<int> best_v;
    // Optional breakdown of the pass, for profiling. Backends report what they can
    // (the CPU backend: accumulate / coverage+max per thread, averaged; merge).
    std::vector<StageTime> timings;
};

class SearchBackend {
public:
    virtual ~SearchBackend() = default;
    // Short identifier used on the command line, e.g. "cpu".
    virtual std::string name() const = 0;
    // Human-readable description, e.g. device name and memory.
    virtual std::string description() const = 0;
    // Run one pass. The progress callback receives the fraction of velocities
    // done (0..1); returning false cancels, after which Cancelled is thrown.
    virtual StackResult run(const StackPass& pass, const ProgressFn& progress = {}) = 0;
};

// Names of the backends compiled into this build, best first, followed by the
// registered ones.
std::vector<std::string> backend_names();

// Register an additional backend under `name` (experimental kernels, studies,
// out-of-tree GPU backends). "auto" never selects a registered backend. A
// registered backend that does not satisfy the contract will fail check-backend.
using BackendFactory = std::function<std::unique_ptr<SearchBackend>()>;
void register_backend(const std::string& name, BackendFactory factory);

// Create a backend by name; "auto" (or "") picks the best available one.
// Throws std::runtime_error for unknown or unavailable backends.
std::unique_ptr<SearchBackend> make_backend(const std::string& name = "auto");

// Validate the inputs of a pass against the contract; throws std::runtime_error.
void validate_pass(const StackPass& pass);

// Plain, unoptimised implementation of the contract: one velocity, one pixel,
// one frame at a time. Far too slow for real data; it exists so the optimised
// backends can be tested against the specification itself.
StackResult run_pass_naive(const StackPass& pass);

// Simple full-image implementation of the contract (one motion at a time over
// the whole image), independent of the optimised CPU backend; bit-identical to
// the naive one and fast enough to check other backends on realistic sizes.
StackResult run_pass_simple(const StackPass& pass);

struct BackendCheck {
    bool ok = false;
    double max_abs_diff = 0.0;      // largest |d max_snr|
    std::size_t pixels = 0;
    std::size_t best_v_mismatch = 0;     // pixels where best_v differs
    std::size_t best_v_unexplained = 0;  // ... and the S/N values are not tied
    std::string message;
};

// Compare a candidate backend's result with a reference result for the same
// pass, following the tolerance rules of the contract.
BackendCheck check_backend(const StackPass& pass, const StackResult& reference,
                           const StackResult& candidate);

}  // namespace fainttrail
