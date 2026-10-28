// SPDX-License-Identifier: GPL-3.0-or-later
// Copyright 2026 Jurgen Kobierczynski
//
// Host-side planning for the CUDA backend: the shift table, 1/sqrt(c), the row
// bands that fit into device memory and the motion chunks per launch. Compiled
// into every build (plain C++), so the unit tests can run the complete GPU
// algorithm on the CPU through emulate_cuda_pass(), even without a GPU or CUDA.
#pragma once

#include <string>
#include <vector>

#include "cuda_stack_core.hpp"
#include "fainttrail/backend.hpp"

namespace fainttrail::detail {

struct CudaPlan {
    int w = 0, h = 0, n = 0, nv = 0, min_count = 1;
    std::vector<int> dx, dy;        // motion-major: j * n + i
    std::vector<float> inv_sqrt;    // c = 0..n
    int min_dy = 0, max_dy = 0;
    std::vector<cuda_core::Band> bands;
    int chunk = 0;                  // motions per kernel launch (multiple of kMotionBlock)
    double frame_bytes = 0;         // device bytes of the largest band's frame rows
};

// Options of the CUDA backend, from a name such as "cuda:device=1,max_mb=64".
struct CudaOptions {
    int device = -1;           // -1: FAINTTRAIL_CUDA_DEVICE or 0
    double max_mb = 0;         // cap on device memory for frames, 0 = 80 % of free memory
    double launch_adds = 4e9;  // additions per kernel launch (progress / cancel granularity)
};
CudaOptions parse_cuda_options(const std::string& spec);  // the part after "cuda:"

// budget_bytes: device memory available for the frame rows of one band.
CudaPlan make_cuda_plan(const StackPass& pass, double budget_bytes, double launch_adds);

// The CUDA algorithm executed on the host, band by band and launch by launch,
// with the same per-pixel code (cuda_core::stack_pixel) the device runs.
StackResult emulate_cuda_pass(const StackPass& pass, double budget_bytes, double launch_adds = 4e9);

// The CUDA backend proper (src/backend_cuda.cu), present when built with
// FAINTTRAIL_CUDA. cuda_device_count() returns 0 when no usable device or driver
// is found; `why` receives the reason.
#ifdef FAINTTRAIL_HAVE_CUDA
int cuda_device_count(std::string* why = nullptr);
std::unique_ptr<SearchBackend> make_cuda_backend(const CudaOptions& options);
#endif

}  // namespace fainttrail::detail
