// SPDX-License-Identifier: GPL-3.0-or-later
// Copyright 2026 Jurgen Kobierczynski
//
// The arithmetic of the CUDA backend, shared by the device kernel and a host
// emulation (src/cuda_emulation.cpp) that the unit tests run on every machine.
// Keeping one copy means the tests check the exact code the GPU executes, minus
// the CUDA plumbing (memory transfers, launches).
//
// Layout on the device: the output is processed in horizontal bands. For a band
// of output rows [y0, y1) the frames are present only for the rows the shifts can
// reach, [r0, r0 + rows); every frame is stored as `rows` x `w` floats, frame after
// frame. One thread owns one output pixel and visits the trial motions in
// ascending order, in chunks; for each block of kMotionBlock motions it walks the
// frames in ascending order, so every sum is formed in exactly the order the
// contract prescribes, and the result is bit-identical to the CPU backend.
#pragma once

#include <algorithm>
#include <cstddef>
#include <vector>

#ifdef __CUDACC__
#define FT_HD __host__ __device__
#define FT_UNROLL _Pragma("unroll")
#else
#define FT_HD
#define FT_UNROLL
#endif

namespace fainttrail::detail::cuda_core {

constexpr int kMotionBlock = 8;  // motions accumulated together in registers

struct BandView {
    const float* frames = nullptr;  // n x rows x w, read-only during a pass
    int w = 0, h = 0;               // full frame size
    int r0 = 0, rows = 0;           // frame rows present: [r0, r0 + rows)
    int n = 0;                      // frames
    const int* dx = nullptr;        // motion j, frame i at j * n + i
    const int* dy = nullptr;
    const float* inv_sqrt = nullptr;  // 1 / sqrt(c), c = 0..n
    int min_count = 1;
};

// Running maximum over the motions [j0, j1) for output pixel (x, y).
FT_HD inline void stack_pixel(const BandView& b, int x, int y, int j0, int j1, float& best, int& arg) {
    for (int j = j0; j < j1; j += kMotionBlock) {
        const int m = (j1 - j) < kMotionBlock ? (j1 - j) : kMotionBlock;
        // fixed-size loops with a guard, fully unrolled on the device, keep acc and
        // cnt in registers
        float acc[kMotionBlock];
        int cnt[kMotionBlock];
        FT_UNROLL
        for (int k = 0; k < kMotionBlock; ++k) {
            acc[k] = 0.0f;
            cnt[k] = 0;
        }
        for (int i = 0; i < b.n; ++i) {  // ascending frame order, as the contract requires
            const float* f = b.frames + static_cast<std::size_t>(i) * b.rows * b.w;
            FT_UNROLL
            for (int k = 0; k < kMotionBlock; ++k) {
                if (k >= m) continue;
                const std::size_t s = static_cast<std::size_t>(j + k) * b.n + i;
                const int sx = x + b.dx[s], sy = y + b.dy[s];
                if (sx < 0 || sx >= b.w || sy < 0 || sy >= b.h) continue;  // frame does not cover
                acc[k] += f[static_cast<std::size_t>(sy - b.r0) * b.w + sx];
                ++cnt[k];
            }
        }
        FT_UNROLL
        for (int k = 0; k < kMotionBlock; ++k) {
            if (k >= m || cnt[k] < b.min_count) continue;
            const float snr = acc[k] * b.inv_sqrt[cnt[k]];
            if (snr > best) {  // motions ascend: ties keep the lowest index
                best = snr;
                arg = j + k;
            }
        }
    }
}

// ---------------------------------------------------------------- host planning
struct Band {
    int y0 = 0, y1 = 0;  // output rows
    int r0 = 0, r1 = 0;  // frame rows needed, [r0, r1)
};

// Cut the output into bands whose frame rows fit into `budget_bytes` on the
// device. min_dy / max_dy: the extreme vertical shifts of all motions and frames.
inline std::vector<Band> plan_bands(int w, int h, int n, int min_dy, int max_dy, double budget_bytes) {
    std::vector<Band> bands;
    const double per_row = 4.0 * w * n;
    const int halo = (max_dy > 0 ? max_dy : 0) + (min_dy < 0 ? -min_dy : 0);
    int out_rows = static_cast<int>(budget_bytes / per_row) - halo;
    out_rows = std::max(1, std::min(out_rows, h));
    for (int y0 = 0; y0 < h; y0 += out_rows) {
        Band b;
        b.y0 = y0;
        b.y1 = std::min(h, y0 + out_rows);
        b.r0 = std::max(0, std::min(h, b.y0 + min_dy));
        b.r1 = std::min(h, std::max(0, b.y1 - 1 + max_dy + 1));
        if (b.r1 <= b.r0) b.r1 = std::min(h, b.r0 + 1);  // every shift leaves the frame: keep a valid buffer
        bands.push_back(b);
    }
    return bands;
}

}  // namespace fainttrail::detail::cuda_core
