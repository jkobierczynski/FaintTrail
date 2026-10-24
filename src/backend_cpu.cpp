// SPDX-License-Identifier: GPL-3.0-or-later
// Copyright 2026 Jurgen Kobierczynski
//
// CPU reference implementation of the stacking contract (see backend.hpp).
//
// Velocities are distributed over threads. For one velocity the shifted frames
// are added row by row into a stack buffer (contiguous, vectorisable), and the
// number of covering frames per pixel comes from a 2-D difference array: each
// frame covers one rectangle, so coverage costs O(n) per velocity instead of
// O(n * pixels).
#include "backend_cpu.hpp"

#include <algorithm>
#include <atomic>
#include <cmath>
#include <limits>
#include <string>
#include <thread>

#include "common.hpp"

namespace fainttrail::detail {

std::string CpuBackend::description() const {
#ifdef FAINTTRAIL_HAVE_OPENMP
    return "CPU reference, OpenMP, " + std::to_string(thread_count(0)) + " threads";
#else
    return "CPU reference, single thread (built without OpenMP)";
#endif
}

StackResult CpuBackend::run(const StackPass& pass, const ProgressFn& progress) {
    validate_pass(pass);
    const auto& frames = *pass.frames;
    const auto& t = *pass.t_h;
    const auto& grid = *pass.velocities;
    const int w = frames[0].width, h = frames[0].height;
    const int n = static_cast<int>(frames.size());
    const int nv = static_cast<int>(grid.size());
    const int min_count = pass.min_count;
    const int nt = thread_count(pass.threads);

    std::vector<float> inv_sqrt(static_cast<std::size_t>(n) + 1, 0.0f);
    for (int c = 1; c <= n; ++c) inv_sqrt[static_cast<std::size_t>(c)] = static_cast<float>(1.0 / std::sqrt(static_cast<double>(c)));

    const float neg = -std::numeric_limits<float>::infinity();
    std::vector<Image> lmax(static_cast<std::size_t>(nt));
    std::vector<std::vector<int>> larg(static_cast<std::size_t>(nt));
    std::atomic<int> done{0};
    std::atomic<bool> cancel{false};

#pragma omp parallel num_threads(nt)
    {
        const int tid = thread_id();
        Image& M = lmax[static_cast<std::size_t>(tid)];
        std::vector<int>& A = larg[static_cast<std::size_t>(tid)];
        M = Image(w, h, neg);
        A.assign(static_cast<std::size_t>(w) * h, -1);
        std::vector<float> stack(static_cast<std::size_t>(w) * h);
        std::vector<int> diff(static_cast<std::size_t>(w + 1) * (h + 1));
        std::vector<int> colacc(static_cast<std::size_t>(w));

#pragma omp for schedule(dynamic, 2)
        for (int vi = 0; vi < nv; ++vi) {
            if (cancel.load(std::memory_order_relaxed)) continue;
            const double vx = grid[static_cast<std::size_t>(vi)].first;
            const double vy = grid[static_cast<std::size_t>(vi)].second;
            std::fill(stack.begin(), stack.end(), 0.0f);
            std::fill(diff.begin(), diff.end(), 0);
            for (int i = 0; i < n; ++i) {  // ascending frame order, as the contract requires
                const int dx = static_cast<int>(std::lround(vx * t[static_cast<std::size_t>(i)]));
                const int dy = static_cast<int>(std::lround(vy * t[static_cast<std::size_t>(i)]));
                const int x0 = std::max(0, -dx), x1 = std::min(w, w - dx);
                const int y0 = std::max(0, -dy), y1 = std::min(h, h - dy);
                if (x0 >= x1 || y0 >= y1) continue;
                diff[static_cast<std::size_t>(y0) * (w + 1) + x0] += 1;
                diff[static_cast<std::size_t>(y0) * (w + 1) + x1] -= 1;
                diff[static_cast<std::size_t>(y1) * (w + 1) + x0] -= 1;
                diff[static_cast<std::size_t>(y1) * (w + 1) + x1] += 1;
                const float* g = frames[static_cast<std::size_t>(i)].data.data();
                for (int y = y0; y < y1; ++y) {
                    float* o = &stack[static_cast<std::size_t>(y) * w];
                    const float* s = g + static_cast<std::ptrdiff_t>(y + dy) * w + dx;
                    for (int x = x0; x < x1; ++x) o[x] += s[x];
                }
            }
            std::fill(colacc.begin(), colacc.end(), 0);
            for (int y = 0; y < h; ++y) {
                int running = 0;
                const int* drow = &diff[static_cast<std::size_t>(y) * (w + 1)];
                const float* srow = &stack[static_cast<std::size_t>(y) * w];
                float* mrow = &M.data[static_cast<std::size_t>(y) * w];
                int* arow = &A[static_cast<std::size_t>(y) * w];
                for (int x = 0; x < w; ++x) {
                    running += drow[x];
                    const int c = (colacc[static_cast<std::size_t>(x)] += running);
                    if (c < min_count) continue;
                    const float snr = srow[x] * inv_sqrt[static_cast<std::size_t>(c)];
                    // ties go to the lower velocity index, independent of the thread schedule
                    if (snr > mrow[x] || (snr == mrow[x] && vi < arow[x])) {
                        mrow[x] = snr;
                        arow[x] = vi;
                    }
                }
            }
            const int d = ++done;
            if (tid == 0 && progress && (d % 8 == 0 || d == nv)) {
                if (!progress(static_cast<double>(d) / nv, "stack")) cancel = true;
            }
        }
    }
    if (cancel) throw Cancelled();

    StackResult out;
    out.max_snr = Image(w, h, neg);
    out.best_v.assign(static_cast<std::size_t>(w) * h, -1);
    for (std::size_t k = 0; k < lmax.size(); ++k) {
        if (lmax[k].empty()) continue;
        for (std::size_t q = 0; q < out.max_snr.data.size(); ++q) {
            const float v = lmax[k].data[q];
            const int a = larg[k][q];
            if (a < 0) continue;
            if (v > out.max_snr.data[q] || (v == out.max_snr.data[q] && (out.best_v[q] < 0 || a < out.best_v[q]))) {
                out.max_snr.data[q] = v;
                out.best_v[q] = a;
            }
        }
    }
    for (auto& v : out.max_snr.data)
        if (!std::isfinite(v)) v = 0.0f;
    if (progress) progress(1.0, "stack");
    return out;
}

}  // namespace fainttrail::detail
