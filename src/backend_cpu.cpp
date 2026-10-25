// SPDX-License-Identifier: GPL-3.0-or-later
// Copyright 2026 Jurgen Kobierczynski
//
// CPU reference implementation of the stacking contract (see backend.hpp).
//
// The kernel is limited by memory traffic, not arithmetic (docs/PROFILING.md).
// A straightforward implementation adds all frames into a full-size buffer for
// one motion at a time and then reads that buffer again: about 12 bytes of
// traffic per addition. This implementation is cache-blocked instead:
//
//   * the output is cut into bands of a few rows; one thread owns a band and
//     writes its rows of the result directly (no merge step);
//   * within a band the motions are processed in blocks of V; for each frame
//     the shifted rows are added into V small accumulators, which stay in L1/L2
//     while the frame's band is reused for all V motions;
//   * the coverage count of a row comes from the frames' column intervals: the
//     row splits into segments of constant coverage, so no count image is
//     needed and the inner loop is a plain, vectorisable maximum.
//
// Per pixel and motion the frames are still summed in ascending order, and the
// motions are visited in ascending order with a strict comparison, so the result
// is bit-identical to the plain definition (run_pass_naive) whatever the block
// sizes or the number of threads.
#include "backend_cpu.hpp"

#include <algorithm>
#include <atomic>
#include <chrono>
#include <cmath>
#include <limits>
#include <string>

#include "common.hpp"

namespace fainttrail::detail {

namespace {

struct Shift {
    int dx, dy;
    int x0, x1, y0, y1;  // output rectangle covered by the frame: [x0,x1) x [y0,y1)
};

}  // namespace

CpuBlocking choose_blocking(int w, int h, int nv, int nt, const CpuBlocking& req) {
    CpuBlocking b;
    // Accumulators of about 600 kB stay in L2 together with the frame rows being
    // read; 8 rows x 8 motions was the best setting on a 2328 px wide image
    // (tools/kernel_bench.cpp, docs/PROFILING.md).
    const double budget = 640.0 * 1024;
    b.band_rows = req.band_rows > 0 ? req.band_rows : 8;
    // enough bands to keep every thread busy with dynamic scheduling
    if (req.band_rows <= 0)
        while (b.band_rows > 1 && (h + b.band_rows - 1) / b.band_rows < 4 * nt) b.band_rows = std::max(1, b.band_rows / 2);
    b.band_rows = std::clamp(b.band_rows, 1, std::max(1, h));
    if (req.vel_block > 0) {
        b.vel_block = req.vel_block;
    } else {
        const double per_motion = 4.0 * b.band_rows * std::max(1, w);
        b.vel_block = static_cast<int>(std::clamp(budget / per_motion, 1.0, 16.0));
    }
    b.vel_block = std::clamp(b.vel_block, 1, std::max(1, nv));
    return b;
}

std::string CpuBackend::description() const {
#ifdef FAINTTRAIL_HAVE_OPENMP
    return "CPU reference, cache-blocked, OpenMP, " + std::to_string(thread_count(0)) + " threads";
#else
    return "CPU reference, cache-blocked, single thread (built without OpenMP)";
#endif
}

StackResult CpuBackend::run(const StackPass& pass, const ProgressFn& progress) {
    validate_pass(pass);
    using clk = std::chrono::steady_clock;
    const auto t_start = clk::now();
    const auto& frames = *pass.frames;
    const int w = frames[0].width, h = frames[0].height;
    const int n = static_cast<int>(frames.size());
    const int nv = pass_motions(pass);
    const int min_count = pass.min_count;
    const int nt = thread_count(pass.threads);
    const CpuBlocking blk = choose_blocking(w, h, nv, nt, blocking_);
    const int band = blk.band_rows, vblock = blk.vel_block;
    const int nbands = (h + band - 1) / band;

    // shifts and covered rectangles, per motion and frame (motion-major)
    std::vector<Shift> sh(static_cast<std::size_t>(nv) * n);
    for (int v = 0; v < nv; ++v)
        for (int i = 0; i < n; ++i) {
            Shift& s = sh[static_cast<std::size_t>(v) * n + i];
            pass_shift(pass, v, i, s.dx, s.dy);
            s.x0 = std::max(0, -s.dx);
            s.x1 = std::min(w, w - s.dx);
            s.y0 = std::max(0, -s.dy);
            s.y1 = std::min(h, h - s.dy);
            if (s.x0 >= s.x1 || s.y0 >= s.y1) s.x0 = s.x1 = s.y0 = s.y1 = 0;  // frame shifted out entirely
        }
    std::vector<float> inv_sqrt(static_cast<std::size_t>(n) + 1, 0.0f);
    for (int c = 1; c <= n; ++c) inv_sqrt[static_cast<std::size_t>(c)] = static_cast<float>(1.0 / std::sqrt(static_cast<double>(c)));

    const float neg = -std::numeric_limits<float>::infinity();
    StackResult out;
    out.max_snr = Image(w, h, neg);
    out.best_v.assign(static_cast<std::size_t>(w) * h, -1);
    std::vector<double> t_acc(static_cast<std::size_t>(nt), 0.0), t_max(static_cast<std::size_t>(nt), 0.0);
    std::atomic<int> done{0};
    std::atomic<bool> cancel{false};
    const auto t_par = clk::now();

#pragma omp parallel num_threads(nt)
    {
        const int tid = thread_id();
        std::vector<float> acc(static_cast<std::size_t>(vblock) * band * w);
        std::vector<std::pair<int, int>> edges(static_cast<std::size_t>(2 * n));
        double ta_sum = 0.0, tm_sum = 0.0;

#pragma omp for schedule(dynamic, 1)
        for (int b = 0; b < nbands; ++b) {
            if (cancel.load(std::memory_order_relaxed)) continue;
            const int ya = b * band, yb = std::min(h, ya + band), bh = yb - ya;
            for (int v0 = 0; v0 < nv; v0 += vblock) {
                const int vn = std::min(vblock, nv - v0);
                const auto t0 = clk::now();
                std::fill(acc.begin(), acc.begin() + static_cast<std::ptrdiff_t>(vn) * bh * w, 0.0f);
                for (int i = 0; i < n; ++i) {  // ascending frame order, as the contract requires
                    const float* g = frames[static_cast<std::size_t>(i)].data.data();
                    for (int k = 0; k < vn; ++k) {
                        const Shift& s = sh[static_cast<std::size_t>(v0 + k) * n + i];
                        const int y0 = std::max(ya, s.y0), y1 = std::min(yb, s.y1);
                        for (int y = y0; y < y1; ++y) {
                            float* o = &acc[(static_cast<std::size_t>(k) * bh + (y - ya)) * w];
                            const float* src = g + static_cast<std::ptrdiff_t>(y + s.dy) * w + s.dx;
                            for (int x = s.x0; x < s.x1; ++x) o[x] += src[x];
                        }
                    }
                }
                const auto t1 = clk::now();
                // coverage per row from the column intervals, then the running maximum;
                // motions ascend, so the strict comparison keeps the lowest index on ties
                for (int k = 0; k < vn; ++k) {
                    const int vi = v0 + k;
                    const Shift* S = &sh[static_cast<std::size_t>(vi) * n];
                    for (int y = ya; y < yb; ++y) {
                        // column intervals of the frames covering this row -> segments of
                        // constant coverage; only the edges are touched, not every column
                        int nedge = 0;
                        for (int i = 0; i < n; ++i)
                            if (y >= S[i].y0 && y < S[i].y1) {
                                edges[static_cast<std::size_t>(nedge++)] = {S[i].x0, +1};
                                edges[static_cast<std::size_t>(nedge++)] = {S[i].x1, -1};
                            }
                        if (nedge / 2 < min_count) continue;
                        std::sort(edges.begin(), edges.begin() + nedge);
                        const float* arow = &acc[(static_cast<std::size_t>(k) * bh + (y - ya)) * w];
                        float* mrow = &out.max_snr.data[static_cast<std::size_t>(y) * w];
                        int* brow = &out.best_v[static_cast<std::size_t>(y) * w];
                        int c = 0;
                        for (int e = 0; e < nedge;) {
                            const int xa = edges[static_cast<std::size_t>(e)].first;
                            while (e < nedge && edges[static_cast<std::size_t>(e)].first == xa) c += edges[static_cast<std::size_t>(e++)].second;
                            const int xb = e < nedge ? edges[static_cast<std::size_t>(e)].first : w;
                            if (c < min_count || xa >= xb) continue;
                            const float f = inv_sqrt[static_cast<std::size_t>(c)];
                            for (int x = xa; x < xb; ++x) {
                                const float snr = arow[x] * f;
                                const bool better = snr > mrow[x];
                                mrow[x] = better ? snr : mrow[x];
                                brow[x] = better ? vi : brow[x];
                            }
                        }
                    }
                }
                const auto t2 = clk::now();
                ta_sum += std::chrono::duration<double>(t1 - t0).count();
                tm_sum += std::chrono::duration<double>(t2 - t1).count();
            }
            const int d = ++done;
            if (tid == 0 && progress && !progress(static_cast<double>(d) / nbands, "stack")) cancel = true;
        }
        t_acc[static_cast<std::size_t>(tid)] = ta_sum;
        t_max[static_cast<std::size_t>(tid)] = tm_sum;
    }
    if (cancel) throw Cancelled();

    const auto t_fin = clk::now();
    for (auto& v : out.max_snr.data)
        if (!std::isfinite(v)) v = 0.0f;
    double sa = 0.0, sm = 0.0;
    for (int k = 0; k < nt; ++k) {
        sa += t_acc[static_cast<std::size_t>(k)];
        sm += t_max[static_cast<std::size_t>(k)];
    }
    out.timings.push_back({"accumulate", sa / nt});
    out.timings.push_back({"coverage + max", sm / nt});
    out.timings.push_back({"setup + finish", std::chrono::duration<double>(t_par - t_start).count() +
                                                 std::chrono::duration<double>(clk::now() - t_fin).count()});
    if (progress) progress(1.0, "stack");
    return out;
}

}  // namespace fainttrail::detail
