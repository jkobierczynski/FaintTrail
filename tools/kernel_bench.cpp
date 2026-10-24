// SPDX-License-Identifier: GPL-3.0-or-later
// Copyright 2026 Jurgen Kobierczynski
//
// Micro-benchmark of the stacking kernel, for profiling and for choosing a
// tiling strategy (see docs/PROFILING.md).
//
//   fainttrail_kernel_bench [--width 2328] [--height 1760] [--frames 20]
//                           [--velocities 48] [--threads N]
//
// It measures
//   1. memory bandwidth of this machine (stream triad and read-only sum),
//   2. the simple full-image kernel (run_pass_simple, FaintTrail <= 0.4),
//   3. the cache-blocked CPU backend with its automatic block sizes and with a
//      range of fixed ones, to tune choose_blocking() for other machines.
// Every result is checked bit for bit against the simple kernel.
#include <algorithm>
#include <chrono>
#include <cmath>
#include <cstdio>
#include <cstdlib>
#include <cstring>
#include <limits>
#include <string>
#include <vector>

#include "fainttrail/backend.hpp"
#include "../src/backend_cpu.hpp"

#ifdef _OPENMP
#include <omp.h>
#endif

using namespace fainttrail;
using clk = std::chrono::steady_clock;

namespace {

double secs(clk::time_point t0) { return std::chrono::duration<double>(clk::now() - t0).count(); }

int max_threads() {
#ifdef _OPENMP
    return omp_get_max_threads();
#else
    return 1;
#endif
}

// ---------------------------------------------------------------- bandwidth
double triad_gbs(int threads) {
    const std::size_t n = 32u << 20;  // 32 M floats = 128 MB per array
    std::vector<float> a(n), b(n, 1.0f), c(n, 2.0f);
    double best = 0.0;
    for (int rep = 0; rep < 4; ++rep) {
        const auto t0 = clk::now();
#pragma omp parallel for num_threads(threads) schedule(static)
        for (long long i = 0; i < static_cast<long long>(n); ++i) a[static_cast<std::size_t>(i)] = b[static_cast<std::size_t>(i)] + 0.5f * c[static_cast<std::size_t>(i)];
        const double s = secs(t0);
        best = std::max(best, 3.0 * 4.0 * n / s / 1e9);  // 2 reads + 1 write
    }
    volatile float sink = a[n / 3];
    (void)sink;
    return best;
}

double read_gbs(int threads) {
    const std::size_t n = 64u << 20;  // 256 MB
    std::vector<float> b(n, 1.0f);
    double best = 0.0;
    for (int rep = 0; rep < 4; ++rep) {
        double sum = 0.0;
        const auto t0 = clk::now();
#pragma omp parallel for num_threads(threads) schedule(static) reduction(+ : sum)
        for (long long i = 0; i < static_cast<long long>(n); ++i) sum += b[static_cast<std::size_t>(i)];
        const double s = secs(t0);
        best = std::max(best, 4.0 * n / s / 1e9);
        if (sum < 0) std::puts("");
    }
    return best;
}

}  // namespace

int main(int argc, char** argv) {
    int w = 2328, h = 1760, n = 20, nv = 48, threads = max_threads();
    for (int i = 1; i + 1 < argc; i += 2) {
        const std::string k = argv[i];
        const int v = std::atoi(argv[i + 1]);
        if (k == "--width") w = v;
        else if (k == "--height") h = v;
        else if (k == "--frames") n = v;
        else if (k == "--velocities") nv = v;
        else if (k == "--threads") threads = v;
    }
    std::printf("machine: %d threads\n", threads);
    std::printf("bandwidth: triad %.1f GB/s (1 thread), %.1f GB/s (%d threads); read %.1f / %.1f GB/s\n",
                triad_gbs(1), triad_gbs(threads), threads, read_gbs(1), read_gbs(threads));

    // random unit-variance frames over 1.6 h, motions spread like a real grid
    std::vector<Image> frames(static_cast<std::size_t>(n), Image(w, h));
    unsigned long long s = 12345;
    for (auto& f : frames)
        for (auto& v : f.data) {
            s ^= s << 13;
            s ^= s >> 7;
            s ^= s << 17;
            v = static_cast<float>((s >> 11) * (1.0 / 9007199254740992.0)) - 0.5f;
        }
    std::vector<double> t(static_cast<std::size_t>(n));
    for (int i = 0; i < n; ++i) t[static_cast<std::size_t>(i)] = -0.8 + 1.6 * i / std::max(1, n - 1);
    std::vector<std::pair<double, double>> vel;
    for (int k = 0; k < nv; ++k) {
        const double a = 2.0 * 3.14159265358979 * k / nv, r = 5.0 + 40.0 * (k % 7) / 6.0;
        vel.emplace_back(r * std::cos(a), r * std::sin(a));
    }
    StackPass pass;
    pass.frames = &frames;
    pass.t_h = &t;
    pass.velocities = &vel;
    pass.min_count = std::max(3, n / 2);
    pass.threads = threads;
    const double adds = static_cast<double>(nv) * n * w * h;
    std::printf("problem: %d x %d px, %d frames, %d motions = %.2f G pixel-adds per run, frames %.0f MB\n\n", w, h, n,
                nv, adds / 1e9, 4.0 * w * h * n / 1e6);

    pass.threads = threads;
    run_pass_simple(pass);  // warm-up
    auto t0 = clk::now();
    const StackResult ref = run_pass_simple(pass);
    double sec = secs(t0);
    std::printf("%-30s %2d thr %7.2f s  %6.2f G adds/s  (%.1f GB/s at 12 B/add)\n", "simple (one motion, whole image)",
                threads, sec, adds / sec / 1e9, 12.0 * adds / sec / 1e9);

    const detail::CpuBlocking automatic = detail::choose_blocking(w, h, nv, threads);
    std::vector<detail::CpuBlocking> cfgs = {automatic, {4, 16}, {8, 8}, {8, 16}, {12, 8}, {16, 4}, {16, 8}, {32, 4}, {64, 1}};
    for (std::size_t c = 0; c < cfgs.size(); ++c) {
        for (int th : {1, threads}) {
            if (th == 1 && c > 0) continue;
            pass.threads = th;
            detail::CpuBackend be(cfgs[c]);
            be.run(pass);
            t0 = clk::now();
            const StackResult r = be.run(pass);
            sec = secs(t0);
            const BackendCheck chk = check_backend(pass, ref, r);
            const bool same = chk.ok && r.best_v == ref.best_v && r.max_snr.data == ref.max_snr.data;
            char name[64];
            std::snprintf(name, sizeof name, "%s%d rows x %d motions", c == 0 ? "cpu (auto) " : "cpu ", cfgs[c].band_rows,
                          cfgs[c].vel_block);
            std::printf("%-30s %2d thr %7.2f s  %6.2f G adds/s  %s\n", name, th, sec, adds / sec / 1e9,
                        same ? "identical" : chk.message.c_str());
        }
    }
    return 0;
}
