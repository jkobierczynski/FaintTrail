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
//   2. the production CPU backend on random frames of the given size,
//   3. prototype cache-blocked variants of the same computation:
//      "band":       rows in bands of B rows, all velocities per band, so the
//                    accumulation buffer stays in L1/L2;
//      "band+vblock": additionally V velocities at a time per frame, so each
//                    frame band is read from memory once for V velocities.
// All variants are checked against the production backend with check_backend()
// (they sum in the same frame order, so they agree bit for bit).
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

// ---------------------------------------------------------------- variants
struct Shift {
    int dx, dy, x0, x1, y0, y1;  // valid output rectangle [x0,x1) x [y0,y1)
};

std::vector<Shift> shifts_for(const StackPass& p, std::size_t vi) {
    const int w = (*p.frames)[0].width, h = (*p.frames)[0].height;
    std::vector<Shift> s;
    for (std::size_t i = 0; i < p.frames->size(); ++i) {
        const double t = (*p.t_h)[i];
        Shift q;
        q.dx = static_cast<int>(std::lround((*p.velocities)[vi].first * t));
        q.dy = static_cast<int>(std::lround((*p.velocities)[vi].second * t));
        q.x0 = std::max(0, -q.dx);
        q.x1 = std::min(w, w - q.dx);
        q.y0 = std::max(0, -q.dy);
        q.y1 = std::min(h, h - q.dy);
        s.push_back(q);
    }
    return s;
}

// Rows in bands of `band` rows; velocities processed `vblock` at a time; within
// a velocity block the frames are the outer loop, so each frame band is reused
// for all velocities of the block while it is in cache.
StackResult run_blocked(const StackPass& p, int band, int vblock, int threads) {
    const auto& frames = *p.frames;
    const int w = frames[0].width, h = frames[0].height;
    const int n = static_cast<int>(frames.size());
    const int nv = static_cast<int>(p.velocities->size());
    std::vector<std::vector<Shift>> sh(static_cast<std::size_t>(nv));
    for (int v = 0; v < nv; ++v) sh[static_cast<std::size_t>(v)] = shifts_for(p, static_cast<std::size_t>(v));
    std::vector<float> inv(static_cast<std::size_t>(n) + 1, 0.0f);
    for (int c = 1; c <= n; ++c) inv[static_cast<std::size_t>(c)] = static_cast<float>(1.0 / std::sqrt(static_cast<double>(c)));

    StackResult out;
    out.max_snr = Image(w, h, 0.0f);
    out.best_v.assign(static_cast<std::size_t>(w) * h, -1);
    const int nbands = (h + band - 1) / band;
    const float neg = -std::numeric_limits<float>::infinity();
#pragma omp parallel num_threads(threads)
    {
        std::vector<float> acc(static_cast<std::size_t>(vblock) * band * w);
        std::vector<int> cov(static_cast<std::size_t>(w) + 1);
        std::vector<float> bmax(static_cast<std::size_t>(band) * w);
        std::vector<int> barg(static_cast<std::size_t>(band) * w);
#pragma omp for schedule(dynamic, 1)
        for (int b = 0; b < nbands; ++b) {
            const int ya = b * band, yb = std::min(h, ya + band), bh = yb - ya;
            std::fill(bmax.begin(), bmax.end(), neg);
            std::fill(barg.begin(), barg.end(), -1);
            for (int v0 = 0; v0 < nv; v0 += vblock) {
                const int vn = std::min(vblock, nv - v0);
                std::fill(acc.begin(), acc.begin() + static_cast<std::ptrdiff_t>(vn) * bh * w, 0.0f);
                for (int i = 0; i < n; ++i) {
                    const float* g = frames[static_cast<std::size_t>(i)].data.data();
                    for (int k = 0; k < vn; ++k) {
                        const Shift& s = sh[static_cast<std::size_t>(v0 + k)][static_cast<std::size_t>(i)];
                        const int y0 = std::max(ya, s.y0), y1 = std::min(yb, s.y1);
                        if (s.x0 >= s.x1) continue;
                        for (int y = y0; y < y1; ++y) {
                            float* o = &acc[(static_cast<std::size_t>(k) * bh + (y - ya)) * w];
                            const float* src = g + static_cast<std::ptrdiff_t>(y + s.dy) * w + s.dx;
                            for (int x = s.x0; x < s.x1; ++x) o[x] += src[x];
                        }
                    }
                }
                // coverage per row from the frames' x intervals, then the running max
                for (int k = 0; k < vn; ++k) {
                    const int vi = v0 + k;
                    const auto& S = sh[static_cast<std::size_t>(vi)];
                    for (int y = ya; y < yb; ++y) {
                        std::fill(cov.begin(), cov.end(), 0);
                        for (const auto& s : S)
                            if (y >= s.y0 && y < s.y1 && s.x0 < s.x1) {
                                ++cov[static_cast<std::size_t>(s.x0)];
                                --cov[static_cast<std::size_t>(s.x1)];
                            }
                        const float* arow = &acc[(static_cast<std::size_t>(k) * bh + (y - ya)) * w];
                        float* mrow = &bmax[static_cast<std::size_t>(y - ya) * w];
                        int* brow = &barg[static_cast<std::size_t>(y - ya) * w];
                        int c = 0;
                        for (int x = 0; x < w; ++x) {
                            c += cov[static_cast<std::size_t>(x)];
                            if (c < p.min_count) continue;
                            const float snr = arow[x] * inv[static_cast<std::size_t>(c)];
                            if (snr > mrow[x]) {  // velocities ascend within a band: ties keep the lower index
                                mrow[x] = snr;
                                brow[x] = vi;
                            }
                        }
                    }
                }
            }
            for (int y = ya; y < yb; ++y)
                for (int x = 0; x < w; ++x) {
                    const std::size_t q = static_cast<std::size_t>(y) * w + x;
                    const int a = barg[static_cast<std::size_t>(y - ya) * w + x];
                    if (a >= 0) {
                        out.max_snr.data[q] = bmax[static_cast<std::size_t>(y - ya) * w + x];
                        out.best_v[q] = a;
                    }
                }
        }
    }
    return out;
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

    auto cpu = make_backend("cpu");
    for (int th : {1, threads}) {
        pass.threads = th;
        cpu->run(pass);  // warm-up
        const auto t0 = clk::now();
        const StackResult r = cpu->run(pass);
        const double sec = secs(t0);
        std::printf("%-26s %2d thr %7.2f s  %6.2f G adds/s  (%.1f GB/s at 12 B/add)", "production (per velocity)", th,
                    sec, adds / sec / 1e9, 12.0 * adds / sec / 1e9);
        for (const auto& st : r.timings) std::printf("  %s %.2f", st.name.c_str(), st.seconds);
        std::printf("\n");
    }
    pass.threads = threads;
    const StackResult ref = cpu->run(pass);

    struct Cfg {
        int band, vblock;
    };
    for (const Cfg c : {Cfg{8, 1}, Cfg{32, 1}, Cfg{128, 1}, Cfg{16, 4}, Cfg{16, 8}, Cfg{8, 16}, Cfg{32, 8}}) {
        for (int th : {1, threads}) {
            run_blocked(pass, c.band, c.vblock, th);
            const auto t0 = clk::now();
            const StackResult r = run_blocked(pass, c.band, c.vblock, th);
            const double sec = secs(t0);
            const BackendCheck chk = check_backend(pass, ref, r);
            char name[64];
            std::snprintf(name, sizeof name, "band %d rows, %d vel/block", c.band, c.vblock);
            std::printf("%-26s %2d thr %7.2f s  %6.2f G adds/s  %s\n", name, th, sec, adds / sec / 1e9,
                        chk.ok && chk.max_abs_diff == 0.0 ? "identical" : chk.message.c_str());
        }
    }
    return 0;
}
