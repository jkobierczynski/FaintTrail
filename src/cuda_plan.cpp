// SPDX-License-Identifier: GPL-3.0-or-later
// Copyright 2026 Jurgen Kobierczynski
#include "cuda_plan.hpp"

#include <algorithm>
#include <cmath>
#include <cstdlib>
#include <limits>
#include <sstream>
#include <stdexcept>

namespace fainttrail::detail {

CudaOptions parse_cuda_options(const std::string& spec) {
    CudaOptions o;
    std::stringstream ss(spec);
    std::string item;
    while (std::getline(ss, item, ',')) {
        if (item.empty()) continue;
        const auto eq = item.find('=');
        const std::string key = item.substr(0, eq), val = eq == std::string::npos ? "" : item.substr(eq + 1);
        char* end = nullptr;
        const double d = std::strtod(val.c_str(), &end);
        if (val.empty() || end == val.c_str() || *end != '\0')
            throw std::runtime_error("cuda backend option '" + item + "': expected key=number");
        if (key == "device") o.device = static_cast<int>(d);
        else if (key == "max_mb") o.max_mb = d;
        else if (key == "launch_adds") o.launch_adds = d;
        else throw std::runtime_error("unknown cuda backend option '" + key + "' (device, max_mb, launch_adds)");
    }
    if (o.max_mb < 0 || o.launch_adds <= 0) throw std::runtime_error("cuda backend options must be positive");
    return o;
}

CudaPlan make_cuda_plan(const StackPass& pass, double budget_bytes, double launch_adds) {
    validate_pass(pass);
    const auto& frames = *pass.frames;
    CudaPlan p;
    p.w = frames[0].width;
    p.h = frames[0].height;
    p.n = static_cast<int>(frames.size());
    p.nv = pass_motions(pass);
    p.min_count = pass.min_count;
    p.dx.resize(static_cast<std::size_t>(p.nv) * p.n);
    p.dy.resize(p.dx.size());
    p.min_dy = std::numeric_limits<int>::max();
    p.max_dy = std::numeric_limits<int>::min();
    for (int j = 0; j < p.nv; ++j)
        for (int i = 0; i < p.n; ++i) {
            const std::size_t s = static_cast<std::size_t>(j) * p.n + i;
            pass_shift(pass, j, i, p.dx[s], p.dy[s]);
            // shifts beyond the frame never cover anything; clamping them keeps
            // the band arithmetic small and does not change any result
            p.dy[s] = std::clamp(p.dy[s], -p.h, p.h);
            p.dx[s] = std::clamp(p.dx[s], -p.w, p.w);
            p.min_dy = std::min(p.min_dy, p.dy[s]);
            p.max_dy = std::max(p.max_dy, p.dy[s]);
        }
    p.inv_sqrt.assign(static_cast<std::size_t>(p.n) + 1, 0.0f);
    for (int c = 1; c <= p.n; ++c) p.inv_sqrt[static_cast<std::size_t>(c)] = static_cast<float>(1.0 / std::sqrt(static_cast<double>(c)));

    p.bands = cuda_core::plan_bands(p.w, p.h, p.n, p.min_dy, p.max_dy, budget_bytes);
    for (const auto& b : p.bands) p.frame_bytes = std::max(p.frame_bytes, 4.0 * p.w * p.n * (b.r1 - b.r0));

    int max_rows = 1;
    for (const auto& b : p.bands) max_rows = std::max(max_rows, b.y1 - b.y0);
    const double per_motion = static_cast<double>(p.w) * max_rows * p.n;
    const double m = std::clamp(launch_adds / std::max(1.0, per_motion), 1.0, static_cast<double>(p.nv));
    const int k = cuda_core::kMotionBlock;
    p.chunk = std::max(k, (static_cast<int>(m) + k - 1) / k * k);
    return p;
}

StackResult emulate_cuda_pass(const StackPass& pass, double budget_bytes, double launch_adds) {
    const CudaPlan p = make_cuda_plan(pass, budget_bytes, launch_adds);
    const auto& frames = *pass.frames;
    StackResult out;
    out.max_snr = Image(p.w, p.h, 0.0f);
    out.best_v.assign(static_cast<std::size_t>(p.w) * p.h, -1);
    std::vector<float> buf;
    for (const auto& b : p.bands) {
        const int rows = b.r1 - b.r0;
        buf.resize(static_cast<std::size_t>(p.n) * rows * p.w);  // the "upload"
        for (int i = 0; i < p.n; ++i)
            std::copy_n(frames[static_cast<std::size_t>(i)].data.begin() + static_cast<std::ptrdiff_t>(b.r0) * p.w,
                        static_cast<std::ptrdiff_t>(rows) * p.w, buf.begin() + static_cast<std::ptrdiff_t>(i) * rows * p.w);
        cuda_core::BandView v;
        v.frames = buf.data();
        v.w = p.w;
        v.h = p.h;
        v.r0 = b.r0;
        v.rows = rows;
        v.n = p.n;
        v.dx = p.dx.data();
        v.dy = p.dy.data();
        v.inv_sqrt = p.inv_sqrt.data();
        v.min_count = p.min_count;
        const std::size_t np = static_cast<std::size_t>(b.y1 - b.y0) * p.w;
        std::vector<float> best(np, -std::numeric_limits<float>::infinity());
        std::vector<int> arg(np, -1);
        for (int j0 = 0; j0 < p.nv; j0 += p.chunk) {  // the "launches"
            const int j1 = std::min(p.nv, j0 + p.chunk);
#pragma omp parallel for schedule(static)
            for (int y = b.y0; y < b.y1; ++y)
                for (int x = 0; x < p.w; ++x) {
                    const std::size_t o = static_cast<std::size_t>(y - b.y0) * p.w + x;
                    cuda_core::stack_pixel(v, x, y, j0, j1, best[o], arg[o]);
                }
        }
        for (std::size_t o = 0; o < np; ++o)  // the "download"
            if (arg[o] >= 0) {
                out.max_snr.data[static_cast<std::size_t>(b.y0) * p.w + o] = best[o];
                out.best_v[static_cast<std::size_t>(b.y0) * p.w + o] = arg[o];
            }
    }
    return out;
}

}  // namespace fainttrail::detail
