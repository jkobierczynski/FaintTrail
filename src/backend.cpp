// SPDX-License-Identifier: GPL-3.0-or-later
// Copyright 2026 Jurgen Kobierczynski
//
// Backend registry, the naive reference implementation of the stacking
// contract, and the conformance check.
#include "fainttrail/backend.hpp"

#include <algorithm>
#include <atomic>
#include <cmath>
#include <limits>
#include <stdexcept>

#include "backend_cpu.hpp"
#include "common.hpp"

namespace fainttrail {

namespace {
std::vector<std::pair<std::string, BackendFactory>>& registry() {
    static std::vector<std::pair<std::string, BackendFactory>> r;
    return r;
}
}  // namespace

void register_backend(const std::string& name, BackendFactory factory) {
    if (name.empty() || name == "auto" || name == "cpu") throw std::runtime_error("reserved backend name '" + name + "'");
    for (auto& e : registry())
        if (e.first == name) {
            e.second = std::move(factory);
            return;
        }
    registry().emplace_back(name, std::move(factory));
}

std::vector<std::string> backend_names() {
    std::vector<std::string> names;
    // GPU backends will be listed here first, when compiled in and a device is present.
    names.push_back("cpu");
    for (const auto& e : registry()) names.push_back(e.first);
    return names;
}

std::unique_ptr<SearchBackend> make_backend(const std::string& name) {
    if (name.empty() || name == "auto" || name == "cpu") return std::make_unique<detail::CpuBackend>();
    for (const auto& e : registry())
        if (e.first == name) return e.second();
    std::string known;
    for (const auto& n : backend_names()) known += (known.empty() ? "" : ", ") + n;
    throw std::runtime_error("unknown or unavailable search backend '" + name + "' (available: " + known + ", auto)");
}

void validate_pass(const StackPass& pass) {
    if (!pass.frames || !pass.t_h || !pass.velocities) throw std::runtime_error("stack pass: missing input");
    const auto& f = *pass.frames;
    if (f.empty()) throw std::runtime_error("stack pass: no frames");
    if (f.size() != pass.t_h->size()) throw std::runtime_error("stack pass: frame and time counts differ");
    for (const auto& im : f)
        if (im.width != f[0].width || im.height != f[0].height || im.size() != static_cast<std::size_t>(im.width) * im.height)
            throw std::runtime_error("stack pass: frames differ in size");
    if (pass.min_count < 1) throw std::runtime_error("stack pass: min_count must be at least 1");
}

namespace {

inline int shift_of(double v, double t) { return static_cast<int>(std::lround(v * t)); }

// S_v at one pixel, exactly as the contract defines it. Returns false if undefined.
bool snr_at(const StackPass& pass, std::size_t vi, int x, int y, float& out) {
    const auto& frames = *pass.frames;
    const auto& t = *pass.t_h;
    const auto [vx, vy] = (*pass.velocities)[vi];
    const int w = frames[0].width, h = frames[0].height;
    float sum = 0.0f;
    int c = 0;
    for (std::size_t i = 0; i < frames.size(); ++i) {
        const int sx = x + shift_of(vx, t[i]), sy = y + shift_of(vy, t[i]);
        if (sx < 0 || sy < 0 || sx >= w || sy >= h) continue;
        sum += frames[i].data[static_cast<std::size_t>(sy) * w + sx];
        ++c;
    }
    if (c < pass.min_count) return false;
    out = sum * static_cast<float>(1.0 / std::sqrt(static_cast<double>(c)));
    return true;
}

}  // namespace

StackResult run_pass_naive(const StackPass& pass) {
    validate_pass(pass);
    const int w = (*pass.frames)[0].width, h = (*pass.frames)[0].height;
    StackResult r;
    r.max_snr = Image(w, h, 0.0f);
    r.best_v.assign(static_cast<std::size_t>(w) * h, -1);
    for (int y = 0; y < h; ++y)
        for (int x = 0; x < w; ++x) {
            float best = -std::numeric_limits<float>::infinity();
            int arg = -1;
            for (std::size_t vi = 0; vi < pass.velocities->size(); ++vi) {
                float s;
                if (!snr_at(pass, vi, x, y, s)) continue;
                if (s > best) {  // strict: ties keep the lower index
                    best = s;
                    arg = static_cast<int>(vi);
                }
            }
            if (arg >= 0) {
                r.max_snr.at(x, y) = best;
                r.best_v[static_cast<std::size_t>(y) * w + x] = arg;
            }
        }
    return r;
}

// Straightforward full-image implementation (FaintTrail 0.2-0.4's CPU kernel):
// one motion at a time, all frames added into a full-size buffer, coverage from
// a 2-D difference array. Independent of the cache-blocked backend and fast
// enough to serve as its reference on realistic image sizes.
StackResult run_pass_simple(const StackPass& pass) {
    const ProgressFn progress;
    validate_pass(pass);
    const auto& frames = *pass.frames;
    const auto& t = *pass.t_h;
    const auto& grid = *pass.velocities;
    const int w = frames[0].width, h = frames[0].height;
    const int n = static_cast<int>(frames.size());
    const int nv = static_cast<int>(grid.size());
    const int min_count = pass.min_count;
    const int nt = detail::thread_count(pass.threads);

    std::vector<float> inv_sqrt(static_cast<std::size_t>(n) + 1, 0.0f);
    for (int c = 1; c <= n; ++c) inv_sqrt[static_cast<std::size_t>(c)] = static_cast<float>(1.0 / std::sqrt(static_cast<double>(c)));

    const float neg = -std::numeric_limits<float>::infinity();
    std::vector<Image> lmax(static_cast<std::size_t>(nt));
    std::vector<std::vector<int>> larg(static_cast<std::size_t>(nt));
    std::atomic<int> done{0};
    std::atomic<bool> cancel{false};

#pragma omp parallel num_threads(nt)
    {
        const int tid = detail::thread_id();
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
    for (int k = 0; k < nt; ++k) {
    }
    if (progress) progress(1.0, "stack");
    return out;
}


BackendCheck check_backend(const StackPass& pass, const StackResult& ref, const StackResult& cand) {
    BackendCheck c;
    validate_pass(pass);
    const int w = (*pass.frames)[0].width, h = (*pass.frames)[0].height;
    if (ref.max_snr.width != w || ref.max_snr.height != h || cand.max_snr.width != w || cand.max_snr.height != h ||
        ref.best_v.size() != ref.max_snr.size() || cand.best_v.size() != cand.max_snr.size()) {
        c.message = "result size does not match the pass";
        return c;
    }
    c.pixels = ref.max_snr.size();
    std::size_t bad_value = 0;
    for (int y = 0; y < h; ++y)
        for (int x = 0; x < w; ++x) {
            const std::size_t k = static_cast<std::size_t>(y) * w + x;
            const double a = ref.max_snr.data[k], b = cand.max_snr.data[k];
            const double tol = 1e-4 * std::max(1.0, std::fabs(a));
            const double d = std::fabs(a - b);
            if (!(d <= tol)) ++bad_value;  // also catches NaN
            if (std::isfinite(d)) c.max_abs_diff = std::max(c.max_abs_diff, d);
            const int va = ref.best_v[k], vb = cand.best_v[k];
            if (va == vb) continue;
            ++c.best_v_mismatch;
            // allowed only if both velocities are defined here and give the same S/N within tolerance
            float sa = 0, sb = 0;
            const bool da = va >= 0 && snr_at(pass, static_cast<std::size_t>(va), x, y, sa);
            const bool db = vb >= 0 && static_cast<std::size_t>(vb) < pass.velocities->size() &&
                            snr_at(pass, static_cast<std::size_t>(vb), x, y, sb);
            if (!(da && db && std::fabs(sa - sb) <= tol)) ++c.best_v_unexplained;
        }
    c.ok = bad_value == 0 && c.best_v_unexplained == 0;
    c.message = c.ok ? "conforms" : std::to_string(bad_value) + " pixels outside the S/N tolerance, " +
                                        std::to_string(c.best_v_unexplained) + " unexplained velocity differences";
    return c;
}

}  // namespace fainttrail
