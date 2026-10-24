// SPDX-License-Identifier: GPL-3.0-or-later
// Copyright 2026 Jurgen Kobierczynski
//
// Backend registry, the naive reference implementation of the stacking
// contract, and the conformance check.
#include "fainttrail/backend.hpp"

#include <cmath>
#include <limits>
#include <stdexcept>

#include "backend_cpu.hpp"

namespace fainttrail {

std::vector<std::string> backend_names() {
    std::vector<std::string> names;
    // GPU backends will be listed here first, when compiled in and a device is present.
    names.push_back("cpu");
    return names;
}

std::unique_ptr<SearchBackend> make_backend(const std::string& name) {
    if (name.empty() || name == "auto" || name == "cpu") return std::make_unique<detail::CpuBackend>();
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
