// SPDX-License-Identifier: GPL-3.0-or-later
// Copyright 2026 Jurgen Kobierczynski
//
// Portable, reproducible random numbers (xoshiro256**). The standard library
// distributions are implementation-defined, so we provide our own to make
// simulations bit-identical across compilers and platforms.
#pragma once

#include <cmath>
#include <cstdint>

namespace fainttrail::detail {

class Rng {
public:
    explicit Rng(std::uint64_t seed) {
        std::uint64_t z = seed;
        for (auto& s : s_) s = splitmix(z);
    }

    std::uint64_t next() {
        const std::uint64_t result = rotl(s_[1] * 5, 7) * 9;
        const std::uint64_t t = s_[1] << 17;
        s_[2] ^= s_[0];
        s_[3] ^= s_[1];
        s_[1] ^= s_[2];
        s_[0] ^= s_[3];
        s_[2] ^= t;
        s_[3] = rotl(s_[3], 45);
        return result;
    }

    // uniform in [0,1)
    double uniform() { return static_cast<double>(next() >> 11) * (1.0 / 9007199254740992.0); }
    double uniform(double a, double b) { return a + (b - a) * uniform(); }
    int uniform_int(int lo, int hi_inclusive) {
        return lo + static_cast<int>(uniform() * (hi_inclusive - lo + 1));
    }

    double normal() {
        if (have_spare_) {
            have_spare_ = false;
            return spare_;
        }
        double u, v, s;
        do {
            u = uniform() * 2.0 - 1.0;
            v = uniform() * 2.0 - 1.0;
            s = u * u + v * v;
        } while (s >= 1.0 || s == 0.0);
        const double m = std::sqrt(-2.0 * std::log(s) / s);
        spare_ = v * m;
        have_spare_ = true;
        return u * m;
    }

    // Poisson deviate; exact for small means, normal approximation above 64.
    double poisson(double mean) {
        if (mean <= 0.0) return 0.0;
        if (mean > 64.0) {
            const double x = mean + std::sqrt(mean) * normal();
            return x < 0.0 ? 0.0 : std::round(x);
        }
        const double l = std::exp(-mean);
        double p = 1.0;
        int k = 0;
        do {
            ++k;
            p *= uniform();
        } while (p > l);
        return static_cast<double>(k - 1);
    }

private:
    static std::uint64_t rotl(std::uint64_t x, int k) { return (x << k) | (x >> (64 - k)); }
    static std::uint64_t splitmix(std::uint64_t& x) {
        std::uint64_t z = (x += 0x9e3779b97f4a7c15ULL);
        z = (z ^ (z >> 30)) * 0xbf58476d1ce4e5b9ULL;
        z = (z ^ (z >> 27)) * 0x94d049bb133111ebULL;
        return z ^ (z >> 31);
    }
    std::uint64_t s_[4]{};
    bool have_spare_ = false;
    double spare_ = 0.0;
};

}  // namespace fainttrail::detail
