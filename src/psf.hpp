// SPDX-License-Identifier: GPL-3.0-or-later
// Copyright 2026 Jurgen Kobierczynski
//
// Pixel-integrated Gaussian PSF rendering, shared by the simulator and the
// injection-recovery tests.
#pragma once

#include <algorithm>
#include <cmath>
#include <vector>

#include "fainttrail/types.hpp"

namespace fainttrail::detail {

constexpr double kPi = 3.14159265358979323846;
constexpr double kFwhmToSigma = 1.0 / 2.354820045030949;

// Add `flux` (in image units) as a pixel-integrated Gaussian centred at (x0,y0).
inline void add_gaussian(Image& img, double x0, double y0, double flux, double sigma,
                         int max_radius = 30) {
    if (flux <= 0.0 || sigma <= 0.0) return;
    const double peak = flux / (2.0 * kPi * sigma * sigma);
    // radius where the profile drops below ~0.05 units
    int r = static_cast<int>(std::ceil(sigma * std::sqrt(2.0 * std::log(std::max(peak / 0.05, 2.0)))));
    r = std::clamp(r, 2, max_radius);
    const int xc = static_cast<int>(std::floor(x0 + 0.5));
    const int yc = static_cast<int>(std::floor(y0 + 0.5));
    const int x_lo = std::max(0, xc - r), x_hi = std::min(img.width - 1, xc + r);
    const int y_lo = std::max(0, yc - r), y_hi = std::min(img.height - 1, yc + r);
    if (x_lo > x_hi || y_lo > y_hi) return;
    const double inv = 1.0 / (sigma * std::sqrt(2.0));
    std::vector<double> wx(static_cast<std::size_t>(x_hi - x_lo + 1));
    std::vector<double> wy(static_cast<std::size_t>(y_hi - y_lo + 1));
    for (int x = x_lo; x <= x_hi; ++x)
        wx[static_cast<std::size_t>(x - x_lo)] =
            0.5 * (std::erf((x + 0.5 - x0) * inv) - std::erf((x - 0.5 - x0) * inv));
    for (int y = y_lo; y <= y_hi; ++y)
        wy[static_cast<std::size_t>(y - y_lo)] =
            0.5 * (std::erf((y + 0.5 - y0) * inv) - std::erf((y - 0.5 - y0) * inv));
    for (int y = y_lo; y <= y_hi; ++y) {
        const double fy = flux * wy[static_cast<std::size_t>(y - y_lo)];
        float* row = &img.data[static_cast<std::size_t>(y) * img.width];
        for (int x = x_lo; x <= x_hi; ++x) row[x] += static_cast<float>(fy * wx[static_cast<std::size_t>(x - x_lo)]);
    }
}

// Add a source that moves from (x0,y0) to (x1,y1) during the exposure.
inline void add_trailed(Image& img, double x0, double y0, double x1, double y1, double flux,
                        double sigma) {
    const double len = std::hypot(x1 - x0, y1 - y0);
    const int steps = std::max(1, static_cast<int>(std::ceil(len / 0.25)));
    const double f = flux / steps;
    for (int k = 0; k < steps; ++k) {
        const double u = (k + 0.5) / steps;
        add_gaussian(img, x0 + u * (x1 - x0), y0 + u * (y1 - y0), f, sigma);
    }
}

}  // namespace fainttrail::detail
