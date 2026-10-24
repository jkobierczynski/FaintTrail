// SPDX-License-Identifier: GPL-3.0-or-later
// Copyright 2026 Jurgen Kobierczynski
//
// Small numeric helpers shared by the implementation files.
#pragma once

#include <algorithm>
#include <cmath>
#include <vector>

#include "fainttrail/types.hpp"

#ifdef FAINTTRAIL_HAVE_OPENMP
#include <omp.h>
#endif

namespace fainttrail::detail {

inline int thread_count(int requested) {
#ifdef FAINTTRAIL_HAVE_OPENMP
    return requested > 0 ? requested : omp_get_max_threads();
#else
    (void)requested;
    return 1;
#endif
}

inline int thread_id() {
#ifdef FAINTTRAIL_HAVE_OPENMP
    return omp_get_thread_num();
#else
    return 0;
#endif
}

// Median of a vector (modified in place). Mean of the two middle values for even sizes.
inline double median_inplace(std::vector<float>& v) {
    if (v.empty()) return 0.0;
    const std::size_t n = v.size(), mid = n / 2;
    std::nth_element(v.begin(), v.begin() + static_cast<std::ptrdiff_t>(mid), v.end());
    double m = v[mid];
    if (n % 2 == 0) {
        const float lower = *std::max_element(v.begin(), v.begin() + static_cast<std::ptrdiff_t>(mid));
        m = 0.5 * (m + lower);
    }
    return m;
}

// Robust standard deviation via the median absolute deviation.
inline double robust_sigma(std::vector<float>& v) {
    if (v.size() < 3) return 0.0;
    const double med = median_inplace(v);
    for (auto& x : v) x = static_cast<float>(std::fabs(x - med));
    return 1.4826 * median_inplace(v);
}

// Regularly sampled finite pixel values (at most ~max_samples of them).
inline std::vector<float> sample_pixels(const Image& img, std::size_t max_samples,
                                        const Image* mask = nullptr, int border = 0) {
    std::vector<float> out;
    const std::size_t usable = static_cast<std::size_t>(std::max(0, img.width - 2 * border)) *
                               static_cast<std::size_t>(std::max(0, img.height - 2 * border));
    if (usable == 0) return out;
    const int step = std::max(1, static_cast<int>(std::sqrt(static_cast<double>(usable) / max_samples)));
    out.reserve(usable / (static_cast<std::size_t>(step) * step) + 16);
    for (int y = border; y < img.height - border; y += step)
        for (int x = border; x < img.width - border; x += step) {
            if (mask && mask->at(x, y) == 0.0f) continue;
            const float v = img.at(x, y);
            if (std::isfinite(v)) out.push_back(v);
        }
    return out;
}

// Quantile of |z| for a standard normal z (inverse of the half-normal CDF),
// q in [0,1). Acklam's rational approximation of the normal quantile.
inline double half_normal_quantile(double q) {
    const double p = 0.5 * (1.0 + std::clamp(q, 0.0, 0.999999));
    static const double a[] = {-3.969683028665376e+01, 2.209460984245205e+02, -2.759285104469687e+02,
                               1.383577518672690e+02, -3.066479806614716e+01, 2.506628277459239e+00};
    static const double b[] = {-5.447609879822406e+01, 1.615858368580409e+02, -1.556989798598866e+02,
                               6.680131188771972e+01, -1.328068155288572e+01};
    static const double c[] = {-7.784894002430293e-03, -3.223964580411365e-01, -2.400758277161838e+00,
                               -2.549732539343734e+00, 4.374664141464968e+00, 2.938163982698783e+00};
    static const double d[] = {7.784695709041462e-03, 3.224671290700398e-01, 2.445134137142996e+00,
                               3.754408661907416e+00};
    if (p > 0.97575) {
        const double r = std::sqrt(-2.0 * std::log(1.0 - p));
        return -(((((c[0] * r + c[1]) * r + c[2]) * r + c[3]) * r + c[4]) * r + c[5]) /
               ((((d[0] * r + d[1]) * r + d[2]) * r + d[3]) * r + 1.0);
    }
    const double r0 = p - 0.5, r = r0 * r0;
    return (((((a[0] * r + a[1]) * r + a[2]) * r + a[3]) * r + a[4]) * r + a[5]) * r0 /
           (((((b[0] * r + b[1]) * r + b[2]) * r + b[3]) * r + b[4]) * r + 1.0);
}

// Bilinear sample; outside pixels return NaN.
inline float bilinear(const Image& img, double x, double y) {
    const int x0 = static_cast<int>(std::floor(x)), y0 = static_cast<int>(std::floor(y));
    if (x0 < 0 || y0 < 0 || x0 + 1 >= img.width || y0 + 1 >= img.height) {
        // allow sampling exactly on the last row/column
        const int xr = static_cast<int>(std::lround(x)), yr = static_cast<int>(std::lround(y));
        if (img.inside(xr, yr) && std::fabs(x - xr) < 1e-9 && std::fabs(y - yr) < 1e-9) return img.at(xr, yr);
        return NAN;
    }
    const double fx = x - x0, fy = y - y0;
    const float* r0 = &img.data[static_cast<std::size_t>(y0) * img.width + x0];
    const float* r1 = r0 + img.width;
    return static_cast<float>((1 - fy) * ((1 - fx) * r0[0] + fx * r0[1]) + fy * ((1 - fx) * r1[0] + fx * r1[1]));
}

}  // namespace fainttrail::detail
