// SPDX-License-Identifier: GPL-3.0-or-later
// Copyright 2026 Jurgen Kobierczynski
#include <cmath>

#include "fainttrail/types.hpp"

namespace fainttrail {

namespace {
constexpr double kPi = 3.14159265358979323846;
constexpr double kDeg = kPi / 180.0;
}  // namespace

bool Wcs::pix2sky(double x, double y, double& ra_deg, double& dec_deg) const {
    if (!valid) return false;
    // FITS pixel coordinates are 1-based
    const double dx = (x + 1.0) - crpix1;
    const double dy = (y + 1.0) - crpix2;
    // intermediate world coordinates in radians on the tangent plane
    const double xi = (cd11 * dx + cd12 * dy) * kDeg;
    const double eta = (cd21 * dx + cd22 * dy) * kDeg;
    const double a0 = crval1 * kDeg, d0 = crval2 * kDeg;
    const double denom = std::cos(d0) - eta * std::sin(d0);
    double ra = a0 + std::atan2(xi, denom);
    const double dec = std::atan2(std::sin(d0) + eta * std::cos(d0),
                                  std::sqrt(xi * xi + denom * denom));
    ra = std::fmod(ra, 2.0 * kPi);
    if (ra < 0) ra += 2.0 * kPi;
    ra_deg = ra / kDeg;
    dec_deg = dec / kDeg;
    return true;
}

bool Wcs::sky2pix(double ra_deg, double dec_deg, double& x, double& y) const {
    if (!valid) return false;
    const double a = ra_deg * kDeg, d = dec_deg * kDeg;
    const double a0 = crval1 * kDeg, d0 = crval2 * kDeg;
    const double cosc = std::sin(d0) * std::sin(d) + std::cos(d0) * std::cos(d) * std::cos(a - a0);
    if (cosc <= 0.0) return false;  // other hemisphere of the projection
    const double xi = std::cos(d) * std::sin(a - a0) / cosc / kDeg;
    const double eta = (std::cos(d0) * std::sin(d) - std::sin(d0) * std::cos(d) * std::cos(a - a0)) /
                       cosc / kDeg;
    const double det = cd11 * cd22 - cd12 * cd21;
    if (det == 0.0) return false;
    const double dx = (cd22 * xi - cd12 * eta) / det;
    const double dy = (-cd21 * xi + cd11 * eta) / det;
    x = dx + crpix1 - 1.0;
    y = dy + crpix2 - 1.0;
    return true;
}

double Wcs::pixel_scale_arcsec() const {
    if (!valid) return 0.0;
    return std::sqrt(std::fabs(cd11 * cd22 - cd12 * cd21)) * 3600.0;
}

Wcs Wcs::make_tan(double ra_deg, double dec_deg, double crpix1_, double crpix2_,
                  double scale_arcsec, double rotation_deg) {
    Wcs w;
    const double s = scale_arcsec / 3600.0;
    const double r = rotation_deg * kDeg;
    // north up, east left (RA increases to the left) for rotation 0
    w.cd11 = -s * std::cos(r);
    w.cd12 = s * std::sin(r);
    w.cd21 = s * std::sin(r);
    w.cd22 = s * std::cos(r);
    w.crval1 = ra_deg;
    w.crval2 = dec_deg;
    w.crpix1 = crpix1_;
    w.crpix2 = crpix2_;
    w.valid = true;
    return w;
}

}  // namespace fainttrail
