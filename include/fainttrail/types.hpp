// SPDX-License-Identifier: GPL-3.0-or-later
// Copyright 2026 Jurgen Kobierczynski
//
// Core data types shared by every part of fainttrail.
//
// Pixel convention for the C++ API: 0-based, the centre of the first pixel
// is (0,0). CSV output and the C API use the FITS convention (1-based).
#pragma once

#include <cstddef>
#include <cstdint>
#include <functional>
#include <stdexcept>
#include <string>
#include <vector>

namespace fainttrail {

// A single-plane float image, row-major, y increasing downwards in memory.
struct Image {
    int width = 0;
    int height = 0;
    std::vector<float> data;

    Image() = default;
    Image(int w, int h, float fill = 0.0f)
        : width(w), height(h), data(static_cast<std::size_t>(w) * h, fill) {}

    bool empty() const { return data.empty(); }
    std::size_t size() const { return data.size(); }
    float& at(int x, int y) { return data[static_cast<std::size_t>(y) * width + x]; }
    float at(int x, int y) const { return data[static_cast<std::size_t>(y) * width + x]; }
    bool inside(int x, int y) const { return x >= 0 && y >= 0 && x < width && y < height; }
};

// Gnomonic (TAN) world coordinate system, enough for registered amateur frames.
// crpix is stored 1-based exactly as in the FITS header.
struct Wcs {
    bool valid = false;
    double crval1 = 0.0, crval2 = 0.0;  // RA, Dec of reference point (deg)
    double crpix1 = 0.0, crpix2 = 0.0;  // reference pixel (1-based FITS)
    double cd11 = 0.0, cd12 = 0.0, cd21 = 0.0, cd22 = 0.0;  // deg / pixel

    // 0-based pixel -> RA/Dec in degrees. Returns false if the WCS is invalid.
    bool pix2sky(double x, double y, double& ra_deg, double& dec_deg) const;
    // RA/Dec in degrees -> 0-based pixel.
    bool sky2pix(double ra_deg, double dec_deg, double& x, double& y) const;
    // Mean pixel scale in arcseconds per pixel (0 if invalid).
    double pixel_scale_arcsec() const;

    static Wcs make_tan(double ra_deg, double dec_deg, double crpix1, double crpix2,
                        double scale_arcsec, double rotation_deg = 0.0);
};

// One registered exposure.
struct Frame {
    Image image;
    double mjd_mid = 0.0;    // UTC MJD at mid-exposure
    double exptime_s = 0.0;  // exposure time in seconds
    Wcs wcs;                 // optional, used for sky coordinates in the output
    // Optional: offset that maps this registered frame back onto the detector,
    // detector = registered + (reg_dx, reg_dy). From the FITS keywords REGDX/REGDY.
    // When present, objects fixed on the detector (warm pixels, defects) are flagged.
    bool has_reg_offset = false;
    double reg_dx = 0.0, reg_dy = 0.0;
    std::string source;      // file name or other provenance
};

// Thrown when a progress callback asked to stop.
struct Cancelled : std::runtime_error {
    Cancelled() : std::runtime_error("cancelled") {}
};

// Progress callback: fraction in [0,1] and a short stage name.
// Return false to cancel the running operation (a Cancelled exception follows).
using ProgressFn = std::function<bool(double fraction, const std::string& stage)>;

}  // namespace fainttrail
