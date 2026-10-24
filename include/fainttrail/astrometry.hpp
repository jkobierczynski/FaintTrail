// SPDX-License-Identifier: GPL-3.0-or-later
// Copyright 2026 Jurgen Kobierczynski
//
// Borrowing a plate solution from another image of the same field.
//
// Typical use: the observer plate-solves one frame or a stack with ASTAP,
// astrometry.net or Siril. That image may differ from fainttrail's registered
// frames in binning, cropping, orientation (rotated by 90/180/270 degrees or
// mirrored) and by a small offset. transfer_wcs() detects stars in both images,
// matches them and maps the solution onto the target pixel grid.
#pragma once

#include <string>
#include <vector>

#include "fainttrail/types.hpp"

namespace fainttrail {

struct Star {
    double x = 0.0, y = 0.0;  // 0-based centroid
    double flux = 0.0;
};

// Detect point sources, brightest first. Saturated or very extended sources are kept;
// they only affect the ordering.
std::vector<Star> detect_stars(const Image& img, int max_stars = 300, double nsigma = 8.0);

// Linear map from source pixels (0-based) to target pixels:
//   target = A * source + b, with A a 2x2 matrix (row major).
struct PixelMap {
    double a11 = 1, a12 = 0, a21 = 0, a22 = 1;
    double bx = 0, by = 0;
    void apply(double x, double y, double& tx, double& ty) const {
        tx = a11 * x + a12 * y + bx;
        ty = a21 * x + a22 * y + by;
    }
    PixelMap inverse() const;
};

struct TransferReport {
    bool ok = false;
    int matched = 0;            // stars used in the final fit
    double rms_px = 0.0;        // residual in target pixels
    double scale = 0.0;         // target pixels per source pixel
    double rotation_deg = 0.0;  // rotation of the map
    bool mirrored = false;
    bool source_has_distortion = false;  // SIP or similar terms were present and are ignored
    std::string message;
    PixelMap map;               // source -> target
};

// Match the stars of `source` (the plate-solved image) and `target` and return the
// pixel map between them. `scale_hint` is target pixels per source pixel; 0 = derive
// it from the image sizes (e.g. 0.5 for a full-resolution stack and 2x2-binned frames).
TransferReport match_images(const Image& source, const Image& target, double scale_hint = 0.0);

// Transfer `source_wcs` (valid for `source`) to the pixel grid of `target`.
// Throws std::runtime_error if the images cannot be matched.
Wcs transfer_wcs(const Wcs& source_wcs, const Image& source, const Image& target,
                 TransferReport* report = nullptr, double scale_hint = 0.0);

// Convenience: read a plate-solved FITS file and transfer its WCS onto every frame.
// The frames must be registered; the middle frame is used for the star match.
TransferReport apply_wcs_from_file(const std::string& solved_fits, std::vector<Frame>& frames,
                                   double scale_hint = 0.0);

}  // namespace fainttrail
