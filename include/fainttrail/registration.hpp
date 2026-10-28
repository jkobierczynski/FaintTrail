// SPDX-License-Identifier: GPL-3.0-or-later
// Copyright 2026 Jurgen Kobierczynski
//
// Frame registration. Synthetic tracking needs the stars on the same pixels in
// every frame. check_registration() measures whether they are; register_frames()
// aligns a sequence on its middle frame: stars are matched, a similarity
// transform (shift, rotation, scale; meridian flips included) is fitted per
// frame, and the frame is resampled with bicubic interpolation. The offset back
// to the detector is recorded (Frame::reg_dx/reg_dy), so objects fixed on the
// detector (hot pixels, defects) are still recognised after alignment.
#pragma once

#include <string>
#include <vector>

#include "fainttrail/types.hpp"

namespace fainttrail {

// How far the star field moves from the first frame (in time) to the middle and
// last ones.
struct DriftCheck {
    bool measured = false;     // enough stars matched to tell
    bool registered = true;    // true when not measured
    double dx = 0.0, dy = 0.0;          // first -> last frame, px
    double mid_dx = 0.0, mid_dy = 0.0;  // first -> middle frame, px
    double rate_pxh = 0.0;     // |first -> last| / time between them
    int stars = 0, matched = 0;
    std::string message;       // explanation when not registered
};

// Frames registered to better than 1.5 px pass. t_h: frame times in hours.
DriftCheck check_registration(const std::vector<const Image*>& frames, const std::vector<double>& t_h);

// Alignment of one frame onto the reference: frame pixel = scale * R(rot) * ref pixel + shift.
struct FrameAlignment {
    std::string source;
    bool aligned = false;      // false: dropped (too few matched stars)
    double dx = 0.0, dy = 0.0; // where the reference frame's centre lies in this frame, minus the centre (px)
    double rot_deg = 0.0, scale = 1.0;
    double rms_px = 0.0;       // residual of the matched stars
    int stars = 0;             // stars used in the fit
};

struct AlignmentReport {
    bool applied = false;      // frames were resampled
    int reference = -1;        // index (in the input order) of the reference frame
    int dropped = 0;
    std::vector<FrameAlignment> frames;  // in the input order
    std::string summary() const;         // one line
};

// Align `frames` in place on the middle frame in time. Frames that cannot be
// aligned (clouds, too few stars) are removed and reported. Pixels outside the
// original frame become NaN. The WCS of the first frame that has one is carried
// over to the reference grid. Throws if fewer than 5 frames remain.
AlignmentReport register_frames(std::vector<Frame>& frames, int threads = 0);

}  // namespace fainttrail
