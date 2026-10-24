// SPDX-License-Identifier: GPL-3.0-or-later
// Copyright 2026 Jurgen Kobierczynski
//
// Tiled search for sequences that do not fit in memory.
//
// The field is cut into square tiles. Every tile is prepared and stacked with a
// halo wide enough that nothing near its edge depends on pixels outside it: the
// longest shift of any trial motion over the sequence, plus the reach of the
// matched filter, cosmic-ray test and mask dilation. Quantities that belong to
// a whole frame (background level, residual noise, matched-filter normalisation,
// PSF width) are measured on the whole frame first, from the same sample
// positions the untiled search uses, and the detection maps of the tiles are
// assembled before peaks are found. The threshold calibration runs the same way
// with the same time permutations.
//
// The result is therefore the same as search(prepare(frames)) bit for bit, while
// only a few tiles' worth of frames are in memory: with a FITS frame source the
// frames are read from disk tile by tile and never loaded whole.
#pragma once

#include <memory>
#include <string>
#include <vector>

#include "fainttrail/astrometry.hpp"
#include "fainttrail/search.hpp"
#include "fainttrail/types.hpp"

namespace fainttrail {

// Metadata of one frame (everything in Frame except the pixels).
struct FrameInfo {
    double mjd_mid = 0.0;
    double exptime_s = 0.0;
    Wcs wcs;
    bool has_reg_offset = false;
    double reg_dx = 0.0, reg_dy = 0.0;
    std::string source;
};

// A registered sequence whose pixels can be read a rectangle at a time.
class FrameSource {
public:
    virtual ~FrameSource() = default;
    virtual int count() const = 0;
    virtual int width() const = 0;
    virtual int height() const = 0;
    virtual const FrameInfo& info(int i) const = 0;
    // Read the rectangle (0-based x0, y0, size w x h) of frame i. Must be safe to
    // call from several threads at once.
    virtual void read(int i, int x0, int y0, int w, int h, Image& out) const = 0;
    // Give every frame this WCS (after a plate-solution transfer).
    virtual void set_wcs(const Wcs& wcs) = 0;
};

// Frames already in memory (not copied; they must outlive the source).
std::unique_ptr<FrameSource> memory_frame_source(const std::vector<Frame>& frames);
// FITS files, read from disk as needed. Only the headers are read here.
std::unique_ptr<FrameSource> fits_frame_source(const std::vector<std::string>& paths);

// Plate-solution transfer (see astrometry.hpp) for a frame source; reads the middle frame.
TransferReport apply_wcs_from_file(const std::string& solved_fits, FrameSource& frames, double scale_hint = 0.0);

struct TileOptions {
    int tile = 0;               // tile size (core, pixels); 0 = largest that fits max_memory_mb
    double max_memory_mb = 0;   // memory budget for choosing the tile size; 0 = 2048 MB
};

struct TilePlan {
    int tile = 0;               // core size
    int stack_halo = 0;         // longest shift of any trial motion, pixels
    int filter_halo = 0;        // reach of matched filter, cosmic-ray test and mask
    int tiles_x = 0, tiles_y = 0;
    double memory_mb = 0.0;     // estimated peak memory of the search
};

// Installed memory of this computer in MB (0 if unknown).
double physical_memory_mb();

// Estimated peak memory (MB) of the untiled search of n frames of w x h pixels,
// including the loaded frames.
double estimate_untiled_memory_mb(int width, int height, int nframes);

// The tiling search_tiled() would use (needs the PSF width; 0 = assume 3 px).
TilePlan plan_tiles(const FrameSource& frames, const SearchParams& params, const TileOptions& options,
                    double fwhm_px = 0.0);

// Search a sequence tile by tile. `meta`, if given, receives the frame times,
// WCS, background levels, PSF width and timings (no images), which is enough for
// write_tracklets_csv() and profile_text().
SearchResult search_tiled(const FrameSource& frames, const SearchParams& params, const TileOptions& options = {},
                          Prepared* meta = nullptr, const ProgressFn& progress = {});

}  // namespace fainttrail
