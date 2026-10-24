// SPDX-License-Identifier: GPL-3.0-or-later
// Copyright 2026 Jurgen Kobierczynski
//
// Minimal, dependency-free FITS reader/writer for single-HDU 2-D images.
// Reads BITPIX 8/16/32/64/-32/-64 (with BSCALE/BZERO); writes BITPIX -32.
#pragma once

#include <string>
#include <utility>
#include <vector>

#include "fainttrail/types.hpp"

namespace fainttrail {

// An ordered list of header cards (keyword, raw value text, comment).
struct FitsHeader {
    struct Card {
        std::string key;
        std::string value;    // raw value field, strings keep their quotes stripped
        std::string comment;
        bool is_string = false;
    };
    std::vector<Card> cards;

    const Card* find(const std::string& key) const;
    bool has(const std::string& key) const { return find(key) != nullptr; }
    double get_double(const std::string& key, double fallback) const;
    long long get_int(const std::string& key, long long fallback) const;
    std::string get_string(const std::string& key, const std::string& fallback = {}) const;

    void set(const std::string& key, double value, const std::string& comment = {});
    void set_int(const std::string& key, long long value, const std::string& comment = {});
    void set(const std::string& key, const std::string& value, const std::string& comment = {});
    void add_history(const std::string& text);
};

struct FitsImage {
    Image image;
    FitsHeader header;
};

// Read the primary HDU (or the first image extension if the primary is empty).
// For 3-D cubes the first plane is returned. Throws std::runtime_error on failure.
FitsImage read_fits(const std::string& path);

// Write a float image as BITPIX = -32. Structural keywords are generated;
// any SIMPLE/BITPIX/NAXIS* cards in `header` are ignored.
void write_fits(const std::string& path, const Image& image, const FitsHeader& header = {});

// TAN plate solution from a header (CD, PC+CDELT or CDELT+CROTA2); invalid if absent.
Wcs wcs_from_header(const FitsHeader& h);

// Convert a FITS header into frame metadata: mid-exposure time and TAN WCS.
// Understands MJD-MID, MJD-OBS, DATE-OBS (+EXPTIME/EXPOSURE), JD, and
// CDi_j / PCi_j+CDELTi / CDELTi+CROTA2 WCS forms.
Frame frame_from_fits(FitsImage&& fits, const std::string& source);

// Load a set of FITS files as frames, sorted by mid-exposure time.
// Throws if sizes differ or a time cannot be determined.
std::vector<Frame> load_frames(const std::vector<std::string>& paths);

// Write a frame back out with its time and WCS keywords.
void write_frame(const std::string& path, const Frame& frame,
                 const FitsHeader& extra = {});

// Calendar helpers (UTC, proleptic Gregorian).
double mjd_from_iso(const std::string& iso);  // "2026-10-24T01:02:03.4"
std::string iso_from_mjd(double mjd, int decimals = 3);

}  // namespace fainttrail
