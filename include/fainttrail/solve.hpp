// SPDX-License-Identifier: GPL-3.0-or-later
// Copyright 2026 Jurgen Kobierczynski
//
// Plate solving through ASTAP (https://www.hnsky.org/astap.htm), and sky
// positions of candidates in the formats astronomers use (sexagesimal, MPC
// 80-column observation lines).
//
// FaintTrail does not ship a star catalogue. When the frames carry no plate
// solution, the static sky (the median of the aligned frames: deep, free of
// cosmic rays and moving objects) is written to a temporary FITS file and solved
// by ASTAP with the star database the user installed for it. The resulting TAN
// solution applies directly to the aligned pixel grid, so every candidate gets
// RA/Dec, a rate in arcsec/min and a position angle.
#pragma once

#include <string>
#include <vector>

#include "fainttrail/fits.hpp"
#include "fainttrail/search.hpp"
#include "fainttrail/types.hpp"

namespace fainttrail {

// Where the telescope pointed and how large the field is, from FITS keywords:
// RA/DEC (degrees), OBJCTRA/OBJCTDEC (sexagesimal), CRVAL1/2; FOCALLEN (mm) with
// XPIXSZ (micron, times XBINNING unless the pixel size already includes it) for
// the field size. Anything missing stays unset and the solver searches wider.
struct SolveHints {
    bool has_center = false;
    double ra_deg = 0.0, dec_deg = 0.0;
    double fov_deg = 0.0;            // image height in degrees, 0 = unknown
    double search_radius_deg = 0.0;  // 0: 30 deg with a centre, whole sky without
};
SolveHints hints_from_header(const FitsHeader& header, int image_height);

struct SolveOptions {
    std::string astap;         // path of astap_cli / astap; empty = search (FAINTTRAIL_ASTAP, PATH, usual places)
    std::string database_dir;  // ASTAP -d (empty = ASTAP's default)
    std::string work_dir;      // temporary files (empty = system temp directory)
    int downsample = 0;        // ASTAP -z (0 = automatic)
    int timeout_s = 300;
    bool keep_files = false;
};

struct SolveResult {
    bool ok = false;
    Wcs wcs;
    double scale_arcsec = 0.0;  // per pixel
    double rotation_deg = 0.0;
    double center_ra_deg = 0.0, center_dec_deg = 0.0;
    std::string solver;         // program that was run
    std::string message;        // why it failed, or a one-line summary
};

// The ASTAP executable, or "" if none is found. `hint`: an explicit path.
std::string find_astap(const std::string& hint = {});

// Solve one image (any pixel values; NaN allowed).
SolveResult solve_with_astap(const Image& image, const SolveHints& hints, const SolveOptions& options = {});

// Solve the static sky of a prepared sequence and install the solution
// (p.wcs). Candidates of a search that already ran get their sky coordinates
// with update_sky_coordinates().
SolveResult solve_prepared(Prepared& p, const SolveHints& hints, const SolveOptions& options = {});

// Recompute RA/Dec, rate and position angle of all candidates from p.wcs.
void update_sky_coordinates(const Prepared& p, SearchResult& res);

// "11 18 26.60" and "+13 45 01.0" (separator " " or ":"; decimals of the seconds).
std::string format_ra(double ra_deg, int decimals = 2, char sep = ' ');
std::string format_dec(double dec_deg, int decimals = 1, char sep = ' ');
// Parse "11 18 26.6", "11:18:26.6", "11h18m26.6s" (hours) or decimal degrees.
bool parse_ra(const std::string& s, double& ra_deg);
bool parse_dec(const std::string& s, double& dec_deg);

// MPC 80-column observation lines (optical, CCD) for the passing candidates:
// `per_object` positions spread over the sequence (first, middle, last frame for
// 3). Designations are temporary ("FT00001" for candidate 1). Requires p.wcs.
// Magnitude and band (columns 66-71) are written when the result was calibrated
// with apply_photometry(), else left blank.
std::vector<std::string> mpc_observations(const Prepared& p, const SearchResult& res, const std::string& obscode = "500",
                                          int per_object = 3, bool include_flagged = false);

}  // namespace fainttrail
