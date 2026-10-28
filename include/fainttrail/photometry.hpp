// SPDX-License-Identifier: GPL-3.0-or-later
// Copyright 2026 Jurgen Kobierczynski
//
// Photometric calibration against Gaia DR3.
//
// Candidates are measured with a PSF-weighted flux on the stacked residuals
// (Candidate::flux). Field stars are measured with exactly the same estimator on
// the static sky (the median of the aligned frames), matched to Gaia DR3 through
// the plate solution, and the zero point is the clipped median of
// catalogue magnitude + 2.5 log10(flux). Only stars of roughly solar colour
// (0.5 < BP-RP < 1.3, like asteroids) are used when the catalogue has colours,
// which keeps the colour term of an unfiltered or luminance camera small.
// Magnitudes: Gaia G directly, Johnson V and Cousins R from G and BP-RP with
// the Gaia EDR3 relations (Riello et al. 2021). Trailing during an exposure
// lowers a PSF flux; candidate magnitudes are corrected for it.
#pragma once

#include <string>
#include <vector>

#include "fainttrail/search.hpp"

namespace fainttrail {

struct CatalogStar {
    double ra_deg = 0.0, dec_deg = 0.0;
    double g = 0.0;               // Gaia G
    double bp_rp = 0.0;           // Gaia BP - RP
    bool has_colour = false;
};

// Magnitude of a catalogue star in band 'G', 'V' or 'R'. NaN for V/R without colour.
double catalog_mag(const CatalogStar& s, char band);

// Read a catalogue: whitespace-, tab- or comma-separated columns RA (deg),
// Dec (deg), G, and optionally BP and RP (or BP-RP); '#' comments and header lines
// are skipped. VizieR's asu-tsv output of I/355/gaiadr3 works as it is.
std::vector<CatalogStar> read_catalog(const std::string& path);

// VizieR query for Gaia DR3 around (ra, dec): radius and faint limit.
std::string gaia_query_url(double ra_deg, double dec_deg, double radius_deg, double g_max);

// Download Gaia DR3 for the field with curl (shipped with Windows 10+, Linux and
// macOS), into cache_dir (empty: the user cache folder). A cached file for the
// same field is reused. Throws with an explanation if neither works.
std::vector<CatalogStar> fetch_gaia(double ra_deg, double dec_deg, double radius_deg, double g_max,
                                    const std::string& cache_dir = {}, std::string* cache_file = nullptr);

// The circle (centre, radius) covering a solved prepared sequence.
void field_circle(const Prepared& p, double& ra_deg, double& dec_deg, double& radius_deg);

struct PhotometryResult {
    bool ok = false;
    char band = 'V';
    double zero_point = 0.0;   // mag = zero_point - 2.5 log10(flux per frame)
    double zp_error = 0.0;     // standard error of the zero point
    double rms = 0.0;          // scatter of the calibration stars
    int stars = 0;             // used after clipping
    int matched = 0;           // catalogue stars matched to measured stars
    std::string message;
};

// Measure the field stars on p.static_sky and fit the zero point. Needs p.wcs.
PhotometryResult calibrate_photometry(const Prepared& p, const std::vector<CatalogStar>& catalog, char band = 'V');

// Calibrated magnitudes of the candidates (Candidate::mag, mag_err), including
// the trailing correction; records the band and zero point in the result.
void apply_photometry(const Prepared& p, const PhotometryResult& cal, SearchResult& res);

// PSF flux (same estimator as Candidate::flux) of a point source at (x, y) on `img`.
double psf_flux(const Image& img, double x, double y, double fwhm_px);

// Fraction of a point source's PSF flux that the estimator recovers for a source
// trailed over `length_px` during the exposure.
double trail_flux_factor(double length_px, double fwhm_px);

}  // namespace fainttrail
