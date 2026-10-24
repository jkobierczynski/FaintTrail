// SPDX-License-Identifier: GPL-3.0-or-later
// Copyright 2026 Jurgen Kobierczynski
//
// Synthetic observing-sequence generator. Produces registered frames with a
// star field, sky background, photon + read noise, hot pixels, cosmic rays,
// seeing and transparency variations, and trailed moving objects whose true
// positions and expected signal-to-noise ratios are known.
#pragma once

#include <cstdint>
#include <string>
#include <vector>

#include "fainttrail/types.hpp"

namespace fainttrail {

struct SimMover {
    std::string id;
    double x_ref = 0.0, y_ref = 0.0;  // 0-based position at the reference epoch
    double vx_pxh = 0.0, vy_pxh = 0.0; // motion in pixels per hour
    double mag = 20.0;
};

struct SimParams {
    int width = 512;
    int height = 512;
    int nframes = 60;
    double exptime_s = 45.0;
    double cadence_s = 60.0;           // start-to-start interval
    double mjd_start = 61337.9;        // 2026-10-24 21:36 UTC
    double ra_deg = 30.0, dec_deg = 12.0;  // field centre (near the ecliptic)
    double pixscale_arcsec = 1.5;
    double fwhm_px = 3.0;
    double fwhm_jitter = 0.08;         // relative frame-to-frame seeing variation (1 sigma)
    double transparency_jitter = 0.02; // relative frame-to-frame throughput variation
    double sky_e_per_s = 40.0;         // sky electrons per pixel per second
    double read_noise_e = 3.0;
    double gain_e_per_adu = 1.0;
    double bias_adu = 500.0;
    double saturation_adu = 65535.0;
    double zero_point = 22.0;          // magnitude giving 1 e-/s
    int nstars = 800;
    double star_mag_bright = 11.0, star_mag_faint = 21.0;
    int nmovers = 8;                   // random movers (ignored if movers are given)
    double mover_mag_bright = 18.5, mover_mag_faint = 21.0;
    double rate_min_arcsec_min = 0.2, rate_max_arcsec_min = 1.0;
    int hot_pixels = 50;
    double cosmic_rays_per_frame = 4.0;
    double sky_gradient = 0.03;        // fractional sky change across the field
    std::uint64_t seed = 1;
};

struct SimTruth {
    SimMover mover;
    double flux_e = 0.0;              // electrons per exposure
    double expected_snr_frame = 0.0;  // matched-filter S/N in one frame
    double expected_snr_stack = 0.0;  // S/N after stacking all frames on the true motion
};

struct SimResult {
    std::vector<Frame> frames;
    std::vector<SimTruth> truth;
    double t_ref_mjd = 0.0;           // reference epoch for x_ref / y_ref
};

// Generate a sequence. If `movers` is non-empty it is used instead of random movers.
SimResult simulate(const SimParams& params, const std::vector<SimMover>& movers = {});

// Matched-filter S/N of a point source of `flux_e` electrons on a background of
// `bg_var_e2` electrons^2 per pixel, with a Gaussian PSF of the given FWHM.
double point_source_snr(double flux_e, double bg_var_e2, double fwhm_px);

// Rate conversion helpers.
inline double arcsec_per_min_to_px_per_h(double rate, double pixscale) { return rate * 60.0 / pixscale; }
inline double px_per_h_to_arcsec_per_min(double v, double pixscale) { return v * pixscale / 60.0; }

}  // namespace fainttrail
