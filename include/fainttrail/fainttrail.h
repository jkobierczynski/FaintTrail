/* SPDX-License-Identifier: GPL-3.0-or-later
 * Copyright 2026 Jurgen Kobierczynski
 *
 * FaintTrail C API — a stable, plain-C interface for host applications such as
 * AstroFind (C++/Qt) or Siril (C/GTK).
 *
 * Typical use:
 *   ft_context* ctx = ft_create();
 *   for each registered frame: ft_add_frame(ctx, w, h, pixels, mjd_mid, exptime);
 *   ft_set_wcs(ctx, crval, crpix, cd);              // optional, enables RA/Dec
 *   ft_params p; ft_params_default(&p); p.vmax_pxh = 60;
 *   if (ft_search(ctx, &p, my_progress, my_data) == FT_OK)
 *       for (i = 0; i < ft_candidate_count(ctx); ++i) ft_get_candidate(ctx, i, &c);
 *   ft_destroy(ctx);
 *
 * Pixel coordinates in this API use the FITS convention (1-based).
 * Velocities are in pixels per hour, times in UTC MJD.
 */
#ifndef FAINTTRAIL_H
#define FAINTTRAIL_H

#ifdef __cplusplus
extern "C" {
#endif

#define FT_OK 0
#define FT_ERROR 1
#define FT_CANCELLED 2

typedef struct ft_context ft_context;

typedef struct ft_params {
    double vmin_pxh;              /* slowest motion searched, px/hour            */
    double vmax_pxh;              /* fastest motion searched, px/hour            */
    double grid_tol_px;           /* allowed drift of a wrong grid velocity, px  */
    double fwhm_px;               /* PSF FWHM, 0 = estimate from the stars        */
    double clip_sigma;
    double cr_sigma;              /* cosmic-ray rejection level                  */
    double mask_scatter;
    double weight_scatter;
    double min_coverage;
    double threshold_snr;         /* 0 = calibrate with time-shuffled nulls       */
    double min_snr_floor;
    double target_false_per_run;  /* accepted noise detections per search         */
    double max_frame_fraction;
    double half_tolerance;
    int static_exclude;           /* frames each side left out of the static sky, -1 = auto */
    int null_runs;
    int max_candidates;
    int edge_px;
    int threads;                  /* 0 = all cores                                */
    unsigned long long seed;
} ft_params;

typedef struct ft_candidate {
    int id;
    int ok;                       /* 1 = passed every consistency test            */
    double x, y;                  /* FITS pixel position at the reference epoch   */
    double vx_pxh, vy_pxh;
    double pos_err_px, vel_err_pxh; /* 1-sigma per axis                           */
    double snr, snr_coarse, snr_half1, snr_half2;
    double max_frame_fraction;
    double false_alarms_per_run;
    double flux, mag_inst;
    double ra_deg, dec_deg;       /* 0 when no WCS was given                      */
    double rate_arcsec_min, pa_deg;
    char flags[64];
} ft_candidate;

/* Return 0 to cancel the running search. */
typedef int (*ft_progress_fn)(double fraction, const char* stage, void* user);

const char* ft_version(void);
void ft_params_default(ft_params* p);

ft_context* ft_create(void);
void ft_destroy(ft_context* ctx);
void ft_clear_frames(ft_context* ctx);

/* Copies the pixels. Frames must be registered to a common grid. */
int ft_add_frame(ft_context* ctx, int width, int height, const float* pixels, double mjd_mid,
                 double exptime_s);
/* Load a FITS file (time from MJD-MID / MJD-OBS / DATE-OBS, TAN WCS if present). */
int ft_add_fits(ft_context* ctx, const char* path);
int ft_frame_count(const ft_context* ctx);

/* TAN WCS of the registered grid: crval (deg), crpix (1-based), cd (deg/px, row major). */
int ft_set_wcs(ft_context* ctx, const double crval[2], const double crpix[2], const double cd[4]);

/* Borrow the plate solution of a solved FITS image of the same field (stack or frame,
 * any binning or 90-degree orientation). Call after adding all frames; returns FT_OK
 * when the star patterns match. `matched` (may be NULL) receives the star count. */
int ft_set_wcs_from_fits(ft_context* ctx, const char* path, int* matched);

int ft_search(ft_context* ctx, const ft_params* params, ft_progress_fn progress, void* user);

double ft_reference_mjd(const ft_context* ctx);
double ft_threshold(const ft_context* ctx);
int ft_candidate_count(const ft_context* ctx);
int ft_get_candidate(const ft_context* ctx, int index, ft_candidate* out);
/* Predicted position of candidate `index` in frame `frame` (0-based, time-sorted). */
int ft_get_track_point(const ft_context* ctx, int index, int frame, double* mjd, double* x, double* y,
                       double* ra_deg, double* dec_deg);

const char* ft_last_error(const ft_context* ctx);

#ifdef __cplusplus
}
#endif

#endif /* FAINTTRAIL_H */
