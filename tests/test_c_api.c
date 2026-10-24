/* SPDX-License-Identifier: GPL-3.0-or-later
 * Copyright 2026 Jurgen Kobierczynski
 *
 * Exercises the C API from plain C: frames with Gaussian noise and one moving
 * point source are handed over pixel by pixel, as a host application would.
 */
#include <math.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>

#include "fainttrail/fainttrail.h"

static unsigned long long rng_state = 88172645463325252ULL;
static double uniform01(void) {
    rng_state ^= rng_state << 13;
    rng_state ^= rng_state >> 7;
    rng_state ^= rng_state << 17;
    return (double)(rng_state >> 11) * (1.0 / 9007199254740992.0);
}
static double gauss(void) {
    double u = uniform01(), v = uniform01();
    if (u < 1e-300) u = 1e-300;
    return sqrt(-2.0 * log(u)) * cos(6.283185307179586 * v);
}

static int progress_calls = 0;
static int on_progress(double f, const char* stage, void* user) {
    (void)f;
    (void)stage;
    (void)user;
    ++progress_calls;
    return 1;
}

int run_c_api_test(void) {
    const int w = 128, h = 128, n = 30;
    const double sigma = 10.0, fwhm = 3.0, s = fwhm / 2.3548;
    const double x0 = 60.0, y0 = 64.0, vx = 18.0, vy = -7.0; /* px/h, 0-based at t_ref */
    const double flux = 9.0 / sqrt((double)n) * sigma * sqrt(4.0 * 3.141592653589793) * s; /* stack S/N 9 */
    const double mjd0 = 61337.9, cadence_d = 60.0 / 86400.0;
    const double tref = mjd0 + 0.5 * (n - 1) * cadence_d;
    float* img = (float*)malloc(sizeof(float) * (size_t)(w * h));
    ft_context* ctx = ft_create();
    ft_params p;
    ft_candidate c;
    int i, x, y, k, found = 0, status;
    if (!ctx || !img) return 1;

    for (i = 0; i < n; ++i) {
        const double mjd = mjd0 + i * cadence_d;
        const double t = (mjd - tref) * 24.0;
        const double cx = x0 + vx * t, cy = y0 + vy * t;
        for (y = 0; y < h; ++y)
            for (x = 0; x < w; ++x) {
                const double r2 = (x - cx) * (x - cx) + (y - cy) * (y - cy);
                img[y * w + x] = (float)(1000.0 + sigma * gauss() + flux * exp(-0.5 * r2 / (s * s)) / (2 * 3.141592653589793 * s * s));
            }
        if (ft_add_frame(ctx, w, h, img, mjd, 30.0) != FT_OK) return 2;
    }
    if (ft_frame_count(ctx) != n) return 3;

    {
        const double crval[2] = {45.0, 10.0}, crpix[2] = {64.5, 64.5};
        const double cd[4] = {-1.5 / 3600.0, 0.0, 0.0, 1.5 / 3600.0};
        if (ft_set_wcs(ctx, crval, crpix, cd) != FT_OK) return 4;
    }

    ft_params_default(&p);
    p.vmin_pxh = 6.0;
    p.vmax_pxh = 25.0;
    p.fwhm_px = fwhm;
    p.threshold_snr = 6.0; /* fixed threshold: no null runs needed */
    status = ft_search(ctx, &p, on_progress, NULL);
    if (status != FT_OK) {
        printf("    ft_search failed: %s\n", ft_last_error(ctx));
        return 5;
    }
    if (progress_calls == 0) return 6;
    for (k = 0; k < ft_candidate_count(ctx); ++k) {
        if (ft_get_candidate(ctx, k, &c) != FT_OK) return 7;
        /* the C API reports FITS 1-based pixel coordinates */
        if (c.ok && fabs(c.x - (x0 + 1.0)) < 1.0 && fabs(c.y - (y0 + 1.0)) < 1.0 && fabs(c.vx_pxh - vx) < 3.0 &&
            fabs(c.vy_pxh - vy) < 3.0) {
            double mjd, tx, ty, ra, dec;
            found = 1;
            if (ft_get_track_point(ctx, k, 0, &mjd, &tx, &ty, &ra, &dec) != FT_OK) return 8;
            if (fabs(mjd - mjd0) > 1e-9 || fabs(ra - 45.0) > 0.05 || fabs(dec - 10.0) > 0.05) return 9;
            if (c.rate_arcsec_min < 0.4 || c.rate_arcsec_min > 0.6) return 10; /* 19.3 px/h * 1.5"/60 = 0.48 */
        }
    }
    if (!found) {
        printf("    injected object not recovered (%d candidates)\n", ft_candidate_count(ctx));
        return 11;
    }
    if (strlen(ft_version()) == 0) return 12;
    ft_destroy(ctx);
    free(img);
    return 0;
}
