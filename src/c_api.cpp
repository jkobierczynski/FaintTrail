// SPDX-License-Identifier: GPL-3.0-or-later
// Copyright 2026 Jurgen Kobierczynski
#include "fainttrail/fainttrail.h"

#include <algorithm>
#include <cstring>
#include <memory>
#include <string>

#include "fainttrail/fainttrail.hpp"

#ifndef FAINTTRAIL_VERSION_STRING
#define FAINTTRAIL_VERSION_STRING "0.0.0"
#endif

namespace fainttrail {
const char* version() { return FAINTTRAIL_VERSION_STRING; }
}  // namespace fainttrail

struct ft_context {
    std::vector<fainttrail::Frame> frames;
    fainttrail::Wcs wcs;
    std::unique_ptr<fainttrail::Prepared> prep;
    fainttrail::SearchResult result;
    bool have_result = false;
    std::string error;
};

namespace {

fainttrail::SearchParams to_cpp(const ft_params& p) {
    fainttrail::SearchParams s;
    s.vmin_pxh = p.vmin_pxh;
    s.vmax_pxh = p.vmax_pxh;
    s.grid_tol_px = p.grid_tol_px;
    s.fwhm_px = p.fwhm_px;
    s.clip_sigma = p.clip_sigma;
    s.cr_sigma = p.cr_sigma;
    s.mask_scatter = p.mask_scatter;
    s.weight_scatter = p.weight_scatter;
    s.static_exclude = p.static_exclude;
    s.min_coverage = p.min_coverage;
    s.threshold_snr = p.threshold_snr;
    s.min_snr_floor = p.min_snr_floor;
    s.target_false_per_run = p.target_false_per_run;
    s.max_frame_fraction = p.max_frame_fraction;
    s.half_tolerance = p.half_tolerance;
    s.null_runs = p.null_runs;
    s.max_candidates = p.max_candidates;
    s.edge_px = p.edge_px;
    s.threads = p.threads;
    s.seed = p.seed;
    s.backend = std::string(p.backend, strnlen(p.backend, sizeof p.backend));
    if (s.backend.empty()) s.backend = "auto";
    return s;
}

}  // namespace

extern "C" {

const char* ft_version(void) { return fainttrail::version(); }

void ft_params_default(ft_params* p) {
    if (!p) return;
    const fainttrail::SearchParams s;
    p->vmin_pxh = s.vmin_pxh;
    p->vmax_pxh = s.vmax_pxh;
    p->grid_tol_px = s.grid_tol_px;
    p->fwhm_px = s.fwhm_px;
    p->clip_sigma = s.clip_sigma;
    p->cr_sigma = s.cr_sigma;
    p->mask_scatter = s.mask_scatter;
    p->weight_scatter = s.weight_scatter;
    p->static_exclude = s.static_exclude;
    p->min_coverage = s.min_coverage;
    p->threshold_snr = s.threshold_snr;
    p->min_snr_floor = s.min_snr_floor;
    p->target_false_per_run = s.target_false_per_run;
    p->max_frame_fraction = s.max_frame_fraction;
    p->half_tolerance = s.half_tolerance;
    p->null_runs = s.null_runs;
    p->max_candidates = s.max_candidates;
    p->edge_px = s.edge_px;
    p->threads = s.threads;
    p->seed = s.seed;
    std::memset(p->backend, 0, sizeof p->backend);
    std::strncpy(p->backend, s.backend.c_str(), sizeof p->backend - 1);
}

ft_context* ft_create(void) {
    try {
        return new ft_context();
    } catch (...) {
        return nullptr;
    }
}

void ft_destroy(ft_context* ctx) { delete ctx; }

void ft_clear_frames(ft_context* ctx) {
    if (!ctx) return;
    ctx->frames.clear();
    ctx->prep.reset();
    ctx->have_result = false;
}

int ft_add_frame(ft_context* ctx, int width, int height, const float* pixels, double mjd_mid, double exptime_s) {
    if (!ctx || !pixels || width <= 0 || height <= 0) return FT_ERROR;
    try {
        fainttrail::Frame f;
        f.image = fainttrail::Image(width, height);
        std::memcpy(f.image.data.data(), pixels, f.image.data.size() * sizeof(float));
        f.mjd_mid = mjd_mid;
        f.exptime_s = exptime_s;
        ctx->frames.push_back(std::move(f));
        ctx->have_result = false;
        return FT_OK;
    } catch (const std::exception& e) {
        ctx->error = e.what();
        return FT_ERROR;
    }
}

int ft_add_fits(ft_context* ctx, const char* path) {
    if (!ctx || !path) return FT_ERROR;
    try {
        ctx->frames.push_back(fainttrail::frame_from_fits(fainttrail::read_fits(path), path));
        ctx->have_result = false;
        return FT_OK;
    } catch (const std::exception& e) {
        ctx->error = e.what();
        return FT_ERROR;
    }
}

int ft_frame_count(const ft_context* ctx) { return ctx ? static_cast<int>(ctx->frames.size()) : 0; }

int ft_set_wcs(ft_context* ctx, const double crval[2], const double crpix[2], const double cd[4]) {
    if (!ctx || !crval || !crpix || !cd) return FT_ERROR;
    fainttrail::Wcs w;
    w.crval1 = crval[0];
    w.crval2 = crval[1];
    w.crpix1 = crpix[0];
    w.crpix2 = crpix[1];
    w.cd11 = cd[0];
    w.cd12 = cd[1];
    w.cd21 = cd[2];
    w.cd22 = cd[3];
    w.valid = (w.cd11 * w.cd22 - w.cd12 * w.cd21) != 0.0;
    ctx->wcs = w;
    return w.valid ? FT_OK : FT_ERROR;
}

int ft_set_wcs_from_fits(ft_context* ctx, const char* path, int* matched) {
    if (!ctx || !path || ctx->frames.empty()) return FT_ERROR;
    try {
        std::vector<fainttrail::Frame> frames = ctx->frames;
        std::stable_sort(frames.begin(), frames.end(),
                         [](const fainttrail::Frame& a, const fainttrail::Frame& b) { return a.mjd_mid < b.mjd_mid; });
        const auto rep = fainttrail::apply_wcs_from_file(path, frames, 0.0);
        if (matched) *matched = rep.matched;
        ctx->wcs = frames.front().wcs;
        ctx->have_result = false;
        return FT_OK;
    } catch (const std::exception& e) {
        ctx->error = e.what();
        return FT_ERROR;
    }
}

int ft_search(ft_context* ctx, const ft_params* params, ft_progress_fn progress, void* user) {
    if (!ctx) return FT_ERROR;
    try {
        ft_params defaults;
        ft_params_default(&defaults);
        const fainttrail::SearchParams sp = to_cpp(params ? *params : defaults);
        std::vector<fainttrail::Frame> frames = ctx->frames;
        std::stable_sort(frames.begin(), frames.end(),
                         [](const fainttrail::Frame& a, const fainttrail::Frame& b) { return a.mjd_mid < b.mjd_mid; });
        if (ctx->wcs.valid)
            for (auto& f : frames) f.wcs = ctx->wcs;
        fainttrail::ProgressFn cb;
        if (progress)
            cb = [progress, user](double f, const std::string& s) { return progress(f, s.c_str(), user) != 0; };
        fainttrail::ProgressFn prep_cb, search_cb;
        if (cb) {
            prep_cb = [&cb](double f, const std::string& s) { return cb(0.1 * f, s); };
            search_cb = [&cb](double f, const std::string& s) { return cb(0.1 + 0.9 * f, s); };
        }
        ctx->prep = std::make_unique<fainttrail::Prepared>(fainttrail::prepare(frames, sp, prep_cb));
        ctx->result = fainttrail::search(*ctx->prep, sp, search_cb);
        ctx->have_result = true;
        return FT_OK;
    } catch (const fainttrail::Cancelled&) {
        ctx->error = "cancelled";
        ctx->have_result = false;
        return FT_CANCELLED;
    } catch (const std::exception& e) {
        ctx->error = e.what();
        ctx->have_result = false;
        return FT_ERROR;
    }
}

double ft_reference_mjd(const ft_context* ctx) {
    return (ctx && ctx->have_result) ? ctx->result.t_ref_mjd : 0.0;
}

double ft_threshold(const ft_context* ctx) {
    return (ctx && ctx->have_result) ? ctx->result.threshold_snr : 0.0;
}

int ft_candidate_count(const ft_context* ctx) {
    return (ctx && ctx->have_result) ? static_cast<int>(ctx->result.candidates.size()) : 0;
}

int ft_get_candidate(const ft_context* ctx, int index, ft_candidate* out) {
    if (!ctx || !out || !ctx->have_result || index < 0 ||
        index >= static_cast<int>(ctx->result.candidates.size()))
        return FT_ERROR;
    const auto& c = ctx->result.candidates[static_cast<std::size_t>(index)];
    std::memset(out, 0, sizeof *out);
    out->id = c.id;
    out->ok = c.ok ? 1 : 0;
    out->x = c.x + 1.0;
    out->y = c.y + 1.0;
    out->vx_pxh = c.vx_pxh;
    out->vy_pxh = c.vy_pxh;
    out->pos_err_px = c.pos_err_px;
    out->vel_err_pxh = c.vel_err_pxh;
    out->snr = c.snr;
    out->snr_coarse = c.snr_coarse;
    out->snr_half1 = c.snr_half1;
    out->snr_half2 = c.snr_half2;
    out->max_frame_fraction = c.max_frame_fraction;
    out->false_alarms_per_run = c.false_alarms_per_run;
    out->flux = c.flux;
    out->mag_inst = c.mag_inst;
    out->ra_deg = c.ra_deg;
    out->dec_deg = c.dec_deg;
    out->rate_arcsec_min = c.rate_arcsec_min;
    out->pa_deg = c.pa_deg;
    std::strncpy(out->flags, c.flags.c_str(), sizeof(out->flags) - 1);
    return FT_OK;
}

int ft_get_track_point(const ft_context* ctx, int index, int frame, double* mjd, double* x, double* y,
                       double* ra_deg, double* dec_deg) {
    if (!ctx || !ctx->have_result || !ctx->prep || index < 0 ||
        index >= static_cast<int>(ctx->result.candidates.size()) || frame < 0 || frame >= ctx->prep->nframes)
        return FT_ERROR;
    const auto pts = fainttrail::track_points(*ctx->prep, ctx->result.candidates[static_cast<std::size_t>(index)]);
    const auto& p = pts[static_cast<std::size_t>(frame)];
    if (mjd) *mjd = p.mjd;
    if (x) *x = p.x + 1.0;
    if (y) *y = p.y + 1.0;
    if (ra_deg) *ra_deg = p.ra_deg;
    if (dec_deg) *dec_deg = p.dec_deg;
    return FT_OK;
}

const char* ft_last_error(const ft_context* ctx) { return ctx ? ctx->error.c_str() : "no context"; }

}  // extern "C"
