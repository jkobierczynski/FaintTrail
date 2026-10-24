// SPDX-License-Identifier: GPL-3.0-or-later
// Copyright 2026 Jurgen Kobierczynski
#include "fainttrail/simulate.hpp"

#include <algorithm>
#include <cmath>
#include <string>

#include "psf.hpp"
#include "random.hpp"

namespace fainttrail {

using detail::kFwhmToSigma;
using detail::kPi;

double point_source_snr(double flux_e, double bg_var_e2, double fwhm_px) {
    const double s = fwhm_px * kFwhmToSigma;
    const double noise2 = flux_e + 4.0 * kPi * s * s * bg_var_e2;
    return noise2 > 0.0 ? flux_e / std::sqrt(noise2) : 0.0;
}

namespace {

struct Star {
    double x, y, flux_e_s;
};

// Magnitudes drawn from N(<m) ~ 10^(0.3 m), a typical star-count slope.
double draw_star_mag(detail::Rng& rng, double bright, double faint) {
    const double a = std::pow(10.0, 0.3 * bright), b = std::pow(10.0, 0.3 * faint);
    return std::log10(a + (b - a) * rng.uniform()) / 0.3;
}

}  // namespace

SimResult simulate(const SimParams& p, const std::vector<SimMover>& movers_in) {
    detail::Rng rng(p.seed);
    SimResult out;
    const int w = p.width, h = p.height;
    const int n = std::max(1, p.nframes);

    // timing
    std::vector<double> mid(static_cast<std::size_t>(n));
    for (int i = 0; i < n; ++i)
        mid[static_cast<std::size_t>(i)] = p.mjd_start + (i * p.cadence_s + 0.5 * p.exptime_s) / 86400.0;
    out.t_ref_mjd = 0.5 * (mid.front() + mid.back());

    // star field
    std::vector<Star> stars;
    stars.reserve(static_cast<std::size_t>(p.nstars));
    for (int i = 0; i < p.nstars; ++i) {
        const double m = draw_star_mag(rng, p.star_mag_bright, p.star_mag_faint);
        stars.push_back({rng.uniform(-5.0, w + 5.0), rng.uniform(-5.0, h + 5.0),
                         std::pow(10.0, -0.4 * (m - p.zero_point))});
    }

    // moving objects
    std::vector<SimMover> movers = movers_in;
    if (movers.empty()) {
        const int margin = 24;
        for (int i = 0; i < p.nmovers; ++i) {
            SimMover m;
            m.id = "S" + std::to_string(i + 1);
            for (int attempt = 0; attempt < 200; ++attempt) {
                m.x_ref = rng.uniform(margin, w - margin);
                m.y_ref = rng.uniform(margin, h - margin);
                bool clear = true;
                for (const auto& o : movers)
                    if (std::hypot(o.x_ref - m.x_ref, o.y_ref - m.y_ref) < 40.0) clear = false;
                if (clear) break;
            }
            const double rate = rng.uniform(p.rate_min_arcsec_min, p.rate_max_arcsec_min);
            const double v = arcsec_per_min_to_px_per_h(rate, p.pixscale_arcsec);
            const double ang = rng.uniform(0.0, 2.0 * kPi);
            m.vx_pxh = v * std::cos(ang);
            m.vy_pxh = v * std::sin(ang);
            m.mag = rng.uniform(p.mover_mag_bright, p.mover_mag_faint);
            movers.push_back(m);
        }
    }

    const double sky_e = p.sky_e_per_s * p.exptime_s;
    const double bg_var = sky_e + p.read_noise_e * p.read_noise_e;
    for (const auto& m : movers) {
        SimTruth t;
        t.mover = m;
        t.flux_e = std::pow(10.0, -0.4 * (m.mag - p.zero_point)) * p.exptime_s;
        t.expected_snr_frame = point_source_snr(t.flux_e, bg_var, p.fwhm_px);
        t.expected_snr_stack = t.expected_snr_frame * std::sqrt(static_cast<double>(n));
        out.truth.push_back(t);
    }

    // static detector defects
    struct Hot { int x, y; float adu; };
    std::vector<Hot> hot;
    for (int i = 0; i < p.hot_pixels; ++i)
        hot.push_back({rng.uniform_int(0, w - 1), rng.uniform_int(0, h - 1),
                       static_cast<float>(rng.uniform(200.0, 6000.0))});

    const Wcs wcs = Wcs::make_tan(p.ra_deg, p.dec_deg, 0.5 * w + 0.5, 0.5 * h + 0.5,
                                  p.pixscale_arcsec, 0.0);

    out.frames.reserve(static_cast<std::size_t>(n));
    for (int i = 0; i < n; ++i) {
        const double fwhm = std::max(0.6 * p.fwhm_px, p.fwhm_px * (1.0 + p.fwhm_jitter * rng.normal()));
        const double sigma = fwhm * kFwhmToSigma;
        const double transp = std::max(0.5, 1.0 + p.transparency_jitter * rng.normal());
        const double sky_i = sky_e * (1.0 + 0.01 * rng.normal());

        Image e(w, h, 0.0f);  // electrons from sources
        for (const auto& s : stars)
            detail::add_gaussian(e, s.x, s.y, s.flux_e_s * p.exptime_s * transp, sigma);

        const double t_start_h = (mid[static_cast<std::size_t>(i)] - out.t_ref_mjd) * 24.0 - 0.5 * p.exptime_s / 3600.0;
        const double t_end_h = t_start_h + p.exptime_s / 3600.0;
        for (const auto& t : out.truth) {
            const auto& m = t.mover;
            detail::add_trailed(e, m.x_ref + m.vx_pxh * t_start_h, m.y_ref + m.vy_pxh * t_start_h,
                                m.x_ref + m.vx_pxh * t_end_h, m.y_ref + m.vy_pxh * t_end_h,
                                t.flux_e * transp, sigma);
        }

        Frame f;
        f.image = Image(w, h);
        for (int y = 0; y < h; ++y) {
            const double gy = 0.5 * p.sky_gradient * (static_cast<double>(y) / h - 0.5);
            for (int x = 0; x < w; ++x) {
                const double g = 1.0 + p.sky_gradient * (static_cast<double>(x) / w - 0.5) + gy;
                const double mean_e = sky_i * g + e.at(x, y);
                const double electrons = rng.poisson(mean_e) + p.read_noise_e * rng.normal();
                f.image.at(x, y) = static_cast<float>(electrons / p.gain_e_per_adu + p.bias_adu);
            }
        }
        for (const auto& hp : hot) f.image.at(hp.x, hp.y) += hp.adu;

        // cosmic rays: short random hits of 1-4 pixels
        const int ncr = static_cast<int>(rng.poisson(p.cosmic_rays_per_frame));
        for (int c = 0; c < ncr; ++c) {
            int x = rng.uniform_int(0, w - 1), y = rng.uniform_int(0, h - 1);
            const int len = rng.uniform_int(1, 4);
            const double e_hit = rng.uniform(800.0, 15000.0);
            const int dx = rng.uniform_int(-1, 1), dy = rng.uniform_int(-1, 1);
            for (int k = 0; k < len && f.image.inside(x, y); ++k, x += dx, y += dy)
                f.image.at(x, y) += static_cast<float>(e_hit / p.gain_e_per_adu / (k + 1));
        }
        for (auto& v : f.image.data) v = std::min(v, static_cast<float>(p.saturation_adu));

        f.mjd_mid = mid[static_cast<std::size_t>(i)];
        f.exptime_s = p.exptime_s;
        f.wcs = wcs;
        f.source = "sim_" + std::to_string(i + 1);
        out.frames.push_back(std::move(f));
    }
    return out;
}

}  // namespace fainttrail
