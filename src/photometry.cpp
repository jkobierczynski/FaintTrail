// SPDX-License-Identifier: GPL-3.0-or-later
// Copyright 2026 Jurgen Kobierczynski
//
// Photometric calibration against Gaia DR3 (photometry.hpp).
#include "fainttrail/photometry.hpp"

#include <algorithm>
#include <cmath>
#include <cstdio>
#include <cstdlib>
#include <filesystem>
#include <fstream>
#include <limits>
#include <sstream>
#include <stdexcept>

#include "common.hpp"
#include "fainttrail/astrometry.hpp"
#include "process.hpp"
#include "psf.hpp"

namespace fs = std::filesystem;

namespace fainttrail {

namespace {

constexpr double kPiP = 3.14159265358979323846;

bool to_double(const std::string& t, double& v) {
    std::size_t a = 0, b = t.size();
    while (a < b && std::isspace(static_cast<unsigned char>(t[a]))) ++a;
    while (b > a && std::isspace(static_cast<unsigned char>(t[b - 1]))) --b;
    if (a == b) return false;
    const std::string s = t.substr(a, b - a);
    char* end = nullptr;
    v = std::strtod(s.c_str(), &end);
    return end != s.c_str() && *end == '\0' && std::isfinite(v);
}

std::vector<std::string> split_fields(const std::string& line) {
    std::vector<std::string> out;
    const char sep = line.find('\t') != std::string::npos ? '\t' : (line.find(',') != std::string::npos ? ',' : 0);
    if (sep) {  // keeps empty fields (VizieR leaves missing colours empty)
        std::size_t s = 0;
        for (;;) {
            const std::size_t e = line.find(sep, s);
            out.push_back(line.substr(s, e == std::string::npos ? std::string::npos : e - s));
            if (e == std::string::npos) break;
            s = e + 1;
        }
    } else {
        std::istringstream is(line);
        std::string t;
        while (is >> t) out.push_back(t);
    }
    return out;
}

double median_of(std::vector<double> v) {
    if (v.empty()) return 0.0;
    std::nth_element(v.begin(), v.begin() + static_cast<std::ptrdiff_t>(v.size() / 2), v.end());
    double m = v[v.size() / 2];
    if (v.size() % 2 == 0) m = 0.5 * (m + *std::max_element(v.begin(), v.begin() + static_cast<std::ptrdiff_t>(v.size() / 2)));
    return m;
}

std::string cache_folder(const std::string& requested) {
    if (!requested.empty()) return requested;
#ifdef _WIN32
    if (const char* d = std::getenv("LOCALAPPDATA")) return (fs::path(d) / "FaintTrail" / "gaia").string();
#else
    if (const char* d = std::getenv("XDG_CACHE_HOME")) return (fs::path(d) / "fainttrail" / "gaia").string();
    if (const char* h = std::getenv("HOME")) return (fs::path(h) / ".cache" / "fainttrail" / "gaia").string();
#endif
    std::error_code ec;
    return (fs::temp_directory_path(ec) / "fainttrail_gaia").string();
}

}  // namespace

double catalog_mag(const CatalogStar& s, char band) {
    if (band == 'G') return s.g;
    if (!s.has_colour) return std::numeric_limits<double>::quiet_NaN();
    const double x = s.bp_rp;
    // Gaia EDR3 documentation 5.5.1 (Riello et al. 2021): G - V and G - R(Cousins)
    if (band == 'V') {
        if (x < -0.5 || x > 5.0) return std::numeric_limits<double>::quiet_NaN();
        return s.g - (-0.02704 + 0.01424 * x - 0.2156 * x * x + 0.01426 * x * x * x);
    }
    if (band == 'R') {
        if (x < -0.5 || x > 4.0) return std::numeric_limits<double>::quiet_NaN();
        return s.g - (-0.02275 + 0.3961 * x - 0.1243 * x * x - 0.01396 * x * x * x + 0.003775 * x * x * x * x);
    }
    throw std::runtime_error(std::string("unknown photometric band '") + band + "' (use V, R or G)");
}

std::vector<CatalogStar> read_catalog(const std::string& path) {
    std::ifstream f(path);
    if (!f) throw std::runtime_error("cannot read catalogue " + path);
    std::vector<CatalogStar> out;
    std::string line;
    while (std::getline(f, line)) {
        if (!line.empty() && line.back() == '\r') line.pop_back();
        if (line.empty() || line[0] == '#') continue;
        const auto fl = split_fields(line);
        if (fl.size() < 3) continue;
        CatalogStar s;
        if (!to_double(fl[0], s.ra_deg) || !to_double(fl[1], s.dec_deg) || !to_double(fl[2], s.g)) continue;  // header, units
        double a = 0, b = 0;
        if (fl.size() >= 5 && to_double(fl[3], a) && to_double(fl[4], b)) {
            s.bp_rp = a - b;
            s.has_colour = true;
        } else if (fl.size() == 4 && to_double(fl[3], a)) {
            s.bp_rp = a;
            s.has_colour = true;
        }
        out.push_back(s);
    }
    return out;
}

std::string gaia_query_url(double ra, double dec, double radius_deg, double g_max) {
    char b[512];
    std::snprintf(b, sizeof b,
                  "https://vizier.cds.unistra.fr/viz-bin/asu-tsv?-source=I/355/gaiadr3"
                  "&-c=%.5f%%20%s%.5f&-c.rd=%.4f&-out=RA_ICRS,DE_ICRS,Gmag,BPmag,RPmag&-out.max=100000&Gmag=%%3C%.1f",
                  ra, dec < 0 ? "-" : "%2B", std::fabs(dec), radius_deg, g_max);
    return b;
}

std::vector<CatalogStar> fetch_gaia(double ra, double dec, double radius_deg, double g_max, const std::string& cache_dir,
                                    std::string* cache_file) {
    const fs::path dir = cache_folder(cache_dir);
    std::error_code ec;
    fs::create_directories(dir, ec);
    char name[160];
    std::snprintf(name, sizeof name, "gaia_%.3f_%+.3f_%.3f_%.1f.tsv", ra, dec, radius_deg, g_max);
    const fs::path file = dir / name;
    if (cache_file) *cache_file = file.string();
    if (fs::is_regular_file(file, ec) && fs::file_size(file, ec) > 200) {
        auto c = read_catalog(file.string());
        if (!c.empty()) return c;
    }
#ifdef _WIN32
    std::string sysdir;
    if (const char* w = std::getenv("SystemRoot")) sysdir = (fs::path(w) / "System32").string();
    const std::string curl = detail::find_program({"curl.exe"}, {sysdir});
#else
    const std::string curl = detail::find_program({"curl"}, {"/usr/bin", "/usr/local/bin"});
#endif
    if (curl.empty()) throw std::runtime_error("curl not found: cannot download Gaia DR3 (use --catalog FILE instead)");
    std::string url = gaia_query_url(ra, dec, radius_deg, g_max);
    const fs::path part = dir / (std::string(name) + ".part"), log = dir / (std::string(name) + ".log");
    std::string why;
    for (const char* host : {"vizier.cds.unistra.fr", "vizier.cfa.harvard.edu"}) {
        const std::string u = std::string("https://") + host + url.substr(url.find("/viz-bin"));
        const int code = detail::run_program(curl, {"-sS", "-f", "-L", "--max-time", "180", "-o", part.string(), u}, 200, log.string());
        if (code == 0) {
            auto c = read_catalog(part.string());
            if (!c.empty()) {
                fs::rename(part, file, ec);
                if (ec) fs::copy_file(part, file, fs::copy_options::overwrite_existing, ec);
                fs::remove(log, ec);
                return c;
            }
            why = "no Gaia stars returned";
        } else {
            std::ifstream lf(log);
            std::getline(lf, why);
            if (why.empty()) why = "curl exit code " + std::to_string(code);
        }
    }
    fs::remove(part, ec);
    throw std::runtime_error("cannot download Gaia DR3 from VizieR (" + why + "); use --catalog FILE with a saved query: " + url);
}

void field_circle(const Prepared& p, double& ra, double& dec, double& radius) {
    if (!p.wcs.valid) throw std::runtime_error("field_circle: no plate solution");
    p.wcs.pix2sky(0.5 * (p.width - 1), 0.5 * (p.height - 1), ra, dec);
    radius = 0.5 * std::hypot(p.width, p.height) * p.wcs.pixel_scale_arcsec() / 3600.0 * 1.02;
}

double psf_flux(const Image& img, double x, double y, double fwhm_px) {
    const double s = fwhm_px * detail::kFwhmToSigma;
    const int R = static_cast<int>(std::ceil(3.0 * s));
    double num = 0.0, den = 0.0;
    for (int v = -R; v <= R; ++v)
        for (int u = -R; u <= R; ++u) {
            const double px = x + u, py = y + v;
            if (px < 0 || py < 0 || px > img.width - 1 || py > img.height - 1) continue;
            const float d = detail::bilinear(img, px, py);
            if (!std::isfinite(d)) continue;
            const double k = std::exp(-0.5 * (u * u + v * v) / (s * s));
            num += d * k;
            den += k * k;
        }
    return den > 0.0 ? num / den * 2.0 * kPiP * s * s : 0.0;
}

double trail_flux_factor(double length_px, double fwhm_px) {
    const double s = fwhm_px * detail::kFwhmToSigma;
    if (length_px <= 1e-6 || s <= 0) return 1.0;
    // a point-source template on a trail of length L: the mean over the trail of
    // exp(-t^2 / (4 s^2)), the overlap of two Gaussians of width s
    return 2.0 * s * std::sqrt(kPiP) / length_px * std::erf(length_px / (4.0 * s));
}

PhotometryResult calibrate_photometry(const Prepared& p, const std::vector<CatalogStar>& catalog, char band) {
    PhotometryResult r;
    r.band = band;
    if (!p.wcs.valid) {
        r.message = "photometry needs a plate solution";
        return r;
    }
    const Image& sky = p.static_sky;
    std::vector<float> smp = detail::sample_pixels(sky, 300000);
    if (smp.size() < 100) {
        r.message = "no static sky";
        return r;
    }
    std::vector<float> tmp = smp;
    const double sig_pix = detail::robust_sigma(tmp);
    const float sky_max = *std::max_element(smp.begin(), smp.end());
    const double fwhm = p.fwhm_px > 0 ? p.fwhm_px : 3.0;
    const double s = fwhm * detail::kFwhmToSigma;
    // noise of the PSF flux estimator: sigma * 2 pi s^2 / sqrt(sum k^2) = sigma * sqrt(4 pi) s
    const double flux_noise = sig_pix * std::sqrt(4.0 * kPiP) * s;

    const auto det = detect_stars(sky, 3000, 10.0);
    const double r_match = std::max(1.5, 0.8 * fwhm);
    const double r_blend = 2.5 * fwhm;

    struct Cal {
        double zp, x;
        bool col_ok;
    };
    std::vector<Cal> cal;
    // project the catalogue once
    struct Proj {
        double x, y;
        const CatalogStar* s;
    };
    std::vector<Proj> proj;
    for (const auto& c : catalog) {
        double x, y;
        if (!p.wcs.sky2pix(c.ra_deg, c.dec_deg, x, y)) continue;
        if (x < -10 || y < -10 || x > p.width + 10 || y > p.height + 10) continue;
        proj.push_back({x, y, &c});
    }
    const int edge = static_cast<int>(std::ceil(3.0 * s)) + 2;
    for (const auto& pr : proj) {
        const CatalogStar& c = *pr.s;
        const double m = catalog_mag(c, band);
        if (!std::isfinite(m)) continue;
        if (pr.x < edge || pr.y < edge || pr.x > p.width - 1 - edge || pr.y > p.height - 1 - edge) continue;
        // a neighbour that adds more than ~10 % of the light: skip the blend
        bool blended = false;
        for (const auto& q : proj)
            if (&q != &pr && std::hypot(q.x - pr.x, q.y - pr.y) < r_blend && q.s->g < c.g + 2.5) {
                blended = true;
                break;
            }
        if (blended) continue;
        const Star* best = nullptr;
        double bd = r_match;
        for (const auto& d : det) {
            const double dd = std::hypot(d.x - pr.x, d.y - pr.y);
            if (dd < bd) {
                bd = dd;
                best = &d;
            }
        }
        if (!best) continue;
        ++r.matched;
        // saturated or non-linear cores: flat tops near the brightest level in the field
        float peak = -1e30f;
        for (int v = -1; v <= 1; ++v)
            for (int u = -1; u <= 1; ++u) {
                const int ix = static_cast<int>(std::lround(best->x)) + u, iy = static_cast<int>(std::lround(best->y)) + v;
                if (sky.inside(ix, iy) && std::isfinite(sky.at(ix, iy))) peak = std::max(peak, sky.at(ix, iy));
            }
        if (peak > 0.6f * sky_max) continue;
        const double f = psf_flux(sky, best->x, best->y, fwhm);
        if (!(f > 20.0 * flux_noise)) continue;
        cal.push_back({m + 2.5 * std::log10(f), c.bp_rp, c.has_colour});
    }
    // asteroid-like (solar) colours first, wider if too few
    std::vector<double> zp;
    for (const auto& range : {std::pair<double, double>{0.5, 1.3}, {0.2, 2.0}, {-1.0, 6.0}}) {
        zp.clear();
        for (const auto& k : cal)
            if (!k.col_ok ? band == 'G' : (k.x >= range.first && k.x <= range.second)) zp.push_back(k.zp);
        if (zp.size() >= 8) break;
    }
    if (zp.size() < 3) {
        r.message = "too few calibration stars (" + std::to_string(zp.size()) + " of " + std::to_string(r.matched) +
                    " matched catalogue stars usable)";
        return r;
    }
    double med = median_of(zp), sd = 0.0;
    for (int it = 0; it < 5; ++it) {
        std::vector<double> dev;
        for (double z : zp) dev.push_back(std::fabs(z - med));
        sd = std::max(0.02, 1.4826 * median_of(dev));
        std::vector<double> keep;
        for (double z : zp)
            if (std::fabs(z - med) <= 3.0 * sd) keep.push_back(z);
        if (keep.size() < 3) break;
        const bool same = keep.size() == zp.size();
        zp = keep;
        med = median_of(zp);
        if (same) break;
    }
    double ss = 0.0;
    for (double z : zp) ss += (z - med) * (z - med);
    r.rms = std::sqrt(ss / std::max<std::size_t>(1, zp.size() - 1));
    r.zero_point = med;
    r.zp_error = 1.253 * r.rms / std::sqrt(static_cast<double>(zp.size()));
    r.stars = static_cast<int>(zp.size());
    r.ok = true;
    char b[256];
    std::snprintf(b, sizeof b, "photometry: %s zero point %.3f +- %.3f from %d Gaia DR3 stars (rms %.3f mag)", std::string(1, band).c_str(),
                  r.zero_point, r.zp_error, r.stars, r.rms);
    r.message = b;
    return r;
}

void apply_photometry(const Prepared& p, const PhotometryResult& cal, SearchResult& res) {
    if (!cal.ok) return;
    res.mag_band = cal.band;
    res.mag_zero_point = cal.zero_point;
    res.mag_zp_error = cal.zp_error;
    res.mag_stars = cal.stars;
    double exp_h = 0.0;
    for (double e : p.exptime_s) exp_h += e;
    exp_h = p.exptime_s.empty() ? 0.0 : exp_h / p.exptime_s.size() / 3600.0;
    const double fwhm = p.fwhm_px > 0 ? p.fwhm_px : res.fwhm_px;
    for (auto& c : res.candidates) {
        if (!(c.flux > 0.0)) {
            c.mag = 99.0;
            continue;
        }
        const double trail = std::hypot(c.vx_pxh, c.vy_pxh) * exp_h;
        const double f = c.flux / trail_flux_factor(trail, fwhm);
        c.mag = cal.zero_point - 2.5 * std::log10(f);
        const double phot = c.snr > 0 ? 1.0857 / c.snr : 1.0;
        c.mag_err = std::sqrt(phot * phot + cal.zp_error * cal.zp_error);
    }
}

}  // namespace fainttrail
