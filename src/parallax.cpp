// SPDX-License-Identifier: GPL-3.0-or-later
// Copyright 2026 Jurgen Kobierczynski
//
// Diurnal parallax model for the near-Earth search (parallax.hpp).
#include "fainttrail/parallax.hpp"

#include <algorithm>
#include <cctype>
#include <cmath>
#include <cstdlib>
#include <sstream>
#include <stdexcept>

namespace fainttrail {

namespace {

constexpr double kPi = 3.14159265358979323846;
constexpr double kDeg = kPi / 180.0;
constexpr double kAuKm = 149597870.7;
constexpr double kEarthA = 6378.137;           // WGS84 equatorial radius, km
constexpr double kEarthF = 1.0 / 298.257223563;

// Degrees from "51.5", "+51 30 00", "51:30:00", "5d18m"; false if unreadable.
bool parse_degrees(std::string s, double& out) {
    for (char& c : s)
        if (!(std::isdigit(static_cast<unsigned char>(c)) || c == '.' || c == '-' || c == '+')) c = ' ';
    std::istringstream is(s);
    std::vector<double> v;
    std::string tok;
    bool neg = false;
    while (is >> tok) {
        if (v.empty() && (tok[0] == '-' || tok[0] == '+')) {
            neg = tok[0] == '-';
            tok.erase(0, 1);
            if (tok.empty()) continue;
        }
        char* end = nullptr;
        const double d = std::strtod(tok.c_str(), &end);
        if (end == tok.c_str() || *end) return false;
        v.push_back(d);
    }
    if (v.empty() || v.size() > 3) return false;
    double x = v[0];
    if (v.size() > 1) x += v[1] / 60.0;
    if (v.size() > 2) x += v[2] / 3600.0;
    out = neg ? -x : x;
    return std::isfinite(out);
}

bool header_angle(const FitsHeader& h, const char* key, double& out) {
    const auto* c = h.find(key);
    return c && parse_degrees(c->value, out);
}

}  // namespace

bool site_from_header(const FitsHeader& h, ObserverSite& site) {
    ObserverSite s;
    const char* lat_keys[] = {"OBSGEO-B", "SITELAT", "LAT-OBS", "OBSLAT"};
    const char* lon_keys[] = {"OBSGEO-L", "SITELONG", "LONG-OBS", "OBSLONG"};
    const char* alt_keys[] = {"OBSGEO-H", "SITEELEV", "ALT-OBS", "OBSALT"};
    for (int k = 0; k < 4 && !s.valid; ++k)
        if (header_angle(h, lat_keys[k], s.lat_deg) && header_angle(h, lon_keys[k], s.lon_deg)) {
            s.valid = true;
            double alt = 0.0;
            if (header_angle(h, alt_keys[k], alt)) s.height_m = alt;
        }
    if (!s.valid || std::fabs(s.lat_deg) > 90.0) return false;
    if (s.lon_deg > 180.0) s.lon_deg -= 360.0;
    site = s;
    return true;
}

ObserverSite parse_site(const std::string& text) {
    std::vector<std::string> parts;
    std::stringstream ss(text);
    std::string item;
    while (std::getline(ss, item, ',')) parts.push_back(item);
    ObserverSite s;
    if (parts.size() < 2 || parts.size() > 3 || !parse_degrees(parts[0], s.lat_deg) || !parse_degrees(parts[1], s.lon_deg) ||
        std::fabs(s.lat_deg) > 90.0)
        throw std::runtime_error("site: give latitude,longitude[,height] in degrees, longitude east positive (e.g. 50.93,5.39,40)");
    if (parts.size() == 3 && !parse_degrees(parts[2], s.height_m)) throw std::runtime_error("site: unreadable height");
    if (s.lon_deg > 180.0) s.lon_deg -= 360.0;
    s.valid = true;
    return s;
}

void observer_position_au(const ObserverSite& site, double mjd_utc, double xyz[3]) {
    const double jd = mjd_utc + 2400000.5;
    const double d = jd - 2451545.0, T = d / 36525.0;
    // Greenwich mean sidereal time (IAU 1982), degrees
    double gmst = 280.46061837 + 360.98564736629 * d + 0.000387933 * T * T - T * T * T / 38710000.0;
    gmst = std::fmod(gmst, 360.0);
    const double theta = (gmst + site.lon_deg) * kDeg;
    const double phi = site.lat_deg * kDeg;
    const double C = 1.0 / std::sqrt(std::cos(phi) * std::cos(phi) + (1 - kEarthF) * (1 - kEarthF) * std::sin(phi) * std::sin(phi));
    const double S = (1 - kEarthF) * (1 - kEarthF) * C;
    const double h = site.height_m / 1000.0;
    const double rc = (kEarthA * C + h) * std::cos(phi), rz = (kEarthA * S + h) * std::sin(phi);
    xyz[0] = rc * std::cos(theta) / kAuKm;
    xyz[1] = rc * std::sin(theta) / kAuKm;
    xyz[2] = rz / kAuKm;
}

void parallax_pixels(const Wcs& wcs, double x, double y, const ObserverSite& site, const std::vector<double>& mjd,
                     std::vector<double>& px, std::vector<double>& py) {
    if (!wcs.valid) throw std::runtime_error("the parallax model needs a plate solution");
    if (!site.valid) throw std::runtime_error("the parallax model needs the observer's site (latitude, longitude)");
    double ra, dec;
    wcs.pix2sky(x, y, ra, dec);
    const double a = ra * kDeg, d = dec * kDeg;
    // east and north unit vectors at the target
    const double e[3] = {-std::sin(a), std::cos(a), 0.0};
    const double n[3] = {-std::sin(d) * std::cos(a), -std::sin(d) * std::sin(a), std::cos(d)};
    // inverse of the CD matrix: intermediate world coordinates (deg) -> pixels
    const double det = wcs.cd11 * wcs.cd22 - wcs.cd12 * wcs.cd21;
    if (det == 0.0) throw std::runtime_error("singular plate solution");
    px.resize(mjd.size());
    py.resize(mjd.size());
    for (std::size_t i = 0; i < mjd.size(); ++i) {
        double r[3];
        observer_position_au(site, mjd[i], r);
        // the object appears displaced opposite to the observer's offset: -r_perp / D
        const double xi = -(r[0] * e[0] + r[1] * e[1] + r[2] * e[2]) / kDeg;   // degrees, east
        const double eta = -(r[0] * n[0] + r[1] * n[1] + r[2] * n[2]) / kDeg;  // degrees, north
        px[i] = (wcs.cd22 * xi - wcs.cd12 * eta) / det;
        py[i] = (-wcs.cd21 * xi + wcs.cd11 * eta) / det;
    }
}

void set_parallax_model(Prepared& p, const ObserverSite& site) {
    std::vector<double> mjd(p.t_h.size());
    for (std::size_t i = 0; i < mjd.size(); ++i) mjd[i] = p.t_ref_mjd + p.t_h[i] / 24.0;
    std::vector<double> px, py;
    parallax_pixels(p.wcs, 0.5 * (p.width - 1), 0.5 * (p.height - 1), site, mjd, px, py);
    // remove the least-squares straight line a + b t: the velocity grid absorbs it
    const double n = static_cast<double>(mjd.size());
    double st = 0, stt = 0;
    for (double t : p.t_h) {
        st += t;
        stt += t * t;
    }
    const double den = n * stt - st * st;
    auto detrend = [&](std::vector<double>& v) {
        double sv = 0, stv = 0;
        for (std::size_t i = 0; i < v.size(); ++i) {
            sv += v[i];
            stv += p.t_h[i] * v[i];
        }
        const double b = den != 0.0 ? (n * stv - st * sv) / den : 0.0, a = (sv - b * st) / n;
        for (std::size_t i = 0; i < v.size(); ++i) v[i] -= a + b * p.t_h[i];
    };
    detrend(px);
    detrend(py);
    p.par_x = std::move(px);
    p.par_y = std::move(py);
}

double parallax_amplitude_px(const Prepared& p) {
    double m = 0.0;
    for (std::size_t i = 0; i < p.par_x.size(); ++i) m = std::max(m, std::hypot(p.par_x[i], p.par_y[i]));
    return m;
}

std::vector<double> inverse_distance_grid(const Prepared& p, const SearchParams& params, double* dq_out) {
    if (dq_out) *dq_out = 0.0;
    if (params.neo_min_dist_au <= 0.0) return {0.0};
    if (static_cast<int>(p.par_x.size()) != p.nframes)
        throw std::runtime_error("the near-Earth search needs the parallax model: a plate solution and the observer's "
                                 "site (set_parallax_model)");
    const double amp = parallax_amplitude_px(p);
    const double qmax = 1.0 / params.neo_min_dist_au;
    const double tol = std::max(0.05, params.grid_tol_px);
    if (amp * qmax <= tol) return {0.0};  // curvature below the grid tolerance: straight tracks suffice
    // neighbouring trials differ by at most 2 tol at the frame of largest displacement,
    // so any distance is within tol of a trial (as for the velocity grid)
    int steps = static_cast<int>(std::ceil(amp * qmax / (2.0 * tol)));
    steps = std::min(steps, std::max(1, params.neo_max_trials - 1));
    std::vector<double> q(static_cast<std::size_t>(steps) + 1);
    for (int k = 0; k <= steps; ++k) q[static_cast<std::size_t>(k)] = qmax * k / steps;
    if (dq_out) *dq_out = qmax / steps;
    return q;
}

}  // namespace fainttrail
