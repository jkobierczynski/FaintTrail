// SPDX-License-Identifier: GPL-3.0-or-later
// Copyright 2026 Jurgen Kobierczynski
//
// Plate solving with ASTAP and sky-position output (solve.hpp).
#include "fainttrail/solve.hpp"

#include <algorithm>
#include <cctype>
#include <chrono>
#include <cmath>
#include <cstdio>
#include <cstdlib>
#include <filesystem>
#include <fstream>
#include <iomanip>
#include <random>
#include <sstream>
#include <stdexcept>
#include <thread>

#include "common.hpp"
#include "process.hpp"


namespace fs = std::filesystem;

namespace fainttrail {

namespace {

constexpr double kDeg = 3.14159265358979323846 / 180.0;

std::string trim(const std::string& s) {
    std::size_t a = 0, b = s.size();
    while (a < b && std::isspace(static_cast<unsigned char>(s[a]))) ++a;
    while (b > a && std::isspace(static_cast<unsigned char>(s[b - 1]))) --b;
    return s.substr(a, b - a);
}

// "11 18 26.6", "11:18:26.6", "11h18m26.6s", "-05 07 30" -> three numbers and the sign.
bool parse_sexagesimal(const std::string& in, double& value) {
    std::string s;
    for (char c : in) s += (std::isdigit(static_cast<unsigned char>(c)) || c == '.' || c == '-' || c == '+') ? c : ' ';
    std::istringstream is(s);
    std::vector<double> v;
    std::string tok;
    bool neg = false;
    while (is >> tok) {
        if (v.empty() && (tok[0] == '-' || tok[0] == '+')) {
            neg = tok[0] == '-';
            tok = tok.substr(1);
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
    value = neg ? -x : x;
    return true;
}

// Header cards from an ASTAP .wcs file (FITS header: 80-character cards, with or
// without line breaks) or .ini file (KEY=VALUE lines).
FitsHeader read_solution_file(const std::string& path) {
    FitsHeader h;
    std::ifstream f(path, std::ios::binary);
    if (!f) return h;
    std::string all((std::istreambuf_iterator<char>(f)), std::istreambuf_iterator<char>());
    std::vector<std::string> lines;
    if (all.find('\n') != std::string::npos) {
        std::istringstream is(all);
        std::string l;
        while (std::getline(is, l)) lines.push_back(l);
    } else {
        for (std::size_t i = 0; i < all.size(); i += 80) lines.push_back(all.substr(i, 80));
    }
    for (std::string l : lines) {
        if (!l.empty() && l.back() == '\r') l.pop_back();
        std::size_t eq = l.find('=');
        if (eq == std::string::npos) continue;
        const std::string key = trim(l.substr(0, eq));
        if (key.empty() || key.size() > 8 || key == "COMMENT" || key == "HISTORY") continue;
        std::string val = l.substr(eq + 1);
        FitsHeader::Card c;
        c.key = key;
        const std::string tv = trim(val);
        if (!tv.empty() && tv[0] == '\'') {
            const std::size_t e = tv.find('\'', 1);
            c.value = trim(tv.substr(1, e == std::string::npos ? std::string::npos : e - 1));
            c.is_string = true;
        } else {
            const std::size_t slash = tv.find('/');
            c.value = trim(slash == std::string::npos ? tv : tv.substr(0, slash));
        }
        h.cards.push_back(c);
    }
    return h;
}

bool is_file(const std::string& p) {
    std::error_code ec;
    return !p.empty() && fs::is_regular_file(p, ec);
}

}  // namespace

bool parse_ra(const std::string& s, double& ra_deg) {
    const std::string t = trim(s);
    const bool sexa = t.find_first_of(" :hHmM") != std::string::npos;
    double v;
    if (!parse_sexagesimal(t, v)) return false;
    ra_deg = sexa ? v * 15.0 : v;
    return ra_deg >= 0 && ra_deg < 360.0 + 1e-9;
}

bool parse_dec(const std::string& s, double& dec_deg) {
    double v;
    if (!parse_sexagesimal(trim(s), v)) return false;
    dec_deg = v;
    return std::fabs(v) <= 90.0;
}

std::string format_ra(double ra_deg, int decimals, char sep) {
    ra_deg = std::fmod(ra_deg, 360.0);
    if (ra_deg < 0) ra_deg += 360.0;
    const double scale = std::pow(10.0, decimals);
    double s = std::round(ra_deg / 15.0 * 3600.0 * scale) / scale;
    if (s >= 86400.0) s -= 86400.0;
    const int h = static_cast<int>(s / 3600.0);
    s -= h * 3600.0;
    const int m = static_cast<int>(s / 60.0);
    s -= m * 60.0;
    char buf[64];
    std::snprintf(buf, sizeof buf, "%02d%c%02d%c%0*.*f", h, sep, m, sep, decimals > 0 ? decimals + 3 : 2, decimals, s);
    return buf;
}

std::string format_dec(double dec_deg, int decimals, char sep) {
    const char sign = dec_deg < 0 ? '-' : '+';
    const double scale = std::pow(10.0, decimals);
    double s = std::round(std::fabs(dec_deg) * 3600.0 * scale) / scale;
    const int d = static_cast<int>(s / 3600.0);
    s -= d * 3600.0;
    const int m = static_cast<int>(s / 60.0);
    s -= m * 60.0;
    char buf[64];
    std::snprintf(buf, sizeof buf, "%c%02d%c%02d%c%0*.*f", sign, d, sep, m, sep, decimals > 0 ? decimals + 3 : 2, decimals, s);
    return buf;
}

SolveHints hints_from_header(const FitsHeader& h, int image_height) {
    SolveHints out;
    double ra = 0, dec = 0;
    auto num = [&](const char* k, double& v) {
        const auto* c = h.find(k);
        if (!c) return false;
        char* end = nullptr;
        v = std::strtod(c->value.c_str(), &end);
        return end != c->value.c_str();
    };
    if (h.has("OBJCTRA") && h.has("OBJCTDEC") && parse_sexagesimal(h.get_string("OBJCTRA"), ra) &&
        parse_dec(h.get_string("OBJCTDEC"), dec)) {
        out.has_center = true;
        out.ra_deg = ra * 15.0;
        out.dec_deg = dec;
    } else if (h.has("RA") && h.has("DEC")) {
        const auto* rc = h.find("RA");
        const auto* dc = h.find("DEC");
        if (rc->is_string) {  // some programs write sexagesimal strings
            if (parse_ra(rc->value, ra) && parse_dec(dc->value, dec)) {
                out.has_center = true;
                out.ra_deg = ra;
                out.dec_deg = dec;
            }
        } else if (num("RA", ra) && num("DEC", dec)) {
            out.has_center = true;
            out.ra_deg = ra;
            out.dec_deg = dec;
        }
    } else if (num("CRVAL1", ra) && num("CRVAL2", dec)) {
        out.has_center = true;
        out.ra_deg = ra;
        out.dec_deg = dec;
    }
    double focal = 0, pix = 0;
    if (num("FOCALLEN", focal) && focal > 0 && (num("XPIXSZ", pix) || num("PIXSIZE1", pix)) && pix > 0) {
        // XPIXSZ includes the binning in most capture programs (MaxIm, NINA, Ekos)
        const double scale = 206.265 * pix / focal;  // arcsec per pixel
        out.fov_deg = scale * image_height / 3600.0;
    }
    return out;
}

std::string find_astap(const std::string& hint) {
    if (!hint.empty()) return is_file(hint) ? hint : std::string();
    if (const char* e = std::getenv("FAINTTRAIL_ASTAP"))
        if (is_file(e)) return e;
#ifdef _WIN32
    // the command-line build first: it waits and writes to the console
    return detail::find_program({"astap_cli.exe", "astap.exe"},
                                {"C:\\Program Files\\astap", "C:\\Program Files (x86)\\astap", "C:\\astap"});
#else
    return detail::find_program({"astap_cli", "astap"},
                                {"/usr/bin", "/usr/local/bin", "/opt/astap", "/Applications/ASTAP.app/Contents/MacOS"});
#endif
}

SolveResult solve_with_astap(const Image& image, const SolveHints& hints, const SolveOptions& opt) {
    SolveResult r;
    r.solver = find_astap(opt.astap);
    if (r.solver.empty()) {
        r.message = opt.astap.empty()
                        ? "ASTAP not found (install it from https://www.hnsky.org/astap.htm with a star database, "
                          "or give its path with --astap / FAINTTRAIL_ASTAP)"
                        : "ASTAP not found at " + opt.astap;
        return r;
    }
    // stretch to 0..65535 so that any ASTAP version reads it like a camera frame
    std::vector<float> s = detail::sample_pixels(image, 400000);
    if (s.size() < 100) {
        r.message = "image has no usable pixels";
        return r;
    }
    std::sort(s.begin(), s.end());
    const double lo = s[static_cast<std::size_t>(0.001 * (s.size() - 1))];
    const double med = s[s.size() / 2];
    const double hi = s[static_cast<std::size_t>(0.9995 * (s.size() - 1))];
    const double span = std::max(1e-12, hi - lo);
    Image out(image.width, image.height);
    for (std::size_t i = 0; i < image.data.size(); ++i) {
        const double v = std::isfinite(image.data[i]) ? image.data[i] : med;
        out.data[i] = static_cast<float>(std::clamp(1000.0 + (v - lo) / span * 60000.0, 0.0, 65535.0));
    }
    std::error_code ec;
    const fs::path dir = opt.work_dir.empty() ? fs::temp_directory_path(ec) : fs::path(opt.work_dir);
    std::mt19937_64 rng(static_cast<unsigned long long>(std::chrono::steady_clock::now().time_since_epoch().count()));
    const std::string stem = "fainttrail_solve_" + std::to_string(rng() % 1000000000ULL);
    const fs::path fits = dir / (stem + ".fits"), base = dir / stem, log = dir / (stem + ".log");
    FitsHeader hdr;
    write_fits(fits.string(), out, hdr, -32);

    std::vector<std::string> args = {"-f", fits.string(), "-o", base.string(), "-wcs"};
    double radius = hints.search_radius_deg;
    if (hints.has_center) {
        char b[64];
        std::snprintf(b, sizeof b, "%.6f", hints.ra_deg / 15.0);
        args.insert(args.end(), {"-ra", b});
        std::snprintf(b, sizeof b, "%.6f", hints.dec_deg + 90.0);  // south pole distance
        args.insert(args.end(), {"-spd", b});
        if (radius <= 0) radius = 30.0;
    } else if (radius <= 0) {
        radius = 180.0;  // blind
    }
    {
        char b[64];
        std::snprintf(b, sizeof b, "%.1f", radius);
        args.insert(args.end(), {"-r", b});
        std::snprintf(b, sizeof b, "%.4f", hints.fov_deg);  // 0 = ASTAP tries itself
        args.insert(args.end(), {"-fov", b});
    }
    if (opt.downsample > 0) args.insert(args.end(), {"-z", std::to_string(opt.downsample)});
    if (!opt.database_dir.empty()) args.insert(args.end(), {"-d", opt.database_dir});

    const int code = detail::run_program(r.solver, args, std::max(10, opt.timeout_s), log.string());
    const FitsHeader ini = read_solution_file(base.string() + ".ini");
    FitsHeader sol = read_solution_file(base.string() + ".wcs");
    if (!sol.has("CRVAL1")) sol = ini;
    if (!sol.has("CTYPE1")) {
        sol.set("CTYPE1", std::string("RA---TAN"));
        sol.set("CTYPE2", std::string("DEC--TAN"));
    }
    const std::string solved = ini.get_string("PLTSOLVD", sol.get_string("PLTSOLVD", ""));
    r.wcs = wcs_from_header(sol);
    r.ok = r.wcs.valid && (solved.empty() || solved[0] == 'T');
    if (r.ok) {
        r.scale_arcsec = r.wcs.pixel_scale_arcsec();
        r.wcs.pix2sky(0.5 * (image.width - 1), 0.5 * (image.height - 1), r.center_ra_deg, r.center_dec_deg);
        r.rotation_deg = std::atan2(r.wcs.cd21, r.wcs.cd22) / kDeg;
        std::ostringstream m;
        m.setf(std::ios::fixed);
        m << "solved by " << fs::path(r.solver).filename().string() << ": centre " << format_ra(r.center_ra_deg, 1) << ' '
          << format_dec(r.center_dec_deg, 0) << ", " << std::setprecision(3) << r.scale_arcsec << "\"/px, rotation "
          << std::setprecision(1) << r.rotation_deg << " deg";
        r.message = m.str();
    } else {
        std::string why = ini.get_string("ERROR", ini.get_string("WARNING", ""));
        if (code == -1) why = "could not start " + r.solver;
        else if (code == -2) why = "no solution within " + std::to_string(opt.timeout_s) + " s";
        else if (why.empty()) why = "no solution (ASTAP exit code " + std::to_string(code) + ")";
        r.message = "plate solving failed: " + why +
                    (hints.has_center ? "" : "; a pointing in the FITS header (RA/DEC or OBJCTRA/OBJCTDEC) or --ra/--dec helps") +
                    (hints.fov_deg > 0 ? "" : "; so does the field height (FOCALLEN and XPIXSZ, or --fov)");
    }
    if (!opt.keep_files)
        for (const char* ext : {".fits", ".ini", ".wcs", ".log"}) fs::remove(base.string() + ext, ec);
    return r;
}

void update_sky_coordinates(const Prepared& p, SearchResult& res) {
    if (!p.wcs.valid) return;
    res.has_sky = true;
    for (auto& c : res.candidates) {
        double ra, dec, ra2, dec2;
        p.wcs.pix2sky(c.x, c.y, ra, dec);
        const double dt = 0.1;
        p.wcs.pix2sky(c.x + c.vx_pxh * dt, c.y + c.vy_pxh * dt, ra2, dec2);
        c.ra_deg = ra;
        c.dec_deg = dec;
        double dra = ra2 - ra;
        if (dra > 180.0) dra -= 360.0;
        if (dra < -180.0) dra += 360.0;
        const double east = dra * std::cos(dec * kDeg) * 3600.0, north = (dec2 - dec) * 3600.0;
        c.rate_arcsec_min = std::hypot(east, north) / (dt * 60.0);
        double pa = std::atan2(east, north) / kDeg;
        if (pa < 0) pa += 360.0;
        c.pa_deg = pa;
    }
}

SolveResult solve_prepared(Prepared& p, const SolveHints& hints, const SolveOptions& options) {
    if (p.static_sky.data.empty()) throw std::runtime_error("solve_prepared: no static sky (prepare() first)");
    SolveResult r = solve_with_astap(p.static_sky, hints, options);
    if (r.ok) p.wcs = r.wcs;
    return r;
}

std::vector<std::string> mpc_observations(const Prepared& p, const SearchResult& res, const std::string& obscode,
                                          int per_object, bool include_flagged) {
    std::vector<std::string> out;
    if (!p.wcs.valid || p.t_h.empty()) return out;
    std::vector<double> t = p.t_h;
    std::sort(t.begin(), t.end());
    std::vector<double> when;
    const int k = std::max(1, per_object);
    for (int i = 0; i < k; ++i) {
        const std::size_t idx = k == 1 ? t.size() / 2 : static_cast<std::size_t>(std::lround(i * (t.size() - 1.0) / (k - 1)));
        when.push_back(t[idx]);
    }
    std::string code = obscode.substr(0, 3);
    while (code.size() < 3) code += ' ';
    for (const auto& c : res.candidates) {
        if (!c.ok && !include_flagged) continue;
        char desig[16];
        std::snprintf(desig, sizeof desig, "FT%05d", c.id % 100000);
        for (const double th : when) {
            double ra, dec;
            if (!p.wcs.pix2sky(c.x + c.vx_pxh * th, c.y + c.vy_pxh * th, ra, dec)) continue;
            // date: UTC, fraction of the day to 5 decimals (0.9 s)
            const double mjd = p.t_ref_mjd + th / 24.0;
            const std::string iso = iso_from_mjd(mjd, 3);  // "YYYY-MM-DDThh:mm:ss.sss"
            const int Y = std::stoi(iso.substr(0, 4)), M = std::stoi(iso.substr(5, 2)), D = std::stoi(iso.substr(8, 2));
            const double frac = mjd - std::floor(mjd);
            char date[32];
            std::snprintf(date, sizeof date, "%04d %02d %08.5f", Y, M, D + frac);
            std::string line(80, ' ');
            auto put = [&](int col, const std::string& s) {  // 1-based column
                for (std::size_t i = 0; i < s.size() && col - 1 + i < line.size(); ++i) line[col - 1 + i] = s[i];
            };
            put(6, desig);
            put(15, "C");
            put(16, date);
            put(33, format_ra(ra, 2, ' '));
            put(45, format_dec(dec, 1, ' '));
            if (res.mag_band && c.mag < 90.0) {  // calibrated (photometry.hpp): columns 66-70 magnitude, 71 band
                char mg[16];
                std::snprintf(mg, sizeof mg, "%4.1f", c.mag);
                put(66, mg);
                put(71, std::string(1, res.mag_band));
            }
            put(78, code);
            out.push_back(line);
        }
    }
    return out;
}

}  // namespace fainttrail
