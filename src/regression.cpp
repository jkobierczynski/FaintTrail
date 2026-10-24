// SPDX-License-Identifier: GPL-3.0-or-later
// Copyright 2026 Jurgen Kobierczynski
//
// Real-data regression cases: parsing, running and checking (see regression.hpp).
#include "fainttrail/regression.hpp"

#include <algorithm>
#include <chrono>
#include <cmath>
#include <cstdarg>
#include <cstdio>
#include <cstdlib>
#include <filesystem>
#include <fstream>
#include <sstream>
#include <stdexcept>

#include "fainttrail/astrometry.hpp"
#include "fainttrail/fits.hpp"
#include "fainttrail/io.hpp"

namespace fs = std::filesystem;

namespace fainttrail {

namespace {

std::string trim(const std::string& s) {
    const auto a = s.find_first_not_of(" \t\r\n");
    if (a == std::string::npos) return {};
    const auto b = s.find_last_not_of(" \t\r\n");
    return s.substr(a, b - a + 1);
}

std::string expand_env(const std::string& s) {
    std::string out;
    for (std::size_t i = 0; i < s.size(); ++i) {
        if (s[i] == '$' && i + 1 < s.size() && s[i + 1] == '{') {
            const auto end = s.find('}', i + 2);
            if (end == std::string::npos) throw std::runtime_error("unterminated ${ in '" + s + "'");
            const std::string var = s.substr(i + 2, end - i - 2);
            const char* v = std::getenv(var.c_str());
            if (!v || !*v) return {};  // unset variable: the path is unavailable
            out += v;
            i = end;
        } else {
            out += s[i];
        }
    }
    return out;
}

std::string resolve(const fs::path& base, const std::string& value) {
    const std::string e = expand_env(value);
    if (e.empty()) return {};
    fs::path p(e);
    if (p.is_relative()) p = base / p;
    return p.lexically_normal().string();
}

double to_num(const std::string& key, const std::string& v, const std::string& where) {
    try {
        std::size_t used = 0;
        const double d = std::stod(v, &used);
        if (trim(v.substr(used)).empty()) return d;
    } catch (...) {
    }
    throw std::runtime_error(where + ": " + key + " expects a number, got '" + v + "'");
}

double seconds_since(std::chrono::steady_clock::time_point t0) {
    return std::chrono::duration<double>(std::chrono::steady_clock::now() - t0).count();
}

std::string fmt(const char* f, ...) {
    char buf[512];
    va_list ap;
    va_start(ap, f);
    std::vsnprintf(buf, sizeof buf, f, ap);
    va_end(ap);
    return buf;
}

// angular distance in arcsec
double sep_arcsec(double ra1, double dec1, double ra2, double dec2) {
    const double d2r = 3.14159265358979323846 / 180.0;
    const double s = std::sin(dec1 * d2r) * std::sin(dec2 * d2r) +
                     std::cos(dec1 * d2r) * std::cos(dec2 * d2r) * std::cos((ra1 - ra2) * d2r);
    return std::acos(std::clamp(s, -1.0, 1.0)) / d2r * 3600.0;
}

double read_baseline_threshold(const std::string& path) {
    std::ifstream in(path);
    std::string line;
    while (std::getline(in, line)) {
        if (line.empty() || line[0] != '#') break;
        const auto p = line.find("threshold_snr=");
        if (p != std::string::npos) return std::atof(line.c_str() + p + 14);
    }
    return 0.0;
}

}  // namespace

const char* status_name(RegressionStatus s) {
    switch (s) {
        case RegressionStatus::Pass: return "PASS";
        case RegressionStatus::Changed: return "CHANGED";
        case RegressionStatus::Fail: return "FAIL";
        case RegressionStatus::Skipped: return "SKIP";
    }
    return "?";
}

RegressionCase load_regression_case(const std::string& path) {
    std::ifstream in(path);
    if (!in) throw std::runtime_error("cannot open regression case " + path);
    const fs::path base = fs::path(path).parent_path();
    RegressionCase c;
    c.case_path = path;
    c.name = fs::path(path).stem().string();
    std::string section;
    RegressionObject* obj = nullptr;
    std::string line;
    int lineno = 0;
    bool have_data = false;
    while (std::getline(in, line)) {
        ++lineno;
        const std::string where = path + ":" + std::to_string(lineno);
        if (const auto h = line.find('#'); h != std::string::npos) line.resize(h);
        line = trim(line);
        if (line.empty()) continue;
        if (line.front() == '[') {
            if (line.back() != ']') throw std::runtime_error(where + ": malformed section header");
            const std::string inner = trim(line.substr(1, line.size() - 2));
            if (inner.rfind("object", 0) == 0) {
                section = "object";
                c.objects.emplace_back();
                obj = &c.objects.back();
                obj->name = trim(inner.substr(6));
                if (obj->name.empty()) obj->name = "object " + std::to_string(c.objects.size());
            } else if (inner == "search" || inner == "expect") {
                section = inner;
                obj = nullptr;
            } else {
                throw std::runtime_error(where + ": unknown section [" + inner + "]");
            }
            continue;
        }
        const auto eq = line.find('=');
        if (eq == std::string::npos) throw std::runtime_error(where + ": expected key = value");
        const std::string key = trim(line.substr(0, eq)), val = trim(line.substr(eq + 1));
        auto num = [&] { return to_num(key, val, where); };
        if (section.empty()) {
            if (key == "name") c.name = val;
            else if (key == "description") c.description = val;
            else if (key == "data") { c.data_dir = resolve(base, val); have_data = true; }
            else if (key == "baseline") c.baseline_path = resolve(base, val);
            else throw std::runtime_error(where + ": unknown key '" + key + "'");
        } else if (section == "search") {
            SearchParams& p = c.params;
            if (key == "vmin") p.vmin_pxh = num();
            else if (key == "vmax") p.vmax_pxh = num();
            else if (key == "threshold") p.threshold_snr = num();
            else if (key == "null_runs") p.null_runs = static_cast<int>(num());
            else if (key == "false_per_run") p.target_false_per_run = num();
            else if (key == "fwhm") p.fwhm_px = num();
            else if (key == "grid_tol") p.grid_tol_px = num();
            else if (key == "static_exclude") p.static_exclude = static_cast<int>(num());
            else if (key == "mask_scatter") p.mask_scatter = num();
            else if (key == "weight_scatter") p.weight_scatter = num();
            else if (key == "seed") p.seed = static_cast<std::uint64_t>(num());
            else if (key == "wcs_from") c.wcs_from = resolve(base, val);
            else if (key == "wcs_scale") c.wcs_scale = num();
            else throw std::runtime_error(where + ": unknown search key '" + key + "'");
        } else if (section == "expect") {
            if (key == "max_unexpected") c.max_unexpected = static_cast<int>(num());
            else if (key == "threshold_min") c.threshold_min = num();
            else if (key == "threshold_max") c.threshold_max = num();
            else throw std::runtime_error(where + ": unknown expect key '" + key + "'");
        } else {
            RegressionObject& o = *obj;
            if (key == "mjd") o.mjd = num();
            else if (key == "x") o.x = num();
            else if (key == "y") o.y = num();
            else if (key == "vx") o.vx = num();
            else if (key == "vy") o.vy = num();
            else if (key == "pos_tol") o.pos_tol = num();
            else if (key == "vel_tol") o.vel_tol = num();
            else if (key == "min_snr") o.min_snr = num();
            else if (key == "ra") { o.ra = num(); o.has_sky = true; }
            else if (key == "dec") { o.dec = num(); o.has_sky = true; }
            else if (key == "sky_tol") o.sky_tol_arcsec = num();
            else if (key == "source") o.source = val;
            else throw std::runtime_error(where + ": unknown object key '" + key + "'");
        }
    }
    if (!have_data) throw std::runtime_error(path + ": no 'data =' given");
    for (const auto& o : c.objects)
        if (o.mjd <= 0.0) throw std::runtime_error(path + ": object '" + o.name + "' has no mjd");
    return c;
}

RegressionOutcome run_regression_case(const RegressionCase& c, const RegressionOptions& opt,
                                      const ProgressFn& progress) {
    RegressionOutcome out;
    if (c.data_dir.empty() || !fs::is_directory(c.data_dir)) {
        out.status = RegressionStatus::Skipped;
        out.skip_reason = c.data_dir.empty() ? "data location not set (environment variable missing)"
                                             : "data directory not found: " + c.data_dir;
        return out;
    }
    if (!c.wcs_from.empty() && !fs::exists(c.wcs_from)) {
        out.status = RegressionStatus::Skipped;
        out.skip_reason = "plate solution not found: " + c.wcs_from;
        return out;
    }

    // ---- load
    auto t0 = std::chrono::steady_clock::now();
    std::vector<std::string> files;
    for (const auto& e : fs::directory_iterator(c.data_dir)) {
        if (!e.is_regular_file()) continue;
        std::string ext = e.path().extension().string();
        std::transform(ext.begin(), ext.end(), ext.begin(), [](unsigned char ch) { return std::tolower(ch); });
        if (ext == ".fits" || ext == ".fit" || ext == ".fts") files.push_back(e.path().string());
    }
    std::sort(files.begin(), files.end());
    if (files.empty()) {
        out.status = RegressionStatus::Skipped;
        out.skip_reason = "no FITS files in " + c.data_dir;
        return out;
    }
    auto frames = load_frames(files);
    auto add = [&](bool ok, bool hard, std::string text) { out.checks.push_back({ok, hard, std::move(text)}); };
    if (!c.wcs_from.empty()) {
        const auto rep = apply_wcs_from_file(c.wcs_from, frames, c.wcs_scale);
        add(rep.ok, true, fmt("plate solution transferred: %d stars, rms %.2f px", rep.matched, rep.rms_px));
    }
    out.load_s = seconds_since(t0);
    out.nframes = frames.size();
    out.width = frames.front().image.width;
    out.height = frames.front().image.height;

    // ---- search
    SearchParams p = c.params;
    if (!opt.backend.empty()) p.backend = opt.backend;
    if (opt.threads > 0) p.threads = opt.threads;
    ProgressFn pp, sp;
    if (progress) {
        pp = [&](double f, const std::string& s) { return progress(0.1 * f, s); };
        sp = [&](double f, const std::string& s) { return progress(0.1 + 0.9 * f, s); };
    }
    t0 = std::chrono::steady_clock::now();
    const Prepared prep = prepare(frames, p, pp);
    out.prepare_s = seconds_since(t0);
    t0 = std::chrono::steady_clock::now();
    out.result = search(prep, p, sp);
    out.search_s = seconds_since(t0);
    out.profile = profile_text(prep, out.result, out.load_s);
    const SearchResult& r = out.result;

    // ---- known objects
    std::vector<bool> explained(r.candidates.size(), false);
    for (const auto& o : c.objects) {
        const double dt_h = (r.t_ref_mjd - o.mjd) * 24.0;
        const double ex = o.x - 1.0 + o.vx * dt_h, ey = o.y - 1.0 + o.vy * dt_h;  // 0-based at t_ref
        int best = -1;
        double best_d = 1e300;
        for (std::size_t k = 0; k < r.candidates.size(); ++k) {
            const auto& cd = r.candidates[k];
            const double d = std::hypot(cd.x - ex, cd.y - ey);
            const double dv = std::hypot(cd.vx_pxh - o.vx, cd.vy_pxh - o.vy);
            if (d <= o.pos_tol && dv <= o.vel_tol && d < best_d) {
                best = static_cast<int>(k);
                best_d = d;
            }
        }
        if (best < 0) {
            add(false, true, fmt("%s: not found within %.1f px / %.1f px/h of (%.2f, %.2f)", o.name.c_str(), o.pos_tol,
                                 o.vel_tol, ex + 1.0, ey + 1.0));
            continue;
        }
        explained[static_cast<std::size_t>(best)] = true;
        const auto& cd = r.candidates[static_cast<std::size_t>(best)];
        const double dv = std::hypot(cd.vx_pxh - o.vx, cd.vy_pxh - o.vy);
        add(cd.ok, true, fmt("%s: candidate %d, S/N %.2f, offset %.2f px, %.2f px/h%s%s", o.name.c_str(), cd.id, cd.snr,
                             best_d, dv, cd.ok ? "" : ", REJECTED: ", cd.ok ? "" : cd.flags.c_str()));
        if (o.min_snr > 0.0) add(cd.snr >= o.min_snr, true, fmt("%s: S/N %.2f >= %.2f", o.name.c_str(), cd.snr, o.min_snr));
        if (o.has_sky) {
            if (!prep.wcs.valid) {
                add(false, true, o.name + ": sky position given but the frames have no WCS");
            } else {
                // compare at the object's epoch
                double ra = 0, dec = 0;
                const double dth = (o.mjd - r.t_ref_mjd) * 24.0;
                prep.wcs.pix2sky(cd.x + cd.vx_pxh * dth, cd.y + cd.vy_pxh * dth, ra, dec);
                const double s = sep_arcsec(ra, dec, o.ra, o.dec);
                add(s <= o.sky_tol_arcsec, true,
                    fmt("%s: sky position off by %.1f\" (limit %.1f\")", o.name.c_str(), s, o.sky_tol_arcsec));
            }
        }
    }

    // ---- everything else
    int unexpected = 0;
    std::string list;
    for (std::size_t k = 0; k < r.candidates.size(); ++k) {
        const auto& cd = r.candidates[k];
        if (!cd.ok || explained[k]) continue;
        ++unexpected;
        list += fmt(" #%d (%.1f, %.1f) S/N %.1f;", cd.id, cd.x + 1.0, cd.y + 1.0, cd.snr);
    }
    add(unexpected <= c.max_unexpected, true,
        fmt("unexplained passing candidates: %d (limit %d)%s", unexpected, c.max_unexpected, list.c_str()));
    if (c.threshold_min > 0.0 || c.threshold_max > 0.0) {
        const bool ok = (c.threshold_min <= 0.0 || r.threshold_snr >= c.threshold_min) &&
                        (c.threshold_max <= 0.0 || r.threshold_snr <= c.threshold_max);
        add(ok, true, fmt("threshold S/N %.2f within [%.2f, %.2f]", r.threshold_snr, c.threshold_min, c.threshold_max));
    }

    // ---- baseline
    if (!c.baseline_path.empty()) {
        if (!fs::exists(c.baseline_path)) {
            add(true, false, "no baseline file yet (" + c.baseline_path + "); create one with --write-baseline");
        } else {
            const auto base = read_candidates_csv(c.baseline_path);
            const double bthr = read_baseline_threshold(c.baseline_path);
            if (bthr > 0.0) {
                const double rel = std::fabs(r.threshold_snr - bthr) / bthr;
                add(rel <= 0.02, false, fmt("threshold %.3f vs baseline %.3f (%+.1f%%)", r.threshold_snr, bthr,
                                            100.0 * (r.threshold_snr - bthr) / bthr));
            }
            const double vtol = std::max(1.0, r.dv_pxh);
            auto same = [&](const Candidate& a, const Candidate& b) {
                return std::hypot(a.x - b.x, a.y - b.y) <= 0.5 && std::hypot(a.vx_pxh - b.vx_pxh, a.vy_pxh - b.vy_pxh) <= vtol;
            };
            int lost = 0, gained = 0, drift = 0, verdict = 0;
            std::string detail;
            std::vector<bool> used(r.candidates.size(), false);
            for (const auto& b : base) {
                int m = -1;
                for (std::size_t k = 0; k < r.candidates.size(); ++k)
                    if (!used[k] && same(b, r.candidates[k])) {
                        m = static_cast<int>(k);
                        break;
                    }
                if (m < 0) {
                    if (b.ok) {
                        ++lost;
                        detail += fmt(" lost passing (%.1f, %.1f) S/N %.1f;", b.x + 1, b.y + 1, b.snr);
                    }
                    continue;
                }
                used[static_cast<std::size_t>(m)] = true;
                const auto& n = r.candidates[static_cast<std::size_t>(m)];
                if (n.ok != b.ok) {
                    ++verdict;
                    detail += fmt(" (%.1f, %.1f) %s -> %s;", b.x + 1, b.y + 1, b.ok ? "passed" : "rejected",
                                  n.ok ? "passed" : "rejected");
                }
                if (b.ok && std::fabs(n.snr - b.snr) > 0.02 * b.snr) {
                    ++drift;
                    detail += fmt(" (%.1f, %.1f) S/N %.2f -> %.2f;", b.x + 1, b.y + 1, b.snr, n.snr);
                }
            }
            for (std::size_t k = 0; k < r.candidates.size(); ++k)
                if (!used[k] && r.candidates[k].ok) {
                    ++gained;
                    detail += fmt(" new passing (%.1f, %.1f) S/N %.1f;", r.candidates[k].x + 1, r.candidates[k].y + 1,
                                  r.candidates[k].snr);
                }
            const bool same_all = lost + gained + drift + verdict == 0;
            add(same_all, false,
                fmt("baseline: %zu candidates before, %zu now; passing lost %d, new %d, verdict changes %d, S/N drift %d%s",
                    base.size(), r.candidates.size(), lost, gained, verdict, drift, detail.c_str()));
        }
    }

    out.status = RegressionStatus::Pass;
    for (const auto& ch : out.checks) {
        if (ch.ok) continue;
        if (ch.hard || opt.strict) {
            out.status = RegressionStatus::Fail;
            break;
        }
        out.status = RegressionStatus::Changed;
    }
    return out;
}

void write_regression_baseline(const std::string& path, const RegressionCase&, const RegressionOutcome& o) {
    const fs::path parent = fs::path(path).parent_path();
    if (!parent.empty()) fs::create_directories(parent);
    write_candidates_csv(path, o.result);
}

}  // namespace fainttrail
