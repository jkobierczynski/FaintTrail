// SPDX-License-Identifier: GPL-3.0-or-later
// Copyright 2026 Jurgen Kobierczynski
#include "fainttrail/io.hpp"

#include <cstdio>
#include <fstream>
#include <iomanip>
#include <map>
#include <sstream>
#include <stdexcept>

#include "common.hpp"
#include "fainttrail/fits.hpp"

namespace fainttrail {

namespace {

std::ofstream open_out(const std::string& path) {
    std::ofstream f(path);
    if (!f) throw std::runtime_error("cannot write " + path);
    f.setf(std::ios::fixed);
    return f;
}

std::vector<std::string> split_csv(const std::string& line) {
    std::vector<std::string> out;
    std::string cur;
    bool quoted = false;
    for (char ch : line) {
        if (ch == '"') quoted = !quoted;
        else if (ch == ',' && !quoted) {
            out.push_back(cur);
            cur.clear();
        } else if (ch != '\r') cur += ch;
    }
    out.push_back(cur);
    return out;
}

struct CsvTable {
    std::map<std::string, std::size_t> col;
    std::vector<std::vector<std::string>> rows;
    double num(std::size_t r, const std::string& name, double fallback = 0.0) const {
        auto it = col.find(name);
        if (it == col.end() || it->second >= rows[r].size()) return fallback;
        try {
            return std::stod(rows[r][it->second]);
        } catch (...) {
            return fallback;
        }
    }
    std::string str(std::size_t r, const std::string& name) const {
        auto it = col.find(name);
        if (it == col.end() || it->second >= rows[r].size()) return {};
        return rows[r][it->second];
    }
};

CsvTable read_csv(const std::string& path) {
    std::ifstream f(path);
    if (!f) throw std::runtime_error("cannot read " + path);
    CsvTable t;
    std::string line;
    bool header = true;
    while (std::getline(f, line)) {
        if (line.empty() || line[0] == '#') continue;
        auto cells = split_csv(line);
        if (header) {
            for (std::size_t i = 0; i < cells.size(); ++i) t.col[cells[i]] = i;
            header = false;
        } else {
            t.rows.push_back(cells);
        }
    }
    return t;
}

}  // namespace

void write_candidates_csv(const std::string& path, const SearchResult& res) {
    auto f = open_out(path);
    f << "# fainttrail candidates; pixel coordinates are FITS 1-based at the reference epoch\n";
    f << "# t_ref_mjd=" << std::setprecision(8) << res.t_ref_mjd << " t_ref_utc=" << iso_from_mjd(res.t_ref_mjd, 2)
      << " threshold_snr=" << std::setprecision(3) << res.threshold_snr << " method=" << res.threshold_method
      << " velocities=" << res.n_velocities << " dv_pxh=" << res.dv_pxh << " fwhm_px=" << res.fwhm_px << "\n";
    f << "id,ok,flags,x,y,pos_err_px,vx_pxh,vy_pxh,vel_err_pxh,snr,snr_coarse,snr_half1,snr_half2,max_frame_fraction,"
         "false_alarms_per_run,flux,mag_inst,ra_deg,dec_deg,rate_arcsec_min,pa_deg\n";
    for (const auto& c : res.candidates) {
        f << c.id << ',' << (c.ok ? 1 : 0) << ",\"" << c.flags << "\"," << std::setprecision(3) << c.x + 1.0 << ','
          << c.y + 1.0 << ',' << c.pos_err_px << ',' << c.vx_pxh << ',' << c.vy_pxh << ',' << c.vel_err_pxh << ','
          << std::setprecision(2) << c.snr << ','
          << c.snr_coarse << ',' << c.snr_half1 << ',' << c.snr_half2 << ',' << std::setprecision(3)
          << c.max_frame_fraction << ',' << std::scientific << std::setprecision(2) << c.false_alarms_per_run
          << std::fixed << ',' << std::setprecision(2) << c.flux << ',' << std::setprecision(3) << c.mag_inst << ','
          << std::setprecision(7) << c.ra_deg << ',' << c.dec_deg << ',' << std::setprecision(4)
          << c.rate_arcsec_min << ',' << std::setprecision(2) << c.pa_deg << '\n';
    }
}

std::vector<Candidate> read_candidates_csv(const std::string& path) {
    const CsvTable t = read_csv(path);
    std::vector<Candidate> out;
    for (std::size_t r = 0; r < t.rows.size(); ++r) {
        Candidate c;
        c.id = static_cast<int>(t.num(r, "id"));
        c.ok = t.num(r, "ok") != 0.0;
        c.flags = t.str(r, "flags");
        c.x = t.num(r, "x") - 1.0;
        c.y = t.num(r, "y") - 1.0;
        c.pos_err_px = t.num(r, "pos_err_px");
        c.vel_err_pxh = t.num(r, "vel_err_pxh");
        c.vx_pxh = t.num(r, "vx_pxh");
        c.vy_pxh = t.num(r, "vy_pxh");
        c.snr = t.num(r, "snr");
        c.snr_coarse = t.num(r, "snr_coarse");
        c.snr_half1 = t.num(r, "snr_half1");
        c.snr_half2 = t.num(r, "snr_half2");
        c.max_frame_fraction = t.num(r, "max_frame_fraction");
        c.false_alarms_per_run = t.num(r, "false_alarms_per_run");
        c.flux = t.num(r, "flux");
        c.mag_inst = t.num(r, "mag_inst");
        c.ra_deg = t.num(r, "ra_deg");
        c.dec_deg = t.num(r, "dec_deg");
        c.rate_arcsec_min = t.num(r, "rate_arcsec_min");
        c.pa_deg = t.num(r, "pa_deg");
        out.push_back(c);
    }
    return out;
}

void write_tracklets_csv(const std::string& path, const Prepared& prep, const SearchResult& res, bool ok_only) {
    auto f = open_out(path);
    f << "# fainttrail tracklets: predicted position of each candidate in each frame (FITS 1-based pixels)\n";
    f << "id,frame,mjd_utc,utc,x,y,ra_deg,dec_deg,mag_inst\n";
    for (const auto& c : res.candidates) {
        if (ok_only && !c.ok) continue;
        for (const auto& p : track_points(prep, c)) {
            f << c.id << ',' << p.frame + 1 << ',' << std::setprecision(7) << p.mjd << ',' << iso_from_mjd(p.mjd, 2)
              << ',' << std::setprecision(3) << p.x + 1.0 << ',' << p.y + 1.0 << ',' << std::setprecision(7)
              << p.ra_deg << ',' << p.dec_deg << ',' << std::setprecision(3) << c.mag_inst << '\n';
        }
    }
}

void write_truth_csv(const std::string& path, const SimResult& sim) {
    auto f = open_out(path);
    f << "# fainttrail simulation truth; FITS 1-based pixels at the reference epoch\n";
    f << "# t_ref_mjd=" << std::setprecision(8) << sim.t_ref_mjd << "\n";
    f << "id,x,y,vx_pxh,vy_pxh,mag,flux_e,expected_snr_frame,expected_snr_stack\n";
    for (const auto& t : sim.truth) {
        f << t.mover.id << ',' << std::setprecision(3) << t.mover.x_ref + 1.0 << ',' << t.mover.y_ref + 1.0 << ','
          << t.mover.vx_pxh << ',' << t.mover.vy_pxh << ',' << t.mover.mag << ',' << std::setprecision(1)
          << t.flux_e << ',' << std::setprecision(3) << t.expected_snr_frame << ',' << t.expected_snr_stack << '\n';
    }
}

std::vector<SimTruth> read_truth_csv(const std::string& path) {
    const CsvTable t = read_csv(path);
    std::vector<SimTruth> out;
    for (std::size_t r = 0; r < t.rows.size(); ++r) {
        SimTruth s;
        s.mover.id = t.str(r, "id");
        s.mover.x_ref = t.num(r, "x") - 1.0;
        s.mover.y_ref = t.num(r, "y") - 1.0;
        s.mover.vx_pxh = t.num(r, "vx_pxh");
        s.mover.vy_pxh = t.num(r, "vy_pxh");
        s.mover.mag = t.num(r, "mag");
        s.flux_e = t.num(r, "flux_e");
        s.expected_snr_frame = t.num(r, "expected_snr_frame");
        s.expected_snr_stack = t.num(r, "expected_snr_stack");
        out.push_back(s);
    }
    return out;
}

void write_completeness_csv(const std::string& path, const CompletenessReport& rep) {
    auto f = open_out(path);
    f << "# fainttrail injection-recovery; threshold_snr=" << std::setprecision(2) << rep.threshold_snr
      << " method=" << rep.threshold_method << " baseline_ok=" << rep.baseline_ok << " snr50=" << rep.snr50 << "\n";
    f << "target_snr,injected,recovered,detected,completeness,mean_snr,mean_dpos_px,mean_dv_pxh\n";
    for (const auto& r : rep.rows) {
        f << std::setprecision(2) << r.target_snr << ',' << r.injected << ',' << r.recovered << ',' << r.detected
          << ',' << std::setprecision(3) << (r.injected ? static_cast<double>(r.recovered) / r.injected : 0.0) << ','
          << std::setprecision(2) << r.mean_snr << ',' << std::setprecision(3) << r.mean_dpos << ',' << r.mean_dv
          << '\n';
    }
}

std::string summary_text(const SearchResult& res, bool ok_only) {
    std::ostringstream s;
    s.setf(std::ios::fixed);
    int nok = 0;
    for (const auto& c : res.candidates) nok += c.ok;
    s << "velocities searched : " << res.n_velocities << " (grid step " << std::setprecision(2) << res.dv_pxh
      << " px/h)\n";
    s << "PSF FWHM            : " << res.fwhm_px << " px\n";
    s << "slowest reliable    : " << res.slowest_reliable_pxh << " px/h (2 FWHM over the sequence)\n";
    s << "threshold           : S/N " << res.threshold_snr << " (" << res.threshold_method << ")\n";
    s << "candidates          : " << nok << " passed, " << res.candidates.size() - nok << " flagged\n";
    if (!res.backend.empty()) s << "backend             : " << res.backend_description << "\n";
    s << "elapsed             : " << std::setprecision(1) << res.elapsed_s << " s\n";
    char line[256];
    std::snprintf(line, sizeof line, "\n %3s %2s %8s %8s %8s %8s %6s %6s %6s %5s %8s  %s\n", "id", "ok", "x", "y",
                  "rate\"/m", "PA", "S/N", "half1", "half2", "frac", "FA/run", "flags");
    s << line;
    for (const auto& c : res.candidates) {
        if (ok_only && !c.ok) continue;
        std::snprintf(line, sizeof line, " %3d %2s %8.2f %8.2f %8.3f %8.1f %6.2f %6.2f %6.2f %5.2f %8.1e  %s\n", c.id,
                      c.ok ? "+" : "-", c.x + 1.0, c.y + 1.0, c.rate_arcsec_min, c.pa_deg, c.snr, c.snr_half1,
                      c.snr_half2, c.max_frame_fraction, c.false_alarms_per_run, c.flags.c_str());
        s << line;
    }
    return s.str();
}

std::string profile_text(const Prepared& prep, const SearchResult& res, double load_s) {
    struct Row {
        std::string name;
        double s = 0.0;
        int n = 0;
        bool sub = false;
        double rss = 0.0;  // highest resident memory seen at the end of this stage
    };
    std::vector<Row> rows;
    auto add = [&](const std::string& prefix, const StageTime& t) {
        const bool sub = t.name.rfind("  ", 0) == 0;
        const std::string name = sub ? "    " + t.name.substr(2) : prefix + t.name;
        for (auto& r : rows)
            if (r.name == name) {
                r.s += t.seconds;
                ++r.n;
                r.rss = std::max(r.rss, t.rss_mb);
                return;
            }
        rows.push_back({name, t.seconds, 1, sub, t.rss_mb});
    };
    if (load_s > 0.0) add("", {"load frames", load_s});
    for (const auto& t : prep.timings) add("prepare: ", t);
    for (const auto& t : res.timings) add("search:  ", t);
    double total = 0.0, stack = 0.0;
    int passes = 0;
    for (const auto& r : rows) {
        if (r.sub) continue;
        total += r.s;
        if (r.name.find("stack (") != std::string::npos) {
            stack += r.s;
            passes += r.n;
        }
    }
    std::ostringstream o;
    o << std::fixed;
    o << "profile (wall clock; null-run stages summed over " << std::max(0, passes - 1)
      << " runs; MB = resident memory at the end of the stage)\n";
    for (const auto& r : rows) {
        o << "  " << std::left << std::setw(34) << r.name << std::right << std::setw(9) << std::setprecision(2) << r.s
          << " s";
        if (!r.sub) o << std::setw(7) << std::setprecision(1) << (total > 0 ? 100.0 * r.s / total : 0.0) << " %";
        if (!r.sub && r.rss > 0.0) o << std::setw(8) << std::setprecision(0) << r.rss << " MB";
        if (r.n > 1) o << "   (" << r.n << "x)";
        o << "\n";
    }
    o << "  " << std::left << std::setw(34) << "total" << std::right << std::setw(9) << std::setprecision(2) << total
      << " s\n";
    const double peak = detail::peak_rss_mb();
    if (peak > 0.0) o << "peak memory of the process so far: " << std::setprecision(0) << peak << " MB\n";
    if (stack > 0.0 && passes > 0) {
        const double adds = static_cast<double>(res.n_velocities) * prep.nframes * prep.width * prep.height * passes;
        o << "stacking: " << passes << " passes x " << res.n_velocities << " motions x " << prep.nframes << " frames x "
          << prep.width << "x" << prep.height << " px = " << std::setprecision(3) << adds / 1e12
          << " T pixel-adds, " << std::setprecision(2) << adds / stack / 1e9 << " G/s\n";
    }
    return o.str();
}

}  // namespace fainttrail
