// SPDX-License-Identifier: GPL-3.0-or-later
// Copyright 2026 Jurgen Kobierczynski
//
// fainttrail command-line tool.
#include <algorithm>
#include <cctype>
#include <memory>
#include <cmath>
#include <cstdio>
#include <cstdlib>
#include <filesystem>
#include <iostream>
#include <map>
#include <sstream>
#include <string>
#include <vector>

#include "fainttrail/fainttrail.hpp"

namespace fs = std::filesystem;
using namespace fainttrail;

namespace {

const char* kUsage = R"(FaintTrail — synthetic tracking for faint asteroids and comets

usage: fainttrail <command> [options]

commands
  demo                      simulate a sequence, search it and score the result
  simulate -o DIR           write a synthetic FITS sequence plus truth.csv
  search FILES|DIR          blind search for moving objects
  stack FILES|DIR           shift-and-stack along a known motion
  completeness FILES|DIR    injection-recovery test on any data set
  validate                  compare a candidates file with a truth file
  version                   print the version

motion range (search, completeness)
  --vmin V / --vmax V       pixels per hour (defaults 4 / 40)
  --rate-min R / --rate-max R
                            arcsec per minute; needs a WCS or --pixscale
  --pixscale S              arcsec per pixel when the frames carry no WCS

search options
  -o FILE                   candidates CSV (default candidates.csv)
  --tracklets FILE          per-frame positions of the passing candidates
  --threshold S             fixed S/N threshold (default: calibrate)
  --null-runs N             time-shuffled null searches for calibration (3)
  --false-per-run F         accepted noise detections per search (0.02)
  --fwhm F                  PSF FWHM in pixels (default: measure)
  --grid-tol P              grid drift tolerance in pixels (0.5)
  --static-exclude N        neighbouring frames left out of each frame's
                            static-sky median (default: automatic)
  --snr-map FILE            write the max-S/N map as FITS
  --wcs-from FILE           borrow the plate solution of FILE (e.g. an ASTAP-solved
                            stack or frame of the same field); stars are matched,
                            so binning, cropping, 90-degree rotations and mirroring
                            are handled. Gives RA/Dec, rates and position angles.
  --wcs-scale S             frame pixels per pixel of FILE (default: from sizes)
  --threads N               worker threads (default: all cores)
  --all                     also print flagged candidates
  --quiet                   no progress output

stack options
  --vx V --vy V             motion in px/h, or --rate R --pa A with a WCS
  -o FILE                   output FITS (default stack.fits)
  --wcs-from FILE           as for search

simulate options
  --width W --height H --frames N --exptime S --cadence S --fwhm F
  --pixscale S --stars N --movers N --mag-bright M --mag-faint M
  --rate-min R --rate-max R (arcsec/min) --sky E (e-/px/s) --seed N

completeness options
  --levels a,b,c            target stack S/N levels (4,5,6,7,8,10,12,15)
  --per-level N             objects injected per level (10)
  -o FILE                   CSV report (default completeness.csv)

validate options
  --truth FILE --candidates FILE [--max-dpos P] [--max-dv V]
)";

struct Args {
    std::vector<std::string> positional;
    std::map<std::string, std::string> opts;
    bool flag(const std::string& k) const { return opts.count(k) > 0; }
    std::string str(const std::string& k, const std::string& d = {}) const {
        auto it = opts.find(k);
        return it == opts.end() ? d : it->second;
    }
    double num(const std::string& k, double d) const {
        auto it = opts.find(k);
        if (it == opts.end()) return d;
        try {
            return std::stod(it->second);
        } catch (...) {
            throw std::runtime_error("option " + k + " expects a number, got '" + it->second + "'");
        }
    }
    int integer(const std::string& k, int d) const { return static_cast<int>(std::lround(num(k, d))); }
};

const std::vector<std::string> kFlags = {"--all", "--quiet", "--help", "-h"};

Args parse(int argc, char** argv, int start) {
    Args a;
    for (int i = start; i < argc; ++i) {
        std::string s = argv[i];
        if (s.size() > 1 && s[0] == '-' && !(s.size() > 1 && (std::isdigit(static_cast<unsigned char>(s[1])) || s[1] == '.'))) {
            if (std::find(kFlags.begin(), kFlags.end(), s) != kFlags.end()) {
                a.opts[s] = "1";
            } else {
                if (i + 1 >= argc) throw std::runtime_error("option " + s + " needs a value");
                a.opts[s] = argv[++i];
            }
        } else {
            a.positional.push_back(s);
        }
    }
    return a;
}

std::vector<std::string> expand_inputs(const std::vector<std::string>& in) {
    std::vector<std::string> out;
    for (const auto& p : in) {
        if (fs::is_directory(p)) {
            std::vector<std::string> files;
            for (const auto& e : fs::directory_iterator(p)) {
                if (!e.is_regular_file()) continue;
                std::string ext = e.path().extension().string();
                std::transform(ext.begin(), ext.end(), ext.begin(), [](unsigned char c) { return std::tolower(c); });
                if (ext == ".fits" || ext == ".fit" || ext == ".fts") files.push_back(e.path().string());
            }
            std::sort(files.begin(), files.end());
            out.insert(out.end(), files.begin(), files.end());
        } else {
            out.push_back(p);
        }
    }
    return out;
}

ProgressFn make_progress(bool quiet) {
    if (quiet) return {};
    auto last = std::make_shared<int>(-1);
    auto last_stage = std::make_shared<std::string>();
    return [last, last_stage](double f, const std::string& stage) {
        const int pct = static_cast<int>(f * 100.0);
        if (pct != *last || stage != *last_stage) {
            *last = pct;
            *last_stage = stage;
            const int bar = pct / 4;
            std::fprintf(stderr, "\r  [%-25s] %3d%%  %-22s", std::string(static_cast<std::size_t>(bar), '#').c_str(),
                         pct, stage.c_str());
            if (f >= 1.0) std::fprintf(stderr, "\n");
            std::fflush(stderr);
        }
        return true;
    };
}

double frame_pixscale(const std::vector<Frame>& frames, const Args& a) {
    if (a.flag("--pixscale")) return a.num("--pixscale", 0.0);
    for (const auto& f : frames)
        if (f.wcs.valid) return f.wcs.pixel_scale_arcsec();
    return 0.0;
}

SearchParams search_params(const Args& a, const std::vector<Frame>& frames) {
    SearchParams p;
    const bool rates = a.flag("--rate-min") || a.flag("--rate-max");
    if (rates) {
        const double scale = frame_pixscale(frames, a);
        if (scale <= 0.0) throw std::runtime_error("--rate-min/--rate-max need a WCS in the frames or --pixscale");
        p.vmin_pxh = arcsec_per_min_to_px_per_h(a.num("--rate-min", px_per_h_to_arcsec_per_min(p.vmin_pxh, scale)), scale);
        p.vmax_pxh = arcsec_per_min_to_px_per_h(a.num("--rate-max", px_per_h_to_arcsec_per_min(p.vmax_pxh, scale)), scale);
    }
    p.vmin_pxh = a.num("--vmin", p.vmin_pxh);
    p.vmax_pxh = a.num("--vmax", p.vmax_pxh);
    p.threshold_snr = a.num("--threshold", p.threshold_snr);
    p.null_runs = a.integer("--null-runs", p.null_runs);
    p.target_false_per_run = a.num("--false-per-run", p.target_false_per_run);
    p.fwhm_px = a.num("--fwhm", p.fwhm_px);
    p.grid_tol_px = a.num("--grid-tol", p.grid_tol_px);
    p.static_exclude = a.integer("--static-exclude", p.static_exclude);
    p.mask_scatter = a.num("--mask-scatter", p.mask_scatter);
    p.weight_scatter = a.num("--weight-scatter", p.weight_scatter);
    p.threads = a.integer("--threads", p.threads);
    p.seed = static_cast<std::uint64_t>(a.num("--seed", static_cast<double>(p.seed)));
    return p;
}

std::vector<Frame> load_inputs(const Args& a) {
    const auto files = expand_inputs(a.positional);
    if (files.empty()) throw std::runtime_error("no input FITS files given");
    std::fprintf(stderr, "loading %zu frames...\n", files.size());
    auto frames = load_frames(files);
    if (a.flag("--wcs-from")) {
        const std::string solved = a.str("--wcs-from");
        std::fprintf(stderr, "transferring plate solution from %s...\n", solved.c_str());
        const auto rep = apply_wcs_from_file(solved, frames, a.num("--wcs-scale", 0.0));
        std::fprintf(stderr, "  %d stars matched, rms %.2f px, scale %.4f, rotation %.2f deg%s\n", rep.matched,
                     rep.rms_px, rep.scale, rep.rotation_deg, rep.mirrored ? ", mirrored" : "");
        std::fprintf(stderr, "  image scale %.3f\"/px\n", frames.front().wcs.pixel_scale_arcsec());
        if (rep.source_has_distortion)
            std::fprintf(stderr, "  note: distortion terms (SIP/PV) in the solution are ignored; positions near the\n"
                                 "        field edges can be off by up to the size of that distortion\n");
    }
    return frames;
}

SimParams sim_params(const Args& a) {
    SimParams s;
    s.width = a.integer("--width", s.width);
    s.height = a.integer("--height", s.height);
    s.nframes = a.integer("--frames", s.nframes);
    s.exptime_s = a.num("--exptime", s.exptime_s);
    s.cadence_s = a.num("--cadence", std::max(s.cadence_s, s.exptime_s));
    s.fwhm_px = a.num("--fwhm", s.fwhm_px);
    s.pixscale_arcsec = a.num("--pixscale", s.pixscale_arcsec);
    s.nstars = a.integer("--stars", s.nstars);
    s.nmovers = a.integer("--movers", s.nmovers);
    s.mover_mag_bright = a.num("--mag-bright", s.mover_mag_bright);
    s.mover_mag_faint = a.num("--mag-faint", s.mover_mag_faint);
    s.rate_min_arcsec_min = a.num("--rate-min", s.rate_min_arcsec_min);
    s.rate_max_arcsec_min = a.num("--rate-max", s.rate_max_arcsec_min);
    s.sky_e_per_s = a.num("--sky", s.sky_e_per_s);
    s.seed = static_cast<std::uint64_t>(a.num("--seed", 1));
    if (s.cadence_s < s.exptime_s) throw std::runtime_error("--cadence must be at least --exptime");
    return s;
}

void print_truth_score(const SearchResult& res, const SimResult& sim) {
    std::vector<SimMover> movers;
    for (const auto& t : sim.truth) movers.push_back(t.mover);
    const double max_dv = std::max(3.0, 2.0 * res.dv_pxh);
    const auto m = match(res.candidates, movers, 2.0, max_dv, true);
    int nok = 0;
    for (const auto& c : res.candidates) nok += c.ok;
    std::printf("\ntruth comparison (match: < 2 px at t_ref, < %.1f px/h)\n", max_dv);
    std::printf(" %-4s %6s %8s %8s %8s  %s\n", "id", "mag", "exp S/N", "found", "S/N", "dpos/dv");
    for (std::size_t t = 0; t < sim.truth.size(); ++t) {
        const auto& tr = sim.truth[t];
        auto it = std::find_if(m.begin(), m.end(), [&](const MatchResult& r) { return r.truth_index == static_cast<int>(t); });
        if (it != m.end())
            std::printf(" %-4s %6.2f %8.1f %8s %8.1f  %.2f px / %.2f px/h\n", tr.mover.id.c_str(), tr.mover.mag,
                        tr.expected_snr_stack, "yes", res.candidates[static_cast<std::size_t>(it->candidate_index)].snr,
                        it->dpos_px, it->dv_pxh);
        else
            std::printf(" %-4s %6.2f %8.1f %8s\n", tr.mover.id.c_str(), tr.mover.mag, tr.expected_snr_stack, "no");
    }
    std::printf("recovered %zu of %zu movers; %d passing candidates are not in the truth list\n", m.size(),
                sim.truth.size(), nok - static_cast<int>(m.size()));
}

int cmd_demo(const Args& a) {
    SimParams sp = sim_params(a);
    if (!a.flag("--width")) sp.width = 384;
    if (!a.flag("--height")) sp.height = 384;
    std::fprintf(stderr, "simulating %d frames of %dx%d px with %d movers (seed %llu)...\n", sp.nframes, sp.width,
                 sp.height, sp.nmovers, static_cast<unsigned long long>(sp.seed));
    const SimResult sim = simulate(sp);
    SearchParams p = search_params(a, sim.frames);
    if (!a.flag("--vmax") && !a.flag("--rate-max"))
        p.vmax_pxh = arcsec_per_min_to_px_per_h(sp.rate_max_arcsec_min, sp.pixscale_arcsec) * 1.1;
    const SearchResult res = search(sim.frames, p, make_progress(a.flag("--quiet")));
    std::printf("%s", summary_text(res, !a.flag("--all")).c_str());
    print_truth_score(res, sim);
    return 0;
}

int cmd_simulate(const Args& a) {
    const std::string dir = a.str("-o", a.str("--out"));
    if (dir.empty()) throw std::runtime_error("simulate needs -o DIR");
    fs::create_directories(dir);
    const SimParams sp = sim_params(a);
    const SimResult sim = simulate(sp);
    for (std::size_t i = 0; i < sim.frames.size(); ++i) {
        char name[64];
        std::snprintf(name, sizeof name, "frame_%03zu.fits", i + 1);
        FitsHeader extra;
        extra.set("OBJECT", std::string("fainttrail simulation"));
        extra.set_int("SIMSEED", static_cast<long long>(sp.seed), "fainttrail simulation seed");
        extra.set("GAIN", sp.gain_e_per_adu, "e-/ADU");
        write_frame((fs::path(dir) / name).string(), sim.frames[i], extra);
    }
    write_truth_csv((fs::path(dir) / "truth.csv").string(), sim);
    std::printf("wrote %zu frames and truth.csv to %s\n", sim.frames.size(), dir.c_str());
    for (const auto& t : sim.truth)
        std::printf("  %-4s mag %.2f  rate %.3f\"/min  expected stack S/N %.1f\n", t.mover.id.c_str(), t.mover.mag,
                    px_per_h_to_arcsec_per_min(std::hypot(t.mover.vx_pxh, t.mover.vy_pxh), sp.pixscale_arcsec),
                    t.expected_snr_stack);
    return 0;
}

int cmd_search(const Args& a) {
    const auto frames = load_inputs(a);
    const SearchParams p = search_params(a, frames);
    const ProgressFn prog = make_progress(a.flag("--quiet"));
    ProgressFn prep_prog, search_prog;
    if (prog) {
        prep_prog = [&](double f, const std::string& s) { return prog(0.1 * f, s); };
        search_prog = [&](double f, const std::string& s) { return prog(0.1 + 0.9 * f, s); };
    }
    const Prepared prep = prepare(frames, p, prep_prog);
    const SearchResult res = search(prep, p, search_prog);
    const std::string out = a.str("-o", "candidates.csv");
    write_candidates_csv(out, res);
    if (a.flag("--tracklets")) write_tracklets_csv(a.str("--tracklets"), prep, res, true);
    if (a.flag("--snr-map")) {
        FitsHeader h;
        h.set("BUNIT", std::string("sigma"), "max S/N over all trial motions");
        Frame f;
        f.image = res.max_snr;
        f.mjd_mid = res.t_ref_mjd;
        f.wcs = prep.wcs;
        write_frame(a.str("--snr-map"), f, h);
    }
    std::printf("%s", summary_text(res, !a.flag("--all")).c_str());
    std::printf("\nwrote %s\n", out.c_str());
    return 0;
}

int cmd_stack(const Args& a) {
    const auto frames = load_inputs(a);
    SearchParams p = search_params(a, frames);
    const Prepared prep = prepare(frames, p);
    double vx = a.num("--vx", NAN), vy = a.num("--vy", NAN);
    if (a.flag("--rate")) {
        if (!prep.wcs.valid) throw std::runtime_error("--rate/--pa need a WCS; use --vx/--vy");
        // convert rate + PA (east of north) into pixel motion via the WCS
        const double rate = a.num("--rate", 0.0), pa = a.num("--pa", 0.0) * 3.14159265358979 / 180.0;
        const double cx = 0.5 * prep.width, cy = 0.5 * prep.height;
        double ra, dec, x2, y2;
        prep.wcs.pix2sky(cx, cy, ra, dec);
        const double dt_min = 60.0;
        const double east = rate * dt_min * std::sin(pa) / 3600.0, north = rate * dt_min * std::cos(pa) / 3600.0;
        prep.wcs.sky2pix(ra + east / std::cos(dec * 3.14159265358979 / 180.0), dec + north, x2, y2);
        vx = x2 - cx;
        vy = y2 - cy;
    }
    if (!std::isfinite(vx) || !std::isfinite(vy)) throw std::runtime_error("stack needs --vx/--vy or --rate/--pa");
    Frame out;
    out.image = stack_full(prep, vx, vy);
    out.mjd_mid = prep.t_ref_mjd;
    out.wcs = prep.wcs;
    FitsHeader h;
    h.set("STVX", vx, "stack motion x [px/h]");
    h.set("STVY", vy, "stack motion y [px/h]");
    h.set_int("NCOMBINE", prep.nframes, "frames stacked");
    const std::string path = a.str("-o", "stack.fits");
    write_frame(path, out, h);
    std::printf("stacked %d frames at (%.3f, %.3f) px/h, reference epoch %s UTC -> %s\n", prep.nframes, vx, vy,
                iso_from_mjd(prep.t_ref_mjd, 2).c_str(), path.c_str());
    return 0;
}

int cmd_completeness(const Args& a) {
    const auto frames = load_inputs(a);
    const SearchParams p = search_params(a, frames);
    CompletenessParams cp;
    if (a.flag("--levels")) {
        cp.snr_levels.clear();
        std::stringstream ss(a.str("--levels"));
        std::string item;
        while (std::getline(ss, item, ',')) cp.snr_levels.push_back(std::stod(item));
    }
    cp.per_level = a.integer("--per-level", cp.per_level);
    cp.seed = static_cast<std::uint64_t>(a.num("--seed", static_cast<double>(cp.seed)));
    const auto rep = completeness(frames, p, cp, make_progress(a.flag("--quiet")));
    const std::string out = a.str("-o", "completeness.csv");
    write_completeness_csv(out, rep);
    std::printf("threshold S/N %.2f (%s); %d candidates passed in the data before injection\n\n", rep.threshold_snr,
                rep.threshold_method.c_str(), rep.baseline_ok);
    std::printf(" %8s %9s %10s %9s %9s %9s\n", "S/N", "injected", "recovered", "detected", "complete", "mean S/N");
    for (const auto& r : rep.rows)
        std::printf(" %8.1f %9d %10d %9d %8.0f%% %9.2f\n", r.target_snr, r.injected, r.recovered, r.detected,
                    r.injected ? 100.0 * r.recovered / r.injected : 0.0, r.mean_snr);
    if (rep.snr50 > 0) std::printf("\n50%% completeness at stack S/N %.2f\n", rep.snr50);
    std::printf("wrote %s\n", out.c_str());
    return 0;
}

int cmd_validate(const Args& a) {
    const std::string tpath = a.str("--truth"), cpath = a.str("--candidates", "candidates.csv");
    if (tpath.empty()) throw std::runtime_error("validate needs --truth FILE");
    const auto truth = read_truth_csv(tpath);
    const auto cands = read_candidates_csv(cpath);
    std::vector<SimMover> movers;
    for (const auto& t : truth) movers.push_back(t.mover);
    const double max_dpos = a.num("--max-dpos", 2.0), max_dv = a.num("--max-dv", 3.0);
    const auto m = match(cands, movers, max_dpos, max_dv, true);
    int nok = 0;
    for (const auto& c : cands) nok += c.ok;
    std::printf(" %-6s %8s %8s %10s %10s\n", "id", "exp S/N", "found", "dpos px", "dv px/h");
    for (std::size_t t = 0; t < truth.size(); ++t) {
        auto it = std::find_if(m.begin(), m.end(), [&](const MatchResult& r) { return r.truth_index == static_cast<int>(t); });
        if (it != m.end())
            std::printf(" %-6s %8.1f %8d %10.3f %10.3f\n", truth[t].mover.id.c_str(), truth[t].expected_snr_stack,
                        cands[static_cast<std::size_t>(it->candidate_index)].id, it->dpos_px, it->dv_pxh);
        else
            std::printf(" %-6s %8.1f %8s\n", truth[t].mover.id.c_str(), truth[t].expected_snr_stack, "-");
    }
    const int fp = nok - static_cast<int>(m.size());
    std::printf("\nrecall %zu/%zu, passing candidates without a match: %d\n", m.size(), truth.size(), fp);
    return 0;
}

}  // namespace

int main(int argc, char** argv) {
    if (argc < 2) {
        std::fputs(kUsage, stdout);
        return 1;
    }
    const std::string cmd = argv[1];
    try {
        if (cmd == "-h" || cmd == "--help" || cmd == "help") {
            std::fputs(kUsage, stdout);
            return 0;
        }
        if (cmd == "version" || cmd == "--version") {
            std::printf("FaintTrail %s\nCopyright (C) 2026 Jurgen Kobierczynski\nLicense GPLv3+: GNU GPL version 3 or later <https://gnu.org/licenses/gpl.html>.\nThis is free software: you are free to change and redistribute it.\nThere is NO WARRANTY, to the extent permitted by law.\n", version());
            return 0;
        }
        const Args a = parse(argc, argv, 2);
        if (a.flag("--help") || a.flag("-h")) {
            std::fputs(kUsage, stdout);
            return 0;
        }
        if (cmd == "demo") return cmd_demo(a);
        if (cmd == "simulate") return cmd_simulate(a);
        if (cmd == "search") return cmd_search(a);
        if (cmd == "stack") return cmd_stack(a);
        if (cmd == "completeness") return cmd_completeness(a);
        if (cmd == "validate") return cmd_validate(a);
        std::fprintf(stderr, "unknown command '%s'\n\n%s", cmd.c_str(), kUsage);
        return 1;
    } catch (const Cancelled&) {
        std::fprintf(stderr, "\ncancelled\n");
        return 2;
    } catch (const std::exception& e) {
        std::fprintf(stderr, "\nerror: %s\n", e.what());
        return 1;
    }
}
