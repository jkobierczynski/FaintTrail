// SPDX-License-Identifier: GPL-3.0-or-later
// Copyright 2026 Jurgen Kobierczynski
//
// fainttrail command-line tool.
#include <algorithm>
#include <cctype>
#include <memory>
#include <chrono>
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
  backends                  list the stacking backends in this build
  crop FRAMES... -o DIR --x0 X --y0 Y --w W --h H [--bitpix 16] [--wcs-from F]
                            cut the same rectangle out of every frame
  regress CASES...          run real-data regression cases (.case files or dirs)
          [--backend NAME] [--strict] [--write-baseline] [--report FILE.csv] [--profile]
          [--tile N]        run the cases as tiled searches (must give the same results)
  check-backend [--backend NAME]
                            verify a backend against the stacking contract
  version                   print the version

motion range (search, completeness)
  --vmin V / --vmax V       pixels per hour (defaults 4 / 40)
  --rate-min R / --rate-max R
                            arcsec per minute; needs a WCS or --pixscale
  --pixscale S              arcsec per pixel when the frames carry no WCS

search options
  -o FILE                   candidates CSV (default candidates.csv)
  --tracklets FILE          per-frame positions of the passing candidates
  --profile                 time per pipeline stage and stacking throughput
  --null-model M            tail model of the threshold calibration: gauss
                            (default) or linear (FaintTrail <= 0.7)
  --tile N                  search in tiles of N x N px (same result, less memory)
  --max-memory MB           memory budget; tiles are used automatically when the
                            sequence would not fit (default: 60 % of the RAM)
  --no-tile                 never tile
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
  --backend NAME            stacking backend: auto (default) or one of
                            'fainttrail backends' 
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

const std::vector<std::string> kFlags = {"--all", "--quiet", "--help", "-h", "--strict", "--write-baseline", "--profile", "--no-tile"};

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

double frame_pixscale(const std::vector<Wcs>& wcs, const Args& a) {
    if (a.flag("--pixscale")) return a.num("--pixscale", 0.0);
    for (const auto& w : wcs)
        if (w.valid) return w.pixel_scale_arcsec();
    return 0.0;
}

std::vector<Wcs> wcs_of(const std::vector<Frame>& frames) {
    std::vector<Wcs> v;
    for (const auto& f : frames) v.push_back(f.wcs);
    return v;
}

SearchParams search_params(const Args& a, const std::vector<Wcs>& wcs) {
    SearchParams p;
    const bool rates = a.flag("--rate-min") || a.flag("--rate-max");
    if (rates) {
        const double scale = frame_pixscale(wcs, a);
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
    p.backend = a.str("--backend", p.backend);
    p.null_model = a.str("--null-model", p.null_model);
    if (p.null_model != "gauss" && p.null_model != "linear") throw std::runtime_error("--null-model must be gauss or linear");
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
    SearchParams p = search_params(a, wcs_of(sim.frames));
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

void report_transfer(const TransferReport& rep, double scale) {
    std::fprintf(stderr, "  %d stars matched, rms %.2f px, scale %.4f, rotation %.2f deg%s\n", rep.matched, rep.rms_px,
                 rep.scale, rep.rotation_deg, rep.mirrored ? ", mirrored" : "");
    std::fprintf(stderr, "  image scale %.3f\"/px\n", scale);
    if (rep.source_has_distortion)
        std::fprintf(stderr, "  note: distortion terms (SIP/PV) in the solution are ignored; positions near the\n"
                             "        field edges can be off by up to the size of that distortion\n");
}

int cmd_search(const Args& a) {
    const auto t_load = std::chrono::steady_clock::now();
    const auto files = expand_inputs(a.positional);
    if (files.empty()) throw std::runtime_error("no input FITS files given");
    // headers only; the pixels are read later, whole or tile by tile
    auto src = fits_frame_source(files);
    if (a.flag("--wcs-from")) {
        const std::string solved = a.str("--wcs-from");
        std::fprintf(stderr, "transferring plate solution from %s...\n", solved.c_str());
        const auto rep = apply_wcs_from_file(solved, *src, a.num("--wcs-scale", 0.0));
        report_transfer(rep, src->info(0).wcs.pixel_scale_arcsec());
    }
    std::vector<Wcs> wcs;
    for (int i = 0; i < src->count(); ++i) wcs.push_back(src->info(i).wcs);
    const SearchParams p = search_params(a, wcs);
    const ProgressFn prog = make_progress(a.flag("--quiet"));

    // tile when asked to, or when the whole sequence would not fit the memory budget
    const double physical = physical_memory_mb();
    const double budget = a.flag("--max-memory") ? a.num("--max-memory", 0.0) : (physical > 0 ? 0.6 * physical : 4096.0);
    const double untiled = estimate_untiled_memory_mb(src->width(), src->height(), src->count());
    const bool tiled = !a.flag("--no-tile") && (a.flag("--tile") || untiled > budget);

    Prepared meta;
    SearchResult res;
    double load_s = 0.0;
    if (tiled) {
        TileOptions to;
        to.tile = a.integer("--tile", 0);
        to.max_memory_mb = budget;
        const TilePlan plan = plan_tiles(*src, p, to);
        std::fprintf(stderr,
                     "tiled search: %d x %d tiles of %d px (halo %d + %d px), estimated memory %.0f MB "
                     "(untiled about %.0f MB, budget %.0f MB)\n",
                     plan.tiles_x, plan.tiles_y, plan.tile, plan.stack_halo, plan.filter_halo, plan.memory_mb, untiled,
                     budget);
        load_s = std::chrono::duration<double>(std::chrono::steady_clock::now() - t_load).count();
        res = search_tiled(*src, p, to, &meta, prog);
    } else {
        std::fprintf(stderr, "loading %zu frames...\n", files.size());
        auto frames = load_frames(files);
        for (int i = 0; i < src->count(); ++i) frames[static_cast<std::size_t>(i)].wcs = src->info(i).wcs;
        load_s = std::chrono::duration<double>(std::chrono::steady_clock::now() - t_load).count();
        ProgressFn prep_prog, search_prog;
        if (prog) {
            prep_prog = [&](double f, const std::string& s) { return prog(0.1 * f, s); };
            search_prog = [&](double f, const std::string& s) { return prog(0.1 + 0.9 * f, s); };
        }
        meta = prepare(std::move(frames), p, prep_prog);  // hands over the pixel buffers
        res = search(meta, p, search_prog);
    }
    const std::string out = a.str("-o", "candidates.csv");
    write_candidates_csv(out, res);
    if (a.flag("--tracklets")) write_tracklets_csv(a.str("--tracklets"), meta, res, true);
    if (a.flag("--snr-map")) {
        FitsHeader h;
        h.set("BUNIT", std::string("sigma"), "max S/N over all trial motions");
        Frame f;
        f.image = res.max_snr;
        f.mjd_mid = res.t_ref_mjd;
        f.wcs = meta.wcs;
        write_frame(a.str("--snr-map"), f, h);
    }
    std::printf("%s", summary_text(res, !a.flag("--all")).c_str());
    if (a.flag("--profile")) std::printf("\n%s", profile_text(meta, res, load_s).c_str());
    std::printf("\nwrote %s\n", out.c_str());
    return 0;
}

int cmd_stack(const Args& a) {
    auto frames = load_inputs(a);
    SearchParams p = search_params(a, wcs_of(frames));
    p.keep_residuals = true;  // stack_full works on the residuals
    const Prepared prep = prepare(std::move(frames), p);
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
    const SearchParams p = search_params(a, wcs_of(frames));
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

// Conformance check of a stacking backend against the contract in backend.hpp:
// first against the naive implementation on a small pass, then against the CPU
// reference on a realistic simulated sequence.
int cmd_check_backend(const Args& a) {
    const std::string name = a.str("--backend", "auto");
    auto cand = make_backend(name);
    std::printf("checking backend '%s' (%s)\n", cand->name().c_str(), cand->description().c_str());
    bool all_ok = true;
    auto report = [&](const char* what, const BackendCheck& c, double secs_ref, double secs_cand) {
        std::printf("  %-34s %s  max |dS/N| %.2e, best-v differences %zu (%zu unexplained)  [%.2f s vs %.2f s]\n", what,
                    c.ok ? "PASS" : "FAIL", c.max_abs_diff, c.best_v_mismatch, c.best_v_unexplained, secs_ref, secs_cand);
        if (!c.ok) {
            std::printf("    %s\n", c.message.c_str());
            all_ok = false;
        }
    };
    auto timed = [](auto&& fn) {
        const auto t0 = std::chrono::steady_clock::now();
        auto r = fn();
        return std::make_pair(std::move(r), std::chrono::duration<double>(std::chrono::steady_clock::now() - t0).count());
    };
    for (int round = 0; round < 2; ++round) {
        SimParams sp;
        sp.width = round == 0 ? 72 : 384;
        sp.height = round == 0 ? 52 : 320;
        sp.nframes = round == 0 ? 12 : 40;
        sp.nstars = round == 0 ? 20 : 400;
        sp.seed = 77 + static_cast<std::uint64_t>(round);
        const SimResult sim = simulate(sp);
        SearchParams p;
        p.vmin_pxh = 2;
        p.vmax_pxh = round == 0 ? 20 : 45;
        const Prepared prep = prepare(sim.frames, p);
        const auto grid = velocity_grid(prep, p);
        StackPass pass;
        pass.frames = &prep.filtered;
        pass.t_h = &prep.t_h;
        pass.velocities = &grid;
        pass.min_count = std::max(3, static_cast<int>(std::ceil(p.min_coverage * prep.nframes)));
        auto [got, tc] = timed([&] { return cand->run(pass); });
        if (round == 0) {
            auto [ref, tr] = timed([&] { return run_pass_naive(pass); });
            report("vs naive contract (72x52, 12 fr)", check_backend(pass, ref, got), tr, tc);
        } else {
            auto [ref, tr] = timed([&] { return run_pass_simple(pass); });
            report("vs simple reference (384x320, 40 fr)", check_backend(pass, ref, got), tr, tc);
        }
    }
    std::printf("%s\n", all_ok ? "backend conforms to the stacking contract" : "backend does NOT conform");
    return all_ok ? 0 : 3;
}

// Cut the same rectangle out of every frame of a registered sequence, keeping
// times, WCS and registration offsets consistent. Used to make test fixtures.
int cmd_crop(const Args& a) {
    const std::string dir = a.str("-o", a.str("--out"));
    if (dir.empty()) throw std::runtime_error("crop needs -o DIR");
    if (!a.flag("--x0") || !a.flag("--y0") || !a.flag("--w") || !a.flag("--h"))
        throw std::runtime_error("crop needs --x0 --y0 --w --h (0-based pixels)");
    const int bitpix = a.integer("--bitpix", -32);
    const auto frames = load_inputs(a);
    fs::create_directories(dir);
    const auto files = expand_inputs(a.positional);
    for (std::size_t i = 0; i < frames.size(); ++i) {
        const Frame f = crop_frame(frames[i], a.integer("--x0", 0), a.integer("--y0", 0), a.integer("--w", 0),
                                   a.integer("--h", 0));
        FitsHeader extra;
        extra.add_history("fainttrail crop of " + fs::path(files[i]).filename().string() + " at x0=" +
                          a.str("--x0") + " y0=" + a.str("--y0"));
        write_frame((fs::path(dir) / fs::path(files[i]).filename()).string(), f, extra, bitpix);
    }
    std::printf("wrote %zu cropped frames to %s\n", frames.size(), dir.c_str());
    return 0;
}

// Run real-data regression cases (see regression.hpp and tests/regression/).
int cmd_regress(const Args& a) {
    std::vector<std::string> cases;
    for (const auto& p : a.positional) {
        if (fs::is_directory(p)) {
            std::vector<std::string> v;
            for (const auto& e : fs::directory_iterator(p))
                if (e.path().extension() == ".case") v.push_back(e.path().string());
            std::sort(v.begin(), v.end());
            cases.insert(cases.end(), v.begin(), v.end());
        } else {
            cases.push_back(p);
        }
    }
    if (cases.empty()) throw std::runtime_error("regress needs one or more .case files or directories");
    RegressionOptions opt;
    opt.backend = a.str("--backend");
    opt.threads = a.integer("--threads", 0);
    opt.strict = a.flag("--strict");
    opt.tile = a.integer("--tile", 0);
    const bool write_base = a.flag("--write-baseline");
    std::FILE* report = nullptr;
    if (a.flag("--report")) {
        report = std::fopen(a.str("--report").c_str(), "w");
        if (!report) throw std::runtime_error("cannot write " + a.str("--report"));
        std::fprintf(report, "case,status,backend,frames,width,height,velocities,threshold,passed,candidates,load_s,prepare_s,search_s,memory_mb\n");
    }
    int n_fail = 0, n_changed = 0, n_skip = 0, n_pass = 0;
    for (const auto& path : cases) {
        const RegressionCase c = load_regression_case(path);
        std::printf("== %s%s%s\n", c.name.c_str(), c.description.empty() ? "" : ": ", c.description.c_str());
        const RegressionOutcome o = run_regression_case(c, opt, make_progress(a.flag("--quiet")));
        if (o.status == RegressionStatus::Skipped) {
            std::printf("   SKIP  %s\n\n", o.skip_reason.c_str());
            ++n_skip;
            continue;
        }
        int passed = 0;
        for (const auto& cd : o.result.candidates) passed += cd.ok;
        std::printf("   %zu frames %dx%d%s, backend %s, %zu motions, threshold %.2f, %d passed / %zu candidates\n",
                    o.nframes, o.width, o.height, o.tiles > 0 ? (", tiles of " + std::to_string(o.tiles) + " px").c_str() : "",
                    o.result.backend.c_str(), o.result.n_velocities,
                    o.result.threshold_snr, passed, o.result.candidates.size());
        std::printf("   time: load %.1f s, prepare %.1f s, search %.1f s; memory %.0f MB\n", o.load_s, o.prepare_s,
                    o.search_s, o.memory_mb);
        if (a.flag("--profile")) std::printf("%s", o.profile.c_str());
        for (const auto& ch : o.checks)
            std::printf("   %-7s %s\n", ch.ok ? "ok" : (ch.hard || opt.strict ? "FAIL" : "changed"), ch.text.c_str());
        if (write_base && !c.baseline_path.empty()) {
            write_regression_baseline(c.baseline_path, c, o);
            std::printf("   wrote baseline %s\n", c.baseline_path.c_str());
        }
        std::printf("   => %s\n\n", status_name(o.status));
        if (report)
            std::fprintf(report, "%s,%s,%s,%zu,%d,%d,%zu,%.3f,%d,%zu,%.2f,%.2f,%.2f,%.0f\n", c.name.c_str(),
                         status_name(o.status), o.result.backend.c_str(), o.nframes, o.width, o.height,
                         o.result.n_velocities, o.result.threshold_snr, passed, o.result.candidates.size(), o.load_s,
                         o.prepare_s, o.search_s, o.memory_mb);
        switch (o.status) {
            case RegressionStatus::Fail: ++n_fail; break;
            case RegressionStatus::Changed: ++n_changed; break;
            default: ++n_pass; break;
        }
    }
    if (report) std::fclose(report);
    std::printf("regression: %d passed, %d changed, %d failed, %d skipped\n", n_pass, n_changed, n_fail, n_skip);
    if (n_fail > 0) return 3;
    if (n_pass + n_changed == 0 && n_skip > 0) return 77;  // everything skipped (ctest SKIP_RETURN_CODE)
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
    fainttrail::tune_allocator_for_images();
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
        if (cmd == "backends") {
            for (const auto& n : backend_names()) {
                const auto b = make_backend(n);
                std::printf("%-8s %s\n", n.c_str(), b->description().c_str());
            }
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
        if (cmd == "check-backend") return cmd_check_backend(a);
        if (cmd == "crop") return cmd_crop(a);
        if (cmd == "regress") return cmd_regress(a);
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
