// SPDX-License-Identifier: GPL-3.0-or-later
// Copyright 2026 Jurgen Kobierczynski
//
// How should the detection threshold be calibrated? (docs/CALIBRATION.md)
//
// The threshold comes from null runs: the time-dependent pipeline is rerun with
// the frame times permuted, the null detection maps are searched for peaks, and
// log N(>S) = a - b S is fitted to the pooled peaks. Three full null runs cost
// three times the real stacking. This tool measures how well cheaper estimators
// reproduce the threshold of a large reference pool of null runs:
//
//   full k       k full null runs (k = 1, 2, 3; 3 is the current default)
//   area k x f   k null runs whose peaks are counted only in a random fraction f
//                of the field (square calibration tiles); stacking costs k * f runs
//   vel k x f    k null runs on a random fraction f of the trial motions
//
//   fainttrail_calibration_study (--data DIR [--wcs-from F] | --sim) [--vmin V] [--vmax V]
//       [--perms N] [--resamples R] [--name NAME] [--csv FILE]
#include <algorithm>
#include <chrono>
#include <cmath>
#include <cstdio>
#include <cstdlib>
#include <filesystem>
#include <fstream>
#include <map>
#include <numeric>
#include <string>
#include <vector>

#include "../src/internal.hpp"
#include "../src/random.hpp"
#include "fainttrail/fainttrail.hpp"

using namespace fainttrail;
namespace fs = std::filesystem;
using clk = std::chrono::steady_clock;

namespace {

double secs(clk::time_point t) { return std::chrono::duration<double>(clk::now() - t).count(); }

struct NullRun {
    std::vector<detail::Peak> peaks;  // full grid, whole field
    std::vector<float> vel_peaks;     // a quarter of the motions, whole field
    std::vector<float> vel2_peaks;    // half of the motions
    double prep_s = 0, stack_s = 0, vel_stack_s = 0, vel2_stack_s = 0;
};

bool threshold_from(const std::vector<float>& pooled, double runs, double target, double floor, double& T, double& b) {
    double a = 0;
    if (!detail::fit_null_tail(pooled, runs, a, b)) return false;
    T = std::max((a - std::log(target)) / b, floor);
    return true;
}

double quantile(std::vector<double> v, double q) {
    if (v.empty()) return NAN;
    std::sort(v.begin(), v.end());
    const double pos = q * (v.size() - 1);
    const std::size_t i = static_cast<std::size_t>(pos);
    const double f = pos - i;
    return i + 1 < v.size() ? v[i] * (1 - f) + v[i + 1] * f : v[i];
}

}  // namespace

int main(int argc, char** argv) {
    std::string data, wcs_from, name = "data", csv, dump;
    bool sim = false;
    int perms = 24, resamples = 400;
    SearchParams params;
    for (int i = 1; i < argc; ++i) {
        const std::string k = argv[i];
        auto next = [&] { return std::string(i + 1 < argc ? argv[++i] : ""); };
        if (k == "--data") data = next();
        else if (k == "--wcs-from") wcs_from = next();
        else if (k == "--sim") sim = true;
        else if (k == "--vmin") params.vmin_pxh = std::atof(next().c_str());
        else if (k == "--vmax") params.vmax_pxh = std::atof(next().c_str());
        else if (k == "--perms") perms = std::atoi(next().c_str());
        else if (k == "--resamples") resamples = std::atoi(next().c_str());
        else if (k == "--name") name = next();
        else if (k == "--csv") csv = next();
        else if (k == "--dump") dump = next();
    }
    tune_allocator_for_images();

    // ---- data
    std::vector<Frame> frames;
    if (sim) {
        SimParams sp;
        sp.width = 512;
        sp.height = 512;
        sp.nframes = 60;
        sp.nstars = 800;
        sp.nmovers = 0;
        sp.seed = 4711;
        frames = simulate(sp).frames;
        if (params.vmax_pxh <= 0 || params.vmax_pxh == SearchParams{}.vmax_pxh) params.vmax_pxh = 30;
        params.vmin_pxh = std::max(params.vmin_pxh, 5.0);
    } else {
        std::vector<std::string> files;
        for (const auto& e : fs::directory_iterator(data))
            if (e.path().extension() == ".fits" || e.path().extension() == ".fit") files.push_back(e.path().string());
        std::sort(files.begin(), files.end());
        frames = load_frames(files);
        if (!wcs_from.empty()) apply_wcs_from_file(wcs_from, frames);
    }
    const Prepared p = prepare(frames, params);
    frames.clear();
    double dv = 0;
    const auto grid = velocity_grid(p, params, &dv);
    const auto backend = make_backend(params.backend);
    const int nms = std::max(2, static_cast<int>(std::lround(p.fwhm_px)));
    const int edge = std::max(params.edge_px, 1);
    std::printf("%s: %d frames %dx%d, %zu motions, FWHM %.2f px, %d null runs\n", name.c_str(), p.nframes, p.width,
                p.height, grid.size(), p.fwhm_px, perms);

    // a random quarter of the motions per run (fixed per run)
    std::vector<NullRun> runs(static_cast<std::size_t>(perms));
    for (int k = 0; k < perms; ++k) {
        NullRun& r = runs[static_cast<std::size_t>(k)];
        const auto t = detail::null_time_orders(p.t_h, 1, 1000 + static_cast<std::uint64_t>(k)).front();
        auto t0 = clk::now();
        const Prepared q = detail::with_times(p, t, params);
        r.prep_s = secs(t0);
        t0 = clk::now();
        const auto full = detail::run_pass(*backend, q, grid, t, params, {}, 0, 1, "null");
        r.stack_s = secs(t0);
        r.peaks = detail::find_peaks(full, 3.0f, nms, edge, 10000000);
        std::vector<std::pair<double, double>> sub;
        detail::Rng rng(5000 + static_cast<std::uint64_t>(k));
        for (const auto& v : grid)
            if (rng.uniform() < 0.25) sub.push_back(v);
        t0 = clk::now();
        const auto vel = detail::run_pass(*backend, q, sub, t, params, {}, 0, 1, "null");
        r.vel_stack_s = secs(t0);
        for (const auto& pk : detail::find_peaks(vel, 3.0f, nms, edge, 10000000)) r.vel_peaks.push_back(pk.snr);
        if (!dump.empty()) {
            std::vector<std::pair<double, double>> half;
            detail::Rng rng2(7000 + static_cast<std::uint64_t>(k));
            for (const auto& v : grid)
                if (rng2.uniform() < 0.5) half.push_back(v);
            t0 = clk::now();
            const auto vel2 = detail::run_pass(*backend, q, half, t, params, {}, 0, 1, "null");
            r.vel2_stack_s = secs(t0);
            for (const auto& pk : detail::find_peaks(vel2, 3.0f, nms, edge, 10000000)) r.vel2_peaks.push_back(pk.snr);
        }
        std::printf("  null run %2d: %zu peaks > 3 (prep %.1f s, stack %.1f s)\n", k + 1, r.peaks.size(), r.prep_s,
                    r.stack_s);
        std::fflush(stdout);
    }
    if (!dump.empty()) {
        // kind: f = full grid (x, y), q = quarter of the motions, h = half of the motions
        std::ofstream d(dump);
        d << "# dataset=" << name << " width=" << p.width << " height=" << p.height << " frames=" << p.nframes
          << " motions=" << grid.size() << " fwhm=" << p.fwhm_px << " target=" << params.target_false_per_run
          << " floor=" << params.min_snr_floor << "\n";
        d << "run,kind,x,y,snr\n";
        for (int k = 0; k < perms; ++k) {
            const NullRun& r = runs[static_cast<std::size_t>(k)];
            d << k << ",T,0,0," << r.prep_s << "\n" << k << ",S,0,0," << r.stack_s << "\n" << k << ",Q,0,0,"
              << r.vel_stack_s << "\n" << k << ",H,0,0," << r.vel2_stack_s << "\n";
            for (const auto& pk : r.peaks) d << k << ",f," << pk.x << ',' << pk.y << ',' << pk.snr << "\n";
            for (float s : r.vel_peaks) d << k << ",q,0,0," << s << "\n";
            for (float s : r.vel2_peaks) d << k << ",h,0,0," << s << "\n";
        }
        std::printf("wrote %s\n", dump.c_str());
    }
    double prep_mean = 0, stack_mean = 0;
    for (const auto& r : runs) {
        prep_mean += r.prep_s / perms;
        stack_mean += r.stack_s / perms;
    }

    // ---- reference: all runs pooled
    std::vector<float> all;
    for (const auto& r : runs)
        for (const auto& pk : r.peaks) all.push_back(pk.snr);
    double T_ref = 0, b_ref = 0;
    if (!threshold_from(all, perms, params.target_false_per_run, params.min_snr_floor, T_ref, b_ref)) {
        std::printf("reference fit failed\n");
        return 1;
    }
    std::printf("reference (%d runs): threshold %.3f, tail slope b = %.2f per sigma\n", perms, T_ref, b_ref);

    // calibration tiles: an 8 x 8 grid over the field
    const int gt = 8;
    auto tile_of = [&](const detail::Peak& pk) {
        const int tx = std::min(gt - 1, pk.x * gt / p.width), ty = std::min(gt - 1, pk.y * gt / p.height);
        return ty * gt + tx;
    };

    struct Method {
        std::string label;
        char kind;  // 'f' full, 'a' area, 'v' velocity
        int k;
        double f;
    };
    const std::vector<Method> methods = {
        {"full 1", 'f', 1, 1.0},        {"full 2", 'f', 2, 1.0},         {"full 3 (current)", 'f', 3, 1.0},
        {"area 3 x 1/2", 'a', 3, 0.5},  {"area 3 x 1/4", 'a', 3, 0.25},  {"area 6 x 1/4", 'a', 6, 0.25},
        {"area 4 x 1/8", 'a', 4, 0.125}, {"area 8 x 1/8", 'a', 8, 0.125}, {"vel 3 x 1/4", 'v', 3, 0.25},
        {"vel 6 x 1/4", 'v', 6, 0.25},
    };
    std::ofstream out;
    if (!csv.empty()) {
        const bool fresh = !fs::exists(csv);
        out.open(csv, std::ios::app);
        if (fresh)
            out << "dataset,method,stack_runs,prep_runs,cost_rel,dT_median,dT_p10,dT_p90,dT_sd,fa_ratio_median,"
                   "fa_ratio_p90,fails\n";
    }
    std::printf("\n%-18s %6s %6s %8s %9s %18s %8s %14s\n", "method", "stack", "prep", "cost", "dT med", "dT 10-90 %",
                "dT sd", "FA ratio med/90%");
    const double cost_full3 = 3 * (prep_mean + stack_mean);
    detail::Rng rng(99);
    for (const auto& m : methods) {
        if (m.k > perms) continue;
        std::vector<double> dts;
        int fails = 0;
        for (int rep = 0; rep < resamples; ++rep) {
            // k distinct runs
            std::vector<int> idx(static_cast<std::size_t>(perms));
            std::iota(idx.begin(), idx.end(), 0);
            for (int i = 0; i < m.k; ++i) {
                const int j = i + static_cast<int>(rng.uniform() * (perms - i));
                std::swap(idx[static_cast<std::size_t>(i)], idx[static_cast<std::size_t>(std::min(j, perms - 1))]);
            }
            std::vector<float> pooled;
            double runs_eff = 0;
            for (int i = 0; i < m.k; ++i) {
                const NullRun& r = runs[static_cast<std::size_t>(idx[static_cast<std::size_t>(i)])];
                if (m.kind == 'f') {
                    for (const auto& pk : r.peaks) pooled.push_back(pk.snr);
                    runs_eff += 1;
                } else if (m.kind == 'v') {
                    pooled.insert(pooled.end(), r.vel_peaks.begin(), r.vel_peaks.end());
                    runs_eff += m.f;  // naive scaling with the fraction of motions
                } else {
                    const int ntiles = gt * gt, nsel = std::max(1, static_cast<int>(std::lround(m.f * ntiles)));
                    std::vector<int> t(static_cast<std::size_t>(ntiles));
                    std::iota(t.begin(), t.end(), 0);
                    for (int a = 0; a < nsel; ++a) {
                        const int j = a + static_cast<int>(rng.uniform() * (ntiles - a));
                        std::swap(t[static_cast<std::size_t>(a)], t[static_cast<std::size_t>(std::min(j, ntiles - 1))]);
                    }
                    std::vector<char> chosen(static_cast<std::size_t>(ntiles), 0);
                    for (int a = 0; a < nsel; ++a) chosen[static_cast<std::size_t>(t[static_cast<std::size_t>(a)])] = 1;
                    for (const auto& pk : r.peaks)
                        if (chosen[static_cast<std::size_t>(tile_of(pk))]) pooled.push_back(pk.snr);
                    runs_eff += static_cast<double>(nsel) / ntiles;
                }
            }
            double T = 0, b = 0;
            if (!threshold_from(pooled, runs_eff, params.target_false_per_run, params.min_snr_floor, T, b)) {
                ++fails;
                continue;
            }
            dts.push_back(T - T_ref);
        }
        double stack_runs = m.kind == 'f' ? m.k : m.k * m.f;
        const double prep_runs = m.k;
        if (m.kind == 'v') stack_runs = m.k * (runs[0].vel_stack_s / std::max(1e-9, runs[0].stack_s));
        const double cost = (prep_runs * prep_mean + stack_runs * stack_mean) / cost_full3;
        const double med = quantile(dts, 0.5), p10 = quantile(dts, 0.1), p90 = quantile(dts, 0.9);
        double sd = 0;
        for (double d : dts) sd += (d - med) * (d - med);
        sd = dts.size() > 1 ? std::sqrt(sd / (dts.size() - 1)) : NAN;
        // true false-alarm rate at the estimated threshold, relative to the requested one
        const double fa_med = std::exp(-b_ref * med), fa_p90 = std::exp(-b_ref * p10);
        std::printf("%-18s %6.2f %6.0f %7.0f%% %+9.3f %+8.3f..%+7.3f %8.3f %7.2f / %5.2f%s\n", m.label.c_str(), stack_runs,
                    prep_runs, 100 * cost, med, p10, p90, sd, fa_med, fa_p90,
                    fails ? ("  (" + std::to_string(fails) + " fits failed)").c_str() : "");
        if (out)
            out << name << ',' << m.label << ',' << stack_runs << ',' << prep_runs << ',' << cost << ',' << med << ','
                << p10 << ',' << p90 << ',' << sd << ',' << fa_med << ',' << fa_p90 << ',' << fails << '\n';
    }
    std::printf("\ncost: relative to 3 full null runs (prep %.1f s + stack %.1f s per run)\n", prep_mean, stack_mean);
    std::printf("FA ratio: false-alarm rate actually obtained / requested, under the reference tail model\n");
    return 0;
}
