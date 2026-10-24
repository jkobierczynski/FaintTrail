// SPDX-License-Identifier: GPL-3.0-or-later
// Copyright 2026 Jurgen Kobierczynski
//
// fainttrail_benchmark: validation of the engine on synthetic data.
//
//   1. false alarms: many pure-noise sequences, passing candidates counted and
//      compared with the requested false-alarm rate;
//   2. completeness: injection-recovery as a function of stack S/N;
//   3. astrometry: position and motion errors against the reported uncertainties.
//
// usage: fainttrail_benchmark OUTDIR [--quick]
#include <cmath>
#include <cstdio>
#include <filesystem>
#include <fstream>
#include <string>
#include <vector>

#include "fainttrail/fainttrail.hpp"

using namespace fainttrail;
namespace fs = std::filesystem;

static SimParams base_sim(std::uint64_t seed, int size) {
    SimParams s;
    s.width = size;
    s.height = size;
    s.nframes = 60;
    s.exptime_s = 45;
    s.cadence_s = 60;
    s.nstars = static_cast<int>(800.0 * size * size / (512.0 * 512.0));
    s.hot_pixels = static_cast<int>(50.0 * size * size / (512.0 * 512.0));
    s.seed = seed;
    return s;
}

static SearchParams base_search() {
    SearchParams p;
    p.vmin_pxh = 7;
    p.vmax_pxh = 44;
    return p;
}

int main(int argc, char** argv) {
    if (argc < 2) {
        std::fprintf(stderr, "usage: %s OUTDIR [--quick] [--astrometry]\n", argv[0]);
        return 1;
    }
    const fs::path out = argv[1];
    bool quick = false, only_astrometry = false;
    for (int i = 2; i < argc; ++i) {
        if (std::string(argv[i]) == "--quick") quick = true;
        if (std::string(argv[i]) == "--astrometry") only_astrometry = true;
    }
    fs::create_directories(out);

    // ---------------------------------------------------------- 1. false alarms
    if (!only_astrometry) {
        std::ofstream f(out / "false_alarms.csv");
        f << "target_false_per_run,runs,passing_candidates,observed_per_run,threshold_mean\n";
        const int runs = quick ? 8 : 40;
        for (double target : {0.02, 0.2, 1.0}) {
            int total = 0;
            double thr_sum = 0;
            for (int k = 0; k < runs; ++k) {
                SimParams s = base_sim(1000 + static_cast<std::uint64_t>(k), 256);
                s.nmovers = 0;
                const SimResult sim = simulate(s);
                SearchParams p = base_search();
                p.target_false_per_run = target;
                p.seed = 777 + static_cast<std::uint64_t>(k);
                const SearchResult r = search(sim.frames, p);
                for (const auto& c : r.candidates) total += c.ok;
                thr_sum += r.threshold_snr;
            }
            f << target << ',' << runs << ',' << total << ',' << static_cast<double>(total) / runs << ','
              << thr_sum / runs << '\n';
            std::printf("false alarms: target %.2f/run -> observed %d in %d runs (%.3f/run), mean threshold %.2f\n",
                        target, total, runs, static_cast<double>(total) / runs, thr_sum / runs);
            std::fflush(stdout);
        }
    }

    // ---------------------------------------------------------- 2. completeness
    if (!only_astrometry) {
        SimParams s = base_sim(5, quick ? 256 : 512);
        s.nmovers = 0;
        const SimResult sim = simulate(s);
        CompletenessParams cp;
        cp.snr_levels = {5, 6, 7, 8, 9, 10, 12, 15, 20};
        cp.per_level = quick ? 6 : 25;
        const auto rep = completeness(sim.frames, base_search(), cp);
        write_completeness_csv((out / "completeness.csv").string(), rep);
        std::printf("completeness: threshold %.2f (%s), 50%% at S/N %.2f\n", rep.threshold_snr,
                    rep.threshold_method.c_str(), rep.snr50);
        for (const auto& r : rep.rows)
            std::printf("  S/N %5.1f  %2d/%2d recovered  mean measured S/N %.2f\n", r.target_snr, r.recovered,
                        r.injected, r.mean_snr);
        std::fflush(stdout);
    }

    // ---------------------------------------------------------- 3. astrometry
    {
        std::ofstream f(out / "astrometry.csv");
        f << "seed,expected_snr,snr,dx_px,dy_px,pos_err_px,dvx_pxh,dvy_pxh,vel_err_pxh,rate_pxh\n";
        const int seeds = quick ? 2 : 8;
        int n = 0;
        double chi_pos = 0, chi_vel = 0;
        for (int k = 0; k < seeds; ++k) {
            SimParams s = base_sim(500 + static_cast<std::uint64_t>(k), 384);
            s.nmovers = 10;
            s.mover_mag_bright = 18.5;
            s.mover_mag_faint = 20.3;
            s.rate_min_arcsec_min = 0.25;
            s.rate_max_arcsec_min = 1.0;
            const SimResult sim = simulate(s);
            const SearchResult r = search(sim.frames, base_search());
            std::vector<SimMover> movers;
            for (const auto& t : sim.truth) movers.push_back(t.mover);
            for (const auto& m : match(r.candidates, movers, 2.0, 6.0, true)) {
                const auto& c = r.candidates[static_cast<std::size_t>(m.candidate_index)];
                const auto& t = sim.truth[static_cast<std::size_t>(m.truth_index)];
                const double dx = c.x - t.mover.x_ref, dy = c.y - t.mover.y_ref;
                const double dvx = c.vx_pxh - t.mover.vx_pxh, dvy = c.vy_pxh - t.mover.vy_pxh;
                f << k << ',' << t.expected_snr_stack << ',' << c.snr << ',' << dx << ',' << dy << ',' << c.pos_err_px
                  << ',' << dvx << ',' << dvy << ',' << c.vel_err_pxh << ','
                  << std::hypot(t.mover.vx_pxh, t.mover.vy_pxh) << '\n';
                chi_pos += (dx * dx + dy * dy) / (2 * c.pos_err_px * c.pos_err_px);
                chi_vel += (dvx * dvx + dvy * dvy) / (2 * c.vel_err_pxh * c.vel_err_pxh);
                ++n;
            }
        }
        std::printf("astrometry: %d matched movers, rms normalised position error %.2f, motion error %.2f "
                    "(1.0 = uncertainties are exact)\n",
                    n, std::sqrt(chi_pos / std::max(1, n)), std::sqrt(chi_vel / std::max(1, n)));
    }
    return 0;
}
