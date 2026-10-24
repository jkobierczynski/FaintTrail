// SPDX-License-Identifier: GPL-3.0-or-later
// Copyright 2026 Jurgen Kobierczynski
//
// End-to-end false alarms of the threshold calibration (docs/CALIBRATION.md).
//
// Many simulated sequences without moving objects (stars, hot pixels, cosmic
// rays, noise) are searched with the whole pipeline. Each search is run once
// with a deliberately low threshold, keeping its null peaks; the thresholds of
// both tail models for several requested rates are then computed from those
// null peaks, and the passing candidates above each threshold are counted.
// (Equivalent to separate searches except for the rare candidate whose ghost
// test depends on a weaker one.)
//
//   fainttrail_false_alarm_study [--runs N] [--size PX] [--frames N] [--csv FILE]
#include <cmath>
#include <cstdio>
#include <cstdlib>
#include <fstream>
#include <string>
#include <vector>

#include "../src/internal.hpp"
#include "fainttrail/fainttrail.hpp"

using namespace fainttrail;

int main(int argc, char** argv) {
    int runs = 100, size = 256, nframes = 60;
    std::string csv;
    for (int i = 1; i + 1 < argc; i += 2) {
        const std::string k = argv[i];
        if (k == "--runs") runs = std::atoi(argv[i + 1]);
        else if (k == "--size") size = std::atoi(argv[i + 1]);
        else if (k == "--frames") nframes = std::atoi(argv[i + 1]);
        else if (k == "--csv") csv = argv[i + 1];
    }
    const std::vector<double> targets = {0.02, 0.1, 0.2, 1.0};
    struct Acc {
        int passing = 0;
        double thr = 0.0;
    };
    std::vector<Acc> lin(targets.size()), gau(targets.size());
    std::ofstream out;
    if (!csv.empty()) {
        out.open(csv);
        out << "run,model,target,threshold,passing\n";
    }
    for (int r = 0; r < runs; ++r) {
        SimParams s;
        s.width = size;
        s.height = size;
        s.nframes = nframes;
        s.exptime_s = 45;
        s.cadence_s = 60;
        s.nstars = static_cast<int>(800.0 * size * size / (512.0 * 512.0));
        s.hot_pixels = static_cast<int>(50.0 * size * size / (512.0 * 512.0));
        s.nmovers = 0;
        s.seed = 20000 + static_cast<std::uint64_t>(r);
        const SimResult sim = simulate(s);
        SearchParams p;
        p.vmin_pxh = 7;
        p.vmax_pxh = 44;
        p.seed = 900 + static_cast<std::uint64_t>(r);
        p.target_false_per_run = 50.0;  // low threshold: keep every candidate that might pass
        p.max_candidates = 5000;
        const SearchResult res = search(sim.frames, p);
        std::vector<float> pooled(res.null_peaks.begin(), res.null_peaks.end());
        for (std::size_t t = 0; t < targets.size(); ++t) {
            double a = 0, b = 0, k = 0, is2 = 0;
            double T_lin = NAN, T_gau = NAN;
            if (detail::fit_null_tail(pooled, p.null_runs, a, b))
                T_lin = std::max((a - std::log(targets[t])) / b, p.min_snr_floor);
            if (detail::fit_null_tail_gauss(pooled, p.null_runs, a, k, is2, 3000.0))
                T_gau = std::max(detail::gauss_tail_level(a, k, is2, targets[t]), p.min_snr_floor);
            int pl = 0, pg = 0;
            for (const auto& c : res.candidates) {
                if (!c.ok) continue;
                pl += c.snr_coarse >= T_lin;
                pg += c.snr_coarse >= T_gau;
            }
            lin[t].passing += pl;
            lin[t].thr += T_lin;
            gau[t].passing += pg;
            gau[t].thr += T_gau;
            if (out) {
                out << r << ",linear," << targets[t] << ',' << T_lin << ',' << pl << '\n';
                out << r << ",gauss," << targets[t] << ',' << T_gau << ',' << pg << '\n';
            }
        }
        if ((r + 1) % 10 == 0) {
            std::printf("after %d runs:", r + 1);
            for (std::size_t t = 0; t < targets.size(); ++t)
                std::printf("  %.2f: linear %d, gauss %d", targets[t], lin[t].passing, gau[t].passing);
            std::printf("\n");
            std::fflush(stdout);
        }
    }
    std::printf("\n%d sequences of %d x %d px, %d frames, no moving objects\n", runs, size, size, nframes);
    std::printf("%-10s %-8s %10s %10s %12s\n", "requested", "model", "threshold", "passing", "per search");
    for (std::size_t t = 0; t < targets.size(); ++t) {
        std::printf("%-10.2f %-8s %10.3f %10d %12.3f\n", targets[t], "linear", lin[t].thr / runs, lin[t].passing,
                    static_cast<double>(lin[t].passing) / runs);
        std::printf("%-10.2f %-8s %10.3f %10d %12.3f\n", targets[t], "gauss", gau[t].thr / runs, gau[t].passing,
                    static_cast<double>(gau[t].passing) / runs);
    }
    return 0;
}
