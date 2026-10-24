// SPDX-License-Identifier: GPL-3.0-or-later
// Copyright 2026 Jurgen Kobierczynski
//
// Which null peaks would the candidate tests reject? (docs/CALIBRATION.md)
//
// For a few null runs, the top peaks are refined and tested exactly like real
// candidates. Prints the share of rejected peaks per S/N band and the threshold
// of each tail model fitted to all null peaks and to the passing ones only.
//
//   fainttrail_null_flags_study (--data DIR [--wcs-from F] | --sim) [--vmin V] [--vmax V]
//       [--runs N] [--top K]
#include <algorithm>
#include <atomic>
#include <cmath>
#include <cstdio>
#include <cstdlib>
#include <filesystem>
#include <map>
#include <string>
#include <vector>

#include "../src/common.hpp"
#include "../src/internal.hpp"
#include "fainttrail/fainttrail.hpp"

using namespace fainttrail;
namespace fs = std::filesystem;

int main(int argc, char** argv) {
    std::string data, wcs_from;
    bool sim = false;
    int runs = 3, top = 2000;
    SearchParams params;
    for (int i = 1; i < argc; ++i) {
        const std::string k = argv[i];
        auto next = [&] { return std::string(i + 1 < argc ? argv[++i] : ""); };
        if (k == "--data") data = next();
        else if (k == "--wcs-from") wcs_from = next();
        else if (k == "--sim") sim = true;
        else if (k == "--vmin") params.vmin_pxh = std::atof(next().c_str());
        else if (k == "--vmax") params.vmax_pxh = std::atof(next().c_str());
        else if (k == "--runs") runs = std::atoi(next().c_str());
        else if (k == "--top") top = std::atoi(next().c_str());
    }
    tune_allocator_for_images();
    std::vector<Frame> frames;
    if (sim) {
        SimParams sp;
        sp.width = 256;
        sp.height = 256;
        sp.nframes = 60;
        sp.nstars = 200;
        sp.hot_pixels = 12;
        sp.nmovers = 0;
        sp.seed = 20001;
        frames = simulate(sp).frames;
        params.vmin_pxh = 7;
        params.vmax_pxh = 44;
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
    const int nt = detail::thread_count(params.threads);
    std::vector<float> all, passing, noise;  // noise: all peaks except artefacts (edge, detector, single)
    std::map<std::string, int> flag_count;
    const double bands[] = {3, 4, 5, 6, 7, 8, 10, 15, 1e9};
    int band_all[8] = {0}, band_bad[8] = {0};
    const auto orders = detail::null_time_orders(p.t_h, runs, params.seed);
    for (int r = 0; r < runs; ++r) {
        const auto& t = orders[static_cast<std::size_t>(r)];
        const Prepared q = detail::with_times(p, t, params);
        const auto null = detail::run_pass(*backend, q, grid, t, params, {}, 0, 1, "null");
        auto peaks = detail::find_peaks(null, 3.0f, nms, edge, 20000);
        for (const auto& pk : peaks) all.push_back(pk.snr);
        for (std::size_t k = static_cast<std::size_t>(top); k < peaks.size(); ++k) noise.push_back(peaks[k].snr);
        if (static_cast<int>(peaks.size()) > top) peaks.resize(static_cast<std::size_t>(top));
        std::vector<detail::RefinedCandidate> rc(peaks.size());
#pragma omp parallel for num_threads(nt) schedule(dynamic)
        for (int k = 0; k < static_cast<int>(peaks.size()); ++k)
            rc[static_cast<std::size_t>(k)] =
                detail::refine_with_contrib(q, peaks[static_cast<std::size_t>(k)],
                                            grid[static_cast<std::size_t>(peaks[static_cast<std::size_t>(k)].vidx)], dv, params);
        SearchResult tmp;
        detail::finish_candidates(q, rc, tmp);  // ghost test and final flags, as for real candidates
        // finish_candidates drops duplicates and reorders; match back by coarse S/N
        // artefact test on every tested peak (finish_candidates drops duplicates; use rc directly)
        for (const auto& r2 : rc) {
            const std::string& f = r2.c.flags;
            const bool artefact = f.find("edge") != std::string::npos || f.find("detector") != std::string::npos ||
                                  f.find("single") != std::string::npos;
            if (!artefact) noise.push_back(static_cast<float>(r2.c.snr_coarse));
        }
        for (const auto& c : tmp.candidates) {
            int b = 0;
            while (c.snr_coarse >= bands[b + 1]) ++b;
            ++band_all[b];
            if (!c.ok) {
                ++band_bad[b];
                std::size_t s = 0;
                while (s <= c.flags.size()) {
                    const std::size_t e = c.flags.find(',', s);
                    ++flag_count[c.flags.substr(s, e == std::string::npos ? std::string::npos : e - s)];
                    if (e == std::string::npos) break;
                    s = e + 1;
                }
            } else {
                passing.push_back(static_cast<float>(c.snr_coarse));
            }
        }
        std::printf("null run %d: top %zu peaks tested\n", r + 1, rc.size());
        std::fflush(stdout);
    }
    std::printf("\nrejected share of the tested null peaks by coarse S/N:\n");
    for (int b = 0; b < 8; ++b)
        if (band_all[b])
            std::printf("  S/N %5.1f-%-5.1f  %6d tested  %5.1f %% rejected\n", bands[b], std::min(bands[b + 1], 99.0),
                        band_all[b], 100.0 * band_bad[b] / band_all[b]);
    std::printf("flags of rejected null peaks:");
    for (const auto& [f, c] : flag_count) std::printf("  %s %d", f.c_str(), c);
    std::printf("\n");
    // thresholds: all null peaks vs passing ones (passing = tested and ok; for the
    // gauss fit only the top `top` per run are used, which are all tested)
    const double target = params.target_false_per_run;
    double a, b, k, is2;
    auto lin = [&](const std::vector<float>& v) {
        return detail::fit_null_tail(v, runs, a, b) ? (a - std::log(target)) / b : NAN;
    };
    auto gau = [&](const std::vector<float>& v, double topn) {
        return detail::fit_null_tail_gauss(v, runs, a, k, is2, topn) ? detail::gauss_tail_level(a, k, is2, target) : NAN;
    };
    std::printf("threshold, linear model : all peaks %.3f, passing peaks %.3f\n", lin(all), lin(passing));
    std::printf("threshold, gauss model  : all peaks %.3f (top 3000/run), passing peaks %.3f (top %d/run tested)\n",
                gau(all, 3000.0), gau(passing, 0.9 * top), top);
    std::printf("threshold, without artefacts (edge/detector/single): gauss %.3f, linear %.3f\n", gau(noise, 3000.0),
                lin(noise));
    return 0;
}
