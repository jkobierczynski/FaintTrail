// SPDX-License-Identifier: GPL-3.0-or-later
// Copyright 2026 Jurgen Kobierczynski
//
// Unit and end-to-end tests. Plain C++, no framework needed: each test is a
// function; a failed CHECK prints the location and marks the test as failed.
#include <algorithm>
#include <array>
#include <cmath>
#include <cstdio>
#include <filesystem>
#include <fstream>
#include <functional>
#include <iomanip>
#include <string>
#include <vector>

#include "fainttrail/fainttrail.hpp"
#include "backend_cpu.hpp"
#include "cuda_plan.hpp"

extern "C" int run_c_api_test(void);

using namespace fainttrail;
namespace fs = std::filesystem;

static int g_failed_checks = 0;

static unsigned long long g_rng = 1;
static void detail_rng_seed(unsigned long long s) { g_rng = s; }
static double detail_rng_uniform() {
    g_rng ^= g_rng << 13; g_rng ^= g_rng >> 7; g_rng ^= g_rng << 17;
    return (static_cast<double>(g_rng >> 11) + 0.5) * (1.0 / 9007199254740992.0);
}
static double detail_rng_normal() {
    return std::sqrt(-2.0 * std::log(detail_rng_uniform())) * std::cos(6.283185307179586 * detail_rng_uniform());
}

#define CHECK(cond)                                                                 \
    do {                                                                            \
        if (!(cond)) {                                                              \
            std::printf("    CHECK failed: %s  (%s:%d)\n", #cond, __FILE__, __LINE__); \
            ++g_failed_checks;                                                      \
        }                                                                           \
    } while (0)

#define CHECK_NEAR(a, b, tol)                                                                     \
    do {                                                                                          \
        const double va_ = (a), vb_ = (b);                                                        \
        if (!(std::fabs(va_ - vb_) <= (tol))) {                                                   \
            std::printf("    CHECK_NEAR failed: %s = %g, %s = %g, tol %g  (%s:%d)\n", #a, va_, #b, vb_, \
                        static_cast<double>(tol), __FILE__, __LINE__);                            \
            ++g_failed_checks;                                                                    \
        }                                                                                         \
    } while (0)

static fs::path temp_dir() {
    fs::path d = fs::temp_directory_path() / "fainttrail_tests";
    fs::create_directories(d);
    return d;
}

// A small, quick simulation used by several tests.
static SimParams small_sim(std::uint64_t seed) {
    SimParams s;
    s.width = 192;
    s.height = 192;
    s.nframes = 40;
    s.exptime_s = 45;
    s.cadence_s = 60;
    s.nstars = 200;
    s.hot_pixels = 15;
    s.cosmic_rays_per_frame = 2;
    s.seed = seed;
    return s;
}

static SearchParams quick_search() {
    SearchParams p;
    p.vmin_pxh = 8;
    p.vmax_pxh = 30;
    p.null_runs = 2;
    return p;
}

// ------------------------------------------------------------------ tests
static void test_fits_roundtrip() {
    Image img(37, 23);
    for (int y = 0; y < img.height; ++y)
        for (int x = 0; x < img.width; ++x) img.at(x, y) = static_cast<float>(x * 1000 + y) * 0.25f - 7.0f;
    Frame f;
    f.image = img;
    f.mjd_mid = 61337.912345678;
    f.exptime_s = 30;
    f.wcs = Wcs::make_tan(123.4, -45.6, 18.5, 12.0, 1.25, 10.0);
    FitsHeader extra;
    extra.set("OBJECT", std::string("O'Brien test"));
    extra.set_int("FOO", 42);
    extra.add_history("written by the fainttrail unit tests");
    const std::string path = (temp_dir() / "roundtrip.fits").string();
    write_frame(path, f, extra);

    FitsImage back = read_fits(path);
    CHECK(back.image.width == 37 && back.image.height == 23);
    bool same = back.image.data.size() == img.data.size();
    for (std::size_t i = 0; same && i < img.data.size(); ++i) same = back.image.data[i] == img.data[i];
    CHECK(same);
    CHECK(back.header.get_string("OBJECT") == "O'Brien test");
    CHECK(back.header.get_int("FOO", 0) == 42);
    CHECK(fs::file_size(path) % 2880 == 0);

    const Frame g = frame_from_fits(std::move(back), path);
    CHECK_NEAR(g.mjd_mid, f.mjd_mid, 1e-9);
    CHECK_NEAR(g.exptime_s, 30.0, 1e-12);
    CHECK(g.wcs.valid);
    CHECK_NEAR(g.wcs.cd12, f.wcs.cd12, 1e-15);
    CHECK_NEAR(g.wcs.crpix1, 18.5, 1e-12);
}

static void test_time_conversion() {
    CHECK_NEAR(mjd_from_iso("1858-11-17T00:00:00"), 0.0, 1e-12);
    CHECK_NEAR(mjd_from_iso("2000-01-01T12:00:00"), 51544.5, 1e-12);
    CHECK_NEAR(mjd_from_iso("2026-10-24"), 61337.0, 1e-12);
    const double m = 61337.9123456;
    CHECK_NEAR(mjd_from_iso(iso_from_mjd(m, 3)), m, 1e-8);
    CHECK(iso_from_mjd(51544.5, 0) == "2000-01-01T12:00:00");
}

static void test_wcs_roundtrip() {
    const Wcs w = Wcs::make_tan(359.9, 89.0, 100.5, 80.5, 1.5, 33.0);
    double worst = 0.0;
    for (double y = -50; y <= 250; y += 37)
        for (double x = -50; x <= 250; x += 41) {
            double ra, dec, x2, y2;
            CHECK(w.pix2sky(x, y, ra, dec));
            CHECK(w.sky2pix(ra, dec, x2, y2));
            worst = std::max(worst, std::hypot(x2 - x, y2 - y));
        }
    CHECK(worst < 1e-6);
    CHECK_NEAR(w.pixel_scale_arcsec(), 1.5, 1e-9);
    double ra, dec;
    w.pix2sky(99.5, 79.5, ra, dec);  // the reference pixel, 0-based
    CHECK_NEAR(ra, 359.9, 1e-9);
    CHECK_NEAR(dec, 89.0, 1e-9);
}

static void test_noise_normalisation() {
    SimParams s = small_sim(11);
    s.nmovers = 0;
    s.nmovers = 0;
    const SimResult sim = simulate(s);
    const Prepared p = prepare(sim.frames, quick_search());
    for (int i = 0; i < p.nframes; i += 7) {
        double sum = 0, sum2 = 0;
        int n = 0;
        for (int y = 10; y < p.height - 10; ++y)
            for (int x = 10; x < p.width - 10; ++x) {
                if (p.mask.at(x, y) == 0.0f) continue;
                const double v = p.filtered[static_cast<std::size_t>(i)].at(x, y);
                sum += v;
                sum2 += v * v;
                ++n;
            }
        const double mean = sum / n, sd = std::sqrt(sum2 / n - mean * mean);
        CHECK(std::fabs(mean) < 0.1);
        CHECK(sd > 0.9 && sd < 1.15);
    }
}

static void test_fwhm_estimate() {
    for (double fwhm : {2.2, 3.0, 4.5}) {
        SimParams s = small_sim(5);
        s.fwhm_px = fwhm;
        s.fwhm_jitter = 0.0;
        s.nmovers = 0;
        s.nmovers = 0;
        s.nframes = 12;
        const SimResult sim = simulate(s);
        const Prepared p = prepare(sim.frames, quick_search());
        CHECK(std::fabs(p.fwhm_px / fwhm - 1.0) < 0.12);
    }
}

static void test_known_motion_stack() {
    SimParams s = small_sim(21);
    SimMover m;
    m.id = "K";
    m.x_ref = 96;
    m.y_ref = 90;
    m.vx_pxh = 20;
    m.vy_pxh = -9;
    m.mag = 19.8;
    const SimResult sim = simulate(s, {m});
    const Prepared p = prepare(sim.frames, quick_search());
    const Image cut = stack_cutout(p, m.x_ref, m.y_ref, m.vx_pxh, m.vy_pxh, 6, true, 0);
    const double snr = cut.at(6, 6);
    const double expected = sim.truth[0].expected_snr_stack;
    CHECK(snr > 0.7 * expected);
    CHECK(snr < 1.2 * expected);
    // a wrong motion must not stack up
    const Image wrong = stack_cutout(p, m.x_ref, m.y_ref, -m.vx_pxh, m.vy_pxh, 6, true, 0);
    CHECK(wrong.at(6, 6) < 0.4 * snr);
}

static void test_blind_search_recovers() {
    SimParams s = small_sim(3);
    std::vector<SimMover> movers;
    const double pos[4][4] = {{40, 50, 22, 10}, {150, 60, -15, 12}, {60, 140, 9, -14}, {140, 150, -12, -18}};
    const double mags[4] = {18.8, 19.1, 19.3, 19.5};
    for (int i = 0; i < 4; ++i) {
        SimMover m;
        m.id = "M" + std::to_string(i + 1);
        m.x_ref = pos[i][0];
        m.y_ref = pos[i][1];
        m.vx_pxh = pos[i][2];
        m.vy_pxh = pos[i][3];
        m.mag = mags[i];
        movers.push_back(m);
    }
    const SimResult sim = simulate(s, movers);
    const SearchResult res = search(sim.frames, quick_search());
    CHECK(res.threshold_method == "null-fit gauss" || res.threshold_method == "null-fit linear" || res.threshold_method == "null-quantile");
    const auto matched = match(res.candidates, movers, 1.0, 2.0 * res.dv_pxh + 1.0, true);
    CHECK(matched.size() == movers.size());
    for (const auto& mr : matched) {
        // errors within 4 sigma of the reported uncertainties
        const auto& c = res.candidates[static_cast<std::size_t>(mr.candidate_index)];
        CHECK(mr.dpos_px < 4.0 * std::sqrt(2.0) * c.pos_err_px);
        CHECK(mr.dv_pxh < 4.0 * std::sqrt(2.0) * c.vel_err_pxh + 0.3);
    }
    int ok = 0;
    for (const auto& c : res.candidates) ok += c.ok;
    CHECK(ok == static_cast<int>(movers.size()));  // no false positives
    // sky coordinates and rates follow from the WCS
    for (const auto& c : res.candidates) {
        if (!c.ok) continue;
        CHECK(c.ra_deg > 29.0 && c.ra_deg < 31.0);
        CHECK(c.rate_arcsec_min > 0.1);
    }
}

static void test_null_has_no_detections() {
    for (std::uint64_t seed : {101ULL, 202ULL}) {
        SimParams s = small_sim(seed);
        s.nmovers = 0;
        const SimResult sim = simulate(s);
        const SearchResult res = search(sim.frames, quick_search());
        int ok = 0;
        for (const auto& c : res.candidates) ok += c.ok;
        CHECK(ok == 0);
        CHECK(res.threshold_snr >= quick_search().min_snr_floor);
    }
}

static void test_completeness() {
    SimParams s = small_sim(7);
    s.nmovers = 0;
    const SimResult sim = simulate(s);
    CompletenessParams cp;
    cp.snr_levels = {4.0, 12.0};
    cp.per_level = 5;
    const auto rep = completeness(sim.frames, quick_search(), cp);
    CHECK(rep.rows.size() == 2);
    CHECK(rep.rows[1].injected >= 4);
    CHECK(rep.rows[1].recovered >= rep.rows[1].injected - 1);
    CHECK(rep.rows[0].recovered <= rep.rows[1].recovered);
    CHECK(rep.baseline_ok == 0);
}

static void test_csv_roundtrip() {
    SearchResult r;
    Candidate c;
    c.id = 1;
    c.ok = true;
    c.flags = "vedge";
    c.x = 10.25;
    c.y = 20.5;
    c.vx_pxh = 3.5;
    c.vy_pxh = -4.25;
    c.snr = 9.5;
    r.candidates.push_back(c);
    const std::string path = (temp_dir() / "cands.csv").string();
    write_candidates_csv(path, r);
    const auto back = read_candidates_csv(path);
    CHECK(back.size() == 1);
    if (!back.empty()) {
        CHECK_NEAR(back[0].x, 10.25, 1e-3);
        CHECK_NEAR(back[0].vy_pxh, -4.25, 1e-3);
        CHECK(back[0].ok);
        CHECK(back[0].flags == "vedge");
    }
}

static void test_cancellation() {
    SimParams s = small_sim(9);
    s.nframes = 10;
    const SimResult sim = simulate(s);
    bool threw = false;
    try {
        search(sim.frames, quick_search(), [](double f, const std::string&) { return f < 0.3; });
    } catch (const Cancelled&) {
        threw = true;
    }
    CHECK(threw);
}


// A solved full-resolution image and 2x2-binned, shifted, rotated and/or mirrored
// frames of the same field: the transferred WCS must reproduce the true sky positions.
static void test_wcs_transfer() {
    SimParams s = small_sim(17);
    s.width = 400;
    s.height = 300;
    s.nframes = 1;
    s.nmovers = 0;
    s.nstars = 400;
    s.pixscale_arcsec = 0.8;
    s.cosmic_rays_per_frame = 0;
    s.hot_pixels = 0;
    const SimResult full = simulate(s);
    const Image& src = full.frames[0].image;
    const Wcs truth = full.frames[0].wcs;  // the "plate solution" of the source image

    struct Case { const char* name; int rot; bool mirror; };
    for (const Case cs : {Case{"binned", 0, false}, Case{"binned+rot90", 1, false}, Case{"binned+rot180+mirror", 2, true}}) {
        // build the target: 2x2 bin, then rotate/mirror, then shift by (3,2) pixels
        const int bw = src.width / 2, bh = src.height / 2;
        Image b(bw, bh);
        for (int y = 0; y < bh; ++y)
            for (int x = 0; x < bw; ++x)
                b.at(x, y) = src.at(2 * x, 2 * y) + src.at(2 * x + 1, 2 * y) + src.at(2 * x, 2 * y + 1) + src.at(2 * x + 1, 2 * y + 1);
        auto xform = [&](double x, double y, double& tx, double& ty, int W, int H) {
            if (cs.mirror) y = H - 1 - y;
            for (int r = 0; r < cs.rot; ++r) { const double nx = H - 1 - y, ny = x; x = nx; y = ny; std::swap(W, H); }
            tx = x + 3; ty = y + 2;
        };
        const int TW = (cs.rot % 2 ? bh : bw) + 6, TH = (cs.rot % 2 ? bw : bh) + 4;
        Image t(TW, TH, static_cast<float>(4 * 540.0));
        for (int y = 0; y < bh; ++y)
            for (int x = 0; x < bw; ++x) {
                double tx, ty;
                xform(x, y, tx, ty, bw, bh);
                t.at(static_cast<int>(std::lround(tx)), static_cast<int>(std::lround(ty))) = b.at(x, y);
            }
        TransferReport rep;
        Wcs w;
        try {
            w = transfer_wcs(truth, src, t, &rep);
        } catch (const std::exception& e) {
            std::printf("    %s: %s\n", cs.name, e.what());
            CHECK(false);
            continue;
        }
        CHECK(rep.ok);
        CHECK(rep.matched > 20);
        CHECK(rep.mirrored == cs.mirror);
        CHECK_NEAR(rep.scale, 0.5, 0.01);
        // compare sky positions of a few source pixels mapped into the target
        double worst = 0;
        for (double sx : {40.5, 200.5, 360.5})
            for (double sy : {30.5, 150.5, 270.5}) {
                double ra0, de0, ra1, de1, tx, ty;
                truth.pix2sky(sx, sy, ra0, de0);
                // binned pixel centre of source pixel (sx,sy) is ((sx-0.5)/2, (sy-0.5)/2)
                xform((sx - 0.5) / 2.0, (sy - 0.5) / 2.0, tx, ty, bw, bh);
                w.pix2sky(tx, ty, ra1, de1);
                worst = std::max(worst, std::hypot((ra1 - ra0) * std::cos(de0 * 3.14159265358979 / 180), de1 - de0) * 3600.0);
            }
        if (worst > 0.5) std::printf("    %s: worst sky error %.3f arcsec\n", cs.name, worst);
        CHECK(worst < 0.5);  // well below one binned pixel (1.6")
    }
}


// The CPU backend must follow the stacking contract exactly: compare it with the
// naive implementation on awkward inputs (odd sizes, shifts beyond the frame,
// half-integer shift values, partial coverage, all-equal frames for tie-breaking).
static void test_backend_contract() {
    detail_rng_seed(4242);
    auto make_pass_data = [](int w, int h, int n, std::vector<Image>& frames, std::vector<double>& t,
                             std::vector<std::pair<double, double>>& vel, bool constant) {
        frames.assign(static_cast<std::size_t>(n), Image(w, h));
        t.resize(static_cast<std::size_t>(n));
        for (int i = 0; i < n; ++i) {
            t[static_cast<std::size_t>(i)] = (i - (n - 1) / 2.0) * 0.25;  // hours, symmetric, includes x.5 shifts
            for (auto& v : frames[static_cast<std::size_t>(i)].data) v = constant ? 1.0f : static_cast<float>(detail_rng_normal());
        }
        vel.clear();
        for (int j = -6; j <= 6; ++j)
            for (int k = -6; k <= 6; ++k) vel.emplace_back(j * 2.0, k * 2.0);  // up to 12 px/h -> shifts up to 9 px
        vel.emplace_back(400.0, 0.0);  // shifts far outside the frame: never covered enough
    };
    for (const bool constant : {false, true}) {
        std::vector<Image> frames;
        std::vector<double> t;
        std::vector<std::pair<double, double>> vel;
        make_pass_data(37, 23, 7, frames, t, vel, constant);
        StackPass pass;
        pass.frames = &frames;
        pass.t_h = &t;
        pass.velocities = &vel;
        pass.min_count = 5;
        const StackResult ref = run_pass_naive(pass);
        for (int threads : {1, 2, 3}) {
            pass.threads = threads;
            auto be = make_backend("cpu");
            const StackResult got = be->run(pass);
            const BackendCheck chk = check_backend(pass, ref, got);
            if (!chk.ok) std::printf("    cpu (%d threads, constant=%d): %s\n", threads, constant, chk.message.c_str());
            CHECK(chk.ok);
            // the CPU reference must agree bit for bit, ties included
            CHECK(got.best_v == ref.best_v);
            CHECK(got.max_snr.data == ref.max_snr.data);
        }
        // block sizes must not change anything: odd bands, single rows, bands
        // taller than the image, motion blocks that do not divide the grid
        for (const auto& blk : {detail::CpuBlocking{1, 1}, detail::CpuBlocking{3, 5}, detail::CpuBlocking{64, 64},
                                detail::CpuBlocking{7, 2}}) {
            pass.threads = 2;
            detail::CpuBackend be(blk);
            const StackResult got = be.run(pass);
            CHECK(got.best_v == ref.best_v);
            CHECK(got.max_snr.data == ref.max_snr.data);
        }
        const StackResult simple = run_pass_simple(pass);
        CHECK(simple.best_v == ref.best_v);
        CHECK(simple.max_snr.data == ref.max_snr.data);
        if (constant) {
            // all velocities tie wherever they are fully covered: lowest index wins
            // pixel (18,11) is covered by all 7 frames even for velocity 0 (shifts up to 9 px)
            CHECK(ref.best_v[static_cast<std::size_t>(11) * 37 + 18] == 0);
            CHECK_NEAR(ref.max_snr.at(18, 11), std::sqrt(7.0), 1e-5);
        }
    }
    // the checker itself must notice a wrong result
    {
        std::vector<Image> frames;
        std::vector<double> t;
        std::vector<std::pair<double, double>> vel;
        make_pass_data(20, 15, 6, frames, t, vel, false);
        StackPass pass;
        pass.frames = &frames;
        pass.t_h = &t;
        pass.velocities = &vel;
        pass.min_count = 4;
        const StackResult ref = run_pass_naive(pass);
        StackResult bad = ref;
        bad.max_snr.data[37] += 0.01f;
        CHECK(!check_backend(pass, ref, bad).ok);
        StackResult bad2 = ref;
        for (auto& v : bad2.best_v)
            if (v >= 0) { v = (v + 1) % static_cast<int>(vel.size()); break; }
        CHECK(!check_backend(pass, ref, bad2).ok);
    }
    CHECK(!backend_names().empty() && backend_names().front() == (cuda_backend_available() ? "cuda" : "cpu"));
    bool threw = false;
    try { make_backend("no-such-backend"); } catch (const std::runtime_error&) { threw = true; }
    CHECK(threw);
}

// Explicit shift tables (D38): the straight-line table reproduces the default
// exactly; arbitrary (curved, jumping, far out of the frame) tables give
// bit-identical results in every implementation and block size.
static void test_shift_tables() {
    detail_rng_seed(777);
    const int w = 41, h = 29, n = 9;
    std::vector<Image> frames(static_cast<std::size_t>(n), Image(w, h));
    std::vector<double> t(static_cast<std::size_t>(n));
    for (int i = 0; i < n; ++i) {
        t[static_cast<std::size_t>(i)] = (i - 4) * 0.3;
        for (auto& v : frames[static_cast<std::size_t>(i)].data) v = static_cast<float>(detail_rng_normal());
    }
    std::vector<std::pair<double, double>> vel;
    for (int j = -5; j <= 5; ++j)
        for (int k = -4; k <= 4; ++k) vel.emplace_back(j * 3.0, k * 2.5);
    StackPass base;
    base.frames = &frames;
    base.t_h = &t;
    base.velocities = &vel;
    base.min_count = 6;
    const StackResult ref = run_pass_naive(base);
    // 1. the straight-line table is the default, bit for bit
    const ShiftTable lin = linear_shift_table(t, vel);
    StackPass tab = base;
    tab.shifts = &lin;
    for (const StackResult& r : {run_pass_naive(tab), run_pass_simple(tab), make_backend("cpu")->run(tab)}) {
        CHECK(r.max_snr.data == ref.max_snr.data);
        CHECK(r.best_v == ref.best_v);
    }
    // 2. an arbitrary table: curved tracks, random jumps, shifts far outside the frame
    ShiftTable cur;
    cur.resize(150, n);
    for (int j = 0; j < cur.motions; ++j)
        for (int i = 0; i < n; ++i) {
            const double ti = t[static_cast<std::size_t>(i)];
            const double ax = (j % 7) - 3.0, ay = (j % 5) - 2.0;
            cur.x(j, i) = static_cast<int>(std::lround((j % 13 - 6) * ti + ax * ti * ti)) + (j == 77 ? 500 : 0);
            cur.y(j, i) = static_cast<int>(std::lround((j % 11 - 5) * ti + ay * ti * ti)) + (j % 17 == 0 && i == 3 ? 4 : 0);
        }
    StackPass free = base;
    free.velocities = nullptr;  // not needed with a table
    free.t_h = nullptr;
    free.shifts = &cur;
    const StackResult naive = run_pass_naive(free);
    const StackResult simple = run_pass_simple(free);
    CHECK(simple.max_snr.data == naive.max_snr.data && simple.best_v == naive.best_v);
    for (const auto& blk : {detail::CpuBlocking{1, 1}, detail::CpuBlocking{5, 7}, detail::CpuBlocking{64, 64}})
        for (int threads : {1, 3}) {
            free.threads = threads;
            detail::CpuBackend be(blk);
            const StackResult r = be.run(free);
            CHECK(r.max_snr.data == naive.max_snr.data);
            CHECK(r.best_v == naive.best_v);
        }
    CHECK(check_backend(free, naive, make_backend("cpu")->run(free)).ok);
    // 3. malformed passes are refused
    ShiftTable bad = cur;
    bad.frames = n - 1;
    StackPass wrong = free;
    wrong.shifts = &bad;
    bool threw = false;
    try {
        run_pass_naive(wrong);
    } catch (const std::runtime_error&) {
        threw = true;
    }
    CHECK(threw);
    StackPass mismatch = free;
    mismatch.velocities = &vel;  // 99 velocities, 150 motions in the table
    threw = false;
    try {
        make_backend("cpu")->run(mismatch);
    } catch (const std::runtime_error&) {
        threw = true;
    }
    CHECK(threw);
}

// The CUDA algorithm (row bands, motion chunks, per-pixel kernel code) run on
// the host: bit-identical to the naive contract for every band and chunk size.
// When a CUDA device is usable, the real backend is checked the same way.
static void test_cuda_algorithm() {
    detail_rng_seed(99);
    const int w = 37, h = 23, n = 7;
    std::vector<Image> frames(static_cast<std::size_t>(n), Image(w, h));
    std::vector<double> t(static_cast<std::size_t>(n));
    for (int i = 0; i < n; ++i) {
        t[static_cast<std::size_t>(i)] = (i - 3) * 0.25;
        for (auto& v : frames[static_cast<std::size_t>(i)].data) v = static_cast<float>(detail_rng_normal());
    }
    frames[2].data[100] = std::nanf("");  // NaN propagates identically
    std::vector<std::pair<double, double>> vel;
    for (int j = -6; j <= 6; ++j)
        for (int k = -6; k <= 6; ++k) vel.emplace_back(j * 2.0, k * 2.0);
    vel.emplace_back(400.0, 0.0);
    vel.emplace_back(0.0, -400.0);  // 171 motions: not a multiple of the motion block
    ShiftTable cur;
    cur.resize(150, n);
    for (int j = 0; j < cur.motions; ++j)
        for (int i = 0; i < n; ++i) {
            const double ti = t[static_cast<std::size_t>(i)];
            cur.x(j, i) = static_cast<int>(std::lround((j % 13 - 6) * ti + ((j % 7) - 3.0) * ti * ti)) + (j == 77 ? 500 : 0);
            cur.y(j, i) = static_cast<int>(std::lround((j % 11 - 5) * ti + ((j % 5) - 2.0) * ti * ti)) - (j == 91 ? 900 : 0);
        }
    const double row = 4.0 * w * n;  // device bytes per frame row
    for (const bool table : {false, true}) {
        StackPass pass;
        pass.frames = &frames;
        pass.t_h = &t;
        pass.velocities = table ? nullptr : &vel;
        pass.shifts = table ? &cur : nullptr;
        pass.min_count = 5;
        const StackResult ref = run_pass_naive(pass);
        for (const double budget : {1.0, 20 * row, 30 * row, 1e9})       // 1-row bands ... one band
            for (const double adds : {1.0, 3.0 * w * h * n, 1e12}) {      // 8 motions per launch ... all
                const detail::CudaPlan plan = detail::make_cuda_plan(pass, budget, adds);
                CHECK(plan.chunk % detail::cuda_core::kMotionBlock == 0);
                for (const auto& b : plan.bands) CHECK(b.r0 >= 0 && b.r1 <= h && b.r1 > b.r0 && b.y1 > b.y0);
                const StackResult got = detail::emulate_cuda_pass(pass, budget, adds);
                CHECK(got.max_snr.data.size() == ref.max_snr.data.size());
                bool same = got.best_v == ref.best_v;
                for (std::size_t k = 0; same && k < ref.max_snr.data.size(); ++k) {
                    const float a = got.max_snr.data[k], b = ref.max_snr.data[k];
                    same = (a == b) || (std::isnan(a) && std::isnan(b));
                }
                if (!same) std::printf("    emulation differs: table=%d budget=%g adds=%g\n", table, budget, adds);
                CHECK(same);
            }
        CHECK(detail::make_cuda_plan(pass, 1.0, 1.0).bands.size() == static_cast<std::size_t>(h));
        CHECK(detail::make_cuda_plan(pass, 1e9, 1e12).bands.size() == 1);
        if (cuda_backend_available()) {
            for (const char* name : {"cuda", "cuda:max_mb=0.01,launch_adds=1"}) {
                const StackResult got = make_backend(name)->run(pass);
                if (!check_backend(pass, ref, got).ok) std::printf("    %s: %s\n", name, check_backend(pass, ref, got).message.c_str());
                CHECK(got.best_v == ref.best_v);
                CHECK(check_backend(pass, ref, got).ok);
            }
        }
    }
    // options
    const detail::CudaOptions o = detail::parse_cuda_options("device=1,max_mb=64");
    CHECK(o.device == 1 && o.max_mb == 64);
    bool threw = false;
    try { detail::parse_cuda_options("speed=11"); } catch (const std::runtime_error&) { threw = true; }
    CHECK(threw);
    if (!cuda_backend_available()) {
        threw = false;
        try { make_backend("cuda"); } catch (const std::runtime_error&) { threw = true; }
        CHECK(threw);
    }
}

// Unregistered input (a drifting mount) is recognised and refused; registered
// frames pass and keep working.
static void test_registration_check() {
    SimParams sp = small_sim(31);
    sp.width = 256;
    sp.height = 192;
    sp.nframes = 20;
    sp.nstars = 150;
    const SimResult sim = simulate(sp);
    std::vector<const Image*> ims;
    std::vector<double> t;
    for (const auto& f : sim.frames) {
        ims.push_back(&f.image);
        t.push_back((f.mjd_mid - sim.frames[0].mjd_mid) * 24.0);
    }
    const DriftCheck ok = check_registration(ims, t);
    CHECK(ok.measured);
    CHECK(ok.registered);
    // shift frame i by round(0.9 i) px in x and round(0.3 i) in y, as a drifting mount would
    std::vector<Frame> drift = sim.frames;
    for (std::size_t i = 0; i < drift.size(); ++i) {
        const int dx = static_cast<int>(std::lround(0.9 * i)), dy = static_cast<int>(std::lround(0.3 * i));
        Image out(sp.width, sp.height, 0.0f);
        for (int y = 0; y < sp.height; ++y)
            for (int x = 0; x < sp.width; ++x) {
                const int sx = x - dx, sy = y - dy;
                out.at(x, y) = (sx >= 0 && sy >= 0 && sx < sp.width && sy < sp.height) ? sim.frames[i].image.at(sx, sy) : 300.0f;
            }
        drift[i].image = std::move(out);
    }
    ims.clear();
    for (const auto& f : drift) ims.push_back(&f.image);
    const DriftCheck bad = check_registration(ims, t);
    CHECK(bad.measured && !bad.registered);
    CHECK_NEAR(bad.dx, std::lround(0.9 * 19), 0.6);
    CHECK_NEAR(bad.dy, std::lround(0.3 * 19), 0.6);
    CHECK(!bad.message.empty());
    SearchParams p = quick_search();
    p.align = "check";
    bool threw = false;
    try { prepare(drift, p); } catch (const std::runtime_error&) { threw = true; }
    CHECK(threw);
    p.align = "off";
    {
        const Prepared q = prepare(drift, p);
        CHECK(!q.alignment.applied);
    }
    // aligned: back on the middle frame's grid, with the offsets to the detector recorded
    std::vector<Frame> al = drift;
    const AlignmentReport rep = register_frames(al);
    CHECK(rep.applied && rep.dropped == 0 && al.size() == drift.size());
    const std::size_t ref = static_cast<std::size_t>(rep.reference);
    CHECK(ref == drift.size() / 2);
    for (std::size_t i = 0; i < al.size(); ++i) {
        const double ex = std::lround(0.9 * i) - std::lround(0.9 * ref), ey = std::lround(0.3 * i) - std::lround(0.3 * ref);
        CHECK_NEAR(rep.frames[i].dx, ex, 0.2);  // frame i is displaced by (ex, ey) from the reference
        CHECK_NEAR(rep.frames[i].dy, ey, 0.2);
        CHECK(rep.frames[i].rms_px < 0.3);
        CHECK(al[i].has_reg_offset);
        CHECK_NEAR(al[i].reg_dx, ex, 0.2);  // detector = registered + offset
    }
    ims.clear();
    for (const auto& f : al) ims.push_back(&f.image);
    CHECK(check_registration(ims, t).registered);
    // pixel values: an aligned frame is its simulated original, displaced like the reference
    {
        const std::size_t i = 3;
        const int rx = static_cast<int>(std::lround(0.9 * ref)), ry = static_cast<int>(std::lround(0.3 * ref));
        double maxdiff = 0, peak = 0;
        for (int y = 30; y < sp.height - 30; ++y)
            for (int x = 30; x < sp.width - 30; ++x) {
                const float v = al[i].image.at(x, y), o = sim.frames[i].image.at(x - rx, y - ry);
                if (!std::isfinite(v)) continue;
                maxdiff = std::max(maxdiff, std::fabs(static_cast<double>(v) - o));
                peak = std::max(peak, static_cast<double>(o));
            }
        CHECK(maxdiff < 0.1 * peak);  // a few hundredths of a pixel of fit error on the brightest star
    }
    // "auto" aligns inside prepare() and the search then works on aligned data
    p.align = "auto";
    const Prepared q = prepare(drift, p);
    CHECK(q.alignment.applied);
    CHECK(q.reg_dx.size() == drift.size());
}

// End to end: a drifting sequence, aligned automatically, finds the same mover
// with nearly the same S/N as the registered original.
static void test_aligned_search() {
    SimParams sp = small_sim(57);
    sp.width = 320;
    sp.height = 256;
    sp.nstars = 600;
    SimMover m;
    m.id = "M";
    m.x_ref = 150;
    m.y_ref = 120;
    m.vx_pxh = 18;
    m.vy_pxh = 7;
    m.mag = 19.6;
    const SimResult sim = simulate(sp, {m});
    SearchParams p = quick_search();
    p.align = "off";
    const SearchResult r0 = search(prepare(sim.frames, p), p);
    std::vector<Frame> drift = sim.frames;
    const std::size_t n = drift.size();
    auto shift_of = [](std::size_t i, int& dx, int& dy) {
        dx = static_cast<int>(std::lround(0.9 * i));
        dy = static_cast<int>(std::lround(-0.35 * i));
    };
    for (std::size_t i = 0; i < n; ++i) {
        int dx, dy;
        shift_of(i, dx, dy);
        Image out(sp.width, sp.height, 0.0f);
        for (int y = 0; y < sp.height; ++y)
            for (int x = 0; x < sp.width; ++x) {
                const int sx = std::clamp(x - dx, 0, sp.width - 1), sy = std::clamp(y - dy, 0, sp.height - 1);
                out.at(x, y) = sim.frames[i].image.at(sx, sy);
            }
        drift[i].image = std::move(out);
    }
    p.align = "auto";
    const Prepared q = prepare(drift, p);
    CHECK(q.alignment.applied);
    const SearchResult r1 = search(q, p);
    int rx, ry;
    shift_of(n / 2, rx, ry);  // the reference (middle) frame's displacement
    auto find = [](const SearchResult& r, double x, double y) -> const Candidate* {
        for (const auto& c : r.candidates)
            if (c.ok && std::hypot(c.x - x, c.y - y) < 1.5) return &c;
        return nullptr;
    };
    const Candidate* c0 = find(r0, m.x_ref, m.y_ref);
    const Candidate* c1 = find(r1, m.x_ref + rx, m.y_ref + ry);
    CHECK(c0 != nullptr);
    CHECK(c1 != nullptr);
    if (c0 && c1) {
        CHECK(c1->snr > 0.85 * c0->snr);
        CHECK_NEAR(c1->vx_pxh, m.vx_pxh, 2.0);
        CHECK_NEAR(c1->vy_pxh, m.vy_pxh, 2.0);
    }
    int ok0 = 0, ok1 = 0;
    for (const auto& c : r0.candidates) ok0 += c.ok;
    for (const auto& c : r1.candidates) ok1 += c.ok;
    CHECK(ok1 <= ok0 + 1);  // no stars or hot pixels turned into movers
    CHECK(r1.threshold_snr < r0.threshold_snr + 1.0);
}

// Sky positions: sexagesimal formatting and parsing, pointing hints from FITS
// headers, MPC 80-column lines, and the ASTAP hand-over (with a stand-in solver
// on POSIX systems).
static void test_sky_positions() {
    CHECK(format_ra(169.610833) == "11 18 26.60");
    CHECK(format_dec(13.750278) == "+13 45 01.0");
    CHECK(format_dec(-0.5, 0, ':') == "-00:30:00");
    CHECK(format_ra(359.99999999) == "00 00 00.00");
    double v = 0;
    CHECK(parse_ra("11:18:26.6", v) && std::fabs(v - 169.610833) < 1e-4);
    CHECK(parse_ra("169.5", v) && std::fabs(v - 169.5) < 1e-12);
    CHECK(parse_dec("-05 07 30", v) && std::fabs(v + 5.125) < 1e-12);
    CHECK(parse_dec("-00 30 00", v) && std::fabs(v + 0.5) < 1e-12);
    CHECK(!parse_dec("95", v));

    FitsHeader h;
    h.set("OBJCTRA", std::string("11 18 30"));
    h.set("OBJCTDEC", std::string("+13 40 00"));
    h.set("FOCALLEN", 530.0);
    h.set("XPIXSZ", 7.6);
    const SolveHints hi = hints_from_header(h, 1760);
    CHECK(hi.has_center);
    CHECK_NEAR(hi.ra_deg, 169.625, 1e-9);
    CHECK_NEAR(hi.dec_deg, 13.666667, 1e-6);
    CHECK_NEAR(hi.fov_deg, 206.265 * 7.6 / 530.0 * 1760 / 3600.0, 1e-9);

    // a sequence with a known plate solution and one mover
    SimParams sp = small_sim(5);
    SimMover m;
    m.id = "A";
    m.x_ref = 80;
    m.y_ref = 100;
    m.vx_pxh = 15;
    m.vy_pxh = -10;
    m.mag = 18.8;
    const SimResult sim = simulate(sp, {m});
    SearchParams p = quick_search();
    Prepared prep = prepare(sim.frames, p);
    CHECK(prep.wcs.valid);
    const SearchResult res = search(prep, p);
    CHECK(res.has_sky);
    const auto lines = mpc_observations(prep, res, "C95");
    int ok = 0;
    for (const auto& c : res.candidates) ok += c.ok;
    CHECK(ok >= 1);
    CHECK(static_cast<int>(lines.size()) == 3 * ok);
    if (!lines.empty()) {
        const std::string& l = lines[1];
        CHECK(l.size() == 80);
        CHECK(l.substr(5, 7) == "FT00001");
        CHECK(l[14] == 'C');
        CHECK(l.substr(77, 3) == "C95");
        // the middle line is at the middle frame: RA/Dec of the track there
        const Candidate& c = res.candidates[0];
        std::vector<double> t = prep.t_h;
        std::sort(t.begin(), t.end());
        double ra, dec;
        prep.wcs.pix2sky(c.x + c.vx_pxh * t[t.size() / 2], c.y + c.vy_pxh * t[t.size() / 2], ra, dec);
        CHECK(l.substr(32, 11) == format_ra(ra, 2));
        CHECK(l.substr(44, 11) == format_dec(dec, 1));
    }

    // without a solution, positions stay in pixels; update_sky_coordinates adds them
    Prepared bare = prep;
    bare.wcs = Wcs{};
    SearchResult r2 = search(bare, p);
    CHECK(!r2.has_sky);
    CHECK(mpc_observations(bare, r2).empty());
    bare.wcs = prep.wcs;
    update_sky_coordinates(bare, r2);
    CHECK(r2.has_sky);
    if (!r2.candidates.empty() && !res.candidates.empty()) {
        CHECK_NEAR(r2.candidates[0].ra_deg, res.candidates[0].ra_deg, 1e-9);
        CHECK_NEAR(r2.candidates[0].rate_arcsec_min, res.candidates[0].rate_arcsec_min, 1e-9);
    }

    SolveOptions none;
    none.astap = (temp_dir() / "no-such-astap").string();
    CHECK(!solve_with_astap(prep.static_sky, {}, none).ok);
#ifndef _WIN32
    // stand-in for ASTAP: answers with the simulation's own solution as a .ini file
    const Wcs& w = prep.wcs;
    const fs::path fake = temp_dir() / "fake_astap.sh";
    {
        std::ofstream f(fake);
        f << "#!/bin/sh\n"
          << "while [ $# -gt 0 ]; do if [ \"$1\" = -o ]; then base=$2; fi; shift; done\n"
          << "printf 'PLTSOLVD=T\\nCRVAL1=%s\\nCRVAL2=%s\\nCRPIX1=%s\\nCRPIX2=%s\\nCD1_1=%s\\nCD1_2=%s\\nCD2_1=%s\\nCD2_2=%s\\n' "
          << std::setprecision(17) << w.crval1 << ' ' << w.crval2 << ' ' << w.crpix1 << ' ' << w.crpix2 << ' ' << w.cd11
          << ' ' << w.cd12 << ' ' << w.cd21 << ' ' << w.cd22 << " > \"$base.ini\"\n";
    }
    fs::permissions(fake, fs::perms::owner_all);
    SolveOptions o;
    o.astap = fake.string();
    Prepared q = prep;
    q.wcs = Wcs{};
    const SolveResult sr = solve_prepared(q, hi, o);
    if (!sr.ok) std::printf("    %s\n", sr.message.c_str());
    CHECK(sr.ok && q.wcs.valid);
    CHECK_NEAR(q.wcs.crval1, w.crval1, 1e-9);
    CHECK_NEAR(sr.scale_arcsec, w.pixel_scale_arcsec(), 1e-6);
#endif
}

// Photometric calibration: field stars of a simulation as the "catalogue", movers
// of known magnitude (one fast enough to trail noticeably) come out right.
static void test_photometry() {
    // band conversions (Gaia EDR3 relations): a solar-coloured star
    CatalogStar sun;
    sun.g = 15.0;
    sun.bp_rp = 0.82;
    sun.has_colour = true;
    CHECK_NEAR(catalog_mag(sun, 'G'), 15.0, 1e-12);
    CHECK_NEAR(catalog_mag(sun, 'V'), 15.0 + 0.1525, 0.002);
    CHECK(catalog_mag(sun, 'R') < 15.0);
    CatalogStar nocol = sun;
    nocol.has_colour = false;
    CHECK(std::isnan(catalog_mag(nocol, 'V')));
    CHECK_NEAR(trail_flux_factor(0.0, 3.0), 1.0, 1e-12);
    CHECK(trail_flux_factor(3.0, 3.0) < 0.95 && trail_flux_factor(3.0, 3.0) > 0.85);

    SimParams sp = small_sim(77);
    sp.width = 320;
    sp.height = 256;
    sp.nstars = 900;
    sp.exptime_s = 180;
    sp.cadence_s = 190;
    sp.nframes = 30;
    std::vector<SimMover> mv(2);
    mv[0] = {"slow", 90, 80, 16, 9, 18.6};
    mv[1] = {"fast", 220, 170, -38, 14, 18.9};  // trails over 2 px in 180 s
    const SimResult sim = simulate(sp, mv);
    SearchParams p = quick_search();
    p.vmax_pxh = 45;
    const Prepared prep = prepare(sim.frames, p);
    SearchResult res = search(prep, p);
    // the simulation's stars as a Gaia-like catalogue (G = simulation magnitude)
    std::vector<CatalogStar> cat;
    for (const auto& st : sim.stars) {
        CatalogStar c;
        prep.wcs.pix2sky(st.x, st.y, c.ra_deg, c.dec_deg);
        c.g = st.mag;
        c.bp_rp = 0.8;
        c.has_colour = true;
        cat.push_back(c);
    }
    const PhotometryResult cal = calibrate_photometry(prep, cat, 'G');
    if (!cal.ok) std::printf("    %s\n", cal.message.c_str());
    CHECK(cal.ok);
    CHECK(cal.stars >= 20);
    CHECK(cal.rms < 0.08);
    // the true zero point of a frame: mag = zero_point - 2.5 log10(flux_e / exptime) with gain 1
    CHECK_NEAR(cal.zero_point, sp.zero_point + 2.5 * std::log10(sp.exptime_s), 0.06);
    apply_photometry(prep, cal, res);
    CHECK(res.mag_band == 'G');
    for (const auto& m : mv) {
        // the photometry is tested here, not the candidate tests: take the strongest
        // candidate on the mover's track whether or not it passed (the noise of the
        // simulation differs slightly between compilers' maths libraries)
        const Candidate* hit = nullptr;
        for (const auto& c : res.candidates)
            if (std::hypot(c.x - m.x_ref, c.y - m.y_ref) < 2.0 && std::hypot(c.vx_pxh - m.vx_pxh, c.vy_pxh - m.vy_pxh) < 3.0 &&
                (!hit || c.snr > hit->snr))
                hit = &c;
        if (!hit) {
            std::printf("    mover %s not found; nearby candidates:\n", m.id.c_str());
            for (const auto& c : res.candidates)
                if (std::hypot(c.x - m.x_ref, c.y - m.y_ref) < 10.0)
                    std::printf("      ok %d (%.1f, %.1f) v (%.1f, %.1f) S/N %.1f %s\n", c.ok, c.x, c.y, c.vx_pxh, c.vy_pxh, c.snr,
                                c.flags.c_str());
        }
        CHECK(hit != nullptr);
        if (hit) {
            CHECK_NEAR(hit->mag, m.mag, std::max(0.12, 3.0 * hit->mag_err));
            CHECK(hit->mag_err > 0.0 && hit->mag_err < 0.2);
        }
    }
    // MPC lines carry the magnitude and band in columns 66-71
    const auto lines = mpc_observations(prep, res, "C95", 3, true);
    CHECK(!lines.empty());
    if (!lines.empty()) {
        CHECK(lines[0][70] == 'G');
        CHECK(std::isdigit(static_cast<unsigned char>(lines[0][66])));
    }
}

// 16-bit fixture files and cropping: quantisation well below the noise, NaN kept,
// WCS and registration offsets consistent after the crop.
static void test_crop_and_int16() {
    detail_rng_seed(5);
    Frame f;
    f.image = Image(64, 48);
    for (auto& v : f.image.data) v = 1000.0f + 10.0f * static_cast<float>(detail_rng_normal());
    f.image.at(40, 30) = 60000.0f;  // saturated star: may be clipped, nothing else
    f.image.at(0, 0) = NAN;
    f.mjd_mid = 60000.5;
    f.exptime_s = 60;
    f.wcs = Wcs::make_tan(10.0, 20.0, 32.5, 24.5, 1.5, 30.0);
    f.has_reg_offset = true;
    f.reg_dx = 3.25;
    f.reg_dy = -1.5;
    const Frame c = crop_frame(f, 10, 5, 40, 30);
    CHECK(c.image.width == 40 && c.image.height == 30);
    CHECK(c.image.at(0, 0) == f.image.at(10, 5));
    double ra1 = 0, de1 = 0, ra2 = 0, de2 = 0;
    f.wcs.pix2sky(25.0, 17.0, ra1, de1);
    c.wcs.pix2sky(15.0, 12.0, ra2, de2);
    CHECK_NEAR(ra1, ra2, 1e-12);
    CHECK_NEAR(de1, de2, 1e-12);
    CHECK_NEAR(c.reg_dx, 13.25, 1e-12);  // detector x = 15 + 13.25 = (15 + 10) + 3.25
    CHECK_NEAR(c.reg_dy, 3.5, 1e-12);
    bool threw = false;
    try {
        crop_frame(f, 30, 0, 40, 10);
    } catch (const std::exception&) {
        threw = true;
    }
    CHECK(threw);

    const std::string path = (temp_dir() / "int16.fits").string();
    write_frame(path, f, {}, 16);
    CHECK(fs::file_size(path) < 64 * 48 * 4);
    FitsImage back = read_fits(path);
    CHECK(back.header.get_int("BITPIX", 0) == 16);
    CHECK(std::isnan(back.image.at(0, 0)));
    double worst = 0.0;
    for (int y = 0; y < 48; ++y)
        for (int x = 0; x < 64; ++x)
            if ((x || y) && !(x == 40 && y == 30))
                worst = std::max(worst, std::fabs(static_cast<double>(back.image.at(x, y)) - f.image.at(x, y)));
    CHECK(worst <= 10.0 / 4.0 * 0.5 + 1e-3);  // half a quantisation step, step <= noise / 4
    CHECK(back.image.at(40, 30) > 1100.0f);   // the star stays far above the sky
    const Frame g = frame_from_fits(std::move(back), path);
    CHECK(g.has_reg_offset);
    CHECK_NEAR(g.reg_dx, 3.25, 1e-9);
}

// Case files parse, resolve paths relative to the case, and skip missing data.
static void test_regression_case() {
    const fs::path dir = temp_dir() / "regcase";
    fs::create_directories(dir / "frames");
    {
        std::ofstream o(dir / "x.case");
        o << "# comment\nname = demo\ndata = frames   # trailing comment\nbaseline = b/x.csv\n"
             "[search]\nvmin = 4\nvmax = 20.5\n[expect]\nmax_unexpected = 1\nthreshold_min = 5\n"
             "[object A one]\nmjd = 60000.25\nx = 10\ny = 11.5\nvx = 1\nvy = -2\nra = 10.5\ndec = -3\n";
    }
    const RegressionCase c = load_regression_case((dir / "x.case").string());
    CHECK(c.name == "demo");
    CHECK(fs::equivalent(c.data_dir, dir / "frames"));
    CHECK(c.baseline_path == (dir / "b" / "x.csv").lexically_normal().string());
    CHECK_NEAR(c.params.vmax_pxh, 20.5, 1e-12);
    CHECK(c.max_unexpected == 1);
    CHECK(c.objects.size() == 1 && c.objects[0].name == "A one" && c.objects[0].has_sky);
    CHECK_NEAR(c.objects[0].y, 11.5, 1e-12);
    // empty frames directory and an unset variable are skips, not failures
    CHECK(run_regression_case(c).status == RegressionStatus::Skipped);
    {
        std::ofstream o(dir / "y.case");
        o << "data = ${FAINTTRAIL_SURELY_UNSET_VARIABLE}/x\n";
    }
    CHECK(run_regression_case(load_regression_case((dir / "y.case").string())).status == RegressionStatus::Skipped);
    {
        std::ofstream o(dir / "z.case");
        o << "data = frames\n[search]\nvmax = fast\n";
    }
    bool threw = false;
    try {
        load_regression_case((dir / "z.case").string());
    } catch (const std::exception&) {
        threw = true;
    }
    CHECK(threw);
}

// Memory layout: residuals are dropped by default and recomputed exactly where
// needed; handing the frames over gives the same result as copying them.
static void test_memory_layout() {
    const SimResult sim = simulate(small_sim(21));
    SearchParams p = quick_search();
    const Prepared a = prepare(sim.frames, p);
    CHECK(a.residual.empty());
    CHECK(static_cast<int>(a.sub.size()) == a.nframes);
    SearchParams pk = p;
    pk.keep_residuals = true;
    const Prepared b = prepare(sim.frames, pk);
    CHECK(static_cast<int>(b.residual.size()) == b.nframes);
    bool same = true;
    for (int i = 0; i < a.nframes && same; ++i) same = a.filtered[static_cast<std::size_t>(i)].data == b.filtered[static_cast<std::size_t>(i)].data;
    CHECK(same);
    // residual cutouts: recomputed (a) equal kept (b), bit for bit, also at the border
    for (const auto& pos : {std::array<double, 4>{60.3, 70.7, 12.0, -7.5}, std::array<double, 4>{1.2, 2.5, 20.0, 3.0},
                            std::array<double, 4>{190.4, 150.0, -15.0, 9.0}}) {
        const Image ca = stack_cutout(a, pos[0], pos[1], pos[2], pos[3], 6, false, 0);
        const Image cb = stack_cutout(b, pos[0], pos[1], pos[2], pos[3], 6, false, 0);
        bool eq = ca.data.size() == cb.data.size();
        for (std::size_t k = 0; eq && k < ca.data.size(); ++k)
            eq = (std::isnan(ca.data[k]) && std::isnan(cb.data[k])) || ca.data[k] == cb.data[k];
        CHECK(eq);
    }
    bool threw = false;
    try {
        stack_full(a, 10.0, 0.0);
    } catch (const std::runtime_error&) {
        threw = true;
    }
    CHECK(threw);
    // moving the frames in gives the same search result
    std::vector<Frame> copy = sim.frames;
    const Prepared m = prepare(std::move(copy), p);
    CHECK(copy.empty() || copy.front().image.data.empty());
    const SearchResult r1 = search(a, p), r2 = search(m, p);
    CHECK(r1.candidates.size() == r2.candidates.size());
    CHECK(r1.threshold_snr == r2.threshold_snr);
    for (std::size_t k = 0; k < std::min(r1.candidates.size(), r2.candidates.size()); ++k) {
        CHECK(r1.candidates[k].x == r2.candidates[k].x);
        CHECK(r1.candidates[k].flux == r2.candidates[k].flux);
    }
    CHECK(!r1.timings.empty() && r1.timings.front().rss_mb >= 0.0);
}

// Tiled search: the same result as the untiled search, bit for bit, for tile
// sizes that do or do not divide the image, from memory and from FITS files,
// with blank registration borders and registration offsets in the data.
static void test_tiled_equals_untiled() {
    SimParams sp = small_sim(31);
    sp.nmovers = 4;
    SimResult sim = simulate(sp);
    // registration-like blank borders that move with the frame, and offsets
    for (std::size_t i = 0; i < sim.frames.size(); ++i) {
        Frame& f = sim.frames[i];
        const int cut = static_cast<int>(i % 5);
        for (int y = 0; y < f.image.height; ++y)
            for (int x = 0; x < cut; ++x) f.image.at(x, y) = NAN;
        f.has_reg_offset = true;
        f.reg_dx = -0.5 * static_cast<double>(i % 5);
        f.reg_dy = 0.25 * static_cast<double>(i % 3);
    }
    SearchParams p = quick_search();
    const SearchResult ref = search(prepare(sim.frames, p), p);
    CHECK(!ref.candidates.empty());
    auto same = [&](const SearchResult& r, const char* what) {
        bool ok = r.threshold_snr == ref.threshold_snr && r.candidates.size() == ref.candidates.size() &&
                  r.max_snr.data == ref.max_snr.data && r.fwhm_px == ref.fwhm_px;
        for (std::size_t k = 0; ok && k < r.candidates.size(); ++k) {
            const Candidate &a = r.candidates[k], &b = ref.candidates[k];
            ok = a.x == b.x && a.y == b.y && a.vx_pxh == b.vx_pxh && a.vy_pxh == b.vy_pxh && a.snr == b.snr &&
                 a.snr_coarse == b.snr_coarse && a.snr_half1 == b.snr_half1 && a.snr_half2 == b.snr_half2 &&
                 a.flux == b.flux && a.flags == b.flags && a.ok == b.ok && a.max_frame_fraction == b.max_frame_fraction;
        }
        if (!ok) std::printf("    %s: threshold %.6f vs %.6f, %zu vs %zu candidates\n", what, r.threshold_snr,
                             ref.threshold_snr, r.candidates.size(), ref.candidates.size());
        CHECK(ok);
    };
    const auto mem = memory_frame_source(sim.frames);
    for (int tile : {48, 100, 1000}) {
        TileOptions to;
        to.tile = tile;
        same(search_tiled(*mem, p, to), tile == 48 ? "tile 48" : tile == 100 ? "tile 100" : "tile 1000");
    }
    // from FITS files on disk
    const fs::path dir = temp_dir() / "tiled";
    fs::create_directories(dir);
    std::vector<std::string> files;
    for (std::size_t i = 0; i < sim.frames.size(); ++i) {
        files.push_back((dir / ("f" + std::to_string(100 + i) + ".fits")).string());
        write_frame(files.back(), sim.frames[i]);
    }
    const SearchResult from_disk_untiled = search(prepare(load_frames(files), p), p);
    const auto disk = fits_frame_source(files);
    TileOptions to;
    to.tile = 64;
    Prepared meta;
    const SearchResult tiled_disk = search_tiled(*disk, p, to, &meta);
    CHECK(tiled_disk.candidates.size() == from_disk_untiled.candidates.size());
    CHECK(tiled_disk.max_snr.data == from_disk_untiled.max_snr.data);
    CHECK(meta.nframes == static_cast<int>(files.size()) && meta.fwhm_px == ref.fwhm_px);
    // a planned tiling under a small budget uses several tiles
    const TilePlan plan = plan_tiles(*mem, p, TileOptions{0, 40.0});
    CHECK(plan.tiles_x * plan.tiles_y > 1);
}

static void test_c_api() { CHECK(run_c_api_test() == 0); }

int main() {
    struct Test {
        const char* name;
        std::function<void()> fn;
    };
    const std::vector<Test> tests = {
        {"fits round trip", test_fits_roundtrip},
        {"time conversion", test_time_conversion},
        {"wcs round trip", test_wcs_roundtrip},
        {"noise normalisation", test_noise_normalisation},
        {"fwhm estimate", test_fwhm_estimate},
        {"known-motion stack", test_known_motion_stack},
        {"blind search recovers movers", test_blind_search_recovers},
        {"no detections in pure noise", test_null_has_no_detections},
        {"injection-recovery", test_completeness},
        {"csv round trip", test_csv_roundtrip},
        {"cancellation", test_cancellation},
        {"plate-solution transfer", test_wcs_transfer},
        {"backend contract", test_backend_contract},
        {"shift tables", test_shift_tables},
        {"CUDA algorithm (emulated)", test_cuda_algorithm},
        {"registration check", test_registration_check},
        {"aligned search of a drifting sequence", test_aligned_search},
        {"sky positions and plate solving", test_sky_positions},
        {"photometric calibration", test_photometry},
        {"crop and 16-bit FITS", test_crop_and_int16},
        {"regression case files", test_regression_case},
        {"memory layout", test_memory_layout},
        {"tiled search equals untiled", test_tiled_equals_untiled},
        {"C API", test_c_api},
    };
    int failed = 0;
    for (const auto& t : tests) {
        const int before = g_failed_checks;
        std::printf("[ RUN  ] %s\n", t.name);
        std::fflush(stdout);
        try {
            t.fn();
        } catch (const std::exception& e) {
            std::printf("    exception: %s\n", e.what());
            ++g_failed_checks;
        }
        const bool ok = g_failed_checks == before;
        std::printf("[ %s ] %s\n", ok ? " OK " : "FAIL", t.name);
        failed += ok ? 0 : 1;
    }
    std::printf("\n%zu tests, %d failed\n", tests.size(), failed);
    return failed == 0 ? 0 : 1;
}
