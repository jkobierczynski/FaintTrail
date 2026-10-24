// SPDX-License-Identifier: GPL-3.0-or-later
// Copyright 2026 Jurgen Kobierczynski
//
// Unit and end-to-end tests. Plain C++, no framework needed: each test is a
// function; a failed CHECK prints the location and marks the test as failed.
#include <cmath>
#include <cstdio>
#include <filesystem>
#include <functional>
#include <string>
#include <vector>

#include "fainttrail/fainttrail.hpp"

extern "C" int run_c_api_test(void);

using namespace fainttrail;
namespace fs = std::filesystem;

static int g_failed_checks = 0;

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
    CHECK(res.threshold_method == "null-fit" || res.threshold_method == "null-quantile");
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
