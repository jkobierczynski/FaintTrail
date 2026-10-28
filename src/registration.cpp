// SPDX-License-Identifier: GPL-3.0-or-later
// Copyright 2026 Jurgen Kobierczynski
//
// Registration check and frame alignment (registration.hpp).
//
// Stars: bright, extended local maxima, so hot pixels and cosmic rays (single
// sharp pixels; hot pixels also do not move and would vote for a zero offset)
// are not used. Matching: every pair of a bright reference star and a frame star
// votes for its offset; the most popular offset (and, for meridian flips, the
// same with the reference rotated by 180 degrees) gives the starting transform,
// which is refined by fitting a similarity transform to nearest-neighbour matches
// with outlier clipping. Resampling: bicubic convolution (Keys, a = -0.5); output
// pixels whose 4 x 4 source window leaves the frame or contains NaN become NaN.
#include "fainttrail/registration.hpp"

#include <algorithm>
#include <cmath>
#include <iomanip>
#include <limits>
#include <map>
#include <sstream>
#include <stdexcept>
#include <utility>
#include <vector>

#include "common.hpp"

namespace fainttrail {

namespace {

struct Star {
    double x, y, peak;
};

std::vector<Star> find_stars(const Image& im, int keep) {
    std::vector<float> s = detail::sample_pixels(im, 200000);
    if (s.size() < 100) return {};
    const double bg = detail::median_inplace(s);
    const double sig = detail::robust_sigma(s);
    if (!(sig > 0)) return {};
    const int w = im.width, h = im.height, b = 4;
    const float thr = static_cast<float>(bg + 5.0 * sig), nthr = static_cast<float>(bg + 2.5 * sig);
    std::vector<Star> out;
    for (int y = b; y < h - b; ++y)
        for (int x = b; x < w - b; ++x) {
            const float v = im.at(x, y);
            if (!(v > thr)) continue;
            bool is_max = true;
            int lit = 0;
            for (int dy = -2; dy <= 2 && is_max; ++dy)
                for (int dx = -2; dx <= 2; ++dx) {
                    if (!dx && !dy) continue;
                    const float u = im.at(x + dx, y + dy);
                    if (!std::isfinite(u) || u > v || (u == v && (dy < 0 || (dy == 0 && dx < 0)))) {
                        is_max = false;
                        break;
                    }
                    if (std::abs(dx) <= 1 && std::abs(dy) <= 1 && u > nthr) ++lit;
                }
            if (!is_max || lit < 5) continue;
            // centroid in a 7 x 7 window, recentred once
            double cx = x, cy = y;
            for (int it = 0; it < 2; ++it) {
                const int ix = static_cast<int>(std::lround(cx)), iy = static_cast<int>(std::lround(cy));
                if (ix < 3 || iy < 3 || ix >= w - 3 || iy >= h - 3) break;
                double sw = 0, sx = 0, sy = 0;
                for (int dy = -3; dy <= 3; ++dy)
                    for (int dx = -3; dx <= 3; ++dx) {
                        const float u = im.at(ix + dx, iy + dy);
                        if (!std::isfinite(u)) continue;
                        const double wgt = std::max(0.0, u - bg - 2.0 * sig);
                        sw += wgt;
                        sx += wgt * (ix + dx);
                        sy += wgt * (iy + dy);
                    }
                if (sw <= 0) break;
                cx = sx / sw;
                cy = sy / sw;
            }
            out.push_back({cx, cy, v - bg});
        }
    std::sort(out.begin(), out.end(), [](const Star& a, const Star& c) { return a.peak > c.peak; });
    if (static_cast<int>(out.size()) > keep) out.resize(static_cast<std::size_t>(keep));
    return out;
}

// frame pixel = (a, -b; b, a) * ref pixel + (tx, ty)
struct Similarity {
    double a = 1, b = 0, tx = 0, ty = 0;
    void apply(double x, double y, double& u, double& v) const {
        u = a * x - b * y + tx;
        v = b * x + a * y + ty;
    }
    double scale() const { return std::hypot(a, b); }
    double rot_deg() const { return std::atan2(b, a) * 180.0 / 3.14159265358979323846; }
};

// Coarse offset by voting; flip: reference rotated by 180 degrees about (cx, cy).
int vote(const std::vector<Star>& ref, const std::vector<Star>& fr, bool flip, double cx, double cy, double max_shift,
         double& ox, double& oy) {
    const std::size_t nr = std::min<std::size_t>(ref.size(), 30), nf = std::min<std::size_t>(fr.size(), 80);
    std::map<std::pair<int, int>, int> votes;  // 2-px bins
    for (std::size_t i = 0; i < nr; ++i) {
        const double px = flip ? 2 * cx - ref[i].x : ref[i].x, py = flip ? 2 * cy - ref[i].y : ref[i].y;
        for (std::size_t k = 0; k < nf; ++k) {
            const double dx = fr[k].x - px, dy = fr[k].y - py;
            if (std::fabs(dx) > max_shift || std::fabs(dy) > max_shift) continue;
            ++votes[{static_cast<int>(std::floor(dx / 2.0)), static_cast<int>(std::floor(dy / 2.0))}];
        }
    }
    int best = 0;
    std::pair<int, int> key{0, 0};
    for (const auto& kv : votes) {  // best 2 x 2 block of bins: an offset on a bin edge is not split
        int sum = 0;
        for (int u = 0; u <= 1; ++u)
            for (int v = 0; v <= 1; ++v) {
                const auto it = votes.find({kv.first.first + u, kv.first.second + v});
                if (it != votes.end()) sum += it->second;
            }
        if (sum > best) {
            best = sum;
            key = kv.first;
        }
    }
    ox = 2.0 * key.first + 2.0;
    oy = 2.0 * key.second + 2.0;
    return best;
}

// Fit the transform mapping reference stars onto frame stars. Returns false when
// fewer than `min_match` stars match.
bool fit_transform(const std::vector<Star>& ref, const std::vector<Star>& fr, double cx, double cy, double max_shift,
                   bool allow_flip, Similarity& T, double& rms, int& used) {
    const std::size_t nr = std::min<std::size_t>(ref.size(), 30);
    if (nr < 6 || fr.size() < 6) return false;
    double ox = 0, oy = 0, fx = 0, fy = 0;
    const int v0 = vote(ref, fr, false, cx, cy, max_shift, ox, oy);
    const int v1 = allow_flip ? vote(ref, fr, true, cx, cy, max_shift, fx, fy) : 0;
    const bool flip = v1 > v0;
    const int votes = std::max(v0, v1);
    if (votes < std::max<int>(5, static_cast<int>(nr) / 3)) return false;
    T = flip ? Similarity{-1, 0, 2 * cx + fx, 2 * cy + fy} : Similarity{1, 0, ox, oy};
    const double radii[] = {4.0, 2.5, 1.5, 1.5};
    std::vector<std::pair<const Star*, const Star*>> pairs;
    for (const double r0 : radii) {
        pairs.clear();
        for (const Star& s : ref) {
            double u, v;
            T.apply(s.x, s.y, u, v);
            const Star* best = nullptr;
            double bd = r0;
            for (const Star& f : fr) {
                const double d = std::hypot(f.x - u, f.y - v);
                if (d < bd) {
                    bd = d;
                    best = &f;
                }
            }
            if (best) pairs.emplace_back(&s, best);
        }
        if (pairs.size() < 5) return false;
        // least squares on centred coordinates, twice with 3-sigma clipping
        for (int pass = 0; pass < 2; ++pass) {
            double mpx = 0, mpy = 0, mqx = 0, mqy = 0;
            for (const auto& pr : pairs) {
                mpx += pr.first->x;
                mpy += pr.first->y;
                mqx += pr.second->x;
                mqy += pr.second->y;
            }
            const double m = static_cast<double>(pairs.size());
            mpx /= m; mpy /= m; mqx /= m; mqy /= m;
            double sab = 0, sbb = 0, spp = 0;
            for (const auto& pr : pairs) {
                const double px = pr.first->x - mpx, py = pr.first->y - mpy;
                const double qx = pr.second->x - mqx, qy = pr.second->y - mqy;
                sab += px * qx + py * qy;
                sbb += px * qy - py * qx;
                spp += px * px + py * py;
            }
            if (spp <= 0) return false;
            if (pairs.size() >= 12) {
                T.a = sab / spp;
                T.b = sbb / spp;
            } else {
                // few stars: rotation and scale are poorly determined; fit the shift only
                // (keeping the 180-degree flip if that is what matched)
                T.a = T.a < 0 ? -1.0 : 1.0;
                T.b = 0.0;
            }
            T.tx = mqx - (T.a * mpx - T.b * mpy);
            T.ty = mqy - (T.b * mpx + T.a * mpy);
            std::vector<double> res;
            double ss = 0;
            for (const auto& pr : pairs) {
                double u, v;
                T.apply(pr.first->x, pr.first->y, u, v);
                const double d = std::hypot(pr.second->x - u, pr.second->y - v);
                res.push_back(d);
                ss += d * d;
            }
            rms = std::sqrt(ss / m);
            const double cut = std::max(0.5, 3.0 * rms);
            std::vector<std::pair<const Star*, const Star*>> keep;
            for (std::size_t k = 0; k < pairs.size(); ++k)
                if (res[k] <= cut) keep.push_back(pairs[k]);
            if (keep.size() < 5) return false;
            pairs.swap(keep);
        }
    }
    used = static_cast<int>(pairs.size());
    return used >= 5;
}

// Keys cubic convolution kernel, a = -0.5.
inline void cubic_weights(double t, double w[4]) {
    const double t2 = t * t, t3 = t2 * t;
    w[0] = -0.5 * t3 + t2 - 0.5 * t;
    w[1] = 1.5 * t3 - 2.5 * t2 + 1.0;
    w[2] = -1.5 * t3 + 2.0 * t2 + 0.5 * t;
    w[3] = 0.5 * t3 - 0.5 * t2;
}

Image resample(const Image& src, int ow, int oh, const Similarity& T, int nt) {
    Image out(ow, oh, std::numeric_limits<float>::quiet_NaN());
    const int w = src.width, h = src.height;
    (void)nt;
#pragma omp parallel for num_threads(nt) schedule(static)
    for (int y = 0; y < oh; ++y)
        for (int x = 0; x < ow; ++x) {
            double u, v;
            T.apply(x, y, u, v);
            const int ix = static_cast<int>(std::floor(u)), iy = static_cast<int>(std::floor(v));
            if (ix < 1 || iy < 1 || ix + 2 >= w || iy + 2 >= h) continue;
            double wx[4], wy[4];
            cubic_weights(u - ix, wx);
            cubic_weights(v - iy, wy);
            double acc = 0;
            bool bad = false;
            for (int j = 0; j < 4 && !bad; ++j) {
                const float* row = &src.data[static_cast<std::size_t>(iy - 1 + j) * w + (ix - 1)];
                double r = 0;
                for (int i = 0; i < 4; ++i) {
                    if (!std::isfinite(row[i])) {
                        bad = true;
                        break;
                    }
                    r += wx[i] * row[i];
                }
                acc += wy[j] * r;
            }
            if (!bad) out.data[static_cast<std::size_t>(y) * ow + x] = static_cast<float>(acc);
        }
    return out;
}

}  // namespace

DriftCheck check_registration(const std::vector<const Image*>& frames, const std::vector<double>& t_h) {
    DriftCheck r;
    const std::size_t n = frames.size();
    if (n < 3 || t_h.size() != n) return r;
    std::vector<std::size_t> order(n);
    for (std::size_t i = 0; i < n; ++i) order[i] = i;
    std::sort(order.begin(), order.end(), [&](std::size_t a, std::size_t b) { return t_h[a] < t_h[b]; });
    const std::size_t first = order.front(), last = order.back(), mid = order[n / 2];
    const Image& a = *frames[first];
    const double cx = 0.5 * (a.width - 1), cy = 0.5 * (a.height - 1);
    const double max_shift = std::min(a.width, a.height) / 3.0;
    const auto sa = find_stars(a, 60);
    const auto sm = find_stars(*frames[mid], 120);
    const auto sl = find_stars(*frames[last], 120);
    r.stars = static_cast<int>(sa.size());
    Similarity tl, tm;
    double rms = 0;
    int nl = 0, nm = 0;
    if (!fit_transform(sa, sl, cx, cy, max_shift, true, tl, rms, nl)) return r;
    r.measured = true;
    r.matched = nl;
    // displacement of the field centre; with a rotation, also the corners count
    auto worst = [&](const Similarity& T, double& dx, double& dy) {
        T.apply(cx, cy, dx, dy);
        dx -= cx;
        dy -= cy;
        double m = 0;
        for (const double x : {0.0, 2 * cx})
            for (const double y : {0.0, 2 * cy}) {
                double u, v;
                T.apply(x, y, u, v);
                m = std::max(m, std::hypot(u - x, v - y));
            }
        return m;
    };
    const double ml = worst(tl, r.dx, r.dy);
    double mm = 0;
    if (fit_transform(sa, sm, cx, cy, max_shift, true, tm, rms, nm)) mm = worst(tm, r.mid_dx, r.mid_dy);
    const double dt = t_h[last] - t_h[first];
    if (dt > 0) r.rate_pxh = std::hypot(r.dx, r.dy) / dt;
    r.registered = ml <= 1.5 && mm <= 1.5;
    if (!r.registered) {
        std::ostringstream s;
        s.setf(std::ios::fixed);
        s.precision(1);
        s << "the frames are not registered: the stars move by (" << r.dx << ", " << r.dy
          << ") px from the first to the last frame (" << r.rate_pxh << " px/h, " << nl << " stars matched";
        if (std::fabs(tl.rot_deg()) > 0.05) s << ", rotation " << tl.rot_deg() << " deg";
        s << "). Without registration every star is detected as a moving object.";
        r.message = s.str();
    }
    return r;
}

std::string AlignmentReport::summary() const {
    std::ostringstream s;
    s.setf(std::ios::fixed);
    s.precision(2);
    if (!applied) return "not needed";
    double maxd = 0, maxrot = 0, sr = 0;
    int k = 0;
    for (const auto& f : frames)
        if (f.aligned) {
            maxd = std::max(maxd, std::hypot(f.dx, f.dy));
            maxrot = std::max(maxrot, std::fabs(f.rot_deg));
            sr += f.rms_px * f.rms_px;
            ++k;
        }
    s << k << " frames aligned on frame " << reference + 1 << ", largest shift " << std::setprecision(1) << maxd << " px";
    if (maxrot > 0.05) s << ", rotation up to " << maxrot << " deg";
    s << std::setprecision(2) << ", star residual " << (k ? std::sqrt(sr / k) : 0.0) << " px rms";
    if (dropped) s << ", " << dropped << " frame(s) dropped (too few stars matched)";
    return s.str();
}

AlignmentReport register_frames(std::vector<Frame>& frames, int threads) {
    AlignmentReport rep;
    const std::size_t n = frames.size();
    if (n == 0) return rep;
    const int nt = detail::thread_count(threads);
    std::vector<std::size_t> order(n);
    for (std::size_t i = 0; i < n; ++i) order[i] = i;
    std::sort(order.begin(), order.end(), [&](std::size_t a, std::size_t b) { return frames[a].mjd_mid < frames[b].mjd_mid; });
    const std::size_t ref = order[n / 2];
    rep.reference = static_cast<int>(ref);
    const Image& R = frames[ref].image;
    const int ow = R.width, oh = R.height;
    const double cx = 0.5 * (ow - 1), cy = 0.5 * (oh - 1);
    const auto rstars = find_stars(R, 250);
    if (rstars.size() < 6)
        throw std::runtime_error("cannot register the frames: only " + std::to_string(rstars.size()) +
                                 " usable stars in the reference frame " + frames[ref].source);

    std::vector<Similarity> T(n);
    rep.frames.resize(n);
    std::vector<int> ok(n, 0);
#pragma omp parallel for num_threads(nt) schedule(dynamic)
    for (int i = 0; i < static_cast<int>(n); ++i) {
        FrameAlignment& fa = rep.frames[static_cast<std::size_t>(i)];
        fa.source = frames[static_cast<std::size_t>(i)].source;
        if (static_cast<std::size_t>(i) == ref) {
            fa.aligned = true;
            fa.stars = static_cast<int>(rstars.size());
            ok[static_cast<std::size_t>(i)] = 1;
            continue;
        }
        const Image& im = frames[static_cast<std::size_t>(i)].image;
        const auto fs = find_stars(im, 350);
        const double max_shift = std::max(im.width, im.height) / 2.0;
        double rms = 0;
        int used = 0;
        Similarity t;
        if (fit_transform(rstars, fs, cx, cy, max_shift, true, t, rms, used) && std::fabs(t.scale() - 1.0) < 0.02) {
            T[static_cast<std::size_t>(i)] = t;
            double u, v;
            t.apply(cx, cy, u, v);
            fa.aligned = true;
            fa.dx = u - cx;
            fa.dy = v - cy;
            fa.rot_deg = t.rot_deg();
            fa.scale = t.scale();
            fa.rms_px = rms;
            fa.stars = used;
            ok[static_cast<std::size_t>(i)] = 1;
        }
    }

    // carry a plate solution over to the reference grid
    Wcs wcs;
    for (std::size_t k = 0; k <= n && !wcs.valid; ++k) {
        const std::size_t i = k == 0 ? ref : k - 1;
        if (i >= n || !ok[i] || !frames[i].wcs.valid) continue;
        const Wcs& w = frames[i].wcs;
        const Similarity& t = T[i];
        // sky = CD (A p0 + t + 1 - crpix)  =>  CD' = CD A, crpix' = 1 + A^-1 (crpix - 1 - t)
        wcs = w;
        wcs.cd11 = w.cd11 * t.a + w.cd12 * t.b;
        wcs.cd12 = -w.cd11 * t.b + w.cd12 * t.a;
        wcs.cd21 = w.cd21 * t.a + w.cd22 * t.b;
        wcs.cd22 = -w.cd21 * t.b + w.cd22 * t.a;
        const double qx = w.crpix1 - 1 - t.tx, qy = w.crpix2 - 1 - t.ty, s2 = t.a * t.a + t.b * t.b;
        wcs.crpix1 = 1 + (t.a * qx + t.b * qy) / s2;
        wcs.crpix2 = 1 + (-t.b * qx + t.a * qy) / s2;
    }

    std::vector<Frame> out;
    out.reserve(n);
    std::vector<FrameAlignment> fas;
    for (std::size_t i = 0; i < n; ++i) {
        if (!ok[i]) {
            ++rep.dropped;
            fas.push_back(rep.frames[i]);
            continue;
        }
        Frame f = std::move(frames[i]);
        if (i != ref) f.image = resample(f.image, ow, oh, T[i], nt);
        const FrameAlignment& fa = rep.frames[i];
        // detector = registered + offset (the frame's own offset, if it had one, adds up)
        f.reg_dx = (f.has_reg_offset ? f.reg_dx : 0.0) + fa.dx;
        f.reg_dy = (f.has_reg_offset ? f.reg_dy : 0.0) + fa.dy;
        f.has_reg_offset = true;
        if (wcs.valid) f.wcs = wcs;
        out.push_back(std::move(f));
        fas.push_back(fa);
    }
    rep.frames = std::move(fas);  // input order; dropped frames have aligned = false
    frames = std::move(out);
    rep.applied = true;
    if (frames.size() < 5)
        throw std::runtime_error("cannot register the frames: only " + std::to_string(frames.size()) +
                                 " frames could be aligned (at least 5 needed)");
    return rep;
}

}  // namespace fainttrail
