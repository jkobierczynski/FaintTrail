// SPDX-License-Identifier: GPL-3.0-or-later
// Copyright 2026 Jurgen Kobierczynski
#include "fainttrail/astrometry.hpp"

#include <algorithm>
#include <cmath>
#include <cstdint>
#include <stdexcept>
#include <unordered_map>

#include "common.hpp"
#include "fainttrail/fits.hpp"

namespace fainttrail {

namespace {
constexpr double kPi = 3.14159265358979323846;

// 3x3 box-free Gaussian smoothing with sigma ~1 px (separable 5-tap kernel).
Image smooth(const Image& in) {
    static const float k[5] = {0.0545f, 0.2442f, 0.4026f, 0.2442f, 0.0545f};
    const int w = in.width, h = in.height;
    Image tmp(w, h), out(w, h);
    for (int y = 0; y < h; ++y)
        for (int x = 0; x < w; ++x) {
            float acc = 0, wsum = 0;
            for (int d = -2; d <= 2; ++d) {
                const int xx = x + d;
                if (xx < 0 || xx >= w) continue;
                const float v = in.at(xx, y);
                if (!std::isfinite(v)) continue;
                acc += k[d + 2] * v;
                wsum += k[d + 2];
            }
            tmp.at(x, y) = wsum > 0 ? acc / wsum : NAN;
        }
    for (int y = 0; y < h; ++y)
        for (int x = 0; x < w; ++x) {
            float acc = 0, wsum = 0;
            for (int d = -2; d <= 2; ++d) {
                const int yy = y + d;
                if (yy < 0 || yy >= h) continue;
                const float v = tmp.at(x, yy);
                if (!std::isfinite(v)) continue;
                acc += k[d + 2] * v;
                wsum += k[d + 2];
            }
            out.at(x, y) = wsum > 0 ? acc / wsum : NAN;
        }
    return out;
}

struct Fit {
    // target = A * source + b (similarity, optionally mirrored)
    PixelMap map;
    double rms = 0.0;
    int n = 0;
};

// Least-squares similarity (4 parameters) from matched pairs; `mirror` flips source y first.
Fit fit_similarity(const std::vector<Star>& src, const std::vector<Star>& dst,
                   const std::vector<std::pair<int, int>>& pairs, bool mirror) {
    Fit f;
    const std::size_t n = pairs.size();
    if (n < 3) return f;
    double sx = 0, sy = 0, tx = 0, ty = 0;
    for (const auto& p : pairs) {
        const double y = mirror ? -src[static_cast<std::size_t>(p.first)].y : src[static_cast<std::size_t>(p.first)].y;
        sx += src[static_cast<std::size_t>(p.first)].x;
        sy += y;
        tx += dst[static_cast<std::size_t>(p.second)].x;
        ty += dst[static_cast<std::size_t>(p.second)].y;
    }
    sx /= n; sy /= n; tx /= n; ty /= n;
    double num_re = 0, num_im = 0, den = 0;
    for (const auto& p : pairs) {
        const double ax = src[static_cast<std::size_t>(p.first)].x - sx;
        const double ay = (mirror ? -src[static_cast<std::size_t>(p.first)].y : src[static_cast<std::size_t>(p.first)].y) - sy;
        const double bx = dst[static_cast<std::size_t>(p.second)].x - tx;
        const double by = dst[static_cast<std::size_t>(p.second)].y - ty;
        // k = sum(b * conj(a)) / sum(|a|^2)
        num_re += bx * ax + by * ay;
        num_im += by * ax - bx * ay;
        den += ax * ax + ay * ay;
    }
    if (den <= 0) return f;
    const double kr = num_re / den, ki = num_im / den;
    // complex multiply: (kr + i ki)(x + i y') with y' = (mirror ? -y : y)
    const double my = mirror ? -1.0 : 1.0;
    f.map.a11 = kr;
    f.map.a12 = -ki * my;
    f.map.a21 = ki;
    f.map.a22 = kr * my;
    f.map.bx = tx - (kr * sx - ki * sy);
    f.map.by = ty - (ki * sx + kr * sy);
    double ss = 0;
    for (const auto& p : pairs) {
        double px, py;
        f.map.apply(src[static_cast<std::size_t>(p.first)].x, src[static_cast<std::size_t>(p.first)].y, px, py);
        ss += std::pow(px - dst[static_cast<std::size_t>(p.second)].x, 2) + std::pow(py - dst[static_cast<std::size_t>(p.second)].y, 2);
    }
    f.n = static_cast<int>(n);
    f.rms = std::sqrt(ss / n);
    return f;
}

// Nearest-neighbour pairs under a map, within `radius` target pixels (one-to-one, greedy).
std::vector<std::pair<int, int>> pair_up(const std::vector<Star>& src, const std::vector<Star>& dst,
                                         const PixelMap& m, double radius) {
    std::vector<std::pair<int, int>> pairs;
    std::vector<char> used(dst.size(), 0);
    for (std::size_t i = 0; i < src.size(); ++i) {
        double px, py;
        m.apply(src[i].x, src[i].y, px, py);
        int best = -1;
        double bd = radius * radius;
        for (std::size_t j = 0; j < dst.size(); ++j) {
            if (used[j]) continue;
            const double d = (dst[j].x - px) * (dst[j].x - px) + (dst[j].y - py) * (dst[j].y - py);
            if (d < bd) {
                bd = d;
                best = static_cast<int>(j);
            }
        }
        if (best >= 0) {
            used[static_cast<std::size_t>(best)] = 1;
            pairs.emplace_back(static_cast<int>(i), best);
        }
    }
    return pairs;
}

}  // namespace

PixelMap PixelMap::inverse() const {
    const double det = a11 * a22 - a12 * a21;
    if (det == 0.0) throw std::runtime_error("singular pixel map");
    PixelMap r;
    r.a11 = a22 / det;
    r.a12 = -a12 / det;
    r.a21 = -a21 / det;
    r.a22 = a11 / det;
    r.bx = -(r.a11 * bx + r.a12 * by);
    r.by = -(r.a21 * bx + r.a22 * by);
    return r;
}

std::vector<Star> detect_stars(const Image& img, int max_stars, double nsigma) {
    std::vector<Star> stars;
    if (img.empty()) return stars;
    auto samples = detail::sample_pixels(img, 200000);
    if (samples.size() < 100) return stars;
    std::vector<float> s2 = samples;
    const double bg = detail::median_inplace(samples);
    const double sigma = std::max(1e-9, detail::robust_sigma(s2));
    const Image sm = smooth(img);
    // smoothing with a ~1 px Gaussian lowers white noise by about a factor 2.8
    const double thr = bg + nsigma * sigma / 2.8;
    const int w = img.width, h = img.height, R = 4;
    for (int y = R + 1; y < h - R - 1; ++y) {
        for (int x = R + 1; x < w - R - 1; ++x) {
            const float v = sm.at(x, y);
            if (!(v > thr)) continue;
            bool peak = true;
            for (int dy = -2; dy <= 2 && peak; ++dy)
                for (int dx = -2; dx <= 2; ++dx) {
                    if (!dx && !dy) continue;
                    const float u = sm.at(x + dx, y + dy);
                    if (u > v || (u == v && (dy < 0 || (dy == 0 && dx < 0)))) { peak = false; break; }
                }
            if (!peak) continue;
            double f = 0, fx = 0, fy = 0;
            bool bad = false;
            for (int dy = -R; dy <= R && !bad; ++dy)
                for (int dx = -R; dx <= R; ++dx) {
                    const float p = img.at(x + dx, y + dy);
                    if (!std::isfinite(p)) { bad = true; break; }
                    const double q = std::max(0.0, p - bg);
                    f += q; fx += q * dx; fy += q * dy;
                }
            if (bad || f <= 0) continue;
            stars.push_back({x + fx / f, y + fy / f, f});
        }
    }
    std::sort(stars.begin(), stars.end(), [](const Star& a, const Star& b) { return a.flux > b.flux; });
    if (static_cast<int>(stars.size()) > max_stars) stars.resize(static_cast<std::size_t>(max_stars));
    return stars;
}

TransferReport match_images(const Image& source, const Image& target, double scale_hint) {
    TransferReport rep;
    if (scale_hint <= 0.0) {
        const double sx = static_cast<double>(target.width) / source.width;
        const double sy = static_cast<double>(target.height) / source.height;
        // the images may be rotated by 90 degrees relative to each other
        const double sx2 = static_cast<double>(target.width) / source.height;
        const double sy2 = static_cast<double>(target.height) / source.width;
        scale_hint = std::fabs(sx - sy) <= std::fabs(sx2 - sy2) ? 0.5 * (sx + sy) : 0.5 * (sx2 + sy2);
    }
    const auto src = detect_stars(source, 200);
    const auto dst = detect_stars(target, 200);
    if (src.size() < 8 || dst.size() < 8) {
        rep.message = "too few stars detected (" + std::to_string(src.size()) + " / " + std::to_string(dst.size()) + ")";
        return rep;
    }
    const std::size_t nv = std::min<std::size_t>(80, std::min(src.size(), dst.size()));

    // Vote for a translation under each of the 8 axis-aligned orientations and a few
    // scale factors around the hint, then refine the best one.
    struct Trial { int votes = 0; PixelMap map; bool mirror = false; };
    Trial best;
    const double bin = 3.0;
    for (int mirror = 0; mirror < 2; ++mirror) {
        for (int rot = 0; rot < 4; ++rot) {
            for (double sf : {1.0, 0.97, 1.03, 0.94, 1.06}) {
                const double s = scale_hint * sf;
                const double c = std::cos(rot * kPi / 2), sn = std::sin(rot * kPi / 2);
                PixelMap m;
                const double my = mirror ? -1.0 : 1.0;
                m.a11 = s * c; m.a12 = -s * sn * my;
                m.a21 = s * sn; m.a22 = s * c * my;
                std::unordered_map<std::int64_t, int> votes;
                int top = 0;
                std::int64_t top_key = 0;
                for (std::size_t i = 0; i < nv; ++i) {
                    double px, py;
                    m.apply(src[i].x, src[i].y, px, py);
                    for (std::size_t j = 0; j < nv; ++j) {
                        const auto kx = static_cast<std::int64_t>(std::floor((dst[j].x - px) / bin));
                        const auto ky = static_cast<std::int64_t>(std::floor((dst[j].y - py) / bin));
                        const std::int64_t key = (kx + 100000) * 200003 + (ky + 100000);
                        const int v = ++votes[key];
                        if (v > top) { top = v; top_key = key; }
                    }
                }
                if (top > best.votes) {
                    best.votes = top;
                    best.mirror = mirror != 0;
                    best.map = m;
                    const std::int64_t kx = top_key / 200003 - 100000, ky = top_key % 200003 - 100000;
                    best.map.bx = (kx + 0.5) * bin;
                    best.map.by = (ky + 0.5) * bin;
                }
            }
        }
    }
    if (best.votes < 6) {
        rep.message = "no consistent star pattern found between the images";
        return rep;
    }
    // refine: pair up all stars, fit, clip, repeat
    Fit fit;
    fit.map = best.map;
    double radius = 2.0 * bin;
    for (int it = 0; it < 6; ++it) {
        auto pairs = pair_up(src, dst, fit.map, radius);
        if (pairs.size() < 6) break;
        Fit f = fit_similarity(src, dst, pairs, best.mirror);
        if (f.n < 6) break;
        fit = f;
        radius = std::max(1.0, 4.0 * fit.rms);
    }
    if (fit.n < 6) {
        rep.message = "too few matched stars (" + std::to_string(fit.n) + ")";
        return rep;
    }
    rep.ok = true;
    rep.matched = fit.n;
    rep.rms_px = fit.rms;
    rep.map = fit.map;
    rep.mirrored = best.mirror;
    rep.scale = std::sqrt(std::fabs(fit.map.a11 * fit.map.a22 - fit.map.a12 * fit.map.a21));
    rep.rotation_deg = std::atan2(fit.map.a21, fit.map.a11) * 180.0 / kPi;
    rep.message = "matched " + std::to_string(fit.n) + " stars";
    return rep;
}

Wcs transfer_wcs(const Wcs& sw, const Image& source, const Image& target, TransferReport* report,
                 double scale_hint) {
    if (!sw.valid) throw std::runtime_error("the source image has no valid TAN plate solution");
    TransferReport rep = match_images(source, target, scale_hint);
    if (report) *report = rep;
    if (!rep.ok) throw std::runtime_error("cannot match the plate-solved image to the frames: " + rep.message);
    // target -> source (0-based)
    const PixelMap inv = rep.map.inverse();
    // TAN: world offset = CD * (source_1based - crpix_s). Source = inv(target), so
    // CD_t = CD_s * J with J = d(source)/d(target), and crpix_t = map(crpix_s - 1) + 1.
    Wcs w = sw;
    w.cd11 = sw.cd11 * inv.a11 + sw.cd12 * inv.a21;
    w.cd12 = sw.cd11 * inv.a12 + sw.cd12 * inv.a22;
    w.cd21 = sw.cd21 * inv.a11 + sw.cd22 * inv.a21;
    w.cd22 = sw.cd21 * inv.a12 + sw.cd22 * inv.a22;
    double tx, ty;
    rep.map.apply(sw.crpix1 - 1.0, sw.crpix2 - 1.0, tx, ty);
    w.crpix1 = tx + 1.0;
    w.crpix2 = ty + 1.0;
    w.valid = true;
    return w;
}

TransferReport apply_wcs_from_file(const std::string& solved_fits, std::vector<Frame>& frames,
                                   double scale_hint) {
    if (frames.empty()) throw std::runtime_error("no frames");
    FitsImage solved = read_fits(solved_fits);
    const bool distortion = solved.header.has("A_ORDER") || solved.header.has("PV1_1") || solved.header.has("PV2_1");
    const Wcs sw = wcs_from_header(solved.header);
    if (!sw.valid) throw std::runtime_error(solved_fits + " has no TAN plate solution (CTYPE RA---TAN with CD or CDELT keywords)");
    const Frame& ref = frames[frames.size() / 2];
    TransferReport rep;
    const Wcs w = transfer_wcs(sw, solved.image, ref.image, &rep, scale_hint);
    rep.source_has_distortion = distortion;
    for (auto& f : frames) f.wcs = w;
    return rep;
}

}  // namespace fainttrail
