// SPDX-License-Identifier: GPL-3.0-or-later
// Copyright 2026 Jurgen Kobierczynski
#include "fainttrail/fits.hpp"

#include <algorithm>
#include <cctype>
#include <cmath>
#include <cstdint>
#include <cstdio>
#include <cstring>
#include <fstream>
#include <iomanip>
#include <limits>
#include <sstream>
#include <stdexcept>

namespace fainttrail {

namespace {

constexpr std::size_t kBlock = 2880;
constexpr std::size_t kCard = 80;

std::string trim(const std::string& s) {
    std::size_t a = 0, b = s.size();
    while (a < b && std::isspace(static_cast<unsigned char>(s[a]))) ++a;
    while (b > a && std::isspace(static_cast<unsigned char>(s[b - 1]))) --b;
    return s.substr(a, b - a);
}

std::string upper(std::string s) {
    for (auto& c : s) c = static_cast<char>(std::toupper(static_cast<unsigned char>(c)));
    return s;
}

// Parse one 80-character card.
bool parse_card(const char* raw, FitsHeader::Card& card, bool& is_end) {
    std::string line(raw, kCard);
    std::string key = trim(line.substr(0, 8));
    is_end = (key == "END");
    if (is_end) return false;
    card = {};
    card.key = key;
    if (key.empty()) return false;
    if (key == "HISTORY" || key == "COMMENT") {
        card.value = trim(line.substr(8));
        card.is_string = true;
        return true;
    }
    if (line.size() < 10 || line[8] != '=' ) {
        card.value = trim(line.substr(8));
        card.is_string = true;
        return true;
    }
    std::string rest = line.substr(10);
    std::size_t i = 0;
    while (i < rest.size() && rest[i] == ' ') ++i;
    if (i < rest.size() && rest[i] == '\'') {
        // quoted string, '' is an escaped quote
        std::string v;
        ++i;
        while (i < rest.size()) {
            if (rest[i] == '\'') {
                if (i + 1 < rest.size() && rest[i + 1] == '\'') {
                    v += '\'';
                    i += 2;
                    continue;
                }
                ++i;
                break;
            }
            v += rest[i++];
        }
        // trailing spaces inside strings are not significant
        while (!v.empty() && v.back() == ' ') v.pop_back();
        card.value = v;
        card.is_string = true;
        std::size_t slash = rest.find('/', i);
        if (slash != std::string::npos) card.comment = trim(rest.substr(slash + 1));
    } else {
        std::size_t slash = rest.find('/');
        card.value = trim(rest.substr(0, slash));
        if (slash != std::string::npos) card.comment = trim(rest.substr(slash + 1));
    }
    return true;
}

std::string format_card(const std::string& key, const std::string& value_field,
                        const std::string& comment) {
    std::string k = upper(key).substr(0, 8);
    std::string line = k;
    line.resize(8, ' ');
    if (k == "HISTORY" || k == "COMMENT") {
        line += value_field;
    } else {
        line += "= ";
        line += value_field;
        if (!comment.empty()) {
            if (line.size() < 31) line.resize(31, ' ');
            line += " / ";
            line += comment;
        }
    }
    line.resize(kCard, ' ');
    return line;
}

std::string format_number(double v) {
    char buf[64];
    if (std::isfinite(v) && v == std::floor(v) && std::fabs(v) < 1e15) {
        std::snprintf(buf, sizeof buf, "%.1f", v);
    } else {
        std::snprintf(buf, sizeof buf, "%.15G", v);
        // FITS wants an E exponent and a decimal point in reals
        std::string s = buf;
        if (s.find('.') == std::string::npos && s.find('E') != std::string::npos) {
            s.insert(s.find('E'), ".0");
        }
        std::snprintf(buf, sizeof buf, "%s", s.c_str());
    }
    std::string s = buf;
    if (s.size() < 20) s = std::string(20 - s.size(), ' ') + s;
    return s;
}

std::string format_string(const std::string& v) {
    std::string esc;
    for (char c : v) {
        esc += c;
        if (c == '\'') esc += '\'';
    }
    if (esc.size() < 8) esc.resize(8, ' ');
    return "'" + esc + "'";
}

template <typename T>
T read_be(const unsigned char* p) {
    std::uint64_t u = 0;
    for (std::size_t i = 0; i < sizeof(T); ++i) u = (u << 8) | p[i];
    T out;
    if constexpr (sizeof(T) == 1) {
        std::uint8_t v = static_cast<std::uint8_t>(u);
        std::memcpy(&out, &v, 1);
    } else if constexpr (sizeof(T) == 2) {
        std::uint16_t v = static_cast<std::uint16_t>(u);
        std::memcpy(&out, &v, 2);
    } else if constexpr (sizeof(T) == 4) {
        std::uint32_t v = static_cast<std::uint32_t>(u);
        std::memcpy(&out, &v, 4);
    } else {
        std::memcpy(&out, &u, 8);
    }
    return out;
}

struct Hdu {
    FitsHeader header;
    std::size_t data_offset = 0;
    std::size_t data_bytes = 0;
};

// Read a header starting at `offset`; returns false at end of file.
bool read_hdu_header(std::istream& in, std::size_t offset, Hdu& hdu) {
    in.clear();
    in.seekg(static_cast<std::streamoff>(offset));
    hdu = {};
    std::vector<char> block(kBlock);
    std::size_t pos = offset;
    bool done = false;
    while (!done) {
        if (!in.read(block.data(), kBlock)) return false;
        pos += kBlock;
        for (std::size_t c = 0; c < kBlock / kCard; ++c) {
            FitsHeader::Card card;
            bool is_end = false;
            if (parse_card(block.data() + c * kCard, card, is_end)) hdu.header.cards.push_back(card);
            if (is_end) {
                done = true;
                break;
            }
        }
    }
    hdu.data_offset = pos;
    const long long bitpix = hdu.header.get_int("BITPIX", 0);
    const long long naxis = hdu.header.get_int("NAXIS", 0);
    std::size_t n = naxis > 0 ? 1 : 0;
    for (long long i = 1; i <= naxis; ++i)
        n *= static_cast<std::size_t>(hdu.header.get_int("NAXIS" + std::to_string(i), 0));
    const long long pcount = hdu.header.get_int("PCOUNT", 0);
    const long long gcount = hdu.header.get_int("GCOUNT", 1);
    std::size_t bytes = static_cast<std::size_t>(std::llabs(bitpix) / 8) *
                        static_cast<std::size_t>(gcount) * (n + static_cast<std::size_t>(pcount));
    hdu.data_bytes = bytes;
    return true;
}

}  // namespace

// ------------------------------------------------------------------ header
const FitsHeader::Card* FitsHeader::find(const std::string& key) const {
    const std::string k = upper(key);
    for (const auto& c : cards)
        if (c.key == k) return &c;
    return nullptr;
}

double FitsHeader::get_double(const std::string& key, double fallback) const {
    const Card* c = find(key);
    if (!c || c->value.empty()) return fallback;
    std::string v = c->value;
    std::replace(v.begin(), v.end(), 'D', 'E');  // Fortran exponents
    std::replace(v.begin(), v.end(), 'd', 'e');
    char* end = nullptr;
    double d = std::strtod(v.c_str(), &end);
    if (end == v.c_str()) return fallback;
    return d;
}

long long FitsHeader::get_int(const std::string& key, long long fallback) const {
    const Card* c = find(key);
    if (!c || c->value.empty()) return fallback;
    char* end = nullptr;
    long long v = std::strtoll(c->value.c_str(), &end, 10);
    if (end == c->value.c_str()) return fallback;
    return v;
}

std::string FitsHeader::get_string(const std::string& key, const std::string& fallback) const {
    const Card* c = find(key);
    return c ? c->value : fallback;
}

static void upsert(FitsHeader& h, const FitsHeader::Card& card) {
    for (auto& c : h.cards) {
        if (c.key == card.key) {
            c = card;
            return;
        }
    }
    h.cards.push_back(card);
}

void FitsHeader::set(const std::string& key, double value, const std::string& comment) {
    upsert(*this, {upper(key), format_number(value), comment, false});
}

void FitsHeader::set_int(const std::string& key, long long value, const std::string& comment) {
    std::string s = std::to_string(value);
    if (s.size() < 20) s = std::string(20 - s.size(), ' ') + s;
    upsert(*this, {upper(key), s, comment, false});
}

void FitsHeader::set(const std::string& key, const std::string& value, const std::string& comment) {
    upsert(*this, {upper(key), value, comment, true});
}

void FitsHeader::add_history(const std::string& text) {
    cards.push_back({"HISTORY", text, {}, true});
}

// ------------------------------------------------------------------ read
FitsImage read_fits(const std::string& path) {
    std::ifstream in(path, std::ios::binary);
    if (!in) throw std::runtime_error("cannot open FITS file: " + path);
    in.seekg(0, std::ios::end);
    const std::size_t file_size = static_cast<std::size_t>(in.tellg());

    std::size_t offset = 0;
    for (int hdu_index = 0; hdu_index < 64 && offset < file_size; ++hdu_index) {
        Hdu hdu;
        if (!read_hdu_header(in, offset, hdu)) break;
        if (hdu_index == 0 && hdu.header.get_string("SIMPLE") != "T")
            throw std::runtime_error("not a FITS file (SIMPLE != T): " + path);
        const long long naxis = hdu.header.get_int("NAXIS", 0);
        const std::string xt = upper(hdu.header.get_string("XTENSION", hdu_index == 0 ? "PRIMARY" : ""));
        const bool image_hdu = (hdu_index == 0 || xt == "IMAGE");
        if (image_hdu && naxis >= 2) {
            const int w = static_cast<int>(hdu.header.get_int("NAXIS1", 0));
            const int h = static_cast<int>(hdu.header.get_int("NAXIS2", 0));
            const long long bitpix = hdu.header.get_int("BITPIX", 0);
            if (w <= 0 || h <= 0) throw std::runtime_error("invalid image size in " + path);
            const std::size_t bpp = static_cast<std::size_t>(std::llabs(bitpix) / 8);
            const std::size_t npix = static_cast<std::size_t>(w) * h;
            if (hdu.data_offset + npix * bpp > file_size)
                throw std::runtime_error("truncated FITS data in " + path);
            std::vector<unsigned char> raw(npix * bpp);
            in.clear();
            in.seekg(static_cast<std::streamoff>(hdu.data_offset));
            in.read(reinterpret_cast<char*>(raw.data()), static_cast<std::streamsize>(raw.size()));
            if (!in) throw std::runtime_error("failed reading FITS data in " + path);
            const double bscale = hdu.header.get_double("BSCALE", 1.0);
            const double bzero = hdu.header.get_double("BZERO", 0.0);
            const bool has_blank = hdu.header.has("BLANK") && bitpix > 0;
            const long long blank = hdu.header.get_int("BLANK", 0);

            FitsImage out;
            out.image = Image(w, h);
            out.header = std::move(hdu.header);
            for (std::size_t i = 0; i < npix; ++i) {
                const unsigned char* p = raw.data() + i * bpp;
                double v = 0.0;
                bool is_blank = false;
                switch (bitpix) {
                    case 8: v = p[0]; is_blank = has_blank && p[0] == blank; break;
                    case 16: { auto s = read_be<std::int16_t>(p); v = s; is_blank = has_blank && s == blank; } break;
                    case 32: { auto s = read_be<std::int32_t>(p); v = s; is_blank = has_blank && s == blank; } break;
                    case 64: { auto s = read_be<std::int64_t>(p); v = static_cast<double>(s); is_blank = has_blank && s == blank; } break;
                    case -32: v = read_be<float>(p); break;
                    case -64: v = read_be<double>(p); break;
                    default: throw std::runtime_error("unsupported BITPIX in " + path);
                }
                out.image.data[i] = is_blank ? NAN : static_cast<float>(bzero + bscale * v);
            }
            return out;
        }
        std::size_t data_padded = (hdu.data_bytes + kBlock - 1) / kBlock * kBlock;
        offset = hdu.data_offset + data_padded;
    }
    throw std::runtime_error("no 2-D image found in " + path);
}

// ------------------------------------------------------------------ write
void write_fits(const std::string& path, const Image& image, const FitsHeader& header, int bitpix) {
    if (bitpix != -32 && bitpix != 16) throw std::runtime_error("write_fits: only BITPIX -32 and 16 are supported");
    std::ofstream out(path, std::ios::binary);
    if (!out) throw std::runtime_error("cannot create FITS file: " + path);
    std::string hdr;
    hdr += format_card("SIMPLE", std::string(19, ' ') + "T", "conforms to FITS standard");
    // 16-bit scaling: robust noise from a subsample, step <= noise / 4
    double bscale = 1.0, bzero = 0.0;
    if (bitpix == 16) {
        std::vector<float> v;
        float lo = std::numeric_limits<float>::infinity(), hi = -lo;
        const std::size_t step = std::max<std::size_t>(1, image.size() / 200000);
        for (std::size_t i = 0; i < image.size(); ++i) {
            const float x = image.data[i];
            if (!std::isfinite(x)) continue;
            lo = std::min(lo, x);
            hi = std::max(hi, x);
            if (i % step == 0) v.push_back(x);
        }
        if (v.empty()) {
            lo = 0.0f;
            hi = 1.0f;
            v.push_back(0.0f);
        }
        auto mid = v.begin() + static_cast<std::ptrdiff_t>(v.size() / 2);
        std::nth_element(v.begin(), mid, v.end());
        const double med = *mid;
        for (auto& x : v) x = std::fabs(x - static_cast<float>(med));
        std::nth_element(v.begin(), mid, v.end());
        const double sigma = 1.4826 * *mid;
        bscale = (static_cast<double>(hi) - lo) / 65530.0;
        if (sigma > 0.0) bscale = std::min(bscale, sigma / 4.0);
        if (!(bscale > 0.0)) bscale = 1.0;
        bzero = lo + 32767.0 * bscale;  // lowest value maps to -32767; -32768 is BLANK
    }
    hdr += format_card("BITPIX", std::string(bitpix == 16 ? 18 : 17, ' ') + std::to_string(bitpix),
                       bitpix == 16 ? "16-bit scaled integers" : "IEEE single precision");
    hdr += format_card("NAXIS", std::string(19, ' ') + "2", "");
    {
        std::string s = std::to_string(image.width);
        hdr += format_card("NAXIS1", std::string(20 - s.size(), ' ') + s, "");
        s = std::to_string(image.height);
        hdr += format_card("NAXIS2", std::string(20 - s.size(), ' ') + s, "");
    }
    for (const auto& c : header.cards) {
        static const char* skip[] = {"SIMPLE", "BITPIX", "NAXIS", "NAXIS1", "NAXIS2",
                                     "NAXIS3", "BSCALE", "BZERO", "BLANK", "END", "EXTEND"};
        bool skipit = false;
        for (const char* s : skip)
            if (c.key == s) skipit = true;
        if (skipit) continue;
        if (c.key == "HISTORY" || c.key == "COMMENT") {
            std::string text = c.value;
            do {
                hdr += format_card(c.key, text.substr(0, 72), "");
                text = text.size() > 72 ? text.substr(72) : std::string();
            } while (!text.empty());
        } else {
            hdr += format_card(c.key, c.is_string ? format_string(c.value) : c.value, c.comment);
        }
    }
    if (bitpix == 16) {
        char b[32];
        std::snprintf(b, sizeof b, "%20.12E", bscale);
        hdr += format_card("BSCALE", b, "physical = BZERO + BSCALE * stored");
        std::snprintf(b, sizeof b, "%20.12E", bzero);
        hdr += format_card("BZERO", b, "");
        hdr += format_card("BLANK", std::string(14, ' ') + "-32768", "undefined pixels");
    }
    hdr += std::string("END") + std::string(kCard - 3, ' ');
    hdr.resize((hdr.size() + kBlock - 1) / kBlock * kBlock, ' ');
    out.write(hdr.data(), static_cast<std::streamsize>(hdr.size()));
    if (bitpix == 16) {
        std::vector<unsigned char> buf(image.size() * 2);
        for (std::size_t i = 0; i < image.size(); ++i) {
            const float x = image.data[i];
            long q = -32768;
            if (std::isfinite(x)) q = std::clamp(std::lround((x - bzero) / bscale), -32767L, 32767L);
            const auto u = static_cast<std::uint16_t>(static_cast<std::int16_t>(q));
            buf[i * 2 + 0] = static_cast<unsigned char>(u >> 8);
            buf[i * 2 + 1] = static_cast<unsigned char>(u);
        }
        buf.resize((buf.size() + kBlock - 1) / kBlock * kBlock, 0);
        out.write(reinterpret_cast<const char*>(buf.data()), static_cast<std::streamsize>(buf.size()));
        if (!out) throw std::runtime_error("failed writing FITS file: " + path);
        return;
    }

    std::vector<unsigned char> buf(image.size() * 4);
    for (std::size_t i = 0; i < image.size(); ++i) {
        std::uint32_t u;
        std::memcpy(&u, &image.data[i], 4);
        buf[i * 4 + 0] = static_cast<unsigned char>(u >> 24);
        buf[i * 4 + 1] = static_cast<unsigned char>(u >> 16);
        buf[i * 4 + 2] = static_cast<unsigned char>(u >> 8);
        buf[i * 4 + 3] = static_cast<unsigned char>(u);
    }
    buf.resize((buf.size() + kBlock - 1) / kBlock * kBlock, 0);
    out.write(reinterpret_cast<const char*>(buf.data()), static_cast<std::streamsize>(buf.size()));
    if (!out) throw std::runtime_error("failed writing FITS file: " + path);
}

// ------------------------------------------------------------------ time
static long long days_from_civil(long long y, unsigned m, unsigned d) {
    y -= m <= 2;
    const long long era = (y >= 0 ? y : y - 399) / 400;
    const unsigned yoe = static_cast<unsigned>(y - era * 400);
    const unsigned doy = (153 * (m + (m > 2 ? -3 : 9)) + 2) / 5 + d - 1;
    const unsigned doe = yoe * 365 + yoe / 4 - yoe / 100 + doy;
    return era * 146097 + static_cast<long long>(doe) - 719468;  // days since 1970-01-01
}

double mjd_from_iso(const std::string& iso) {
    int Y = 0, M = 0, D = 0, h = 0, m = 0;
    double s = 0.0;
    std::string t = iso;
    std::replace(t.begin(), t.end(), 'T', ' ');
    int n = std::sscanf(t.c_str(), "%d-%d-%d %d:%d:%lf", &Y, &M, &D, &h, &m, &s);
    if (n < 3) throw std::runtime_error("cannot parse date: " + iso);
    const long long days = days_from_civil(Y, static_cast<unsigned>(M), static_cast<unsigned>(D));
    return 40587.0 + static_cast<double>(days) + (h + (m + s / 60.0) / 60.0) / 24.0;
}

std::string iso_from_mjd(double mjd, int decimals) {
    double days = std::floor(mjd);
    double frac_s = (mjd - days) * 86400.0;
    const double scale = std::pow(10.0, decimals);
    frac_s = std::round(frac_s * scale) / scale;
    if (frac_s >= 86400.0) {
        frac_s -= 86400.0;
        days += 1.0;
    }
    long long z = static_cast<long long>(days) - 40587 + 719468;
    const long long era = (z >= 0 ? z : z - 146096) / 146097;
    const unsigned doe = static_cast<unsigned>(z - era * 146097);
    const unsigned yoe = (doe - doe / 1460 + doe / 36524 - doe / 146096) / 365;
    long long y = static_cast<long long>(yoe) + era * 400;
    const unsigned doy = doe - (365 * yoe + yoe / 4 - yoe / 100);
    const unsigned mp = (5 * doy + 2) / 153;
    const unsigned d = doy - (153 * mp + 2) / 5 + 1;
    const unsigned mo = mp + (mp < 10 ? 3 : -9);
    y += (mo <= 2);
    const int hh = static_cast<int>(frac_s / 3600.0);
    const int mm = static_cast<int>((frac_s - hh * 3600.0) / 60.0);
    const double ss = frac_s - hh * 3600.0 - mm * 60.0;
    char buf[64];
    if (decimals > 0)
        std::snprintf(buf, sizeof buf, "%04lld-%02u-%02uT%02d:%02d:%0*.*f", y, mo, d, hh, mm,
                      decimals + 3, decimals, ss);
    else
        std::snprintf(buf, sizeof buf, "%04lld-%02u-%02uT%02d:%02d:%02d", y, mo, d, hh, mm,
                      static_cast<int>(std::lround(ss)));
    return buf;
}

// ------------------------------------------------------------------ frames
Wcs wcs_from_header(const FitsHeader& h) {
    Wcs w;
    const std::string ct1 = upper(h.get_string("CTYPE1"));
    if (ct1.find("TAN") == std::string::npos || !h.has("CRVAL1") || !h.has("CRPIX1")) return w;
    w.crval1 = h.get_double("CRVAL1", 0.0);
    w.crval2 = h.get_double("CRVAL2", 0.0);
    w.crpix1 = h.get_double("CRPIX1", 0.0);
    w.crpix2 = h.get_double("CRPIX2", 0.0);
    if (h.has("CD1_1")) {
        w.cd11 = h.get_double("CD1_1", 0.0);
        w.cd12 = h.get_double("CD1_2", 0.0);
        w.cd21 = h.get_double("CD2_1", 0.0);
        w.cd22 = h.get_double("CD2_2", 0.0);
    } else if (h.has("CDELT1")) {
        const double c1 = h.get_double("CDELT1", 0.0), c2 = h.get_double("CDELT2", 0.0);
        if (h.has("PC1_1")) {
            w.cd11 = c1 * h.get_double("PC1_1", 1.0);
            w.cd12 = c1 * h.get_double("PC1_2", 0.0);
            w.cd21 = c2 * h.get_double("PC2_1", 0.0);
            w.cd22 = c2 * h.get_double("PC2_2", 1.0);
        } else {
            const double rot = h.get_double("CROTA2", 0.0) * 3.14159265358979323846 / 180.0;
            w.cd11 = c1 * std::cos(rot);
            w.cd12 = -c2 * std::sin(rot);
            w.cd21 = c1 * std::sin(rot);
            w.cd22 = c2 * std::cos(rot);
        }
    }
    w.valid = (w.cd11 * w.cd22 - w.cd12 * w.cd21) != 0.0;
    return w;
}

Frame frame_from_fits(FitsImage&& fits, const std::string& source) {
    Frame f;
    const FitsHeader& h = fits.header;
    f.exptime_s = h.get_double("EXPTIME", h.get_double("EXPOSURE", h.get_double("EXP_TIME", 0.0)));
    double mjd = NAN;
    if (h.has("MJD-MID")) {
        mjd = h.get_double("MJD-MID", NAN);
    } else if (h.has("MJD-OBS")) {
        mjd = h.get_double("MJD-OBS", NAN) + f.exptime_s / 172800.0;
    } else if (h.has("DATE-OBS")) {
        std::string d = h.get_string("DATE-OBS");
        if (d.find('T') == std::string::npos && h.has("TIME-OBS")) d += "T" + h.get_string("TIME-OBS");
        mjd = mjd_from_iso(d) + f.exptime_s / 172800.0;
    } else if (h.has("JD")) {
        mjd = h.get_double("JD", NAN) - 2400000.5 + f.exptime_s / 172800.0;
    }
    if (!std::isfinite(mjd)) throw std::runtime_error("no usable observation time in " + source);
    f.mjd_mid = mjd;

    f.wcs = wcs_from_header(h);
    if (h.has("REGDX") && h.has("REGDY")) {
        f.has_reg_offset = true;
        f.reg_dx = h.get_double("REGDX", 0.0);
        f.reg_dy = h.get_double("REGDY", 0.0);
    }
    f.image = std::move(fits.image);
    f.source = source;
    return f;
}

std::vector<Frame> load_frames(const std::vector<std::string>& paths) {
    std::vector<Frame> frames;
    frames.reserve(paths.size());
    for (const auto& p : paths) frames.push_back(frame_from_fits(read_fits(p), p));
    if (frames.empty()) return frames;
    for (const auto& f : frames)
        if (f.image.width != frames[0].image.width || f.image.height != frames[0].image.height)
            throw std::runtime_error("frame size mismatch: " + f.source);
    std::stable_sort(frames.begin(), frames.end(),
                     [](const Frame& a, const Frame& b) { return a.mjd_mid < b.mjd_mid; });
    return frames;
}

void write_frame(const std::string& path, const Frame& frame, const FitsHeader& extra, int bitpix) {
    FitsHeader h;
    h.set("DATE-OBS", iso_from_mjd(frame.mjd_mid - frame.exptime_s / 172800.0, 3),
          "UTC start of exposure");
    h.set("MJD-MID", frame.mjd_mid, "UTC MJD at mid-exposure");
    h.set("EXPTIME", frame.exptime_s, "exposure time [s]");
    if (frame.wcs.valid) {
        h.set("CTYPE1", std::string("RA---TAN"));
        h.set("CTYPE2", std::string("DEC--TAN"));
        h.set("CRVAL1", frame.wcs.crval1);
        h.set("CRVAL2", frame.wcs.crval2);
        h.set("CRPIX1", frame.wcs.crpix1);
        h.set("CRPIX2", frame.wcs.crpix2);
        h.set("CD1_1", frame.wcs.cd11);
        h.set("CD1_2", frame.wcs.cd12);
        h.set("CD2_1", frame.wcs.cd21);
        h.set("CD2_2", frame.wcs.cd22);
        h.set("RADESYS", std::string("ICRS"));
    }
    if (frame.has_reg_offset) {
        h.set("REGDX", frame.reg_dx, "detector x = registered x + REGDX");
        h.set("REGDY", frame.reg_dy, "detector y = registered y + REGDY");
    }
    for (const auto& c : extra.cards) {
        if (c.key == "HISTORY" || c.key == "COMMENT") h.cards.push_back(c);
        else upsert(h, c);
    }
    write_fits(path, frame.image, h, bitpix);
}

Frame crop_frame(const Frame& frame, int x0, int y0, int w, int h) {
    const Image& src = frame.image;
    if (w <= 0 || h <= 0 || x0 < 0 || y0 < 0 || x0 + w > src.width || y0 + h > src.height)
        throw std::runtime_error("crop rectangle " + std::to_string(w) + "x" + std::to_string(h) + "+" +
                                 std::to_string(x0) + "+" + std::to_string(y0) + " leaves the " +
                                 std::to_string(src.width) + "x" + std::to_string(src.height) + " image");
    Frame out = frame;
    out.image = Image(w, h);
    for (int y = 0; y < h; ++y)
        for (int x = 0; x < w; ++x) out.image.at(x, y) = src.at(x0 + x, y0 + y);
    if (out.wcs.valid) {
        out.wcs.crpix1 -= x0;
        out.wcs.crpix2 -= y0;
    }
    if (out.has_reg_offset) {  // detector = registered + reg; registered = cropped + x0
        out.reg_dx += x0;
        out.reg_dy += y0;
    }
    return out;
}

}  // namespace fainttrail
