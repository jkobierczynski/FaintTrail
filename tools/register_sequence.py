#!/usr/bin/env python3
# SPDX-License-Identifier: GPL-3.0-or-later
# Copyright 2026 Jurgen Kobierczynski
"""Register a sequence of raw frames for FaintTrail.

FaintTrail expects registered frames (stars on the same pixels). This script is
the simple registration used to prepare the real-data regression sets:

  * optional N x N binning (sum),
  * frames whose pixel data are identical to an earlier frame are dropped
    (some capture programs repeat the last frame when the camera stalls),
  * stars are found in every frame, a coarse offset comes from FFT
    cross-correlation, and a similarity transform (shift, rotation, scale) is
    fitted to the matched stars with outlier clipping,
  * frames are resampled (cubic spline) onto the reference frame; pixels from
    outside the original frame become NaN,
  * DATE-OBS and EXPTIME are copied, and REGDX/REGDY record the offset back to
    the detector, so that FaintTrail can flag objects fixed on the detector.

No dark or flat correction is done. Needs numpy, scipy and astropy.

    python3 tools/register_sequence.py RAW_DIR -o OUT_DIR --bin 2 \
        [--exclude 'name_008.fit'] [--min-exptime 30] [--ref INDEX]
"""
import argparse
import fnmatch
import glob
import hashlib
import os
import sys

import numpy as np
from astropy.io import fits
from scipy import ndimage as ndi


def list_inputs(paths):
    out = []
    for p in paths:
        if os.path.isdir(p):
            fs = [f for f in glob.glob(os.path.join(p, "*")) if f.lower().endswith((".fit", ".fits", ".fts"))]
            out += sorted(fs)
        else:
            out.append(p)
    return out


def binned(a, n):
    if n == 1:
        return a.astype(np.float32)
    h, w = a.shape[0] // n * n, a.shape[1] // n * n
    return a[:h, :w].reshape(h // n, n, w // n, n).sum(axis=(1, 3), dtype=np.float64).astype(np.float32)


def robust(im):
    s = im[::4, ::4]
    s = s[np.isfinite(s)]
    med = np.median(s)
    return med, 1.4826 * np.median(np.abs(s - med))


def find_stars(im, nsig=10.0, nmax=400, saturation=None):
    med, sd = robust(im)
    sm = ndi.gaussian_filter(im - med, 1.2)
    mx = ndi.maximum_filter(sm, 11)
    ys, xs = np.where((sm == mx) & (sm > nsig * sd / 2.0))
    out = []
    for y, x in zip(ys, xs):
        if 8 < y < im.shape[0] - 8 and 8 < x < im.shape[1] - 8:
            if saturation is not None and im[y, x] >= saturation:
                continue
            cut = np.clip(im[y - 5:y + 6, x - 5:x + 6] - med, 0, None)
            yy, xx = np.mgrid[-5:6, -5:6]
            f = cut.sum()
            if f > 0:
                out.append((x + (cut * xx).sum() / f, y + (cut * yy).sum() / f, f))
    out = np.array(out).reshape(-1, 3)
    return out[np.argsort(-out[:, 2])][:nmax]


def normalised(im):
    med, sd = robust(im)
    return np.clip((im - med) / max(sd, 1e-6), 0, 100)


def main():
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("inputs", nargs="+", help="raw FITS files or directories")
    ap.add_argument("-o", "--out", required=True)
    ap.add_argument("--bin", type=int, default=1)
    ap.add_argument("--exclude", action="append", default=[], help="file name pattern to skip (repeatable)")
    ap.add_argument("--min-exptime", type=float, default=0.0, help="skip shorter exposures (test frames)")
    ap.add_argument("--ref", type=int, default=-1, help="index of the reference frame (default: middle)")
    ap.add_argument("--prefix", default="", help="output name prefix (default: keep input names)")
    args = ap.parse_args()

    files = [f for f in list_inputs(args.inputs)
             if not any(fnmatch.fnmatch(os.path.basename(f), pat) for pat in args.exclude)]
    frames, meta, seen = [], [], {}
    for f in files:
        with fits.open(f, memmap=False) as hd:
            h = hd[0].header
            data = hd[0].data
        if data is None or data.ndim != 2:
            print(f"skip {f}: no 2-D image", file=sys.stderr)
            continue
        exp = float(h.get("EXPTIME", h.get("EXPOSURE", 0.0)))
        if exp < args.min_exptime:
            print(f"skip {os.path.basename(f)}: exposure {exp} s", file=sys.stderr)
            continue
        digest = hashlib.sha1(np.ascontiguousarray(data).tobytes()).hexdigest()
        if digest in seen:
            print(f"skip {os.path.basename(f)}: same pixels as {seen[digest]}", file=sys.stderr)
            continue
        seen[digest] = os.path.basename(f)
        if "DATE-OBS" not in h:
            print(f"skip {f}: no DATE-OBS", file=sys.stderr)
            continue
        frames.append(binned(np.asarray(data, dtype=np.float32), args.bin))
        meta.append((f, h["DATE-OBS"], exp))
    if len(frames) < 3:
        sys.exit("fewer than 3 usable frames")
    print(f"{len(frames)} frames, {frames[0].shape[1]} x {frames[0].shape[0]} px after {args.bin}x{args.bin} binning")

    ref_i = args.ref if args.ref >= 0 else len(frames) // 2
    ref = find_stars(frames[ref_i])
    if len(ref) < 10:
        sys.exit(f"only {len(ref)} stars in the reference frame")
    H, W = frames[0].shape
    F0 = np.fft.rfft2(normalised(frames[ref_i]))
    yy, xx = np.mgrid[0:H, 0:W]
    os.makedirs(args.out, exist_ok=True)
    for i, (im, (path, date, exp)) in enumerate(zip(frames, meta)):
        X = np.fft.irfft2(np.fft.rfft2(normalised(im)) * np.conj(F0), s=(H, W))
        cy, cx = np.unravel_index(np.argmax(X), X.shape)
        cy = cy - H if cy > H / 2 else cy
        cx = cx - W if cx > W / 2 else cx
        src = find_stars(im)
        d = np.hypot(src[:, None, 0] - cx - ref[None, :, 0], src[:, None, 1] - cy - ref[None, :, 1])
        j = np.argmin(d, axis=1)
        ok = d[np.arange(len(src)), j] < 3
        A, B = src[ok, :2], ref[j[ok], :2]
        if len(A) < 6:
            print(f"skip {os.path.basename(path)}: only {len(A)} matched stars", file=sys.stderr)
            continue
        keep = np.ones(len(A), bool)
        for _ in range(4):  # B = k A + t in complex form, with clipping
            za, zb = A[keep, 0] + 1j * A[keep, 1], B[keep, 0] + 1j * B[keep, 1]
            ma, mb = za.mean(), zb.mean()
            k = np.sum((zb - mb) * np.conj(za - ma)) / np.sum(np.abs(za - ma) ** 2)
            t = mb - k * ma
            res = np.abs((A[:, 0] + 1j * A[:, 1]) * k + t - (B[:, 0] + 1j * B[:, 1]))
            keep = res < max(3 * 1.4826 * np.median(res[keep]), 0.2)
        rms = float(np.sqrt(np.mean(res[keep] ** 2)))
        zin = (xx + 1j * yy - t) / k  # registered pixel -> input pixel
        out = ndi.map_coordinates(im, [zin.imag, zin.real], order=3, mode="constant", cval=np.nan).astype(np.float32)
        c = W / 2 + 1j * H / 2
        off = (c - t) / k - c  # detector = registered + off, at the image centre
        hdr = fits.Header()
        hdr["DATE-OBS"] = (date, "UTC start of exposure")
        hdr["EXPTIME"] = (exp, "exposure time [s]")
        hdr["REGDX"] = (round(float(off.real), 4), "detector x = registered x + REGDX")
        hdr["REGDY"] = (round(float(off.imag), 4), "detector y = registered y + REGDY")
        hdr["HISTORY"] = f"register_sequence.py: {os.path.basename(path)} bin {args.bin}, {int(keep.sum())} stars, rms {rms:.3f} px"
        name = args.prefix + os.path.splitext(os.path.basename(path))[0] + ".fits" if args.prefix else \
            os.path.splitext(os.path.basename(path))[0] + ".fits"
        fits.PrimaryHDU(out, hdr).writeto(os.path.join(args.out, name), overwrite=True)
        print(f"{os.path.basename(path)}  {date}  shift ({off.real:7.2f}, {off.imag:7.2f})  "
              f"rot {np.degrees(np.angle(k)):7.4f} deg  scale {abs(k):.5f}  {int(keep.sum()):3d} stars  rms {rms:.3f} px")


if __name__ == "__main__":
    main()
