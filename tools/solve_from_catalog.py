#!/usr/bin/env python3
# SPDX-License-Identifier: GPL-3.0-or-later
# Copyright 2026 Jurgen Kobierczynski
"""Plate-solve a stack against a small star list, when the field centre is known.

For fields where a blind solver is not at hand: give an approximate centre (for
example from a planetarium program), an approximate scale, and a catalogue
extract (RA, Dec in degrees and magnitude, whitespace separated; e.g. Gaia DR3
from VizieR, 15-20' radius, brightest 50-100 stars). The script

  1. median-stacks the registered frames (or takes one image),
  2. finds the pattern by brute force: mirror, rotation in 1 deg steps and scale
     within +-10 %, voting for the translation,
  3. refines with an affine fit and outlier clipping, then fits a TAN WCS,
  4. writes the stack with that WCS, ready for `fainttrail --wcs-from`.

    python3 tools/solve_from_catalog.py FRAMES_DIR_OR_FITS catalogue.tsv \
        --ra 234.984 --dec -6.618 --scale 3.0 -o stack_solved.fits

Needs numpy, scipy and astropy.
"""
import argparse
import glob
import os
import sys

import numpy as np
from astropy.io import fits
from astropy.coordinates import SkyCoord
from astropy.wcs.utils import fit_wcs_from_points
import astropy.units as u
from scipy import ndimage as ndi

D2R = np.pi / 180.0


def load_stack(path):
    if os.path.isdir(path):
        fs = sorted(f for f in glob.glob(os.path.join(path, "*")) if f.lower().endswith((".fit", ".fits", ".fts")))
        step = max(1, len(fs) // 25)
        cube = np.array([fits.getdata(f).astype(np.float32) for f in fs[::step]])
        with np.errstate(all="ignore"):
            st = np.nanmedian(cube, axis=0)
        hdr = fits.getheader(fs[len(fs) // 2])
    else:
        st = fits.getdata(path).astype(np.float32)
        hdr = fits.getheader(path)
    fill = np.nanmedian(st)
    return np.where(np.isfinite(st), st, fill).astype(np.float32), hdr


def stars(im, nsig=6.0):
    med = np.median(im)
    sd = 1.4826 * np.median(np.abs(im - med))
    sm = ndi.gaussian_filter(im - med, 1.0)
    mx = ndi.maximum_filter(sm, 9)
    ys, xs = np.where((sm == mx) & (sm > nsig * sd / 2))
    out = []
    for y, x in zip(ys, xs):
        if 5 < y < im.shape[0] - 5 and 5 < x < im.shape[1] - 5:
            cut = np.clip(im[y - 4:y + 5, x - 4:x + 5] - med, 0, None)
            yy, xx = np.mgrid[-4:5, -4:5]
            f = cut.sum()
            out.append((x + (cut * xx).sum() / f, y + (cut * yy).sum() / f, sm[y, x]))
    out = np.array(out).reshape(-1, 3)
    return out[np.argsort(-out[:, 2])]


def tangent(ra, dec, ra0, dec0):
    ra, dec, a0, d0 = ra * D2R, dec * D2R, ra0 * D2R, dec0 * D2R
    cc = np.sin(d0) * np.sin(dec) + np.cos(d0) * np.cos(dec) * np.cos(ra - a0)
    xi = np.cos(dec) * np.sin(ra - a0) / cc
    eta = (np.cos(d0) * np.sin(dec) - np.sin(d0) * np.cos(dec) * np.cos(ra - a0)) / cc
    return np.c_[xi, eta] / D2R * 3600.0  # arcsec


def main():
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("image", help="directory of registered frames, or one FITS image")
    ap.add_argument("catalogue", help="text file: ra_deg dec_deg mag per line")
    ap.add_argument("--ra", type=float, required=True, help="approximate centre RA (deg)")
    ap.add_argument("--dec", type=float, required=True, help="approximate centre Dec (deg)")
    ap.add_argument("--scale", type=float, required=True, help="approximate scale (arcsec/px)")
    ap.add_argument("-o", "--out", required=True)
    args = ap.parse_args()

    im, src_hdr = load_stack(args.image)
    P = stars(im)
    cat = np.loadtxt(args.catalogue, ndmin=2)
    cat = cat[np.argsort(cat[:, 2])]
    C = tangent(cat[:, 0], cat[:, 1], args.ra, args.dec)
    if len(P) < 8 or len(C) < 8:
        sys.exit(f"too few stars: {len(P)} in the image, {len(C)} in the catalogue")
    Pb, Cb = P[:80, :2], C[:40]
    H, W = im.shape
    best = (0,)
    for flip in (1, -1):
        for scale in np.arange(0.9, 1.1, 0.01) * args.scale:
            for rot in np.arange(0.0, 360.0, 1.0):
                t = rot * D2R
                R = np.array([[np.cos(t), -np.sin(t)], [np.sin(t), np.cos(t)]])
                Q = (Cb * [flip, 1]) @ R.T / scale
                D = (Pb[:, None, :] - Q[None, :, :]).reshape(-1, 2)
                hist, xe, ye = np.histogram2d(D[:, 0], D[:, 1], bins=[np.arange(-W, 2 * W, 4), np.arange(-H, 2 * H, 4)])
                k = hist.max()
                if k > best[0]:
                    i, j = np.unravel_index(hist.argmax(), hist.shape)
                    best = (k, flip, scale, rot, xe[i] + 2, ye[j] + 2)
    k, flip, scale, rot, tx, ty = best
    print(f"pattern: {int(k)} votes, {'mirrored, ' if flip < 0 else ''}rotation {rot:.0f} deg, scale {scale:.3f}\"/px")
    if k < 6:
        sys.exit("no convincing match; check the centre, the scale and the catalogue")
    t = rot * D2R
    R = np.array([[np.cos(t), -np.sin(t)], [np.sin(t), np.cos(t)]])
    pred = (C * [flip, 1]) @ R.T / scale + [tx, ty]
    for tol in (4.0, 2.0, 1.5):
        d = np.hypot(pred[:, None, 0] - P[None, :, 0], pred[:, None, 1] - P[None, :, 1])
        j = d.argmin(1)
        ok = d[np.arange(len(C)), j] < tol
        X = np.c_[C[ok], np.ones(ok.sum())]
        coef, *_ = np.linalg.lstsq(X, P[j[ok], :2], rcond=None)
        pred = np.c_[C, np.ones(len(C))] @ coef
    sky = SkyCoord(cat[ok, 0] * u.deg, cat[ok, 1] * u.deg)
    xy = P[j[ok], :2]
    w = fit_wcs_from_points((xy[:, 0], xy[:, 1]), sky, projection="TAN")
    px = w.world_to_pixel(sky)
    rms = float(np.sqrt(np.mean((px[0] - xy[:, 0]) ** 2 + (px[1] - xy[:, 1]) ** 2)))
    sc = float(np.sqrt(abs(np.linalg.det(w.pixel_scale_matrix))) * 3600)
    print(f"WCS: {ok.sum()} stars, rms {rms:.3f} px ({rms * sc:.2f}\"), scale {sc:.4f}\"/px")
    hdr = w.to_header()
    for key in ("DATE-OBS", "EXPTIME"):
        if key in src_hdr:
            hdr[key] = src_hdr[key]
    hdr["HISTORY"] = f"solve_from_catalog.py: TAN fit to {int(ok.sum())} catalogue stars, rms {rms * sc:.2f} arcsec"
    fits.PrimaryHDU(im, hdr).writeto(args.out, overwrite=True)
    print("wrote", args.out)


if __name__ == "__main__":
    main()
