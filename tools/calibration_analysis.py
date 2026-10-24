#!/usr/bin/env python3
# SPDX-License-Identifier: GPL-3.0-or-later
# Copyright 2026 Jurgen Kobierczynski
"""Evaluate threshold-calibration strategies on dumped null runs.

Input: peak dumps written by `fainttrail_calibration_study --dump FILE` (many
independent null runs of one dataset). For every strategy the script repeatedly

  1. takes the null runs the strategy would compute (k runs, possibly only part
     of the field or of the motion grid),
  2. fits the tail and sets the threshold for the requested false-alarm rate,
  3. measures the false-alarm rate that threshold actually gives on all the
     *other* (held-out) null runs, which are full-field, full-grid searches.

The mean of (3) over many repetitions is the false-alarm rate a user gets on
average; the mean threshold is the price paid in sensitivity.

    python3 tools/calibration_analysis.py sim512.csv mariekebaan.csv ... [--csv out.csv]
"""
import argparse
import math
import random
import sys

import numpy as np


def load(path):
    meta = {}
    runs = {}
    with open(path) as f:
        for line in f:
            if line.startswith("#"):
                for tok in line[1:].split():
                    k, v = tok.split("=")
                    meta[k] = v
                continue
            if line.startswith("run,"):
                continue
            r, kind, x, y, s = line.strip().split(",")
            r = int(r)
            d = runs.setdefault(r, {"f": [], "q": [], "h": [], "T": 0, "S": 0, "Q": 0, "H": 0})
            if kind in "TSQH":
                d[kind] = float(s)
            elif kind == "f":
                d["f"].append((int(x), int(y), float(s)))
            else:
                d[kind].append(float(s))
    out = []
    for r in sorted(runs):
        d = runs[r]
        f = np.array(d["f"], dtype=float).reshape(-1, 3)
        out.append({"x": f[:, 0], "y": f[:, 1], "s": f[:, 2], "q": np.array(d["q"]), "h": np.array(d["h"]),
                    "prep": d["T"], "stack": d["S"], "stack_q": d["Q"], "stack_h": d["H"]})
    return meta, out


# ---------------------------------------------------------------- tail fits
def fit_ls(pooled, runs, lo=0.67, hi=30.0):
    """FaintTrail's fit: least squares of log N(>s) per run against s, for the
    sorted peaks with lo <= N per run <= hi."""
    s = np.sort(pooled)[::-1]
    n = (np.arange(len(s)) + 1) / runs
    sel = (n >= lo) & (n <= hi)
    if sel.sum() < 8:
        return None
    x, y = s[sel], np.log(n[sel])
    slope = np.polyfit(x, y, 1)[0]
    if slope >= 0:
        return None
    b = -slope
    a = y.mean() + b * x.mean()
    return a, b


def fit_pot(pooled, runs, n_u=30.0):
    """Peaks over threshold, exponential tail, maximum likelihood: u is the level
    exceeded n_u times per run; b = 1 / mean excess over u."""
    s = np.sort(pooled)[::-1]
    m = int(round(n_u * runs))
    if m < 8 or m >= len(s):
        return None
    u = s[m]
    exc = s[:m] - u
    if exc.mean() <= 0:
        return None
    b = 1.0 / exc.mean()
    a = math.log(m / runs) + b * u
    return a, b


def above(sorted_values, t):
    return len(sorted_values) - int(np.searchsorted(sorted_values, t, side="right"))


def fit_gauss(pooled, runs, lo=0.67, hi=300.0, k=None):
    """Tail of the peaks of a smooth Gaussian random field: the expected number of
    local maxima above u falls like u^k exp(-u^2 / 2 s^2) (Euler-characteristic
    heuristic; k = 3 for a field over position and motion, 4 dimensions). Fit
    log N(>u) = a + k log u - u^2 / (2 s^2) to the sorted peaks with lo <= N/run <= hi;
    k fixed (default 3) or free."""
    srt = np.sort(pooled)[::-1]
    n = (np.arange(len(srt)) + 1) / runs
    sel = (n >= lo) & (n <= hi)
    if sel.sum() < 8:
        return None
    u, y = srt[sel], np.log(n[sel])
    if k is None:
        A = np.column_stack([np.ones_like(u), np.log(u), -0.5 * u * u])
    else:
        A = np.column_stack([np.ones_like(u), -0.5 * u * u])
        y = y - k * np.log(u)
    coef, *_ = np.linalg.lstsq(A, y, rcond=None)
    if k is None:
        a, kk, inv_s2 = coef
    else:
        (a, inv_s2), kk = coef, k
    if inv_s2 <= 0:
        return None
    return ("gauss", a, kk, inv_s2)


def fit_quad(pooled, runs, lo=0.67, hi=300.0):
    srt = np.sort(pooled)[::-1]
    n = (np.arange(len(srt)) + 1) / runs
    sel = (n >= lo) & (n <= hi)
    if sel.sum() < 8:
        return None
    c2, c1, c0 = np.polyfit(srt[sel], np.log(n[sel]), 2)
    if c2 >= 0:
        return None
    return ("quad", c0, c1, c2)


def log_rate(fit, u):
    if fit[0] == "gauss":
        _, a, k, inv_s2 = fit
        return a + k * math.log(u) - 0.5 * inv_s2 * u * u
    _, c0, c1, c2 = fit
    return c0 + c1 * u + c2 * u * u


def threshold(fit, target, floor):
    if isinstance(fit[0], str):
        # solve log_rate(u) = log(target) on the falling branch by bisection
        lo, hi = 3.0, 30.0
        lt = math.log(target)
        if log_rate(fit, lo) < lt:
            return floor
        for _ in range(100):
            mid = 0.5 * (lo + hi)
            if log_rate(fit, mid) > lt:
                lo = mid
            else:
                hi = mid
        return max(hi, floor)
    a, b = fit
    return max((a - math.log(target)) / b, floor)


# ---------------------------------------------------------------- strategies
def select(runs, idx, kind, frac, rng, width, height, grid=8):
    pooled, eff = [], 0.0
    for i in idx:
        r = runs[i]
        if kind == "full":
            pooled.append(r["s"])
            eff += 1
        elif kind == "vel":
            pooled.append(r["q"] if frac == 0.25 else r["h"])
            eff += frac
        else:  # area: random calibration tiles of an 8 x 8 grid
            ntiles = grid * grid
            nsel = max(1, round(frac * ntiles))
            chosen = np.zeros(ntiles, bool)
            chosen[rng.sample(range(ntiles), nsel)] = True
            tx = np.minimum(grid - 1, (r["x"] * grid / width).astype(int))
            ty = np.minimum(grid - 1, (r["y"] * grid / height).astype(int))
            pooled.append(r["s"][chosen[ty * grid + tx]])
            eff += nsel / ntiles
    return np.concatenate(pooled), eff


def cost(runs, kind, k, frac):
    prep = np.mean([r["prep"] for r in runs])
    stack = np.mean([r["stack"] for r in runs])
    if kind == "full":
        st = stack
    elif kind == "vel":
        st = np.mean([r["stack_q"] if frac == 0.25 else r["stack_h"] for r in runs])
    else:
        st = frac * stack
    return k * (prep + st) / (3 * (prep + stack))


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("dumps", nargs="+")
    ap.add_argument("--resamples", type=int, default=600)
    ap.add_argument("--target", type=float, default=0.0)
    ap.add_argument("--csv", default="")
    args = ap.parse_args()

    strategies = [("full", 1, 1.0), ("full", 2, 1.0), ("full", 3, 1.0), ("full", 6, 1.0),
                  ("area", 3, 0.5), ("area", 3, 0.25), ("area", 6, 0.25),
                  ("vel", 3, 0.25), ("vel", 6, 0.25), ("vel", 6, 0.5)]
    fits = {"LS 0.67-30 (current)": lambda p, e: fit_ls(p, e),
            "LS 0.67-300": lambda p, e: fit_ls(p, e, hi=300.0),
            "POT n_u=100": lambda p, e: fit_pot(p, e, 100.0),
            "quadratic 0.67-1000": lambda p, e: fit_quad(p, e, 0.67, 1000.0),
            "Gauss k=3, top 3000": lambda p, e: fit_gauss(p, e, 0.0, 3000.0, 3),
            "Gauss (new), top 3000": lambda p, e: fit_gauss(p, e, 0.0, 3000.0, None)}
    rows = []
    for path in args.dumps:
        meta, runs = load(path)
        name = meta.get("dataset", path)
        target = args.target or float(meta.get("target", 0.02))
        floor = float(meta.get("floor", 5.0))
        W, H = int(meta["width"]), int(meta["height"])
        n = len(runs)
        print(f"\n=== {name}: {n} null runs, {meta.get('frames')} frames {W}x{H}, {meta.get('motions')} motions, "
              f"requested {target}/run")
        # empirical reference: rate at the threshold from all runs, as a sanity line
        allp = np.concatenate([r["s"] for r in runs])
        all_sorted = np.sort(allp)
        run_sorted = [np.sort(r["s"]) for r in runs]
        ref = fit_ls(allp, n)
        print(f"    all {n} runs, current fit: threshold {threshold(ref, target, floor):.3f}, slope b {ref[1]:.2f}")
        srt = all_sorted[::-1]
        emp = []
        for rate in (0.02, 0.05, 0.1, 0.3, 1.0, 3.0, 10.0, 30.0, 100.0):
            j = int(round(rate * n))
            if 1 <= j <= len(srt):
                emp.append(f"{rate:g}/run: {srt[j - 1]:.2f}")
        print("    empirical levels: " + ", ".join(emp))
        print(f"    {'strategy':<14} {'fit':<22} {'cost':>5} {'mean T':>7} {'sd T':>6} {'FA/req':>7}")
        rng = random.Random(1)
        for kind, k, frac in strategies:
            if k >= n:
                continue
            c = cost(runs, kind, k, frac)
            label = f"{kind} {k}" + ("" if kind == "full" else f" x {frac:g}")
            for fname, fit in fits.items():
                Ts, fas = [], []
                for rep in range(args.resamples):
                    idx = rng.sample(range(n), k)
                    pooled, eff = select(runs, idx, kind, frac, rng, W, H)
                    f = fit(pooled, eff)
                    if f is None:
                        continue
                    T = threshold(f, target, floor)
                    cnt = above(all_sorted, T) - sum(above(run_sorted[i], T) for i in idx)
                    Ts.append(T)
                    fas.append(cnt / (n - k))
                if not Ts:
                    continue
                Ts, fas = np.array(Ts), np.array(fas)
                print(f"    {label:<14} {fname:<22} {c:5.2f} {Ts.mean():7.3f} {Ts.std():6.3f} {fas.mean() / target:7.2f}")
                rows.append((name, label, fname, c, Ts.mean(), Ts.std(), fas.mean() / target))
    if args.csv:
        with open(args.csv, "w") as f:
            f.write("dataset,strategy,fit,cost,mean_T,sd_T,fa_ratio\n")
            for r in rows:
                f.write(",".join(str(v) for v in r) + "\n")


if __name__ == "__main__":
    main()
