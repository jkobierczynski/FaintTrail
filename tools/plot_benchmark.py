#!/usr/bin/env python3
# SPDX-License-Identifier: GPL-3.0-or-later
# Copyright 2026 Jurgen Kobierczynski
"""Plot the output of fainttrail_benchmark.

usage: plot_benchmark.py BENCH_DIR OUTPUT.png
Requires numpy and matplotlib.
"""
import csv
import sys

import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np

BLUE = "#2a78d6"
ORANGE = "#eb6834"
INK = "#2b2b29"
MUTED = "#8a897f"
GRID = "#e4e3dd"


def read(path):
    with open(path) as f:
        rows = [r for r in csv.reader(line for line in f if not line.startswith("#"))]
    head, body = rows[0], rows[1:]
    return {h: np.array([float(r[i]) for r in body]) for i, h in enumerate(head)}


def style(ax, title, xlabel, ylabel):
    ax.set_title(title, loc="left", fontsize=11, color=INK, pad=10)
    ax.set_xlabel(xlabel, color=MUTED)
    ax.set_ylabel(ylabel, color=MUTED)
    ax.grid(True, color=GRID, linewidth=0.8)
    ax.set_axisbelow(True)
    for s in ("top", "right"):
        ax.spines[s].set_visible(False)
    for s in ("left", "bottom"):
        ax.spines[s].set_color(GRID)
    ax.tick_params(colors=MUTED)


def main():
    bench, out = sys.argv[1], sys.argv[2]
    comp = read(f"{bench}/completeness.csv")
    fa = read(f"{bench}/false_alarms.csv")
    ast = read(f"{bench}/astrometry.csv")

    fig, axes = plt.subplots(1, 3, figsize=(13, 3.9))

    ax = axes[0]
    ax.plot(comp["target_snr"], 100 * comp["completeness"], color=BLUE, linewidth=2, marker="o", markersize=6)
    ax.set_ylim(-3, 103)
    style(ax, "Recovered fraction of injected objects", "injected stack S/N", "recovered (%)")

    ax = axes[1]
    lim = [0.005, 2.0]
    ax.plot(lim, lim, color=MUTED, linewidth=1, linestyle="--")
    ax.text(0.012, 0.018, "observed = requested", color=MUTED, fontsize=8, rotation=33)
    # a zero count is drawn as an upper limit (one count) with a downward marker
    for t, o, n, k in zip(fa["target_false_per_run"], fa["observed_per_run"], fa["runs"], fa["passing_candidates"]):
        if k > 0:
            ax.plot([t], [o], color=BLUE, linewidth=0, marker="o", markersize=8)
            label = f"{int(k)} in {int(n)}"
        else:
            o = 1.0 / n
            ax.plot([t], [o], color=BLUE, linewidth=0, marker="v", markersize=9, markerfacecolor="white")
            label = f"0 in {int(n)} (upper limit)"
        ax.annotate(label, (t, o), textcoords="offset points", xytext=(8, -4), fontsize=8, color=INK)
    ax.set_xscale("log")
    ax.set_yscale("log")
    ax.set_xlim(lim)
    ax.set_ylim(lim)
    style(ax, "False alarms in pure-noise sequences", "requested per search", "observed per search")

    ax = axes[2]
    norm = np.hypot(ast["dx_px"], ast["dy_px"]) / (np.sqrt(2) * ast["pos_err_px"])
    ax.scatter(ast["snr"], norm, s=22, color=BLUE, edgecolor="white", linewidth=0.6)
    ax.axhline(1.0, color=MUTED, linewidth=1, linestyle="--")
    ax.set_ylim(0, max(3.5, norm.max() * 1.1))
    style(ax, "Position error / reported uncertainty", "measured S/N", "normalised error")

    fig.tight_layout(w_pad=3)
    fig.savefig(out, dpi=140)
    print(f"wrote {out}")


if __name__ == "__main__":
    main()
