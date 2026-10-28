# FaintTrail

**Synthetic tracking for faint asteroids and comets — an open-source library, CLI and minimal Qt 6 GUI.**

FaintTrail finds moving objects that are too faint to see in any single frame. It
shift-and-stacks a sequence of short exposures along thousands of trial motions,
keeps what adds up, and — the part that matters most — tells you honestly how
likely each detection is to be noise.

It is meant as an engine: a C++17 library with a stable C API, so that astrometry
applications such as AstroFind or Siril can use it, plus a command-line tool and a
small GUI for working with it directly.

> Status: **0.9.0, CPU reference backend and a CUDA backend (bit-identical by design, not yet run on a GPU) behind a pluggable backend interface; validated on synthetic data and a real-data regression suite.**
> It has not yet been checked end to end against real observations with known
> results. Do not submit astrometry to the Minor Planet Center based on it without
> verifying each object independently.

![FaintTrail GUI](docs/screenshot.png)

## Why

Amateur synthetic-tracking work currently depends on closed tools. FaintTrail aims
to be the open, inspectable alternative, with false-alarm control built in rather
than left to the observer:

* **Calibrated threshold.** The search is repeated on copies of the data with the
  frame times shuffled, which removes every real moving object but keeps the noise,
  cosmic rays and star residuals. The threshold is set from the noise peaks of these
  null searches for a chosen false-alarm rate (default 0.02 per search).
* **Per-candidate tests.** Each detection must appear in both halves of the
  sequence, must not depend on a single frame, and must not be a sidelobe of a
  brighter object.
* **Measured completeness.** An injection-recovery test works on any data set,
  real or simulated, and reports the recovered fraction against S/N.
* **Uncertainties.** Positions and motions come with 1-σ errors.

## Download

Binaries for Windows (CLI and GUI, with the CUDA backend), Linux (AppImage GUI and
command-line tool) and macOS (Apple silicon) are attached to each
[release](../../releases). See
[BUILDING.md](BUILDING.md#releases) for what each file contains.

## Building

See [BUILDING.md](BUILDING.md) for Windows (Visual Studio 2022, CUDA, Qt), Linux
and macOS, the build options, and the optional run-time helpers (ASTAP for plate
solving, curl for the Gaia photometry). In short, on Linux:

```sh
cmake -S . -B build -DCMAKE_BUILD_TYPE=Release
cmake --build build --config Release -j
ctest --test-dir build -C Release --output-on-failure
```

and on Windows with the GPU backend and the GUI:

```bat
cmake -S . -B build -G "Visual Studio 17 2022" -A x64 -T cuda=12.8 -DFAINTTRAIL_CUDA=ON -DCMAKE_CUDA_ARCHITECTURES=120 -DCMAKE_PREFIX_PATH=C:\Qt\6.11.2\msvc2022_64
cmake --build build --config Release -j
cmake --build build --config Release --target deploy_gui
ctest --test-dir build -C Release --output-on-failure
```

## Quick start

```sh
# simulate a sequence, search it and score the result against the truth
./build/fainttrail demo

# the same in steps, with files
./build/fainttrail simulate -o sim --seed 7
./build/fainttrail search sim --rate-min 0.2 --rate-max 1.1 -o candidates.csv --tracklets tracklets.csv
./build/fainttrail validate --truth sim/truth.csv --candidates candidates.csv

# how deep does the search go on this data?
./build/fainttrail completeness sim --rate-min 0.2 --rate-max 1.1

# real data: borrow the plate solution of an ASTAP-solved stack or frame
./build/fainttrail search frames/ --wcs-from solved_stack.fits --rate-min 0.1 --rate-max 2

# stack along a known motion (e.g. from an ephemeris)
./build/fainttrail stack sim --rate 0.62 --pa 282 -o stack.fits

# GUI
./build/fainttrail-gui            # File > Open frames, or File > Simulate sequence
./build/fainttrail-gui --demo     # simulate and search at start-up
```

## Input

* **Calibrated frames**, one FITS image each, registered or not. Since 0.9
  FaintTrail checks whether the stars stay on the same pixels; if they move (a
  drifting mount, raw frames), it aligns the frames on the middle one: stars are
  matched, a shift, rotation and scale (meridian flips included) is fitted per
  frame, and the frame is resampled. The offset back to the detector is kept, so
  hot pixels and other detector defects are still recognised. `--align` chooses:
  `auto` (default), `always`, `check` (refuse unregistered frames) or `off`.
  `fainttrail register FRAMES -o DIR` aligns and writes the frames without
  searching. Registration needs at least about 6 stars per frame; frames without
  enough stars (clouds) are dropped and reported.
* **Time** from `MJD-MID`, `MJD-OBS`, `DATE-OBS` (+`EXPTIME`) or `JD`.
* **Plate solution** (optional; without one, positions are in pixels): taken from
  the frames when they have one; borrowed with `--wcs-from FILE` (GUI: *Import
  plate solution*) from any solved image of the same field; or solved on request
  with [ASTAP](https://www.hnsky.org/astap.htm): the **Plate solve (ASTAP)** button
  in the GUI toolbar (Ctrl+P) or `--solve` on the command line. FaintTrail stacks
  the aligned frames and lets ASTAP solve the stack, using the pointing
  (`RA`/`DEC` or `OBJCTRA`/`OBJCTDEC`) and field size (`FOCALLEN`, `XPIXSZ`) of the
  FITS header as hints (the GUI shows them for editing; `--ra`/`--dec`/`--fov` on
  the command line). Solving after a search adds RA/Dec, rates and calibrated
  magnitudes to the candidates without searching again.
* **Optional TAN WCS** (`CD` matrix, `PC`+`CDELT` or `CDELT`+`CROTA2`) for sky
  coordinates, rates in arcsec/min and position angles.
* **Optional** `REGDX`/`REGDY`: the offset from the registered grid back to the
  detector (`detector = registered + offset`). With it, objects that stay fixed on
  the sensor (warm pixels carried along by the registration) are rejected.
* At least 5 frames; synthetic tracking pays off from about 20. The sequence should
  last long enough for targets to move at least two PSF widths.
* **Large sequences** are searched in tiles automatically when they would not fit in
  60 % of the RAM (`--max-memory MB` sets the budget, `--tile N` forces a tile size,
  `--no-tile` turns it off). Frames are then read from disk a tile at a time and the
  result is identical to the untiled search; it takes longer (about 1.4× for the
  Leo Triplet case).

## Output

`candidates.csv` — one row per candidate, FITS 1-based pixels at the reference
epoch (middle of the sequence, written in the header line):

| column | meaning |
|---|---|
| `ok` | 1 if the candidate passed every test |
| `flags` | `half`, `single`, `ghost`, `edge`, `detector` (rejecting); `masked`, `vedge`, `slow` (informational) |
| `x`, `y`, `pos_err_px` | position at the reference epoch and its 1-σ error per axis |
| `vx_pxh`, `vy_pxh`, `vel_err_pxh` | motion in pixels per hour and its error |
| `snr`, `snr_coarse` | refined S/N, and the grid S/N the threshold applies to |
| `snr_half1`, `snr_half2` | S/N in the first and second half of the frames |
| `max_frame_fraction` | largest share of the signal from one frame |
| `false_alarms_per_run` | expected number of noise peaks this strong per search |
| `flux`, `mag_inst` | PSF flux per frame (data units), instrumental magnitude |
| `ra_deg`, `dec_deg`, `rate_arcsec_min`, `pa_deg` | with a WCS |

With a plate solution the terminal table shows RA and Dec (J2000) of every
candidate at the reference epoch, and `--mpc FILE --obscode CODE` (GUI: *File ›
Export MPC observations*) writes MPC 80-column observation lines, three per
passing candidate (first, middle and last frame), with temporary designations
`FT00001`... Paste them into MPChecker or Find_Orb to identify the objects.
Magnitudes are calibrated against Gaia DR3 (since 0.9): FaintTrail downloads the
Gaia stars of the field from VizieR with `curl` (cached per field), measures them
on the stacked sky with the same PSF photometry as the candidates, and fits the
zero point on stars of solar colour. `--band V` (default), `R` or `G`; V and R
come from Gaia G and BP−RP. Fast movers are corrected for trailing. Offline, give
a saved extract with `--catalog FILE`; `--photometry off` skips it. The MPC lines
then carry magnitude and band. Check every object before submitting anything to
the Minor Planet Center.

`tracklets.csv` — predicted position of each passing candidate in every frame
(time, pixel, RA/Dec). This is what a host application turns into ADES records
after its own astrometric and photometric calibration.

## Using the library

C++:

```cpp
#include "fainttrail/fainttrail.hpp"

auto frames = fainttrail::load_frames(paths);          // or fill Frame objects yourself
fainttrail::SearchParams p;
p.vmin_pxh = 8; p.vmax_pxh = 60;
auto prep = fainttrail::prepare(frames, p);
auto res  = fainttrail::search(prep, p, [](double f, const std::string& stage) {
    std::printf("%s %.0f%%\n", stage.c_str(), 100 * f);
    return true;                                       // false cancels
});
for (const auto& c : res.candidates)
    if (c.ok) { /* c.x, c.y, c.vx_pxh, c.snr, c.ra_deg, ... */ }
```

C (`fainttrail/fainttrail.h`), for hosts written in C or wanting a stable ABI:

```c
ft_context* ctx = ft_create();
ft_add_frame(ctx, w, h, pixels, mjd_mid, exptime);    /* per frame, or ft_add_fits() */
ft_set_wcs(ctx, crval, crpix, cd);                    /* optional */
ft_params p; ft_params_default(&p); p.vmax_pxh = 60;
if (ft_search(ctx, &p, progress_cb, user) == FT_OK)
    for (int i = 0; i < ft_candidate_count(ctx); ++i) { ft_candidate c; ft_get_candidate(ctx, i, &c); }
ft_destroy(ctx);
```

Pixel coordinates are 0-based in the C++ API and FITS 1-based in files and the C API.

## How it works

See [docs/ALGORITHM.md](docs/ALGORITHM.md): static-sky subtraction with a leave-out
median, matched filtering, the shift-and-stack grid search, the time-shuffled null
calibration, refinement and the consistency tests.

## Backends

The stacking loop runs behind an interchangeable backend interface with an exact,
testable contract; the cache-blocked CPU (OpenMP) backend is the reference. Since
0.9 there is a CUDA backend (`-DFAINTTRAIL_CUDA=ON`), used automatically when a GPU
is present and designed to give bit-identical results. It has been compiled and its
algorithm tested on the CPU, but not yet run on a GPU: check it with
`fainttrail check-backend --backend cuda` and `ctest`. See
[docs/BACKENDS.md](docs/BACKENDS.md).

## Profiling

`--profile` on `fainttrail search` and `fainttrail regress` prints the time of every
pipeline stage and the stacking throughput; `fainttrail_kernel_bench` measures the
stacking kernel itself. [docs/PROFILING.md](docs/PROFILING.md) has the results for
the full Leo Triplet run: stacking is 88 % of the time (three quarters of it in the
null calibration) and the kernel is limited by memory bandwidth. Since 0.5 the CPU
backend is cache-blocked: the full Leo Triplet run went from 578 s to 175 s on the
same 2-core machine, with bit-identical results. Since 0.6 a search keeps about
3 copies of the sequence in memory at its peak instead of 6 (Leo Triplet: 2.1 GB →
1.1 GB). Since 0.7 sequences that do not fit are searched in tiles, read from disk
a tile at a time, with identical results: the same run in 2 × 2 tiles needs 380 MB.

## Validation

See [docs/VALIDATION.md](docs/VALIDATION.md) for the false-alarm, completeness and
astrometry benchmarks on simulated data, and how to reproduce them with
`fainttrail_benchmark`.

The detection threshold is calibrated with time-shuffled null searches. Since 0.8
artefacts are screened out of the null runs and a curved Gaussian-field tail is
fitted ([docs/CALIBRATION.md](docs/CALIBRATION.md)): the requested false-alarm rate
is now delivered instead of a data-dependent fraction of it, the threshold is
lower (Leo Triplet: S/N 11.2 → 7.1), and a second asteroid appeared in the Leo
Triplet field: (4456) Mawson at V 19.5, confirmed with MPChecker.

Real data are covered by a regression suite ([tests/regression](tests/regression/README.md)):
three asteroids identified in the author's captures, (12631) Mariekebaan,
(21996) 1999 XP31 and (4456) Mawson, checked against their catalogue/ephemeris sky
positions, and a null field
around Ceres. Small 16-bit crops run with every `ctest`; the full frames run
when `FAINTTRAIL_REGRESSION_DATA` is set.

```sh
fainttrail regress tests/regression [--backend NAME] [--strict]
```

## Roadmap

* CUDA backend: verify and tune on the RTX 5080; then decide on GPU preprocessing
* curved motion for fast near-Earth objects
* more real-data cases (fainter and faster objects), then integration with AstroFind and Siril

## Kernel features

[docs/KERNEL_FEATURES.md](docs/KERNEL_FEATURES.md) compares whole-pixel and sub-pixel
shifts, measures what straight-line stacking loses on curved tracks, and tests
per-frame weights, all at equal false-alarm rates. Whole-pixel shifts stay; since
0.8.2 the kernel accepts an optional shift table (curved, ephemeris-guided or
parallax-corrected tracks); weights stay out of the kernel.

## Development

Design decisions, their reasons and the open questions are recorded in
[DEVELOPMENT.md](DEVELOPMENT.md).

## Licence

Copyright © 2026 Jurgen Kobierczynski.

FaintTrail is free software: you can redistribute it and/or modify it under the terms
of the GNU General Public License as published by the Free Software Foundation, either
version 3 of the License, or (at your option) any later version. It is distributed in
the hope that it will be useful, but WITHOUT ANY WARRANTY; without even the implied
warranty of MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See [LICENSE](LICENSE)
for the full text (GPL-3.0-or-later).

Because the library is GPL as well, programs that link it (the CLI, the GUI, or host
applications through the C API) must be distributed under GPL-compatible terms. AstroFind
(AGPL-3.0) and Siril (GPL-3.0) can both use it.
