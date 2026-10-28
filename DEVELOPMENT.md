# FaintTrail development decisions

This file records the decisions that shaped FaintTrail: what was decided, why, what
was rejected, and where the evidence is. New decisions are added at the end of
their section, with the version that introduced them. The detailed studies live in
`docs/`; this file is the index and the reasoning.

Status words: **in force**, **superseded** (with the replacement), **open**.

---

## 1. Scope and product

### D1. Fill the open-source gap for faint-asteroid detection (0.1)
**Decision.** Build an open-source tool for detecting and measuring asteroids and
comets too faint to see in single frames, by synthetic tracking (shift-and-stack
along trial motions).
**Why.** Astrometrica is closed and Windows-only; Tycho Tracker, which does
synthetic tracking well, is closed and commercial. There was no open, scriptable
implementation with calibrated false-alarm control.
**Status.** In force.

### D2. A library with a CLI and a minimal GUI, integrating with existing tools (0.1)
**Decision.** A C++ library is the product; a command-line tool and a minimal Qt 6
GUI sit on top of it. Integration with AstroFind and Siril rather than a new
full-featured astronomy application. A C API (`fainttrail.h`) makes the library
usable from other languages and hosts.
**Why.** The value is in the detection engine. Calibration, stacking, plate solving
and reporting are already covered by tools people use; FaintTrail should plug into
them.
**Rejected.** A complete desktop application of its own.
**Status.** In force.

### D3. Prove it on synthetic data first, then on real data (0.1)
**Decision.** The first version was built and validated against a simulator with
known truth (Gaussian PSF, photon and read noise, sky gradients, seeing and
transparency jitter, hot pixels, cosmic rays), then tested on the author's own
captures.
**Status.** In force (see section 6).

### D4. Name: FaintTrail (0.1.2)
Working name *synthtrack*; renamed to **FaintTrail**.

### D5. Licence: GPL-3.0-or-later (0.1.2)
**Decision.** GPL-3.0-or-later.
**Superseded.** The first builds were Apache-2.0.

---

## 2. Language, build and platforms

### D6. C++17, CMake ≥ 3.20, OpenMP, Qt 6 (0.1)
**Decision.** C++17 for the library and CLI, CMake build, OpenMP for parallelism,
Qt 6 Widgets + Concurrent for the GUI (optional). Portable code by default
(`FAINTTRAIL_NATIVE=OFF`); `-march=native` is opt-in.
**Why.** Fast, dependency-light, and the same language as Siril's ecosystem and a
future CUDA backend. Portable default builds because `-march=native` gained only
about 15 % on the memory-bound kernel (PROFILING.md).
**Status.** In force.

### D7. CI on Linux, macOS and Windows (0.1)
GitHub Actions builds and runs the unit tests and the real-data fixtures on all
three platforms; the GUI is smoke-tested offscreen on Linux.

### D7a. Windows builds with the Visual Studio 2022 toolset (0.9)
**Decision.** Windows builds use the Visual Studio 2022 generator (MSVC v143), also
when Visual Studio 2026 is installed; CUDA via `-T cuda=<version>`, the CUDA
architecture given explicitly, Qt's MSVC 2022 kit, the Qt runtime copied by the
`deploy_gui` target, and the CUDA runtime linked statically. BUILDING.md has the
steps and the errors met on the way.
**Why.** CUDA 12.8 accepts host compilers only up to VS 2022 and installs its
MSBuild integration only into the Visual Studio versions it knows; the first
build on the author's PC (VS 2026) failed on both. Found and fixed while bringing
up the CUDA backend: an import-library name clash (LNK1149), the static/DLL C
runtime mismatch (LNK4098), `/O2` breaking Debug builds (D8016), MSVC's
`INFINITY` in device code, and a GUI without its Qt DLLs.

### D7b. Tagged releases with binaries from CI (0.9)
**Decision.** A tag `vX.Y.Z` (checked against the CMake version) triggers
`.github/workflows/release.yml`: Windows x64 (VS 2022, CUDA 12.8, Qt 6.8 MSVC,
`windeployqt` with the compiler runtime, OpenMP DLL), Linux x86_64 on Ubuntu 22.04
(CUDA 12.8; CLI as tar.gz, GUI as AppImage via linuxdeploy), macOS arm64 (Homebrew
Qt and libomp, `macdeployqt`, ad-hoc signed, not notarised). Every job runs the
tests before packaging; the release carries SHA-256 checksums. CUDA code for
compute capability 7.5–12.0 plus PTX. The packaging lives in `packaging/` scripts
so it can be run and debugged locally.
**Status.** The Linux packaging was run locally (AppImage started and searched the
demo sequence offscreen); the Windows and macOS jobs run for the first time on
GitHub.

### D8. No external astronomy libraries in the core (0.1)
FITS reading and writing, TAN WCS, time conversion and star matching are
implemented in the library. Python (numpy, scipy, astropy) is used only by tools
outside the library (registration, catalogue plate solving, analysis scripts).

---

## 3. Input assumptions

### D9. Input is calibrated, registered frames (0.1)
**Decision.** FaintTrail does not calibrate or align frames. Times come from
`MJD-MID`, `MJD-OBS`, `DATE-OBS` (+`EXPTIME`) or `JD`.
**Why.** Stacking software already does this well; duplicating it adds risk.
`tools/register_sequence.py` exists only to prepare test data reproducibly.
**Status.** Superseded in 0.9: the input is checked (D42) and aligned when needed (D43).

### D42. Refuse unregistered input (0.9)
**Decision.** `prepare()` finds stars in the first, middle and last frame (by time;
extended local maxima, so hot pixels and cosmic rays do not count), matches them by
voting on the offsets between star pairs, and stops with an explanation when the
star field moves by more than 1.5 px. `allow_unregistered` / `--allow-unregistered`
overrides it; the GUI asks. The check takes about 0.03 s.
**Why.** Raw frames from a drifting mount (the author's Crescent Nebula H-alpha
sequence, 33 frames, drifting 54 px/h) gave 4 passing and 129 flagged "candidates",
all of them stars, with a calibrated threshold of S/N 55. Registered, the same
frames gave none, at S/N 7.6. The measured drift, (58.9, 3.2) px, matched the
registration script's (58.6, 3.2) px.
**Status.** The refusal became the `check` mode of D43; the default now aligns.

### D43. Built-in registration (0.9)
**Decision.** FaintTrail aligns frames itself when the check (D42) finds them
unregistered (`align = "auto"`, the default; also `always`, `check`, `off`;
`fainttrail register` writes aligned frames). Reference: the middle frame in time.
Stars as in D42, centroided in a 7 × 7 window; the starting offset by voting on
star-pair offsets, with and without a 180° rotation (meridian flip); then a
similarity transform (shift, rotation, scale; shift only with fewer than 12
stars) fitted to nearest-neighbour matches with 3σ clipping. Resampling with
bicubic convolution (Keys, a = −0.5); pixels whose source window leaves the frame
become NaN, so the existing border handling applies. The offset to the detector is
recorded (`REGDX`/`REGDY` semantics, added to any offset the frame already had),
which keeps the detector flag working for hot pixels in raw frames; the plate
solution of the first solved frame is transformed onto the reference grid. Frames
with too few matched stars are dropped and reported. The GUI aligns at load time,
so the viewer and the candidates share one pixel grid.
**Why.** The author's sequences come straight from the camera; asking for an
external registration step made every unguided or drifting sequence produce
star-filled candidate lists. The original reason for D9 (duplicating stacking
software adds risk) is met by testing: the aligner is checked against known shifts
and against the regression pipeline.
**Evidence.** Crescent Nebula, 33 raw frames drifting 54 px/h: aligned in 1.2 s on
2 cores, star residual 0.09 px rms (the Python script: 0.25 px), threshold 7.53
(Python-registered: 7.61), no candidates in either. Simulated drifting sequence
with a mover: S/N 12.18 registered vs 12.48 aligned, threshold 6.43 vs 6.44, no
extra candidates (unit test).
**Not covered.** Tiled searches cannot align (they read a tile at a time): the CLI
checks three frames first and searches untiled when alignment is needed. Frames of
different sizes are aligned by the GUI and the library, but the CLI's tiled frame
reader still requires equal sizes. Field rotation beyond a similarity (distortion)
is not modelled. The small regression crops have too few stars (2–6) to register.

### D10. Registration offsets `REGDX`/`REGDY` (0.1.x)
**Decision.** Frames may carry the offset from the registered grid back to the
detector (`detector = registered + offset`). With it, objects fixed on the sensor
are rejected (`detector` flag).
**Why.** On the Mariekebaan data a warm pixel, carried along by the registration,
looked like a moving object. A simple scatter test failed; a hypothesis test
("does a fixed detector position explain the signal about as well?", ≥ 0.8 of the
track S/N) works.
**Status.** In force.

### D11. Plate solutions can be borrowed (`--wcs-from`) (0.1.2)
**Decision.** A TAN WCS can be transferred from any plate-solved image of the same
field (for example an ASTAP-solved stack): stars are matched by translation voting
over 8 orientations × 5 scales, then a similarity transform is fitted. SIP/PV
distortion terms are ignored (reported).
**Why.** Many amateur sequences are not solved frame by frame, but a stack is easy
to solve.
**Status.** In force. Distortion support is open.

### D12. Catalogue plate solving as a tool, not a library feature (0.3)
`tools/solve_from_catalog.py` solves a stack from an approximate centre and a short
catalogue extract (Gaia DR3) when no blind solver is at hand. Used for the
Mariekebaan sequence (67 stars, 0.41″ rms).
**Status.** Still in force for the catalogue itself; blind solving comes from ASTAP
since 0.9 (D44).

### D44. Plate solving through ASTAP; sky positions and MPC lines (0.9)
**Decision.** When the frames carry no plate solution, FaintTrail writes the static
sky (median of the aligned frames: deep, no cosmic rays, no movers) to a temporary
FITS file and runs ASTAP on it (`astap_cli` preferred, else `astap`; found through
`--astap`, `FAINTTRAIL_ASTAP`, the PATH or the usual install folders), passing the
header's pointing and field size as hints. The TAN solution is read from ASTAP's
`.wcs` or `.ini` output and applies directly to the aligned grid. `--solve auto`
(default) solves only when needed and ASTAP is present; the GUI does the same.
Candidates get RA/Dec (terminal table, GUI table, CSV) and `--mpc` writes MPC
80-column lines (3 per passing candidate, temporary designations, no magnitude).
**Why.** No star catalogue in FaintTrail (D8); ASTAP is free, runs on Windows,
Linux and macOS, solves blind in seconds and is what the author already uses.
Solving the stack instead of single frames is faster and more reliable.
**Tested.** Unit test with a stand-in solver (POSIX), the leo-triplet-cut crop with
its WCS removed and a stand-in returning the known solution: (21996) 1999 XP31 at
11 18 26.56 +13 45 03.8 for 03:52:23 UTC, against MPChecker's geocentric 11 18 26.6
+13 45 01 (2.8″, within the case's 5″ tolerance; topocentric parallax is of that
order). The hand-over to a real ASTAP is still to be checked on the author's PC.

### D45. Photometric calibration against Gaia DR3 (0.9)
**Decision.** With a plate solution, the Gaia DR3 stars of the field (G < 19 by
default) are downloaded from VizieR with curl (no HTTP library in the core, D8;
curl ships with Windows 10+, Linux and macOS) and cached per field; `--catalog`
reads a saved extract instead. Field stars are measured on the static sky with
the candidates' own PSF-weighted flux estimator (same Gaussian template, same
bilinear sampling); blends, saturated cores (above 60 % of the brightest level),
stars below S/N 20 and stars too close to the edge are skipped. The zero point is
the 3σ-clipped median over stars of solar colour (0.5 < BP−RP < 1.3, widened
when fewer than 8). Bands: Gaia G, Johnson V and Cousins R from G and BP−RP
(Riello et al. 2021). Candidate fluxes are corrected for trailing within an
exposure (the PSF estimator recovers (2s√π/L)·erf(L/4s) of a trail of length L).
Magnitudes go to the terminal and GUI tables, the CSV (`mag`, `mag_err`,
`mag_band`) and the MPC lines (columns 66–71).
**Why.** MPChecker and the MPC compare magnitudes in V (or R/G); an instrumental
magnitude cannot be compared. Calibrating on solar-coloured stars matches an
asteroid's colour, so an unfiltered or luminance camera needs no colour term to
first order.
**Tested.** Simulation with the field stars as the catalogue: zero point within
0.06 mag of the truth, star scatter < 0.08 mag; a mover trailing 2.2 px per
exposure measured 18.88 against a true 18.90. Not yet run on real data (VizieR is
not reachable from the development machine); first check: (4456) Mawson against
its predicted V 19.5.

---

## 4. Algorithm

### D13. Static sky by leave-out median (0.1)
**Decision.** The static sky of frame *i* is the per-pixel median over time that
leaves out the frames within `static_exclude` ranks of *i* (default N/8).
**Why.** A plain median subtracts part of a slow mover from its own frames; the
first version lost signal that way.
**Status.** In force. Cost is O(n²) per pixel; a faster form is open (PROFILING.md
item 4).

### D14. Robust per-frame noise, temporal-scatter mask and weights, cosmic-ray rejection (0.1)
MAD noise per frame; pixels with excess temporal scatter are masked or down-weighted;
single-frame defects sharper than the PSF are removed before filtering.
Cosmic-ray rejection was added after cosmic rays produced false candidates.

### D15. Matched filter normalised to unit variance; adaptive-moment FWHM (0.1)
**Why.** A unit-variance S/N image makes stacks of any length directly comparable
to a threshold. The first FWHM estimator underestimated wide PSFs; adaptive
(Gaussian-weighted) moments of isolated stars fixed it.

### D16. Blank registration borders are grown by 3 px (0.1.x)
**Why.** Interpolation damages the pixels next to the blank border of registered
frames; on Mariekebaan these produced edge tracks.

### D17. Integer shifts, linear motion, velocity grid with spacing 4·grid_tol/span (0.1)
**Decision.** Shifts are `round(v·t)` (half away from zero); the grid spacing keeps
the worst-case mismatch at `grid_tol_px` (default 0.5 px) over the sequence.
**Why.** Integer shifts are exact, cheap and make the kernel a contiguous row add.
**Status.** In force; sub-pixel shifts are an **open** kernel decision (section 9).

### D18. Candidate tests (0.1, extended 0.1.x)
Rejecting flags: **half** (both halves of the sequence must show the object),
**single** (no frame may dominate), **ghost** (sidelobes and duplicates of a
stronger track), **edge** (track touches blank borders), **detector** (fixed on the
sensor). Informational: masked, vedge, slow.

### D19. Position uncertainty with an empirical factor 1.25 (0.1)
`pos_err = 1.25 · √2 · σ_psf / (S/N)`. The ideal matched-filter precision
underestimated the scatter on simulated data by about 25 %. The factor needs
confirming on independent and real data.

### D20. Threshold from time-shuffled null searches (0.1)
**Decision.** The search is repeated `null_runs` times (default 3) with the frame
times permuted; the null peaks give the noise tail and the threshold is set for
`target_false_per_run` (default 0.02).
**Why.** Analytic thresholds underestimate real false alarms; the shuffle keeps
every property of the data except coherent motion.
**Detail that matters.** The complete time-dependent pipeline is rerun for each
permutation. Reusing the real-time preprocessing produced 12 false detections in 30
runs at a requested 0.02 (the leave-out static sky correlates neighbouring frames).
**Status.** In force; the tail model was replaced in 0.8 (D35, D36).

### D21. Linear tail fit for the null peaks (0.1)
Straight line through log N(>S) for 0.67–30 peaks per run.
**Superseded** by D35 (0.8): biased by a data-dependent amount.

---

## 5. Engineering principles

### D22. Determinism and exact reproducibility
Results must not depend on the thread count, block sizes, tiling or platform RNG.
FaintTrail uses its own RNG (no `std::shuffle`/distributions), ties in the stacking
go to the lowest motion index, and every optimisation since 0.5 has been required
to give **bit-identical** results to the plain implementation.

### D23. Every optimisation is checked against an independent implementation
`run_pass_naive` (one pixel at a time) and `run_pass_simple` (the 0.2–0.4 kernel)
are kept as references; unit tests compare bit for bit, including awkward block
sizes and tie-breaking.

### D24. Real-data regression suite with independent truth (0.3)
**Decision.** Case files describe what must be found (known objects with positions,
motions and catalogue sky positions), what may not appear (unexplained passing
candidates), the threshold range, and a baseline candidate list. Truth must be
independent of FaintTrail (MPChecker, ephemeris, observer).
**Why.** Synthetic benchmarks do not show real-camera artefacts (edges, warm pixels,
missing darks), as Mariekebaan and Leo proved.
**Status.** In force. Mutation tests confirmed the suite catches planted faults.

### D25. Small 16-bit crops in the repository, full frames outside (0.3)
**Decision.** Fixtures are 16-bit crops (quantisation ≤ ¼ of the noise; identical
candidates to float crops), about 6 MB, run in every `ctest` and CI. Full-frame
cases are skipped unless `FAINTTRAIL_REGRESSION_DATA` is set.
**Permission.** The author agreed to publishing crops of his own captures
(Mariekebaan, Leo Triplet, Ceres).

### D26. Each case documents how its data were made
Every case file contains the exact `register_sequence.py`, `solve_from_catalog.py`
and `fainttrail crop` commands. Full-frame pixel positions refer to that
registration (reference frame = middle frame).

---

## 6. Validation record

* **Synthetic** (VALIDATION.md): false alarms, injection-recovery completeness,
  astrometric accuracy.
* **(12631) Mariekebaan**, 2013-08-01, 82 × 30 s: found (S/N ≈ 10); sky position
  1.1″ from the ephemeris after the Gaia solution.
* **(21996) 1999 XP31**, Leo Triplet, 2016-11-30, 20 × 240 s: found (S/N 15.5);
  3.0″ from MPChecker.
* **(4456) Mawson**, same Leo Triplet sequence: found in 0.8 (S/N 8.4) once the
  calibration stopped counting artefacts; 2.4″ from MPChecker. It was below the
  0.7 threshold.
* **Ceres field**, 11 × 30 s over 5.3 min: too short for synthetic tracking; null
  case (nothing may pass).

---

## 7. Performance architecture (steps 1–4 of the pre-GPU plan)

### D27. Pre-GPU plan (0.1.2)
Before writing a GPU backend: (1) backend interface, (2) real-data regression
suite, (3) profiling, (4) tiling and memory layout, (5) decide the kernel features,
(6) CUDA. Steps 1–4 are done; step 5 is mostly settled (D34–D39), with half
precision, multi-candidate output, core-only tile output and GPU preprocessing
still open.

### D28. Backend interface with an exact contract (0.2)
**Decision.** The stacking loop sits behind `SearchBackend`. The contract fixes
rounding, coverage, summation order (ascending frames, single precision), outputs
(0 / −1 where undefined), tie-breaking, and the tolerance other backends may use
(1e-4·max(1,|S/N|); motion index may differ only on ties within that tolerance).
`check_backend()` and `fainttrail check-backend` test it.
**Why.** A GPU backend must be verifiable pixel by pixel against the CPU.
See BACKENDS.md.

### D29. Profile before optimising (0.4)
`--profile` (time and memory per stage), `fainttrail_kernel_bench`. Finding: the
stacking was 88 % of the Leo run, three quarters of it in the null runs, and the
kernel was limited by memory traffic (about 12 bytes per addition), not arithmetic.
See PROFILING.md.

### D30. Cache-blocked CPU kernel (0.5)
**Decision.** Bands of 8 rows owned by one thread, motions in blocks sized to keep
about 600 kB of accumulators in L2, coverage from the frames' column intervals as
row segments. Bit-identical to the plain kernel.
**Result.** Leo Triplet 578 s → 175 s on 2 cores (stacking 3.7×).

### D31. Memory layout: about three copies of the sequence at peak instead of six (0.6)
**Decision.** Residuals are freed as soon as each frame is filtered and recomputed
exactly where needed (photometry patches; the edge test uses the background-subtracted
frames and the static sky); null runs keep only their filtered frames;
`prepare(std::move(frames))` takes over the loaded buffers; the CLI and GUI tune
glibc's mmap threshold so freed frames return to the system. The GUI and
`fainttrail stack` keep residuals (`keep_residuals`).
**Result.** Leo 2073 → 1148 MB, Mariekebaan 701 → 365 MB, identical output.

### D32. Spatial tiling with exact equivalence (0.7)
**Decision.** `search_tiled()` reads frames from disk a rectangle at a time
(`FrameSource`), gives every tile a stack halo (longest shift) and a filter halo,
measures whole-frame quantities first from the same sample positions, assembles
the detection maps before peak finding, and refines on tiles with a refinement halo.
Results are byte-identical to the untiled search.
**Why.** Exactness keeps the regression baselines and the backend contract valid
whatever the memory budget.
**Cost.** About twice the preprocessing plus the halos (Leo 1.4–1.6×, Mariekebaan
3.2× because of the O(n²) leave-out median).
**Policy.** The CLI tiles only when the untiled search would need more than 60 % of
the RAM (`--max-memory`, `--tile`, `--no-tile`). The GUI and C API stay untiled.

### D33. Rejected: halving precision, merging passes to save memory
Half-precision storage would halve traffic but break the 1e-4 tolerance (open, see
section 9). Dropping the real filtered frames during null runs and rebuilding them
afterwards would save one copy at about 4 % of run time; not done yet.

---

## 8. Threshold calibration (0.8)

### D34. Keep full null runs; no field or motion subsets
**Why.** The accuracy is set by the number of independent noise peaks, which grows
with field area × runs. Calibrating on part of the field equals fewer runs; on part
of the motion grid it biases the tail (13–23 % of the requested rate). Three full
runs stay the default; `--null-runs 2` is an acceptable speed trade.

### D35. Gaussian-field tail model
**Decision.** `log N(>u) = a + k log u − u²/(2σ²)`, fitted to the top 3000 null
peaks per run (`null_model = "gauss"`).
**Why.** The straight line delivered 4–70 % of the requested false-alarm rate,
depending on the data; the curved model delivers it (72–126 % in held-out tests,
0.020 observed at a requested 0.02 in 200 end-to-end searches) with a third to a
quarter of the threshold scatter. `--null-model linear` keeps the old behaviour.

### D36. Screen artefacts out of the null runs
**Decision.** The top 3000 null peaks of each run are checked with the artefact
tests (single, edge, detector) at their grid position and motion; artefacts are
left out of the fit. Statistical tests (half, ghost) are not applied to nulls.
**Why.** The time shuffle keeps defects that sit still on the frames. On the full
Leo field every null peak above S/N 7 was an artefact the real search rejects;
counting them raised the threshold to 11.2 and hid (4456) Mawson at S/N 8.4. Half
and ghost reject genuine noise peaks in an S/N-dependent way and would bend the tail.
Refined screening gave the same threshold at many times the cost.
See CALIBRATION.md.

---

## 9. Kernel features and GPU (steps 5 and 6)

### Fixed by the calibration study
* Null passes remain full passes over the whole field and motion grid; no subset
  inputs.
* Batching the real and null passes in one call saves no memory traffic (each null
  run has its own filtered frames); all four passes must be fast.
* Null passes need the detection map and the motion index at the peaks.
* Screening and the tail fit stay on the CPU.

### Target hardware
The author's PC: Intel Core Ultra 7 265K (20 threads) and an NVIDIA RTX 5080
(16 GB, about 960 GB/s). The binned Leo search fits in video memory several times
over; GPU tiling is needed only for very large or long sequences. CUDA via
VirtualBox is not an option; the GPU backend targets native Linux/Windows.

### D37. Whole-pixel shifts and the 0.5 px grid stay (0.8.1)
**Why.** At equal false-alarm rate, sub-pixel (bilinear) stacking gained 1.6 % on
undersampled images (FWHM 1.6 px) and lost 0.6 % on well-sampled ones (FWHM 3 px),
at about twice the kernel cost; a finer whole-pixel grid gained 1–2 % for 4× the
motions; a coarser sub-pixel grid lost 2 %. Interpolation smooths signal and noise
alike. The GPU kernel needs no texture interpolation. See KERNEL_FEATURES.md.

### D38. An optional explicit shift table in the backend contract (decided 0.8.1, implemented 0.8.2)
**Why.** Straight-line stacking keeps 95 % of the S/N up to a departure of 0.4 FWHM
(about 2″): all main-belt objects, and near-Earth objects beyond about 0.1 AU for
1–2 h sequences. Closer objects or long sequences lose a quarter or more. The
cache-blocked kernel already stacks from a shift table per motion and frame, so
accepting an explicit table costs nothing per addition and enables curved,
ephemeris-guided and parallax-corrected stacks later. `round(v·t)` stays the default
generator; a curvature search is not planned now.
**Implementation.** `ShiftTable` (integer dx, dy per motion and frame) and
`StackPass::shifts`; `linear_shift_table()` builds the default; backends read shifts
through `pass_shift()`. The naive, simple and cache-blocked implementations agree
bit for bit on arbitrary tables (curved, jumping, far outside the frame), the
straight-line table reproduces the default exactly, and `check-backend` includes a
curved table. Given the curvature, the kernel recovers the full S/N of curved
tracks (KERNEL_FEATURES.md).

### D39. No per-frame weights in the kernel (0.8.1)
**Why.** On the three real sequences frame efficiency varies by 8–12 %, worth 0.3–0.7 %
in S/N. In simulated ordinary conditions weighting raised the calibrated threshold
by 5–7 % and detected slightly fewer objects; only with 40 % seeing and transparency
jitter did it help strongly (S/N/threshold 0.54 → 0.84). Weights can be folded into
the frames before stacking, so the contract stays unchanged; a conditional
down-weighting or rejection of poor frames is an open preprocessing item.

### D40. Backends can be registered (0.8.1)
`register_backend(name, factory)` adds experimental or out-of-tree backends (the
sub-pixel study backend; later possibly a GPU backend built as a plug-in). `auto`
never selects a registered backend.

### D41. CUDA backend: one thread per pixel, exact order, bit-identical (0.9)
**Decision.** The CUDA backend sums in the contract's order (frames ascending for
each motion, motions ascending for the maximum), in single precision without fused
multiply-adds, so it is bit-identical to the CPU backend rather than merely within
tolerance. One thread per output pixel, 8 motions per frame sweep in registers,
coverage tested per frame in the kernel, row bands when the frames do not fit in
device memory, launches of about 4·10⁹ additions with the running maximum kept on
the device. Off by default (`-DFAINTTRAIL_CUDA=ON`); `auto` prefers it when a device
is usable and falls back to the CPU otherwise; the runtime is linked statically.
**Why.** Exactness keeps the whole validation (regression baselines with `--strict`,
tiled = untiled) valid on the GPU without new tolerances, and makes any
discrepancy a bug rather than a judgement call. The per-thread order costs nothing
on a GPU: each thread owns its pixel, so no reduction across threads is needed.
**Testing without a GPU.** The per-pixel code is one header shared by the kernel and
a host emulation of the complete algorithm (bands, launches); unit tests compare the
emulation with the naive contract on every machine. The backend source was compiled
with Clang's CUDA mode against CUDA 12.9 headers (no register spills, 44 registers;
PTX assembles for sm_89 and sm_120), and a build with CUDA support was checked to
fall back to the CPU without a driver. CI compiles it with nvcc 12.8.
**Status.** Verified on the RTX 5080 (Windows, CUDA 12.8, sm_120, 2026-10-26):
`check-backend --backend cuda` bit-identical to the reference in 4 of 4 checks
(naive contract, curved shift table, simple reference, forced row bands). Still to
`ctest`: all regression crops identical to their baselines on the GPU. Full Leo
Triplet (22 × 4656 × 3520 px, 3376 motions, 4 passes): stacking 92.4 s on the 265K →
7.4 s on the GPU (661 G adds/s; kernel alone 770 G adds/s), whole search 105 s →
20.4 s. Preprocessing (static sky in the real and null runs) is now the largest
cost (PROFILING.md); kernel tuning waits until that is addressed.

### Still open
* **Half-precision frames** (halves traffic, breaks the current tolerance).
* **Output**: best motion per pixel only, or several candidates per pixel.
* **Core-only output window** for tiles, to avoid stacking the halos.
* **Preprocessing on the GPU**: after a GPU kernel, preprocessing (27 s for Leo,
  more for long sequences) becomes most of the run; the faster leave-out median is
  part of this decision.
* **Baselines on the target PC**: rerun `fainttrail_kernel_bench` and the profile on
  the 265K, and the same runs with `--backend cuda`.

---

## 10. Other open items

* SIP/PV distortion in transferred plate solutions.
* Curved motion (see section 9).
* Slow movers in short sequences lose signal to the static sky.
* The 1.25 astrometric factor needs independent confirmation.
* Faster leave-out median (O(n log n)).
* Conditional down-weighting or rejection of poor frames (D39), with its own validation.
* Tiling in the GUI and the C API.
* More real cases: fainter and faster objects, unbinned frames.
* AstroFind and Siril integration.

---

## Version history

| version | main change |
|---|---|
| 0.1 | library, CLI, Qt GUI; synthetic validation; null calibration |
| 0.1.x | real data: Mariekebaan, 1999 XP31; edge and detector flags; `--wcs-from` |
| 0.1.2 | renamed FaintTrail; GPL-3.0-or-later |
| 0.2 | backend interface and contract |
| 0.3 | real-data regression suite; crop and registration tools; Gaia solution of Mariekebaan |
| 0.4 | profiling (`--profile`, kernel benchmark) |
| 0.5 | cache-blocked CPU kernel |
| 0.6 | memory layout |
| 0.7 | spatial tiling with exact equivalence |
| 0.8 | threshold calibration: Gaussian-field tail, artefact screening; (4456) Mawson found |
| 0.8.1 | kernel feature study (shifts, curved tracks, weights); backend registration |
| 0.8.2 | explicit shift tables in the backend contract |
| 0.9.0 | CUDA backend (verified bit-identical); built-in frame alignment; ASTAP plate solving, RA/Dec and MPC output |
