# Kernel features: shifts, curved tracks, per-frame weights

Three open questions from DEVELOPMENT.md section 9, answered with
`tools/kernel_feature_study.cpp`. Every comparison that involves a threshold runs
the full null calibration with the method under test, so detections are compared
at equal false-alarm rates. The key figure is **S/N / threshold** for the same
injected objects: higher means easier to detect.

## 1. Whole-pixel or sub-pixel shifts

The current kernel rounds every shift to a whole pixel. The alternative samples
each frame with bilinear interpolation at the exact shift; it reads two rows per
frame instead of one and does four times the arithmetic. It could allow a coarser
motion grid. An experimental sub-pixel backend was registered next to the CPU
backend (`register_backend`, new in 0.8.1) and both were run through the same
pipeline: two 384 × 384 sequences of 60 frames per PSF width, 80 injected objects
near the detection limit.

| FWHM | kernel | motions | threshold | S/N / ideal | **S/N / threshold** | time / pass |
|---:|---|---:|---:|---:|---:|---:|
| 1.6 px | whole pixel, grid 0.5 px (current) | 724 | 6.68 | 0.855 | **0.834** | 3.8 s |
| | whole pixel, grid 0.25 px | 2760 | 6.85 | 0.894 | 0.849 (+1.8 %) | 5.2 s |
| | sub-pixel, grid 0.5 px | 724 | 5.90 | 0.766 | 0.847 (+1.6 %) | 10.7 s |
| | sub-pixel, grid 1.0 px | 192 | 5.79 | 0.723 | 0.815 (−2.3 %) | 3.8 s |
| 3.0 px | whole pixel, grid 0.5 px (current) | 724 | 6.77 | 1.005 | **0.682** | 2.3 s |
| | whole pixel, grid 0.25 px | 2760 | 6.90 | 1.035 | 0.689 (+1.0 %) | 4.9 s |
| | sub-pixel, grid 0.5 px | 724 | 6.46 | 0.951 | 0.678 (−0.6 %) | 10.8 s |
| | sub-pixel, grid 1.0 px | 192 | 6.39 | 0.932 | 0.671 (−1.6 %) | 3.7 s |

(The sub-pixel backend is a plain, unblocked implementation; an optimised one
would be about twice as expensive per motion as the whole-pixel kernel, not three
times.)

* Interpolation smooths the signal and the noise alike: the S/N per trial drops,
  the threshold drops with it, and the net effect is within ±2 %.
* Only for undersampled images (FWHM 1.6 px, like the Mariekebaan frames) does
  sub-pixel stacking gain anything (1.6 %), and a finer whole-pixel grid gains the
  same. For well-sampled images it is slightly worse.
* A coarser sub-pixel grid costs about the same as the current kernel and detects
  worse.

**Decision: whole-pixel shifts and the 0.5 px grid stay.** The GPU kernel does not
need texture interpolation. The backend contract is unchanged.

## 2. Curved tracks

The kernel stacks along straight lines. Real tracks curve because of the diurnal
parallax (the observer moves with the Earth's rotation) and, for close objects,
the changing viewing geometry. How far a track departs from the best straight
line over a sequence of T hours (diurnal parallax only, latitude 51°):

| distance | 1 h | 2 h | 4 h | 6 h |
|---:|---:|---:|---:|---:|
| 0.01 AU | 3.2″ | 12.7″ | 51″ | 114″ |
| 0.02 AU | 1.6″ | 6.4″ | 25″ | 57″ |
| 0.05 AU | 0.6″ | 2.5″ | 10″ | 23″ |
| 0.1 AU | 0.3″ | 1.3″ | 5.1″ | 11.4″ |
| 0.2 AU | 0.2″ | 0.6″ | 2.5″ | 5.7″ |
| 0.5 AU | 0.06″ | 0.25″ | 1.0″ | 2.3″ |
| 1 AU (main belt) | 0.03″ | 0.13″ | 0.5″ | 1.1″ |

How much S/N a straight-line search keeps, from tracks injected with a known
curvature into a simulated sequence (FWHM 2.4 px, 60 frames over 59 min), compared
with stacking along the true curved track:

| departure from a straight line | S/N kept |
|---:|---:|
| 0.2 FWHM | 99 % |
| 0.4 FWHM | 95 % |
| 0.6 FWHM | 90 % |
| 0.8 FWHM | 84 % |
| 1.25 FWHM | 73 % |
| 2.1 FWHM | 56 % |

With 1.5–2″/px and a FWHM of 2–3 px (4–5″), the loss stays below 5 % up to a
departure of about 2″: everything in the main belt even over 6 hours, and near-Earth
objects beyond about 0.1 AU for 1–2 hour sequences. Closer objects, or long
sequences on objects within 0.2–0.5 AU, lose a quarter of their signal or more.

**Decision: generalise the kernel's input, not the search.** The cache-blocked CPU
kernel already works from a table of integer shifts per motion and frame; accepting
an explicit table in place of `round(v·t)` costs nothing per addition. With it, the
same kernel can stack along curved tracks (a curvature search for near-Earth
objects), along an ephemeris (known objects), or with parallax-corrected motions.
The contract will take an optional shift table, with `round(v·t)` as the default
generator. A curvature *search* (two more grid dimensions) is not planned now; it
belongs to a later near-Earth mode.

## 3. Per-frame weights

All frames count equally. With changing seeing or transparency, weighting each
frame by its point-source efficiency (signal-to-noise of a star in that frame) is
optimal in theory: the gain is √(n Σs²) / Σs.

**Real sequences.** Frame efficiency measured from 58–135 stars per sequence:

| sequence | frames | efficiency range | variation (CV) | possible S/N gain |
|---|---:|---:|---:|---:|
| Mariekebaan (unguided) | 82 | 0.72–1.23 | 12 % | 0.7 % |
| Leo Triplet | 20 | 0.80–1.09 | 9 % | 0.4 % |
| Ceres field | 11 | 0.84–1.16 | 8 % | 0.3 % |

**Simulated poor conditions**, weights estimated from the stars as a pipeline would
do, full null calibration for both:

| seeing and transparency jitter | weights | threshold | S/N / threshold | detected |
|---:|---|---:|---:|---:|
| 10 % | equal | 6.99 | **0.765** | 3/40 |
| | by efficiency | 7.45 | 0.724 | 2/40 |
| 25 % | equal | 6.77 | **0.819** | 6/40 |
| | by efficiency | 7.15 | 0.801 | 3/40 |
| 40 % | equal | 10.49 | 0.541 | 0/40 |
| | by efficiency | 7.24 | **0.835** | 6/40 |

* In ordinary conditions, weighting does not pay: the theoretical gain is under
  4 %, and the calibrated threshold rises by 5–7 % (fewer effective frames give a
  heavier noise tail), so detection gets slightly worse.
* In very poor conditions, it helps a lot: with 40 % jitter the equal-weight
  threshold rose to 10.5. This is most likely because the worst frames leave
  residuals the screening does not catch; this explanation has not been checked.
  Weighting suppresses those frames.

**Decision: no weights in the kernel.** Weights can be folded into the frames
before stacking (scaling frame i by w_i·√n/√Σw² gives Σwg/√Σw² wherever all frames
cover a pixel), so the kernel and its contract stay as they are. What poor
conditions need is a preprocessing option: down-weight or reject frames whose
efficiency falls far below the rest, applied only when the frames really differ,
and validated like the calibration. It is recorded as an open item.

## Summary for the GPU kernel

| feature | decision |
|---|---|
| shifts | whole pixel, unchanged |
| motion grid | 0.5 px tolerance, unchanged |
| track shape | optional explicit shift table per motion and frame (free in the kernel) |
| frame weights | none in the kernel; possible preprocessing option for poor nights |

## Reproducing

```sh
cmake -S . -B build -DFAINTTRAIL_BUILD_TOOLS=ON && cmake --build build -j
./build/fainttrail_kernel_feature_study shifts
./build/fainttrail_kernel_feature_study curvature
./build/fainttrail_kernel_feature_study weights DIR...
./build/fainttrail_kernel_feature_study weights-sim
```
