# Threshold calibration

How should FaintTrail decide which detections are real? This is the study behind
the calibration of 0.8, and what it means for the stacking kernel.

## The question

The threshold comes from **null runs**: the time-dependent pipeline is rerun with
the frame times shuffled, which destroys every coherent motion but keeps the noise,
the star residuals and the defects; the peaks of the null detection maps show how
strong noise gets, and the threshold is set where the expected number of noise
peaks per search equals the requested rate (default 0.02, one false detection in
fifty searches).

Three null runs cost three times the real stacking, 75 % of all stacking. Before
fixing the kernel's features we wanted to know whether that can be cheaper: null
runs on part of the field, or on part of the motion grid, or fewer of them. The
study answered a different question first: the calibration itself was wrong.

## Method

* **Held-out null runs.** For three data sets (a synthetic 512 × 512 field with 60
  frames, the full Mariekebaan sequence, the full Ceres field) 96 independent null
  runs were computed and their peaks stored (`tools/calibration_study.cpp`). Each
  strategy was then applied to a few of them, many times over, and the false-alarm
  rate its threshold actually gives was measured on all the other runs
  (`tools/calibration_analysis.py`). The mean of that rate is what a user gets on
  average; its ratio to the requested rate is reported as FA/req.
* **Artefact screening.** For the crops and the full Leo Triplet field the strongest
  null peaks were refined and tested exactly like candidates
  (`tools/null_flags_study.cpp`).
* **End to end.** 200 simulated sequences without moving objects were searched with
  the whole pipeline and the passing candidates counted
  (`tools/false_alarm_study.cpp`).

## Finding 1: the straight-line tail fit was wrong, by a data-dependent amount

Up to 0.7 a straight line was fitted to log N(>S) for 0.67 to 30 null peaks per run
and extrapolated to 0.02 per run. The tail of the noise peaks is curved, so the
extrapolation depends on how far it has to reach, and it uses only the top few
dozen peaks per run, so it scatters. Held-out results for three null runs:

| data | model | mean threshold | spread (sd) | FA/req |
|---|---|---:|---:|---:|
| synthetic 512² | straight line (≤ 0.7) | 7.26 | 0.20 | **0.04** |
| | Gaussian-field tail (0.8) | 6.86 | 0.06 | **1.03** |
| Mariekebaan | straight line | 7.22 | 0.23 | 0.70 |
| | Gaussian-field tail | 6.98 | 0.09 | 1.26 |
| Ceres field | straight line | 7.20 | 0.21 | 0.54 |
| | Gaussian-field tail | 7.05 | 0.11 | 0.72 |

The straight line delivered between 4 % and 70 % of the requested rate depending on
the data, so the threshold was too high by up to 0.4 σ, by an amount nobody could
predict. Wider fit windows, a peaks-over-threshold estimator and a quadratic were
also tried; only the curved models were calibrated, and the Gaussian-field form
was the most stable.

**The 0.8 model.** By the Euler-characteristic heuristic, the expected number of
local maxima above u of a smooth Gaussian random field falls like u^k exp(−u²/2σ²).
The null detection map is the maximum over a family of correlated Gaussian fields
(position and motion), so

    log N(>u) = a + k log u − u² / (2 σ²)

with a, k and σ fitted by least squares to the top 3000 null peaks per run. It uses
a hundred times more peaks than the straight line, follows the curvature, and gives
a threshold with a third to a quarter of the scatter.

## Finding 2: artefacts in the null runs inflated the threshold on real data

On the Leo Triplet data the new model gave S/N 15–16, above the straight line's
already high 11.2, and would have missed 1999 XP31 (S/N 15.5). The reason: the
strongest null peaks there are not noise at all. The time shuffle keeps defects
that sit still on the frames, so the null maps contain the same edge artefacts and
single-frame defects (no darks were used) that the candidate tests reject in the
real search:

| full Leo field, top 3000 null peaks per run | rejected by the tests |
|---|---:|
| S/N 4–5 | 29 % |
| S/N 5–6 | 44 % |
| S/N 6–7 | 95 % |
| S/N above 7 | **100 %** |

Flags: single-frame defect 2931, edge 465, detector 24 (plus statistical flags).
Counting these as noise asked the threshold to rise above artefacts that can never
pass. **0.8 screens the top 3000 null peaks of every run with the artefact tests
(single, edge, detector) at their grid position and motion and leaves artefacts
out of the fit.** The statistical tests (half, ghost) are not applied: they also
reject genuine noise peaks, in an S/N-dependent way that would bend the fitted tail.
The screening costs 0.3–1 s per null run (Leo Triplet 1.0 s for all three runs,
Mariekebaan with 82 frames 2.8 s). Screening with fully refined
candidates gave the same threshold (Leo: 7.11 against 7.10) at many times the cost.

| case | threshold 0.7 | threshold 0.8 | passing candidates |
|---|---:|---:|---|
| Leo Triplet, full field | 11.15 | **7.10** | 1999 XP31 (S/N 15.5) and **(4456) Mawson (S/N 8.4)** |
| Leo Triplet crop | 11.18 | 6.40 | 1999 XP31 |
| Mariekebaan, full | 7.44 | 7.02 | Mariekebaan |
| Mariekebaan crop | 8.11 | 7.60 | Mariekebaan |
| Ceres field, full | 6.90 | 6.86 | none |
| Ceres crop | 6.81 | 6.42 | none |

The new Leo candidate (11h19m33.5s +13°19′21″ at 2016-11-30 03:52:23 UTC, 26″/h at
position angle 113°, halves S/N 6.1 and 5.8, about 0.6 mag fainter than 1999 XP31)
lies 26′ from 1999 XP31, outside the 15′ MPChecker query made for that object. A new
query identified it as **(4456) Mawson**: MPChecker gives 11h19m35.0s +13°19′14″ at
2016-11-30.20 UT, V 19.5, motion +25″/h, −10″/h; the candidate's track at that time
is 2.4″ away (geocentric against topocentric) and moves at +24″/h, −10″/h. With the
0.7 calibration this main-belt asteroid was below the threshold.

## Finding 3: subsets do not make the calibration cheaper

What limits the accuracy is the number of independent noise peaks in the fit. That
grows with field area × number of runs, so:

* **Part of the field** (random calibration tiles) is statistically the same as
  fewer runs over the whole field: 3 runs on a quarter of the field scatter like one
  full run. It saves stacking but not accuracy per unit of work.
* **Part of the motion grid** changes the shape of the tail (fewer, less correlated
  trials), so the scaled counts are biased: 6 runs on a quarter of the grid
  delivered 13–23 % of the requested rate.
* **Fewer full runs** are the honest trade: with the new model two runs scatter by
  0.07–0.14 σ instead of 0.06–0.11 σ for three.

The default stays at three full null runs; `--null-runs 2` (or
`SearchParams::null_runs = 2`) is a reasonable choice when speed matters more.

## End-to-end check

200 simulated sequences of 256 × 256 px with stars, hot pixels and cosmic rays but
no moving objects, searched with the complete pipeline:

| requested per search | 0.7: passing per search | 0.8: passing per search | 0.8 threshold |
|---:|---:|---:|---:|
| 0.02 (default) | 0.005 (1 in 200) | **0.020** (4 in 200) | 6.65 (0.7: 7.03) |
| 0.1 | 0.060 | 0.160 | 6.36 |
| 0.2 | 0.135 | 0.225 | 6.23 |
| 1.0 | 0.68 | 1.02 | 5.90 |

The 0.8 calibration delivers the requested rate; the 0.7 one delivered a quarter
of it at the default. At 0.1 per search the 32 observed false alarms are 1.6 times
the 20 requested, so the fit is slightly optimistic at that intermediate level;
at the default and at 0.2 and 1.0 it agrees within the counting noise.

On the simulated completeness benchmark (512² field, 25 objects per level) the
threshold fell from 7.61 to 6.94 and the 50 % completeness point from S/N 9.2 to
about 7.4.

## What this means for the kernel

* No motion-subset or area-subset input is needed: null runs stay full passes over
  the whole field and the whole motion grid.
* Batching the real pass and the null passes in one kernel call would save nothing
  in memory traffic: every null run has its own filtered frames (the static sky
  depends on the time order). The passes stay separate, and the null passes are
  75 % of the stacking, so the GPU must accelerate all four equally.
* The null passes need the detection map and the motion index at the peaks (the
  artefact screen evaluates the track of each peak), as the real pass does.
* The screening and the tail fit are cheap and stay on the CPU.

## Reproducing

```sh
cmake -S . -B build -DFAINTTRAIL_BUILD_TOOLS=ON && cmake --build build -j
./build/fainttrail_calibration_study --sim --perms 96 --name sim512 --dump sim512.csv
python3 tools/calibration_analysis.py sim512.csv
./build/fainttrail_null_flags_study --data DIR [--wcs-from F] --top 3000
./build/fainttrail_false_alarm_study --runs 200
```

`--null-model linear` restores the 0.7 calibration for comparison.
