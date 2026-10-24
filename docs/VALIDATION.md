# Validation on synthetic data

All numbers below come from `fainttrail_benchmark` (version 0.1.0) on simulated
sequences. Simulated data are not real data: the simulator models a Gaussian PSF,
photon and read noise, sky gradients, seeing and transparency jitter, hot pixels
and cosmic rays, but not every defect of a real camera, mount or sky. Validation
on real sequences with known asteroids is the next step.

![benchmark results](validation.png)

Common settings: 60 frames of 45 s at a 60 s cadence (59 min span), 1.5"/px,
PSF FWHM 3 px, sky 40 e⁻/px/s, read noise 3 e⁻, motion range 7–44 px/h
(0.18–1.1"/min), default search parameters.

## False alarms

40 pure-noise sequences of 256 × 256 px (stars, hot pixels and cosmic rays, no
moving objects) were searched at three requested false-alarm rates:

| requested per search | passing candidates | observed per search | mean threshold S/N |
|---:|---:|---:|---:|
| 0.02 (default) | 0 in 40 runs | < 0.025 | 7.02 |
| 0.2 | 10 in 40 runs | 0.25 | 6.42 |
| 1.0 | 56 in 40 runs | 1.40 | 5.99 |

At the default the calibration holds: about 0.8 false detections were expected
and none occurred. At loose settings the observed rate is 25–40 % above the
requested one. The tail fit is slightly optimistic there, so treat a looser
setting as a lower bound on the real false-alarm rate.

An earlier version reused the real-time preprocessing for the time-shuffled null
searches. That produced 12 false detections in 30 runs at a requested 0.02 per
search: the leave-out static sky correlates the noise of frames that are close in
time, and the shuffled null did not see that correlation. Rerunning the complete
time-dependent pipeline for every null search fixed it.

## Completeness

Injection-recovery on one 512 × 512 px sequence: 25 objects per level at random
positions and motions; the threshold was calibrated once (S/N 7.61) and kept fixed.

| injected stack S/N | recovered | measured S/N |
|---:|---:|---:|
| 5 | 0 / 25 | – |
| 6 | 3 / 25 | 7.96 |
| 7 | 7 / 25 | 8.17 |
| 8 | 11 / 25 | 8.19 |
| 9 | 15 / 25 | 8.68 |
| 10 | 19 / 25 | 9.44 |
| 12 | 25 / 25 | 10.83 |
| 15 | 25 / 25 | 12.59 |
| 20 | 24 / 25 | 18.14 |

50 % completeness is reached at an injected S/N of 8.4. The "injected S/N" is the
ideal background-limited value; the measured S/N is about 10–15 % lower because the
static-sky model absorbs part of the signal of slow objects and the integer-shift
grid loses a little. At S/N 20 one object was detected but flagged.

## Astrometry

76 matched movers (S/N 9–40) in 8 sequences of 384 × 384 px. The reported 1-σ
uncertainties were compared with the true errors:

| quantity | rms (true error / reported uncertainty) |
|---|---:|
| position at the reference epoch | 1.00 |
| motion | 0.97 |

The ideal matched-filter precision, `sqrt(2) σ_psf / (S/N)`, underestimated the
scatter by 25 %; the reported uncertainties include an empirical factor of 1.25 for
that. The factor was determined on these same simulations, so the 1.00 above is
in-sample; it needs confirming on independent and on real data. No systematic
offset was found (mean errors below 0.02 px).

## First real data: (12631) Mariekebaan, 2013-08-01

82 unguided 30 s frames (ATIK 314L, 21:06–21:50 UTC, 44 min) of the field of the
main-belt asteroid (12631) Mariekebaan, registered with a similarity transform,
binned 2 × 2, no darks or flats. A further 17 frames (178–194) in the capture were
byte-identical copies of frame 177 and were removed; duplicated frames would
otherwise look like a stationary chunk of the sequence.

* **One passing candidate**, S/N 10.4, present in both halves (6.9 / 7.8) and in no
  single frame dominantly (largest share 3 %). Motion 17.0 px/h (unbinned) at 31.4°
  from the image x axis.
* **Identified as (12631) Mariekebaan.** The orbit was integrated back from the JPL
  osculating elements (epoch 2026) to the observing time with the Sun, planets and
  Moon as perturbers. Predicted motion: 25.9"/h at 30.0° from the RA direction
  (dRA·cosδ +22.5"/h, dDec −13.0"/h), position 15h39m54s −6°37′. The candidate moves
  at 31.4° from the image x axis, which is the RA direction for this mount: the
  directions agree within 1.5°. The rate implies an image scale of 1.52"/px,
  consistent with the 750 mm f/5 telescope used with a focal-reducing corrector.
  The observer confirmed that the candidate is the asteroid he identified in these
  frames in 2013.
* **Artefacts found and fixed.** The first runs passed (a) tracks along the moving
  blank border of the registered frames and (b) a warm pixel, which the
  registration turns into an object that moves with the mount's drift. FaintTrail
  now grows blank borders, requires enough frames for the static sky, flags `edge`
  tracks, and, when frames carry their registration offsets (`REGDX`/`REGDY`), tests
  every candidate against the hypothesis "fixed on the detector" (`detector` flag).
  After these fixes the warm pixel and all edge tracks are rejected.

## Second real identification: (21996) 1999 XP31, Leo Triplet, 2016-11-30

20 guided 240 s frames (ZWO ASI1600MM, 03:03–04:41 UTC, 1 h 34 min) of the Leo
Triplet, no darks or flats, registered and binned 2 × 2 (2328 × 1760 px). The search
covered 7084 trial motions; the null-calibrated threshold was S/N 11.0.

* **One passing candidate**, S/N 15.5 (halves 10.0 / 12.0, largest single-frame share
  8 %), invisible in the individual frames and in the ordinary stack. The two flagged
  detections were edge artefacts (`half,single,edge`).
* **Astrometry.** The observer's ASTAP plate solution of the stack (0.852"/px) was tied
  to the search frames with 165 stars (0.16 px rms). Candidate position at
  2016-11-30 03:52:23 UTC: RA 11h18m26.56s, Dec +13°45′03.6″ (±0.2″), motion
  +33.9"/h in RA, −0.8"/h in Dec; about magnitude 19 from ASTAP's zero point.
* **Built-in transfer.** The same result comes out of a single command,
  `fainttrail search bin2 --wcs-from stack_solved.fits`: 142 stars matched between the
  full-resolution stack and the 2 × 2-binned frames (0.25 px rms, scale 0.5000), and
  the reported position differs by 0.3″ from the hand-made mapping.
* **MPChecker** (geocentric, 15′ radius, V < 21) returns (21996) 1999 XP31 at
  11h18m26.6s +13°45′01″, V 18.6, motion +34"/h, 0"/h: offset below 3″ (within the
  geocentric-versus-topocentric parallax) and the motion agrees within 1"/h.
* The only other known object in the region, (182337) 2001 PV54 at V 20.6, is about
  four times fainter, which puts it near S/N 4 here, below the threshold; it was
  correctly not reported.

A third real set (Ceres field, 2014-05-16, 11 frames over 5.3 min) is too short
for synthetic tracking: nothing in it moves two PSF widths. It gave no false
detections, and injection-recovery on its real noise reached 50 % completeness at
S/N 10.6 for fast motions (200–600 px/h). It also showed that objects moving only
2–7 PSF widths over a short sequence lose much of their signal to the static-sky
model; that limitation is still open.

## Reproducing

```sh
cmake -S . -B build -DFAINTTRAIL_BUILD_TOOLS=ON && cmake --build build -j
./build/fainttrail_benchmark bench          # about 15 minutes on 2 cores
python3 tools/plot_benchmark.py bench docs/validation.png
```

`--quick` runs a reduced version in a minute or two. The raw CSV files of the run
above are in [`benchmark/`](benchmark/).
