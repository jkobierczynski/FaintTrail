# Real-data regression suite

Each `.case` file names a sequence of registered frames, the search settings, and
what the search must find there. `fainttrail regress` runs them:

```sh
fainttrail regress tests/regression                     # all cases
fainttrail regress tests/regression --backend cuda      # same cases, another backend
fainttrail regress tests/regression --strict            # baseline changes fail too
fainttrail regress tests/regression --report times.csv  # one line per case, with timings
ctest --test-dir build -L regression                    # as part of the test suite (strict)
```

## What is checked

| check | kind | fails the case |
|---|---|---|
| every known object is found within its position and motion tolerance | truth | always |
| ... passes every test, reaches `min_snr`, and (if a WCS is present) lies within `sky_tol` of its catalogue position | truth | always |
| no more than `max_unexpected` other passing candidates | limit | always |
| threshold inside `[threshold_min, threshold_max]` | limit | always |
| same candidates and verdicts as the baseline; threshold and S/N of passing candidates within 2 % | baseline | with `--strict` |

The truth checks hold the pipeline to independent facts: an asteroid identified
by MPChecker or confirmed by the observer must be found, and nothing else may pass.
The baseline catches everything those checks would let through: a new rejected
candidate, a verdict that flips, a drift of the threshold. After an intentional
change of the algorithm, review the differences and rewrite the baselines with
`--write-baseline`.

## The cases

| case | data | what it covers |
|---|---|---|
| `mariekebaan-cut` | 82 × 128 × 128 px, in the repository | (12631) Mariekebaan at S/N ≈ 11, a slow mover (9 px/h binned) in unguided data, checked against its ephemeris position; WCS from a Gaia DR3 solution |
| `leo-triplet-cut` | 20 × 200 × 160 px, in the repository, with WCS | (21996) 1999 XP31 at S/N ≈ 15, checked against its MPChecker position; two edge artefacts on the left border must be rejected |
| `ceres-cut` | 11 × 256 × 256 px, in the repository | null case: real noise around saturated, nearly static Ceres; nothing may pass |
| `mariekebaan-full` | 82 × 695 × 519 px, external | the same asteroid among moving blank borders and a warm pixel (edge and detector tests) |
| `leo-triplet-full` | 20 × 2328 × 1760 px, external | the same asteroid in the full field, and (4456) Mawson at S/N 8.4 (found once the 0.8 calibration screened artefacts out); plate solution transferred from the solved stack (about 6 min on 2 cores) |
| `ceres-full` | 11 × 1391 × 1039 px, external | null case on the full field |

The fixtures in `data/` are 16-bit crops of the full sets (about 6 MB together);
quantisation is at most a quarter of the noise and gives the same candidates as
the floating-point crops. They come from the author's own captures.

## Full-frame data

The `-full` cases are skipped unless `FAINTTRAIL_REGRESSION_DATA` points at a
directory with this layout:

```
$FAINTTRAIL_REGRESSION_DATA/
  mariekebaan/                      registered frames (see mariekebaan-cut.case)
  leo-triplet/                      registered frames (see leo-triplet-cut.case)
  mariekebaan-stack-solved.fits     the stack solved with tools/solve_from_catalog.py
  leo-triplet-stack-solved.fits     the plate-solved stack (ASTAP)
  ceres/                            registered frames (see ceres-cut.case)
  baselines/                        written by --write-baseline
```

Each case file gives the exact `tools/register_sequence.py` command that produced
its frames from the raw captures. Pixel positions in the `-full` cases refer to
that registration (reference frame = middle frame); registering differently moves
them.

## Adding a case

1. Register the frames (`tools/register_sequence.py`, or any tool that writes
   `DATE-OBS`, `EXPTIME` and ideally `REGDX`/`REGDY`).
2. Identify the objects independently (MPChecker, an ephemeris, the observer)
   and write the `[object]` sections; positions are FITS 1-based pixels at `mjd`.
3. For the repository, crop around the objects with
   `fainttrail crop ... --bitpix 16` and keep the data small.
4. Run `fainttrail regress NEW.case --write-baseline`, then commit case, data and
   baseline.
