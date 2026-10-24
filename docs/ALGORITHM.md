# How FaintTrail works

FaintTrail looks for objects that are too faint to see in any single frame by
adding many short exposures along trial motions. If a trial motion matches the
motion of a real object, the object's light adds up in one place while the noise
averages down; the signal-to-noise ratio grows roughly as the square root of the
number of frames.

The hard part is not the stacking but knowing when a peak is real. A search over
thousands of trial motions and millions of pixels produces strong noise peaks by
chance. FaintTrail therefore treats false-alarm control as a first-class part of
the pipeline.

## 1. Preprocessing (`prepare`)

Input: registered (aligned) frames with mid-exposure times.

1. **Background.** The median of each frame is subtracted, which absorbs changes
   in sky level and transparency.
2. **Static sky.** For every pixel the median over time is the static sky (stars,
   galaxies, hot pixels). The residual of frame *i* uses a median that leaves out
   the frames closest in time to *i* (`static_exclude`, default N/8 on each side),
   so a moving object does not subtract part of itself from its own frames.
3. **Noise.** Each residual frame gets a robust noise estimate (MAD).
4. **Mask and weights.** The temporal scatter of every pixel, in units of the frame
   noise, is measured. Pixels that scatter far more than noise (stellar cores
   under changing seeing, bad columns) are masked; moderately noisy ones are
   down-weighted.
5. **Cosmic rays.** Pixels more than `cr_sigma` above the noise whose neighbours are
   much fainter than a stellar profile allows are removed from that frame.
6. **Clipping and matched filter.** Residuals are clipped (`clip_sigma`) and
   convolved with a Gaussian of the measured PSF width, normalised so that the
   result is a unit-variance S/N image. The PSF width is measured from the static
   sky with adaptive (Gaussian-weighted) moments of isolated stars.

## 2. Search (`search`)

For every trial velocity **v** on a grid (spacing `4 * grid_tol_px / span`), each
matched-filtered frame is shifted by `round(v * t_i)` and the frames are summed:

    S_v(p) = sum_i g_i(p + v t_i) / sqrt(n(p))

`n(p)` is the number of frames that overlap pixel `p` at this velocity. Only the
per-pixel maximum over all velocities and the velocity that produced it are kept.
Local maxima of this max-S/N map are the raw detections.

The computation is limited by memory traffic, so the CPU backend is cache-blocked
(see [PROFILING.md](PROFILING.md)): the image is processed in bands of 8 rows, each
band by one thread; within a band the velocities are taken 8 or more at a time, and
each frame's rows are added into small per-velocity accumulators that stay in the
L2 cache. `n(p)` follows from the frames' column intervals, which split every row
into segments of constant coverage. The summation order per pixel is unchanged, so
the result is bit-identical to the plain definition. A GPU backend would use the
same layout, with a band in shared memory.

## 3. Threshold calibration

The search is repeated `null_runs` times with the frame **times randomly permuted**.
Permuting the times destroys every coherent linear motion but keeps everything else
about the data: the noise, the cosmic rays, the star residuals, the hot pixels. The
complete time-dependent pipeline (static sky, residuals, mask, filter) is rerun for
each permutation, so the null searches see exactly the same estimator as the real
one. (Reusing the real-time preprocessing would be wrong: the leave-out static sky
correlates the noise of neighbouring frames, and only a full rerun reproduces that
correlation in the null.)

The peaks of the null searches give the tail of the noise distribution. A straight
line is fitted to `log N(>S)` against `S` for 0.67 to 30 noise peaks per search,
and the threshold is set where the line predicts `target_false_per_run` peaks per
search (default 0.02, one false detection in fifty searches). Every candidate also
reports the expected number of noise peaks at least as strong as itself.

## 4. Refinement and tests

Each peak above the threshold is refined on a finer grid of sub-pixel positions and
motions using bilinear interpolation, followed by parabolic interpolation. Then:

* **half**: the first and second half of the sequence are stacked separately; a real
  object shows up in both with about S/N/√2 each. Artefacts that live in one half
  fail.
* **single**: the share of the signal coming from the strongest single frame. A
  real object spreads its signal evenly; a defect does not.
* **ghost**: weaker detections whose tracks overlap a brighter object's track and
  get much of their signal from the overlapping frames are sidelobes of that object.
  Detections that coincide with a brighter one in more than half of the frames are
  dropped as duplicates.
* informational flags: **masked** (track mostly over masked pixels), **vedge**
  (at the fastest searched motion, the object may be faster), **slow** (moves less
  than two PSF widths during the sequence).

A candidate passes (`ok`) when it is above the threshold and has none of the
rejecting flags (half, single, ghost).

## 5. Uncertainties

Position uncertainty per axis is `sqrt(2) * sigma_psf / (S/N)`, the precision of a
matched-filter centroid. The motion uncertainty follows from the time sampling:
`pos_err / rms(t)`. The benchmark (`tools/benchmark.cpp`) checks these against the
actual errors on simulated data.

## 6. Tiled search (`search_tiled`)

For sequences that do not fit in memory the field is cut into square tiles and
searched tile by tile, with frames read from the FITS files a rectangle at a time
(`FrameSource`). The result is identical to the untiled search, bit for bit:

* Every tile is prepared and stacked with a halo: the longest shift of any trial
  motion over the sequence (the *stack halo*), plus the reach of the matched
  filter, the cosmic-ray test and the mask dilation (the *filter halo*), plus the
  3 px growth of blank borders. Within the core of a tile nothing then depends on
  pixels outside the tile, and the coverage count sees the true image edges.
* Quantities that belong to a whole frame are measured on the whole frame first,
  from exactly the sample positions the untiled search uses: the background level,
  the residual noise (residuals are computed only at the sample positions), the PSF
  width (from the static sky of the whole field) and the matched-filter
  normalisation (one extra filtering pass over the tiles).
* The detection maps of the tile cores are assembled into whole-field maps before
  peaks are found, so peaks on tile boundaries are found once.
* Each null run uses the same time permutation as the untiled search and goes
  through the same three passes.
* Candidates are refined on a tile around them with a refinement halo that covers
  the refined track, the photometry and edge boxes and the detector hypothesis
  (which follows the registration offsets).

The price is the extra preprocessing (about twice the untiled amount) and the
stacking of the halos. The tile size is chosen from a memory budget
(`TileOptions::max_memory_mb`) or given directly; the command-line tool tiles
automatically when the untiled search would need more than 60 % of the RAM.

## 7. Known limitations

* Motion is linear. Fast near-Earth objects observed over long sequences curve;
  splitting the sequence or adding a curvature term is planned.
* Objects that move less than about two PSF widths during the sequence are partly
  absorbed into the static sky. Longer sequences help.
* Memory: the untiled search holds about three float copies of the sequence at its
  peak. Larger sequences are searched in tiles (section 6), which costs time.
* Long sequences (many frames) spend much time in the leave-out median of the
  static sky, which costs O(n²) per pixel; in a tiled search it is computed about
  twice.
* The frames must already be calibrated and registered. FaintTrail does not align.
