# Profiling

Where does the time go in a real search, and what should be optimised before a
GPU backend is written? Measured on the full Leo Triplet regression case
(`tests/regression/leo-triplet-full.case`), with the other two full-frame cases for
comparison.

**Machine.** A cloud VM: 2 cores of an Intel Xeon at 2.1 GHz, 2 MB L2 per core,
260 MB shared L3, 7 GB RAM; GCC, Release build (`-O3`, portable SSE2 code unless
noted). The VM is noisy: repeated runs differ by up to 15 %. A desktop with more
cores and a smaller L3 will give different absolute numbers, but the shape of the
profile is what matters here.

**Tools.** Two instruments were added for this:

* `--profile` on `fainttrail search` and `fainttrail regress` prints the wall time
  of every pipeline stage (null-run stages summed), the CPU backend's breakdown of
  each stacking pass, and the stacking throughput. The data are also available in
  the library (`Prepared::timings`, `SearchResult::timings`, `StackResult::timings`).
* `fainttrail_kernel_bench` (built with `-DFAINTTRAIL_BUILD_TOOLS=ON`) measures the
  machine's memory bandwidth, the production stacking kernel, and prototype
  cache-blocked variants of it on random frames of any size.

## The full Leo Triplet run

20 frames of 2328 × 1760 px (2 × 2 binned), 3160 trial motions, 3 null runs.

| stage | time | share |
|---|---:|---:|
| load frames + plate-solution transfer | 33.3 s (cold disk cache; 1.1 s warm) | 7.2 % |
| preprocessing (background, static sky, mask, matched filter) | 6.2 s | 1.3 % |
| **stacking, real data** | **101.3 s** | **21.9 %** |
| null-run preprocessing (3×) | 17.0 s | 3.7 % |
| **stacking, null runs (3×)** | **305.6 s** | **65.9 %** |
| peaks, threshold fit, refinement, ghost tests | 0.3 s | 0.1 % |
| **total** | **463.6 s** | |

Inside the stacking passes (per-thread time, averaged over threads):
accumulating the shifted frames 341 s (84 %), coverage count and running maximum
66 s (16 %), merging the per-thread results 0.2 s.

The stacking throughput is **2.55 G pixel-additions per second** on 2 cores
(1.04 T additions in total). Peak memory is **2.1 GB**. User CPU time was 843 s
for 464 s wall, so both cores were busy for almost the whole run.

The other full-frame cases have the same shape with two differences:

| case | stacking | null rebuild | other notable stage |
|---|---:|---:|---|
| Leo Triplet (20 × 4.1 Mpx) | 88 % | 4 % | |
| Mariekebaan (82 × 0.36 Mpx) | 61 % | **28 %** | static sky 7 % |
| Ceres (11 × 1.4 Mpx, fast motions) | 83 % | 5 % | coverage + max **30 %** of stacking |

## The stacking kernel

`fainttrail_kernel_bench` on frames of the Leo size (2328 × 1760, 20 frames,
48 motions):

| variant | 1 thread | 2 threads | result |
|---|---:|---:|---|
| memory bandwidth (stream triad) | 11.8 GB/s | 24.1 GB/s | |
| production backend (one velocity at a time, whole image) | 1.26 G adds/s | 2.13 G adds/s | reference |
| production backend, `-march=native` | 1.41 | 2.45 | identical |
| bands of 8 rows | 3.12 | 5.62 | identical |
| bands of 16 rows, 8 motions per block | 3.18 | 6.07 | identical |
| bands of 8 rows, 16 motions per block | 3.48 | **6.27** | identical |
| bands of 128 rows | 2.33 | 5.25 | identical |

What this shows:

1. **The production kernel is limited by memory traffic, not arithmetic.** For
   every velocity it adds 20 shifted frames into a 16 MB accumulation image, then
   reads that image again for the coverage and maximum. That is about 12 bytes of
   memory traffic per addition; at 2.13 G adds/s this is 25 GB/s, the full
   bandwidth of the machine. Wider SIMD (`-march=native`) gains only 15 %, and
   going from 1 to 2 threads gains 1.7× instead of 2×.
2. **Cache blocking removes most of that traffic.** Processing the image in bands
   of a few rows, with all velocities per band, keeps the accumulation buffer in
   L1/L2; only the frames are streamed. On 2 threads this is **2.6–2.9× faster**
   and scales better with threads (1.8–2.0× from 1 to 2), so the gain should grow
   on a desktop with more cores. The sums are formed in the same order, so the
   results are **bit-identical** (checked with `check_backend`), and the backend
   contract does not change.
3. Band height matters more than the velocity block; 8–16 rows was best here.

## Conclusions for the next steps

In order of expected gain for this run:

| # | change | where it helps | expected effect on the Leo run |
|---|---|---|---|
| 1 | **cache-blocked CPU kernel** (bands of 8–16 rows, velocity blocks, per-row coverage) | all stacking | stacking 407 s → about 150 s; total about 465 s → 210 s |
| 2 | **cheaper null calibration**: the 3 null runs are 75 % of all stacking. Run them on a random subset of the velocity grid and correct the tail fit for the number of trials. This changes the statistics, so it must be re-validated with the false-alarm benchmark first. | all runs | up to a further 2× |
| 3 | **memory**: `Prepared` keeps background-subtracted frames, residuals and filtered frames, and each null run builds a second set, alongside the frames loaded by the caller (2.1 GB here). Dropping residuals after filtering and keeping only what the null runs need cuts this roughly in half. | large sensors | the same sequence unbinned (4656 × 3520 px) would otherwise need about 8 GB |
| 4 | **leave-out median**: the static sky with leave-out costs O(n²) per pixel; with 82 frames it is a third of the Mariekebaan run. A sliding-window form brings it to O(n log n). | long sequences | none here; Mariekebaan −30 % |
| 5 | wider SIMD by default | kernel | ≤ 15 %, not worth losing portability |

The disk-read time (33 s cold, 1 s warm) is the cloud machine, not the code.

**Scaling.** Without binning, the Leo sequence has 4× the pixels and, because the
motion grid is in pixels, about 4× the motions: 16× the work, about 2 hours on
these 2 cores with the current kernel, about 50 minutes with item 1.

**For the GPU backend (step 6).** The kernel is memory-bound, so a GPU should be
judged by memory bandwidth and by how well it keeps the accumulation on chip. The
banded CPU layout is the same structure a CUDA kernel would use (a band in shared
memory, frames streamed through), so item 1 is also the design study for the GPU
kernel. A mid-range GPU with 300–500 GB/s is 10–20× the bandwidth of this VM.

## Reproducing

```sh
cmake -S . -B build -DFAINTTRAIL_BUILD_TOOLS=ON && cmake --build build -j
FAINTTRAIL_REGRESSION_DATA=/path/to/data \
  ./build/fainttrail regress tests/regression/leo-triplet-full.case --profile
./build/fainttrail_kernel_bench --width 2328 --height 1760 --frames 20 --velocities 48
```
