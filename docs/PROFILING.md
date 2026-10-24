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
  machine's memory bandwidth, the simple whole-image kernel and the cache-blocked
  backend with a range of block sizes, on random frames of any size. (In 0.4 it
  compared prototypes; the tables below come from that version.)

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

## After step 4: the cache-blocked CPU backend (0.5.0)

Item 1 is implemented: the CPU backend works in bands of 8 rows and blocks of
8–16 motions (chosen from the image width so the accumulators take about 600 kB),
and computes coverage from the frames' column intervals, so a row becomes a few
segments of constant coverage with a plain, vectorisable maximum inside. Each band
belongs to one thread, so the per-thread result images and their merge are gone.

The old (0.4) and new (0.5) builds were run back to back on the same machine,
with `--strict` against the 0.4 baselines. All candidates, verdicts and
thresholds are identical, as are the unit tests and `check-backend`, which
compare bit for bit with two independent implementations.

| case | stacking 0.4 | stacking 0.5 | speed-up | throughput 0.5 | total 0.4 → 0.5 |
|---|---:|---:|---:|---:|---:|
| Leo Triplet (20 × 2328 × 1760) | 551 s | 147 s | **3.7×** | 7.0 G adds/s | 578 s → **175 s** |
| Ceres (11 × 1391 × 1039) | 76.8 s | 26.3 s | 2.9× | 5.6 G adds/s | 81 s → 31 s |
| Mariekebaan (82 × 695 × 519) | 27.9 s | 13.4 s | 2.1× | 6.4 G adds/s | 48 s → 35 s |

(The machine was slower on this day than in the first measurement above: the 0.4
build reached 1.9 G adds/s instead of 2.55.)

The kernel benchmark gives the same picture (2328 × 1760, 20 frames): 1.57 G
adds/s for the simple kernel and 4.5 G adds/s for the blocked one on 2 threads.
Block sizes from 4 to 16 rows and 8 to 16 motions are within about 15 % of each
other; `fainttrail_kernel_bench` lists them, so the defaults can be checked on
other machines.

What is left, for the Leo run:

| stage | share after 0.5 |
|---|---:|
| stacking, null runs | 62 % |
| stacking, real data | 21 % |
| null-run preprocessing | 11 % |
| preprocessing and the rest | 6 % |

Inside the stacking, accumulation is now 82 % and coverage + maximum 18 %. The
next gains are the cheaper null calibration (item 2) and, for long sequences, the
leave-out median (item 4); memory (item 3) is unchanged at 2.1 GB.

## After step 4: memory layout (0.6.0)

Item 3. In 0.5 a search held six copies of the sequence at its peak: the frames
as loaded; the background-subtracted frames; the residuals (frame − background −
static sky); the matched-filtered frames; and, during each null run, a second set
of residuals and filtered frames. Now:

* **Residuals are dropped as soon as each frame is filtered** (the noise and
  scatter maps are computed first). Only three things used them later: the edge
  test, the photometry cut-outs and the GUI's residual view. The edge test now
  reads blank pixels from the background-subtracted frames and the static sky,
  which is exactly equivalent. The photometry cut-outs recompute the residuals of
  the few pixels they need with the same leave-out median, so the results are
  bit-identical. The GUI and `fainttrail stack` keep the residuals
  (`SearchParams::keep_residuals`).
* **The null runs keep only their filtered frames**, and their residuals are also
  released one frame at a time.
* **`prepare(std::move(frames))` takes over the loaded pixel buffers** instead of
  copying them; the CLI, the GUI and `fainttrail regress` use it.
* **Freed memory goes back to the system.** With glibc, large freed blocks stayed
  with the process: glibc raises its mmap threshold after a large block is freed,
  so later frame-sized buffers came from the heap and were never returned. The CLI
  and GUI call `tune_allocator_for_images()` at start-up, which fixes the threshold
  at 1 MB; the library returns freed heap memory after preparation and between
  null runs. A host application embedding the library can call the same function.
* `--profile` shows the resident memory after every stage, and `fainttrail
  regress` reports each case's memory (also as a `memory_mb` column in `--report`).

Peak resident memory, 0.5 → 0.6, with identical candidate lists (compared
byte for byte, including fluxes) and all regression cases passing in strict mode:

| case | 0.5 | 0.6 | change |
|---|---:|---:|---:|
| Leo Triplet (20 × 2328 × 1760, 328 MB of frames) | 2073 MB | 1148 MB | −45 % |
| Mariekebaan (82 × 695 × 519, 118 MB) | 701 MB | 365 MB | −48 % |
| Ceres (11 × 1391 × 1039, 64 MB) | 443 MB | 261 MB | −41 % |

Run time is unchanged (Leo: 188 s for both builds on the same day).

The peak is now during a null run's preprocessing: the background-subtracted
and filtered frames of the real data (2 copies), plus that run's residuals while
they are being filtered (1 copy), plus working buffers. After the search, 2
copies remain. The unbinned Leo sequence (4656 × 3520) would now need about
4.6 GB instead of 8.3 GB. Going below 3 copies needs spatial tiling (part 3 of
step 4).

## After step 4: spatial tiling (0.7.0)

Item 3, part 2: `search_tiled()` searches the field in square tiles, reading the
frames from the FITS files one rectangle at a time (see ALGORITHM.md, section 6).
Whole-frame quantities are measured first from the same sample positions, every
tile carries the halos it needs, and the detection maps are assembled before the
peaks are found, so **the candidate lists are byte-identical to the untiled
search**. This is checked by a unit test (several tile sizes, from memory and
from FITS files, with blank borders and registration offsets), by the three
fixtures run again as tiled searches in `ctest`, and on the full data below.

| case | untiled | tiled | result |
|---|---|---|---|
| Leo Triplet, `--max-memory 500` → 2 × 2 tiles of 1280 px | 1148 MB, 188 s | **381 MB**, 284 s | identical |
| Leo Triplet, `--tile 512` → 5 × 4 tiles | 1148 MB, 188 s | **209 MB**, 311 s | identical |
| Mariekebaan, `--tile 256` → 3 × 3 tiles | 365 MB, 35 s | **177 MB**, 112 s | identical |

Where the extra time goes in the 2 × 2 run of the Leo Triplet (256 s in the profiled
run against 188 s untiled):

| stage | untiled | tiled |
|---|---:|---:|
| preprocessing, all four time orders | 27 s | 70 s (noise pass 6 s, normalisation pass 32 s, tile preparation 32 s) |
| stacking | 147–160 s | 172 s (the stack halo adds about 7 % of area here) |

Tiling costs most where preprocessing dominates: Mariekebaan has 82 frames, and
its leave-out median (O(n²) per pixel) is computed about twice per time order and
over the halos too. Item 4 (a faster leave-out median) would remove most of that
difference. With smaller tiles the halos grow relative to the cores, so the
command-line tool picks the largest tile that fits the memory budget and does not
tile at all when the whole sequence fits.

The unbinned Leo sequence (4656 × 3520, about 4.6 GB untiled) can now be searched
within any budget down to a few hundred MB.

## Reproducing

```sh
cmake -S . -B build -DFAINTTRAIL_BUILD_TOOLS=ON && cmake --build build -j
FAINTTRAIL_REGRESSION_DATA=/path/to/data \
  ./build/fainttrail regress tests/regression/leo-triplet-full.case --profile
./build/fainttrail_kernel_bench --width 2328 --height 1760 --frames 20 --velocities 48
```
