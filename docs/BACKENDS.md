# Search backends

Almost all of FaintTrail's search time is spent in one operation: stacking the
matched-filtered frames along every trial motion and keeping, per pixel, the best
S/N and the motion that produced it. That operation sits behind a small interface,
`fainttrail::SearchBackend` in [`include/fainttrail/backend.hpp`](../include/fainttrail/backend.hpp),
so that faster implementations can be added without touching the rest of the
pipeline. Since 0.9 there are two: `cpu` (the reference) and `cuda`.

```
fainttrail backends                       # list backends in this build
fainttrail search frames/ --backend cpu   # choose one (default: auto = best available)
fainttrail check-backend --backend NAME   # verify one against the contract
```

The GUI has a *Backend* selector in the Detection group; the C API has
`ft_params.backend`.

## The contract

The header defines the computation exactly; in short:

* Shifts are integers, `round(v * t_i)` with halves rounded away from zero,
  evaluated in double precision, or, since 0.8.2, taken from an explicit
  `ShiftTable` (one integer shift per trial motion and frame) when the pass carries
  one. Tables allow curved, ephemeris-guided or parallax-corrected tracks;
  `linear_shift_table()` reproduces the straight-line rule exactly. A backend reads
  every shift through `pass_shift()` and must support both forms.
* A frame covers an output pixel when the shifted source pixel lies inside it.
  Pixels covered by fewer than `min_count` frames are undefined for that motion.
* `S_v = sum(g_i) / sqrt(c)`, summed in single precision in ascending frame order.
* Output: the maximum over all motions and the index of the maximising motion;
  0 and -1 where no motion is defined; ties go to the lowest motion index.

The CPU backend (cache-blocked, see [PROFILING.md](PROFILING.md)) follows this bit
for bit. It is tested against two independent implementations: a deliberately
naive one, `run_pass_naive` (one pixel at a time, including tie-breaking and
partial coverage), and `run_pass_simple`, the straightforward whole-image kernel
of FaintTrail 0.2–0.4, which is fast enough to check other backends on realistic
image sizes. Every block size gives the same result.
Other backends may sum in a different order: their S/N must agree within
1e-4 · max(1, |S/N|), and the chosen motion may differ only where the two
motions' S/N values are within that same tolerance. `check_backend()` and the
`check-backend` command test exactly these rules.

## The CUDA backend (0.9)

### Building

Needs the CUDA Toolkit; for an RTX 50xx card (Blackwell, compute capability 12.0)
version 12.8 or newer. The backend is off by default, so builds without CUDA are
unchanged. Full instructions, including the Visual Studio 2022 requirement on
Windows, are in [BUILDING.md](../BUILDING.md).

Linux:

```sh
cmake -S . -B build -DCMAKE_BUILD_TYPE=Release -DFAINTTRAIL_CUDA=ON
cmake --build build -j
./build/fainttrail backends                       # "cuda" should be listed first
./build/fainttrail check-backend --backend cuda   # must say: bit-identical in 4 of 4 checks
ctest --test-dir build --output-on-failure        # regression cases run on the GPU ("auto")
```

Windows (Visual Studio 2022 build tools, CUDA 12.8 with its Visual Studio
integration):

```bat
cmake -S . -B build -G "Visual Studio 17 2022" -A x64 -T cuda=12.8 -DFAINTTRAIL_CUDA=ON -DCMAKE_CUDA_ARCHITECTURES=120
cmake --build build --config Release -j
build\Release\fainttrail.exe check-backend --backend cuda
ctest --test-dir build -C Release --output-on-failure
```

Without `CMAKE_CUDA_ARCHITECTURES` the code is compiled for the GPU in the build
machine (`native`, CMake ≥ 3.24). For a binary that runs on other cards, list them,
e.g. `-DCMAKE_CUDA_ARCHITECTURES="86;89;120"`. The CUDA runtime is linked
statically; only the NVIDIA driver is needed at run time. A binary built with CUDA
runs on machines without a GPU or driver: `auto` then falls back to the CPU and
`fainttrail backends` says why. Verified on an RTX 5080 (Windows 11, CUDA 12.8,
MSVC 14.44).

### Choosing and options

`auto` (the default everywhere) picks `cuda` when a device can be used, else `cpu`.
`--backend cpu` forces the CPU. Options follow a colon:

| name | meaning |
|---|---|
| `cuda:device=1` | use device 1 (default: `FAINTTRAIL_CUDA_DEVICE` or 0) |
| `cuda:max_mb=256` | cap on device memory for frames; smaller values force more row bands |
| `cuda:launch_adds=1e9` | additions per kernel launch (granularity of progress and cancel) |

### How it works

* **One thread per output pixel**, blocks of 32 × 8 threads, so a warp reads 32
  consecutive pixels of a frame row: coalesced loads, and the shifts are the same
  for the whole warp (broadcast from cache).
* Each thread walks the trial motions in ascending order, 8 at a time with the 8
  sums and coverage counts in registers, and for each group of 8 visits the frames
  in ascending order. Every sum is therefore formed in exactly the order the
  contract prescribes, and the result is **bit-identical to the CPU backend**, not
  just within tolerance. (Compiled with `--fmad=false`; the kernel has no
  multiply-add, the flag protects against future edits.)
* Coverage is decided per frame inside the kernel (is the shifted pixel inside the
  frame?), so no coverage map is needed; S/N = sum × 1/√c from a table computed on
  the host in double precision, as on the CPU.
* **Row bands.** The frames are uploaded only for the rows a band of output rows
  can reach (band ± the extreme vertical shifts). Normally the whole sequence fits
  in one band (the binned Leo Triplet sequence uses under 1 GB of the 16 GB of an
  RTX 5080); larger ones are processed band by band within 80 % of free device
  memory.
* **Launches** cover a chunk of motions (about 4·10⁹ additions each, a few
  milliseconds); every thread keeps its running maximum in device memory between
  launches. Progress and cancellation are checked between launches.
* Timings reported to `--profile`: upload, kernel, download, setup.

The per-pixel code (`src/cuda_stack_core.hpp`) is shared with a host emulation
(`emulate_cuda_pass` in `src/cuda_plan.cpp`) that executes the complete GPU
algorithm, bands and launches included, on the CPU. The unit tests run it on every
machine against the naive contract with 1-row bands, single-group launches, NaN
input and shift tables far outside the frame, so the band and launch logic is
tested even where no GPU is present. On a machine with a GPU the same test also
runs the real backend.

### Not done yet

The kernel is a straightforward first version: shifts are read from global memory
(cached), frames are uploaded for each pass (each null run has its own frames
anyway), and nothing has been tuned on real hardware. Measurements on the RTX 5080
come first; the candidates for tuning are shifts in shared memory, more motions
per thread, and keeping frames on the device between the passes of one search.

## Adding a backend

1. Implement `SearchBackend` (`name()`, `description()`, `run()`), normally in
   `src/backend_<name>.cpp`, compiled only when its toolkit is found by CMake.
2. Register it in `backend_names()` / `make_backend()` in `src/backend.cpp`,
   before `"cpu"` so that `auto` prefers it, and only when a usable device is present.
3. Run `fainttrail check-backend --backend <name>` (it includes a curved shift
   table) and `ctest`; add the backend to
   the conformance loop in `tests/test_main.cpp`.
4. Run the real-data regression cases with the new backend:
   `fainttrail regress tests/regression --backend <name> --strict`. The candidates,
   verdicts and thresholds must match the CPU baselines (see
   [tests/regression/README.md](../tests/regression/README.md)).

Progress is reported as the fraction of motions processed; returning `false`
from the callback must stop the pass promptly, after which `run()` throws
`fainttrail::Cancelled`.
