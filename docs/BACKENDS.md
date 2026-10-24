# Search backends

Almost all of FaintTrail's search time is spent in one operation: stacking the
matched-filtered frames along every trial motion and keeping, per pixel, the best
S/N and the motion that produced it. That operation sits behind a small interface,
`fainttrail::SearchBackend` in [`include/fainttrail/backend.hpp`](../include/fainttrail/backend.hpp),
so that faster implementations (CUDA, Vulkan compute, ...) can be added without
touching the rest of the pipeline.

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
  evaluated in double precision.
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

## Adding a backend

1. Implement `SearchBackend` (`name()`, `description()`, `run()`), normally in
   `src/backend_<name>.cpp`, compiled only when its toolkit is found by CMake.
2. Register it in `backend_names()` / `make_backend()` in `src/backend.cpp`,
   before `"cpu"` so that `auto` prefers it, and only when a usable device is present.
3. Run `fainttrail check-backend --backend <name>` and `ctest`; add the backend to
   the conformance loop in `tests/test_main.cpp`.
4. Run the real-data regression cases with the new backend:
   `fainttrail regress tests/regression --backend <name> --strict`. The candidates,
   verdicts and thresholds must match the CPU baselines (see
   [tests/regression/README.md](../tests/regression/README.md)).

Progress is reported as the fraction of motions processed; returning `false`
from the callback must stop the pass promptly, after which `run()` throws
`fainttrail::Cancelled`.
