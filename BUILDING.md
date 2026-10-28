# Building FaintTrail

FaintTrail is a C++17 library with a command-line tool (`fainttrail`), a Qt 6 GUI
(`fainttrail-gui`) and an optional CUDA backend. It builds with CMake on Windows,
Linux and macOS.

| needed for | what |
|---|---|
| everything | CMake ≥ 3.20 (≥ 3.24 for automatic GPU detection), a C++17 compiler |
| multithreading | OpenMP (included with MSVC and GCC; `libomp` on macOS) |
| the GUI | Qt 6 (Widgets, Concurrent), 6.2 or newer |
| the GPU backend | an NVIDIA GPU and the CUDA Toolkit (12.8 or newer for RTX 50xx cards) |

Nothing else is needed to build: FITS reading and writing, registration and the
search are built in.

Always build the **Release** configuration. A Debug build is many times slower.

---

## Windows

### 1. Tools

* **Visual Studio 2022 C++ build tools (MSVC v143).** Install either Visual Studio
  2022 or the free *Build Tools for Visual Studio 2022*:

  ```bat
  winget install Microsoft.VisualStudio.2022.BuildTools --override "--wait --passive --add Microsoft.VisualStudio.Workload.VCTools --includeRecommended"
  ```

  Visual Studio 2026 can stay installed next to it. If you only have VS 2026, CUDA
  12.8 cannot be used: its compiler only accepts MSVC up to VS 2022, and its
  installer does not integrate with VS 2026. A CPU-only build works with VS 2026
  (`-G "Visual Studio 18 2026"`, without `-DFAINTTRAIL_CUDA=ON`).
* **CMake** from cmake.org, or the one included with Visual Studio.
* **Qt 6 for MSVC 2022 64-bit** (for the GUI), e.g. `C:\Qt\6.11.2\msvc2022_64`. Use
  the Qt Online Installer, or without an account:
  `pip install aqtinstall` then `aqt install-qt windows desktop 6.8.3 win64_msvc2022_64 -O C:\Qt`.
  A `mingw_64` Qt kit does not work with MSVC.
* **CUDA Toolkit 12.8 or newer** (for the GPU backend). Install it **after** the
  VS 2022 build tools, so that its *Visual Studio Integration* is installed into
  VS 2022. If CUDA was installed first, run the CUDA installer again, choose
  *Custom* and select only *Visual Studio Integration* (untick the driver and
  PhysX components so that newer versions are not replaced). Check that
  `CUDA 12.8.props` and `CUDA 12.8.targets` are in
  `C:\Program Files (x86)\Microsoft Visual Studio\2022\BuildTools\MSBuild\Microsoft\VC\v170\BuildCustomizations`
  (or the same folder under `2022\Community`); if not, copy them there from
  `C:\Program Files\NVIDIA GPU Computing Toolkit\CUDA\v12.8\extras\visual_studio_integration\MSBuildExtensions`.

### 2. Configure, build, test

From a normal Command Prompt in the source folder (adjust the Qt path and the CUDA
version; leave out the CUDA options for a CPU-only build, the Qt path for a build
without GUI):

```bat
cmake -S . -B build -G "Visual Studio 17 2022" -A x64 -T cuda=12.8 ^
      -DFAINTTRAIL_CUDA=ON -DCMAKE_CUDA_ARCHITECTURES=120 ^
      -DCMAKE_PREFIX_PATH=C:\Qt\6.11.2\msvc2022_64
cmake --build build --config Release -j
ctest --test-dir build -C Release --output-on-failure
```

Check the configure output for:

```
-- FaintTrail: CUDA backend enabled (CUDA 12.8..., architectures 120)
-- fainttrail: building Qt 6.11.2 GUI
-- fainttrail: deploy_gui target uses ...windeployqt.exe
```

`CMAKE_CUDA_ARCHITECTURES` is the compute capability of the card without the dot:
120 for RTX 50xx, 89 for RTX 40xx, 86 for RTX 30xx. Several can be listed for a
binary that runs on other cards: `-DCMAKE_CUDA_ARCHITECTURES="86;89;120"`.

`--config Release` (and `-C Release` for ctest) matters: the Visual Studio
generator builds Debug by default. After changing generator, Visual Studio
version or CUDA version, delete the `build` folder before configuring again.

### 3. Run and deploy

The programs are in `build\Release`:

```bat
build\Release\fainttrail.exe backends
build\Release\fainttrail.exe check-backend --backend cuda
build\Release\fainttrail.exe demo
```

`check-backend --backend cuda` must end with *bit-identical to the reference in 4
of 4 checks*.

The GUI needs the Qt DLLs next to it. Copy them once after every build:

```bat
cmake --build build --config Release --target deploy_gui
build\Release\fainttrail-gui.exe
```

(or run `C:\Qt\6.11.2\msvc2022_64\bin\windeployqt.exe --release build\Release\fainttrail-gui.exe`).
The folder `build\Release` then runs on any 64-bit Windows PC with the Microsoft
Visual C++ Redistributable (x64); the CUDA runtime is built in, so other PCs need
only an NVIDIA driver, and fall back to the CPU without one.

### Common errors

| message | cause and fix |
|---|---|
| `No CUDA toolset found` | the CUDA Visual Studio integration is missing for the selected Visual Studio; see step 1, and use `-G "Visual Studio 17 2022" -T cuda=12.8` |
| `unsupported Microsoft Visual Studio version` (nvcc) | the host compiler is VS 2026; use the VS 2022 build tools |
| `D8016: '/O2' and '/RTC1' are incompatible` | a Debug build with older sources; update, and build with `--config Release` |
| `LNK1149: output filename matches input filename` | older sources with CUDA enabled; update |
| `Qt6Widgets.dll was not found` when starting the GUI | run the `deploy_gui` target (step 3) |
| `fainttrail: Qt 6 not found, GUI disabled` | `CMAKE_PREFIX_PATH` must be the Qt kit folder that contains `lib\cmake\Qt6` |
| `backend: CPU reference ... (CUDA not used: ...)` | the reason is given; `fainttrail backends` shows it too |

---

## Linux

```sh
sudo apt install cmake g++ qt6-base-dev          # Debian/Ubuntu; Fedora: cmake gcc-c++ qt6-qtbase-devel
cmake -S . -B build -DCMAKE_BUILD_TYPE=Release
cmake --build build -j
ctest --test-dir build --output-on-failure
./build/fainttrail demo
```

With CUDA (toolkit from NVIDIA's repository, e.g. packages `cuda-nvcc-12-8` and
`cuda-cudart-dev-12-8`):

```sh
cmake -S . -B build -DCMAKE_BUILD_TYPE=Release -DFAINTTRAIL_CUDA=ON \
      -DCMAKE_CUDA_COMPILER=/usr/local/cuda-12.8/bin/nvcc
cmake --build build -j
./build/fainttrail check-backend --backend cuda
```

Without `CMAKE_CUDA_ARCHITECTURES` the code is compiled for the GPU in the build
machine (`native`, CMake ≥ 3.24). Clang can be used as the CUDA compiler as well
(`-DCMAKE_CUDA_COMPILER=clang++`).

## macOS

```sh
brew install cmake qt libomp
cmake -S . -B build -DCMAKE_BUILD_TYPE=Release \
      -DCMAKE_PREFIX_PATH="$(brew --prefix qt);$(brew --prefix libomp)" \
      -DOpenMP_ROOT="$(brew --prefix libomp)"
cmake --build build -j
ctest --test-dir build --output-on-failure
```

There is no CUDA on macOS; the CPU backend is used.

---

## Build options

| option | default | |
|---|---|---|
| `FAINTTRAIL_BUILD_CLI` | ON | the `fainttrail` command-line tool |
| `FAINTTRAIL_BUILD_GUI` | ON | the Qt GUI (skipped with a message when Qt is not found) |
| `FAINTTRAIL_BUILD_TESTS` | ON | unit tests and the regression cases for `ctest` |
| `FAINTTRAIL_BUILD_TOOLS` | OFF | benchmarks and studies (`fainttrail_benchmark`, `fainttrail_kernel_bench`, ...) |
| `FAINTTRAIL_CUDA` | OFF | the CUDA stacking backend |
| `FAINTTRAIL_NATIVE` | OFF | `-march=native` (GCC/Clang) |
| `CMAKE_CUDA_ARCHITECTURES` | native | GPU architectures to compile for |

## Run-time helpers (optional, not needed to build)

| program | used for | without it |
|---|---|---|
| [ASTAP](https://www.hnsky.org/astap.htm) with a star database (e.g. D50) | plate solving on request (GUI button *Plate solve (ASTAP)*, `--solve`) | positions in pixels only, unless the frames are solved or `--wcs-from` is used |
| `curl` (included in Windows 10 and 11, Linux, macOS) | downloading Gaia DR3 for the photometric calibration | instrumental magnitudes only, unless `--catalog FILE` is given |
| Python 3 with numpy, scipy, astropy | the scripts in `tools/` (registration for the regression data, catalogue solving) | not needed for searching |

## Releases

Ready-made binaries are built by GitHub Actions (`.github/workflows/release.yml`)
whenever a version tag is pushed:

| file | contents |
|---|---|
| `fainttrail-X.Y.Z-windows-x64.zip` | `fainttrail-gui.exe`, `fainttrail.exe`, CUDA backend, Qt and runtimes |
| `FaintTrail-X.Y.Z-x86_64.AppImage` | Linux GUI with Qt inside, CUDA backend |
| `fainttrail-X.Y.Z-linux-x86_64.tar.gz` | Linux command-line tool, CUDA backend |
| `fainttrail-X.Y.Z-macos-arm64.zip` | `FaintTrail.app` and the command-line tool (CPU backend) |
| `SHA256SUMS.txt` | checksums |

The CUDA code in the release builds covers compute capability 7.5 to 12.0 (GTX 16xx
and RTX 20xx up to RTX 50xx) and runs with NVIDIA driver 570 or newer; without
one, the CPU backend is used.

To make a release:

1. Set the version in `CMakeLists.txt` (`project(fainttrail VERSION X.Y.Z ...)`)
   and add it to the version history in DEVELOPMENT.md.
2. Commit, then tag and push the tag; the tag must be `v` plus that version
   (the workflow checks it):
   ```sh
   git tag vX.Y.Z
   git push origin vX.Y.Z
   ```
3. The workflow builds and tests on Windows, Linux and macOS, packages, and
   publishes the release with the files above and generated release notes.

*Actions › release › Run workflow* builds the same packages as downloadable
artifacts without publishing anything, to try the packaging. The packaging
scripts are in `packaging/` and also work locally, e.g.
`bash packaging/package-linux.sh 0.9.0 build dist`.

## Paths in the other documents

The documents in `docs/` write programs as `./build/fainttrail` (Linux, macOS).
On Windows with the Visual Studio generator they are in `build\Release\`, e.g.
`build\Release\fainttrail.exe`, `build\Release\fainttrail_kernel_bench.exe`.

## Regression data

The small fixtures in `tests/regression/data` run with every `ctest`. The full
sequences are optional; see [tests/regression/README.md](tests/regression/README.md)
for the folder layout and `FAINTTRAIL_REGRESSION_DATA`.
