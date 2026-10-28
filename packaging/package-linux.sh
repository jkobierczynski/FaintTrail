#!/usr/bin/env bash
# Usage: packaging/package-linux.sh VERSION BUILD_DIR OUT_DIR
# Makes OUT_DIR/fainttrail-VERSION-linux-x86_64.tar.gz (command-line tool and
# documents) and OUT_DIR/FaintTrail-VERSION-x86_64.AppImage (the GUI with its Qt
# libraries). Needs qmake of the Qt the GUI was built with in QMAKE or on the PATH,
# and downloads linuxdeploy and its Qt plug-in from GitHub.
set -euo pipefail
ver="${1:?version}"; build="$(realpath "${2:?build dir}")"; out="$(realpath -m "${3:?out dir}")"
src="$(cd "$(dirname "$0")/.." && pwd)"
mkdir -p "$out"
work="$(mktemp -d)"
trap 'rm -rf "$work"' EXIT

# --- command-line tool
name="fainttrail-$ver-linux-x86_64"
mkdir -p "$work/$name"
install -m 755 "$build/fainttrail" "$work/$name/"
cp "$src/LICENSE" "$src/README.md" "$src/BUILDING.md" "$src/DEVELOPMENT.md" "$work/$name/"
cp -r "$src/docs" "$work/$name/docs"
cat > "$work/$name/README-BINARY.txt" <<TXT
FaintTrail $ver for Linux x86_64

fainttrail                 command-line tool (run ./fainttrail --help)
FaintTrail-*.AppImage      the GUI, released next to this archive

Needs glibc 2.35 or newer (Ubuntu 22.04, Debian 12, Fedora 36 and later) and the
GNU OpenMP runtime (libgomp1, installed on almost every system). The CUDA backend
is built in and used when an NVIDIA driver and GPU are present; otherwise the CPU
is used. Optional: ASTAP for plate solving, curl for the Gaia photometry.
TXT
tar -C "$work" -czf "$out/$name.tar.gz" "$name"

# --- GUI as AppImage
if [ -x "$build/fainttrail-gui" ]; then
  cd "$work"
  for t in linuxdeploy-x86_64.AppImage linuxdeploy-plugin-qt-x86_64.AppImage; do
    url="https://github.com/linuxdeploy/${t%%-x86_64.AppImage}/releases/download/continuous/$t"
    curl -fsSL -o "$t" "$url"
    chmod +x "$t"
  done
  export APPIMAGE_EXTRACT_AND_RUN=1 NO_STRIP=1
  export LINUXDEPLOY_OUTPUT_VERSION="$ver"
  export EXTRA_PLATFORM_PLUGINS="libqoffscreen.so"  # lets CI start the GUI without a display
  export QMAKE="${QMAKE:-$(command -v qmake6 || command -v qmake)}"
  ./linuxdeploy-x86_64.AppImage --appdir AppDir \
    --executable "$build/fainttrail-gui" \
    --desktop-file "$src/packaging/linux/fainttrail-gui.desktop" \
    --icon-file "$src/packaging/fainttrail.svg" \
    --plugin qt --output appimage
  mv FaintTrail-*.AppImage "$out/"
fi
ls -l "$out"
