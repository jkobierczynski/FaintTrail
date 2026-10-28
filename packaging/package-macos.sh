#!/usr/bin/env bash
# Usage: packaging/package-macos.sh VERSION BUILD_DIR OUT_DIR
# Makes OUT_DIR/fainttrail-VERSION-macos-arm64.zip: FaintTrail.app (GUI, with Qt and
# OpenMP inside) and the command-line tool with its OpenMP library next to it.
set -euo pipefail
ver="${1:?version}"; build="$(cd "${2:?build dir}" && pwd)"; out="${3:?out dir}"
src="$(cd "$(dirname "$0")/.." && pwd)"
mkdir -p "$out"; out="$(cd "$out" && pwd)"
name="fainttrail-$ver-macos-arm64"
work="$(mktemp -d)"; trap 'rm -rf "$work"' EXIT
dst="$work/$name"; mkdir -p "$dst"

# GUI bundle
cp -R "$build/fainttrail-gui.app" "$dst/FaintTrail.app"
"$(brew --prefix qt)/bin/macdeployqt" "$dst/FaintTrail.app" -always-overwrite

# command-line tool: carry libomp along and point the binary at it
cp "$build/fainttrail" "$dst/"
omp=$(otool -L "$dst/fainttrail" | awk '/libomp/ {print $1}' | head -1)
if [ -n "$omp" ]; then
  cp "$(echo "$omp" | sed 's|@rpath|'"$(brew --prefix libomp)"'/lib|')" "$dst/libomp.dylib"
  install_name_tool -change "$omp" "@executable_path/libomp.dylib" "$dst/fainttrail"
fi
codesign --force --sign - "$dst/fainttrail" "$dst/libomp.dylib" 2>/dev/null || true
codesign --force --deep --sign - "$dst/FaintTrail.app" 2>/dev/null || true

cp "$src/LICENSE" "$src/README.md" "$src/BUILDING.md" "$dst/"
cat > "$dst/README-BINARY.txt" <<TXT
FaintTrail $ver for macOS (Apple silicon)

FaintTrail.app   the GUI
fainttrail       command-line tool (run ./fainttrail --help in Terminal)

The binaries are not notarised. On first start macOS refuses them; allow them in
System Settings > Privacy & Security, or run once in Terminal:
  xattr -dr com.apple.quarantine FaintTrail.app fainttrail libomp.dylib
There is no CUDA on macOS; the CPU backend is used.
TXT
(cd "$work" && ditto -c -k --keepParent "$name" "$out/$name.zip")
ls -l "$out"
