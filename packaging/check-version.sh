#!/usr/bin/env bash
# Usage: packaging/check-version.sh vX.Y.Z
# Fails unless the tag matches project(fainttrail VERSION X.Y.Z) in CMakeLists.txt.
# Prints the version.
set -euo pipefail
tag="${1:?tag}"
here="$(cd "$(dirname "$0")/.." && pwd)"
ver=$(sed -n 's/^project(fainttrail VERSION \([0-9.]*\).*/\1/p' "$here/CMakeLists.txt")
if [ -z "$ver" ]; then echo "cannot read the version from CMakeLists.txt" >&2; exit 1; fi
case "$tag" in
  v*) want="${tag#v}" ;;
  *) echo "not a release tag: $tag (expected v$ver)"; echo "$ver"; exit 0 ;;
esac
if [ "$want" != "$ver" ]; then
  echo "tag $tag does not match the version $ver in CMakeLists.txt" >&2
  exit 1
fi
echo "$ver"
