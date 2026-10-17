#!/usr/bin/env bash
# Create a source tarball of the repository INCLUDING the submodule sources
# (git archive leaves submodules out), so the engine builds offline:
#   dist/skarstrom-<version>-full-source.tar.gz  ->  skarstrom-<version>/...
#
# Usage: bash scripts/make_source_tarball.sh [version] [outdir]
# The checkout must have its submodules populated (git submodule update --init --recursive).
set -euo pipefail

SCRIPT_DIR="$(cd "$(dirname "$0")" && pwd)"
ROOT="$(cd "$SCRIPT_DIR/.." && pwd)"
cd "$ROOT"

VERSION="${1:-$(sed -n 's/^set(ADSIM_PROJECT_VERSION \([0-9.]*\)).*/\1/p' CMakeLists.txt)$(sed -n 's/^set(ADSIM_VERSION_SUFFIX "\([a-z0-9]*\)").*/\1/p' CMakeLists.txt)}"
OUTDIR="${2:-dist}"
NAME="skarstrom-${VERSION}"
mkdir -p "$OUTDIR"
OUT="$(cd "$OUTDIR" && pwd)/${NAME}-full-source.tar.gz"

LIST="$(mktemp)"
trap 'rm -f "$LIST"' EXIT
# Every tracked file of the repository and of each (recursive) submodule.
git ls-files --recurse-submodules > "$LIST"
tar --create --gzip --file "$OUT" \
    --transform "s,^,${NAME}/," \
    --owner=0 --group=0 --numeric-owner --mode='u+rwX,go+rX,go-w' \
    --files-from "$LIST"
echo "wrote $OUT ($(du -h "$OUT" | cut -f1), $(wc -l < "$LIST") files)"
