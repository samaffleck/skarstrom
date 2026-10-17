#!/usr/bin/env bash
# Prove the adsim install/export works end to end:
#   1. configure + build the library (no tests) with ADSIM_INSTALL=ON,
#   2. install it into a throw-away prefix,
#   3. configure + build adsim/install_consumer against ONLY that prefix
#      (find_package(Adsim)), and
#   4. run the consumer on an example deck.
#
# Usage: bash scripts/test_install_consumer.sh [deck.json]
# Env knobs (optional):
#   BUILD_DIR=...    reuse an existing build tree instead of build_install/
#                    (it must be configured with ADSIM_INSTALL=ON)
#   PREFIX=...       install prefix (default: a fresh mktemp -d, removed on exit)
#   JOBS=...         parallel build jobs (default 2)
set -euo pipefail

SCRIPT_DIR="$(cd "$(dirname "$0")" && pwd)"
ROOT="$(cd "$SCRIPT_DIR/.." && pwd)"
cd "$ROOT"

BUILD_DIR="${BUILD_DIR:-build_install}"
JOBS="${JOBS:-2}"
DECK="${1:-examples/decks/Ward_2022_Q=50[ccm].json}"
WORK="$(mktemp -d)"
PREFIX="${PREFIX:-$WORK/prefix}"
trap 'rm -rf "$WORK"' EXIT

GEN=()
command -v ninja >/dev/null 2>&1 && GEN=(-G Ninja)

echo "== [1/4] Build adsim ($BUILD_DIR) =="
if [ ! -f "$BUILD_DIR/CMakeCache.txt" ]; then
  cmake -S . -B "$BUILD_DIR" "${GEN[@]}" \
    -DCMAKE_BUILD_TYPE=Release \
    -DADSIM_BUILD_TESTS=OFF -DADSIM_INSTALL=ON
fi
cmake --build "$BUILD_DIR" --config Release -j "$JOBS"

echo "== [2/4] Install into $PREFIX =="
cmake --install "$BUILD_DIR" --config Release --prefix "$PREFIX" > "$WORK/install.log"
test -f "$PREFIX/lib/cmake/Adsim/AdsimConfig.cmake" \
  || test -f "$PREFIX/lib64/cmake/Adsim/AdsimConfig.cmake" \
  || { echo "AdsimConfig.cmake was not installed"; exit 1; }

echo "== [3/4] Configure + build the consumer against the install tree only =="
# CMAKE_FIND_USE_PACKAGE_REGISTRY etc. off: the package must come from PREFIX.
cmake -S adsim/install_consumer -B "$WORK/consumer" "${GEN[@]}" \
  -DCMAKE_BUILD_TYPE=Release \
  -DCMAKE_PREFIX_PATH="$PREFIX" \
  -DCMAKE_FIND_USE_PACKAGE_REGISTRY=OFF \
  -DCMAKE_FIND_USE_SYSTEM_PACKAGE_REGISTRY=OFF
cmake --build "$WORK/consumer" --config Release -j "$JOBS"

echo "== [4/4] Run the consumer on $DECK =="
CONSUMER="$WORK/consumer/adsim_consumer"
[ -x "$CONSUMER" ] || CONSUMER="$WORK/consumer/Release/adsim_consumer.exe"
"$CONSUMER" "$DECK"
echo "install consumer: OK"
