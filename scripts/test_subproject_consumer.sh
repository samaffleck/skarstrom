#!/usr/bin/env bash
# Prove that a front end can build the engine as a subproject:
#   add_subdirectory(<this checkout>) + target_link_libraries(... Adsim::adsim)
# and run it on an example deck. With --emscripten the same consumer is built
# with emcmake (the engine compiled to WebAssembly) and run under Node.js.
#
# Usage: bash scripts/test_subproject_consumer.sh [--emscripten] [deck.json]
# Env knobs: BUILD_DIR (default build_subproject[_emscripten]), JOBS (default 2),
# CMAKE_ARGS (extra configure arguments, e.g. -DADSIM_EMSCRIPTEN_EXCEPTIONS=js).
# For --emscripten, emcmake/emcc and node must be on PATH (source emsdk_env.sh).
set -euo pipefail

SCRIPT_DIR="$(cd "$(dirname "$0")" && pwd)"
ROOT="$(cd "$SCRIPT_DIR/.." && pwd)"
cd "$ROOT"

EMSCRIPTEN=0
if [ "${1:-}" = "--emscripten" ]; then EMSCRIPTEN=1; shift; fi
DECK="${1:-examples/decks/Ward2022_Q100ccm.json}"
JOBS="${JOBS:-2}"
GEN=()
command -v ninja >/dev/null 2>&1 && GEN=(-G Ninja)

if [ "$EMSCRIPTEN" = 1 ]; then
  BUILD_DIR="${BUILD_DIR:-build_subproject_emscripten}"
  # shellcheck disable=SC2086 # CMAKE_ARGS is a list of arguments
  emcmake cmake -S core/subproject_consumer -B "$BUILD_DIR" "${GEN[@]}" -DCMAKE_BUILD_TYPE=Release ${CMAKE_ARGS:-}
  cmake --build "$BUILD_DIR" -j "$JOBS" --target frontend
  node "$BUILD_DIR/frontend.js" "$DECK"
else
  BUILD_DIR="${BUILD_DIR:-build_subproject}"
  cmake -S core/subproject_consumer -B "$BUILD_DIR" "${GEN[@]}" -DCMAKE_BUILD_TYPE=Release
  cmake --build "$BUILD_DIR" -j "$JOBS" --target frontend --config Release
  EXE="$BUILD_DIR/frontend"
  [ -x "$EXE" ] || EXE="$BUILD_DIR/Release/frontend.exe"
  "$EXE" "$DECK"
fi
echo "subproject consumer: OK"
