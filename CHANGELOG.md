# Changelog

All notable changes to Skarstrom are recorded here. The format follows
[Keep a Changelog](https://keepachangelog.com/en/1.1.0/) and the project uses
[Semantic Versioning](https://semver.org/). The single source of the version is
`ADSIM_PROJECT_VERSION` (plus `ADSIM_VERSION_SUFFIX` for pre-releases) in the
root `CMakeLists.txt`.

## [Unreleased]

### Fixed
- Competitive isotherms summed their per-component terms in `std::unordered_map`
  order, which depends on the order of the entries in the deck file and on the
  standard library. A deck and its re-saved copy could therefore run to
  slightly different numbers (seen on macOS arm64, where FMA contraction makes
  even a two-component sum order-dependent). Isotherm parameters are now kept
  in key order (`std::map`); the deck format is unchanged.
- CI: the "example decks unchanged" check no longer uses `A && B || C`
  (shellcheck SC2015) and also fails when the test suite deletes a file.

### Added
- `DeckRoundTrip` tests: every example deck's numbers parse exactly and
  survive save/reload bit for bit, and a run is unchanged when every map in
  the deck is listed in reverse order.

## [0.1.0] — initial open-source release (unreleased)

First public release of the Skarstrom engine under the Apache License 2.0.

### Added
- The `adsim` C++20 simulation engine: 1D finite-volume packed-bed model
  (mass, momentum, energy and wall energy balances; 20 isotherm models;
  linear driving force kinetics with several coefficient correlations;
  dispersion and heat-transfer correlations; CoolProp fluid properties),
  multi-step and multi-bed cycles with valves, pumps, buffer tanks and step
  connections, SUNDIALS IDA integration, cyclic-steady-state detection and
  Aitken acceleration, KPIs, a parameter registry and a dlib optimiser.
- Public headers for front ends: `adsim/objects/deck_io.h` (load decks with
  either cereal root), `adsim/objects/deck_edit.h` (component, adsorbent and
  connection edits that keep a deck consistent),
  `adsim/equations/isotherm_kinds.h` (isotherm kind tags, parameter access,
  evaluation) and `adsim/solver/isotherm_fit.h` (the isotherm fitter).
- CMake build with presets for Linux, macOS and Windows; `find_package(Adsim)`
  install/export; use as a subproject (`add_subdirectory`) with front-end
  friendly defaults; builds with Emscripten.
- GoogleTest suite (~670 tests): unit, integration, analytic verification and
  quantitative validation against ten published studies.
- The `skarstrom` Python package (scikit-build-core + pybind11) and CLI
  (`run`, `params`, `examples list|run|path|copy`, JSON output), bundling all
  example decks and the adsorbent library; structural edits from Python
  (`add_component`, `remove_component`, `add_library_adsorbent`).
- 33 example decks, a runnable Python script per deck family in `examples/`,
  and three executed notebooks.
- A curated, cited adsorbent library (`data/adsorbents/`) and cycle templates
  (`data/cycles/`).
- MkDocs documentation with generated validation and adsorbent-library pages.
- GitHub Actions: CI on Linux, macOS and Windows, the Python wheel, the
  Emscripten build and the docs; release workflow with cibuildwheel wheels,
  sdist, source tarball and gated (Test)PyPI trusted publishing; docs
  deployment to GitHub Pages.
