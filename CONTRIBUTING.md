# Contributing to Skarstrom

Thank you for helping. Bug reports, validation cases, adsorbent data, docs and
code are all welcome. Please read the [Code of Conduct](CODE_OF_CONDUCT.md)
first; for security problems see [SECURITY.md](SECURITY.md).

## Getting set up

```bash
git clone --recursive https://github.com/samaffleck/skarstrom.git
cd skarstrom
cmake --preset linux-release         # macos-release, windows (see CMakePresets.json)
cmake --build --preset linux-release
```

You need CMake 3.25+ for the presets (3.15+ without them), a C++20 compiler
(GCC 11+, Clang 14+, Apple clang from Xcode 15+, MSVC 2022), Python 3 (CoolProp
runs it at build time) and, recommended, Ninja. Forgot `--recursive`? Run
`git submodule update --init --recursive`.

## Running the tests

```bash
# C++ (working directory core/; ~670 tests, ~10 minutes on 4 cores)
ctest --preset linux-release -j3
# fast tier: everything but the paper replications (label `paper`) and the slow
# verification tests (label `slow`); CI runs it on Linux, macOS and Windows for
# every pull request, and every tier on main, nightly and manual runs
ctest --preset linux-release -j3 -LE "paper|slow"
rm -rf core/Results                 # the suite writes ~4 GB of output there

# Python package, CLI and example smoke tests
pip install . pytest
pytest python/tests

# Docs (and the generated pages)
pip install -r requirements-docs.txt
python scripts/docs/gen_validation_page.py --check
python scripts/docs/gen_adsorbent_library.py --check
mkdocs build --strict
```

Front-end compatibility: `bash scripts/test_subproject_consumer.sh` builds the
engine with `add_subdirectory`, and `--emscripten` (with emsdk activated) builds
it with Emscripten; `bash scripts/test_install_consumer.sh` checks
`find_package(Adsim)`.

## Making a change

- Open an issue first for anything large, so we can agree on the approach.
- Keep pull requests focused; add or update tests with every behaviour change.
- CI runs on pull requests (a draft is fine) and on main, not on a plain push
  to a branch; to run it on a branch without a pull request, use Actions > CI >
  Run workflow. A new commit on a pull request cancels the run for the last one.
  Physics changes need a test whose expected value is derived from physics
  (`core/test/validation/`) or a published reference (`core/test/papers/`),
  not one recorded from the current output.
- Changes to a serialised class must bump its `cereal_class_version` and keep
  older decks loading.
- The test suite re-saves the example decks (`examples/decks/`) and some
  adsorbent files; small changes to their embedded numbers are expected —
  commit them with your change.
- Update the docs (`docs/`) and add a changelog entry as a **new file**
  `changelog.d/<topic>.md` (see `changelog.d/README.md`), not an edit to
  `CHANGELOG.md`, so pull requests in flight at the same time do not
  conflict. It holds one or more `### Added`, `Changed`, `Deprecated`,
  `Removed`, `Fixed` or `Security` subsections (the kinds choose the next
  version). Do not change the version in `CMakeLists.txt`; the release pull
  request does (see RELEASING.md).
- Keep the engine front-end free: no GUI, web or JSON-schema code in `core/`;
  it must keep building as a subproject and with Emscripten.

## Data integrity (adsorbent and validation data)

The value of Skarstrom's data lies in being traceable. These rules are not
negotiable:

- **Every adsorbent parameter is cited.** A library entry
  (`data/adsorbents/library.json` + `<key>.json`) records, for each isotherm,
  heat of adsorption and material property, the source (DOI, table/figure) and
  its provenance (`source`, `derived` with the formula, `repo`, `default`,
  `not_stated`). A value the source does not state is never filled in silently:
  its reference begins `Not stated in source (checked: ...)`.
  `python scripts/docs/gen_adsorbent_library.py` validates the catalogue against
  the files and must pass; the `AdsorbentLibrary*` tests must pass.
- **One recommended entry per material.** Every entry names its material
  (`material_id`, a key of `materials`), and every material has exactly one
  `recommended` entry with a one-line `recommended_reason`. The recommended
  entry has the best data quality (each isotherm's `data_quality`:
  `measured_tabulated` > `measured_fit` > `refit` > `secondary`), then the best
  coverage of the material's usual components, then the widest range. A
  material gets a second entry only when it is a different product, form or
  data range a user can name: that entry's `variant`. Otherwise the better
  entry replaces the weaker one, whose key goes into the material's `aliases`
  so it keeps resolving. An entry's label is the material name, plus the
  variant in brackets; no author names.
- **Reference values are printed numbers.** `reference_values.csv` in a paper
  folder holds only numbers printed in the paper (tables, values stated in the
  text), each with its table/figure/page. Anything derived is marked as such.
  Never eyeball values off a plot; digitised curves go in `core/data/<paper>/`
  and say so.
- **Tolerances are regression bounds.** Set each from the measured error,
  rounded up modestly, with a comment giving the measured value. Never loosen a
  tolerance to hide a disagreement: document it in the paper's README.
- **No copyrighted paper material.** Do not commit PDFs, supplementary
  documents, figure images or full-text extractions (`.gitignore` blocks the
  usual culprits under `core/test/papers/`).
- **Generated pages stay generated.** After running the paper tests, update the
  validation page with
  `python scripts/docs/gen_validation_page.py --merge build/linux-release/validation`;
  never edit `docs/validation.md` or `docs/adsorbent-library.md` by hand.
- File names in `examples/decks/`, `data/adsorbents/` and `data/cycles/` must be
  shell-safe (`DataFileNamesAreShellSafe` test).

## Developer Certificate of Origin (DCO)

Every commit must be signed off, certifying that you wrote the change or have
the right to submit it under the project's licence, as described by the
[Developer Certificate of Origin 1.1](https://developercertificate.org/):

```bash
git commit -s -m "Fix the Toth isotherm derivative"
```

This adds a `Signed-off-by: Your Name <you@example.com>` line. To sign off
commits you already made: `git rebase --signoff main`. Contributions are
licensed under the Apache License 2.0 (see [LICENSE](LICENSE)).

## Releases

Releases are automatic: once a day, the Daily release workflow releases `main`
to PyPI when it has changelog fragments and its CI is green. Between releases
`main` carries a development version (`0.9.1.dev0`). See [RELEASING.md](RELEASING.md).
