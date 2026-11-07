### Changed — one definition of the gas constant and the standard atmosphere
- `kIdealGasConstant` and the new `kStandardAtmosphere` (101325 Pa) in `adsim/math/constants.h` replace the literals in the expression parser, the boundary-condition defaults and the initial-condition pressure default. Results are unchanged.

### Added — merged adsorbent library file
- `scripts/docs/gen_adsorbent_library.py` also writes `data/adsorbents/library.merged.json`, the validated catalogue with every entry's provenance and serialised adsorbent in one file, in the order of the docs page; `--check` (run in CI) fails when it is out of date.
