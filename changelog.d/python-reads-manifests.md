### Changed — the Python package takes component formulas from the engine's component manifest
- The formula → CoolProp name table behind `sim.steps[0].bottom.y["CO2"]` and the other typed mappings is read from the component manifest (once per process) instead of a copy in `model.py`. More formulas and CoolProp aliases are accepted than before.

### Added — `skarstrom.component_aliases()` and reference-client tests
- `skarstrom.component_aliases()` returns the formula and alias → CoolProp name table.
- Tests that act as a reference client of the engine's views: for every bundled deck the deck view validates against its schema and every value sets back by its path; every `sim.edit` cascade and kind (read from `model_kinds()`) leaves a deck that passes `check()` and round-trips; the manifests are non-empty and consistent.
- Docs: "Views and manifests" in the Python guide.
