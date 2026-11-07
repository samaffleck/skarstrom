### Added — deck view, view schema and display manifests for integrators
- `adsim::DeckView` (`adsim/objects/deck_view.h`; Python `Simulation.deck_view()`): a read-only JSON view of a deck's structure. Each node carries the parameter-registry entries that belong to it (path, value, text, unit, bounds, choices, label). `DeckViewSchema()` (Python `skarstrom.deck_view_schema()`) is its JSON Schema, with the builder's kind names as enumerations.
- Display manifests in `adsim/objects/manifests.h` and the Python package: KPIs (`kpi_manifest`), unit categories with suggested units resolved by the unit parser, plus `check_unit` (`unit_manifest`, `check_unit`), advisory warnings (`warning_manifest`), example decks (`example_manifest`, declared in `examples/example_decks.json`) and components (`component_manifest`). The unit parser exposes `UnitCategoryOf` and `UnitOverrideSymbols`.
- Docs page "Deck View and Manifests".

### Changed — CLI reads KPI columns and example families from the engine
- The `run` KPI table and `kpis.csv` take their columns from the KPI manifest. The `kpis.csv` headers and order are unchanged; the terminal table now shows production rate before the energies. `skarstrom examples --json` takes each deck's family from the example manifest.
