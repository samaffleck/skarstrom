### Added — molecular diffusion has its own switch in the macropore term
- The resistance model has a new switch, `molecular` (registry leaf
  `kinetics[C].resistances.molecular`), that refines the macropore term next to
  `knudsen`. The pore diffusivity is D_m (molecular only), D_K (Knudsen only) or
  the Bosanquet mean 1/D = 1/D_K + 1/D_m (both). One of the two stays on while
  macropore is on; switching macropore on with neither on switches molecular on.
- `molecular_diffusion_coefficient` belongs to the `molecular` switch.
- `MassTransferResistances` is at cereal version 2. A deck from an earlier version
  loads with `molecular` on wherever `macropore` is on, so its results do not change.
- Python: `set_rate_model(resistances=[...])` keeps molecular diffusion on with
  macropore; switch it off with `set_rate_resistance(component, "molecular", False)`.
  The deck-to-script writer does this for a Knudsen-only site.

### Changed — the coefficient's display equation has one definition per line
- `MassTransferEquationLatex` puts each definition on its own line, aligned on
  the equals sign (it joined them on one line before). The Bosanquet mean and
  D_K are separate lines.
- With the film term on in a packed layer, the deck view's `equation` also shows
  the Sherwood number of the layer's film correlation (Wakao–Funazkri, Gunn or
  Ranz–Marshall) and the Reynolds and Schmidt numbers it uses
  (`FilmCorrelationLatex`).
- The model-kind manifest's resistances have `needs_child` (true for macropore).
