### Added — registry leaves, descriptors and a display equation for the mass-transfer resistances
- Each resistance of the coefficient model is a Boolean registry leaf,
  `adsorbent[A].kinetics[C].resistances.<tag>` (`film`, `macropore`, `knudsen` while macropore is on,
  `micropore`, `lumped`), and the capacity-factor mode is the leaf `kinetics[C].capacity_factor_mode`
  (`none`, `slope`, `secant`, `fixed`; while film or macropore is on). Setting them by path goes through
  `SetMassTransferResistance` and `SetCapacityFactorMode`, with their rules. A dual-rate slow site has
  `slow_resistances.<tag>` and `slow_capacity_factor_mode`.
- The parameters keep their flat names (`kinetics[C].k`, `.tortuosity`, `.slow_k`, ...) and are listed only
  while their resistance (or mode) is on. The former arrhenius leaf `kinetics[C].pre_exp_factor`
  (`slow_pre_exp_factor`) still resolves: it is the lumped `k` (`slow_k`).
- The `series_resistance` kind descriptor lists the resistances (tag, label, parent, whether the capacity
  factor scales it, parameters, LaTeX fragments for spheres and slabs, sources), the capacity-factor modes
  and the presets; the six former kinds are marked `legacy`. `ModelKindManifest()` emits the model kinds
  as JSON (also `skarstrom.model_kinds()`); `build::RateCoefficientPresets()` and
  `build_kinds()["rate_coefficient_preset"]` list the preset tags.
- `MassTransferEquationLatex(model, shape)` composes the coefficient's equation from the descriptor's
  fragments. The deck view's kinetics entries have `resistances`, `capacity_factor_mode`, `preset` (the
  derived preset tag, or `custom`), `equation` and `equation_slab`, and the same with `slow_` for a
  dual-rate slow site.
- Python: `set_rate_model(..., resistances=[...], capacity_factor_mode=...)`,
  `Adsorbent.set_rate_resistance(component, tag, on, site="fast")`, and the edits
  `sim.edit.set_mass_transfer_resistance` and `sim.edit.set_capacity_factor_mode`. `to_script` writes the
  switches of a custom model.

### Changed — the coefficient kind is always series_resistance
- A coefficient model reports the kind `series_resistance` (the deck view's `rate_coefficient`,
  `describe()`, `build::RateCoefficientKinds()`); the preset it matches is a separate, derived value.
  The former kind tags are still accepted by the builder, `set_rate_model` and `SetMassTransferKind`.
- A model created from the tag `series_resistance` has the class defaults (the lumped term only,
  k = 0.5 1/s); its terms no longer switch on when a parameter becomes non-zero. Set the switches instead
  (`resistances=` or the registry leaves). `SetMassTransferKind(..., series_resistance)` on an existing
  model keeps its switches, mode and values. Version-1 decks load with the switches their zero
  parameters meant, and their results do not change.
