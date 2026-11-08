### Changed — one mass-transfer coefficient model with a switch per resistance
- The `series_resistance` coefficient (`MassTransferCoefficientSeriesResistance`, cereal version 2) is now
  the one coefficient model: external film, macropore (with a Knudsen sub-switch), micropore and a lumped
  fitted coefficient (k with an optional activation energy), each switched on or off explicitly and added in
  series, and an explicit capacity-factor mode (`none`, `slope`, `secant`, `fixed`). Every parameter is kept
  whatever its switch, so switching a term off and on restores it.
- The older kinds (`constant`, `arrhenius`, `macropore`, `macropore_knudsen`, `haynes`, `macropore_capacity`)
  are presets of it. Decks and library entries that store them still load (the old classes are kept as
  load-only types) and are converted on load; each preset is evaluated in the old class's operation order,
  so results are bit for bit those of the old kind. Every old tag is still accepted by the builder,
  `set_rate_model` and `SetMassTransferKind`, and the reported kind is the preset the switches match.
  Version-1 `series_resistance` archives load with the switches their zero sentinels meant.
- New deck edits `SetMassTransferResistance` (one switch, fast or slow site; the last term cannot be
  switched off, Knudsen needs macropore) and `SetCapacityFactorMode`. Each dual-rate site has its own
  switches.
- New deck checks: no resistance on, or an enabled term with a parameter out of range, is an error; a film
  term with a zero film resistance factor, or a second-order LDF with a physical resistance, is a run-start
  warning.
