### Added — DAC breakthrough example and kinetics documentation
- Example deck `Illustrative_DAC_Breakthrough` (script `examples/illustrative_dac_breakthrough.py`): an
  isothermal 400 ppm CO₂ breakthrough on the Lewatit VP OC 1065 library adsorbent using the dual-rate LDF
  and the series-resistance coefficient together. An integration test checks its mole balances, the early
  front and slow tail, and that split = 1 with unit capacity factor and no micropore term reproduces the
  LDF + Haynes run bit for bit.
- `docs/kinetic-models.md` has a "Which to use for direct air capture" section and one notation for both
  kinetic sections; `docs/python.md` and the deck-builder kind list name `dual_rate_ldf`,
  `series_resistance` and `macropore_capacity`.
