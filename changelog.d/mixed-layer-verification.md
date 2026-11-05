### Added — verification of mixed adsorbent layers
- `core/test/validation/test_mixed_layers.cpp` (`MixedLayers.*`): solved tests for
  layers packed with several solids at set mass fractions. A blend of an adsorbent
  and an identical copy equals the pure layer, and a blend with an inert solid
  equals the pure adsorbent at the blend's density, voidage, conductivity and heat
  capacity, field by field in LTE and LTNE. The stoichiometric breakthrough time of
  a two-adsorbent layer on linear isotherms matches its mass-weighted capacity, and
  the mole and energy balances close on mixed beds. No defect was found.
