### Added — series-resistance LDF coefficient
- A **series-resistance (Glueckauf) LDF coefficient** (`MassTransferCoefficientSeriesResistance`,
  kind `series_resistance`): 1/k = f K r_p/(3 k_f) + K τ r_p²/(15 ε_p D_m) + r_c²/(15 D_c), external
  film (k_f from the layer's film correlation, as `haynes`), macropore and micropore resistances in
  series, with the local capacity factor K = ρ_p ∂q*/∂c from the isotherm slope (secant where the
  slope is not positive, floored at 1; or a fixed `capacity_factor`) and an Arrhenius micropore
  diffusivity D_c = D_c0 exp(−E_a/RT). Each term is switched off by its parameter. With the
  micropore term off and K = 1 it equals `haynes` bit for bit. Defined once in `model_kinds`, so the
  registry (`kinetics[C].crystal_radius`, `.micropore_diffusivity`, ...), the builder and Python
  `model_kinds()` pick it up. Existing decks are unchanged.
