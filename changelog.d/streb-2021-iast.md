### Fixed — IAST with steep Sips terms and vanishing mole fractions
- The IAST mixing rule failed at the start of a breakthrough into a bed of inert gas. Two causes, both fixed in the IAST solve:
  - The integrator's decaying front toes reach subnormal mole fractions (about 5 × 10⁻³²¹). There the spreading pressure underflowed and the solve returned NaN loadings. Mole fractions at or below 10⁻¹⁵⁰ now take IAST's first-order dilute limit.
  - A Sips-family term with exponent n < 1 has an infinite Henry slope, so the dilute loading of an absent component diverged (−1.3 × 10⁷ mol/kg at y = −10⁻¹⁴ in the Streb 2021 CO₂/CH₄ decks). Inside the IAST solve, such a term is now continued below 10⁻⁴ of the total pressure by a monotone C¹ form whose spreading pressure is exact there. So IAST is unchanged for every component with a mole fraction of 10⁻⁴ or more.
- No bundled deck uses IAST, so no example deck changes.

### Added — validation: Streb (2021) under IAST
- Every Streb 2021 breakthrough case is also run under IAST over the pure Sips terms of the same parameters, and recorded without assertion (`Experiment<case>_IAST` on the Validation page). This is the competition model the paper prefers.
- IAST moves every CO₂ front 0.4–14 % later and every CH₄ front 1–3 % earlier. It is better on 14 of the 47 asserted metrics and worse on 33. Extended Sips therefore stays the decks' model, and no bound changes. The suite's README has the table and the attribution.
