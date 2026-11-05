### Fixed — extended Langmuir–Freundlich trace loadings
- `extended_langmuir_freundlich` no longer multiplies each competitive term by
  the smooth step y² / (y² + 1e-10) (with a softplus on y). That step biased
  trace loadings by about 1e-10 / y²: 1e-4 at y = 1e-3, 1 % at 1e-4, about
  50 % at 1e-5. The power (P y / P_ref)^n is now continued below y = 1e-5 by
  the same C1 quadratic the competitive Sips kinds use, so the isotherm is
  exact for y ≥ 1e-5 in both the solver and the single-component evaluation,
  with a finite Henry slope below and exactly zero loading at y = 0. Results
  of runs with trace components on this kind move; no bundled deck uses it.
