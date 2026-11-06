### Added — regression test for CSS acceleration on a thermal cycle
- A slow-tier test checks that Aitken acceleration shortens the approach to
  cyclic steady state on the thermal four-step LPP cycle (87 cycles without
  it, 39 with it) and that both runs stop at the same KPIs. The jump limit
  that replaced the `|mu| <= 5` cap already gives the reduction; earlier,
  both runs took the same number of cycles.
