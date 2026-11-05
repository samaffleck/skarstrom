### Added — Skarstrom-cycle equilibrium-theory verification
- `test_skarstrom_cycle_equilibrium_theory.cpp` checks a full four-step
  Skarstrom cycle (pressurisation with light product, feed, blowdown, purge)
  on linear isotherms against the Knaebel-Hill equilibrium theory: light-product
  purity 1, the closed-form recovery at two purge-to-feed ratios and the
  heavy-product composition. A slow-tier test refines the grid and the LDF
  rate and shows the gap falling at first order and extrapolating to the theory.
