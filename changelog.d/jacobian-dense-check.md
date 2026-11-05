### Added — dense finite-difference check of the solver Jacobian
- New integration test `JacobianDenseCheck`: at states in the middle of each
  step of real decks, the Jacobian IDA receives (coloured, and the band
  evaluations that learn and audit its pattern) is compared entry by entry
  with a dense forward difference over the whole matrix, with no colouring
  and no band assumption. It covers upwind, WENO and van Leer, thermal and
  isothermal runs, Kozeny-Carman and Ergun momentum, a wall energy balance,
  a breakthrough, a PSA and a TSA cycle, and each deck's grid as well as 10
  and 20 cells. Every entry agrees to rounding (relative error about 2e-16);
  a guard case shows that the comparison detects an incomplete pattern.
- `BandedJacobian::SetEvaluationObserverForTesting`: a read-only test hook
  that sees every assembled Jacobian and its column increments.
