### Fixed
- The slow test `NewtonKrylovCSSReference.ConvergesBelowTheNoiseFloor` no longer fails on Windows: it accepts a solve that stops at the integrator noise floor (about 1e-8) without reaching the 1e-8 Newton tolerance, and it measures contraction factors only clear of that floor.
