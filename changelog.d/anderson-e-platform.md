### Fixed — Anderson deck test failed on macOS arm64
- `CSSAndersonDecks` (slow tier) now holds the Anderson stop to the plain stop's KPI distance from
  the reference (the larger of the purity and recovery errors), not to each KPI. The stop rule
  bounds the bed state, not each KPI, and the split between purity and recovery at an Anderson stop
  changes with the compiler and CPU. On macOS arm64 the Perez2019_E case stopped at cycle 42 with a
  recovery error of 5.8e-5 (allowance 5e-5) while its distance stayed inside the plain stop's 6.6e-5.
  No engine change.
