### Changed — faster test builds: a precompiled header for the test target
- `AdsorptionSimulatorTests` gets its own precompiled header (the engine's
  `core/src/pch.h` plus GoogleTest), so its ~130 sources no longer parse
  cereal, Eigen, CoolProp, SUNDIALS and GoogleTest from scratch.
  `_USE_MATH_DEFINES` becomes a definition on the test target, since a
  `#define` above `<cmath>` in a source file no longer reaches a header that
  the precompiled header already includes.
