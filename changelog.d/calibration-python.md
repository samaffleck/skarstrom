### Added — calibration against measured data from Python and the command line
- `skarstrom.calibrate(cases, parameters, ...)` and `Simulation.calibrate(parameters, datasets=None, ...)`
  fit deck parameters to the measured datasets stored in one or more decks, with the engine's
  calibration object (`adsim/solver/calibration.h`). A case is a Simulation, a
  `(Simulation, [dataset names], weight)` tuple or a `skarstrom.calibration.Case`; a parameter is a
  registry path (bounds from the registry) or a `skarstrom.calibration.Parameter` with bounds, a
  linear or log scale, a start and shared or per-case values. Bounds and starts take units like the
  rest of the package. Least squares (the default) or the global search, central or forward
  differences, a run budget, threads, a progress callback, Ctrl-C and a `CancelToken`.
- `CalibrationResult`: the estimates with standard errors and bound flags, the correlation matrix,
  the fit of each dataset (RMSE, NRMSE, maximum error, residuals, the simulated series), the
  objective trace, the runs, why it stopped, `table()`, `apply()`, `calibrated()`, `plot()` and a
  JSON round trip (`to_json()` / `from_json()`).
- `skarstrom calibrate deck.json --param PATH[:LO:HI[:log]] --dataset NAME[:WEIGHT] [--case OTHER.json]
  [--method least_squares|global] [--threads N] [--max-simulations N] [--out FILE] [--json FILE]`.
- The "Fitting parameters to measured data" guide (`docs/parameterisation.md`) is rewritten around
  the calibration API; `AGENTS.md` points agents to it.
