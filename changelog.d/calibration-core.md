### Added — calibration against measured data (C++)
- `adsim/solver/calibration.h`: fit model parameters to the datasets stored in one or more decks.
  A problem is a set of cases (a deck, the bound datasets to reproduce, a run policy) and a set of
  parameters (registry paths with bounds, a linear or log scale, and shared across cases or one
  value per case), validated against every deck.
- The objective is a weighted least-squares sum over every measured point, scaled by the point's
  uncertainty or the dataset's range and by the dataset's weight per point. The default solver is a
  bounded Levenberg–Marquardt method with a finite-difference Jacobian whose columns run in
  parallel threads (warm-started from cached cyclic steady states where the deck has them); the
  global MaxLIPO search is also available, optionally refined by least squares. Progress and cancel
  callbacks; a failed run is counted, never thrown.
- The result gives each estimate with its standard error and whether it is at a bound, the
  correlation matrix, per-dataset RMSE, NRMSE, maximum error, residuals and simulated series, the
  objective trace, the runs made, why the solver stopped, and a calibrated deck per case.
- `RunFromState` (`adsim/solver/state_cache.h`): the run half of `RunWithStateCache`, for callers
  that choose the warm start and cache the final state themselves.
