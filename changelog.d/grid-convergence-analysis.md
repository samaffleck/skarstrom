### Added — grid convergence analysis in the engine
- New public header `adsim/solver/grid_convergence.h`: observed order of
  accuracy, Richardson extrapolation, Roache's grid convergence index and the
  cell count for a target error, for a result (a cycle KPI, say) on three or
  more grids at any constant refinement ratio. It reports whether the grids
  are in the asymptotic range: monotone, oscillatory (no order reported),
  divergent, grid independent or too few levels. The cycle-KPI grid
  convergence test now uses it (same numbers).
