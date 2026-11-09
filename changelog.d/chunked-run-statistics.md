### Added — integrator statistics for chunked runs
- `CoupledSolverChuncked::GetStatistics()` returns the integrator counters of a run in slices (steps,
  residual evaluations, error-test and convergence failures, Jacobian evaluations and `recovered_steps`).
  They are the same counters, from the same code, that a blocking run reports in `RunResult::stats`.
