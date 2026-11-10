### Added — run-health warnings
- Cheap checks after every cycle flag a run that completes but should not be trusted, as structured
  warnings `{code, value, threshold, cycle, step, suspects, hint}`: `MOLE_FRACTION_SUM` (the gas mole
  fractions drift from summing to 1), `MASS_BALANCE_STUCK` (a component's cycle mole balance stays open
  while the end-of-cycle state has settled: a model error, not slow CSS) and `STEP_RETRIES` (many rejected
  integrator attempts per accepted step). Thresholds leave every example deck and paper case silent
  with margin; the measured distribution is in the new docs page *Run health*.
- `RunResult.healthy` and `RunResult.health_warnings` (C++), `result.healthy` and `result.warnings`
  (Python), `Simulation::RunHealth()` / `sim.run_health()` with the per-cycle measures. The run status,
  exit code and results are unchanged. The solver logs the warnings once at the end of the run (LOW).
- CLI: the `--json` summary has `healthy` and `run_health` (structured; their messages are also appended
  to `warnings`), the text summary a `Health:` line, and run-health warnings go to stderr even with `-q`.
- Switches `settings.warn_mole_fraction_sum`, `warn_mass_balance_stuck` and `warn_step_retries` (on by
  default; `WarningSettings` class version 2), listed in the warning
  manifest with their `code`.
