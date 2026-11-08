### Changed — recovered integrator failures are quiet and counted
- SUNDIALS error messages go to the engine log instead of stderr: the solver replaces the SUNDIALS
  context's default error handler, so nothing from IDA is printed directly any more. Every IDA message
  appears in the log at level *high* (`IDA: ...`).
- A step that IDA fails but the engine's retry recovers is logged at *medium* as a "recovered step" with
  IDA's message (the chunked solver now logs these too). Only a step whose retries run out is logged at
  *low*, with IDA's last message. Results and retry behaviour are unchanged.

### Added — `recovered_steps` run statistic
- `IntegrationStats::recovered_steps` (C++ `RunResult.stats`), `result.stats["recovered_steps"]` (Python)
  and `stats.recovered_steps` in the CLI `--json` summary and `run.json`: the integrator failures the
  engine recovered by retrying. The run-health `STEP_RETRIES` measure reads the same count.
- Docs: starting a flow inlet into a bed at low pressure (ramp it over about 1 s) and what a recovered step
  means in the log (*Boundary conditions*).
