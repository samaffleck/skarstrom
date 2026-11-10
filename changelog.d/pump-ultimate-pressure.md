### Added — vacuum pump ultimate pressure and step-end pressures
- A vacuum pump stream (a *Pump* device on an outlet) takes an ultimate (blank-off) pressure,
  `endpoint[N].pump.ultimate_pressure` [Pa] (`sk.Pump(..., ultimate_pressure=...)` on a `sk.PumpStream`).
  The flow from the power is multiplied by the rotary-pump speed factor `1 - p_ult/P1`, rounded to a C1
  zero at `p_ult`, so the bed pressure approaches the ultimate pressure from above and is held there.
  The default, 0, is off: no shipped deck changes its results. `PumpBoundaryStream` class version 2.
- The pressure each step ends at: `Simulation::StepEndPressures(cycle)` / `sim.step_end_pressures(cycle)`,
  one record per step with the bed-end pressures (`bottom`, `top`, the last value of the step's bed-end
  series) and `bed_min`, `bed_mean`, `bed_max` over the packed cells. Recorded with the KPIs when each
  cycle completes, so every cycle of a run is available; the CLI's `--json` output lists the last cycle's
  as `step_end_pressures`.
- Sensitivity studies and the optimiser read them through a step reference in place of a component:
  outputs `step[<ID or name>].end_pressure_bottom` (`_top`, `_min`, `_mean`, `_max`;
  `sa.step_output_names()`), and objective and constraint terms such as
  `opt.at_least("step[Evacuation]", "end_pressure_bottom", 1e4)` (C++ `KPIType::END_PRESSURE_*`, appended
  to the enum).
