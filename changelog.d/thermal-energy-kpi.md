### Added — thermal energy and cooling duty KPIs
- Two KPIs from the step energy sources, per kg of net product (J/kg),
  defined and undefined exactly as the electrical energy:
  `thermal_energy`, the heat supplied over the cycle (sources with a
  positive power), and `cooling_duty`, the heat removed (negative power),
  reported as a positive magnitude and not subtracted from the thermal
  energy. Each `KPIRecord` also has `step_thermal_energy` and
  `step_cooling_duty`, every step's share in the deck's step order. They are
  in `KPIRecord` (C++ and Python), `as_dict()`, `kpi_table()`, the CLI's
  `--json` and its printed table, two columns appended to `kpis.csv`, the
  sensitivity outputs (`Component.thermal_energy`, `Component.cooling_duty`)
  and the `KPIType` / `KPIVariableType` enumerators `THERMAL_ENERGY` and
  `COOLING_DUTY`, appended so no stored index moves.
- They use the power the solver applied in each step, logged once per step
  with the step's log (`Log::heat_duty`, evicted with its cycle): a disabled
  source, or one left on a removed layer, counts nothing, and both are 0 in an
  isothermal run. Results do not change.
