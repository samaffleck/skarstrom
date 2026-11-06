### Added — heating-medium duty
- Each step's log records the heat its heating medium passed into the wall
  and took back out over the step (`StepHeatDuty::medium_heat_in` and
  `medium_heat_out`, J). It is integrated with the trapezoidal rule over
  every logged point, separately from the energy-source powers.
- The thermal-energy and cooling-duty KPIs, total and per step, now include
  the medium's heat in and heat out. `CalculateStepThermalDuty` also reports
  the medium's share (`medium_heating`, `medium_cooling`). Results do not
  change for decks without a medium.
