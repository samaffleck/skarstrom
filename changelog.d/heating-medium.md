### Added — heating medium on the wall's outer face
- Each step can heat or cool the column through a jacket: a heating medium
  at a set temperature (constant, or a ramp with the inlet-temperature
  profiles) and a heat transfer coefficient on the wall's outer face. While
  the step runs they replace the ambient temperature and every layer's
  external heat transfer coefficient; steps without a medium keep the
  ambient exchange. The medium temperature is set, not solved.
- Registry parameters `step[ID].medium.enabled`, `.temperature` (or
  `.temperature_linear_start`, `_end`, `_time`, ... for a ramp) and
  `.heat_transfer_coefficient`, listed for every step, so a medium can be
  switched on from the CLI (`--set`), swept, optimised and used in
  sensitivity studies; in Python, `sim.steps["Heat"].medium`.
- `Step` is at class version 6 and carries a new `HeatingMedium` (class
  version 0). Older decks load with no medium and give the same results;
  re-saved decks gain a disabled `Heating_Medium` in every step.
