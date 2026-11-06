### Added — internal heat exchanger
- A step's heating medium can also heat or cool the bed directly, through an
  internal heat exchanger: `a_x h_x (T_m - T)` per bed volume in every bed
  cell, at the same medium temperature as the jacket. Under LTNE the gas
  takes the `eb` share and the solid the `1 - eb` share. It is off by default
  (`a_x = 0`).
- Registry parameters `step[ID].medium.exchanger_area_per_volume` (m2/m3)
  and `step[ID].medium.exchanger_heat_transfer_coefficient` (W/(m2 K)).
- The exchanger's heat is part of the logged medium duty
  (`StepHeatDuty::medium_heat_in` / `medium_heat_out`), and so of the
  thermal-energy and cooling-duty KPIs. Its share is also logged on its own
  (`exchanger_heat_in` / `exchanger_heat_out`).
- `HeatingMedium` is at class version 1. Older decks load with no exchanger
  and give the same results. Re-saved decks gain the two exchanger fields.
