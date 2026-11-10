### Added — valve opening ramp, and the case study's pressurisation valves open on it
- Valve streams take an opening ramp (`endpoint[ID].valve.opening_start`,
  `endpoint[ID].valve.opening_ramp_time`; Python `sk.ValveStream(..., opening_start=,
  opening_ramp_time=)` and `endpoint.valve.opening_start`; deck `"Opening_Start"`,
  `"Opening_Ramp_Time"`). The valve law's flow is multiplied by an opening fraction that rises
  from the start opening to 1 over the ramp time from the start of the step, along a cubic with a
  zero slope at both ends (no kink for the solver's Newton iteration). It applies to connected
  valves too. The defaults (1 and 0 s) give a valve that is fully open for the whole step, so
  existing decks keep their results; older decks load with the defaults (`Valve` version 2).
- The case-study decks `Perez2019_CaseStudy_Base` and `_Base_LPP` open their pressurisation valves
  on a ramp (variables `valve_ramp_start`, `valve_ramp_time`: 0.1 over 20 s for the feed-header
  valve, 0.2 over 8 s for both ends of the LPP valve). The four-step deck's approach to
  fluidisation at pressurisation falls from 3.6 to 0.95. The vacuum pump power and blowdown share
  were solved again, so the step-end pressures still match the validated decks to 2 %.
- `examples/perez_2019_case_study_optimise.py` searches the two ramp variables, and its
  `--summary DIR` option writes the fronts, the cost of each practical constraint and the
  comparison with the paper from a results file (committed under `examples/case_study_results/`).
  A new notebook, `examples/notebooks/case_study_vsa.ipynb`, runs the base case, reads the
  operating-limit and KPI reports, runs a small optimisation and plots the three fronts.
