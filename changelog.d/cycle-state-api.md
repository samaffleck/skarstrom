### Added — cycle state: restart, warm start and one-cycle runs
- New public header `adsim/solver/cycle_state.h`: `CycleState`, the full state
  at a cycle boundary. It holds the DAE unknowns and their derivatives in the
  solver's order, a field map, and the differential entries as a flat vector
  (`Differential` / `SetDifferential`). It also holds every endpoint's logged
  history, which connected inlets replay, and the buffer-tank summaries. It
  records the deck structure it belongs to: components, layers with cell
  counts, steps, endpoints, step times and the number of cycles behind it.
- `Simulation::FinalState()` is the state at the end of the last completed
  run. `RunCycle(state)` runs exactly one cycle from a state, so `FinalState()`
  then gives the cycle map applied to it; `RunCycle()` with no state steps on
  from the last final state. `SetInitialState(state)` sets a warm start that
  later runs start from instead of the initial conditions. `CheckState(state)`
  gives the reason a state cannot start the deck: other components, layers or
  cell counts, another step structure, a changed state layout, or changed step
  times when connections or buffer tanks read the histories. A warm-started
  run counts the state's cycles towards the CSS minimum number of cycles.
- Python: `skarstrom.CycleState`; `Simulation.final_state()`, `run_cycle()`,
  `set_initial_state()`, `warm_start()`, `clear_initial_state()`,
  `initial_state` and `check_state()`.
- The state is kept for the session only and is never written to the deck
  file. Runs without a state are unchanged.
