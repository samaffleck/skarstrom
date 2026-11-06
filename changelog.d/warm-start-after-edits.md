### Added — warm starts after an edit
- `RunOptions::initial_state` / Python `run(initial_state=state)` start one
  run from a previous run's state, for example its `final_state()`. A state
  set with `SetInitialState` is used for every later run.
- A state now starts a run after the step times change, where it was refused
  before. Each changed step's logged histories are stretched onto the new
  duration, and its buffer-tank summary moles are scaled with it
  (`ResampleCycleStateHistories`).
- `PrepareCycleState` / `Simulation::PrepareState` / Python
  `prepare_state()` return the state as a run would start from it (remapped
  onto other cells per layer, histories stretched) with a note of what was
  done. The run logs the same note. Other components, layers, void regions or
  step structure are still refused with the reason.
- Fewer cycles to the new CSS after typical edits (a step time, a feed
  composition, the cell count), with the same purity and recovery as a cold
  start. The table is in `docs/cyclic-steady-state.md`.

### Fixed — grid remap overshot the pressure at temperature fronts
- `RemapCycleState` now remaps the fluid temperature weighted by the gas
  concentration. The remapped pressure is then the plain cell average of the
  old pressure. Before, it could overshoot both neighbours where the
  temperature changed sharply (by 20 Pa on the TSA example, 30 → 45 cells),
  and the consistent-initial-condition solve then failed.
