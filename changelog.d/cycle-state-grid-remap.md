### Added — moving a cycle state to another grid
- `RemapCycleState` (`adsim/solver/cycle_state.h`) and
  `Simulation::RemapState` transfer a `CycleState` onto a deck with other
  cells per layer. Inside each layer, every per-cell field is remapped by
  overlap-weighted cell averages. The gas is remapped as component
  concentrations, so every component's gas and adsorbed inventory in every
  layer is kept to rounding. Connection histories and buffer-tank summaries
  are carried over unchanged. Refining by a whole-number ratio and coarsening
  back gives the original state.
- A warm start or `RunCycle` from a state on another grid now remaps it,
  where it was refused before. Different components, layers, void regions or
  step structure are still refused with the reason.
