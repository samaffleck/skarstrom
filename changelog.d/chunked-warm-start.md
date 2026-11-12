### Added — warm starts for the chunked run
- `Simulation::Start` (the chunked run driven with `Advance`) can start from a cycle state: a new `initial_state` argument, or the `SetInitialState` state. As for a blocking run, the state is adapted to the deck first, its cycles count towards the CSS minimum number of cycles, and a state that does not fit makes `Start` return false with the reason.
- A chunked run that completes or reaches cyclic steady state now keeps its `FinalState()`, so a caller that steps the run itself can use the state cache (`adsim/solver/state_cache.h` describes the pattern).
- A run without a start state gives the same results as before.
