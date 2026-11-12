### Added — calibration examples
- `examples/calibrate_ward_2022.py` fits the three Ward 2022 breakthrough runs together with
  `skarstrom.calibrate` (shared q_s,2 and dispersion factor, one inside heat transfer coefficient
  per flow, the hand calibration's weights) and reproduces the calibration of the shipped decks;
  `--quick` runs it in seconds. `examples/notebooks/calibrate_ward_2022.ipynb` goes through it
  step by step with the overlay and residual plots.
- `examples/calibrate_cycle_kpis.py` fits two parameters of `Perez2019_A` to per-cycle CO₂ purity
  and recovery (datasets bound to KPI observables): a twin experiment that recovers the values that
  made the data.
