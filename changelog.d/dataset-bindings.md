### Added — experimental datasets record what they measure
- A dataset stored in a deck can now name the simulated quantity it measures (its observable): a
  bulk or component variable at the top, the bottom or a position, over one cycle, or a per-cycle
  KPI of a tracked component. The binding also gives the units the points are in (a time unit for
  x; any unit of the observable's dimension for y, offsets such as degC included), an optional
  per-point uncertainty and a weight. Decks saved before load with their datasets unbound.
- Python: `sk.Observable`; `sim.experimental_data.add(..., observable=, x_unit=, y_unit=,
  uncertainty=, weight=)` and the same in `load_csv` (plus `uncertainty_column=`); writable
  `dataset.observable`, `x_unit`, `y_unit`, `uncertainty`, `weight`; `dataset.simulated(sim)` (the
  simulated counterpart in the dataset's units) and `dataset.compare(sim)` (RMSE, NRMSE, maximum
  error, points compared and the weighted sum of squares). `compare(x, y)` is unchanged.
- C++: `XYData::observable` (`adsim/solver/experimental_data.h`), `SimulatedSeries` and
  `CompareWithSimulation` (`adsim/solver/parameterisation.h`), and the edits `AddExperimentalData`,
  `SetExperimentalDataBinding` and `RemoveExperimentalData` (`adsim/objects/deck_edit.h`, also
  `sim.edit` operations).
- The deck view has a root `experimental_data` list with each dataset's points and binding.
- The bundled decks with data (Ward 2022, Streb 2021, Nguyen 2022, the diluted-bed example) are bound.
