# Changelog

All notable changes to Skarstrom are recorded here. The format follows
[Keep a Changelog](https://keepachangelog.com/en/1.1.0/) and the project uses
[Semantic Versioning](https://semver.org/). The single source of the version is
`ADSIM_PROJECT_VERSION` (plus `ADSIM_VERSION_SUFFIX` for pre-releases) in the
root `CMakeLists.txt`.

## [Unreleased]

## [0.12.0] — 2026-11-13

### Fixed — per-component CSS test with an absent component
- The Mass Balance (Per-Component) CSS test now skips a component whose inflow over the cycle is below 10⁻¹² times the moles of all components fed in the cycle (before: below 10⁻³⁰ mol). Under the IAST mixing rule, a component that is in the deck but never fed takes round-off loadings, so its flows were round-off-sized but not zero, its relative balance error was O(1) on about half the cycles, and runs reached CSS only by chance or went to the cycle cap. Such a run now stops where the same deck without the component stops. The solver logs each skipped component once per run ("negligible throughput, mass balance not tested"). Fed components are tested exactly as before.

### Added — water on CALF-20: a type V water isotherm and water factors on the dual-site Langmuir kinds
- New isotherm kind `cubic_langmuir`: the type V (Cubic-Langmuir) water isotherm of Nguyen et al., Ind. Eng. Chem. Res. 63 (2024) 3265–3281, in the gas concentration, with one energy for the three cooperative coefficients and one for the Langmuir term. The defaults are the paper's adsorption-branch values for CALF-20. It has a closed IAST spreading pressure, a fitter set-up, and a deck check that its denominator stays positive up to saturation at the deck's temperatures. The loading is held at 0 where a negative first coefficient would make it negative.
- The dual-site Langmuir kinds (`dual_site_langmuir_p`/`_c` and their extended forms) take four optional water factors: per-site capacity factors 1/(1 + λ_k q_w^λ3) (`lambda1`, `lambda2`, `lambda3`) and an affinity factor exp(−γ2 q_w) (`gamma2`), with q_w the adsorbed water loading of the same adsorbent. All default to 0, which is the dry isotherm exactly. The fields are written to a deck only when set, so existing decks and their results do not change. On an extended kind, the other components keep the shared denominator. With a factor set, the IAST mixing rule rejects the isotherm. `skarstrom.evaluate_isotherm(..., water_loading=)` and `fit_isotherm(..., water_loading=)` support the factors.

### Added — molecular diffusion has its own switch in the macropore term
- The resistance model has a new switch, `molecular` (registry leaf
  `kinetics[C].resistances.molecular`), that refines the macropore term next to
  `knudsen`. The pore diffusivity is D_m (molecular only), D_K (Knudsen only) or
  the Bosanquet mean 1/D = 1/D_K + 1/D_m (both). One of the two stays on while
  macropore is on; switching macropore on with neither on switches molecular on.
- `molecular_diffusion_coefficient` belongs to the `molecular` switch.
- `MassTransferResistances` is at cereal version 2. A deck from an earlier version
  loads with `molecular` on wherever `macropore` is on, so its results do not change.
- Python: `set_rate_model(resistances=[...])` keeps molecular diffusion on with
  macropore; switch it off with `set_rate_resistance(component, "molecular", False)`.
  The deck-to-script writer does this for a Knudsen-only site.

### Changed — the coefficient's display equation has one definition per line
- `MassTransferEquationLatex` puts each definition on its own line, aligned on
  the equals sign (it joined them on one line before). The Bosanquet mean and
  D_K are separate lines.
- With the film term on in a packed layer, the deck view's `equation` also shows
  the Sherwood number of the layer's film correlation (Wakao–Funazkri, Gunn or
  Ranz–Marshall) and the Reynolds and Schmidt numbers it uses
  (`FilmCorrelationLatex`).
- The model-kind manifest's resistances have `needs_child` (true for macropore).

### Added — warm starts for the chunked run
- `Simulation::Start` (the chunked run driven with `Advance`) can start from a cycle state: a new `initial_state` argument, or the `SetInitialState` state. As for a blocking run, the state is adapted to the deck first, its cycles count towards the CSS minimum number of cycles, and a state that does not fit makes `Start` return false with the reason.
- A chunked run that completes or reaches cyclic steady state now keeps its `FinalState()`, so a caller that steps the run itself can use the state cache (`adsim/solver/state_cache.h` describes the pattern).
- A run without a start state gives the same results as before.

### Fixed — Nguyen 2022: compare with the measured data, not the paper's model
- The files in `core/data/nguyen_2022/Fig4/` and `Fig5/` were the paper's model curves, digitised; the paper test and the `Nguyen2022_Fig4` and `Nguyen2022_Fig5` decks used them as the experiment. They are now `core/data/nguyen_2022/paper_model/` (the test's "vs paper model" check). The measured markers are re-digitised from the paper's figures (`core/data/nguyen_2022/digitise_figures.py`), including the outlet N2 of Fig. 5.
- The decks follow the paper and the author's thesis (Table A2.1): Fig. 5 feeds CO2 in nitrogen (the decks fed CO2 in helium), each run has its own temperature and pressure, and the decks use the constant gas heat capacities of Table 2, the capacity-factor LDF of Eq. 5 with D_m = 1.496e-5 m2/s, Eq. 6 dispersion with the same D_m, the Blake–Kozeny pressure drop and no gas holdup in the pores.
- The measured thermocouple 1.82 cm from the outlet is now compared too.
- Nothing is fitted. Against the measurement, the CO2 fronts are within 3–11 % (the paper's own model: 5–15 %). The Fig. 5 thermocouple RMSE is 1.0–2.7 K. The Fig. 4 temperature peaks are 6–8 K too high; this is an open item. `core/test/papers/nguyen_2022/README.md` has the provenance, the attribution and the open items.

### Added — Nguyen 2022: all six breakthrough runs as example decks
- `Nguyen2022_Fig4_50pct`, `Nguyen2022_Fig4_100pct`, `Nguyen2022_Fig5_5pct` and `Nguyen2022_Fig5_30pct` join `Nguyen2022_Fig4` (70 % CO2) and `Nguyen2022_Fig5` (15 % CO2). `examples/nguyen_2022.py` builds all six.
- The test cases are renamed `Nguyen2022Test.Fig4_100pct` … `Fig5_5pct`.

### Fixed — registry coverage test across map entries
- `RegistryCoverage` intersects its three passes per map entry (one component's value at a time), not per generic path. A deck whose adsorbent mixed a capacity-factor component with a lumped one was reported as uncovered, although each value is covered in one of the passes.

### Added — Streb 2021: all nine breakthrough experiments as example decks
- `Streb2021_A2`, `Streb2021_A4` … `Streb2021_A7`, `Streb2021_B1` and `Streb2021_B2` join `Streb2021_A1` and `Streb2021_A3`. Each deck is the paper test's case: the calibration on A3, the experiment's flow factor and its digitised data, bound to the deck.
- `examples/streb_2021.py` builds all nine, including the ternary H2/CO2/CH4 cases. It runs A1 and A3 by default; `--case` selects the others.

### Added — calibration examples
- `examples/calibrate_ward_2022.py` fits the three Ward 2022 breakthrough runs together with
  `skarstrom.calibrate` (shared q_s,2 and dispersion factor, one inside heat transfer coefficient
  per flow, the hand calibration's weights) and reproduces the calibration of the shipped decks;
  `--quick` runs it in seconds. `examples/notebooks/calibrate_ward_2022.ipynb` goes through it
  step by step with the overlay and residual plots.
- `examples/calibrate_cycle_kpis.py` fits two parameters of `Perez2019_A` to per-cycle CO₂ purity
  and recovery (datasets bound to KPI observables): a twin experiment that recovers the values that
  made the data.

### Changed — Streb 2021 decks calibrated on experiment A3
- The `Streb2021_A1` and `Streb2021_A3` decks, and the paper test's nine cases, are calibrated on the reference experiment A3 alone:
  - the CO2 and CH4 LDF coefficients: 0.0592 and 0.137 1/s (Table 4 gives 0.06 and 0.2, also fitted to A3);
  - the wall-to-thermofluid coefficient: 49.7 W/(m² K) (the paper takes 220 from a reference);
  - the sampling line to the mass spectrometer: 0.142 L of gas at column pressure. It is an inert, nearly empty layer above the bed (bed voidage 0.99, 40 cells), so the line carries the outlet gas as plug flow. A void region is one well-mixed cell and smears the CH4 front.
- Each case's feed flow is its Table 2 set point times a flow factor, fitted to the speed of its own CO2 temperature front between the thermocouples (A1 1.051, A3 1.020). The paper also fits the feed velocity per experiment. Six of the nine factors are more than 2 % from the set point (up to −11.3 % for B2).
- The other eight cases are predictions. Their CO2 fronts are now within 1.4 % of the measurement (before: 0.3–11 %) and their CH4 fronts within 2.4 %, except B1 at 4.7 % (before: 0.2–18 %). The thermocouple RMSE is now 5.5–11.4 K (before: 8.7–28.5 K). `core/test/papers/streb_2021/README.md`, "Calibration on A3", has the fit, the flow factors, the result held to 2 % of the set point and the limits. `streb_calibration.py` in the same folder repeats the calibration.
- `CellVariableAt` (test fixture `validation_metrics.h`) takes the log index of the layer's first cell, for beds of more than one layer.

### Added — calibration against measured data from Python and the command line
- `skarstrom.calibrate(cases, parameters, ...)` and `Simulation.calibrate(parameters, datasets=None, ...)`
  fit deck parameters to the measured datasets stored in one or more decks, with the engine's
  calibration object (`adsim/solver/calibration.h`). A case is a Simulation, a
  `(Simulation, [dataset names], weight)` tuple or a `skarstrom.calibration.Case`; a parameter is a
  registry path (bounds from the registry) or a `skarstrom.calibration.Parameter` with bounds, a
  linear or log scale, a start and shared or per-case values. Bounds and starts take units like the
  rest of the package. Least squares (the default) or the global search, central or forward
  differences, a run budget, threads, a progress callback, Ctrl-C and a `CancelToken`.
- `CalibrationResult`: the estimates with standard errors and bound flags, the correlation matrix,
  the fit of each dataset (RMSE, NRMSE, maximum error, residuals, the simulated series), the
  objective trace, the runs, why it stopped, `table()`, `apply()`, `calibrated()`, `plot()` and a
  JSON round trip (`to_json()` / `from_json()`).
- `skarstrom calibrate deck.json --param PATH[:LO:HI[:log]] --dataset NAME[:WEIGHT] [--case OTHER.json]
  [--method least_squares|global] [--threads N] [--max-simulations N] [--out FILE] [--json FILE]`.
- The "Fitting parameters to measured data" guide (`docs/parameterisation.md`) is rewritten around
  the calibration API; `AGENTS.md` points agents to it.

### Added — calibration against measured data (C++)
- `adsim/solver/calibration.h`: fit model parameters to the datasets stored in one or more decks.
  A problem is a set of cases (a deck, the bound datasets to reproduce, a run policy) and a set of
  parameters (registry paths with bounds, a linear or log scale, and shared across cases or one
  value per case), validated against every deck.
- The objective is a weighted least-squares sum over every measured point, scaled by the point's
  uncertainty or the dataset's range and by the dataset's weight per point. The default solver is a
  bounded Levenberg–Marquardt method with a finite-difference Jacobian whose columns run in
  parallel threads (warm-started from cached cyclic steady states where the deck has them); the
  global MaxLIPO search is also available, optionally refined by least squares. Progress and cancel
  callbacks; a failed run is counted, never thrown.
- The result gives each estimate with its standard error and whether it is at a bound, the
  correlation matrix, per-dataset RMSE, NRMSE, maximum error, residuals and simulated series, the
  objective trace, the runs made, why the solver stopped, and a calibrated deck per case.
- `RunFromState` (`adsim/solver/state_cache.h`): the run half of `RunWithStateCache`, for callers
  that choose the warm start and cache the final state themselves.

## [0.11.0] — 2026-11-12

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

## [0.10.0] — 2026-11-12

### Changed — Ward 2022 decks calibrated against the experiment
- The `Ward2022_Q50ccm`, `Ward2022_Q100ccm` and `Ward2022_Q200ccm` decks are calibrated against the digitised outlet and thermocouple curves of all three flows together. The calibrated values are the ambient and wall-end temperature (293 K; Table 1 gives 291.4 K, but the measured thermocouple settles at 293 K after the front), the site-2 capacity q_s,2 (5.89 mol/kg; Table 2 gives 6.04, inside the paper's 95 % interval), the factor on the paper's dispersion correlation (3.75; the paper fits 4) and the inside heat transfer coefficient (97.3 / 81.5 / 69.0 W/(m² K); the paper fits 77.5 / 57.9 / 45.2). The decks also now use Table 1's particle density, ρ_b / (1 − ε) = 799.6 kg/m³ (799 before).
- The adsorbent library entry `activated_carbon_ward2022` keeps the published Table 2 isotherm and is no longer deck-checked against `Ward2022_Q50ccm`.
- The CO2 front time (y = 0.5) is now within 0.4 % of the measurement at all three flows (it was 1.7-4.5 % late), the outlet-curve NRMSE is 2.7-3.9 % (it was 14-35 %) and the thermocouple RMSE 0.18-1.40 K (it was 1.15-2.56 K). See `core/test/papers/ward_2022/README.md`, "Calibration". `examples/ward_2022.py` builds the calibrated decks.

### Changed — heat capacity seeds
- A component's constant gas heat capacity (`fluid.gas_heat_capacity[<component>]`, used in the `constant` mode) now starts at its ideal-gas heat capacity at 25 °C, from its NASA polynomial at 298.15 K rounded to 0.01 J/(mol·K), for example 37.14 for CO2 and 29.12 for N2. It was 29.1 J/(mol·K) (7/2 R) for every component. `Fluid::DefaultGasHeatCapacity(component)` gives the seed. The example decks in the `ideal_gas` mode carry the new seeds. Their results do not change, because that mode does not use the constants.
- A component added with `add_component` (`AddComponentToDeck`, which the web workspace also uses) now starts with the same value as its adsorbed-phase heat capacity (`fluid.adsorbed_heat_capacity[<component>]`). It was 0. In the default `constant` adsorbed mode this changes the results of a new deck that does not set the value. A deck that has no value for a component still uses 0, so a saved deck gives the same results as before. The `arora_2021` and `illustrative_dac_breakthrough` example scripts now set their 0 values explicitly.
- Both values stay editable.

### Removed — six adsorbent library entries
- The adsorbent library no longer ships these entries: `cms_3k` (Carbon molecular sieve, CMS-3K, source unknown), `aluminium_ho1972`, `ceramic_fibre_paper_headley2019`, `cordierite_robie1995`, `stainless_steel_304l_kim1975` and `stainless_steel_316l_kim1975`. Their blobs, their library records, the materials `aluminium`, `ceramic-fibre-paper`, `cordierite`, `stainless-steel-304l` and `stainless-steel-316l`, and the sources only they cited are gone. `sim.add_library_adsorbent("cordierite")` (and the others) now fails with an unknown-material error. A deck that already holds one of these adsorbents keeps it, because a deck stores its adsorbents in full.
- The library now has no `support` material. The `support` family and its validation rules stay, so a new support solid can be added later. The `test_support_solids` integration test is removed, because it had no entries to run.
- `cms_3kt_zhang2021` is now the only carbon molecular sieve entry. Its label is now `Carbon molecular sieve`, and it has no variant. Its blob and its isotherm data do not change.
- The wall-material presets (`ss316l`, `ss304l`, `aluminium_6061` and the others) are not library entries and do not change.

### Fixed — Anderson deck test failed on macOS arm64
- `CSSAndersonDecks` (slow tier) now holds the Anderson stop to the plain stop's KPI distance from
  the reference (the larger of the purity and recovery errors), not to each KPI. The stop rule
  bounds the bed state, not each KPI, and the split between purity and recovery at an Anderson stop
  changes with the compiler and CPU. On macOS arm64 the Perez2019_E case stopped at cycle 42 with a
  recovery error of 5.8e-5 (allowance 5e-5) while its distance stayed inside the plain stop's 6.6e-5.
  No engine change.

### Added — Nguyen et al. (2025) dry VSA cases on CALF-20
- A paper validation suite `core/test/papers/nguyen_2024/` (`Nguyen2024Test.CaseA0`, `CaseB0`) for the two dry cases of Table 1 of Nguyen, Shimizu & Rajendran, Adsorption 31 (2025) 45 (doi:10.1007/s10450-025-00600-z). The decks are the `Nguyen2023_Exp_VII` four-step light-product-pressurisation deck at this paper's conditions (0.96 bar, 22 °C, each case's adsorption time, feed velocity and blowdown and evacuation pressures); the light-product connection is unchanged. The suite asserts the purity of both cases and the productivity of B0, and records the recovery, which the engine over-predicts (README).
- Two example decks, `Nguyen2024_A0` and `Nguyen2024_B0`, saved by the suite.
- The humid cases of the paper are not reproduced: the paper prints no water isotherm for CALF-20, and the engine has no competitive CO₂/H₂O kind for a physisorbent.

### Removed — CoolProp from the build
- The engine no longer builds or links CoolProp: the `external/CoolProp` submodule, its `add_subdirectory` and cache settings, the `FMT_CONSTEVAL=` definition, the CoolProp link, include and install rules, and the `Adsim::CoolProp` imported target of the installed CMake package are gone. The component constants come from the generated table `component_properties.h` (since the previous release), so no result changes. A clean Release build of `adsim` (GCC 13, 4 jobs) drops from 3 min 50 s / 881 CPU-s to 2 min 45 s / 632 CPU-s, and a WebAssembly build no longer compiles CoolProp's fluid data and backends into the module.
- `CMAKE_POLICY_VERSION_MINIMUM` is no longer set: it was there for CoolProp's `cmake_minimum_required` calls, and every remaining dependency configures with CMake 4 without it.
- For C++ consumers: `Adsim::adsim` no longer carries `FMT_CONSTEVAL=` or the CoolProp include path, and the install tree has no `libCoolProp.a`, CoolProp headers or fmt. A program that used CoolProp through the engine must now bring its own.
- `scripts/gen_component_properties.py` reads the CoolProp fluid files it is pinned to (CoolProp 6.7.0, commit `e84ec7a8ab16`, each file by SHA-256): it downloads them, or reads a local `dev/fluids` directory with `--fluids`.

### Changed — component constants come from a generated table
- The engine no longer calls CoolProp at run time. The molar mass, formula, CAS number and aliases of each component, and the critical point, triple point and vapour-pressure equation that give the saturation pressure of the GAB isotherms, are in the generated table `core/include/adsim/objects/component_properties.h`. `scripts/gen_component_properties.py` writes it from the CoolProp fluid files (`--check` fails when the committed table differs; CI runs it). The values are the ones CoolProp gave, so no result changes.
- `fluid.h` and the precompiled header no longer include `CoolProp.h` (and CoolProp's fmt). CoolProp is still built and linked in this release.
- The component manifest's `properties_source` names the table ("CoolProp 6.7.0 fluid files (generated table)") instead of "CoolProp". Its `heat_capacity_held` field is always false; it stays for one release.

### Removed — Hexane
- Hexane is no longer a component. It was the only component without an ideal-gas heat capacity polynomial (NASA TM-4513), so its heat capacity was CoolProp's value held at the initial conditions. No example deck, adsorbent library entry or cycle template used it. A deck that has Hexane does not load: the error names the component and says that it was removed.
- With it go the held heat-capacity path (`GetGasHeatCapacity`, `Fluid::ComponentsWithHeldHeatCapacity`), the `gas_heat_capacity` run-start warning, and the unused `GetHeatCapacity`, `GetViscosity` and `GetThermalConductivity` CoolProp getters of `fluid.h`. `skarstrom.components_with_held_heat_capacity()` now returns an empty list.

### Fixed — light-product decks amplified a trace component through the vent
- When a step draws gas from a bed end through one stream (for example a light-product draw to
  another step) and another stream at the same end flows back in (a vent valve), the streams that
  leave now carry a mix: the bed's net outflow at the bed-end composition, plus the back-fed gas at
  the vent's own composition (its `y`, normalised). Before, the back-fed gas took the bed-end
  composition, so the draw carried about twice the gas that left the bed, all at the top composition.
  In `Nguyen2023_Exp_VI`–`X`, a component present only at round-off (for example water added to the
  deck) grew about 1.9 times per cycle and changed the cyclic steady state (#270). A vent whose `y`
  is all zero keeps the old behaviour.
- The vents of the light-product decks (`Nguyen2023_Exp_VI`–`X`, `Perez2019_D1`/`E`/`F`,
  `FiveStepProductPurge`, `Perez2019_CaseStudy_Base_LPP`) now carry N2, the feed's inert, as their
  backflow composition. Their results change slightly (see the pull request).

## [0.9.1] — 2026-11-10

### Fixed — the Ward 2022 example builds its decks to the last bit
- `examples/ward_2022.py` computed the paper's dispersion coefficient with `math.pi` and rounded it to 15 digits. The shipped decks use the engine's 15-digit pi at full precision. So the 100 and 200 ccm cases differed from their decks in the last bit. This changed IDA's step sequence, and `test_example_runs_like_its_deck[Ward2022_Q100ccm]` failed on macOS (the v0.9.0 release run). The builder now uses the same formula as the deck, so each built case runs bit for bit like its deck. A new test checks every registry value of the three cases exactly.

## [0.9.0] — 2026-11-10

### Fixed — the Daily release workflow failed with exit code 120 while writing the job summary
- CI: `release.py notes | head -n 200` made Python fail with `BrokenPipeError` once `head` closed the pipe, and `pipefail` stopped the step before the release commit was pushed. The workflow now writes the notes and the release commit stat to files and runs `head` on the files, and `scripts/release.py` restores the default `SIGPIPE` handling so a closed pipe ends it quietly.

### Added — vacuum pump ultimate pressure and step-end pressures
- A vacuum pump stream (a *Pump* device on an outlet) takes an ultimate (blank-off) pressure,
  `endpoint[N].pump.ultimate_pressure` [Pa] (`sk.Pump(..., ultimate_pressure=...)` on a `sk.PumpStream`).
  The flow from the power is multiplied by the rotary-pump speed factor `1 - p_ult/P1`, rounded to a C1
  zero at `p_ult`, so the bed pressure approaches the ultimate pressure from above and is held there.
  The default, 0, is off: no shipped deck changes its results. `PumpBoundaryStream` class version 2.
- The pressure each step ends at: `Simulation::StepEndPressures(cycle)` / `sim.step_end_pressures(cycle)`,
  one record per step with the bed-end pressures (`bottom`, `top`, the last value of the step's bed-end
  series) and `bed_min`, `bed_mean`, `bed_max` over the packed cells. Recorded with the KPIs when each
  cycle completes, so every cycle of a run is available; the CLI's `--json` output lists the last cycle's
  as `step_end_pressures`.
- Sensitivity studies and the optimiser read them through a step reference in place of a component:
  outputs `step[<ID or name>].end_pressure_bottom` (`_top`, `_min`, `_mean`, `_max`;
  `sa.step_output_names()`), and objective and constraint terms such as
  `opt.at_least("step[Evacuation]", "end_pressure_bottom", 1e4)` (C++ `KPIType::END_PRESSURE_*`, appended
  to the enum).

### Fixed — IAST with steep Sips terms and vanishing mole fractions
- The IAST mixing rule failed at the start of a breakthrough into a bed of inert gas. Two causes, both fixed in the IAST solve:
  - The integrator's decaying front toes reach subnormal mole fractions (about 5 × 10⁻³²¹). There the spreading pressure underflowed and the solve returned NaN loadings. Mole fractions at or below 10⁻¹⁵⁰ now take IAST's first-order dilute limit.
  - A Sips-family term with exponent n < 1 has an infinite Henry slope, so the dilute loading of an absent component diverged (−1.3 × 10⁷ mol/kg at y = −10⁻¹⁴ in the Streb 2021 CO₂/CH₄ decks). Inside the IAST solve, such a term is now continued below 10⁻⁴ of the total pressure by a monotone C¹ form whose spreading pressure is exact there. So IAST is unchanged for every component with a mole fraction of 10⁻⁴ or more.
- No bundled deck uses IAST, so no example deck changes.

### Added — validation: Streb (2021) under IAST
- Every Streb 2021 breakthrough case is also run under IAST over the pure Sips terms of the same parameters, and recorded without assertion (`Experiment<case>_IAST` on the Validation page). This is the competition model the paper prefers.
- IAST moves every CO₂ front 0.4–14 % later and every CH₄ front 1–3 % earlier. It is better on 14 of the 47 asserted metrics and worse on 33. Extended Sips therefore stays the decks' model, and no bound changes. The suite's README has the table and the attribution.

### Fixed — NSGA-II example uses the double-exponential evacuation pressure
- `examples/optimise_nsga2.py` and the `skarstrom.optimisation` docstring varied `pressure_exponential_end`, which the evacuation step of `Perez2019_A` no longer has. They now vary `pressure_double_exponential_end`. The Python wheel CI job installs `pymoo`, so this example test runs there too.

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

### Changed — CI: Windows all-tiers job has a longer time limit
- The C++ job on Windows now has a 240-minute limit (Linux and macOS keep 150). The all-tiers run on Windows took more than 100 minutes of tests after a 47-minute cold build, and the old 150-minute limit cancelled it before it finished. Pull-request fast-tier runs are not affected.

### Added — integrator statistics for chunked runs
- `CoupledSolverChuncked::GetStatistics()` returns the integrator counters of a run in slices (steps,
  residual evaluations, error-test and convergence failures, Jacobian evaluations and `recovered_steps`).
  They are the same counters, from the same code, that a blocking run reports in `RunResult::stats`.

### Fixed — Halabi (2011) paper test fails on the start-up REACTION_STIFF warning
- `Halabi2011Test.SorptionEnhanced_RecordOnly` failed its `healthy` check since the `REACTION_STIFF` run-health code was added: five of its twelve runs start from a gas whose CO₂ is replaced by Ar, out of reaction equilibrium,, and the integrator's step collapses in the first 3 ms. The Halabi tests now print every run-health warning (code, cycle, step, value, message, hint) and accept `REACTION_STIFF` only when all steps shorter than 1e-4 s end by t = 0.1 s; any other code, or short steps later in the run, still fails. The results and the tolerances are unchanged (#255).

### Changed — Perez 2019 tutorial shows the current numbers
- The Perez 2019 tutorial now states the results of the current engine for `Perez2019_D1`: CSS at cycle 92, CO₂ purity 94.08 % and recovery 89.48 %, and the Aitken comparison (CSS at cycle 36). Docs only: no result moves.

### Changed
- A sensitivity study that names no outputs now also watches `Component.thermal_energy` (heat supplied per kg of product, listed last) when the deck supplies heat: a step has an enabled energy source with positive power or an enabled heating medium. PSA and VSA decks keep the four defaults (purity, recovery, productivity_mass, energy_consumption). `adsim::sensitivity::DeckSuppliesHeat` exposes the rule.

### Fixed
- The slow test `NewtonKrylovCSSReference.ConvergesBelowTheNoiseFloor` no longer fails on Windows: it accepts a solve that stops at the integrator noise floor (about 1e-8) without reaching the 1e-8 Newton tolerance, and it measures contraction factors only clear of that floor.

### Fixed — automatic releases skipped because main's CI runs were cancelled
- Every push to `main` now gets its own CI run (they no longer share one concurrency group, where GitHub cancelled all but one pending run). The *Daily release* workflow decides "is main green?" from the newest completed CI run for main's newest non-release commit, or else for an ancestor of it from the last 36 hours, and ignores cancelled runs. See `RELEASING.md`.

### Added — re-optimising the CO₂ VSA case study under plant constraints
- `examples/perez_2019_case_study_optimise.py`: NSGA-II (pymoo) over the device settings and schedule of `Perez2019_CaseStudy_Base` and `_Base_LPP` (feed time, blower flow, vacuum-pump power, blowdown share, valve coefficients, the LPP feed share), each candidate run to cyclic steady state in a process pool and warm-started from the nearest cyclic steady state its worker has (`skarstrom.optimisation.run_warm`). Two objective sets: specific electrical energy (kWh per tonne CO₂) against productivity, subject to the purity and recovery targets and the plant limits; and purity against recovery, with or without the plant limits, for the comparison with the paper. Constraints, as NSGA-II constraint violations normalised by their bounds: purity and recovery targets, approach to fluidisation ≤ 1 and bed pressure drop per metre (operating-limit report and the pressure profile, every step), the lowest blowdown and evacuation pressure ≥ `variable[p_evac_floor]`, at least 1 % CO₂ recovery (a run that makes nothing is not a design), cyclic steady state reached with a healthy run, and optionally no backflow through an unconnected outlet. The pump power and the valve coefficients are searched in log10. Every evaluation is appended to `results.csv` (variables, objectives, each constraint's value and margin, CSS cycle, cycles run, warm start, wall time); running the same command again resumes the search from the file.
- "Re-optimising under constraints" in the case-study documentation: how to run it, what each constraint does, and the constrained, unconstrained and paper fronts.

### Added — Perez 2019 case-study base decks
- Example decks `Perez2019_CaseStudy_Base` (four-step cycle of configuration A) and `Perez2019_CaseStudy_Base_LPP` (configuration D1 with light-product pressurisation), built by `examples/perez_2019_case_study.py`. They run the paper's cycles the way a plant runs them: an n-bed schedule in fraction step-time mode (continuous feed and vacuum pump; the step fractions are expressions of `variable[n_beds]`), a feed blower, a vent valve, pressurisation through a valve (from the feed header, or from the partner bed for LPP) and one vacuum pump with a pressure-dependent power and efficiency curve for blowdown and evacuation, with the Ergun momentum balance. Every device setting and the case-study constraint values (evacuation floor, purity and recovery targets) are user variables. The devices reach the validated decks' step-end pressures to within 1.3 %, and purity, recovery and specific energy agree to within 1.3 points and 3 %. See docs/case-study-vsa.md.

### Added — connection pressure warning
- A new end-of-run warning, `connection_pressure` (`settings.warn_connection_pressure`, on by
  default), reports gas that enters a bed through a connection while that bed's end pressure is
  more than 1 kPa above the partner's logged pressure (more than 0.1 % of the cycle feed in one
  step, last cycle; pump streams are not checked). It is a diagnostic: a connected endpoint still
  replays its partner's previous run and does not enforce the order of the pressures. (A cap on
  connected pressure profiles and a pressure-difference scaling of replayed flows were tried and
  withdrawn before release: they made the light-product decks settle more slowly and noisily.)

### Fixed — saved decks and library files list their entries in a fixed order
- Saving a deck or an adsorbent file wrote the entries of every internal hash map (components, isotherms, mass-transfer models, heat capacities, adsorbents) in hash-table order. The order depended on the standard library, so the same deck gave different files on different platforms and the test suite re-saved the committed example decks with spurious differences. These maps are now written in key order. The file format is unchanged and files in either order load. The containers are unchanged, so no result changes.

### Changed — specific surface area of a cylindrical packed layer
- The specific surface area of a packed layer follows the particle shape: `6 (1 - eb) / d_p` for spheres (unchanged) and `4 (1 - eb) / d_p` for cylinders, in the solver (gas–solid heat transfer) and in `MixedAdsorbentLayer::GetSpecificSurfaceArea()` / `LayerGeometrySummary`, with the hydraulic diameter `4 eb / a_v`. Both use one shape factor. The film correlations keep `d_p` as the Reynolds and Sherwood length. Decks with spheres give the same results bit for bit.

### Changed — the log keeps one cycle by default
- `settings.max_cycles_to_keep` defaults to 1 (was 10), and every example deck keeps 1. The log stores
  every solver step of each kept cycle, and with 10 cycles some example decks (Arora2021_LiLSX_PSA,
  Ferre2026_NaYF) ran out of the 2 GB that the web solver has. The KPIs and the CSS history still cover
  every cycle, and the log still keeps the cycles the CSS check reads. Raise the setting only to plot the
  series or profiles of earlier cycles.
- There is no run-health warning when the log keeps fewer cycles than the run: the `HISTORY_TRUNCATED`
  check and its switch were removed before release. `WarningSettings` is class version 5; a deck saved
  with that switch (version 2, 3 or 4) still loads, and the value is ignored.

### Fixed — fast reactions started from a product-free bed
- A reaction that relaxes far faster than the flow, started from a bed that holds none of its products (a steam-reforming catalyst fed methane into a bed of steam and argon), no longer holds the integrator at the reaction's microsecond time scale and fails with "too much work". The cause was the difference-quotient Jacobian: the mole fractions of the reaction's components sat far below their perturbation (the error weight, about 1e-8), so a falling one was perturbed below zero, where the rate law clips the activity, and a power below 1 (H2 at -2.5 in steam-reforming kinetics) was differenced across a secant far flatter than its tangent; the Newton iteration lost the reaction's coupling. The mole-fraction column of every component a reaction in the cell involves is now perturbed away from zero, and by 1e-4 of its own value (or of the activity floor) when the rate raises it to a power below 1. On a Xu-Froment reformer from a steam/argon bed, 1200 s of feed take 0.7 s instead of 52 s (509 steps instead of 11,600) and give the same outlet conversion as a bed started at equilibrium. Decks without reactions are unchanged bit for bit; decks with reactions move by about the integrator tolerance (the reaction tests pass unchanged).

### Added — `REACTION_STIFF` run-health warning
- A new run-health code, `REACTION_STIFF`: in one cycle step, 100 accepted steps shorter than 1e-4 s with the largest local error at a reacting cell (or a failed step after such steps). The warning names the reaction and the component whose near-zero mole fraction drives it, and points to docs/reactions.md, "Start-up of fast reactions". It is counted within the step and also reported for a run that then fails. Switch: `settings.warn_reaction_stiff` (on by default; `WarningSettings` cereal version 4).
- `IReactionRate::MinExponent(component, nu)`: the smallest nonzero power to which a rate law raises a component's activity, and `ActivityFloorMoleFraction(basis, T, P)`.

### Changed — deck view stream types are the stream kind tags
- The deck view (`DeckView`, `sim.deck_view()`) now reports a set-flow stream's `type` as `"flow"`, not `"stream"`. Every stream `type` is now its tag in the `stream` kind list (`flow`, `valve`, `pump`), the tag that `set_stream_kind` takes, so a front end can show the device with the kind list and need not translate. The schema enumerations of a stream `type` and an end `kind` now come from the same kind tables as the edit ops. `kDeckViewVersion` is now 2. `sim.describe()` is unchanged.

### Added — validation: Halabi (2011), sorption-enhanced steam methane reforming
- New paper suite `Halabi2011Test` (`core/test/papers/halabi_2011/`, ctest label
  `paper`): Rh/Ce0.6Zr0.4O2 catalyst mixed with K2CO3-promoted hydrotalcite
  (M. H. Halabi, PhD thesis, TU Eindhoven 2011, doi:10.6100/IR709035). The
  decks are built in C++ from kind names. The thesis's rate model goes into
  three `lhhw` reactions; its two site balances are expanded exactly into one
  denominator of 14 terms.
- Asserted: with an inert sorbent, the outlet CH₄ conversion matches the
  equilibria the thesis prints (Fig. 7 within 0.46 points, bound 1; Table 6
  within 1.95 points, bound 2.5).
- Recorded only: the twelve sorption-enhanced runs of Table 6. CO₂ breaks
  through after 1.6–3.4 min in the model against 4–15 min in the thesis, and
  the conversion has no plateau before breakthrough. The open question is the
  sorbent capacity below 0.05 bar (#245). See the suite's README and the
  Validation page.

### Added — reactions in the Python builder
- `Adsorbent.add_reaction(name, nu, kind=, basis=, effectiveness=, key_component=, adsorption_terms=, **values)` creates a reaction on an adsorbent, and `Adsorbent.reactions` / `remove_reaction` read and remove them. The returned `Reaction` has `set_kind`, `set_stoichiometry`, `set_basis`, `set_effectiveness`, `add_adsorption_term` and `remove_adsorption_term`, and its values (`rx.rate.pre_exp_factor`, `rx.rate.order[C]`, `heat_of_reaction`...) are fields. The `Simulation.build` operations gain `set_reaction_kind`, `set_reaction_stoichiometry`, `set_effectiveness`, `add_adsorption_term`, `remove_adsorption_term` and `remove_reaction`; `skarstrom.build_kinds()` lists `reaction_rate`, `reaction_basis` and `effectiveness`.
- `sim.describe()` lists the reactions of each adsorbent, and `skarstrom.to_script` reproduces them, so a catalytic or sorption-enhanced deck can be built, saved as a script and rebuilt without writing deck JSON. New docs page: Reactions.

### Added — check valve on valve streams, and a warning for backflow through an outlet
- Valve streams take an opt-in check valve (`endpoint[ID].valve.check_valve`, Python
  `sk.ValveStream(..., check_valve=True)`, deck `"Check_Valve"`). An outlet then passes gas only out
  of the bed and an inlet only into it. The closure is smoothed over 10 Pa, so the boundary residual
  stays differentiable as the flow decays through zero. It applies to unconnected streams only;
  connected valves (pressure equalisation, light-product pressurisation) still carry flow both ways.
  Off by default: existing decks keep their results.
- A new end-of-run warning, `outlet_backflow` (`settings.warn_outlet_backflow`, on by default), is
  given when gas flows back into the bed through an outlet with no connection to another step. It
  fires when more than 0.1 % of the cycle feed flows back in through one outlet in one step of the
  last logged cycle. That gas carries the bed-end composition, not the stream's. The warning names
  the step, the endpoint, the moles and the fraction of the feed, and suggests the check valve.

### Fixed — the docs said an outlet's `y` and temperature enter on backflow
- `docs/python.md` said that gas flowing back in through an outlet takes the outlet's `y` and
  temperature. It takes the bed-end composition and temperature (an outlet is zero-gradient), as
  `docs/boundary-conditions.md` says.

### Added — uncounted-inflow diagnostic on the KPI records
- Each per-component KPI record now has `uncounted_inflow` (moles of the component that entered the bed in the cycle through unconnected endpoints that are neither feed nor recycle endpoints of it, from the gross inflow of the KPI quadratures) and `uncounted_inflow_fraction` (that amount over the counted feed). It shows gas that purity and recovery do not count, for example a vent valve that backflows or a stream left out of the feed set. Connected endpoints (pressure equalisation, every connection) are excluded. The fields are in `sim.kpis()`, `sim.kpi_table()`, the `--json` output of `skarstrom run`, and, as `uncounted inflow fraction [mol/mol]`, in a new last column of `kpis.csv`. No existing KPI value changes.

### Fixed — paper fixtures no longer rewrite example decks in the legacy form
- The paper test fixtures build the generic mass-transfer coefficient model, so the full test tier no longer rewrites the example decks in the legacy form.

### Added — cylindrical particles (extrudates, fibres) in packed layers
- A packed layer has a particle shape: `sphere` (the default, and the shape of every existing deck) or
  `cylinder`, an extrudate or fibre of diameter `particle_diameter` = 2 r_f, taken as infinitely long.
  The linear-driving-force diffusion geometry follows the shape: an LDF constant of 8 on the radius
  (k = 8 D / r_f^2) and V/S = r_f / 2 for the film term, from Glueckauf's parabolic profile. The film,
  macropore, Knudsen and capacity-factor forms of the generic coefficient model use it; the micropore
  term stays that of a spherical crystal. Pressure drop, dispersion, the film correlations and the
  specific surface area keep d_p as an equivalent-sphere diameter.
- Registry leaf `layer[L].particle_shape`, deck edit `SetParticleShape` (Python `Layer.particle_shape`,
  `sim.edit` op `set_particle_shape`, `build_kinds()["particle_shape"]`), deck-view layer field
  `particle_shape`, and the `cylinder` LaTeX fragments of the resistance descriptors (the display
  equation takes the layer's shape). `MixedAdsorbentLayer` is at class version 7; older decks load
  as spheres and give bitwise the same results.

### Added — resolved adsorbent library metadata
- `data/adsorbents/library.resolved.json`: the adsorbent library with the engine's own resolution applied (isotherm kind tags and labels, parameters in engine order with units, source points in SI with a summary, completeness, property units and provenance, mass-transfer parameters, DOI URLs). `scripts/docs/gen_adsorbent_library.py` generates it and checks it against `library.resolved.schema.json`; a front end can embed it as it is. `library.merged.json` is unchanged.

### Changed — README: get started with the web interface or an AI agent
- The README's *Get started* section points people to the web interface and AI agents to the repository, with a one-line example prompt.

### Added — registry leaves, descriptors and a display equation for the mass-transfer resistances
- Each resistance of the coefficient model is a Boolean registry leaf,
  `adsorbent[A].kinetics[C].resistances.<tag>` (`film`, `macropore`, `knudsen` while macropore is on,
  `micropore`, `lumped`), and the capacity-factor mode is the leaf `kinetics[C].capacity_factor_mode`
  (`none`, `slope`, `secant`, `fixed`; while film or macropore is on). Setting them by path goes through
  `SetMassTransferResistance` and `SetCapacityFactorMode`, with their rules. A dual-rate slow site has
  `slow_resistances.<tag>` and `slow_capacity_factor_mode`.
- The parameters keep their flat names (`kinetics[C].k`, `.tortuosity`, `.slow_k`, ...) and are listed only
  while their resistance (or mode) is on. The former arrhenius leaf `kinetics[C].pre_exp_factor`
  (`slow_pre_exp_factor`) still resolves: it is the lumped `k` (`slow_k`).
- The `series_resistance` kind descriptor lists the resistances (tag, label, parent, whether the capacity
  factor scales it, parameters, LaTeX fragments for spheres and slabs, sources), the capacity-factor modes
  and the presets; the six former kinds are marked `legacy`. `ModelKindManifest()` emits the model kinds
  as JSON (also `skarstrom.model_kinds()`); `build::RateCoefficientPresets()` and
  `build_kinds()["rate_coefficient_preset"]` list the preset tags.
- `MassTransferEquationLatex(model, shape)` composes the coefficient's equation from the descriptor's
  fragments. The deck view's kinetics entries have `resistances`, `capacity_factor_mode`, `preset` (the
  derived preset tag, or `custom`), `equation` and `equation_slab`, and the same with `slow_` for a
  dual-rate slow site.
- Python: `set_rate_model(..., resistances=[...], capacity_factor_mode=...)`,
  `Adsorbent.set_rate_resistance(component, tag, on, site="fast")`, and the edits
  `sim.edit.set_mass_transfer_resistance` and `sim.edit.set_capacity_factor_mode`. `to_script` writes the
  switches of a custom model.

### Changed — the coefficient kind is always series_resistance
- A coefficient model reports the kind `series_resistance` (the deck view's `rate_coefficient`,
  `describe()`, `build::RateCoefficientKinds()`); the preset it matches is a separate, derived value.
  The former kind tags are still accepted by the builder, `set_rate_model` and `SetMassTransferKind`.
- A model created from the tag `series_resistance` has the class defaults (the lumped term only,
  k = 0.5 1/s); its terms no longer switch on when a parameter becomes non-zero. Set the switches instead
  (`resistances=` or the registry leaves). `SetMassTransferKind(..., series_resistance)` on an existing
  model keeps its switches, mode and values. Version-1 decks load with the switches their zero
  parameters meant, and their results do not change.

### Changed — example decks and library entries stored on the generic mass-transfer coefficient model
- The bundled example decks and the adsorbent library entries (`data/adsorbents/*.json`, `library.merged.json`) now store their mass-transfer coefficients in the `series_resistance` form (explicit resistance switches, capacity-factor mode, every parameter present) instead of the old fixed-combination classes. The file format changes only; results are unchanged, bit for bit. A new test checks this for every example deck.

### Added — water and humid CO2 parameters for three more library entries
- Lewatit VP OC 1065 and Purolite A110 (Low et al. 2023) now include water: the Young-type GAB isotherm (`gab_young`, Table 9, Eqs. 17-22), the limiting water heat of adsorption (Table 10) and the measured water isotherm points (288.25, 298.15 and 308.15 K).
- Triamine-grafted alumina (Grossmann et al. 2025) now includes the GAB water isotherm (Table 5), the water heat of adsorption and the CO2 weighted-average dual-site Toth (`toth_wadst`, Tables 4 and 6, referenced to T0 = 308 K).

### Changed — Lewatit VP OC 1065 CO2 isotherm is now the humid WADST fit
- The CO2 isotherm of `lewatit_vp_oc_1065_low2023` is the `toth_wadst` fit of Young et al. (2021) instead of the dry td-Toth fit of Low et al. (2023). At zero water loading it is the dry site of that fit, which differs from the earlier dry isotherm. The dry td-Toth parameters stay in the entry notes.

### Fixed — removing a component left an empty isotherm object
- Removing a component from an adsorbent (including loading a library adsorbent into a deck that lacks one of its components) now also drops an isotherm object that no component uses any more. Before, a library adsorbent with a Water isotherm loaded into a deck without Water left an empty isotherm in the deck, which the deck-to-script round trip could not reproduce. Results do not change: the empty isotherm was never evaluated.

### Changed — one mass-transfer coefficient model with a switch per resistance
- The `series_resistance` coefficient (`MassTransferCoefficientSeriesResistance`, cereal version 2) is now
  the one coefficient model: external film, macropore (with a Knudsen sub-switch), micropore and a lumped
  fitted coefficient (k with an optional activation energy), each switched on or off explicitly and added in
  series, and an explicit capacity-factor mode (`none`, `slope`, `secant`, `fixed`). Every parameter is kept
  whatever its switch, so switching a term off and on restores it.
- The older kinds (`constant`, `arrhenius`, `macropore`, `macropore_knudsen`, `haynes`, `macropore_capacity`)
  are presets of it. Decks and library entries that store them still load (the old classes are kept as
  load-only types) and are converted on load; each preset is evaluated in the old class's operation order,
  so results are bit for bit those of the old kind. Every old tag is still accepted by the builder,
  `set_rate_model` and `SetMassTransferKind`, and the reported kind is the preset the switches match.
  Version-1 `series_resistance` archives load with the switches their zero sentinels meant.
- New deck edits `SetMassTransferResistance` (one switch, fast or slow site; the last term cannot be
  switched off, Knudsen needs macropore) and `SetCapacityFactorMode`. Each dual-rate site has its own
  switches.
- New deck checks: no resistance on, or an enabled term with a parameter out of range, is an error; a film
  term with a zero film resistance factor, or a second-order LDF with a physical resistance, is a run-start
  warning.

### Changed — Lennard-Jones parameters cited to NASA TR R-132
- The Lennard-Jones parameters are now cited to Svehla (1962), NASA TR R-132, Table I(a). The values of all other components do not change. Octane (452.09 K, 6.362 Å) and toluene (472 K, 5.240 Å) take CoolProp's viscosity-fit pairs. Isopentane uses the n-pentane row (estimate). Octane D_m rises about 11 %, toluene about 2 %.
- C++: `ComponentData::source` (`LJSource`) and `Components::LennardJonesCitation(name)`.

### Added — Newton–Krylov cyclic steady state
- New CSS acceleration method `newton_krylov` (`settings.acceleration_method`, `AccelerationMethod::NEWTON_KRYLOV`): Newton's method on G(x) − x = 0 over the bed's differential state, with GMRES and a finite-difference product of the cycle map's Jacobian (one plain cycle per Krylov vector, run in a scratch copy so the logged run is untouched). Damped by the Glover rule, inexact (Eisenstat–Walker forcing), with a fallback to Anderson cycles when a Newton step does not reduce the residual. It reaches a CSS well below the stopping rule's tolerance, and from a warm start it needs one or two Newton steps. Off by default; the other methods' results are unchanged.
- `sim.css_history()` / `CSSRecord` report each Newton–Krylov solve: `newton_iterations`, `newton_cycles` (cycle equivalents), `newton_residuals`, `newton_krylov_vectors`, `floquet_multiplier` (the dominant Floquet multiplier estimated from the Krylov basis) and `newton_outcome`; the solver log at level *medium* prints the same.

### Changed — automatic releases
- Releases are automatic: the new *Daily release* workflow (`.github/workflows/daily-release.yml`) releases `main` once a day when it has changelog fragments and its CI is green. It cuts the changelog, commits "Release X.Y.Z" to `main` and starts the *Release* workflow, which tags, creates the GitHub Release and publishes to PyPI. The *Prepare release* workflow and its "Release X.Y.Z" pull request on `bot/release` are retired; a hand-picked bump or pre-release is now a manual run of *Daily release* (`bump`, `suffix`, `dry_run`). Releases can be paused with the repository variable `AUTO_RELEASE=false`. *Release* accepts a `ref` input (a commit on `main` to release) and never releases a development version. See `RELEASING.md`.

### Added — honest development versions
- Between releases `main` carries a PEP 440 development version, the next patch with `.dev0` (0.8.1.dev0 now, 0.9.1.dev0 after 0.9.0), so a source build never reports a released version. A development build from a git checkout also names its commit: `skarstrom --version`, `skarstrom.__version__` and `adsim::kVersionString` read e.g. `0.9.1.dev0+g0a19084` (new `ADSIM_VERSION_LOCAL` in `adsim/version.h`; read when CMake configures). Release builds are unchanged (`0.9.0`). New `scripts/release.py develop` and `check-dev`; CI fails if `main` has an already-tagged version without a dev marker.

### Changed — reported KPIs on the exact integrals
- The reported KPIs (purity, recovery, productivity, specific electrical energy, production rate, `mass_balance_error`, and the moles fed behind the reaction KPIs) now use the KPI quadratures the solver integrates alongside the solution, not the trapezoid over the logged samples. Reported KPIs move in about the fourth significant figure or below (1e-7 to 2e-5 relative on Perez 2019 A–C). The trapezoid values are the new `trapezoid_*` fields of `KPIRecord` (Python `sim.kpis()`, `kpi_table()`, CLI `--json`, the KPI manifest). They replace the `quadrature_*` fields and `quadrature_mass_balance_error` (renamed `trapezoid_mass_balance_error`) added in this release's KPI quadratures entry, which never shipped. An endpoint without quadratures (turned off, or a failed step) falls back to the trapezoid. `KPIParams::integration_` now defaults to `KPIIntegration::QUADRATURE`.
- At an outlet, an endpoint's component flow, logged and integrated, is its molar flow times the bed-end face mole fraction the species balance convects, not the bed-end average (#171). Gas leaving through an inlet (a back-flow) is counted the same way. The cycle component balance now closes to 1e-7 to 3e-7 on Perez 2019 A–C, the same as the face fluxes (9e-6 to 2e-5 before). The logged stream composition keeps the bed-end sampling, so connected streams are unchanged. A buffer tank's average composition is built from the component flows, so it moves by the same small amount.
- `productivity_volume` and the bed-volume space-time yield divide by the packed cells' volume, without the two ghost cells (and without void end cells) (#201). `productivity_volume` rises by (N+2)/N for N cells.

### Added — humid DAC breakthrough example
- Example deck `Illustrative_Humid_DAC_Breakthrough`: a non-isothermal breakthrough of 400 ppm CO₂ in humid
  N₂ (y_H₂O = 0.0115, 49 % RH at 20 °C, 1 bar) on the APDES-NFC library adsorbent, with the GAB water
  isotherm, the humid Toth CO₂ isotherm that reads the adsorbed water loading, LDF kinetics and the water
  heat of sorption. An integration test checks that the run is healthy and the mole balances of all three
  components close, that water breaks through first and the outlet loadings reach the GAB and humid Toth
  values at the local state, that a dry feed reaches the lower dry CO₂ capacity, and that the temperature
  peaks before the CO₂ front.
- `docs/isotherms.md` evaluates the APDES-NFC dry and humid loadings at the feed state; the TSA tutorial
  (`docs/tsa-regeneration.md`) has a "Humid feeds" section.

### Added — reaction KPIs and reaction-rate series
- Reaction KPIs, configured by a key reactant, the products to report and a space-time-yield basis (`kpi.reaction.key_reactant`, `kpi.reaction.product[C]`, `kpi.reaction.basis`; `build::SetReactionKPI`; Python `sim.kpi.track_reaction(...)`). Per cycle, from the moles the reactions made in the bed and the key reactant's moles fed: `conversion` X_A = -G_A / F_A,in on the key reactant's record; `selectivity` G_P / (-G_A), `reaction_yield` G_P / F_A,in and `space_time_yield` (kg per kg of catalyst, or per m3 of bed, per second) on each product's record; the cycle's `peak_temperature` with its position and time, and `reaction_generation`, on both. Selectivity and yield are plain molar ratios, without stoichiometric normalisation. The key reactant's moles fed are its KPI feed endpoints, or every endpoint's inflow of it when it has none. The new `KPIRecord` fields are in `Simulation::KPIs()`, `sim.kpis()`, `kpi_table()`, the CLI's `--json` and `kpis.csv`, and the KPI manifest (flagged `reaction`); records of a deck without reactions are unchanged, and so is its `kpis.csv`. The deck view gains an optional `reaction_kpi` section (view version unchanged). A component removed from the deck leaves the reaction KPIs.
- The rate of each reaction is reachable like the other logged series: `Simulation::Series` / `Profile` take a `reaction` (name or ID, with `adsorbent` to pick the solid), Python `sim.series("reaction_rate", reaction="SMR")` / `sim.profile(..., reaction=...)`, `sim.series_reactions()` and `skarstrom.reaction_variables()`.

### Added — water and humid parameters in the adsorbent library
- The APDES-NFC library entry (`apdes_nfc_stampi2020`) now includes water: the GAB isotherm of Stampi-Bombelli et al. (2020), the water rate constant and heat of adsorption, and the humid modified Toth for CO2 (`toth_humid`) with the corrected γ and β of the 2022 correction. At zero water loading it reduces to the dry fit.

### Added — mixed adsorbent layers documented, with an example deck
- New documentation page [Mixed Adsorbent Layers](https://samaffleck.github.io/skarstrom/mixed-layers/): how a layer's composition is stored and checked, the mixing rules with their equations, the shared-temperature assumption, mixture versus thin layers, per-solid results, and how to set a composition from Python and the CLI (`layer[ID].mass_fraction[ID]`).
- New bundled deck `Illustrative_Diluted_Bed`: the Ward 2022 100 ccm CO2 breakthrough on activated carbon diluted with 30 % glass beads by mass in a mixed layer (illustrative, not a validated reproduction).

### Added — reactions enter the balances
- The heterogeneous reactions on a solid (`adsorbent[ID].reaction[ID].*`) now act in a run. Per cell, each reaction's rate (at the gas composition and pressure and the catalyst temperature: the gas's in LTE and isothermal runs, the solid's in LTNE) times its effectiveness factor and the solid's bulk density gives the component sources S_i = sum nu_i rho_b eta r, which enter the component balances; their sum enters the total balance, so a reaction that changes the number of moles moves the velocity field. The heat of reaction (Kirchhoff by default) enters the bed energy balance in LTE and the solid's in LTNE; an isothermal run has no heat term. Reactions add no unknowns, and a deck without reactions gives the same results as before.
- New logged series, only for a deck with reactions: `reaction_heat` (bulk, W/m3), `reaction_source` (per component, mol/(m3 s)) and the rate of each reaction on each solid (`"<reaction> Rate (<solid>)"`, mol/(kg s), `EquationLog::GetReactionRate`), listed in `SeriesManifest()` (`reaction_variables` for the rates). The adsorption source series `component_mass_source` and `energy_source` are unchanged.
- The Thiele effectiveness model is evaluated: the generalised (Aris-Bischoff) modulus of the key reactant for a spherical, isothermal particle (radius from the layer's particle diameter, the solid's density, D_eff = eps_p / tau times the cell's molecular diffusivity), integrated along the stoichiometry with an 8-point Gauss-Legendre rule down to the equilibrium concentration of a reversible rate. For a first-order rate it is eta = 3 / phi^2 (phi coth phi - 1).
- The mass-balance CSS methods, the cycle mole-balance diagnostic and its run-health check count the moles the reactions make, |in + G - out| / in, so a reactive deck can reach cyclic steady state. The reaction warnings are written to the run-start log.

### Changed — inert loadings are no longer solved
- The loading of a component on a solid where its isotherm is Inert is no longer a DAE unknown: it is held at exactly 0 (logged, in the KPIs and the per-adsorbent series as before). A pure catalyst or support solid no longer adds a loading unknown per component per cell. Decks with inert components solve a smaller system, so their results move by about the integrator tolerance (the re-saved example decks' embedded numbers change in the last digits), and a cycle state saved from such a deck before this change no longer fits its layout.

### Added
- Agent guide at the top of the "Using Skarstrom" half of `AGENTS.md`: check `skarstrom --version` and read the matching docs, prefer the wheel, run a few cycles first, treat run-health warnings as wrong answers, read "Matching Another Model" before comparing with another model, and where the machine-readable pointers are.
- The docs site publishes `llms.txt` (an index of the pages by nav section) and `llms-full.txt` (every page's Markdown in nav order), generated from the MkDocs nav by `scripts/docs/llms_txt_hook.py`. URLs follow `site_url`, so the released root and `/dev/` each point at themselves.

### Changed — documentation site shows the latest release
- The GitHub Pages root is built from the latest `v<semver>` tag with that tag's own configuration, so it
  describes the installable wheel. `main` is previewed under `/dev/` with a "not on PyPI yet" announcement bar.
  The site is rebuilt on every push to `main` and after every release.
- Convention: documented new names carry a "New in X.Y" note (AGENTS.md); applied to the run-health page and
  the `recovered_steps` statistic.

### Added — per-adsorbent series for mixed layers
- In a bed with a layer of several solids, the solver log now also holds each solid's own adsorbed loading, equilibrium loading, mass-transfer coefficient and (with dual-rate LDF) slow-site loading, per kg of that solid, named `"<component> Adsorbed Concentration (<adsorbent>)"`. The mass-weighted aggregate series are unchanged. Beds without a mixed layer log nothing new.
- `Simulation::Series` / `Profile` take an `adsorbent` (a name or deck ID); Python `sim.series(..., adsorbent=)` and `sim.profile(..., adsorbent=)` take a name, a deck ID or an adsorbent object. For a bed without a mixed layer the per-adsorbent series are read off the aggregate. `sim.series_adsorbents()` lists the adsorbents, `skarstrom.per_adsorbent_variables()` the variables. `save_logs` and the CLI's `--csv` logs write them as extra CSV files.
- `SeriesManifest()` (`skarstrom.series_manifest()`) lists every bed series variable with its log name pattern, unit and whether it is also logged per adsorbent.

### Added — documentation code blocks run in CI
- `scripts/docs/run_doc_blocks.py` executes every Python and shell code block of `docs/*.md` and
  `AGENTS.md` against the installed package and reports each block's file, line, result and time; CI runs
  it in the Python wheel job (about 35 s, each run capped at 3 cycles). `--full` runs the uncapped set,
  including blocks marked `<!-- long -->`, for a nightly run. `<!-- no-run -->` excludes a block and a
  hidden `<!-- doc-setup -->` block defines what a fragment assumes; the convention is in AGENTS.md.

### Fixed — documentation examples the engine rejected
- A call to `sim.edit.set_pump_kind` on a valve stream, `sim.stream_series` with a wrong argument list and
  an `evaluate_iast` call with a pressure array (AGENTS.md); a `fit_isotherm` example whose data were
  overwritten by the example before it, and one that wrote a fit of one isotherm kind into a deck of
  another (isotherm fitting); `enumerate(sim.steps)` formatted as a name (TSA regeneration); a variable
  and pressure path the example deck does not have (global variables); the evacuation pressure path of the
  NSGA-II example (parameterisation, Perez 2019); and several fragments that used
  names no earlier block defined.

### Changed — recovered integrator failures are quiet and counted
- SUNDIALS error messages go to the engine log instead of stderr: the solver replaces the SUNDIALS
  context's default error handler, so nothing from IDA is printed directly any more. Every IDA message
  appears in the log at level *high* (`IDA: ...`).
- A step that IDA fails but the engine's retry recovers is logged at *medium* as a "recovered step" with
  IDA's message (the chunked solver now logs these too). Only a step whose retries run out is logged at
  *low*, with IDA's last message. Results and retry behaviour are unchanged.

### Added — `recovered_steps` run statistic
- `IntegrationStats::recovered_steps` (C++ `RunResult.stats`), `result.stats["recovered_steps"]` (Python)
  and `stats.recovered_steps` in the CLI `--json` summary and `run.json`: the integrator failures the
  engine recovered by retrying. The run-health `STEP_RETRIES` measure reads the same count.
- Docs: starting a flow inlet into a bed at low pressure (ramp it over about 1 s) and what a recovered step
  means in the log (*Boundary conditions*).

### Changed — every cycle's KPIs and CSS record are kept
- `settings.max_cycles_to_keep` now limits only the logged series and profiles (the memory-heavy data).
  `sim.kpis()`, `sim.kpi_table()`, `sim.css_history()` and the CLI's `kpis.csv` / `css.csv` cover every
  cycle of the run: the solver records each cycle's CSS check (`CSSHistory`, `GetCSSHistory()`) when the
  cycle completes, as it already did for the KPIs. One record is about 200 bytes, so a 2000-cycle run adds
  well under 1 MB. Results are unchanged.

### Added — DAC breakthrough example and kinetics documentation
- Example deck `Illustrative_DAC_Breakthrough` (script `examples/illustrative_dac_breakthrough.py`): an
  isothermal 400 ppm CO₂ breakthrough on the Lewatit VP OC 1065 library adsorbent using the dual-rate LDF
  and the series-resistance coefficient together. An integration test checks its mole balances, the early
  front and slow tail, and that split = 1 with unit capacity factor and no micropore term reproduces the
  LDF + Haynes run bit for bit.
- `docs/kinetic-models.md` has a "Which to use for direct air capture" section and one notation for both
  kinetic sections; `docs/python.md` and the deck-builder kind list name `dual_rate_ldf`,
  `series_resistance` and `macropore_capacity`.

### Added — series-resistance LDF coefficient
- A **series-resistance (Glueckauf) LDF coefficient** (`MassTransferCoefficientSeriesResistance`,
  kind `series_resistance`): 1/k = f K r_p/(3 k_f) + K τ r_p²/(15 ε_p D_m) + r_c²/(15 D_c), external
  film (k_f from the layer's film correlation, as `haynes`), macropore and micropore resistances in
  series, with the local capacity factor K = ρ_p ∂q*/∂c from the isotherm slope (secant where the
  slope is not positive, floored at 1; or a fixed `capacity_factor`) and an Arrhenius micropore
  diffusivity D_c = D_c0 exp(−E_a/RT). Each term is switched off by its parameter. With the
  micropore term off and K = 1 it equals `haynes` bit for bit. Defined once in `model_kinds`, so the
  registry (`kinetics[C].crystal_radius`, `.micropore_diffusivity`, ...), the builder and Python
  `model_kinds()` pick it up. Existing decks are unchanged.

### Added — run-health warnings
- Cheap checks after every cycle flag a run that completes but should not be trusted, as structured
  warnings `{code, value, threshold, cycle, step, suspects, hint}`: `MOLE_FRACTION_SUM` (the gas mole
  fractions drift from summing to 1), `MASS_BALANCE_STUCK` (a component's cycle mole balance stays open
  while the end-of-cycle state has settled: a model error, not slow CSS) and `STEP_RETRIES` (many rejected
  integrator attempts per accepted step). Thresholds leave every example deck and paper case silent
  with margin; the measured distribution is in the new docs page *Run health*.
- `RunResult.healthy` and `RunResult.health_warnings` (C++), `result.healthy` and `result.warnings`
  (Python), `Simulation::RunHealth()` / `sim.run_health()` with the per-cycle measures. The run status,
  exit code and results are unchanged. The solver logs the warnings once at the end of the run (LOW).
- CLI: the `--json` summary has `healthy` and `run_health` (structured; their messages are also appended
  to `warnings`), the text summary a `Health:` line, and run-health warnings go to stderr even with `-q`.
- Switches `settings.warn_mole_fraction_sum`, `warn_mass_balance_stuck` and `warn_step_retries` (on by
  default; `WarningSettings` class version 2), listed in the warning
  manifest with their `code`.

### Added — dual-rate LDF rate model
- New rate model `dual_rate_ldf`: two parallel LDF sites share the component's equilibrium loading. The fast site holds the fraction `split` of q*, the slow site the rest, and each site has its own coefficient model of the same kind (registry leaves `adsorbent[A].kinetics[C].split` and `slow_<name>`, e.g. `slow_k`). The slow-site loading is a new state variable. It has its own initial condition, it is part of the CSS state test and the warm-start remap, and it is logged as `slow_adsorbed_concentration`. With `split = 1` the model is the first-order LDF exactly, with no extra unknown. Use it for the slow uptake tails of direct-air-capture breakthroughs (Stampi-Bombelli et al. 2024, doi:10.1021/acs.iecr.4c01392). Decks that do not use the new model load and run unchanged.

### Changed — the Python package takes component formulas from the engine's component manifest
- The formula → CoolProp name table behind `sim.steps[0].bottom.y["CO2"]` and the other typed mappings is read from the component manifest (once per process) instead of a copy in `model.py`. More formulas and CoolProp aliases are accepted than before.

### Added — `skarstrom.component_aliases()` and reference-client tests
- `skarstrom.component_aliases()` returns the formula and alias → CoolProp name table.
- Tests that act as a reference client of the engine's views: for every bundled deck the deck view validates against its schema and every value sets back by its path; every `sim.edit` cascade and kind (read from `model_kinds()`) leaves a deck that passes `check()` and round-trips; the manifests are non-empty and consistent.
- Docs: "Views and manifests" in the Python guide.

### Changed — one definition per model kind
- The kind descriptors of `adsim/objects/model_kinds.h` are now the single source for each model kind. They provide the builder's kind lists, the registry's entries for each kind's parameters (units, bounds, labels, tunable flag), the isotherm kind tags and flags (`IsothermKindString`, `IsCompetitiveIsothermKind`, `IsWaterCoupledIsothermKind`, `PairedIsothermKind`), the IAST family table and the isotherm fitter's kind list. A test checks that the descriptors agree with the model classes (names, order, defaults), the registry and the fitter's default parameters.
- The registry output of every example deck is unchanged. The one visible difference is a display label: the ramp time of an exponential temperature or flow profile is now named "exponential ramp time", as it already was for pressure profiles.

### Added — KPI diagnostics and Cp provenance in the display manifests
- The KPI manifest lists the `KPIRecord` diagnostics (`quadrature_*` KPIs, `mass_balance_error`, `quadrature_mass_balance_error`) after the reported KPIs, with `diagnostic: true`. The reported set and `sensitivity::KPINames()` are unchanged, and the CLI tables and `kpis.csv` leave the diagnostics out.
- The component manifest gives each component's `cp_species` and `cp_source` and a top-level `cp_reference`. **Changed:** `ComponentManifest()` / `skarstrom.component_manifest()` now return `{cp_reference, components}` instead of a bare list.

### Added — registry paths for the cost estimate, the sensitivity study and KPI tracking
- `kpi[Component].tracked`, `cost.*` (plant, finance and operating values, `cost.vessel.*`, `cost.heat_exchanger.*`, `cost.compressor[N].*`, `cost.vacuum_pump[N].*`, `cost.adsorbent[ID].cost_per_kg`) and `sensitivity.*` (the stored study's seed, Morris and Monte Carlo settings). Python: `sim.kpi_tracked`, `sim.cost`, `sim.study_settings`.
- A test checks, leaf by leaf, that every editable value of every example deck and of a deck with every reaction kind has a registry path; the leaves that are structure (set through `deck_edit.h` / the builder), values of a variant not in use, legacy fields or run results are listed in the test with their reason.

### Added — complete deck-editing API
- `adsim/objects/deck_edit.h` covers every structural edit a front end makes, each with its cascade: `AddLayer` / `RemoveLayer` / `MoveLayer`, `AddStep` / `RemoveStep` / `DuplicateStep` / `MoveStep`, `AddStream` / `RemoveStream`, `ConnectEndpoints`, and kind switches that keep every parameter the old and new kinds share (same name and unit): `SetIsothermKind`, `SetMassTransferKind`, `SetDispersionKind`, `SetHeatDispersionKind`, `SetInternalHeatTransferKind`, `SetExternalHeatTransferKind`, `SetParticleHeatTransferKind`, `SetBoundaryConditionKind`, `SetStreamKind`, `SetFlowKind`, `SetProfileKind`, `SetPumpKind`, `SetPumpEfficiencyKind`, `SetPumpPowerKind`. Also `EnsureIds` (stable IDs by first appearance, connection pair ids 1..N), `RepairConnections`, `PruneBufferTanks`, `SyncKpiConfigToComponents` / `PruneKpiConfig` and `FillDeckGaps`.
- `adsim/objects/model_kinds.h`: one descriptor per model kind (tag, label, ordered parameters with unit, default, bounds and tunable flag; isotherms also carry their IAST family and fitter flags) and `ModelFields`, which binds a live model's parameters to the descriptor names.
- Python: `sim.edit` (`skarstrom.Editor`) exposes the cascades by kind tag, and `skarstrom.model_kinds()` lists the descriptors.

### Changed
- The builder's `Connect` gives a new connection a pair id no other connection of the deck uses.

### Changed — one definition of the gas constant and the standard atmosphere
- `kIdealGasConstant` and the new `kStandardAtmosphere` (101325 Pa) in `adsim/math/constants.h` replace the literals in the expression parser, the boundary-condition defaults and the initial-condition pressure default. Results are unchanged.

### Added — merged adsorbent library file
- `scripts/docs/gen_adsorbent_library.py` also writes `data/adsorbents/library.merged.json`, the validated catalogue with every entry's provenance and serialised adsorbent in one file, in the order of the docs page; `--check` (run in CI) fails when it is out of date.

### Added — deck view, view schema and display manifests for integrators
- `adsim::DeckView` (`adsim/objects/deck_view.h`; Python `Simulation.deck_view()`): a read-only JSON view of a deck's structure. Each node carries the parameter-registry entries that belong to it (path, value, text, unit, bounds, choices, label). `DeckViewSchema()` (Python `skarstrom.deck_view_schema()`) is its JSON Schema, with the builder's kind names as enumerations.
- Display manifests in `adsim/objects/manifests.h` and the Python package: KPIs (`kpi_manifest`), unit categories with suggested units resolved by the unit parser, plus `check_unit` (`unit_manifest`, `check_unit`), advisory warnings (`warning_manifest`), example decks (`example_manifest`, declared in `examples/example_decks.json`) and components (`component_manifest`). The unit parser exposes `UnitCategoryOf` and `UnitOverrideSymbols`.
- Docs page "Deck View and Manifests".

### Changed — CLI reads KPI columns and example families from the engine
- The `run` KPI table and `kpis.csv` take their columns from the KPI manifest. The `kpis.csv` headers and order are unchanged; the terminal table now shows production rate before the energies. `skarstrom examples --json` takes each deck's family from the example manifest.

### Changed — one ideal-gas heat capacity set per component, to 1000 K
- Each component's ideal-gas heat capacity polynomial (`ideal_gas_heat_capacity.h`) now holds one NASA TM-4513 coefficient set, the low-temperature one, valid to 1000 K. The 1000–6000 K set is dropped: adsorption processes stay well below 1000 K. Above 1000 K the value at 1000 K is used (as below the lower limit). Results do not change below 1000 K.
- C++: `NasaPolynomial` is `{t_min, t_max, a}` (a0..a6); `t_mid`, `low` and `high` are gone.
- The display name of the default gas heat capacity mode is now "Ideal gas" (was "Ideal gas Cp(T)"); the registry name is unchanged (`ideal_gas`).

### Changed — paper decks hold gas in the bed voids only
- The Perez 2019 and Ward & Pini 2024 decks hold the gas in the bed voids only, as the papers do (Perez et al. 2019, SI Eqs. S1, S2, S6; Ward & Pini 2022, Eqs. A1, A2, A4): the adsorbent's `particle_void_fraction` is 0, and the particle voidage stays in the Perez decks' macropore LDF coefficient, so the rate is unchanged. Results move by up to 1.1 pp (Perez F purity 75.65 → 76.78 %); five Perez bounds are loosened to the new measured errors.

### Added — provenance of the ideal-gas heat-capacity table as data
- `NasaPolynomial` gains `species` (the TM-4513 species label) and `source` (its data-source code, or Air's mixing note); new `kIdealGasHeatCapacityReference`. No numeric change.
- Python: `skarstrom.ideal_gas_heat_capacity_reference()`, `ideal_gas_heat_capacity_source(component)` and `components_with_held_heat_capacity()`.

### Fixed — slow-tier CSS tests on the Perez 2019 decks
- `CSSAndersonDecks` and `CSSDistanceDecks.Perez2019EReportsTheRemainingDistance`
  use stops and reference KPIs re-measured on the current `Perez2019_A` and
  `Perez2019_E` decks. The decks' constant gas heat capacity moved the plain
  stops (A: 48 → 47, E: 205 → 206) and the long-run reference KPIs, so the
  tests' earlier constants no longer matched. No tolerance changed.

### Added — KPI quadratures
- The solver integrates each endpoint's molar flow in and out per component, its pump or compressor power, the buffer-tank summary integrals and the species flux through the bed faces as IDAS quadrature states alongside the solution. KPIs from these integrals are reported beside the trapezoid KPIs (which stay the reported values) as `quadrature_purity`, `quadrature_recovery`, `quadrature_productivity_mass`, `quadrature_productivity_volume`, `quadrature_energy_consumption` and `quadrature_production_rate`, with the cycle mole-balance error from either (`mass_balance_error`, `quadrature_mass_balance_error`), in `sim.kpis()`, `kpi_table()` and `--json`. Quadrature error control is off, so runs and reported results are unchanged.
- Per-step mole balance (column inventory change against the integrated endpoint and face flows; for a step with one open end, the inventory-difference cross-check of its flow), in the solver log at the High level and in `BoundaryLog::GetStepBalance`.

### Added — CO2 isotherms that depend on adsorbed water
- Two isotherm kinds for CO2 on sorbents where adsorbed water raises the CO2 capacity. In both, the CO2 loading depends on the adsorbed loading of the component named `Water` on the same adsorbent, in the same cell. Water keeps its own isotherm (for example `gab` or `gab_young`) and its own rate model.
  - `toth_wadst` is the weighted-average dual-site Toth of Young et al. (2021, Energy Environ. Sci. 14, 5377), a dry and a wet temperature-dependent Toth site weighted by exp(-A / q_w).
  - `toth_humid` is the humid Toth of Stampi-Bombelli, Storione and Mazzotti (2020, Adsorption 26, 1183; corrected 2022): Toth 2, with the capacity divided by (1 - gamma q_w) and the affinity multiplied by (1 + beta q_w).
  - The defaults are the papers' Lewatit VP OC 1065 and APDES-NFC values. With no water adsorbed, or no `Water` component, each kind is its dry Toth exactly.
- Isotherms can read the adsorbent's actual loadings (`IsothermUpdateParams::q_map`); only the new kinds use them.
- At the start of a run, a deck with these kinds takes one more equilibrium pass, so CO2 starts in equilibrium with the initial water loading.
- The IAST mixing rule rejects these kinds with a deck-check error. A deck that uses them without a `Water` component gets a warning.
- The isotherm fitter supports both kinds. A point's water loading is `isofit::Point::q_w`, or `water_loading=` in `skarstrom.fit_isotherm`. Dry points fit the dry Toth and humid points fit the water terms; a block the data cannot identify is held.
- `skarstrom.evaluate_isotherm` takes `water_loading=`, and `adsim::EvalIsotherm` and `IsostericHeat` take `q_water`.
- Existing kinds and decks are unchanged.

### Changed — the engine links IDAS instead of IDA
- The `adsim` library now links SUNDIALS IDAS only (it was linked against IDA and IDAS, with IDA's copies of the shared entry points winning). IDAS has every IDA entry point, so runs are unchanged (the end state and the integrator counters of six example decks are bit-identical), and its forward sensitivities can act on the solver's own integrator.

### Added — prototype: cycle-map Jacobian-vector products by forward sensitivities
- Internal, experimental `CycleJacobianVector` (`adsim/solver/cycle_sensitivity.h`): one cycle of a single-bed deck without connections with IDAS forward sensitivities along given directions of the initial differential state, returning the end-of-cycle state and Φ_x·v. The sensitivities restart with the integrator at every step switch. `SolverData` gains a re-initialisation observer for it (null by default). Not exposed in Python. The measured cost and accuracy, and where it fails (across the step switches of a PSA cycle), are in the Cyclic Steady State page, "Cycle-map Jacobian by forward sensitivities".

### Added — Haghpanah et al. (2013) four-step 13X PVSA validation
- New paper suite `Haghpanah2013Test` (`core/test/papers/haghpanah_2013/`) and seven bundled decks with their example script `examples/haghpanah_2013.py`: the four-step CO₂/N₂ PVSA on zeolite 13X of Haghpanah et al. (2013, doi:10.1021/ie302658y), with the paper's loading-dependent LDF coefficient, exponential pressure profiles and compressor/vacuum-pump energy. Cases: the Sect. 6.1 reference design, the five minimum-energy designs of Table 6 (compared with Haghpanah's energies and with Glover et al. 2026, arXiv:2606.21523, Table 4), and Glover's CSS benchmark design (dominant multiplier 0.9513 against 0.9509).
- Purity, recovery and productivity of the Table 6 designs agree with Glover's runs within 1.6 pp and 3.4 %; the specific energy is 12 % above Table 6 at evacuation pressures of 0.02–0.03 bar and the Sect. 6.1 purity 9 pp below the paper (open; attribution in the paper README).

### Changed — paper decks use their paper's gas heat capacity
- The Casas 2012, Ferré 2026, Nguyen 2023, Estupiñan Perez 2019, Vysyaraju 2024, Ward & Pini 2022 and Ward & Pini 2024 decks now use the constant gas heat capacity mode with the paper's value: the printed constant where there is one (30.7 J/(mol K); per-kilogram NIST values converted with each component's molar mass; Ward & Pini 2022's single 840 J/(kg K) given to each component per kilogram, so every CO₂/He mixture has 840 J/(kg K)), and the mean of the ideal-gas Cp0 over 298–373 K for Casas 2012, which states that window but prints no value. Nguyen 2023 also takes the adsorbed-phase heat capacity from the gas, as the paper does. Arora & Hasan 2021 and Nguyen 2022 keep the default. The per-paper table is in `core/test/papers/README.md`.
- **Results move** slightly for these decks (at most 0.06 pp in a purity or recovery); every validation bound holds. The example decks and the example scripts are updated.

### Added — packed-bed film correlations
- Each packed layer has one gas-particle film correlation, `layer[ID].film_correlation` (Python `layer.film_correlation`, builder `set_film_correlation`): `wakao_funazkri` (default), `gunn` (Gunn 1978) or `ranz_marshall` (Ranz and Marshall 1952). By the heat-mass analogy the same correlation gives the particle Nusselt number of the Film gas-solid heat transfer model and the Sherwood number of the Haynes rate coefficient's film resistance. The default is the Sherwood number used before, so existing decks give the same results bit for bit. Ranz-Marshall is a single-sphere correlation and the run log warns when a layer selects it. Sources, validity ranges and the points still to be confirmed are in docs/transport-correlations.md.

### Changed — deck format
- `MixedAdsorbentLayer` is at cereal class version 6 (the film correlation); older decks load with the default, Wakao-Funazkri. The packed path of the Film gas-solid heat model now evaluates the Wakao-Kaguei-Funazkri cube root as the Sherwood number always has (`pow(Pr, 0.3333333)`), a relative change below 1e-7.

### Added — warm-started optimisation and sensitivity studies
- `CycleStateCache` and `RunWithStateCache` (`adsim/solver/state_cache.h`): a cache of cyclic steady states keyed by the decision vector, owned by the caller. Each run starts from the state of the **nearest** point evaluated so far (Euclidean distance in the variables scaled to [0, 1] by their bounds), through the ordinary warm-start path (grid remap, stretched histories). A state that does not fit the deck (another structure) is skipped with a note, a warm-started run that fails is repeated cold, and only runs that reached CSS are cached, so decks without a CSS check and fixed-cycle studies run cold as before. Python: `skarstrom.StateCache` and `skarstrom.optimisation.run_warm`.
- `OptimizeSimulation` (`adsim/solver/optimisation.h`) and `skarstrom.optimisation.optimise`: the whole single-objective loop over registry paths (`GlobalOptimizer`, MaxLIPO with a trust region), the cost from the candidates' KPIs. Candidates are warm-started from the cache by default (`warm_start=False` turns it off); the optimum is re-run from the deck's initial conditions and that cold value is the reported cost, next to the search's own value (`search_cost`). Every evaluation reports whether it was warm, its cycles and wall time. `GlobalOptimizer` also takes plain bounds and an objective of the candidate vector, and exposes `CurrentPoint()` and `GetBestValue()`.
- Sensitivity studies warm-start each run from the nearest earlier run's CSS, keyed by the run's probability coordinates (`sa::EvaluateSample` overload with a cache; Python `Study(warm_start=True)` by default, one cache per worker process; CLI `skarstrom sensitivity --no-warm-start`). Samples report `warm`, `cycles_run`, `css_reached` and `note`.
- Measured savings and the agreement with cold runs are in the cyclic-steady-state documentation (*Warm-started optimisation and studies*). Runs outside optimisation and studies are unchanged.

### Changed — the gas heat capacity is the ideal-gas Cp0(T) at the local temperature
- The default gas heat capacity is now each component's ideal-gas heat capacity C_p0,i(T), evaluated at the local gas temperature in every cell. It comes from the NASA 7-coefficient polynomials of McBride, Gordon and Reno (1993, NASA TM-4513), with the source of every species in `ideal_gas_heat_capacity.h`. It replaces CoolProp's value evaluated once at the initial-condition pressure and temperature and held for the run. The initial temperature is now only a start: editing it no longer changes the cyclic steady state (#145), and a TSA cycle's gas heat capacity follows its temperature.
- The polynomials cover 33 of the engine's 34 components; Air is built from N2, O2 and Ar at CoolProp's air composition. Each one agrees with CoolProp's ideal-gas Cp0 within 1 % (light and permanent gases) or 2 % (others) from 200 to 800 K; the table is in the Equation of State page. Hexane, which has no TM-4513 entry, keeps the held CoolProp value, and a run-start warning says so.
- **Results move** for every deck in the default mode. The example decks are re-saved.
- The registry name of the mode is now `ideal_gas` (`fluid.gas_heat_capacity_mode`: `ideal_gas`, `constant`); `coolprop` is still accepted as its former name. In C++, `GasHeatCapacityMode::IDEAL_GAS`, with `COOLPROP` kept as an alias.
- The *Gas phase* adsorbed-phase heat capacity follows the adsorbed phase's temperature. Kirchhoff heats of reaction integrate the heat capacities exactly.
- `Fluid::MolarHeatCapacity`, `MixtureHeatCapacity` and `AdsorbedPhaseHeatCapacity` take the temperature; new `MolarEnthalpyChange`, `HeatCapacityPolynomial` and `ComponentsWithHeldHeatCapacity`. Python: `skarstrom.ideal_gas_heat_capacity(component, T)`.
- `CSSStartIndependence` gains hot-deck and cold-deck starts: the deck's own initial temperatures edited by ±40 K.

### Added — estimated distance to cyclic steady state
- Every CSS check now estimates how far the run still is from the cyclic steady state: the convergence rate ρ from the decay of the last few end-of-cycle state changes, and the remaining distance ≈ last change × ρ/(1−ρ), in the norm of `css.state_tolerance` and as a multiple of it. It is reported in the solver's CSS log lines, as `distance_estimated`, `convergence_rate`, `distance`, `distance_ratio` and `distance_note` on `sim.css_history()` records, in the `css` entries of `skarstrom run --json` and in `css.csv`. The cycles right after an acceleration jump, and windows that are not a monotone decay, report "not estimated" with the reason.
- A run that stops at CSS more than 100 state tolerances from it (by this estimate) logs a warning. The stopping rule is unchanged: no result moves.

### Added — a layer's mass fractions can be expressions
- A layer's mass fractions are expressions over the deck variables, like other deck values: a two-solid blend can be written as `x` and `1 - x` with a deck variable `x`, and a sweep or a sensitivity study varies `variable[x]`. The registry path `layer[ID].mass_fraction[ID]` is now an Expression parameter, so `--set`, `sim[...] = "1 - x"` and `layer.set_composition({"13X": "x", "Glass beads": "1 - x"})` accept text. The composition check uses the evaluated fractions, shows the formulas with their values when the sum is off, and reports a fraction whose formula does not evaluate. The blend edits (`AddSolidToLayer`, `SetLayerMassFraction`, `RemoveSolidFromLayer`) refuse a layer whose fractions include a formula, with nothing changed.
- `MixedAdsorbentLayer` archive version 5: the fractions are written as expressions. Older decks load with their numbers unchanged, bit for bit.

## [0.8.0] — 2026-11-07

### Added — CSS stops verified against a long reference
- A slow-tier test (`CSSReference`) checks plain, Aitken and Anderson stops
  and a perturbed warm start against a long reference run, on eight example
  decks. It compares the full end-of-cycle state (loadings, mole fractions,
  temperatures and pressure), the endpoint histories and the KPIs. Every
  distance must be within c times the stop's state tolerance, with one bound
  (c <= 100) for all decks and methods. The measured c is at most 37.
- `docs/cyclic-steady-state.md` tabulates c per deck and method, gives the
  reference drift, and lists the decks left out and why.

### Added — CSS start independence verified
- A slow-tier test (`CSSStartIndependence`) runs eight example decks to CSS
  from the shipped start, a clean bed, a feed-saturated bed, a hot and a
  cold bed (+-40 K, through the warm-start API) and, where one exists, a
  neighbouring deck's CSS. Every start must reach the shipped start's CSS
  within c <= 100 times the state tolerance, in the full state, the
  histories and the KPIs. The measured c is at most 74 (Perez 2019 E), and
  no deck has a second steady state.
- `docs/cyclic-steady-state.md` tabulates the result, and explains why the
  hot and cold starts do not edit the initial temperature: with CoolProp
  gas heat capacities that also changes the model for the whole run.

### Added — grid sequencing benchmark
- A slow-tier test (`GridSequencing`) compares three ways of reaching a
  fine-grid cyclic steady state on six example decks: the fine grid from
  cold, a coarse grid to CSS followed by a warm-started fine grid, and three
  levels (N/4, N/2, N cells). It asserts that every path reaches the same CSS
  within the CSS reference bound and that the warm-started fine level needs
  fewer cycles than the cold run; timings are printed, not asserted.
- `docs/cyclic-steady-state.md` tabulates cycles and wall time per path and
  recommends when sequencing pays: 1.4-2.1x on the slower VSA decks, nothing
  on a fast TSA deck or the PSA deck; a grid-refinement study should always
  warm-start each level from the previous level's CSS.

### Changed — TSA tutorial heated by a heating medium
- The illustrative TSA deck (`Illustrative_TSA_13X`) and its tutorial now
  heat with a 150 C medium (steam or hot oil) and cool with 25 C water.
  Both pass through a jacket on the wall and an internal heat exchanger in
  the packing, replacing the 30 W heater and the -12 W cooler. The deck's
  results change: CSS at cycle 27, CO2 purity 86.3 %, recovery 96.0 %,
  thermal energy 10.9 MJ/kg, all of it from the medium duty.
- The tutorial explains the choice, with a table showing that a 100 C
  jacket alone is limited by the bed-to-wall heat transfer on this column.

### Added — heating-medium verification
- Analytic heat-up tests against closed forms: the jacket heating the wall
  alone, the jacket heating wall and bed (exact two-node solution), and the
  internal exchanger heating the bed alone.
- `sk.to_script` writes a step's heating medium (its temperature, jacket
  coefficient and internal exchanger). A ramped medium temperature is left
  as a comment, since the builder cannot create the ramp yet.

### Added — inert support solids in the adsorbent library
- Five non-adsorbing support solids for structured contactors, each a
  library material of the new family `support`: cordierite (dense crystal),
  pure aluminium, stainless steel 304L and 316L, and alumina-silica ceramic
  fibre paper. Mix one by mass with a sorbent in a monolith or laminate
  layer's composition (`sim.add_library_adsorbent("cordierite")`). Every
  value is cited (USGS Bulletin 2131 and Open-File Report 88-690, Ho, Powell
  and Liley 1972, ANL-75-55, Headley et al. 2019); the ceramic paper's
  thermal conductivity has no citable value at adsorption temperatures and
  keeps the engine default, as its entry says.
- Library catalogue: a property may carry a `temperature_table` (values
  against temperature, transcribed from the source); the engine uses the
  constant value at 300 K. Support entries need density, heat capacity and
  thermal conductivity instead of a particle size, and their particle
  voidage is 0 by convention (provenance `convention`).

### Added — heat in monolith and laminate layers
- Heat-dispersion model `channel` for structured layers: the gas conducts over the open frontal area, λ_f = OFA·k_g, under LTNE; under LTE the single temperature sees OFA·k_g + (1 − OFA)·k_s.
- Axial conduction through the solid of a structured layer under LTNE, λ_s = (1 − OFA)·k_s, with the face conductance of the two half cells in series. It is zero on faces that touch a packed layer or a void and at the bed ends; packed beds are unchanged.
- Verification tests `StructuredHeat.*`: wall-to-solid and gas-to-solid relaxation against the analytic exponentials, steady Fourier conduction through the solid, energy closure, and packed beds free of solid conduction.

### Changed — structured layers take the column-wall flux through the solid
- Under LTNE, the whole bed-to-wall flux of a structured layer goes into the solid; the gas share is 0 instead of the open frontal area. Packed layers keep the ε_b split.
- The deck checks refuse Ruthven and Yagi–Kunii heat dispersion and Dixon and Leva wall coefficients on a structured layer, and the `channel` model on a packed one. The structured-layer warning about these closures is gone; the LTE warning stays.
- `build::SetLayerKind` to a structured kind installs `channel` heat dispersion and a constant h_i (replacing Dixon or Leva; a constant coefficient keeps its value). Switching back to packed restores constant heat dispersion.
- New run-start warning: a structured layer whose constant h_i is 0 while the column wall exchanges heat (h_o > 0 or a correlation, a fixed-temperature wall end, or a heating medium) is cut off from the wall; the warning names `layer[ID].internal_htc`.

### Added — structured layers: monolith and laminate geometry
- A layer can now be a monolith (square or circular channels) or a laminate
  (sheets with a gap) as well as a packed bed: `MixedAdsorbentLayer::kind_`
  with `StructuredGeometry` (cell density, wall thickness, sheet thickness,
  gap, washcoat thickness). Builder: `build::SetLayerKind`,
  `build::SetChannelShape`, `build::LayerKinds`, `build::ChannelShapeKinds`.
- The flow in a structured layer is laminar Hagen-Poiseuille,
  `K = OFA D_h^2 / (2 fRe)` with fRe 14.23 (square), 16 (circle), 24 (plates),
  and no inertial term; the layer's momentum model is not used. The open
  frontal area takes the bed void fraction's place (gas holdup, bulk density,
  dispersion area, heat split).
- Derived geometry for front ends: `LayerGeometrySummary` (open fraction,
  hydraulic diameter, specific surface, diffusion length, channel count,
  washcoat fraction).
- Parameter registry: a structured layer lists its kind's fields
  (`channel_shape`, `cell_density`, `wall_thickness` for a monolith;
  `sheet_thickness`, `gap` for a laminate; `washcoat_thickness`) and keeps
  `particle_diameter` while the packed-bed transfer closures still read it;
  `bed_void_fraction` and `momentum_model` are not listed for it.
- Deck checks refuse a geometry without an open channel. Run-start warnings:
  LTE with a structured layer (LTNE recommended), and the packed-bed closures
  a structured layer still uses.
- `MixedAdsorbentLayer` is serialised at version 4 (`Kind`,
  `Structured_Geometry`), `StructuredGeometry` at version 1. Older decks load
  as packed beds with unchanged results.

### Removed
- The unused `StructuredLayer` class.

### Added — structured layers: channel heat and mass transfer
- Reynolds number of a structured layer on the channel velocity and hydraulic
  diameter, `rho (|u| / OFA) D_h / mu`.
- Channel film coefficients: the fully developed laminar Nu = Sh of Shah and
  London (2.976 square, 3.657 circle, 7.541 plates), or Hawthorn's
  developing-flow form over the layer length with the new per-layer switch
  `layer[ID].entrance_effects` (`StructuredGeometry` v2).
- The macropore rate coefficients use Glueckauf's slab LDF `3 D_eff / l^2`
  (half the wall or sheet, or the washcoat) on a structured layer; Haynes uses
  the slab's film and pore resistances with the channel Sh. Packed beds keep
  `60 D_eff / d_p^2`, unchanged.
- New gas-solid model `film` (bed-wide): Wakao-Funazkri on packed layers,
  the channel Nu on structured layers. The constant model uses the layer's
  specific surface (`6 (1 - eb) / dp` for packed beds, unchanged).
- New dispersion model `taylor_aris`: `D_m + u_ch^2 D_h^2 / (kappa D_m)`,
  kappa 192 (circle; square as an approximation) and 840 (plates).
- `build::SetLayerKind` installs Taylor-Aris and `film` for a structured
  layer and restores the packed defaults when switching back.
- Deck checks reject Ruthven dispersion on a structured layer and Taylor-Aris
  on a packed one.
- Operating limits: a warning when a structured layer's channel Reynolds
  number passes 2000 (`settings.warn_channel_reynolds`, `WarningSettings` v1).
- Verification tests `StructuredTransfer.*` against each source's own result
  and a parallel-passage breakthrough against the LDF moment analysis.

### Added — packed and structured layers give the same breakthrough
- Verification tests `StructuredMetamorphic.*`: a packed Kozeny–Carman layer, a laminate and a square monolith with the same open fraction (bed void), sorbent mass, constant dispersion, constant LDF coefficient and matched permeability give the same breakthrough curve, isothermal (Henry) and non-isothermal under LTE (Langmuir with a heat of adsorption and wall exchange).

## [0.7.0] — 2026-11-06

### Added — Anderson acceleration of the approach to cyclic steady state
- New CSS acceleration method *Anderson* (`AccelerationMethod::ANDERSON`).
  It moves the slow fields (loadings, solid and wall temperatures, mole
  fractions) along a combination of up to four earlier cycle-to-cycle steps,
  so it can remove several slowly decaying modes in one jump, where Aitken
  removes one. It uses the same projection and the same re-seed check as
  Aitken. Measured CSS stops (plain / Aitken / Anderson): Perez 2019 E
  205 / 65 / 45, Nguyen 2023 X over 800 / 548 / 190, Perez 2019 A 48 / 32 / 35,
  at KPIs as close to a long reference run as the plain stop's. Decks that do
  not select it are unchanged.

### Changed — hint for varying a blend fraction
- The error for layer mass fractions that do not sum to 1 now also says how to
  vary one solid's fraction in a study: change the others with it (for two
  solids, x and 1 - x together at each point). A mass fraction is a number, not
  an expression. `docs/model-conventions.md` gives a Python loop and the CLI
  form, and says that a sensitivity study cannot vary a blend fraction yet.

### Changed — cycle counts in the CSS acceleration docs
- The Numerics, Cyclic Steady State and Perez 2019 tutorial pages now give
  the stopping cycles and KPI differences of the current engine. Among them:
  Nguyen 2023 X on 20 cells reaches CSS at cycle 84 with Aitken instead of
  349, and the accelerated test cycle stops 28 cycles before the plain rule.
  The counts for a clean bed and for the old upwind 30-cell deck are removed.

### Changed — faster builds: only the dlib sources the engine uses
- The build compiles 6 of dlib's ~50 sources: the ones behind
  `dlib::find_min_global` (the global optimiser, its thread pool and dlib's
  build-consistency check), not the
  whole library. A cold dlib build needs about a seventh of the CPU time it did.
  The `dlib::dlib` target, the installed `Adsim::dlib` export and the
  installed dlib headers are unchanged (`cmake/dlib_minimal.cmake`).

### Changed — shorter pull-request CI: the JS-exceptions Emscripten build runs nightly
- The "Engine with Emscripten" CI job builds the engine with native
  WebAssembly exceptions (the default) on every pull request. The second build,
  with JavaScript-based exceptions (`ADSIM_EMSCRIPTEN_EXCEPTIONS=js`), now runs
  only on pushes to main, the nightly run and manual runs.

### Changed — face temperature follows the advection scheme
- The energy balance's convective face temperature is now reconstructed with
  the selected advection scheme (WENO or van Leer), with the same stencils and
  the same upwind faces (bed ends, void cells) as the mole fractions. It used
  to be the upwind value whatever the scheme, which kept the temperature first
  order. In the manufactured-solution test the temperature now converges at
  order 2 with WENO and van Leer, and the other fields of a thermal run keep
  their second order on fine grids. Upwind runs are unchanged. Results of
  thermal runs with WENO (the default) or van Leer move; the re-saved example
  decks and the paper comparisons record by how much.

### Changed — energy KPI documentation
- The performance metrics page lists the three energy figures (electrical
  energy, thermal energy, cooling duty) in one table and shows how to convert
  J/kg to MJ/kg, GJ/t and kWh/t, and says that the thermal energy and cooling
  duty count the energy sources only: a heating medium's duty is not logged
  yet. The TSA tutorial reads the heater's thermal
  energy (6.07 MJ per kg of CO₂) and the cooling duty from the KPIs instead
  of working them out by hand, and changes the cooler power through the
  parameter registry.
- A closed-bed test checks that the logged heat duty of a step equals the
  analytic power × time, for heating and for cooling. Results do not change.

### Fixed — memory grows without bound when a step takes very many solver steps
- The result log records a full snapshot of every variable at every accepted
  solver step (about 100 kB per solver step at 320 cells), so a stalled or
  very stiff step could exhaust memory before it failed. Each (cycle, step)
  now has a memory budget of 2 GiB, sized from the row size: when a step
  reaches it, every other interior row of that step is dropped (the first
  row and the newest row are always kept, so the end-of-step state used by
  the CSS check is exact), and the solver log says so once at the LOW level.
  The budget is above every bundled deck (the largest, Streb2021_A1, logs
  about 1.07 GB in its one step), so no result changes.
- The `max_log_points_per_step` setting is now documented as what it does:
  it decimates the profile CSV files and adds fill-in points in long solver
  steps; it does not limit the in-memory log.

### Fixed — GAB fit set-up and its test
- The default fit set-up for the GAB kinds keeps K between 0.3 and 1.2 and
  the monolayer capacity below ten times the largest loading, the physical
  range for water, so the search stays out of the near-linear corner
  (K -> 0 with a very large q_m) where the parameters are not identifiable.
- The `gab_young` fit test used a non-physical target (K = 0.14 at 298 K),
  where the fit is ill-conditioned and ended in a different corner on macOS
  arm64; it now uses the consistent Young et al. (2021) Lewatit set, and a
  new test checks that the GAB fits recover the same parameters for every
  global-search budget and a last-digit perturbation of the data.

### Added — heating-medium duty
- Each step's log records the heat its heating medium passed into the wall
  and took back out over the step (`StepHeatDuty::medium_heat_in` and
  `medium_heat_out`, J). It is integrated with the trapezoidal rule over
  every logged point, separately from the energy-source powers.
- The thermal-energy and cooling-duty KPIs, total and per step, now include
  the medium's heat in and heat out. `CalculateStepThermalDuty` also reports
  the medium's share (`medium_heating`, `medium_cooling`). Results do not
  change for decks without a medium.

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

### Changed — layer compositions are checked
- **Behaviour change:** a deck whose layer mass fractions do not sum to 1 now
  fails the deck checks instead of running with scaled properties. The
  message names the layer, the sum and the fix, e.g. `Layer "Bottom": mass
  fractions sum to 0.8; they must sum to 1`. No bundled deck is affected.
- A layer whose mass fractions do not sum to 1 (within 1e-6), or that holds a
  negative fraction, is now refused by the deck checks, with a message naming
  the layer and the sum. Every property of a layer's solid scales with the
  sum, so a sum of 0.8 used to mean 80 % of the solid without notice. The
  fractions are never renormalised.
- A solid at mass fraction 0 is allowed but reported as a warning (it takes no
  part in the run).
- The rules live in one engine function, `CheckLayerComposition`
  (`porous_media.h`), which the deck checks (C++ runs, Python `sim.check()`),
  `LayerCompositionWarnings` (`deck_warnings.h`) and front ends share.
- `RemoveAdsorbentFromDeck` refuses to remove a solid that is part of a blend
  (a layer that also holds another solid) and says which layers use it, like
  removing a component that a reaction uses. A layer that holds only that
  solid is still re-seeded with another one. New checked overload
  `RemoveAdsorbentFromDeck(data, id, err)` and `AdsorbentInUseByMixtures`.
- New blend edits in `deck_edit.h` that keep a layer's fractions summing to 1
  (the other solids are rescaled in proportion): `AddSolidToLayer`,
  `SetLayerMassFraction`, `RemoveSolidFromLayer`. `SeedLayerComposition` no
  longer adds a second solid at 1.0 to a layer that already holds one.
- The solver logs the composition warnings at the start of a run, with the
  other deck warnings.

### Added — verification of mixed-layer edges
- `MixedLayers.*`: a solid at mass fraction 0 changes nothing; a layer of an
  inert solid alone holds nothing and balances; splitting a layer into two of
  the same blend gives the same solution bit for bit; the mole balance closes
  across the interface between two layers of different blends.

### Changed — faster test builds: a precompiled header for the test target
- `AdsorptionSimulatorTests` gets its own precompiled header (the engine's
  `core/src/pch.h` plus GoogleTest), so its ~130 sources no longer parse
  cereal, Eigen, CoolProp, SUNDIALS and GoogleTest from scratch.
  `_USE_MATH_DEFINES` becomes a definition on the test target, since a
  `#define` above `<cmath>` in a source file no longer reaches a header that
  the precompiled header already includes.

### Added — thermal energy and cooling duty KPIs
- Two KPIs from the step energy sources, per kg of net product (J/kg),
  defined and undefined exactly as the electrical energy:
  `thermal_energy`, the heat supplied over the cycle (sources with a
  positive power), and `cooling_duty`, the heat removed (negative power),
  reported as a positive magnitude and not subtracted from the thermal
  energy. Each `KPIRecord` also has `step_thermal_energy` and
  `step_cooling_duty`, every step's share in the deck's step order. They are
  in `KPIRecord` (C++ and Python), `as_dict()`, `kpi_table()`, the CLI's
  `--json` and its printed table, two columns appended to `kpis.csv`, the
  sensitivity outputs (`Component.thermal_energy`, `Component.cooling_duty`)
  and the `KPIType` / `KPIVariableType` enumerators `THERMAL_ENERGY` and
  `COOLING_DUTY`, appended so no stored index moves.
- They use the power the solver applied in each step, logged once per step
  with the step's log (`Log::heat_duty`, evicted with its cycle): a disabled
  source, or one left on a removed layer, counts nothing, and both are 0 in an
  isothermal run. Results do not change.

### Fixed — Ward 2022 outlet-flow label and paper README links
- The outlet volumetric flow files of the Ward & Pini (2022) breakthrough data
  (`core/data/ward_2022/Qin=*/outlet_volumetric_flowrate_*.csv`) now carry
  the unit of their values, `[m3/s]`, instead of `[ccm]`; the values are
  unchanged. The bundled `Ward2022_Q50ccm`, `Ward2022_Q100ccm` and
  `Ward2022_Q200ccm` decks show the corrected label. Results do not change.
- The data README records the source figure as "Fig. 2 or 3, to be checked"
  rather than asserting either, and the Nguyen 2022 test README no longer
  links a `reference_values.csv` that does not exist (that suite compares
  only with digitised curves).

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

## [0.6.0] — 2026-11-06

### Added — regression test for CSS acceleration on a thermal cycle
- A slow-tier test checks that Aitken acceleration shortens the approach to
  cyclic steady state on the thermal four-step LPP cycle (87 cycles without
  it, 39 with it) and that both runs stop at the same KPIs. The jump limit
  that replaced the `|mu| <= 5` cap already gives the reduction; earlier,
  both runs took the same number of cycles.

### Fixed — library re-saves under `ctest --parallel`
- The tests that re-save `data/adsorbents/*.json` now write through the
  suite's atomic-write helper (`test::SaveAdsorbentToJsonAtomically`): the
  rename over the target is retried while a reader holds the file open
  (Windows), and a failed write throws instead of leaving the old file in
  place unnoticed. The library files are unchanged.

### Fixed — flaky tests reading example decks under `ctest --parallel`
- The tests that re-save `examples/decks/*.json` now write each deck to a
  temporary file in the same directory and rename it over the deck, so a test
  loading the deck at the same time always reads a complete file instead of
  truncated JSON (`rapidjson internal assertion failure: IsObject()`). On
  Windows, where the replace fails while a reader holds the file open, the
  rename is retried for up to about 5 s. The deck contents are unchanged.

### Fixed — class versions visible wherever a deck is saved or loaded
- Every serialised engine class now declares its `CEREAL_CLASS_VERSION` in the
  header that holds its inline `serialize` template, not in its `.cpp`. A
  translation unit that saved or loaded a deck without seeing the declaration
  used version 0, so which version a save wrote (and which branch a load took)
  depended on which instantiation the linker kept. Saved decks do not change.
  A new unit test (`CerealClassVersions.*`) pins the version each class writes.

### Fixed — class versions of the fluid, connections and CSS acceleration
- `adsim::Fluid`, `adsim::Connection` and `adsim::AccelerateCyclicSteadyState`
  now declare their `CEREAL_CLASS_VERSION` in their headers too, so every
  serialised engine class does. Saved decks do not change; the
  `CerealClassVersions.*` test now pins these three classes.

### Changed — changelog entries are fragments in `changelog.d/`
- A pull request adds its changelog entry as a new file
  `changelog.d/<topic>.md` instead of editing `CHANGELOG.md`, so pull
  requests in flight together no longer conflict on the top of *Unreleased*
  (GitHub ignores the union merge driver). `scripts/release.py` reads the
  fragments as part of *Unreleased*; the release pull request moves them into
  `CHANGELOG.md` and deletes them. New commands: `release.py fragments`
  (check) and `release.py collect`. CONTRIBUTING.md, RELEASING.md and the pull
  request template describe it.

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

### Added — electrical energy KPI
- The specific energy KPI is now also reported under its explicit name,
  `electrical_energy` (compressor and vacuum-pump work over the cycle per kg
  of net product, J/kg): a `KPIRecord` field and Python attribute, a key of
  `KPIRecord.as_dict()` / `kpi_table()` / the CLI's `--json` KPIs, a
  sensitivity output (`Component.electrical_energy`), an appended
  `electrical energy [J/kg]` column of the CLI's `kpis.csv`, and the
  `ELECTRICAL_ENERGY` enumerator of `KPIType` and `KPIVariableType` (the
  same enumerator as `ENERGY_CONSUMPTION`). The KPI's display label is now
  "Electrical Energy", and `skarstrom run` prints it as
  `electrical energy [J/kg]`.

### Changed — `energy_consumption` is documented as the electrical energy
- `energy_consumption` keeps its name, value and meaning in KPI records,
  sensitivity outputs (still the default), optimiser objectives, the cost
  estimate and the existing `kpis.csv` columns; it is documented as the
  specific electrical energy. Step energy sources (heaters, coolers) are
  still not included. Results do not change.

### Fixed — extended Langmuir–Freundlich trace loadings
- `extended_langmuir_freundlich` no longer multiplies each competitive term by
  the smooth step y² / (y² + 1e-10) (with a softplus on y). That step biased
  trace loadings by about 1e-10 / y²: 1e-4 at y = 1e-3, 1 % at 1e-4, about
  50 % at 1e-5. The power (P y / P_ref)^n is now continued below y = 1e-5 by
  the same C1 quadratic the competitive Sips kinds use, so the isotherm is
  exact for y ≥ 1e-5 in both the solver and the single-component evaluation,
  with a finite Henry slope below and exactly zero loading at y = 0. Results
  of runs with trace components on this kind move; no bundled deck uses it.

### Fixed — energy-balance docs equation
- The LTE energy-balance equation in the docs now matches the engine: the adsorption heat term no longer carries a stray $RT/P$ factor, its sign follows the stored heat of adsorption (positive for exothermic adsorption, so adsorption heats the bed), and the equation now shows the adsorbed-phase heat capacity in the accumulation term and the step energy source.

### Added — GAB and BET water isotherms
- Two isotherm kinds for water and other condensable vapours, in the relative
  pressure x = p_i / p_sat(T): `gab` (Guggenheim-Anderson-de Boer with
  C = C0 exp(dH_C/RT), K = K0 exp(dH_K/RT), q_m = q_m0 exp(dH_m/RT); BET is
  K0 = 1, dH_K = 0) and `gab_young` (the layer-energy form of Young et al.
  2021 and Low et al. 2023, with the papers' own parameters). The saturation
  pressure comes from CoolProp's saturation ancillary (water: within 0.02 % of
  IAPWS-95). The loading is exact for 0 <= x <= min(1, 0.95/K) and continued
  smoothly beyond (the Henry line for x < 0, the tangent above), so it stays
  finite at and above saturation. Not competitive: water adsorbs on its own.
- `adsim::IsostericHeat(isotherm, component, T, P, y)`: the isosteric heat
  implied by any isotherm's temperature dependence, to help choose the
  constant heat of adsorption the energy balance uses.
- The GAB kinds work under the IAST mixing rule (closed-form spreading
  pressure).

### Fixed — IAST with an isotherm kind it does not know
- An isotherm kind without an IAST form now raises an error instead of being
  treated as inert (which silently gave that component zero loading).

### Added — grid convergence analysis in the engine
- New public header `adsim/solver/grid_convergence.h`: observed order of
  accuracy, Richardson extrapolation, Roache's grid convergence index and the
  cell count for a target error, for a result (a cycle KPI, say) on three or
  more grids at any constant refinement ratio. It reports whether the grids
  are in the asymptotic range: monotone, oscillatory (no order reported),
  divergent, grid independent or too few levels. The cycle-KPI grid
  convergence test now uses it (same numbers).

### Added — dense finite-difference check of the solver Jacobian
- New integration test `JacobianDenseCheck`: at states in the middle of each
  step of real decks, the Jacobian IDA receives (coloured, and the band
  evaluations that learn and audit its pattern) is compared entry by entry
  with a dense forward difference over the whole matrix, with no colouring
  and no band assumption. It covers upwind, WENO and van Leer, thermal and
  isothermal runs, Kozeny-Carman and Ergun momentum, a wall energy balance,
  a breakthrough, a PSA and a TSA cycle, and each deck's grid as well as 10
  and 20 cells. Every entry agrees to rounding (relative error about 2e-16);
  a guard case shows that the comparison detects an incomplete pattern.
- `BandedJacobian::SetEvaluationObserverForTesting`: a read-only test hook
  that sees every assembled Jacobian and its column increments.

### Fixed — macOS builds in CI and releases
- The macOS arm64 CI job and the macOS arm64 wheel build ran on GitHub's `macos-14` runner, which GitHub has retired ("The configuration 'macos-14-arm64-us-default' is not supported"), so the job was cancelled before it started and a release could not build the macOS arm64 wheel. Both now run on `macos-15`. The wheels still target macOS 13.3 and later (`MACOSX_DEPLOYMENT_TARGET`).

### Added — verification of mixed adsorbent layers
- `core/test/validation/test_mixed_layers.cpp` (`MixedLayers.*`): solved tests for
  layers packed with several solids at set mass fractions. A blend of an adsorbent
  and an identical copy equals the pure layer, and a blend with an inert solid
  equals the pure adsorbent at the blend's density, voidage, conductivity and heat
  capacity, field by field in LTE and LTNE. The stoichiometric breakthrough time of
  a two-adsorbent layer on linear isotherms matches its mass-weighted capacity, and
  the mole and energy balances close on mixed beds. No defect was found.

### Added — documentation for matching another model
- New docs page *Matching Another Model* (`docs/model-conventions.md`): what `solid_density` is and how the bulk density follows from it, the total voidage and how to set it equal to the bed voidage, the heat-of-adsorption sign with a worked example, which composition the capacity-factor rate uses, the formula of every boundary-value profile (including the pressure-ramp start value, ramp time and exponential factor), how velocity feeds are converted to moles, where profile values sit, and a checklist of the settings that most often differ between codes. Linked from the mass- and energy-balance pages.

### Added — GAB isotherms from Python for any vapour
- `skarstrom.evaluate_isotherm`, `fit_isotherm` (`component=`) and
  `evaluate_iast` (`components=`) name the vapour whose saturation pressure
  the GAB kinds use (a CoolProp name; default water). `isotherm_kinds()` now
  lists `gab` and `gab_young`, and the fitter has default set-ups for both
  (`gab`: q_m0, C0, K0, with dH_C and dH_K at several temperatures;
  `gab_young`: q_m, C_E, F, with G at several temperatures). In C++,
  `isofit::FitOptions::component` and `isofit::Evaluate(..., component)`.

### Added — reactions on a catalyst solid (data model and API)
- A solid (adsorbent) can carry heterogeneous reactions
  (`Adsorbent::reaction_source_`, `adsim/objects/reaction.h`): stoichiometry,
  a rate law per kg of the solid, a rate basis (partial pressure in bar, or
  concentration), a heat of reaction (constant, or Kirchhoff with the gas
  heat capacities the energy balance uses) and an effectiveness factor.
  A catalyst mixed into a layer with a sorbent makes a sorption-enhanced
  reactor. The reactions do not enter the balances of a run yet.
- Rate kinds (`adsim/equations/reaction_kinds.h`): `power_law`,
  `reversible_power_law`, `lhhw` (Langmuir-Hinshelwood-Hougen-Watson with any
  number of adsorption terms; reproduces the Xu and Froment steam-reforming
  rates) and `expression` (a formula of T, P and each component's `p_<C>`,
  `c_<C>`, `y_<C>`, which expressions can now use).
- Element formulas for the built-in components and an element-balance check
  of every reaction.
- Builder: `AddReaction`, `SetReactionKind`, `SetReactionStoichiometry`,
  `SetReactionBasis`, `SetEffectiveness`, `AddAdsorptionTerm`,
  `RemoveAdsorptionTerm`, `RemoveReaction`, `ReactionRateKinds`. Registry
  paths `adsorbent[ID].reaction[ID].*` (category `reaction`).
- Deck checks (`ReactionChecks`) and warnings (`ReactionWarnings`).
- `RemoveComponentFromDeck(data, name, err)` refuses to remove a component a
  reaction uses.

### Changed — adsorbent archive version 5
- `Adsorbent` is saved at cereal version 5 with a `Reaction_Source` node.
  Decks and library files at version 4 load unchanged, with no reactions.
- The library files in `data/adsorbents/` and `examples/decks/Illustrative_TSA_13X.json`
  are re-saved by the test suite at version 5 (version stamp and an empty
  `Reaction_Source` only).

### Fixed — saving an adsorbent file fails silently
- `SaveAdsorbentToJSON` now throws `std::runtime_error` naming the file and
  the reason when it cannot write it (missing or read-only directory, full
  disk, a rename that keeps failing), instead of returning without a file.
  The temporary file is removed and the target is left as it was. The
  rename is retried for about 5 s, as on Windows it fails while another
  process has the target open.

### Changed — shorter energy-balance option labels
- The energy-balance and thermal-equilibrium option labels (`EnergyBalanceModeList`, `ThermalEquilibriumModeList`) drop their bracketed notes: **Non-isothermal**, **Isothermal**, **Local thermal equilibrium** and **Local thermal non-equilibrium**. Only the displayed text changes; the stored values are the same.

### Added — Skarstrom-cycle equilibrium-theory verification
- `test_skarstrom_cycle_equilibrium_theory.cpp` checks a full four-step
  Skarstrom cycle (pressurisation with light product, feed, blowdown, purge)
  on linear isotherms against the Knaebel-Hill equilibrium theory: light-product
  purity 1, the closed-form recovery at two purge-to-feed ratios and the
  heavy-product composition. A slow-tier test refines the grid and the LDF
  rate and shows the gap falling at first order and extrapolating to the theory.

### Added — paper-validation protocol and data provenance
- `core/test/papers/README.md` now holds the seven-step validation protocol every paper goes through (source, provenance audit, model-form audit, data, metrics, attribution, report) and a provenance-table template with a filled example.
- `core/data/<paper>/README.md` for every experimental data set (Nguyen et al. 2022, Streb & Mazzotti 2021, Ward & Pini 2022): source figure, points and data range of each file, how it was obtained (all three are digitised figures) and known caveats, with facts the repository does not record marked "not recorded".

### Fixed — gas heat capacity of condensable components
- In the CoolProp gas heat capacity mode, a component that is condensed in its
  pure state at the initial pressure and temperature (water, ethanol, hexane,
  toluene, ... at 1 bar and 298 K) now takes its ideal-gas heat capacity
  instead of CoolProp's liquid value: water vapour is about 1865 J/(kg·K),
  where the liquid's 4181 J/(kg·K) was used before. Components in a gas or
  supercritical state, including the permanent gases, are unchanged, as is
  the user-defined heat capacity mode. Viscosity and thermal conductivity
  already used vapour-phase fits. No bundled example deck has a condensable
  component, so their results do not change.

## [0.5.0] — 2026-11-05

### Changed — WENO is the default advection scheme; example decks on grid-converged cell counts
- `IDASettings::advection_scheme_` (`settings.advection_scheme`) now defaults to
  `AdvectionScheme::WENO` instead of `Upwind`. A grid study of every example
  deck (10-160 cells, Richardson extrapolation) found that recovery within
  0.1 percentage point of its grid limit needs 200-600 or more upwind cells on
  the cyclic decks but 25-110 WENO cells; a WENO run costs about 1.5x an
  upwind run on the same grid, so at equal accuracy it is 3-10x cheaper (see
  `docs/numerics.md`, "Advection schemes"). Upwind and van Leer remain
  options. Decks keep the scheme they were saved with; a deck saved before the
  setting existed (`Solver_Settings` version < 2) still loads as Upwind.
- The example and paper decks now run WENO on grid-converged cell counts:
  Perez2019 A-F, Nguyen2023 I-X (was 100), Vysyaraju2024 and
  FiveStepProductPurge 80; Ward2024 and Arora2021 PVSA 100; Ferre2026 120;
  Arora2021 PSA, Casas2012 and Streb2021 160 (Ward2022 60, Nguyen2022 and
  Illustrative_TSA_13X 30 unchanged). The Python example scripts build the
  same decks.
- The Perez2019 and Nguyen2023 decks start from a feed-loaded bed (the feed gas,
  15 % CO2, with the loadings in equilibrium with it) instead of a clean N2
  bed. The cyclic steady state is the same (Perez2019_F: 75.65 / 97.69 % after
  517 cycles from the feed-loaded bed, 75.58 / 97.69 % after 2276 from a clean
  one), reached in far fewer cycles: Perez2019 A-F after 49, 57, 65, 94, 207 and
  517 cycles, Nguyen2023 I-X after 37-50 (I-V), 80, 326, 447, 1166 and 1558.
  Cycle caps are about 30 % above those (Perez2019 E 270, F 680; Nguyen2023
  VI-X 110, 430, 590, 1520, 2030; the others keep 250). Nguyen2023 X stays slow
  from any start (slowest mode 0.997 per cycle, from the deck's
  composition-only coupling of the light-product pressurisation).
- Every paper metric was re-measured and the moved bounds re-set from the
  grid-converged errors (34 bounds: loosened where the error grew past the
  bound, tightened where it fell well below it), with before/after tables in
  each paper README ("Grid and advection scheme"). Most VSA recoveries rose by
  0.15-2.6 pp; that moves some comparisons towards the reference (Perez B and
  C1, Ferre, Casas Fig. 4 and 25/75 H2 times, the Ward2022 temperature rises,
  Vysyaraju vs the paper model) and some away (Perez A, D1, E, F recovery,
  Ward2024 recovery, Arora PSA purity, Casas 75/25 H2 time, Ward2022 Q200
  outlet NRMSE).
- The Perez2019 D1, E and F tests also re-run each deck under the paper's own
  stopping rule (each species' balance within 1 % over 5 consecutive cycles, no
  state test) and print the result next to the CSS one, report only: they stop
  after 44, 84 and 152 cycles, with recoveries about 0.1 pp and F's purity
  0.5 pp above the CSS values (`core/test/papers/perez_2019/README.md`).
- Tests whose reference is the upwind discretisation (first-order order of
  accuracy, numerical dispersion u dz/2 in an analytic reference, tanks in
  series, the upwind Jacobian band), or whose tolerance is upwind's
  round-off-level monotonicity, now set `AdvectionScheme::Upwind` explicitly,
  each with the WENO number in a comment; `CSSRobustDeck` keeps the 30-cell,
  clean-start upwind form of Perez2019_A its converged KPIs were measured on.
- Run time, measured on the same quiet 4-core machine before and after: paper
  suites (outside CI) 1016 s summed / 512 s wall at -j2 before, 2850 s /
  1431 s after (Ferre2026 66 -> 826 s, Nguyen2023 743 -> 1295 s, Perez2019
  149 -> 403 s including the paper-protocol re-runs). With the CI tiers:
  the fast tier (`ctest -LE "paper|slow"`, 878 tests) takes 98 s summed / 32 s
  wall at -j3 and the slow tier (25 tests) 163 s / 64 s; the facade and
  sensitivity tests that run Perez2019_A / Arora2021_LiLSX_PSA pin a small grid.

### Fixed — Aitken CSS acceleration stopped with WENO / van Leer
- `AccelerateCyclicSteadyState::Accelerate` skipped every three-cycle window
  whose fitted jump exceeded five natural steps (|mu| > 5, a decay rate above
  0.83). Once a run settled on a slow mode (mu about -28 at a rate of 0.965,
  -140 at 0.993) no window passed again: with WENO or van Leer the
  acceleration test cycle got two jumps and reached CSS at the same cycle as a
  plain run. The jump is now limited to twenty steps instead
  (`AccelerateCyclicSteadyState::kMaxStepMultiple`); a limited jump undershoots
  the fixed point along a contracting mode.
- A window whose coefficient has the wrong sign for a geometric decay (mu > 0
  with steps that point the same way, typically a drifting front) is skipped:
  it used to be applied and moved the state backwards along its trajectory.
- Skipped windows are logged as "CSS acceleration skipped (reason)" (MEDIUM),
  like rejected re-seeds, and each window's fit at HIGH; `Accelerate` takes an
  optional `AccelerationNote*` that carries them.
- Effect: the test cycle reaches CSS in 24-26 cycles instead of 39 with every
  scheme (`CSSAcceleration.ReducesCyclesToCSS{Upwind,Weno,VanLeer}`); Perez
  2019 D1 at cycle 144 instead of 189 (upwind, 30 cells; 179 before the fix)
  and 198 instead of 237 with WENO on 80 cells (250 before); Nguyen 2023 X on
  20 cells (one slow mode of rate 0.985) at 111 instead of 435 (435 before;
  `CSSAcceleration.SlowGeometricModeReachesCssSooner`). Purity and recovery
  agree with the plain run's within 3e-5. `CSSConvergence` is unchanged.

### Changed — fast CI: smaller smoke tests and a slow tier
- Tests are in three tiers, set by ctest labels in `core/test/CMakeLists.txt`:
  the fast tier (no label), `paper` (the paper replications) and `slow`
  (verification that needs fine grids or many cycles: the cycle-KPI study of
  the higher-order schemes, the CSS-acceleration studies
  (`ReducesCyclesToCSS*`, `SlowGeometricModeReachesCssSooner`) and slow-CSS
  studies, the Perez2019_A stopping-rule check, the fine-grid and long-run
  regressions and the cross-code comparisons). Pull-request CI runs the fast tier
  (`ctest -LE "paper|slow"`) on Linux, macOS and Windows; pushes to main, the
  nightly run and manual runs run every tier. Linux and macOS no longer run
  the paper suites on pull requests.
- Smoke and behaviour tests run on smaller grids and fewer cycles: the biogas
  8-step cycle (5 cycles on 25 cells), the four-step LPP family (10 or 20
  cycles; the mass-balance CSS tests isothermal), the map-order round trip
  (15 cells), the coloured-Jacobian checks (30 cells), the robustness fuzz test
  (32 decks; 96 on main and nightly). Verification tests lost only waste: the
  Langmuir shock moment test runs on 30 cells (its first moment is a mass
  balance; the tolerance went from 3 % to 0.5 %), the dispersion CSS balance
  12 cycles instead of 60, the equalisation tests 12 cycles instead of 30, the
  advection "less diffusive" tests two grids instead of three, the cycle
  rotation test 20 cells. The mass-balance CSS tests now count the cycles run
  from the last cycle index (the size of the log's cycle set is capped by
  `max_cycles_to_keep`, so their bound was vacuous) and assert that CSS was
  reached.
- Build: CI sets `CCACHE_SLOPPINESS=pch_defines,time_macros`, so ccache caches
  the engine sources built with the precompiled header (before, every engine
  source was recompiled on every run), and Clang builds the precompiled header
  without its timestamp (`-Xclang -fno-pch-timestamp`), which a cached PCH
  needs.

## [0.4.0] — 2026-11-05

### Added — automatic releases
- Merging a version change into `main` now releases it: `release.yml` runs on
  a push to `main` that changes the version in `CMakeLists.txt`, creates the
  tag `v<version>` itself, uses the version's CHANGELOG section as the
  GitHub Release notes and publishes to PyPI (TestPyPI for pre-releases). Before,
  only a hand-pushed tag released, and 0.2.0 and 0.3.0 were never tagged or
  published. A pushed tag still works. A final version is published only if
  CHANGELOG.md has its section with nothing left under Unreleased, and
  CITATION.cff carries it.
- New workflow `prepare-release.yml` keeps one "Release X.Y.Z" pull request
  open (branch `bot/release`) whenever Unreleased has entries, rebuilt from
  `main` on every changelog change. The next version comes from the Unreleased
  subsections: Added, Changed, Removed or Deprecated give a minor release,
  Fixed or Security alone a patch release; major releases and pre-releases are
  chosen by hand (Run workflow with `bump` / `suffix`).
- `scripts/release.py` (`current`, `next`, `prepare`, `check`, `notes`) does the
  version arithmetic and the three-file release change, with tests in
  `scripts/test_release.py` (run by CI's docs job). The version on `main` is
  now the last released version; it is no longer bumped for development.
  RELEASING.md and CONTRIBUTING.md describe the new flow.

### Added — KPIs of every cycle (KPIHistory)
- The solvers record the KPIs of every tracked component when each cycle
  completes (`adsim::KPIHistory`, `adsim/solver/kpi_history.h`;
  `CoupledSolver::GetKPIHistory`, `CoupledSolverChuncked::GetKPIHistory`), so a
  run's KPIs cover all its cycles, not only the last
  `settings.max_cycles_to_keep` the logs hold. `Simulation::KPIs()` (Python
  `sim.kpis()`, the CLI's `kpis` and `kpis.csv`) now lists every cycle:
  cycles still in the log are computed from it as before, evicted ones come
  from the history. New: `Simulation::KPICycles()`,
  `Simulation::KPIHistoryRecords()`, and `adsim::ComputeCycleKPIs`, the
  per-cycle computation both use. `KPIRecord` moved from `simulation.h` to
  `kpi_history.h` (which `simulation.h` includes).

### Fixed — WENO / van Leer slow when a component is absent
- The WENO and van Leer interior face mole fractions were clamped to
  [0, 1]. A component absent from the feed and the bed holds only round-off
  around zero, and the clamp zeroed the face (and its derivative) of every
  cell marginally below zero: negative round-off could not leave the bed,
  the component drifted to about -6e-10, and the Newton iteration kept
  failing to converge. WENO took 25 to 80 times longer than Upwind on such
  decks, van Leer 4 to 15 times (run with each scheme,
  `UserDataTest.IsothermalMassSource`: Upwind 0.08 s, WENO 6.1 s, van Leer
  0.52 s;
  `ClosedBedEnergyTest.OpenSystemEnergyBalanceConvergesUnderGridRefinement`:
  0.34 s, 27 s, 5.5 s). The face is now bounded to [0, 1] widened to the
  values of the two cells that share it: identical for cells inside [0, 1],
  and round-off is carried along and out of the bed as under Upwind. Both
  tests now take Upwind's steps with either scheme (WENO 0.09 s and 0.64 s,
  van Leer 0.09 s and 0.51 s). Upwind is untouched. With WENO or van Leer, results of decks
  move at the solver-tolerance level (KPIs by up to 1e-6 relative on ten
  example decks), and an absent component no longer leaks a spurious trace
  out of the bed (9e-7 to 1e-24 mol-fraction seconds on `Nguyen2022_Fig4`).
  Test: `AdvectionSchemeRobustness.AnAbsentComponentCostsNoExtraSteps`.

### Added — own ramp time for pressure ramps; exponential factor
- Linear, exponential and log-linear pressure profiles (pressure boundaries
  and valve streams) take an optional ramp time of their own:
  `use_step_time_` (deck field `Use_Step_Time`, default on) keeps the ramp
  over the whole step as before; off, the profile's `ramp_time_` /
  `step_time_` is used and must not exceed the step time
  (`CheckPressureProfile`, run by `Cycle::PerformChecks`). Registry:
  `endpoint[i].pressure_<kind>_use_step_time`, and `pressure_<kind>_time` when
  it is off. Python: `sk.linear(end, time=10)` for a pressure.
- The exponential profile's factor is a parameter: `factor_` (deck field
  `Factor`, registry `<head>_exponential_factor`, Python
  `sk.exponential(..., factor=)`), v(t) = v1 + (v0 - v1) F^(-t / t_r).
  The default 100 is the old fixed value (99 % at the ramp time); F must be
  greater than 1 (`CheckProfileFactor`).
- Class versions: `LinearBoundaryValue` 4, `ExponentialBoundaryValue` 4,
  `LogLinearBoundaryValue` 2. Older decks load with the step time and F = 100,
  so their results are unchanged.

### Added — method-of-manufactured-solutions hook and order test
- Verification hook `UserData::manufactured_sources_`
  (`adsim/equations/manufactured_sources.h`): source terms, functions of z
  and t, for the overall balance, each component balance, each loading's
  rate law, the momentum faces (a pressure-gradient source on Darcy faces
  and interior Ergun faces) and the LTE energy balance, added to the
  residuals in the scaling of the neighbouring physical terms. Testing only:
  not serialised and not reachable from decks, the parameter registry,
  Python or the web. With no source set the results are bit-identical
  (`bench_solver` logs of two example decks are byte-identical to the
  previous residual code) and the residual cost is unchanged.
  `UpdateResiduals` records the integrator time in `UserData::residual_time_`
  for it.
- `core/test/validation/test_manufactured_solutions.cpp`: a smooth,
  time-dependent manufactured solution with dispersion, LDF kinetics on a
  Langmuir isotherm, Kozeny-Carman flow and (thermal case) the energy
  balance, sources derived by hand and checked by finite differences.
  Observed orders on 20-160 cells: upwind 1.0 in every field; WENO and
  van Leer 2 (WENO interior 2.5-2.8 before the second-order central terms
  take over); with the energy balance the temperature is first order for
  every scheme (the convective face temperature is always upwind). About 4 s.

### Fixed — the CSS state test never passed for an inert bed
- For a bed of inert packing (an adsorbent with no isotherm) the adsorbed
  loadings are round-off, 1e-26 to 1e-22 mol/kg of either sign. The state
  test of the mass-balance CSS methods
  (`CyclicSteadyState::EndOfCycleStateChange`) scales the loading change by
  the RMS total loading and skipped the field only when that scale was below
  1e-300, so it divided noise by noise and reported a change of 0.77 to 51
  on a bed whose state repeats: a mass-balance method with the state test on
  (the default) never stopped.
  The loading scale is now floored at the solver's absolute tolerance on the
  loading unknowns (1e-8 mol/kg); a loading below it is not resolved. Every
  adsorbing deck has an RMS loading far above the floor, so its state change
  is unchanged to the last bit. Test:
  `CyclicSteadyState.InertBedReachesCSSUnderTheStateTest` (the shared PSA
  cycle with inert packing stops at the minimum-cycle floor, cycle 4). Found
  by the dynamic-tank CSS test of the verification programme (WP19-E2a).

### Added — convergence verification of cycle results (grid, integrator tolerance, CSS)
- `CycleKpiGridConvergence` (`core/test/validation/test_grid_convergence.cpp`):
  purity, recovery and specific energy of the shared PSA cycle at cyclic
  steady state on 10, 20, 40 and 80 cells, for each advection scheme.
  Convergence is monotone; the observed order (Richardson, three finest grids)
  is 1.55 for upwind (inside the `[1, 2]` its error expansion allows) and 2.4
  to 2.5 for WENO and van Leer; Roache's GCI of the 80-cell upwind KPIs is
  0.06 to 0.12 % and contains the higher-order solutions.
- `IntegratorTolerance` (`test_integrator_tolerance.cpp`): the error falls with
  reltol at slope 0.68 (breakthrough) and 0.77 (cycle state), inside the
  `[0.5, 1]` of BDF error-per-step control, and floors at the absolute
  tolerance; `max_step` changes results by less than the tolerance error. The
  KPIs converge with reltol at slope 0.42, because the boundary flows are
  integrated over the logged points and the flow transient after a pressure
  switch (~2e-5 s) is resolved only by the integrator's own steps.
- `CSSConvergence` (`test_css_convergence.cpp`): the approach to CSS is
  geometric (rate 0.106 isothermal, 0.966 with the energy balance), the KPIs
  follow it, the error the stopping rule leaves is `g Delta rho / (1 - rho)`
  (measured 6.11e-5 in recovery at the default state tolerance, 6.21e-5
  predicted), and Aitken acceleration stops 25 cycles earlier on the same KPIs.
- `scripts/verification/grid_convergence.py` writes the grid study of the
  example decks (observed order, GCI, extrapolated KPIs); the table and the grid
  each studied deck needs are on the Numerics page ("Convergence of the
  results").

### Added — extended Langmuir-Freundlich isotherm (toPSAil form)
- New competitive isotherm kind `extended_langmuir_freundlich`
  (`IsothermType::EXTENDED_LANGMUIR_FREUNDLICH`,
  `ExtendedLangmuirFreundlichIsotherm`): q_i = (k1 + k2 T) z_i / (1 + sum z_j),
  z_i = k3 exp(k4 / T) (P y_i / P_ref)^{n_i}, 1/n_i = k5 + k6 / T (Yang 1987,
  toPSAil's `calcIsothermExtLangFreu.m`; P_ref = 1e5 Pa keeps toPSAil's
  bar-based k3). No existing kind maps it unless k2 = k6 = 0. Registry,
  Python `isotherm_kinds()`, IAST (an exact Sips term) and the docs know it;
  the fitter does not offer it yet. Near y = 0 each term is regularised with a
  step at 1e-5 (loading changed by 1e-10 / y^2: 1e-4 at y = 1e-3).

### Fixed — competitive Sips kinds no longer suppress trace-component loadings
- `extended_sips_2` and `arrhenius_sips` multiplied every competitive term
  z_i = (b_i P y_i)^{s_i} by y_i^2 / (y_i^2 + 1e-6) in the solver only, to
  keep the slope of y^s (s < 1) finite at y = 0. That halved a component's
  equilibrium loading at y = 1e-3, cut it by 1 % at 1e-2 and by 99 % at 1e-4
  (a trace-impurity front on these kinds travelled about twice as fast as the
  published isotherm gives), while `evaluate_isotherm` and the fitter used
  the unregularised form. Below y0 = 1e-5 the power is now continued by the
  quadratic in y that is zero at y = 0 and matches the power's value and
  slope at y0; above y0 the published isotherm is used exactly. The
  continuation is concave (no inflection), has a finite Henry slope, and the
  single-component evaluation (`IsothermFunction::ExtendedSips2`, new, and
  `IsothermFunction::ArrheniusSips`, so `EvalIsotherm`, Python
  `evaluate_isotherm` and the fitter) uses the same function. Documented in
  `docs/isotherms.md` ("Competitive Sips kinds near zero mole fraction"),
  which also gains the missing Arrhenius Sips section.
- Paper cases move little (the fronts' plateaux are unchanged): Casas 2012
  CO2 breakthrough times (y > 0.002) are 0.6-0.8 s later (errors against
  Table 5: 0.92 → 0.73 %, 0.64 → 0.45 %, 9.81 → 9.99 %); every Streb 2021
  metric stays within its bound, the largest change being the ternary CH4
  fronts (B1 2.13 → 2.38 %, B2 13.7 → 14.0 %). The breakthrough decks take
  10-40 % more integrator steps because the front's toe is now a shock down
  to y = 0 (wall time still under 1 s for A3 and Casas); the 20-cycle biogas
  PSA, also from a CO2-free bed, takes 2-3 % more steps and needs no step
  retries. Tests: `CompetitiveSipsTraceLoading.*`,
  `CompetitiveSipsBelowThreshold.*` and
  `CompetitiveSipsTraceBreakthrough.FirstMomentMatchesPublishedIsotherm`.

### Fixed — breakthrough decks no longer stall on fine grids
- `Ward2022_Q50ccm` ran on 10-160 cells but failed on 320 with every
  advection scheme ("too much work": 50 000 steps, 24 000 Newton convergence
  failures, 72 000 Jacobians) as the CO2 front left the column. The
  difference-quotient Jacobian perturbed the pressure by IDA's increment
  floor reltol |P| + abstol (0.1 Pa at 1 bar). A pressure outlet imposes the
  outlet-face velocity as the explicit Ergun root of the end-cell pressure
  drop, which is linear only over A² / (4B), a pressure scale proportional to
  the cell width (0.12 Pa on 60 cells, 0.02 Pa on 320): the 0.1 Pa secant
  was half the true slope there, and the modified Newton iteration diverged
  unless the step was tiny. The pressure columns of the `BandDQ`, `Coloured`
  and `ColouredCheck` Jacobians now use sqrt(eps) max(|P|, |h P'|) floored at
  the pressure's absolute tolerance (1.5e-3 Pa at 1 bar); every other column
  keeps IDA's increment. The deck now runs on 320 cells in 2 169 steps with
  24 convergence failures (3 s), and on 640 cells in 2 047 steps. Results of
  every deck move at integrator-tolerance level, since the Jacobian's
  pressure columns change in the last digits: across the paper validation
  metrics, front times by at most 2e-5 relative and curve RMSEs by at most
  1e-5 absolute, and every paper test still passes.
  `JacobianMode::IdaInternal` (IDA's own Jacobian) is unchanged. Test:
  `FineGridBreakthrough.Ward2022Q50ccmConvergesTo320Cells`.

### Fixed — WENO / van Leer outlet-cell mole-fraction sum drifts during fast outflow
- With `AdvectionScheme::WENO` or `VanLeer` the cell at a bed end that gas
  leaves through could see its mole fractions drift from summing to one
  while the outflow lasted (a valve surge, a blowdown): 5.0e-4 with WENO and
  4.5e-4 with van Leer in `ErgunValveSurgeTest` (Kozeny-Carman), against
  1.1e-5 with upwind. The schemes clamped every face mole fraction to [0, 1], the bed-end faces
  included; those take the upwind cell's value and are not rescaled to its
  sum, so an end cell holding pure nitrogen at y = 1 + 1e-9 sent out a face
  that summed to one instead of the cell's sum, and the cell's departure from
  one was no longer carried out of the bed. Only the reconstructed values are
  clamped now; the upwind faces are left as they are, as under Upwind. The
  same clamp made the logged balances disagree with the solver (the stream
  logs integrate the end cell's value, the solver used the clamped face):
  the DAC charge with WENO had its total mole balance off by 3.45e-2 mol
  (now 2.7e-5, the upwind value) and the outlet cell's sum by 4.0e-5 (now
  9e-15); the temperature purge with van Leer logged 2.4e-5 mol of a
  nitrogen that is in neither the feed nor the bed (now 1.3e-8). The drift
  also reached 1.6e-6 (WENO) in the steep-front breakthrough of
  `AdvectionSchemeConservation` at reltol 1e-8; every scheme now keeps about
  1e-12 there. Upwind is unchanged (its faces were never clamped). Tests:
  `AdvectionSchemeConservation.EveryFaceSumsToItsUpwindCellSum` (the
  discrete identity, face by face); the valve-surge test run with every
  scheme (`ErgunValveSurgeTest.*WENO`, `*VanLeer`), which bounds the outlet
  cell's sum at the end of each feed step by reltol; `UserDataTest.DAC` and
  `TemperaturePurge` run with every scheme (`*WENO`, `*VanLeer`), with their
  logged balances bounded at the upwind level; and the steep-front sum
  check, now at reltol 1e-8 with a 1e-9 bound. Found by the WENO-everywhere
  experiment of the verification programme (WP20-E1).

### Fixed — the ColouredCheck Jacobian diagnostic no longer changes the run
- In `JacobianMode::ColouredCheck` the band evaluation built for the
  comparison learned every coupling it found into the pattern without
  recolouring, and its count was discarded: a coupling the coloured pattern
  missed stayed missing from the coloured Jacobian for the rest of the run,
  no audit ever counted it (`NumMisses()` stayed 0), and the check run took a
  different path from the Coloured run it was meant to observe. The
  comparison now only observes (no learning, and no residual-call,
  out-of-band or pattern counts, and a failed comparison evaluation does not
  fail the step), so a ColouredCheck run is the Coloured run exactly and
  reports a missed coupling until an audit adds it. Test:
  `ColouredJacobian.ColouredCheckObservesTheColouredRun` (it starts from a
  pattern learned from one evaluation, the new testing switch
  `IDASettings::jacobian_learn_once_for_testing_`).
- With the WENO and van Leer schemes the pattern learned in the first cycle
  missed couplings of the columns in the three cells at each end (e.g. the
  first cell's mole fractions on the inlet ghost cell's other mole
  fractions); the audits added them later. Those columns now also take the
  interior pattern, so the coloured Jacobian equals the band DQ Jacobian entry
  for entry for all three schemes (`ColouredJacobian.MatchesBandDQ_*_WENO` /
  `_VanLeer`). Results of WENO / van Leer runs in the default Coloured mode
  move at integrator-tolerance level (they are now bit-identical to BandDQ
  runs); upwind runs are unchanged, bit for bit.

### Added — constant gas heat capacity and viscosity; heat of adsorption as −ΔU + RT
- **`fluid.gas_heat_capacity_mode`**: `coolprop` (default, unchanged) or
  `constant`, with a molar value per component in
  `fluid.gas_heat_capacity[<component>]` [J/(mol·K)] (a new component starts
  at 29.1, 7R/2). The constants replace the CoolProp values that
  `Fluid::Initialise` caches; the molar mixing rule is unchanged, so equal
  values give that heat capacity at every composition. The *Gas phase*
  adsorbed-phase heat capacity follows them.
- **`fluid.viscosity_mode`**: `polynomial` (default, unchanged) or `constant`,
  with one mixture viscosity `fluid.viscosity` [Pa·s] in every cell.
- **`fluid.heat_of_adsorption_convention = "isosteric_rt"`**: the stored heat
  of adsorption is −ΔU, the internal energy of adsorption, and the heat
  released per mole adsorbed is −ΔU + R·T at the local adsorbed-phase
  temperature (LTE bed or LTNE solid temperature). This is the source term
  of Yancy-Caballero et al. (2020) / PEESEgroup PSA, −(ΔU − RT); a constant
  isosteric value fixes it at one temperature.
- With the three switches a deck reproduces models that use fixed gas
  properties (C_p,g = 30.7 J/(mol·K), μ = 1.72e-5 Pa·s in Yancy-Caballero
  2020 and Haghpanah 2013). Defaults are unchanged, so every example deck
  and validation result is unchanged. `Fluid` is class version 4; older
  decks load with the defaults. Non-positive constants fail the deck
  checks. Registry, Python (`sim.fluid.gas_heat_capacity_mode`,
  `gas_heat_capacity[...]`, `viscosity_mode`, `viscosity`) and the script
  export carry the new fields. Docs: `equation-of-state.md` (constant gas
  properties), `energy-balance.md` (heat-of-adsorption conventions). Tests:
  `FluidPropertyModes*` and the `IsostericRT` cases of
  `ComponentEnergySourceLog`.

### Added — valves can read the bed at the end cell
- Solver setting `settings.valve_bed_pressure` (`IDASettings::valve_bed_pressure_`,
  `"Valve_Bed_Pressure"`; `IDASettings` cereal version 7): `face` (the default,
  unchanged behaviour) evaluates the bed side of every valve stream's linear
  law at the bed-end face, `end_cell` at the last cell, as toPSAil does. A
  connected valve then reads its partner's end-cell pressure, which every
  stream and pressure boundary now logs (`StreamData::GetEndCellPressure`).
  Decks from before load with `face`.

### Added — sensitivity studies saved with the deck
- A deck now saves its sensitivity study's set-up (`UserData::sensitivity_study`,
  `Sensitivity_Study` in the file; `UserData` cereal version 4): inputs and their
  distributions, outputs, run policy, seed, and the screening and propagation
  sizes. Results are not saved. Decks from before load with no study.
- `adsim/solver/sensitivity_study.h` holds `Distribution`, `Input`, `RunPolicy`
  (moved from `sensitivity.h`, which includes it; the API is unchanged) and the
  new `Study`. `sensitivity.h` adds `StudyOutputs` and `CheckStudy(deck, study)`.
- Python: `Study.from_deck(sim, **overrides)`, `Study.save_to(sim)` and
  `Study.setup()`; the study's sizes are constructor arguments, the defaults of
  `screen()`, `propagate()` and `Screening.important()`.
- CLI: `skarstrom sensitivity deck.json` without `--input` runs the deck's
  saved study; options given override it.

### Added — cross-code regression against the PEESE model and toPSAil
- `core/test/crosscode/`: twelve single-column breakthrough and purge cases
  (CO2/N2 on 13X with concentration- and pressure-form Langmuir and Henry
  isotherms, 1-50 % feeds, purge, adiabatic, a 25 % Ergun pressure drop,
  CO2/CH4 on LiX, a steel wall with ambient loss) run in the PEESEgroup/PSA
  column model (MIT) and in toPSAil (GPL-3.0; numerical outputs only), frozen
  as data with their decks and replayed by `CrossCodePeese.*` and
  `CrossCodeToPSAil.*` (about 15 s). Bounds come from the measured, attributed
  differences: isothermal cases agree within 0.6 s on t50 at the same grid;
  the adiabatic offsets (up to 14.5 s) are the CoolProp gas heat capacity and
  the constant heat of adsorption. The folder README has the attribution.

### Fixed — a connection replays a bottom pressure boundary's outflow with the right sign
- `PressureBoundaryCondition::Log` did not record the bed end it logged at,
  so a pressure boundary's `StreamData` always said TOP. For a pressure
  boundary at the BOTTOM of its bed, `StreamData::GetOutflow` then returned
  gas leaving the bed as a negative outflow, and a no-valve stream connected
  to it (e.g. a rinse or heavy-reflux inlet fed from a blowdown or evacuation)
  drew that gas out of its own bed instead of taking it in: the connection
  destroyed twice the transferred moles, and with a step-change blowdown IDA
  failed at the start of the receiving step. The log now records the
  boundary's location, as a stream's log does. A TOP pressure boundary was not
  affected, and no bundled deck, cycle template or paper case connects a
  stream to a BOTTOM pressure boundary, so no published result moves. Tests:
  `ConnectedStreamTest.PressureBoundaryDonorLogRecordsItsBedEnd`,
  `ConnectionSign.*`.

### Added — verification: extreme operating corners
- `core/test/validation/test_extreme_conditions.cpp`: very high and very low
  pressure ratios, near-zero flow, very short and very long steps and high
  temperatures each converge or fail cleanly, with finite and physical states.

### Added — verification: thermodynamic bound on the separation energy
- `core/test/validation/test_thermodynamic_bound.cpp`: at cyclic steady state
  the logged pump energy of a vacuum-swing and a pressure-swing cycle exceeds
  the reversible pump work, the bed's exergy change and the minimum work of
  separation of the products, in that order.

### Added — verification: parameter limits and degenerate mixtures
- `core/test/validation/test_parameter_limits.cpp`: the rate at which the
  model reaches its limits (LDF k to infinity and to zero, heat of adsorption
  to zero, wall heat capacity and bed-wall coefficient to infinity, Langmuir b
  to zero, pressure ratio to one, trace feed), and two degenerate mixtures
  (identical components keep their feed ratio; zero affinity is inert).

### Added — verification: metamorphic relations
- `core/test/validation/test_metamorphic_relations.cpp`: equivalent decks must
  give the same solution to within the integrator's measured error, for the
  component order, a packing split into identical layers, a step cut in two,
  inputs written in other units, and the step a cyclic steady state is started
  from (with a connection whose lag changes under the rotation).
  `BreakthroughCase` gains `n_layers`.

### Added — verification: evacuation, pressure equalisation and buffer tank
- `core/test/validation/test_unit_operations.cpp` checks constant-speed
  evacuation (`P = P0 exp(-S t / V a R T)`), two-bed pressure equalisation
  through connected valve steps at the cyclic steady state (mole balance and
  the `exp(-2 cv t / 1e5 V a)` decay of the pressure difference; with a
  Langmuir isotherm, conservation of the two beds' inventory), and a
  buffer-tank fill (it delivers the set percentage of the vented moles), each
  without and with Henry adsorption.

### Added — verification: thermal fronts
- `core/test/validation/test_thermal_fronts.cpp` checks the speed of an
  inert-bed thermal front against the energy-balance jump condition
  (including the gas heat-capacity term), and the steady exponential axial
  temperature profile of a flowing bed losing heat through the wall, with
  first-order convergence to it.

### Added — verification: blowdown composition path
- `core/test/validation/test_blowdown_composition_path.cpp` checks that a
  binary on linear isotherms follows the equilibrium-theory (Knaebel-Hill)
  composition path `y(P)` during a valve blowdown, with the departure from
  the path equal to the LDF lag (it halves when `k` doubles).

### Added — verification: binary competitive Langmuir fronts
- `core/test/validation/test_competitive_langmuir_fronts.cpp` compares a
  binary extended-Langmuir breakthrough (constant separation factor) with
  equilibrium theory: the roll-up plateau of the weaker component from the
  coherence condition, and the times of both fronts from the jump conditions.

### Added — verification: Langmuir desorption simple wave
- `core/test/validation/test_desorption_simple_wave.cpp` checks that a
  Langmuir desorption front spreads as the equilibrium-theory simple wave:
  the outlet times of three levels inside the wave converge at first order to
  the method-of-characteristics times, and their Richardson extrapolation
  agrees to 0.1 %.

### Added — verification: LDF-limited front against the Thomas solution
- `core/test/validation/test_ldf_thomas.cpp` compares the whole outlet curve of
  a Henry-isotherm, LDF-limited, dispersion-free front with the
  Anzelius/Thomas closed form, and checks that the grid error is exactly the
  upwind scheme's: the variance minus `L dx R^2/u^2` is the Thomas variance,
  the curve error falls at first order and, after the closed-form `O(dx)`
  correction, at second order.

### Changed — CI skips the builds for docs-only pull requests
- A new `changes` job lists a pull request's files. The C++ (Linux, macOS,
  Windows), Python wheel, Emscripten and subproject jobs run only if one of
  them is code: anything except Markdown, `docs/`, `mkdocs.yml`,
  `requirements-docs.txt`, `.gitattributes`, `CITATION.cff` and issue
  templates. `docs/python.md` (its example is run by `python/tests`),
  `python/README.md`, `THIRD_PARTY_NOTICES.md`, `LICENSE` and `NOTICE` count as
  code. The docs and workflow-lint jobs always run, and pushes to `main`,
  nightly and manual runs build everything.
- `.gitattributes`: `CHANGELOG.md` uses git's union merge driver, so two
  pull requests adding entries at the top of Unreleased no longer conflict in
  local merges.

### Fixed — connected streams conserve moles when the donor's flow reverses
- A connected no-valve stream (e.g. a velocity-0 equalisation consumer)
  replayed the magnitude of the connected stream's logged flow. When the
  donor's flow reversed (a valve outlet drawing gas back into its bed), the
  backflow reached the receiving step as more inflow and the connection
  created gas. The replay now keeps the sign: what left the donor's bed enters
  the receiver and what flowed back into it is drawn out (`StreamData::GetOutflow`;
  the log records the bed end it was taken at). Test:
  `ConnectedStreamTest.ReplayKeepsTheSignOfTheDonorFlow`.
- Vysyaraju 2024 `PSA_Exp1` was the one case affected: PER took in 29 % more
  gas than PED released, so the conservation-based He recovery read 17 % (product
  He / fed He was 69 %). With the fix the two recovery measures agree (63.1 %),
  He purity is 99.36 %, and the tolerances are re-measured (see the paper's
  README); the deck's step cap is 0.5 s. The case still fails on macOS and
  Windows (it did before the fix too): its single-bed equalisation lags PER's
  pressure by one cycle, and which state the cycle settles in depends on
  round-off. Linux passes.
- CI: one run per pull request (pushes to feature branches no longer start a
  second run), superseded runs are cancelled, the paper replications run on
  pull requests (full suite on Linux and macOS, fast subset on Windows), Windows
  builds with Ninja + ccache (`windows-ninja-release` preset), the Emscripten
  builds use ccache, and ctest's per-test timings are cached so the longest
  tests start first.

### Fixed — face values of logged cell variables on non-uniform meshes
- **`EquationLog::GetFaceVariable` / `GetComponentFaceVariable` interpolate
  with the right weights.** The cell-centre-to-face interpolation weighted
  each neighbour by its own width, `(y_L dx_L + y_R dx_R) / (dx_L + dx_R)`,
  which biases the face value towards the larger cell. It is now
  `(y_L dx_R + y_R dx_L) / (dx_L + dx_R)`, exact for a linear profile, so a
  face between cells of different widths (between layers with different cell
  sizes, or between a void cell and the bed) sits where it should. The two
  boundary faces now take the plain mean of the ghost and end cells. That is
  the face value the boundary conditions impose (`P_face = (P_0 + P_1) / 2`),
  and it matters because the ghost cell copies the width of the first or last
  *bed* cell, not of a void end cell. A uniform mesh is unchanged: every
  example deck, the paper validation and `docs/validation.md` are unaffected.
  Tests: `EquationLogFaces.*` (hand-built log, non-uniform mesh, linear
  profile exact at every interior face) and
  `VoidRegionSolverTest.FaceReconstructedBoundaryMolesMatchTheStreamLogs`
  (with voids, the face-rebuilt boundary moles now equal the stream logs to
  round-off; they were 1e-5 to 2e-5 off).
  **What moves:** the feed/product KPIs (purity, recovery, productivity,
  energy) do not use these face values and do not change. The total
  mass-balance CSS check (`css.method` 2), which reads the boundary-face
  pressure and temperature, changes for beds with void regions. Its error had
  a spurious floor. On Perez 2019 case A with both voids enabled (1e-4 m³,
  void cell 0.124 m against 0.01 m bed cells), the floor was about 3.46e-3
  (tolerance 5e-3) and is now 2.4e-5, the same order as without voids
  (1.9e-5); CSS is still reached at cycle 75. With 3e-4 m³ voids the old
  floor of 1.26e-2 never passed the tolerance (no CSS in 200 cycles); CSS is
  now reached at cycle 198. Face profiles plotted by front ends also move at
  layer and void interfaces.

### Fixed — every cycle step starts the integrator afresh
- **A step can no longer step over a short pulse at its start.** Cycle steps
  ran back to back on one IDA history, so the first internal step of a step
  inherited the size of the previous step's last, `tstop`-clipped, step. The
  boundary condition's algebraic jump is out of the error test
  (`suppress_alg_`) and an implicit step samples the forcing only at its end,
  so a step whose forcing is a pulse shorter than that inherited step was
  integrated as if the pulse were not there. On `Vysyaraju2024_PSA_exp1` the
  PER step replays the PED blowdown (2.8 mol/s decaying in about 0.1 s): when
  BLO happened to end on a 0.2 s remainder instead of a 0.015 s one, PER's
  first step jumped the whole pulse, the bed missed its equalisation gas and
  the run restarted its approach to cyclic steady state (PRES uptake
  0.1497 → 0.284 → −0.002 → 0.071 mol …). Which cycle tripped hung on
  round-off: natively at cycle 3 (and at cycle 48 with a reltol 1e-4 away),
  in the WebAssembly build at cycle 51, so the web solver disagreed with the
  native one at the last cycle (relative error 0.72). The engine before
  #36/#37 tripped too, but early (cycle 6 natively), and had re-converged by
  cycle 60. Now
  `SolverData::RestartForNewStep`, called by the blocking and the chunked
  solver at the start of every cycle step, restarts IDA from a consistent
  initial condition at the new boundary conditions with fresh history and
  `hin = init_step` (falling back to the unchanged end-of-step point if the
  consistent IC fails). Native and WebAssembly runs of the deck now agree to
  3e-5 per cycle, and native runs no longer depend on the reltol's last
  digits. The run's IDA statistics stay cumulative across re-initialisations
  (they used to restart from zero at each retry restart). Tests:
  `StepBoundaryRestart.*`. Costs about 6 % run time on the deck. Other
  decks' endpoint totals move by at most 3.4e-4 relative (Nguyen 2023 VII);
  Vysyaraju 2024's PRES uptake at cycle 60 by 1.2e-3.

### Fixed — gas-mixture heat capacity mixed on the molar basis
- **The gas heat capacity of a mixture is mass-consistent.** `UpdateHeatCapacity`
  mixed the per-kilogram pure-component heat capacities that `Fluid` caches
  from CoolProp by mole fraction, `C_p,g = Σ y_i C_p,i`, which is right only
  when every component has the same molar mass. It now uses
  `C_p,g = Σ y_i M_i C_p,i / Σ y_i M_i` [J/(kg·K)] (the molar mixture heat
  capacity over the mixture molar mass), so `ρ_g C_p,g = C Σ y_i M_i C_p,i`.
  For 50/50 H₂/CO₂ the gas heat capacity was about 174 J/(mol·K) and is now
  about 33; for CO₂/N₂ 15/85 it was 1.3 % high. Every user of the gas heat
  capacity takes this one value: the LTE and LTNE energy balances, the
  convective enthalpy terms, the thermal dispersion and heat-transfer
  correlations, and the logged heat capacity. The rule lives in
  `Fluid::MixtureHeatCapacity` / `Fluid::MixtureHeatCapacityFromSums`, with
  `Fluid::MolarHeatCapacity` (`C_p,i M_i`) and
  `Fluid::AdsorbedPhaseHeatCapacity` (the *Gas phase* adsorbed-phase option,
  which was already per component and is unchanged) replacing the two copies
  in the equations. A single-component gas, or a mixture with all but one
  mole fraction zero, is unchanged. The stream "enthalpy flows" and the pump
  and compressor energy do not use the gas heat capacity (they are molar flow
  × temperature and polytropic relations) and are unchanged. Tests:
  `GasHeatCapacityTest.*` and `GasHeatCapacitySolverTest.*`.
  **Results move for decks with light components**: Casas 2012 (H₂/CO₂)
  breakthrough times fall by 0.2–2.2 % (Fig. 4 50/50: H₂ 87.98 → 87.40 s,
  CO₂ 445.1 → 436.0 s; Fig. 6 75/25: H₂ 133.2 → 131.5 s, CO₂ 359.3 →
  351.4 s), the Streb 2021 He-carrier breakthroughs move (thermocouple RMSE
  down 2–18 %, the A-series outlet and front-time errors up), Nguyen 2022
  front times fall by up to 1.5 %, Vysyaraju 2024 He recovery rises from
  17.54 to 17.91 %; CO₂/N₂ and O₂/N₂ cycles move in the fourth digit.

### Fixed — bed-variable results at the right positions and boundary values
- **`Simulation::Series` / `sim.series` and `Simulation::Profile` /
  `sim.profile` read the log at the right places.** They mapped the raw cell
  rows, ghost cells included, onto the face positions, so every profile and
  every `z` query sat half a cell downstream, `bottom` returned the inlet ghost
  cell and `top` the last cell. In an adiabatic CO2/N2 breakthrough (feed
  313.15 K, 30 cells, t = 500 s) `series("temperature", location="bottom")`
  read 281.3 K, the ghost's 2 T_feed − T_1 (as low as 269 K during the run),
  and `top` pressure read 100004 Pa against the 100000 Pa set pressure. Now:
  - cell variables are at the cell centres and face variables (velocity,
    volumetric and molar flow rate, approach to fluidisation) at the faces;
    profiles leave out the ghost cells (`AxialProfile::faces` / `prof.faces`
    says which grid `z` is). `Simulation::Positions` / `sim.positions` are
    the cell centres, as documented; the new `FacePositions` /
    `sim.face_positions` are the faces.
  - `top` and `bottom` are the physical boundary values. Pressure,
    temperature and mole fraction are the logged values of the step's
    endpoint at that end, identical to `stream_series` (several streams are
    mixed by flow): at an inlet the value the boundary imposes (feed
    composition, feed temperature after any compressor, set pressure), at an
    outlet the bed-end face value. A closed end, and every other variable,
    gives the face value (the mean of the ghost and the end cell; a mole
    fraction's ghost, which is the face value itself; the nearest packed cell
    for adsorbent-only variables). A `z` between a bed end and the first cell
    centre interpolates between the two. The time base is unchanged.
  - The same logic is the new public header `adsim/solver/bed_series.h`
    (`ExtractBedTimeSeries`, `CellCentrePositions`, `IsFaceRow`);
    `ExtractBulkTimeSeries` / `ExtractComponentTimeSeries`
    (`parameterisation.h`) use it and take an optional `BedBoundarySources`
    (boundary log and steps) for the endpoint values. Front ends that fit to
    or plot these series move with it.
- **`stream_series(..., "pressure")`**: the bed-end face pressure of an
  endpoint (a pressure boundary's set value) is now readable
  (`BoundaryVariableType::PRESSURE`, appended to `kBoundaryVariableTypeList`),
  and `BoundaryLog::GetStreamData` returns an endpoint's logged data.
- Fitted or validated numbers that compared bed-end or `z` series (the
  example scripts' outlet curves and thermocouples) move slightly: outlet
  compositions and temperatures were already the end cell, so mainly
  thermocouple positions (half a cell) and inlet-end series change.

### Fixed — component mass balance with axial dispersion (Danckwerts inlets)
- **Inlet faces are Danckwerts.** With non-zero axial dispersion the component
  balance did not close at cyclic steady state: the species residual let a
  dispersive flux eb·D·C·(y₁ − y_in)/Δz cross the inlet face on top of the
  convective u·C·y_in, because the inlet ghost cell holds the feed composition.
  Totals still balanced (the dispersive fluxes sum to zero over components) but
  each component did not, by an amount proportional to D and nearly independent
  of the grid, and the stream logs, KPIs and the CSS mass-balance check never
  saw the extra species (reported: CO₂ out/in 1.0035 at D = 3.1e-4 m²/s and
  1.032 at 3e-3 m²/s). No dispersive flux now crosses the two boundary faces of
  the bed (`UpdateMoleFraction`), so an inlet face carries exactly u·C·y_in and
  the first cell may differ from the feed; at outlets and closed ends the ghost
  mirrors the first cell, so nothing changes there. The energy balance had the
  same defect: conduction across the inlet face, and an upwind value at the
  inlet face equal to the mirrored ghost temperature 2·T_in − T₁ rather than
  T_in. Both are fixed (`UpdateTemperature`, `UpdateFluidTemperature`). The WENO
  and van Leer reconstructions take the ghost value at the boundary faces
  (they averaged ghost and first cell, which at an inlet is not the feed).
  `GetTotalComponentMolesIn` (CSS mass-balance check, mass-balance panels)
  reconstructs the boundary flux without the dispersive term, matching the
  residual.
- **Flow boundaries on Darcy layers (Kozeny–Carman, Blake–Kozeny) deliver exactly the stream flow.**
  `FlowBoundaryCondition::Apply` computed the boundary-face flow of a Darcy
  layer with the ghost cell's mobility alone, K/(μ_ghost·Δz), while the face
  velocity that carries the gas into the bed obeys the face's Darcy law with
  the two half-cells in series. Whenever the ghost and first-cell viscosities
  differ (feed composition vs the bed gas, or the mirrored ghost temperature)
  the bed took in a different flow from the one the stream logs record (≈1 %
  in the new test; several per cent in the Streb 2021 decks, whose beds start
  in helium). The condition now matches the stream flow against the face's own
  velocity unknown, as Ergun faces (and the momentum-balance docs) already did.
- New regression tests `DispersionMassBalance.*`
  (`core/test/validation/test_dispersion_mass_balance.cpp`): a dispersive
  (D = 3e-3 m²/s) single-bed cycle with a flow-condition inlet at the bottom and
  a pressure-condition inlet at the top. Per step and component, the moles the
  stream logs record entering the bed equal the bed's inventory change (before
  the fix: errors of 0.4–50 %), and at cyclic steady state out/in of every
  component is 1 within 5e-4, the log quadrature error (before: O₂ 1.053, N₂ 0.983).
- Verification tests adjusted: `OpenSystemEnergyBalanceConvergesUnderGridRefinement`
  asserted first-order convergence of an energy closure error that came from the
  inlet ghost upwinding; the error is now 1.5e-3 → 6.6e-5 over the four grids
  (at least first order is asserted, and < 2e-4 on the finest grid instead of
  < 5e-3). `EndOfCycleProfileConverges` allows an absolute 2e-4 mol/m³ wobble at
  the ~1e-3 mol/m³ floor the gaps level off at.
- **Validation numbers that moved** (example decks re-saved; tolerances
  re-measured where they no longer held, each documented in the paper's
  README):
  - Ward 2022 (species fix; the bed had taken in ~2 % more CO₂ than fed): CO₂
    front time vs experiment +1.9/+1.7/+0.1 % → +4.2/+3.3/+1.5 % at 50/100/200
    cm³/min, outlet-curve NRMSE 0.18/0.13/0.029 → 0.33/0.24/0.12; thermocouple
    RMSE 1.19/2.06/4.01 → 1.18/1.18/2.76 K.
  - Streb 2021 (flow-face fix): CH₄ fronts ~3 % earlier; CH₄ front-time error
    A2 7.6 → 10.5 %, A3 0.3 → 3.1 %, A5 12.1 → 14.2 %, A6 8.7 → 11.7 %,
    A7 7.3 → 10.2 %, B1 0.7 → 2.1 %; A3 y_CH₄ NRMSE 0.21 → 0.50; A6 y_CO₂ NRMSE
    0.048 → 0.10; A2 thermocouple RMSE 14.1 → 18.6 K.
  - Nguyen 2023: CO₂ recovery of IX and X no longer exceeds 100 %
    (100.35/100.78 → 98.27/97.64 %); purity X 80.84 → 79.67 %, IX 84.69 →
    84.28 %; VIII productivity 0.1566 → 0.1545 (−1.4 %).
  - Within their bounds: Perez 2019 A–F recovery down by up to 1.1 pp
    (purity ≤ 0.2 pp); Arora 2021 O₂ recovery at cycle 50 4.24 → 4.53 % (PSA)
    and 19.88 → 20.74 % (PVSA); Ferré 2026 CH₄ recovery up by up to 0.75 pp;
    Casas 2012 H₂ breakthrough times ≤ 1.2 % earlier; Nguyen 2022 outlet-curve
    NRMSE improves (Fig. 4 50 %: 0.061 → 0.042); Vysyaraju 2024 He recovery
    17.32 → 17.54 %; Ward & Pini 2024 ≤ 0.2 %. `docs/validation.md` and
    `validation_summary.json` are regenerated.

### Changed
- **Cyclic steady state: the mass-balance methods no longer stop by accident.**
  On the way to CSS the per-cycle balance error swings through zero (a bed
  alternately gaining and losing a component), so the old single-cycle test
  with its 0.5 % band passed early: on the Yancy-Caballero (2020) 13X cycle at
  10 cells it stopped at the fifth cycle with purity and recovery 4–13 × 10⁻³
  from their converged values, and on `Perez2019_A`…`_F` 26 to 166 cycles
  early (KPIs up to 1.8 × 10⁻³ off). *Mass Balance (Total)* and *Mass Balance
  (Per-Component)* now also require the **end-of-cycle bed state** to have
  stopped changing: the largest relative RMS change since the previous cycle
  of the component loadings (scaled by the total loading), the mole fractions
  and the fluid, solid and wall temperatures
  (`CyclicSteadyState::EndOfCycleStateChange`) must be within the new
  **`css.state_tolerance`** (default 2 × 10⁻⁵; 0 switches it off). Every
  method's test must now hold on **`css.consecutive_cycles`** cycles in a row
  (default 3) before CSS is declared. The balance tolerance is unchanged: a
  tighter one cannot be relied on, because the balance carries a small
  systematic error of its own (up to 1 × 10⁻⁴ at CSS on `Perez2019_E`).
  With the defaults the Yancy cycle stops at cycles 59 and 100 (C1, C2,
  0-based) instead of 4, within 1.4 × 10⁻⁵ of the 150-cycle KPIs, and the
  example decks within 4.7 × 10⁻⁵ of their converged KPIs (`Perez2019_F`, which
  closes only ~1.5 % of its gap per cycle, 1.3 × 10⁻⁴); the stopping cycles
  and errors of 14 decks are tabulated in `docs/cyclic-steady-state.md`.
  `CyclicSteadyState` is class version 2: decks saved before it load with the
  defaults, so **existing decks that use a mass-balance method run longer and
  report converged (different) purity and recovery**, and Molar Profile decks
  stop at least two cycles later. `css.state_tolerance = 0` with
  `css.consecutive_cycles = 1` restores the old test. The solver's logs keep at
  least the cycles the check reads (one more than `consecutive_cycles`), even
  when `settings.max_cycles_to_keep` is smaller. `CSSResult`, `CSSRecord` and
  Python's `sim.css_history()` gain `state_change`, `state_tolerance`,
  `passing_cycles` and `required_cycles`; the solver log at level *medium*
  prints them after every cycle. Every paper case that runs to CSS moved
  (`docs/validation.md`): Perez 2019 A–F by up to 0.17 pp (A's recovery
  against the experiment crossed its bound, 0.294 → 0.417 pp; the bound is now
  0.5, see the Perez README), Nguyen 2023 I–X by up to 0.22 pp, Ferre 2026 by
  up to 0.32 pp, Ward & Pini 2024 by 0.04 pp and Vysyaraju 2024's He recovery
  by 0.64 pp (17.32 → 16.69 %); all within their bounds. Every bundled deck is
  re-saved with the two new fields, and the CSS decks with their new results.
- **Adsorbent library: one recommended entry per material** (catalogue schema 3).
  `library.json` groups its entries into 22 `materials` (id, name, family,
  `recommended` entry and why, `aliases`). A material's other entries are
  *variants*, each naming the product, form or data range that sets it apart;
  labels are the material name plus the variant ("Zeolite 13X (Zeochem Z10-02,
  CO2/N2 to 1 bar)") instead of "Material (Author year)". Each isotherm records
  its `data_quality` (measured with tabulated data, measured fit only, refit of
  other authors' data, copied from an earlier paper, not yet checked). The
  validator, the native library tests and the generated docs page (now grouped
  by material) enforce it. `family` moved from the entry to its material.
- `AddLibraryAdsorbent`, `Simulation.add_library_adsorbent` and
  `skarstrom.adsorbent_path` accept a material id (e.g. `"zeolite-13x"`) and add
  or return that material's recommended entry; `AdsorbentLibrary` gains
  `materials`, `FindMaterial`, `Resolve` and `Variants`. Entry keys are unchanged.

### Fixed
- **WENO and van Leer advection** (`advection_scheme = weno | van_leer`) are fixed; Upwind (the
  default, and every bundled deck) is unchanged. Compared line by line with the WENO of
  Haghpanah et al. (2013) / Yancy-Caballero et al. (2020) (PEESEgroup/PSA `CycleSteps/WENO.m`):
  (1) WENO's upwind-side candidate was the face value one face upstream, `(y[n-1] + y[n])/2`,
  instead of the extrapolation `3/2 y[n] - 1/2 y[n-1]`, which made the scheme first order;
  (2) both schemes gave the bed-end faces the average of the ghost and the first cell, but the
  ghost holds the inlet composition itself, so the bed received `(y_feed + y_1)/2` while the
  stream logs, KPIs and cycle mass balance counted `y_feed` (CO2 balance 0.94-0.98 on the
  Yancy-Caballero 13X cycle at 10 cells with dispersion off, where upwind closes to 1.0000);
  (3) built component by component, the limited faces did not sum to the cell sum, and
  linearised about a steep front they weight the downwind cell, so a departure of the cell sum
  of mole fractions from one grew without bound (van Leer reached a sum of 90 in a CO2/N2
  breakthrough). Now the bed-end faces take the upwind value; the first interior face next to
  an inlet extrapolates from the inlet value half a cell away (`2 y_1 - y_0`, as the reference;
  van Leer doubles that slope likewise); void cells are treated as well mixed; WENO uses the
  cell widths on non-uniform grids; and the interior faces are rescaled so their mole fractions
  sum to the upwind cell's sum, which couples the mole fractions two cells apart, so the
  Jacobian half-bandwidth for these schemes grows from `2 neq` to `2 neq + components - 1`.
  On the Yancy-Caballero 13X cycle at 10 cells, WENO now gives a CO2 recovery of 90.0 % for
  the minimum-energy case (PSASimulator.jl WENO reference 89.7 %; was 77.5 %), and its CO2
  balance with dispersion off closes to 1e-5 at every grid. Tests:
  `AdvectionSchemeConservation.EachSpeciesIsConservedThroughTheBedEnds` and
  `AdvectionSchemeConservation.CellMoleFractionsKeepSummingToOneAtASteepFront`.
- `CreateBoundaryValue` and `CreateInletTemperatureProfile` know the double-exponential
  profile. The web document layer builds profiles through the first factory, so a deck
  with a double-exponential pressure boundary (the Perez 2019 decks) lost it to a null
  profile, i.e. a constant 1 bar, in the browser. A unit test now checks every listed
  kind has a factory entry.

### Added
- A **Blake–Kozeny** momentum model (`MomentumModel::BlakeKozeny` = 2, builder kind and
  registry choice `blake_kozeny`, Python `layer.set_momentum("blake_kozeny")` /
  `add_layer(..., momentum="blake_kozeny")`, display name "Darcy (Blake-Kozeny)"):
  Darcy's law with K = ε³ d_p² / (150 (1 − ε)²), i.e. Ergun's viscous term without the
  inertial one. It is the pressure-drop law of Haghpanah et al. (2013) and Estupiñan
  Perez et al. (2019, SI Eq. S5). It takes the Darcy face path everywhere, like
  Kozeny–Carman. The stored values 0 (Kozeny–Carman) and 1 (Ergun) are unchanged, so
  existing decks load as before; a deck that uses the new value needs this engine.
- An optional **constant molecular diffusivity** for the macropore LDF coefficients
  (`macropore`, `macropore_capacity`, `macropore_knudsen`): registry leaf
  `adsorbent[...].kinetics[<component>].molecular_diffusion_coefficient` [m2/s], Python
  `set_rate_model(..., molecular_diffusion_coefficient=1.6e-5)` / `kinetics[c].molecular_diffusion_coefficient`.
  Above 0 it replaces the binary-diffusion D_m computed from T and P (which scales as
  1/P, so at 3-10 kPa evacuation pressures the coefficient was 10-33 times the value
  of the papers that use a constant D_m, e.g. Haghpanah et al. 2013 and Estupiñan Perez
  et al. 2019, Table 1). The default 0 keeps the previous behaviour bit for bit. Class
  versions bumped (macropore and macropore + Knudsen to 3, capacity factor to 2); older
  decks load with 0. Haynes is unchanged: its film term uses D_m with the Sherwood
  correlation, so a constant there is not a drop-in.
- **Ideal Adsorbed Solution Theory (IAST)** as a per-adsorbent **mixing rule**
  (`AdsorptionSource::mixing_rule_`, `MixingRule::None|IAST`; registry
  `adsorbent[<id>].mixing_rule` with choices `none`, `iast`; builder
  `SetMixingRule` and kind list `MixingRuleKinds()`; Python `ads.mixing_rule`,
  `ads.set_mixing_rule(kind)`). Under `iast` the mixture loadings come from the
  components' pure isotherms (any kinds, mixed freely), solved per cell by a
  bracketed Newton iteration on the reduced spreading pressure, warm-started
  from the cell's last solve. Closed-form spreading pressures for most kinds,
  Gauss–Legendre quadrature for Toth, Toth 2 and Redlich–Peterson. Absent
  components (y ≤ 0) take the infinite-dilution limit, so the loading is
  smooth through zero mole fraction. The evaluation functions are public in
  `adsim/equations/iast.h` and in Python as `skarstrom.evaluate_iast`. The
  solver log warns when IAST is combined with an extended (competitive) kind.
  The capacity-factor macropore coefficient (`macropore_capacity`), which
  evaluates q\* itself, takes the IAST mixture loading under the IAST rule.
  `AdsorptionSource` is at cereal version 1; older decks load with `none` and
  run exactly as before. Docs: the isotherms page.
- `IsCompetitiveIsothermKind(type)` and `PairedIsothermKind(type, out)` in
  `adsim/equations/isotherm_kinds.h`: which kinds couple components through a
  shared denominator, and each one's pure counterpart with the same parameters
  (and back). The IAST warning uses them and names the pure kind to switch to.
- A **double exponential** profile (`DoubleExponentialBoundaryValue`,
  `BoundaryValueType::DOUBLE_EXPONENTIAL`, builder kind `double_exponential`, Python
  `sk.double_exponential(end, rate_1, rate_2, weight, start=)`): a two-term decay
  towards the end value, v(t) = v1 + (v0 − v1)[w e^(−r1 t) + (1 − w) e^(−r2 t)], the
  form Estupiñan Perez et al. (2019) fitted to their rig's evacuation pressure. The
  rates are entered in 1/s. Registry leaves `<head>_double_exponential_start/_end/
  _rate_1/_rate_2/_weight` and `pressure_double_exponential_end/_rate_1/_rate_2/_weight`.
- A **capacity-factor macropore LDF coefficient** (`MassTransferCoefficientMacroporeCapacity`,
  kind `macropore_capacity`): k = (15 ε_p D_m / (τ r_p²)) · c / (ρ_p q\*), the
  Haghpanah et al. (2013) form the Rajendran-group papers use. The ratio is evaluated
  from the isotherm at max(c, c_floor), with the isotherm re-evaluated at the floored
  concentration, so it stays finite where c and q\* both vanish; capped at 1. The plain
  `macropore` kind is unchanged (it never applied the factor; its `ratio = 1.0` stub is
  gone).

### Changed
- **Native WebAssembly exceptions by default.** A new CMake option, `ADSIM_EMSCRIPTEN_EXCEPTIONS`, chooses the C++
  exception scheme of an Emscripten build: `wasm` (the new default, `-fwasm-exceptions`) or `js` (`-fexceptions`,
  the previous behaviour). The whole engine tree is compiled with the chosen flag, and `Adsim::adsim` now passes it to
  the consumer's compile and link. On the Skarstrom web solver, native exceptions ran 10 % faster on six example decks
  (4–17 %) with a 12 % smaller `.wasm`; they need Chrome 95, Firefox 100, Safari 15.2 or Node 17. **Consumers:** drop
  your own exception flags. `-sNO_DISABLE_EXCEPTION_CATCHING=1` is rejected by Emscripten together with
  `-fwasm-exceptions`; keep the old scheme with `-DADSIM_EMSCRIPTEN_EXCEPTIONS=js`. The subproject consumer test now
  checks that an exception thrown inside the engine reaches the consumer's catch, and CI runs it with each scheme.
- CI builds the Emscripten consumer with Emscripten 6.0.13 (was 4.0.10).
- **Perez 2019 example decks** (`Perez2019_A/B/C1/D1/E/F`) now use the paper's model
  forms, after an audit of the deck against the paper and its SI
  (`core/test/papers/perez_2019/README.md`, provenance and attribution tables): the
  two-term exponential evacuation pressure of Sect. 4.3 (0.019 and 0.101 1/s; the split
  of the terms is not printed and an equal split is used, chosen on A), wall ends at
  ambient (SI Table S2) and the capacity-factor LDF of SI Eq. S4. The four-step
  recoveries, 8-10 pp above the experiments before, are now within 0.7-1.1 pp for A and
  B; the attribution shows the evacuation pressure history was the cause. The decks keep
  Kozeny-Carman: the paper calls the pressure drop negligible, Ergun is worth 0.3 pp, and
  it fails the initial-condition calculation of the fast-feed connected configuration in
  `ConfigD1_ConnectedBC` (an open engine item). `examples/perez_2019.py` builds the same
  decks.

### Changed
- **Ward 2022 example decks** (`Ward2022_Q50ccm/Q100ccm/Q200ccm`) now use the
  paper's model forms, after an audit of the deck against the paper text
  (`core/test/papers/ward_2022/README.md`, provenance and attribution tables):
  helium as the carrier (the CO2/He experiment), the paper's dispersion
  coefficient (Eq. 4 with D_m = 0.605e-5 m2/s, times the factor 4 the paper
  fitted), -dU_1 + RT as the heat of adsorption (32.4 kJ/mol), the wall ends
  held at ambient (Table 4) and the Ergun pressure drop. The CO2 front now
  arrives 0.1-1.9 % late instead of 4-6 % early, the curve NRMSE is 3-18 %
  instead of 34-38 %, and the temperature rise is within 0.9-2.9 K of the
  measurement (it was 8-11 K high). The Ward test's temperature-rise
  indicator follows the paper's Table 9 (against the feed temperature).
  `examples/ward_2022.py` builds the same decks.

### Added
- A **log-linear** profile (`LogLinearBoundaryValue`, `BoundaryValueType::LOG_LINEAR`,
  builder kind `log_linear`, Python `sk.log_linear(end, start=, time=)`) for pressures,
  flows and inlet temperatures: v(t) = v0 (v1 / v0)^(t / t_r), a geometric ramp whose
  relative rate of rise d(ln v)/dt is constant. A pressurisation from vacuum with a
  linear ramp puts most of the relative pressure rise, and with it the inlet pressure
  gradient and approach to fluidisation, in the first seconds; the log-linear ramp
  spreads it over the step. It holds its end value after the ramp, and ramps linearly
  when the start or end value is zero or negative. Registry leaves
  `<head>_log_linear_start/_end/_time` and `pressure_log_linear_end`.

### Changed
- Solver residual, byte-identical: the competitive concentration-based
  Langmuir isotherms (`ExtendedLangmuirIsotherm_Concentration`,
  `ExtendedDualSiteLangmuirIsotherm_Concentration`, the kinds most example
  decks use) gain the single-pass `UpdateAll` the pressure-based ones had, so
  their shared denominators and exponentials are computed once per cell
  instead of once per component (Illustrative TSA: 17% fewer instructions for
  the run, the isotherm 2.6x cheaper). Every competitive isotherm's
  per-component `Update` is now that pass for one component, which also
  serves the capacity-factor MTC that re-runs the isotherm. The binary
  diffusion model memoises the Neufield collision integral and `T^1.5` per
  cell on the exact temperature bits (`TemperatureCache`), since most residual
  evaluations repeat the previous temperature in most cells (15% fewer
  instructions in that model).
- Solver residual, byte-identical: the Langmuir and dual-site Langmuir isotherm
  functions (and their concentration-based variants) compute each `exp(q/(R T))`
  factor once instead of twice (GCC's default `-fmath-errno` does not merge the
  duplicate calls, so Linux builds paid 4 `exp` per cell per component), and a
  constant mass-transfer coefficient is evaluated once per layer instead of once
  per cell (`ExpressionParam::IsRunConstant()`). -17% instructions per residual
  call on FiveStepProductPurge; CSV output identical on six example decks.
- Solver: two settings are now on by default (`IDASettings`, transient):
  - `suppress_alg_` with `pressure_differential_`: the pressure is tagged a
    differential unknown (its residual carries dP/dt) and IDA's local error
    test leaves out the algebraic entries, i.e. the ghost cells (mirror
    relations with no time derivative) and the face velocity. Their jumps at
    a boundary-condition switch used to fail the error test and force tiny
    restart steps; the cyclic example decks now take 40-60% fewer steps
    (1.3-2.1x wall) at the same accuracy against a reltol = 1e-8 reference.
  - `alias_lte_solid_temperature_` and `derived_concentration_`: with local
    thermal equilibrium the solid temperature is a copy of the fluid
    temperature instead of an unknown, and the gas concentration is derived
    from P and T (C = P / RT; no residual reads dC/dt). Two unknowns per cell
    fewer, a narrower band, one colour fewer (1.1-1.25x).
  The equations solved are unchanged; the integration path is not, so
  results move in their last digits (the example decks are re-saved). The
  per-step boundary/outlet logs have fewer points than before, so a value read
  from them at an exact time is interpolated over a longer step. Set the
  switches to `false` for the previous behaviour.

### Changed
- Solver: the banded Jacobian matrix is sized from the coupling reach a
  difference-quotient probe measures at the initial state (12 / 16 instead of
  2 neq = 20 / 20 diagonals at ten equations per cell), about 45% less LU
  work; results are byte-identical. If a Jacobian audit later finds a
  coupling outside the band, the solver widens it to the full reach and
  restarts from the current point (logged at medium level).
  `IDASettings::band_from_probe_` switches it off.

### Added
- Transient (not serialised) `IDASettings` switches for solver experiments,
  all off by default: `alias_lte_solid_temperature_` and
  `derived_concentration_` (drop Ts and C from the DAE unknowns),
  `suppress_alg_` (`IDASetSuppressAlg`), `delta_cj_lsetup_`
  (`IDASetDeltaCjLSetup`) and `pressure_differential_` (tag P as a
  differential unknown). `bench_solver` gets `--unknowns`, `--band`,
  `--suppress-alg`, `--dcj`, `--p-diff` and `--reltol`; its report adds the
  band widths, out-of-band count and whether the band was widened.
- **Engine internals de-duplicated** (no behaviour or deck-format change): the
  twenty concrete isotherm classes share one `ParamIsotherm<Params, kType>` base and
  only declare their `Update` (archives are byte-identical: the parameter map keeps
  its name and each class its `cereal_class_version`); the four flow boundary value
  kinds keep their `profile_` and serialisation in `IFlowBoundaryValue`; the
  per-stream `GetComponentMolarFlow` / `GetEnthalpyFlow`, which were copied into
  every stream kind, live once in `IBoundaryStream` with a `BedInletTemperature`
  hook for the pump-driven kinds' compression heating; `EquationLog` / `BoundaryLog`
  getters and face interpolation, the KPI bed-accumulation integrals and the
  compressor / vacuum-pump cost structs each have a single implementation. The
  unused `CalculateAccumulation_Volume` is gone. Test fixtures that were copied
  across files (the validation PSA cycle, the feed/hold integration deck, the
  Ogata-Banks solution) are shared headers now.

## [0.3.0] — 2026-10-25

### Added
- Wall pressure design (`adsim/objects/wall_design.h`): the minimum wall
  thickness of a cylindrical column for ASME BPVC VIII-1 UG-27 (with
  Appendix 1-2 for thick shells) or ASME B31.3 304.1.2, with a corrosion
  allowance, and for full vacuum by the formula method of ASME Code Case
  2286 (VIII-2 4.4.5.1). New `BedWall` fields `design_code_` and
  `corrosion_allowance_`; registry paths `wall.design_code`,
  `wall.design_pressure`, `wall.design_temperature`,
  `wall.weld_joint_efficiency`, `wall.corrosion_allowance`,
  `wall.vacuum_service`. Elastic moduli of the wall materials
  (`GetElasticModulus`). Python: `sim.wall_design()`.
- Operating-limit warnings (`adsim/solver/operating_limits.h`): at the end of
  a completed run the solver checks the last logged cycle and warns, per step,
  when the pressure passes the wall's design pressure or goes below
  atmospheric without vacuum service, when the wall passes its design
  temperature, and when the approach to fluidisation passes 1 (with where and
  when). At the start of a run it warns when the wall is thinner than its
  design requires (`WallDesignWarnings`). Python: `sim.operating_limits()`;
  CLI: `warnings` in `skarstrom run --json`. Docs: `docs/wall-design.md`.
- Warning switches (`adsim/solver/warning_settings.h`, `IDASettings::warnings_`,
  cereal version 6): each advisory warning (pressure step change, wall
  thickness, design pressure, vacuum, design temperature, fluidisation) can be
  switched off in the deck; registry paths `settings.warn_<name>`. Docs: the
  Warnings section of `docs/solver-settings.md`.

### Changed
- Inlet mole fractions (pressure boundaries and every stream kind) and the
  initial-condition mole fractions are expressions (`ExpressionMap`, a
  component-keyed map of `ExpressionParam`s in `adsim/math/expression_map.h`)
  instead of plain numbers, so a composition can use global variables
  (`y_feed`, `1 - y_feed`) or a fraction unit (`21 [%]`, `400 [ppm]`). The
  registry lists `endpoint[ID].y[C]` and `initial.y[C]` as Expression
  parameters (a plain number stays tunable). The sum-to-1 checks use the
  evaluated values. Archive versions: `IBoundaryStream` 5,
  `PressureBoundaryCondition` 5, `InitialConditions` 3; older decks, whose
  mole fractions are numbers, still load. A deck saved by this version does
  not load in 0.2.0 or earlier.
- Vacuum service no longer sizes the wall with 1 atm as an internal pressure
  in UG-27 (not conservative: buckling governs) and no longer ignores the
  design pressure; the larger of the internal and the external requirement
  is used.
- `BedWall` cereal version 4. Its version was declared in `wall.cpp`, where
  no archive saw it, so decks carry version 0 for the wall; they load
  unchanged. The declaration is now in `wall.h`.

### Removed
- `CalculateAsmeMinWallThickness`, `AsmeWallThicknessResult` and
  `BedWall::CalculateMinThickness`; use `DesignWall` / `DesignBedWall`.

## [0.2.0] — 2026-10-24

### Added
- Sensitivity analysis and KPI error bars (`adsim/solver/sensitivity.h`):
  Morris screening (which inputs move the KPIs) and Monte Carlo propagation
  with Latin hypercube sampling (median, 5-95 % band with bootstrap
  intervals, a convergence check, and each input's share of the spread from
  a rank regression). Inputs are registry paths with uniform, triangular or
  normal distributions; plans are seeded and reproducible on every platform.
  The engine plans and analyses; callers run the samples
  (`EvaluateSample`), which must reach cyclic steady state by default.
  Python: `skarstrom.sensitivity` (`Study.screen`, `Study.propagate`,
  parallel `workers`); CLI: `skarstrom sensitivity`. Docs:
  `docs/sensitivity.md`.

### Changed
- `examples/notebooks/sensitivity_analysis.ipynb` now runs the two-stage
  study instead of one-at-a-time finite differences.
- The solver log no longer ends a run with the IDA statistics line
  ("Statistics: N steps, N residual evals, ..."); the completion line with
  the wall time is the summary. The counters stay available through
  `CoupledSolver::PrintStatistics()` and `SolverData::GetStatistics()`.

### Removed
- The unsourced legacy adsorbent library entry `zeolite_13x` ("Zeolite 13X,
  CO2/CH4 (legacy)") and its blob `data/adsorbents/zeolite_13x.json`. The
  sourced 13X entries (`zeolite_13x_ferre2026`, `zeolite_13x_perez2019`,
  `zeolite_13x_streb2021`) remain.
- The `bulk_density` property from `data/adsorbents/library.json`: the bed
  (bulk) density depends on how a layer is packed, so it is a layer property,
  not a property of the adsorbent material. The layer's bed voidage and the
  particle density determine it.

### Fixed
- The polynomial pump efficiency model had no P^2 term: it evaluated
  eta = a + b P + c P^3 + d P^4 + e P^5. It is now the full fifth-degree
  polynomial eta = a + b P + c P^2 + d P^3 + e P^4 + f P^5, with a new
  coefficient `f` (registry and Python `poly_f`). `PolyPumpEfficiencyModel`
  is now class version 3; a version-2 deck's c, d, e load as d, e, f with
  c = 0, so its efficiency and results are unchanged.
- `AddComponentToDeck` now seeds a 0 adsorbed-phase heat capacity for the new
  component (and `RemoveComponentFromDeck` drops it), as the Python
  `AddComponent` already did. A front end that lists the map's entries showed
  no heat-capacity field for a component added through the shared cascade
  until the deck had been run once.

### Security
- Release wheels build with cibuildwheel v4.3.4 (was v4.3.1), whose pip,
  certifi and manylinux images no longer trust the compromised
  eMudhra/emSign certificate authority (October 2026). Every CI job now
  upgrades pip before installing from PyPI. The 0.1.0 files on PyPI were
  checked (hashes, RECORD files, no `.pth` files, Python sources identical
  to the tag) and are unaffected.

## [0.1.0] — 2026-10-20

First public release of the Skarstrom engine under the Apache License 2.0.

### Changed (units)
- Unit strings are parsed and converted with the LLNL units library
  (BSD-3-Clause, new submodule `external/units`, compiled into the engine)
  instead of a hand-kept table. Any unit LLNL units knows works, with SI
  prefixes, and typos are still errors (strict parsing). The engineering
  conventions are kept: everything after `/` is in the denominator
  (`W/m2-K` = W/(m²·K)), a lone `degC` is absolute, `C`, `F` and `A` are
  °C, °F and the ångström, and `ccm`, `SCCM`, `SCFM` and `Nm3` are gas-flow
  units at 0 °C and 1 atm. Only units of length, mass, time, temperature and
  amount of substance are accepted. `LookupUnitSymbol` is removed.
- An inline `[unit]` now gets the unit category of its dimension: heat
  transfer coefficients and viscosities were Custom before.
- Python: pint is a dependency, and `skarstrom.ureg` / `skarstrom.Quantity`
  give its application registry. Scripts set values with pint quantities
  (`sim.wall.thickness = 3.1 * u.mm`) or a number and its unit
  (`node.set("thickness", 3.1, "mm")`); every example script and the docs
  now do, instead of unit strings. Strings still accept deck-syntax
  quantities (`"2 [min]"`) and expressions.

### Added (Python)
- Pump-driven streams, whose flow follows from a shaft power (constant, or a
  polynomial of the boundary pressure): `sk.PumpStream(name, power=, pump=)`,
  `adsim::build::AddPumpStream`, and registry entries `endpoint[N].power`,
  `power_poly_a..f` and the pump's settings.
- Buffer tanks: `inlet.set_buffer_tank(source, percentage=)`,
  `adsim::build::SetBufferTank`, `endpoint[N].buffer_tank.percentage` [%].
- Experimental data stored in the deck: `sim.experimental_data` (`add`,
  `load_csv`, `remove`, `clear`, `[name].compare(x, y)` for RMSE / NRMSE).
- `to_script` and `describe()` cover pump-driven streams and buffer tanks.

### Added (KPIs)
- A **recycle** role in the KPI configuration (`Recycle_Endpoints` in the
  deck, `recycle=` in `sim.kpi.track` and `adsim::build::SetKPI`): inlets that
  return product to the bed (a purge or a light-product pressurisation fed
  from the product). Recovery, productivity, production rate and specific
  energy use the net product (product minus recycle); purity stays the purity
  of the collected product. Before, a recycle inlet could only be ticked as
  product, which subtracted it from the purity as well.

### Changed
- Python: `skarstrom.Pump` is now the pump specification the builder takes
  (`sk.Pump("vacuum_pump", efficiency=0.7)`); the deck object behind
  `endpoint.pump` is `skarstrom.model.Pump`.
- Contact is through LinkedIn instead of email (README "Get involved",
  SECURITY.md, CODE_OF_CONDUCT.md); the README and docs gain the logo, a
  "Why Skarstrom" introduction and a call for steering-committee members.
- README: a "Browser interface" section introducing the Skarstrom website
  (https://www.skarstrom.co.uk/), the browser front end built on this engine.

### Fixed (registry)
- `layer[ID].heat_dispersion_coefficient` (constant model) was declared in
  m²/s; it is the energy balance's effective axial thermal conductivity, in
  W/(m·K). Its unit, name and suggested bounds are corrected (values are
  unchanged).
- The extended-Sips `B1` parameter was declared in 1/bar; the isotherm
  multiplies it by the pressure in Pa, so it is 1/Pa (unit and bounds fixed,
  values unchanged). A value written with a unit, e.g. Table 3's
  `"1.129e-4 1/bar"`, now converts correctly.

### Added (choices)
- Model and mode switches are named: the registry lists `choices` for
  enum-valued integers (`css.method`, `settings.acceleration_method`,
  `settings.advection_scheme`, `energy_balance.energy_balance_mode`,
  `cycle.step_time_mode`, `bed.thermal_equilibrium_mode`, the fluid modes,
  `layer[ID].momentum_model`); `SetText` (Python strings, `--set`) accepts the
  name, e.g. `sim.css.method = "mass_balance_component"`, and the Python
  object layer reads names back. Numbers still work; out-of-range numbers are
  rejected.

### Removed
- `adsim/tools/eigenlldb.py`, Eigen's LLDB pretty-printer (unused).
- Example decks renamed so every name is a Python identifier:
  `Ward_2022_Q=50[ccm]` → `Ward2022_Q50ccm` (and `Q100ccm`, `Q200ccm`),
  `Nguyen_2022_Fig4`/`Fig5` → `Nguyen2022_Fig4`/`Fig5`, `streb_2021_A1`/`A3`
  → `Streb2021_A1`/`A3`. The shell-safe-name allowlist is now empty.
- Python: `Simulation.steps` is a collection of step objects (names:
  `sim.steps.names`); `skarstrom.examples()`, `example_path()` and
  `copy_example()` are replaced by the `skarstrom.examples` module; the CLI's
  `examples` command lists names without `.json`.

### Fixed
- Ergun layers no longer fail at a step switch that opens a flow boundary
  whose streams follow the face pressure (a vent valve, a pump) on a bed the
  previous step left away from the streams' pressure: the Perez 2019 D1 cycle
  with a 0.2 m/s feed, the light product vented through a Cv = 10 valve at
  1 bar and recycled through a connection to a 150 s pressurisation ending at
  101325 Pa failed IDA's corrector at the start of every feed step and the
  consistent-IC restart in the third cycle ("convergence test failure (may
  indicate initial condition calculation failure)"). The face velocity has to
  jump from a standstill to several m/s, and the implicit Ergun residual's
  slope in it changes by an order of magnitude across the jump. The flow
  boundary of an Ergun layer now solves its face problem (the streams' flow
  law and the Ergun law) for the face velocity and ghost pressure inside the
  residual and imposes the root, so the two unknowns' Jacobian block is the
  identity (`FlowBoundaryCondition::ImposeFlowFace`). Results of Ergun decks
  change only within the solver tolerance; Darcy (Kozeny-Carman) faces are
  untouched. Regression test `ErgunValveSurgeTest`;
  `Perez2019Test.ConfigD1_ConnectedBC` is back on Ergun.
- Ergun layers no longer fail at a step switch that moves the boundary
  pressure by a large factor (a constant-pressure feed re-pressurising an
  evacuated bed, 0.05 -> 1 bar, failed IDA's corrector at the start of every
  pressurisation step; the restart from a consistent initial condition that
  rescued it gave out after about ten cycles, and sooner with a bottom void).
  The boundary faces' Ergun coefficients used the ghost cell's density and
  viscosity, a mirror construct (P = 2 Pf - P1) that jumps with the boundary
  pressure; they now use the interior neighbour's. Results of Ergun decks move
  by less than 0.1 % (the density change over half a cell); Darcy
  (Kozeny-Carman) results are unchanged. Regression test
  `ErgunPressureSwitchTest`.
- The KPI configurations of the light-product-pressurisation decks
  (Perez2019_E, _F, Nguyen2023_Exp_VI–X) listed the pressurisation inlet as
  product. That stream is N2-rich gas moving within the cycle, not CO2
  product, and its negative contribution made the CO2 purity NaN. The product
  is now the evacuation stream alone, as in Perez2019_D1 and the papers
  (a defined CO2 purity instead of NaN). The validation numbers were already
  computed this way and do not change.
- Perez2019_E, _F and Nguyen2023_Exp_VIII–X stopped at their 250-cycle cap
  before reaching cyclic steady state. Their CO2 inventory settles slowly:
  the per-component mass balance closes to 0.5 % only after about 450 (E),
  1000 (F) and 300–500 (VIII–X) cycles, long after the gas profiles stop
  changing. The caps are now 1500 (Perez E/F) and 1000 (Nguyen LPP), so the
  decks and the validation report results at CSS. Perez E/F agree better
  with the experiments (CO2 purity 1.9 and 2.3 pp low instead of 5.0 and
  8.9 pp); their tolerances are tightened accordingly.
- A value with an inline unit of the wrong dimension (`step time = "2 [bar]"`)
  was converted silently (to 200 000 s); it is now rejected.
- Competitive isotherms summed their per-component terms in `std::unordered_map`
  order, which depends on the order of the entries in the deck file and on the
  standard library. A deck and its re-saved copy could therefore run to
  slightly different numbers (seen on macOS arm64, where FMA contraction makes
  even a two-component sum order-dependent). Isotherm parameters are now kept
  in key order (`std::map`); the deck format is unchanged.
- CI: the "example decks unchanged" check no longer uses `A && B || C`
  (shellcheck SC2015) and also fails when the test suite deletes a file.

### Added
- The solver log warns, once at the start of a run, about a step whose
  constant pressure boundary differs by more than 25 % from the pressure the
  step is expected to start at (the end value of the previous pressure
  boundary, or the initial pressure for the first step): the bed end then
  sees an instantaneous pressure change the integrator has to resolve, and a
  linear or exponential ramp avoids it. Ramps never trigger it
  (`BoundaryPressureStepChangeWarnings`, `adsim/solver/deck_warnings.h`).
- Building a simulation from scratch. In the engine, `adsim/objects/deck_builder.h`
  creates components, adsorbents (isotherms by kind, rate models), the bed
  shape, wall, layers and their models, steps, pressure and flow boundaries,
  streams and valve streams, profiles, pumps, connections, energy sources,
  KPI endpoints and user variables, each model chosen by a snake_case kind
  name. In Python, `Simulation.new()`, `add_variable`, `add_component`,
  `add_adsorbent`, `add_layer`, `add_step`, `connect`, `set_isotherm`,
  `set_rate_model`, `set_shape`,
  `set_material`, `step.bottom.set_streams` / `set_pressure` / `set_closed`,
  `add_energy_source`, `sim.kpi.track`, with `skarstrom.velocity`,
  `molar_flow`, `linear`, `exponential`, `Stream`, `ValveStream`, `Pump`;
  `skarstrom.build_kinds()`, `sim.describe()` and `skarstrom.to_script(sim)`,
  which writes the builder code for any deck; a test checks that this round
  trip reproduces all 33 example decks and that five run to identical results.
- Every example script (`examples/*.py`) builds its cases from scratch in
  code, every value explicit with its unit, instead of loading a deck;
  `build(case)` returns the simulation and a test checks it reproduces the
  shipped deck.
- Parameter registry: every value of the example decks is now settable —
  initial conditions (`initial.*`), `fluid.*` (adsorbed-phase heat
  capacities), `layer[ID].cells` and `mass_fraction[...]`, profile ramps of
  inlet temperatures and flows, valve coefficients, the other pump
  efficiency models, energy sources (`step[ID].energy[layer].power`), void
  regions, particle heat transfer, `cycle.step_time_mode`, and inlet and
  backflow compositions of every endpoint.
- Python: deck objects held in variables (a step, a layer, `sim.bed`) follow
  structural edits instead of going stale.
- Python: the deck as typed objects (`skarstrom.model`). `sim.steps["Feed"].time`,
  `sim.steps[0].bottom.temperature`, `sim.steps["Feed"].bottom.y["CO2"]`,
  `sim.layers[0].length`, `sim.adsorbents["13X"].isotherm["CO2"].b1`,
  `sim.variables[...]`, `sim.cycle`, `sim.bed`, `sim.wall`, `sim.settings`,
  `sim.css`, `sim.energy_balance`, `sim.info`. They are views onto the
  parameter registry (string paths keep working); misspelt or absent fields
  raise `AttributeError` with suggestions; editors complete the fields and
  show their SI units. A test checks every registry parameter of every
  example deck is reachable this way.
- Python: `skarstrom.examples`, one attribute per bundled deck
  (`skarstrom.examples.Perez2019_A.load()`, `.path`, `.title`, `.copy_to()`),
  generated from `examples/decks/` by `scripts/gen_examples_module.py`
  (checked in CI).
- Units. A bare number is SI, as before; a value can now carry its unit as a
  string (`"2 min"`, `"25 degC"`, `"35 kJ/mol"`, also `"2 [min]"`), as a pint
  Quantity, or with `unit=`; `get(..., unit=...)` reads in any unit;
  `Simulation(deck, require_units=True)` rejects bare numbers for dimensional
  values; `skarstrom.units.convert()`. In the engine: `ParseUnitExpression`
  reads compound units (`W/(m^2*K)`, `W/m2-K`, `kJ/(kg-K)`, `1/bar`,
  `Nm3/h`...), `ConvertUnit` and `SplitQuantity`; `Parameter::SetText` (and
  so the CLI's `--set`) takes `"<number> <unit>"` and checks the unit's
  dimension against the parameter's. `ParseUnitSuffix` falls back to unit
  expressions for symbols missing from its dictionary.
- `DeckRoundTrip` tests: every example deck's numbers parse exactly and
  survive save/reload bit for bit, and a run is unchanged when every map in
  the deck is listed in reverse order.

### Added
- The `adsim` C++20 simulation engine: 1D finite-volume packed-bed model
  (mass, momentum, energy and wall energy balances; 20 isotherm models;
  linear driving force kinetics with several coefficient correlations;
  dispersion and heat-transfer correlations; CoolProp fluid properties),
  multi-step and multi-bed cycles with valves, pumps, buffer tanks and step
  connections, SUNDIALS IDA integration, cyclic-steady-state detection and
  Aitken acceleration, KPIs, a parameter registry and a dlib optimiser.
- Public headers for front ends: `adsim/objects/deck_io.h` (load decks with
  either cereal root), `adsim/objects/deck_edit.h` (component, adsorbent and
  connection edits that keep a deck consistent),
  `adsim/equations/isotherm_kinds.h` (isotherm kind tags, parameter access,
  evaluation) and `adsim/solver/isotherm_fit.h` (the isotherm fitter).
- CMake build with presets for Linux, macOS and Windows; `find_package(Adsim)`
  install/export; use as a subproject (`add_subdirectory`) with front-end
  friendly defaults; builds with Emscripten.
- GoogleTest suite (~670 tests): unit, integration, analytic verification and
  quantitative validation against ten published studies.
- The `skarstrom` Python package (scikit-build-core + pybind11) and CLI
  (`run`, `params`, `examples list|run|path|copy`, JSON output), bundling all
  example decks and the adsorbent library; structural edits from Python
  (`add_component`, `remove_component`, `add_library_adsorbent`).
- 33 example decks, a runnable Python script per deck family in `examples/`,
  and three executed notebooks.
- A curated, cited adsorbent library (`data/adsorbents/`) and cycle templates
  (`data/cycles/`).
- MkDocs documentation with generated validation and adsorbent-library pages.
- GitHub Actions: CI on Linux, macOS and Windows, the Python wheel, the
  Emscripten build and the docs; release workflow with cibuildwheel wheels (CPython 3.10–3.15),
  sdist, source tarball and gated (Test)PyPI trusted publishing; docs
  deployment to GitHub Pages.
