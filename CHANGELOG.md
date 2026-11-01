# Changelog

All notable changes to Skarstrom are recorded here. The format follows
[Keep a Changelog](https://keepachangelog.com/en/1.1.0/) and the project uses
[Semantic Versioning](https://semver.org/). The single source of the version is
`ADSIM_PROJECT_VERSION` (plus `ADSIM_VERSION_SUFFIX` for pre-releases) in the
root `CMakeLists.txt`.

## [Unreleased]

### Fixed
- `CreateBoundaryValue` and `CreateInletTemperatureProfile` know the double-exponential
  profile. The web document layer builds profiles through the first factory, so a deck
  with a double-exponential pressure boundary (the Perez 2019 decks) lost it to a null
  profile, i.e. a constant 1 bar, in the browser. A unit test now checks every listed
  kind has a factory entry.

### Added
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
