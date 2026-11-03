# Changelog

All notable changes to Skarstrom are recorded here. The format follows
[Keep a Changelog](https://keepachangelog.com/en/1.1.0/) and the project uses
[Semantic Versioning](https://semver.org/). The single source of the version is
`ADSIM_PROJECT_VERSION` (plus `ADSIM_VERSION_SUFFIX` for pre-releases) in the
root `CMakeLists.txt`.

## [Unreleased]

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
