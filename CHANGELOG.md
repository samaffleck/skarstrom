# Changelog

All notable changes to Skarstrom are recorded here. The format follows
[Keep a Changelog](https://keepachangelog.com/en/1.1.0/) and the project uses
[Semantic Versioning](https://semver.org/). The single source of the version is
`ADSIM_PROJECT_VERSION` (plus `ADSIM_VERSION_SUFFIX` for pre-releases) in the
root `CMakeLists.txt`.

## [Unreleased]

## [0.9.0] — 2026-11-07

### Added — estimated distance to cyclic steady state
- Every CSS check now estimates how far the run still is from the cyclic steady state: the convergence rate ρ from the decay of the last few end-of-cycle state changes, and the remaining distance ≈ last change × ρ/(1−ρ), in the norm of `css.state_tolerance` and as a multiple of it. It is reported in the solver's CSS log lines, as `distance_estimated`, `convergence_rate`, `distance`, `distance_ratio` and `distance_note` on `sim.css_history()` records, in the `css` entries of `skarstrom run --json` and in `css.csv`. The cycles right after an acceleration jump, and windows that are not a monotone decay, report "not estimated" with the reason.
- A run that stops at CSS more than 100 state tolerances from it (by this estimate) logs a warning. The stopping rule is unchanged: no result moves.

### Changed — the gas heat capacity is the ideal-gas Cp0(T) at the local temperature
- The default gas heat capacity is now each component's ideal-gas heat capacity C_p0,i(T), evaluated at the local gas temperature in every cell. It comes from the NASA 7-coefficient polynomials of McBride, Gordon and Reno (1993, NASA TM-4513), with the source of every species in `ideal_gas_heat_capacity.h`. It replaces CoolProp's value evaluated once at the initial-condition pressure and temperature and held for the run. The initial temperature is now only a start: editing it no longer changes the cyclic steady state (#145), and a TSA cycle's gas heat capacity follows its temperature.
- The polynomials cover 33 of the engine's 34 components; Air is built from N2, O2 and Ar at CoolProp's air composition. Each one agrees with CoolProp's ideal-gas Cp0 within 1 % (light and permanent gases) or 2 % (others) from 200 to 800 K; the table is in the Equation of State page. Hexane, which has no TM-4513 entry, keeps the held CoolProp value, and a run-start warning says so.
- **Results move** for every deck in the default mode. The example decks are re-saved.
- The registry name of the mode is now `ideal_gas` (`fluid.gas_heat_capacity_mode`: `ideal_gas`, `constant`); `coolprop` is still accepted as its former name. In C++, `GasHeatCapacityMode::IDEAL_GAS`, with `COOLPROP` kept as an alias.
- The *Gas phase* adsorbed-phase heat capacity follows the adsorbed phase's temperature. Kirchhoff heats of reaction integrate the heat capacities exactly.
- `Fluid::MolarHeatCapacity`, `MixtureHeatCapacity` and `AdsorbedPhaseHeatCapacity` take the temperature; new `MolarEnthalpyChange`, `HeatCapacityPolynomial` and `ComponentsWithHeldHeatCapacity`. Python: `skarstrom.ideal_gas_heat_capacity(component, T)`.
- `CSSStartIndependence` gains hot-deck and cold-deck starts: the deck's own initial temperatures edited by ±40 K.

### Added — a layer's mass fractions can be expressions
- A layer's mass fractions are expressions over the deck variables, like other deck values: a two-solid blend can be written as `x` and `1 - x` with a deck variable `x`, and a sweep or a sensitivity study varies `variable[x]`. The registry path `layer[ID].mass_fraction[ID]` is now an Expression parameter, so `--set`, `sim[...] = "1 - x"` and `layer.set_composition({"13X": "x", "Glass beads": "1 - x"})` accept text. The composition check uses the evaluated fractions, shows the formulas with their values when the sum is off, and reports a fraction whose formula does not evaluate. The blend edits (`AddSolidToLayer`, `SetLayerMassFraction`, `RemoveSolidFromLayer`) refuse a layer whose fractions include a formula, with nothing changed.
- `MixedAdsorbentLayer` archive version 5: the fractions are written as expressions. Older decks load with their numbers unchanged, bit for bit.

### Changed — paper decks use their paper's gas heat capacity
- The Casas 2012, Ferré 2026, Nguyen 2023, Estupiñan Perez 2019, Vysyaraju 2024, Ward & Pini 2022 and Ward & Pini 2024 decks now use the constant gas heat capacity mode with the paper's value: the printed constant where there is one (30.7 J/(mol K); per-kilogram NIST values converted with each component's molar mass; Ward & Pini 2022's single 840 J/(kg K) given to each component per kilogram, so every CO₂/He mixture has 840 J/(kg K)), and the mean of the ideal-gas Cp0 over 298–373 K for Casas 2012, which states that window but prints no value. Nguyen 2023 also takes the adsorbed-phase heat capacity from the gas, as the paper does. Arora & Hasan 2021 and Nguyen 2022 keep the default. The per-paper table is in `core/test/papers/README.md`.
- **Results move** slightly for these decks (at most 0.06 pp in a purity or recovery); every validation bound holds. The example decks and the example scripts are updated.

### Added — warm-started optimisation and sensitivity studies
- `CycleStateCache` and `RunWithStateCache` (`adsim/solver/state_cache.h`): a cache of cyclic steady states keyed by the decision vector, owned by the caller. Each run starts from the state of the **nearest** point evaluated so far (Euclidean distance in the variables scaled to [0, 1] by their bounds), through the ordinary warm-start path (grid remap, stretched histories). A state that does not fit the deck (another structure) is skipped with a note, a warm-started run that fails is repeated cold, and only runs that reached CSS are cached, so decks without a CSS check and fixed-cycle studies run cold as before. Python: `skarstrom.StateCache` and `skarstrom.optimisation.run_warm`.
- `OptimizeSimulation` (`adsim/solver/optimisation.h`) and `skarstrom.optimisation.optimise`: the whole single-objective loop over registry paths (`GlobalOptimizer`, MaxLIPO with a trust region), the cost from the candidates' KPIs. Candidates are warm-started from the cache by default (`warm_start=False` turns it off); the optimum is re-run from the deck's initial conditions and that cold value is the reported cost, next to the search's own value (`search_cost`). Every evaluation reports whether it was warm, its cycles and wall time. `GlobalOptimizer` also takes plain bounds and an objective of the candidate vector, and exposes `CurrentPoint()` and `GetBestValue()`.
- Sensitivity studies warm-start each run from the nearest earlier run's CSS, keyed by the run's probability coordinates (`sa::EvaluateSample` overload with a cache; Python `Study(warm_start=True)` by default, one cache per worker process; CLI `skarstrom sensitivity --no-warm-start`). Samples report `warm`, `cycles_run`, `css_reached` and `note`.
- Measured savings and the agreement with cold runs are in the cyclic-steady-state documentation (*Warm-started optimisation and studies*). Runs outside optimisation and studies are unchanged.

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
