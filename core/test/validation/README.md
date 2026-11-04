# Analytical validation suite

Tests in this directory compare the solver against **expected values derived
from physics**, not against previously-recorded output. There are no golden
files: if a test here fails, either the solver is wrong or the derivation in the
test comment is wrong, and both are checkable by reading.

This complements the existing coverage rather than duplicating it:

| Layer | Where | What it proves |
|---|---|---|
| Unit tests | `test/unit/` | Individual components behave as specified |
| Integration tests | `test/integration/` | Cycles run, and mass is conserved |
| Paper validation | `test/papers/` | Results match published experiments |
| **Analytical validation** | **`test/validation/`** | **Results match closed-form solutions of the governing equations** |

## Design rules

1. **The reference is written from the physics, not from `core/src`.** The
   helpers in `validation_helpers.h` re-derive cylinder geometry, energy
   inventories and heat fluxes from textbook forms. Do not refactor them to
   share code with the solver — that would make the tests self-confirming.
2. **Collapse the PDE to something solvable.** Most energy tests seal both ends
   of the bed so that convection and axial conduction vanish identically and
   every cell obeys the same ODE. The spatial discretisation then cannot
   influence the answer, so any discrepancy is a defect in the physics.
3. **Assert that the test can fail.** Each test contains a guard confirming the
   effect it targets is actually present and large (`ASSERT_GT(gas_share, 0.1)`,
   `ASSERT_LT(dU_bed, -100.0)`, `EXPECT_TRUE(flow_reversed)`, and so on). A
   validation test that passes because the term under test is negligible is
   worse than no test.
4. **Tolerances come from the derivation.** Where a closed form is only
   asymptotically exact, the tolerance is computed from the neglected term
   (the Henry-limit tests scale theirs by `b*p`; the invariant tests set theirs
   from the relevant equation's `abstol`) rather than picked to make the test
   pass.
5. **Prefer asserting a rate over asserting a number.** Several tests check the
   *observed order of convergence* — that halving `dx` halves the error, or that
   the polytropic work approaches the isothermal limit linearly in `n-1`. That
   is a far stronger statement than any single tolerance, and it cannot be
   satisfied by accident.

## Files

| File | Covers |
|---|---|
| `validation_helpers.h` | Independent reference implementations: energy inventories, cylinder areas, moments |
| `breakthrough_case.h` | One shared, deliberately plain isothermal breakthrough case |
| `psa_cycle_case.h` | One shared minimal O2/N2 pressurise/feed/blowdown cycle on a 13X bed |
| `test_energy_balance.cpp` | Closed and open-system energy balances, wall/ambient coupling, cyclic steady state |
| `test_analytic_transport.cpp` | Residence time, pressure drop, breakthrough time, valve blowdown |
| `test_analytic_dispersion.cpp` | Physical axial dispersion: inert outlet curve against Ogata-Banks |
| `test_retarded_dispersion.cpp` | Linear-isotherm advection-dispersion: front retarded by R, width by R^2 |
| `test_dispersion_mass_balance.cpp` | Component balance with dispersion and Danckwerts inlets, per step and at CSS |
| `test_advection_schemes.cpp` | Upwind order on a smooth solution, WENO / van Leer less diffusive, conservation, sum y = 1 |
| `test_equilibrium_shock.cpp` | Langmuir shock front speed against the equilibrium-theory secant balance |
| `test_ldf_kinetics.cpp` | LDF relaxation at rate k, first- and second-order rate laws |
| `test_extended_langmuir.cpp` | Extended (competitive) Langmuir loadings, single-component limit, selectivity |
| `test_thermodynamic_consistency.cpp` | Isosteric heat vs van 't Hoff, Henry limits, monotonicity |
| `test_adiabatic_adsorption.cpp` | Adiabatic breakthrough temperature plateau against the energy balance |
| `test_ltne_relaxation.cpp` | Two-temperature (LTNE) gas-solid relaxation rate, bed-wall energy conservation |
| `test_physics_parity.cpp` | Ergun / Blake-Kozeny pressure drop, isothermal mode, inlet temperature ramps, serialisation |
| `test_compressor_work.cpp` | Polytropic work, isothermal and reversible-adiabatic limits |
| `test_numerical_verification.cpp` | Order of accuracy, tolerance independence, scale invariance, flow symmetry |
| `test_cyclic_steady_state.cpp` | End-of-cycle profile converges to a cyclic steady state |
| `test_css_acceleration.cpp` | CSS acceleration reduces the cycles to CSS; reseeded ramps start at the step start |
| `test_invariants.cpp` | Maximum principle, composition and loading bounds |
| `test_robustness_guards.cpp` | Steep fronts and thin walls stay finite and bounded |
| `test_robustness_fuzz.cpp` | Randomised physical decks converge or fail cleanly, never silently NaN |
| `test_extreme_conditions.cpp` | Extreme pressure ratios, near-zero flow, very short / long steps, high temperature: converge or fail cleanly |
| `test_thermodynamic_bound.cpp` | Separation energy at CSS against the reversible pump work, the bed exergy change and the minimum work |
| `test_parameter_limits.cpp` | Rates of approach to parameter limits (k_LDF, dH, wall, Langmuir b, pressure ratio, trace feed); degenerate mixtures |
| `test_metamorphic_relations.cpp` | Equivalent decks give the same solution: component order, identical layers, split steps, units, cycle phase at CSS |
| `test_unit_operations.cpp` | Constant-speed evacuation, two-bed pressure equalisation at CSS, buffer-tank fill, without and with Henry adsorption |
| `test_thermal_fronts.cpp` | Inert-bed thermal front speed; steady exponential temperature profile with wall losses |
| `test_blowdown_composition_path.cpp` | Blowdown of a binary on linear isotherms: equilibrium-theory composition path y(P), departure scaling with the LDF lag |
| `test_competitive_langmuir_fronts.cpp` | Binary competitive (constant separation factor) Langmuir: roll-up plateau and both front times against equilibrium theory |
| `test_desorption_simple_wave.cpp` | Langmuir desorption: outlet times of the simple wave against the method of characteristics, first-order convergence and extrapolation |
| `test_ldf_thomas.cpp` | LDF-limited linear-isotherm front: whole outlet curve against the Anzelius/Thomas solution, with the exact upwind error |
| `psa_cycle_case.h` | The shared PSA cycle (pressurise, feed, blowdown) and its O2 product KPIs |
| `test_grid_convergence.cpp` | Cycle-KPI grid convergence at CSS: monotone convergence, observed order, GCI, all three schemes |
| `test_integrator_tolerance.cpp` | Error against reltol (BDF slope, abstol floor) on a breakthrough and a cycle; `max_step` independence |
| `test_css_convergence.cpp` | Geometric approach to CSS, KPI error left by the stopping rule, Aitken lands on the same KPIs |

## `test_energy_balance.cpp`

The energy equation had no conservation check before this suite; the mass
balance had two.

- **`QuiescentClosedBedIsExactlyPreserved`** — a closed, adiabatic bed with no
  source has zero right-hand side. The exact solution is the initial condition.
  Catches spurious source terms.
- **`ElectricHeatingGivesAnalyticTemperatureRamp`** — closed adiabatic bed with a
  known volumetric heat input. Because the vessel is rigid and sealed, the molar
  density is a constant of the motion and `dP/dt = C R dT/dt`, which turns the
  `Cp` in the temperature equation into `Cv`:

      T(t) = T0 + Q t / ( et C0 Cv + rho_b Cps ),   P(t) = C0 R T(t)

  A straight line with an analytically known slope. This pins the accumulation
  term, the P-V work term and the ideal-gas coupling. The test also asserts the
  simulated slope is *inconsistent* with the `Cp` (no work term) prediction, so
  it cannot pass for the wrong reason.
- **`HeatDrivenDesorptionMatchesLumpedModel`** — the same bed, now with a Henry
  loading and a heater that drives desorption. The spatially uniform PDE system
  reduces exactly to three coupled ODEs in `(T, C, q)`, which the test integrates
  independently with RK4 at `dt = 1e-3`. Exercises the heat-of-adsorption source,
  the gas/adsorbed-phase mole exchange and the LDF kinetics simultaneously. The
  test asserts the latent heat is at least 20 % of the energy input.
- **`ClosedBedWallAndAmbientEnergyBalanceCloses`** — first-law bookkeeping over
  bed + wall: `dU_bed + dU_wall + Q_ambient = 0`, plus the requirement that the
  bed-to-wall exchange is a pure transfer. This is the test that found the wall
  flux-normalisation bug described below.
- **`WallAmbientCoolingMatchesTwoBodyExponentialSolution`** — conservation says
  the books balance; it does not say the wall responds at the right *rate*. With
  the bed uniform and the wall axially insulated, bed + wall is a two-state
  linear system whose solution is a sum of two exponentials with analytically
  known eigenvalues. The test asserts the fast (bed/wall equilibration) and slow
  (loss to ambient) modes are separated by at least 5x, so a wrong wall heat
  capacity or area cannot hide.
- **`OpenSystemEnergyBalanceConvergesUnderGridRefinement`** — with flow the
  control volume is open: `dU = int (n_in h_in - n_out h_out) dt`. This is *not*
  expected to hold exactly on a coarse grid, because the temperature equation
  discretises convection in non-conservative form and neighbouring cells' fluxes
  do not telescope. For a single-component isobaric feed the two forms differ by
  `T d(rho Cp u)/dz`, which vanishes identically, so the residual is pure
  discretisation error. The test asserts it converges at **first order** — the
  measured sequence is 1.92 %, 0.95 %, 0.47 %, 0.23 % over 20/40/80/160 cells.
- **`SealedThermalCycleEnergyBalanceClosesAtCyclicSteadyState`** — at CSS the net
  energy accumulation over a cycle must be zero: whatever the heater puts in
  must leave through the wall. The mass balance has had this as a convergence
  criterion for a while; the energy balance had no equivalent. A sealed bed is
  used so the cycle carries no flow enthalpy and the balance is exact rather
  than discretisation-limited.

## `test_analytic_transport.cpp`

- **`InertTracerMeanResidenceTimeMatchesVoidVolume`** — the first moment of a
  step-change breakthrough curve equals `N_void / n_dot` exactly, and is
  *invariant to dispersion*. Physical dispersion, numerical smearing from the
  upwind scheme and cell count all change the shape of the curve; none of them
  move its centre of mass. Sharp test of void volume, area and the
  velocity/flow coupling; immune to everything a shape comparison would trip on.
- **`SteadyFlowPressureDropMatchesCarmanKozeny`** — the only quantitative test of
  the momentum equation. For isothermal compressible Darcy flow,
  `P_in^2 - P_out^2 = 2 N L mu R T / K` with `K = eb^3 dp^2 / (180 (1-eb)^2)`.
  Every other test in the suite is satisfied by any velocity field that conserves
  mass.
- **`LinearIsothermBreakthroughTimeMatchesCapacity`** — the mass-balance tests
  confirm that everything fed in is accounted for; they cannot see *when* it
  arrives. The stoichiometric time
  `t_st = V P [ et/(R T) + rho_b k exp(Q/RT) ] / n_dot` is the column capacity
  divided by the feed rate, and like the inert moment it is invariant to
  dispersion and to mass-transfer resistance. The `exp(Q/RT)` factor is what
  separates it from the inert case, so a sign error in the van 't Hoff exponent
  shows up here and nowhere else.
- **`ValveBlowdownFollowsAnalyticExponentialDecay`** — the valve law is linear in
  the pressure difference, so discharging a rigid isothermal vessel gives a
  first-order decay with `tau = V_void * 1e5 / (R T cv)`. Exercises the valve
  boundary condition, the coupling between a boundary molar flow and the
  interior pressure field, and the accumulation term under a 5:1 pressure swing.
  The test asserts the run really was isothermal before applying the closed form.

## `test_thermodynamic_consistency.cpp`

The heat of adsorption is a **separate user input** from the isotherm
parameters, and nothing in the solver forces them to agree. Thermodynamics does:
`q_st = R T^2 (d ln P / d T)_q`. These tests evaluate that derivative numerically
from the isotherm functions and check it against the parameter a user would
reach for. No solver run, so they are near-free.

- Henry and single-site Langmuir: `q_st` equals the fitted energy parameter
  exactly, at every coverage.
- Dual-site Langmuir: `q_st` is loading-dependent, bounded by the two site
  energies, decreasing as the strong site saturates, and tends to the
  Henry-weighted mean at zero coverage.
- **Concentration-basis Langmuir: `q_st = q1 + R T`, not `q1`.** Worth knowing —
  at 300 K that is an extra 2.5 kJ/mol, roughly 8 % of a typical CO2 heat of
  adsorption, and it lands straight in the energy balance if the heat of
  adsorption field is filled in with the fitted `q1`.
- Low-pressure Henry limits and monotonicity/saturation bounds.

## `test_compressor_work.cpp`

Compression duty is the largest operating cost in most PSA plants, so the power
model feeds every economic KPI. The polytropic work has two limits worth
pinning: `n -> 1` must recover the isothermal work `n_dot R T ln(P2/P1)` (both
the prefactor and the bracket diverge, so this is a genuine L'Hopital limit),
and `n = gamma` must give the reversible adiabatic work, which must equal the
enthalpy rise `n_dot Cp (T2 - T1)` computed from the model's own outlet
temperature. Also checks that `CalculateMolarFlow` inverts `CalculatePower`.

## `test_numerical_verification.cpp`

Verification rather than validation: not "does it solve the right equations?"
but "does it solve the equations right?". These catch defects that are present
on *every* grid and are therefore invisible to any fixed-grid comparison.

- **`UpwindSchemeConvergesAtFirstOrder`** — with dispersion off and no
  adsorption the exact breakthrough curve is a step with zero variance, so the
  variance of the simulated curve is a direct measure of numerical diffusion
  (`D_num = u dx / 2`) with a known exact value of zero. Halving `dx` must halve
  it. The test measures the observed order and requires `1.0 +/- 0.2`, while
  separately requiring the first moment *not* to drift (it is fixed by the
  column capacity and must be grid-invariant).
- **`SolutionIsIndependentOfIntegratorTolerance`** — loosening `reltol` from
  1e-9 to 1e-4 must move the outlet composition by far less than the
  discretisation error. If it does not, every result in the suite is quietly
  tolerance-dependent.
- **`ResultsAreInvariantUnderColumnScaling`** — the equations see the
  cross-section only through the superficial velocity. Doubling the diameter
  while quadrupling the feed rate must reproduce the profiles cell-for-cell.
  Catches an area, volume or diameter that has leaked into the physics with the
  wrong power.
- **`TopFedBedMirrorsBottomFedBed`** — a uniform bed has no preferred direction,
  so a top-fed run must be the mirror image of a bottom-fed one. This is the
  sharpest available test of the upwinding, which switches branch on the sign of
  the velocity; an asymmetry there is invisible in every co-current test.

## `test_invariants.cpp`

Properties the exact solution has for structural reasons, asserted over every
cell and every logged instant rather than through one scalar functional. They
catch localised misbehaviour — an overshoot at a front, a negative loading in a
cell the front has not reached — that an integral comparison averages away.

- **`TemperatureObeysDiscreteMaximumPrinciple`** — a scalar convected with no
  source can never leave the envelope of its initial and boundary values. Upwind
  preserves this; an unlimited higher-order scheme does not, and the symptom is
  a small local overshoot that no energy balance would see. Run both ways (hot
  feed into a cold bed and vice versa).
- **`MoleFractionsAndLoadingsStayBounded`** — `0 <= y <= 1`, `sum y = 1`,
  `q >= 0`, and `q` below the equilibrium value at the feed partial pressure.
- **`BoundsHoldThroughAPressureSwingCycle`** — the same bounds through a
  pressurise/feed/blowdown cycle with two competing adsorbates, which reverses
  the flow. Reversal is where boundedness is most likely to break. The test
  asserts the flow genuinely reverses.

## `test_extreme_conditions.cpp`

Targeted corners an operator reaches by turning one knob to an extreme on an
ordinary deck (the fuzz samples at random): pressure ratio 500 (50 bar to
0.1 bar) and 1.0001, blowdown to 1 mbar, feeds of 1e-9 and 1e-10 mol/s, steps
of 1 ms and of 1e5 s (which packs ~10^4 bar into the closed bed), a 1e6 s
breakthrough with no step cap, a bed at 600, 900 and 1200 K (above the 650 K
end of the viscosity and conductivity fits) and the PSA cycle at 700 K. No
accuracy is claimed. Each run must converge or fail cleanly (a negative IDA
code or an exception, within tightened work caps), every logged state must be
finite, and a converged run must be physical: P > 0, T > 0, 0 <= y <= 1 and
q >= 0 to their absolute tolerance, sum y = 1 to 1e-3 (the fuzz bound;
measured worst 4e-6, at 1 mbar). Each corner has a guard that it really is
extreme (the pressure or temperature reached). All corners converge today.
## `test_thermodynamic_bound.cpp`

At cyclic steady state the separation energy must be at least the minimum
(reversible) work of producing the product streams from the feed. For a
plant at CSS exchanging heat with surroundings at T0 and handling ideal gases
at T0 (the decks are isothermal), the second law gives
W >= sum_out n b - sum_in n b, with b = R T0 [sum y ln y + ln(P/P0)]; with all
streams at P0 this is the textbook minimum work of separation. The test
checks it as a chain of separately meaningful inequalities on the last cycle:

    W_electric  >=  W_rev  >=  dB_bed  >=  W_min

W_electric is the engine's reported energy, the integral of every endpoint's
logged pump power (the E_total of the energy KPI; polytropic, eta = 0.7,
n = 1.4); W_rev the reversible isothermal work of the same pumps on the same
logged flows (pump model consistency); dB_bed the exergy change of the
instantaneous streams at the bed ends (the bed's own second law: dispersion,
pressure drop and finite-rate adsorption can only destroy exergy); W_min the
mixed products at their delivery pressures (the ideal-mixing bound). The
slack for the not-quite-CSS cycle is R T0 (max |ln y| + max |ln P/P0|) times
the measured component imbalance (< 0.15 J).

- **`VacuumSwingEnergyExceedsMinimumSeparationWork`** — feed and product at
  1 bar, blowdown to 0.2 bar through a vacuum pump: 2486 J >= 1513 J >=
  1245 J >= 172 J (pure separation work); second-law efficiency 6.9 %.
- **`PressureSwingEnergyExceedsMinimumWork`** — 1 -> 3 bar with compressors,
  light product delivered at 3 bar (its pressure exergy counts): 15004 J >=
  9050 J >= 8691 J >= 5851 J (of which 266 J separation).
## `test_parameter_limits.cpp`

Limits as a parameter goes to zero or infinity, asserted through the RATE of
approach over a geometric sequence (observed order, tolerance from the next
term of the expansion), plus two degenerate mixtures.

- **`LdfFrontApproachesEquilibriumTheoryAsOneOverK`** — linear isotherm, no
  dispersion. The cumulants of the semi-discrete column are exact: the mean is
  t_st whatever k is, and the variance is 2 (L/u_s) H / k plus a k-free
  upwind part, so var(k) - var(2k) = (L/u_s) H / k. Measured: mean constant
  to 1e-3, order in 1/k 1.002, and the size within 1.1-1.3 % at y = 0.005,
  0.5-0.6 % at y = 0.0025; the error is O(y) (the uptake changes the gas
  velocity) and its extrapolation to y = 0 leaves 0.04 %.
- **`VanishingLdfCoefficientGivesAnInertTracer`** — k -> 0: the outlet curve
  approaches the inert tracer's at order 0.98-0.99 in k, and the deficit behind the
  front, 1 - exp(-k H L/u_s), gives deficit/k = 120.775 s extrapolated to k = 0
  against H L/u_s = 120.766 s.
- **`AdiabaticRunApproachesIsothermalAsHeatOfAdsorptionVanishes`** — the
  temperature rise and the outlet-curve shift are O(dH) (order 1.000 for
  dH = 10, 5, 2.5 kJ/mol); at dH = 0 the run matches the isothermal mode to
  1.3e-5 y_feed, 1.3 % of the smallest dH effect.
- **`WallAndBedBecomeIsothermalAsHeatCapacityAndExchangeGrow`** — the wall's
  excursion falls as 1/C (order 0.998, 1.000 per decade of wall heat
  capacity); with that wall, the bed's falls as 1/h_i once the bed-wall time
  constant is below the uptake time (order 0.99, 1.00 from h_i = 1e5 W/m2K).
- **`LangmuirRunReducesToHenryRunAsAffinityVanishes`** — b -> 0 with m b fixed:
  the simulated outlet curve approaches the Henry run's at order 1.00 in b' p.
- **`PressureSwingSeparationVanishesAsThePressureRatioGoesToOne`** — at CSS
  the heavy product's composition deviation is O(s) (order 0.92, 0.90) and the
  light product's enrichment O(s^2) (order 1.91, 1.87) for P_high/P_low =
  1 + s, s = 0.1 .. 0.025; with no swing the enrichment is 2.3e-7 (the feed
  step's own pressure drop).
- **`TraceFeedReducesToTheSingleComponentLimit`** — a trace of O2 in N2 on a
  competitive Langmuir: the carrier's loading tends to the pure-N2 run's at
  order 1.00 in y, the normalised trace curve to a limit at order ~0.95, and
  that limit's first moment (50.524 s, extrapolated) to linear chromatography
  with the trace's Henry constant at the carrier loading (50.507 s).
- **`DegenerateMixtures.IdenticalComponentsKeepTheirFeedRatio`** — O2 split
  into O2 and Ar with identical isotherm, kinetics and heat of adsorption,
  2:1, through three PSA cycles: y_O2 = 2 y_Ar to 1.3e-10 and q_O2 = 2 q_Ar to
  1.2e-12 mol/kg in every cell at every instant (tolerance one error-test
  weight), while the pair is separated from N2 (y_O2 0.04 .. 0.34).
- **`DegenerateMixtures.ZeroAffinityComponentBehavesAsAnInert`** — a Henry
  constant of zero (and a 40 kJ/mol heat of adsorption that must do nothing)
  reproduces the inert tracer to 3e-11, with the loading at 3e-19 mol/kg.
## `test_metamorphic_relations.cpp`

Pairs of runs that pose the same physical problem in two ways and must give
the same answer, compared field by field (pressure, temperature, every mole
fraction and loading, every cell, the end of every step). No reference
solution is needed. The tolerance comes from the integrator: two equivalent
runs are two adaptive-step integrations of one semi-discrete system, each
within its global error E of it, so they may differ by 2E (E measured in the
test by tightening reltol 100x) plus one error-test weight w = reltol |x| +
atol (atol fixed in the code, `docs/numerics.md`).

- **`ComponentOrderDoesNotChangeTheSolution`** — reversing the fluid's
  component list on the shared PSA cycle (two adsorbates, heat of adsorption,
  wall exchange, flow reversal). Measured: y differs by 1.4e-7 to 2.7e-7,
  against a tolerance of 1.1e-6 to 1.6e-6 (E = 3.5e-7).
- **`IdenticalLayersEqualOneLayer`** — the packing split into two and four
  identical layers (dispersion, heat of adsorption, bed-wall exchange and
  harmonic-mean face permeability all on, front across the interfaces).
  Bitwise identical on Linux.
- **`SplitStepEqualsSingleStep`** — a step cut into two halves with the same
  boundary conditions. Approximate by construction: every step restarts the
  integrator (first order, initial step, re-solved algebraic state), which may
  move a single unknown by up to sqrt(N) w (the WRMS error test over N
  unknowns), so the tolerance is 2E + sqrt(N) w. Measured on y_O2: 3.7e-8
  (0.6 w); on a longer run it followed w down to the absolute-tolerance floor
  as reltol was tightened, the signature of an integrator effect.
- **`InputsInOtherUnitsGiveTheSameSolution`** — the same deck with its inputs
  in cm, mm, bar, mbar, degC, min, g/cm3, kJ/(kg K) and mol/h. The converted
  values equal the SI ones to 1e-14 and the solutions agree (y: 6.9e-8,
  E = 6.5e-8).
- **`RotatedCycleReachesTheSameCyclicSteadyState`** — the engine solves one
  bed, so the beds of a symmetric multi-bed cycle are that bed at other phases
  and are identical by construction; the testable consequence is that the
  CSS cannot depend on which step is listed first. An isothermal
  pressurise / feed / blowdown / purge cycle, the purge fed with the feed's
  product through a connection, is started two steps later: listed before the
  feed, the purge then reads the PREVIOUS cycle's feed instead of the same
  cycle's. After 30 cycles the two CSS agree to 9.4e-8 in y (tolerance: twice
  each run's own cycle-to-cycle change, 1.7e-7), while their first cycles
  differ by 0.21.

Not testable as a relation: "one bed = two beds in series of half length".
Other beds exist only through connections that replay a partner step's
previous run, so two beds working in series at the same time cannot be
written as a deck that should equal one long bed. The layer split is the
exact statement the engine can express (a series connection with no dead
volume).
## `test_unit_operations.cpp`

The unit operations between and around beds, each reduced to a lumped ODE: the
bed is isothermal, the pressure drop negligible, and with the same Henry slope
on both gases the inventory is `V a P`, `a = eps/(R T) + rho_b K`, whatever the
composition.

- **`ConstantSpeedEvacuationDecaysExponentially`**: a volumetric-flow outlet
  (constant pump speed `S`) gives `P = P0 exp(-t/tau_S)`,
  `tau_S = V a R T / S`. Fitted over three time constants: 52.366 s against
  52.360 s inert and 247.03 s against 246.96 s with adsorption, which
  stretches the decay by the Henry retardation factor (4.7). About 0.03 s.
- **`TwoBedEqualisationConservesMolesAndDecaysAtTwiceTheRate`**: equalisation
  through a linear valve, as two connected steps of the single-bed model
  (pressurise, ED, blow down, EU) run to the cyclic steady state (30 cycles).
  The mole balance keeps `P_ED + P_EU = P_H + P_L` (to 1.2e-5 inert, 2.7e-5
  with adsorption) and the difference decays as `exp(-t/tau_E)`,
  `tau_E = 1e5 V a / (2 cv)`: fitted 31.699 s against 31.699 s inert, 149.53 s
  against 149.51 s with adsorption. The factor 2 is the signature of two beds
  exchanging gas rather than one venting to a fixed pressure. About 2.8 s.
- **`TwoBedEqualisationWithLangmuirConservesTheInventory`**: the same with
  a Langmuir isotherm (`K P = 0.5 .. 1.5`). The pressures no longer add up (the
  sum falls 9.0 % short, the guard); what is conserved is the inventory
  `n(P_ED) + n(P_EU) = n(P_H) + n(P_L)`, `n = eps P/(R T) + rho_b q*(P)`,
  measured to 1.1e-4 (tolerance 2e-3, the LDF lag). About 2.3 s.
- **`BufferTankFillDeliversTheVentedMoles`**: a fill step drawing 60 % of a
  vent step's moles from a buffer tank raises the pressure by 60 % of the
  vent's pressure drop (103 760 Pa against 103 759 Pa inert; 103 766 against
  103 757 with adsorption), linearly in time. About 0.04 s.
## `test_thermal_fronts.cpp`

Convective heat transport in a flowing bed, which the sealed-bed energy tests
cannot see.

- **`InertBedFrontTravelsAtEnergyBalanceSpeed`**: a 10 K step in feed
  temperature into an inert, adiabatic bed. At constant pressure the gas
  internal energy per volume `eps Cv P/R` is constant, so the energy and mole
  jump conditions give `w = G1 Cp / (rho_b Cps + eps C1 Cp)`. The outlet
  temperature's first moment is 347.54 s against `L/w = 347.63` s (2.5e-4);
  without the gas heat-capacity term it would be 323.1 s, so the test sees
  that term (7 % of the total) with a 1e-3 tolerance. The jump conditions
  make the characteristic speeds on both sides equal to `w` (a contact
  discontinuity), so the smeared front still travels at `w`. About 0.1 s.
- **`SteadyFlowWithWallLossGivesExponentialProfile`**: hot gas through a bed
  losing heat through the bed-wall and wall-ambient films in series reaches
  `T(z) = T_amb + (T_in - T_amb) exp(-lambda z)`,
  `lambda = U pi d_i / (n_dot Cp)`,
  `1/(U pi d_i) = 1/(h_i pi d_i) + 1/(h_o pi d_o)` (`lambda L = 1.34`). On
  25, 50 and 100 cells the profile and outlet errors fall at first order
  (outlet 9.3e-3, 4.7e-3, 2.35e-3), and the finest outlet error equals the
  leading upwind term `theta (lambda L)(lambda dx)/2 = 2.352e-3` (measured
  2.345e-3). About 0.5 s.
## `test_blowdown_composition_path.cpp`

A pressure-change step. A uniform bed of a binary on linear isotherms
(`K_A = 8 K_B`, isothermal) blows down through a valve with the other end
closed. In local equilibrium the composition stays uniform and the outflow
carries it, so `d(a_A P y)/y = d(a_B P (1-y))/(1-y)` with
`a_i = eps/(R T) + rho_b K_i`, which integrates to the Knaebel-Hill path
`P/P0 = (y0/y)^(1/(1-beta)) ((1-y)/(1-y0))^(beta/(1-beta))`,
`beta = a_B/a_A`. The path depends on the isotherm slopes only, not on the
valve or the time.

- **`FollowsEquilibriumTheoryPath`**: blowing down 5 -> 1 bar enriches the
  strong component from 0.30 to 0.7615. The largest departure of `ln P` from
  the path is 7.1e-4 at `k = 5 1/s` and 3.6e-4 at `k = 10 1/s`: it is the LDF
  lag, halving with `1/k` (ratio 2.00), and its Richardson extrapolation to
  `k -> inf` vanishes. The final composition matches the path at the final
  pressure, and the composition stays uniform to 1.3e-5. Guard: `beta < 0.3`,
  so a model without the adsorbed phase (`beta = 1`, `y` constant) fails by
  0.46 in `y`. About 0.4 s.
## `test_competitive_langmuir_fronts.cpp`

Competitive adsorption fronts. Two trace adsorbates on an extended Langmuir
isotherm with a common capacity (constant separation factor `alpha = 3`) are
fed into a clean bed. Equilibrium theory (Rhee, Aris & Amundson 1970; Ruthven
1984 ch. 9) gives two shocks with the weak component alone, ROLLED UP above
its feed level, between them. The coherence condition on the slow front (both
components' jumps give the same speed) fixes the plateau as the positive root
of `alpha x^2 - [(alpha-1) b0 - alpha + D_f] x - (alpha-1) b0 = 0`
(`x = K_B p_B*`, `b0 = K_B p_B0`, `D_f = 1 + K_A p_A0 + K_B p_B0`), and the
jump conditions give the two front times.

- **`TwoFrontsAndRollUpMatchEquilibriumTheory`**: on 200 cells the plateau of
  the weak component is `1.5363e-3` against `1.5352e-3` (a 1.54x roll-up; error
  0.073 %), the fast front's first moment is 312.39 s against 312.40 s
  (3e-5), and the slow front is 532.60 s in both the strong component and the
  weak one's drop from the plateau, against 532.72 s (2.3e-4). Tolerance
  `y_A + y_B = 2e-3` relative, an `O(y)` bound on the sorption effect, which
  was measured by varying the feed level (plateau error about
  `0.35 (y_A + y_B)`). The quantities are converged in everything numerical:
  10x in `reltol`, `max_step` 0.1 to 2 s, 100 to 400 cells and LDF 5 to
  20 1/s move them by < 1.5e-5. First moments are mass balances, so the
  numerical and LDF smearing of the shocks does not move them; the partial
  integrals interpolate at the window ends, because the solver logs at
  irregular instants. Guards: the roll-up exceeds 1.4x and shows in the
  simulated outlet, and A is absent from the plateau. About 3.4 s.

## `test_desorption_simple_wave.cpp`

The unfavourable direction of the Langmuir front that `test_equilibrium_shock`
pins: a bed saturated at a trace level (`K p0 = 3`) purged with clean carrier.
The front spreads into a simple wave in which each concentration travels at
`u / (1 + (rho_b/eps) f'(c))`, so the outlet sees `c` at
`t(c) = (L/u) (1 + M / (1 + K p)^2)` (method of characteristics; de Vault 1943,
Ruthven 1984 ch. 8), and the whole wave probes the isotherm slope rather than
its secant.

- **`LangmuirDesorptionFollowsCharacteristics`**: on 100 and 200 cells the
  outlet times of the levels `c/c0 = 0.6, 0.5, 0.4` converge at first order
  (measured 0.97, 1.00, 0.98) to `t(c)`, and the Richardson-extrapolated times
  agree with `t(c)` to 3e-4, 1e-5 and 1.2e-4 (tolerance 1e-3: LDF lag
  `1/(k t)`, the `O(y0)` sorption effect and the mole-fraction `abstol`).
  The sorption effect was measured: `1.3e-3` of `t(c)` at `y0 = 1e-3`,
  `1.3e-4` at `y0 = 1e-4`, which is why the test uses `y0 = 1e-4`. Guard: the
  wave spreads over 450 s, more than 5x its head time, and `f'(c)` varies
  16-fold across it. About 0.5 s.
## `test_ldf_thomas.cpp`

A front shaped by the LDF rate rather than by dispersion: a Henry isotherm with
`k = 0.1 1/s` (12 transfer units), dispersion off, isothermal, compared point
by point with the Anzelius/Thomas closed form
`c/c0 = J(x, tau) = 1 - int_0^x exp(-tau - s) I0(2 sqrt(tau s)) ds`
(Ruthven 1984, ch. 8), which the test evaluates by quadrature with its own
Bessel series. The upwind scheme makes the column `n` tanks in series, whose
transfer function `(1 + (dx/u) g(s))^(-L/dx)` is known exactly, so the grid
error is not a nuisance to be tolerated but a second prediction.

- **`ReferenceSolutionSatisfiesItsIdentities`**: the two integral forms of `J`
  agree, `J(a,b) + J(b,a) = 1 + exp(-a-b) I0(2 sqrt(ab))`, the analytic second
  derivatives match finite differences, and the cumulants are `x` and `2x`. No
  solver run.
- **`LinearIsothermOutletCurveMatchesThomasSolution`**: on 50, 100 and 200
  cells, (1) the outlet variance minus the exact upwind term `L dx R^2 / u^2`
  equals the Thomas variance `2 (R-1) L / (u k)` (measured 2415.5 vs
  2416.2 s^2, 0.03 %, on every grid); (2) `max |F - J|` falls at first order
  (0.0246, 0.0128, 0.0065); (3) after subtracting the closed-form leading error
  `(L dx / 2) d2F/dz2`, the residual falls at second order (2.2e-3, 6.1e-4,
  1.7e-4; orders 1.89, 1.84) and is 2.6 % of the error on the finest grid.
  Guard: a Gaussian curve with the same mean and variance misses the
  simulated one by 0.042, 250x the finest residual, so the asymmetric Thomas
  shape is what is being matched. About 1.3 s.

## Convergence of cycle results

`test_grid_convergence.cpp`, `test_integrator_tolerance.cpp` and
`test_css_convergence.cpp`. The tests above verify the discretisation on
single breakthrough curves. What a user reports is a cycle KPI at cyclic steady state, which carries three
numerical errors at once: the grid, the integrator tolerance and the distance
from CSS. These files measure each one on the shared PSA cycle
(`psa_cycle_case.h`, with the O2 product KPIs of `AddO2ProductKpis`: purity,
recovery and the specific energy of a feed compressor). The full grid study
of the example decks, and what each needs, is in `docs/numerics.md`
("Convergence of the results").

- **`CycleKpiGridConvergence.*`** — the isothermal cycle at CSS on 10, 20, 40
  and 80 cells with each advection scheme. The successive differences must
  keep one sign and shrink (monotone convergence), and the observed order from
  the three finest grids (Richardson) must lie in `[1, 2]` for upwind: with
  error `a h + b h^2`, same-signed terms put the observed order between the
  formal order and one more (measured 1.54 to 1.56; offline on 80 to 640 cells
  it falls to 1.37 and 1.26). WENO and van Leer must converge faster than
  upwind (larger order, smaller difference at every grid; measured 2.40 to
  2.50). Their asymptotic order is capped at 1 by the first-order bed-end faces
  (`test_advection_schemes.cpp`), but the KPIs are read at the feed end, where
  the face carries the feed exactly, and at the product end, where the
  gradient is small, so that term has a small coefficient. Finally Roache's
  grid convergence index (safety factor 1.25) of the 80-cell upwind KPI must
  be a true error band: it must contain the WENO and van Leer 80-cell values,
  which are ten times closer to the limit (upwind purity 0.39572 +- 0.00046,
  WENO 0.39611, van Leer 0.39611).
- **`IntegratorTolerance.*`** — IDA's error-per-step control of a BDF method of
  order k gives a global error `~ tol^(k/(k+1))`, so the observed slope of
  log(error) against log(reltol) must lie in `[0.5, 1]`: measured 0.68 on a
  breakthrough (1e-3 to 1e-7, against a 1e-12 reference) and 0.77 on the
  cycle's end-of-cycle loadings. Below `reltol |y| = abstol` the error stops
  falling: in a dilute breakthrough the loading error stays at 0.75 to 1.4 x
  abstol (1e-8 mol/kg) from reltol 1e-8 to 1e-10. The cycle KPIs converge more
  slowly (recovery slope 0.42): they integrate the logged boundary flows by
  the trapezoid rule, and a step that opens a pressure boundary drives a flow
  transient on the pressure diffusion time of one cell (~2e-5 s), resolved
  only by IDA's own steps. At the default reltol 1e-6 the recovery error is
  6e-5 after two cycles, below the 2.4e-4 change from 20 to 40 cells.
  `max_step` caps at or above IDA's own steps leave a run bit-for-bit the same;
  a 0.1 s cap (2 to 2.3 x the steps) moves profiles and KPIs by less than the
  tolerance error of a 100 x looser reltol. Decks that set `max_step` do it
  for robustness, not accuracy (Vysyaraju 2024: IDA stalled at its minimum
  step with a 1 s cap; 0.5 s and reltol 1e-7 at 1 s agree to 0.01 points).
- **`CSSConvergence.*`** — linearised at its fixed point the end-of-cycle map
  contracts by its slowest eigenvalue rho, so the state change decays
  geometrically, every KPI converges at the same rate, and the error left
  when the stopping rule fires at a state change `Delta <= eps` is
  `g Delta rho / (1 - rho)` (g the KPI change per unit state change). On the
  isothermal cycle (rho = 0.106) the successive ratios agree with the fitted
  rate within 10 % over five decades, the KPI error is a fixed multiple of the
  state change, and each decade of state tolerance costs one cycle and buys a
  decade of KPI accuracy. With the energy balance on (rho = 0.966, a slow
  thermal mode) the CSS recovery fitted to cycles 40-79 predicts the KPI of
  cycles 99 to 139 within the noise, and the default rule (state tolerance
  2e-5) stops at cycle 108 with a recovery error of 6.11e-5 against 6.21e-5
  predicted (purity 7.20e-5 against 7.29e-5), inside a tolerance of 3 standard
  deviations computed from the noise measured in the run (11 % and 4 %).
  Aitken acceleration stops at cycle 83 on KPIs 2.6e-5 (recovery) and 4.0e-5
  (purity) from that CSS value, inside the error the rule allows a plain run.
  The KPIs carry cycle-to-cycle noise that depends on the integrator's step
  sequence, and so on the platform; the tests measure it in the run and fit
  every CSS value over many cycles (an Aitken extrapolation from three
  cycles amplifies the noise five-fold: on macOS that moved f* by 15 % of the
  error under test).

## Fixed along the way

### Wall flux normalisation

`ClosedBedWallAndAmbientEnergyBalanceCloses` originally failed by ~60 %. The
bed-side flux was correct — the bed's internal-energy loss matched
`pi di dx hi (T - Tw)` to within 2 % — but `UpdateWallTemperature` normalised its
fluxes by `(d_i + t_w)^2 - d_i^2` instead of the annulus area factor
`d_o^2 - d_i^2 = 4 t_w (d_i + t_w)`, about half as large for a thin wall. The
numerators were also inconsistent with each other: `flux_i` was written in
diameter form (`4 d_i h_i / (d_o^2 - d_i^2)`) and `flux_o` in radius form
(`2 r_o h_o / (r_o^2 - r_i^2)`), which differ by a factor of four.

At `d_i = 0.1 m`, `t_w = 5 mm` the wall absorbed **2.05x** the heat the bed gave
up and rejected **0.51x** the physical Newton-cooling flux to ambient. The
measured wall internal-energy change agreed with those factors to five
significant figures.

Fixed in `src/equations/adsorbent_equations.cpp`: both fluxes now share the
annulus denominator `d_o^2 - d_i^2` and `flux_o` uses the diameter form, so the
wall's gain equals the bed's loss and the ambient loss equals
`pi d_o dx h_o (Tw - T_amb)`.

### Total-mole balance was a no-op

`IsMolesConserved` in `test/fixtures/user_data_test.h` had its `return false`
commented out: it printed a diagnostic and then reported success regardless. The
reason it could not be enabled was the boundary term — it used `GetTotalMolesIn`,
which evaluates `u*C` with `C` from the ghost cell and ignores axial dispersion,
rather than the discrete flux the species equations actually use. It now sums the
per-component fluxes from `GetTotalComponentMolesIn` (upwinded convection plus
the boundary dispersive flux), which closes to the same 0.2 % tolerance as the
component check, and the assertion is live.

## Known model behaviour worth a decision

`ReactorEquations` defaults `wall_temperature_boundary_condition_` to
`WallTemperatureBoundaryConditionFixed` at **293 K**, independent of
`bed.ambient_temperature_`. On a hot column the ghost cells sit far below the bed
temperature (186 K was observed for a 400 K bed) and axial conduction drains a
large amount of energy through the two end cells. The energy-balance tests set
`ZERO_FLUX` explicitly to close the control volume. Worth deciding whether
`FIXED @ 293 K` is the intended default, and whether it should at least track the
ambient temperature.

## Further work

The obvious remaining gap is an **open-system energy balance at cyclic steady
state** with real flow — the CSS test here uses a sealed bed to avoid the flow
enthalpy bookkeeping, and the open-system test runs a single step. Combining them
would need the boundary enthalpy integrals to be tracked across valves and
connections, which is a larger piece of plumbing than the rest of this suite.

Two smaller ones: an **IAST / competitive-isotherm consistency** check (the
extended Langmuir forms should reduce to the single-component isotherm at
`y = 1`, and should satisfy the Gibbs adsorption isotherm), and a **second-order
scheme comparison** if a limiter is ever added — `UpwindSchemeConvergesAtFirstOrder`
is written so that the expected order is a single number to change.
