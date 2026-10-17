# Verification

**Verification** asks whether the solver solves its equations correctly; [Validation](validation.md) asks whether those are the right equations for real columns. Skarstrom's verification tests compare the solver with closed-form solutions of the governing equations, conservation identities and observed orders of convergence, so a failure means a defect, not a modelling choice. They live in [`adsim/test/validation/`](https://github.com/samaffleck/skarstrom/tree/main/adsim/test/validation) and run with the rest of the test suite (`ctest -R "EnergyBalance|Analytic|Numerical|Invariant|Thermodynamic|Compressor|Shock|LDF|LTNE|Adiabatic|ExtendedLangmuir|CyclicSteadyState|Robustness"`).

The rest of this page is the suite's README.

## Analytical verification suite

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

### Design rules

1. **The reference is written from the physics, not from `adsim/src`.** The
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

### Files

| File | Covers |
|---|---|
| `validation_helpers.h` | Independent reference implementations: energy inventories, cylinder areas, moments |
| `breakthrough_case.h` | One shared, deliberately plain isothermal breakthrough case |
| `test_energy_balance.cpp` | Closed and open-system energy balances, wall/ambient coupling, cyclic steady state |
| `test_analytic_transport.cpp` | Residence time, pressure drop, breakthrough time, valve blowdown |
| `test_thermodynamic_consistency.cpp` | Isosteric heat vs van 't Hoff, Henry limits, monotonicity |
| `test_compressor_work.cpp` | Polytropic work, isothermal and reversible-adiabatic limits |
| `test_numerical_verification.cpp` | Order of accuracy, tolerance independence, scale invariance, flow symmetry |
| `test_invariants.cpp` | Maximum principle, composition and loading bounds |

### `test_energy_balance.cpp`

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

### `test_analytic_transport.cpp`

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

### `test_thermodynamic_consistency.cpp`

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

### `test_compressor_work.cpp`

Compression duty is the largest operating cost in most PSA plants, so the power
model feeds every economic KPI. The polytropic work has two limits worth
pinning: `n -> 1` must recover the isothermal work `n_dot R T ln(P2/P1)` (both
the prefactor and the bracket diverge, so this is a genuine L'Hopital limit),
and `n = gamma` must give the reversible adiabatic work, which must equal the
enthalpy rise `n_dot Cp (T2 - T1)` computed from the model's own outlet
temperature. Also checks that `CalculateMolarFlow` inverts `CalculatePower`.

### `test_numerical_verification.cpp`

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

### `test_invariants.cpp`

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
