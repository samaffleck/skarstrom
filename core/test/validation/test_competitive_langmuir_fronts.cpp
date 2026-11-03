/*
 * Skarstrom - pressure and temperature swing adsorption simulator
 * Copyright (c) 2024-2026 Sam Affleck and contributors
 *
 * Licensed under the terms in the LICENSE file at the repository root.
 */

// Analytical validation of COMPETITIVE adsorption fronts: equilibrium theory
// for a binary extended-Langmuir mixture with a constant separation factor.
//
// test_equilibrium_shock pins a single Langmuir shock; test_extended_langmuir
// pins the competitive isotherm in a closed bed. Neither sees what competition
// does to a breakthrough: two fronts, and the weaker component ROLLED UP above
// its feed concentration between them, because the stronger one displaces it.
//
// Theory (Rhee, Aris & Amundson 1970; Ruthven 1984, ch. 9; LeVan & Carta, in
// Perry's, sec. 16). Two adsorbates A (strong) and B (weak) at trace levels in an
// inert carrier, at constant pressure and temperature, so the interstitial
// velocity u is constant. The extended Langmuir isotherm with a COMMON capacity
// m has a constant separation factor alpha = K_A / K_B:
//
//     q_i = m K_i p_i / (1 + K_A p_A + K_B p_B).
//
// Feeding (p_A0, p_B0) into a clean bed gives, in local equilibrium, two
// shocks (both transitions are favourable). Ahead of the fast front the bed is
// clean; between the fronts the gas holds B only, at a plateau p_B*; behind the
// slow front it holds the feed. A shock between states 1 and 2 travels at
//
//     w = u / (1 + (rho_b/eps) [q_i]/[c_i])   for EVERY component i
//
// ([.] is the jump, c_i = p_i / RT). For the slow front the two components'
// jumps must give the same speed (the coherence condition), which fixes the
// plateau. With x = K_B p_B*, b0 = K_B p_B0 and D_f = 1 + K_A p_A0 + K_B p_B0,
//
//     q_A(f)/c_A0 = (q_B(f) - q_B(*)) / (c_B0 - c_B*)
//   =>  alpha x^2 - [ (alpha - 1) b0 - alpha + D_f ] x - (alpha - 1) b0 = 0,
//
// whose positive root is the plateau. The outlet times of the two fronts are
//
//     t_1 = (L/u) (1 + M / (1 + x)),     t_2 = (L/u) (1 + M alpha / D_f),
//
// with M = (rho_b/eps) m K_B R T (the B Henry retardation minus one).
//
// What is measured: the outlet plateau of B (its roll-up), and each front's
// time as the first moment of its outlet curve. The first moment of a shock's
// curve is its equilibrium-theory time however much the front is smeared
// (numerical diffusion, LDF), because it is a mass balance: the moles a front
// has stored behind it. So the moments are exact up to the neglected terms:
//
//   * the sorption effect: the carrier velocity changes by O(y) as the
//     adsorbates leave the gas. Measured by varying the feed level (y_A = y_B
//     = 1e-4, 1e-3, 3e-3), the plateau error is +0.011 %, +0.073 % and +0.21 %,
//     i.e. about 0.35 (y_A + y_B); the front times move by less than 3e-4;
//   * front overlap: the plateau lasts about 220 s, and the test asserts that
//     A is absent from its middle third, so the fronts' tails do not reach the
//     mid-plateau split.
//
// The measured quantities are converged in everything numerical: changing
// reltol 10x either way, max_step from 0.1 to 2 s, the grid from 100 to 400
// cells or the LDF rate from 5 to 20 1/s moves the plateau by < 1e-6 relative
// and the front times by < 1.5e-5. (Integrals over part of the run must
// interpolate at the window ends: the solver logs at irregular instants, every
// IDA step plus a minimum resolution, so clipping a window to whole samples
// drops a sliver that depends on the step sequence; that made an earlier
// version of this test move by 0.3 % between Linux and Windows.)
//
// Tolerance: 2e-3 relative on every quantity, an O(y) bound with coefficient
// 1 on y_A + y_B = 2e-3 (about 3x the measured sorption effect on the
// plateau); measured errors are 7.3e-4 (plateau), 3e-5 (t1) and 2.3e-4 (t2).
//
// GUARDS: the roll-up must be large (plateau > 1.4 p_B0) and visible in the
// simulated outlet, and the two fronts must be well separated, so neither the
// competition nor the coherence condition can be missing.

#include "validation_helpers.h"

#include <algorithm>

namespace adsim {

  using namespace adsim::validation;

  namespace {

    constexpr double kT = 298.0;           // [K]
    constexpr double kP = 1.0e5;           // [Pa]
    constexpr double kEps = 0.4;           // bed voidage (ep = 0)
    constexpr double kRhoS = 1000.0;       // solid density [kg/m3]
    constexpr double kL = 1.0;             // [m]
    constexpr double kD = 0.1;             // [m]
    constexpr double kFlow = 0.005;        // [mol/s]
    constexpr double kYA = 1.0e-3;         // feed mole fractions (trace)
    constexpr double kYB = 1.0e-3;
    constexpr double kKA = 0.015;          // [1/Pa]: K_A p_A0 = 1.5
    constexpr double kKB = 0.005;          // [1/Pa]: K_B p_B0 = 0.5 (alpha = 3)
    constexpr double kM = 20.0;            // B Henry retardation minus one
    constexpr double kLdf = 5.0;           // [1/s], fast
    constexpr double kRunTime = 800.0;     // [s]
    constexpr double kTol = 1.0 * (kYA + kYB);  // O(y) sorption-effect bound, see above

    const char* kCarrier = "Nitrogen";
    const char* kA = "CarbonDioxide";      // strong
    const char* kB = "Methane";            // weak

    double RhoB() { return (1.0 - kEps) * kRhoS; }
    // Common capacity m [mol/kg] that gives (rho_b/eps) m K_B R T = kM.
    double Capacity() { return kM / ((RhoB() / kEps) * kKB * kR * kT); }

    struct Theory {
      double x;        // K_B p_B* at the plateau
      double t1, t2;   // outlet times of the fast (B) and slow (A) fronts [s]
    };

    Theory EquilibriumTheory(double u) {
      const double alpha = kKA / kKB;
      const double b0 = kKB * kP * kYB;
      const double Df = 1.0 + kKA * kP * kYA + b0;
      const double bq = (alpha - 1.0) * b0 - alpha + Df;
      const double x = (bq + std::sqrt(bq * bq + 4.0 * alpha * (alpha - 1.0) * b0)) / (2.0 * alpha);
      const double tL = kL / u;
      return {x, tL * (1.0 + kM / (1.0 + x)), tL * (1.0 + kM * alpha / Df)};
    }

    struct BinaryRun {
      std::vector<double> t, yA, yB;
      double u = 0.0, dP_rel = 0.0;
    };

    BinaryRun RunBinaryBreakthrough(int n_cells) {
      UserData data;
      data.bed.wall_.thickness_ = 0.005;
      data.fluid.AddComponent(kCarrier);
      data.fluid.AddComponent(kA);
      data.fluid.AddComponent(kB);

      auto ads = data.adsorbent_directory.Add("Packing");
      auto lid = data.bed.layers_.Add("Packing");
      auto& layer = data.bed.layers_.Get(lid);
      layer.momentum_model_ = MomentumModel::KarmenCozeny;
      layer.domain_.length = kL;
      layer.domain_.n_cells = n_cells;
      layer.eb_ = kEps;
      layer.particle_diameter_ = 0.002;
      layer.composition_[ads] = 1.0;
      std::get<CylinderBed>(data.bed.shape_).diameter = kD;
      auto disp = std::make_unique<DispersionModelConstant>();
      disp->dispersion_coefficient_ = 0.0;
      layer.equations_.dispersion_model_ = std::move(disp);

      auto& adsorbent = data.adsorbent_directory.Get(ads);
      adsorbent.ep_ = 0.0;
      adsorbent.solid_density_ = kRhoS;
      adsorbent.heat_capacity_ = 900.0;
      adsorbent.thermal_conductivity_ = 0.0;
      auto iso = std::make_unique<ExtendedLangmuirIsotherm>();
      for (const auto& [comp, K] : {std::pair{kA, kKA}, std::pair{kB, kKB}}) {
        iso->params_[comp].m1 = Capacity();
        iso->params_[comp].b1 = K;
        iso->params_[comp].q1 = 0.0;
      }
      adsorbent.adsorption_source_.AddIsotherm({kA, kB}, std::move(iso));
      for (const char* comp : {kA, kB}) {
        auto ldf = std::make_unique<LinearDrivingForce>();
        auto mtc = std::make_unique<MassTransferCoefficientConstant>();
        mtc->constant_mtc_ = kLdf;
        ldf->mtc_model_ = std::move(mtc);
        adsorbent.adsorption_source_.mass_transfer_models_[comp] = std::move(ldf);
        adsorbent.adsorption_source_.heat_of_adsorption_[comp] = 0.0;
      }
      data.bed.equations_.energy_balance_mode_ = EnergyBalanceMode::ISOTHERMAL;

      data.ic.pressure = kP;
      data.ic.temperature = kT;
      data.ic.wall_temperature = kT;
      data.ic.velocity = 0.0;
      data.ic.mole_fraction[kCarrier] = 1.0;
      data.ic.mole_fraction[kA] = 0.0;
      data.ic.mole_fraction[kB] = 0.0;

      data.cycle.AddStep("Feed");
      Step& step = data.cycle.GetStep("Feed");
      step.step_time_ = kRunTime;

      auto inlet = std::make_unique<BoundaryStream>();
      auto flow = std::make_unique<MolarFlowBoundaryValue>();
      static_cast<ConstantBoundaryValue*>(flow->profile_.get())->constant_value_ = kFlow;
      inlet->boundary_flow_ = std::move(flow);
      inlet->mole_fraction[kCarrier] = 1.0 - kYA - kYB;
      inlet->mole_fraction[kA] = kYA;
      inlet->mole_fraction[kB] = kYB;
      inlet->temperature_ = kT;
      auto inlet_bc = std::make_unique<FlowBoundaryCondition>();
      inlet_bc->direction_ = BoundaryDirection::INLET;
      inlet_bc->location_ = BoundaryLocation::BOTTOM;
      inlet_bc->streams_.push_back(std::move(inlet));
      step.bottom_boundary_condition_ = std::move(inlet_bc);

      auto outlet_bc = std::make_unique<PressureBoundaryCondition>();
      outlet_bc->direction_ = BoundaryDirection::OUTLET;
      outlet_bc->location_ = BoundaryLocation::TOP;
      auto outlet_p = std::make_unique<ConstantBoundaryValue>();
      outlet_p->constant_value_ = kP;
      outlet_bc->pressure_ = std::move(outlet_p);
      outlet_bc->mole_fraction_[kCarrier] = 1.0 - kYA - kYB;
      outlet_bc->mole_fraction_[kA] = kYA;
      outlet_bc->mole_fraction_[kB] = kYB;
      step.top_boundary_condition_ = std::move(outlet_bc);

      data.cycle.number_of_cycles_ = 1;
      data.settings.reltol = 1e-9;
      data.settings.max_step = 0.5;
      data.settings.max_log_points_per_step_ = 800;  // at least one snapshot per second

      CoupledSolver solver{data};
      std::string error;
      EXPECT_TRUE(data.PerformChecks(error)) << error;
      EXPECT_GE(solver.IntegrateCycle(data.cycle, false), 0)
        << "Solver did not complete (n_cells = " << n_cells << ")";

      const auto& log = solver.GetLog();
      const auto& dom = log.GetDomain();
      BinaryRun run;
      const auto P = log.GetVariable(0, 0, VariableType::PRESSURE);
      const auto& Pf = P.values.back();
      run.dP_rel = std::abs(Pf[1] - Pf[Pf.size() - 2]) / kP;
      const auto vel = log.GetVariable(0, 0, VariableType::VELOCITY);
      const auto& uf = vel.values.back();
      double u_sum = 0.0;
      int u_cnt = 0;
      for (size_t i = 1; i + 1 < uf.size(); ++i) { u_sum += uf[i]; ++u_cnt; }
      run.u = (u_sum / u_cnt) / dom.et[1];

      for (const auto& [name, series] : {std::pair{kA, &run.yA}, std::pair{kB, &run.yB}}) {
        const auto y = log.GetComponentFaceVariable(0, 0, ComponentVariableType::MOLE_FRACTION, name);
        run.t = y.times;
        for (const auto& snap : y.values) series->push_back(snap.empty() ? 0.0 : snap.back());
      }
      return run;
    }

    // int_{t_a}^{t_b} f dt over a logged series: trapezoids on the samples
    // inside the window plus the partial intervals at both ends, with f
    // interpolated linearly at t_a and t_b. The solver logs at irregular
    // instants, so clipping the window to whole samples would drop a sliver
    // at each end that depends on the integrator's step sequence.
    double Integral(const std::vector<double>& t, const std::vector<double>& f,
                    double t_a, double t_b) {
      double s = 0.0;
      for (size_t i = 1; i < t.size(); ++i) {
        const double lo = std::max(t[i - 1], t_a), hi = std::min(t[i], t_b);
        if (hi <= lo) continue;
        const double slope = (f[i] - f[i - 1]) / (t[i] - t[i - 1]);
        s += 0.5 * ((f[i - 1] + slope * (lo - t[i - 1])) + (f[i - 1] + slope * (hi - t[i - 1]))) * (hi - lo);
      }
      return s;
    }

  } // namespace

  TEST(CompetitiveLangmuirFronts, TwoFrontsAndRollUpMatchEquilibriumTheory) {
    const BinaryRun run = RunBinaryBreakthrough(200);
    ASSERT_FALSE(run.t.empty());
    const Theory th = EquilibriumTheory(run.u);
    const double yB_star = th.x / (kKB * kP);  // plateau mole fraction of B

    // --- the derivation's assumptions ------------------------------------
    ASSERT_LT(run.dP_rel, 2e-3) << "Pressure drop too large for a constant velocity";
    ASSERT_NEAR(run.yA.back(), kYA, 1e-3 * kYA) << "A has not broken through";
    ASSERT_NEAR(run.yB.back(), kYB, 1e-3 * kYB) << "B has not returned to its feed level";

    // GUARDS: a large roll-up, and two well separated fronts.
    ASSERT_GT(yB_star, 1.4 * kYB) << "Roll-up too small to test the coherence condition";
    ASSERT_GT(th.t2, 1.5 * th.t1) << "Fronts too close together";
    const double yB_max = *std::max_element(run.yB.begin(), run.yB.end());
    ASSERT_GT(yB_max, 1.4 * kYB) << "The simulated outlet shows no roll-up of the weak component";

    // --- 1. the plateau (the roll-up) -------------------------------------
    // Mean over the middle third of the plateau.
    const double t_a = th.t1 + (th.t2 - th.t1) / 3.0;
    const double t_b = th.t2 - (th.t2 - th.t1) / 3.0;
    const double yB_plateau = Integral(run.t, run.yB, t_a, t_b) / (t_b - t_a);
    std::cout.precision(8);
    std::cout << "[CompetitiveLangmuir] u=" << run.u << " plateau yB=" << yB_plateau
              << " (theory " << yB_star << ", roll-up " << yB_star / kYB << ")";
    EXPECT_NEAR(yB_plateau, yB_star, kTol * yB_star)
      << "Plateau of the weak component " << yB_plateau << " != equilibrium-theory roll-up "
      << yB_star << " (feed " << kYB << ")";
    // A is absent from the plateau.
    EXPECT_LT(Integral(run.t, run.yA, t_a, t_b) / (t_b - t_a), 1e-3 * kYA);

    // --- 2. the two front times (first moments) ---------------------------
    // Fast front: B rises from 0 to the plateau; split mid-plateau.
    const double t_mid = 0.5 * (th.t1 + th.t2);
    std::vector<double> sB(run.t.size()), sA(run.t.size()), dB(run.t.size());
    for (size_t i = 0; i < run.t.size(); ++i) {
      sB[i] = 1.0 - run.yB[i] / yB_plateau;
      sA[i] = 1.0 - run.yA[i] / kYA;
      dB[i] = (run.yB[i] - kYB) / (yB_plateau - kYB);
    }
    const double t1_meas = Integral(run.t, sB, 0.0, t_mid);
    const double t2_meas = Integral(run.t, sA, 0.0, run.t.back());
    // The slow front seen in B: B falls from the plateau to the feed level.
    const double t2_meas_B = t_mid + Integral(run.t, dB, t_mid, run.t.back());
    std::cout << " t1=" << t1_meas << " (theory " << th.t1 << ") t2=" << t2_meas
              << " / " << t2_meas_B << " (theory " << th.t2 << ")\n";

    EXPECT_NEAR(t1_meas, th.t1, kTol * th.t1)
      << "Fast (weak-component) front at " << t1_meas << " s, equilibrium theory " << th.t1 << " s";
    EXPECT_NEAR(t2_meas, th.t2, kTol * th.t2)
      << "Slow (strong-component) front at " << t2_meas << " s, equilibrium theory " << th.t2 << " s";
    // Coherence: both components see the slow front at the same time.
    EXPECT_NEAR(t2_meas_B, th.t2, kTol * th.t2)
      << "The weak component's drop from the plateau (" << t2_meas_B
      << " s) does not travel with the strong component's front (" << th.t2 << " s)";
  }

} // namespace adsim
