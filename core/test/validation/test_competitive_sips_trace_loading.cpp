/*
 * Skarstrom - pressure and temperature swing adsorption simulator
 * Copyright (c) 2024-2026 Sam Affleck and contributors
 *
 * Licensed under the terms in the LICENSE file at the repository root.
 */

// Trace-component loadings of the competitive Sips kinds (arrhenius_sips,
// extended_sips_2) against the published isotherms.
//
// Both kinds share the denominator of the empirical competitive Sips form
// (Ruthven 1984; Schell et al. 2012, Eqs. 5-7; Streb & Mazzotti 2021, Eqs. 1-4):
//     q_i = q_max,i z_i / (1 + sum_j z_j),   z_i = (b_i P y_i)^{s_i}.
// For s < 1 the power has an infinite slope at y = 0, so the solver path
// (UpdateAll) regularises each z_i below a mole fraction y0 = 1e-5: there it is
// replaced by the quadratic in y that is zero at y = 0 and matches the value and
// slope of the power at y0,
//     z = (b P y0)^s t ((2 - s) - (1 - s) t),   t = y / y0
// (docs/isotherms.md, "Competitive Sips kinds near zero mole fraction"). Above y0
// the published form is used unchanged. Until WP19-1 the kinds instead
// multiplied every z by y^2 / (y^2 + 1e-6), which halved the loading of a
// component at y = 1e-3 and cut it by 1 % at y = 1e-2.
//
// A closed, uniform, isothermal bed starts at equilibrium (Equations seeds
// q = q* from the initial gas), so the loadings logged at t = 0 are the
// solver's equilibrium loadings; no time integration enters the comparison.
//
//   * y >= y0: the reference is the published closed form, written here from
//     the papers' equations (never from the solver). Tolerance 1e-6 relative:
//     the two evaluate the same algebra, so only round-off separates them.
//   * y < y0: the reference is the documented continuation above; at y = 0
//     the loading must be exactly 0 (to 1e-12 mol/kg).
//
// GUARD: at every y >= y0 tested the former step changed the loading by
// delta^2 / (y^2 + delta^2) >= 1e-4 (delta = 1e-3), 99x the tolerance, so the
// test cannot pass with the former regularisation. Measured on main before the
// fix: arrhenius_sips (Casas 2012 CO2 at 15 bar, y_H2 = 0.5) -0.96 % at
// y = 1e-2, -49.9 % at 1e-3, -99.0 % at 1e-4; extended_sips_2 (Streb 2021 CO2
// at 10 bar, y_CH4 = 0.5) -0.65 %, -45.5 %, -98.9 %. After it: below 1e-15.

#include "validation_helpers.h"

#include <cctype>
#include <iomanip>
#include <sstream>

namespace adsim {

  using namespace adsim::validation;

  namespace {

    constexpr double kY0 = 1e-5;  // documented regularisation threshold

    enum class Kind { ArrheniusSips, ExtendedSips2 };

    // Casas 2012 / Schell 2012 parameters (core/test/papers/casas_2012/casas_2012.h).
    struct ArrSips { double a, b, A, B, c_ref, c_amp, c_beta, c_Tref; };
    const ArrSips kCasasCO2{1.38286, 5625.93, 1.6814e-08, 9160.22, 0.832365, 0.116571, 0.0928714, 330.129};
    const ArrSips kCasasH2{6.66, 0.0, 6.6886e-10, 9954.73, 0.96, 0.0, 0.0, 298.15};

    // Streb 2021, Table 3 (core/test/papers/streb_2021/streb_2021.cpp).
    struct Sips2 { double A1, A2, B1, B2, C1, C2, T_ref; };
    const Sips2 kStrebCO2{7.268, -0.6168, 1.129e-9, 28389.0, 0.425, 0.7238, 298.0};
    const Sips2 kStrebCH4{4.473, -0.6569, 3.605e-9, 15561.0, 1.073, 0.0, 298.0};

    // Affinity b [1/Pa], exponent s and capacity q_max [mol/kg] at T, from the
    // papers' temperature laws.
    struct Law { double b, s, q_max; };
    Law Eval(const ArrSips& p, double T) {
      return {p.A * std::exp(p.B / (kR * T)),
              p.c_ref + p.c_amp * (2.0 / M_PI) * std::atan(p.c_beta * (T - p.c_Tref)),
              p.a * std::exp(p.b / (kR * T))};
    }
    Law Eval(const Sips2& p, double T) {
      return {p.B1 * std::exp(p.B2 / (kR * T)),
              p.C1 + p.C2 * (T / p.T_ref - 1.0),
              p.A1 * std::exp(p.A2 * (T / p.T_ref - 1.0))};
    }

    // Competitive term: the published power at y >= y0, the documented
    // continuation below.
    double Z(const Law& l, double P, double y) {
      if (y >= kY0) return std::pow(l.b * P * y, l.s);
      const double t = y / kY0;
      return std::pow(l.b * P * kY0, l.s) * t * ((2.0 - l.s) - (1.0 - l.s) * t);
    }

    struct Case {
      Kind kind;
      double y;  // mole fraction of the trace component
    };

    // The trace component and its competitor for each kind.
    struct TraceSetup {
      const char* trace;
      const char* other;
      double P, T, y_other;
      Law trace_law, other_law;
    };
    TraceSetup SetupFor(Kind k) {
      constexpr double T = 298.15;
      if (k == Kind::ArrheniusSips)  // Casas 2012: CO2 in CO2/H2 at 15 bar
        return {"CarbonDioxide", "Hydrogen", 15.0e5, T, 0.5, Eval(kCasasCO2, T), Eval(kCasasH2, T)};
      // Streb 2021: CO2 in CO2/CH4 at 10 bar (s_CO2 = 0.425, the most singular exponent in the library)
      return {"CarbonDioxide", "Methane", 10.0e5, T, 0.5, Eval(kStrebCO2, T), Eval(kStrebCH4, T)};
    }

    double Reference(const TraceSetup& s, double y) {
      const double z = Z(s.trace_law, s.P, y);
      const double z_other = Z(s.other_law, s.P, s.y_other);
      return s.trace_law.q_max * z / (1.0 + z + z_other);
    }

    void Set(ArrheniusSipsIsothermParameters& to, const ArrSips& p) {
      to.n_inf_a = p.a; to.n_inf_b = p.b; to.k_A = p.A; to.k_B = p.B;
      to.c_ref = p.c_ref; to.c_amp = p.c_amp; to.c_beta = p.c_beta; to.c_Tref = p.c_Tref;
    }
    void Set(Sips2IsothermParameters& to, const Sips2& p) {
      to.A1 = p.A1; to.A2 = p.A2; to.B1 = p.B1; to.B2 = p.B2;
      to.C1 = p.C1; to.C2 = p.C2; to.T_ref = p.T_ref;
    }

    std::string Name(const ::testing::TestParamInfo<Case>& info) {
      std::ostringstream s;
      s << (info.param.kind == Kind::ArrheniusSips ? "ArrheniusSips" : "ExtendedSips2") << "_y"
        << std::scientific << std::setprecision(0) << info.param.y;
      std::string n = s.str();
      for (char& c : n) if (!std::isalnum(static_cast<unsigned char>(c))) c = '_';
      return n;
    }

  } // namespace

  class CompetitiveSipsTraceLoading : public UserDataTest,
                                      public ::testing::WithParamInterface<Case> {
  public:
    // Seeded (t = 0) loading of the trace component, mean over the interior.
    double SeededTraceLoading(Kind kind, const TraceSetup& s, double y) {
      AddComponents({s.trace, s.other, "Helium"});
      auto [layer_id, ads_id] = AddLayer("Packing");
      auto& layer = GetLayer(layer_id);
      layer.domain_.length = 1.0;
      layer.domain_.n_cells = 5;
      layer.eb_ = 0.4;
      layer.particle_diameter_ = 0.002;
      std::get<CylinderBed>(data.bed.shape_).diameter = 0.05;
      auto& ads = GetAdsorbent(ads_id);
      ads.ep_ = 0.0;
      ads.solid_density_ = 1000.0;
      ads.heat_capacity_ = 900.0;
      ads.thermal_conductivity_ = 0.0;
      AddConstantInternalHeatTransferCoefficientModel(0.0, layer_id);
      AddConstantExternalHeatTransferCoefficientModel(0.0, layer_id);
      if (kind == Kind::ArrheniusSips) {
        auto iso = std::make_unique<ArrheniusSipsIsotherm>();
        Set(iso->params_[s.trace], kCasasCO2);
        Set(iso->params_[s.other], kCasasH2);
        ads.adsorption_source_.AddIsotherm({s.trace, s.other}, std::move(iso));
      } else {
        auto iso = std::make_unique<ExtendedSips2Isotherm>();
        Set(iso->params_[s.trace], kStrebCO2);
        Set(iso->params_[s.other], kStrebCH4);
        ads.adsorption_source_.AddIsotherm({s.trace, s.other}, std::move(iso));
      }
      for (const char* c : {s.trace, s.other}) {
        AddLinearDrivingForceModel(ads_id, c, 1.0);
        SetHOA(ads_id, c, 0.0);
      }
      data.bed.equations_.energy_balance_mode_ = EnergyBalanceMode::ISOTHERMAL;
      AddInitialConditions(s.P, s.T, s.T, {{s.trace, y}, {s.other, s.y_other},
                                           {"Helium", 1.0 - y - s.y_other}}, 0.0);
      AddStep("Hold", 1.0);
      AddClosedBC(BoundaryLocation::BOTTOM, "Hold");
      AddClosedBC(BoundaryLocation::TOP, "Hold");
      ConfigureCycle(1);
      data.settings.reltol = 1e-10;
      CoupledSolver solver{data};
      std::string error;
      EXPECT_TRUE(data.PerformChecks(error)) << error;
      EXPECT_GE(solver.IntegrateCycle(data.cycle, false), 0);
      const auto& log = solver.GetLog();
      const auto q = log.GetComponentVariable(0, 0, ComponentVariableType::ADSORBED_CONCENTRATION, s.trace);
      EXPECT_FALSE(q.values.empty());
      return q.values.empty() ? std::nan("") : InteriorMean(q.values.front(), log.GetDomain().volume);
    }
  };

  TEST_P(CompetitiveSipsTraceLoading, MatchesPublishedIsothermAboveThreshold) {
    const Case c = GetParam();
    const TraceSetup s = SetupFor(c.kind);
    const double ref = Reference(s, c.y);
    // GUARD: the former step y^2 / (y^2 + 1e-6) would move this loading by
    // at least 99x the tolerance (1e-4 at y = 0.1).
    constexpr double kTol = 1e-6;
    ASSERT_GE(c.y, kY0);
    ASSERT_GT(1e-6 / (c.y * c.y + 1e-6), 50.0 * kTol);

    const double q = SeededTraceLoading(c.kind, s, c.y);
    std::cout << "[CompetitiveSipsTrace] y = " << c.y << ": solver " << q << ", published "
              << ref << " mol/kg (" << 100.0 * (q / ref - 1.0) << " %)\n";
    EXPECT_NEAR(q, ref, kTol * ref);
  }

  INSTANTIATE_TEST_SUITE_P(TraceFractions, CompetitiveSipsTraceLoading,
      ::testing::Values(Case{Kind::ArrheniusSips, 1e-1}, Case{Kind::ArrheniusSips, 1e-2},
                        Case{Kind::ArrheniusSips, 1e-3}, Case{Kind::ArrheniusSips, 1e-4},
                        Case{Kind::ArrheniusSips, 2e-5},
                        Case{Kind::ExtendedSips2, 1e-1}, Case{Kind::ExtendedSips2, 1e-2},
                        Case{Kind::ExtendedSips2, 1e-3}, Case{Kind::ExtendedSips2, 1e-4},
                        Case{Kind::ExtendedSips2, 2e-5}),
      Name);

  // Below y0 the loading follows the documented continuation, and is exactly 0
  // with the component absent.
  class CompetitiveSipsBelowThreshold : public CompetitiveSipsTraceLoading {};

  TEST_P(CompetitiveSipsBelowThreshold, FollowsDocumentedContinuation) {
    const Case c = GetParam();
    const TraceSetup s = SetupFor(c.kind);
    const double ref = Reference(s, c.y);
    const double q = SeededTraceLoading(c.kind, s, c.y);
    std::cout << "[CompetitiveSipsTrace] y = " << c.y << ": solver " << q << ", documented "
              << ref << " mol/kg\n";
    if (c.y == 0.0) {
      EXPECT_NEAR(q, 0.0, 1e-12);
    } else {
      EXPECT_NEAR(q, ref, 1e-6 * ref);
    }
  }

  INSTANTIATE_TEST_SUITE_P(TraceFractions, CompetitiveSipsBelowThreshold,
      ::testing::Values(Case{Kind::ArrheniusSips, 5e-6}, Case{Kind::ArrheniusSips, 1e-7},
                        Case{Kind::ArrheniusSips, 0.0},
                        Case{Kind::ExtendedSips2, 5e-6}, Case{Kind::ExtendedSips2, 1e-7},
                        Case{Kind::ExtendedSips2, 0.0}),
      Name);

  // ---------------------------------------------------------------------------
  // A trace breakthrough on arrhenius_sips: the front of a component fed at
  // y_f = 1e-3 into a clean bed must travel at the speed the published isotherm
  // gives. Helium carrier, CO2 the only adsorbate, isothermal, P = 1 bar.
  //
  // Theory: for a favourable isotherm (s = 0.74 < 1 at 298 K) the front is a
  // shock, and the first moment of its outlet curve is its equilibrium-theory
  // time however much it is smeared (a mass balance over the moles stored
  // behind it):
  //     t_1 = int_0^inf (1 - y_out / y_f) dt = (L / u) (1 + (rho_b / eps) q*(y_f) / c_f),
  // with c_f = P y_f / (R T). The neglected term is the sorption effect (the
  // velocity changes by O(y_f) across the front); tolerance 2 y_f = 2e-3, as in
  // test_competitive_langmuir_fronts. The toe of the front and the bed ahead of
  // it pass through y < y0, so the run also exercises the regularised region.
  //
  // GUARD: the former step halved q*(y_f) at y_f = 1e-3, which moves t_1 by
  // about 49 %, far outside the tolerance (measured on main: 543 s against
  // 1062 s). After the fix: -1.1e-4 relative. About 1 s.
  TEST(CompetitiveSipsTraceBreakthrough, FirstMomentMatchesPublishedIsotherm) {
    constexpr double kT = 298.15, kP = 1.0e5, kEps = 0.4, kRhoS = 100.0;
    constexpr double kL = 1.0, kD = 0.1, kFlow = 0.005, kYf = 1.0e-3, kLdf = 5.0;
    constexpr double kTol = 2.0 * kYf;
    const char* kCarrier = "Helium";
    const char* kA = "CarbonDioxide";

    const Law law = Eval(kCasasCO2, kT);
    auto q_star = [&](double y, bool former_step) {
      double z = std::pow(law.b * kP * y, law.s);
      if (former_step) z *= y * y / (y * y + 1e-6);
      return law.q_max * z / (1.0 + z);
    };
    const double c_f = kP * kYf / (kR * kT);
    const double rho_b = (1.0 - kEps) * kRhoS;
    const double M = rho_b / kEps * q_star(kYf, false) / c_f;
    const double M_former = rho_b / kEps * q_star(kYf, true) / c_f;
    ASSERT_LT(law.s, 1.0) << "the isotherm must be favourable for a shock";
    ASSERT_GT(std::abs(M_former / M - 1.0), 100.0 * kTol);

    UserData data;
    data.bed.wall_.thickness_ = 0.005;
    data.fluid.AddComponent(kCarrier);
    data.fluid.AddComponent(kA);
    auto ads = data.adsorbent_directory.Add("Packing");
    auto lid = data.bed.layers_.Add("Packing");
    auto& layer = data.bed.layers_.Get(lid);
    layer.momentum_model_ = MomentumModel::KarmenCozeny;
    layer.domain_.length = kL;
    layer.domain_.n_cells = 100;
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
    auto iso = std::make_unique<ArrheniusSipsIsotherm>();
    Set(iso->params_[kA], kCasasCO2);
    adsorbent.adsorption_source_.AddIsotherm({kA}, std::move(iso));
    auto ldf = std::make_unique<LinearDrivingForce>();
    auto mtc = std::make_unique<MassTransferCoefficientConstant>();
    mtc->constant_mtc_ = kLdf;
    ldf->mtc_model_ = std::move(mtc);
    adsorbent.adsorption_source_.mass_transfer_models_[kA] = std::move(ldf);
    adsorbent.adsorption_source_.heat_of_adsorption_[kA] = 0.0;
    data.bed.equations_.energy_balance_mode_ = EnergyBalanceMode::ISOTHERMAL;

    data.ic.pressure = kP;
    data.ic.temperature = kT;
    data.ic.wall_temperature = kT;
    data.ic.velocity = 0.0;
    data.ic.mole_fraction[kCarrier] = 1.0;
    data.ic.mole_fraction[kA] = 0.0;

    // Run to about twice the predicted front time (u ~ 0.039 m/s here).
    const double u_guess = kFlow / (kP / (kR * kT) * M_PI * kD * kD / 4.0 * kEps);
    const double run_time = std::ceil(2.0 * kL / u_guess * (1.0 + M));
    data.cycle.AddStep("Feed");
    Step& step = data.cycle.GetStep("Feed");
    step.step_time_ = run_time;
    auto inlet = std::make_unique<BoundaryStream>();
    auto flow = std::make_unique<MolarFlowBoundaryValue>();
    static_cast<ConstantBoundaryValue*>(flow->profile_.get())->constant_value_ = kFlow;
    inlet->boundary_flow_ = std::move(flow);
    inlet->mole_fraction[kCarrier] = 1.0 - kYf;
    inlet->mole_fraction[kA] = kYf;
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
    outlet_bc->mole_fraction_[kCarrier] = 1.0 - kYf;
    outlet_bc->mole_fraction_[kA] = kYf;
    step.top_boundary_condition_ = std::move(outlet_bc);

    data.cycle.number_of_cycles_ = 1;
    data.settings.reltol = 1e-9;
    data.settings.max_step = 1.0;
    data.settings.max_log_points_per_step_ = static_cast<int>(run_time);  // a snapshot per second

    CoupledSolver solver{data};
    std::string error;
    ASSERT_TRUE(data.PerformChecks(error)) << error;
    ASSERT_GE(solver.IntegrateCycle(data.cycle, false), 0) << "Solver did not complete";

    const auto& log = solver.GetLog();
    const auto& dom = log.GetDomain();
    const auto P = log.GetVariable(0, 0, VariableType::PRESSURE);
    const auto& Pf = P.values.back();
    const double dP_rel = std::abs(Pf[1] - Pf[Pf.size() - 2]) / kP;
    const auto vel = log.GetVariable(0, 0, VariableType::VELOCITY);
    const auto& uf = vel.values.back();
    double u_sum = 0.0;
    int u_cnt = 0;
    for (size_t i = 1; i + 1 < uf.size(); ++i) { u_sum += uf[i]; ++u_cnt; }
    const double u = (u_sum / u_cnt) / dom.et[1];
    const auto y = log.GetComponentFaceVariable(0, 0, ComponentVariableType::MOLE_FRACTION, kA);
    ASSERT_GT(y.times.size(), 2u);
    double t1 = 0.0;
    for (size_t i = 1; i < y.times.size(); ++i) {
      const double s0 = 1.0 - y.values[i - 1].back() / kYf, s1 = 1.0 - y.values[i].back() / kYf;
      t1 += 0.5 * (s0 + s1) * (y.times[i] - y.times[i - 1]);
    }
    const double t1_theory = kL / u * (1.0 + M);

    ASSERT_LT(dP_rel, 2e-3) << "Pressure drop too large for a constant velocity";
    ASSERT_NEAR(y.values.back().back(), kYf, 1e-3 * kYf) << "CO2 has not broken through";
    std::cout.precision(8);
    std::cout << "[CompetitiveSipsTraceBreakthrough] u = " << u << " m/s, t1 = " << t1
              << " s, theory " << t1_theory << " s (" << 100.0 * (t1 / t1_theory - 1.0)
              << " %), former step would give " << kL / u * (1.0 + M_former) << " s\n";
    EXPECT_NEAR(t1, t1_theory, kTol * t1_theory);
  }

} // end adsim namespace
