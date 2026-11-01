/*
 * Skarstrom - pressure and temperature swing adsorption simulator
 * Copyright (c) 2024-2026 Sam Affleck and contributors
 *
 * Licensed under the terms in the LICENSE file at the repository root.
 */

// Roadmap item 8 — physics parity with the literature:
//
//   A. Ergun pressure drop (MomentumModel::Ergun) vs the legacy Darcy /
//      Kozeny-Carman law, checked against the closed-form steady-state Ergun
//      pressure drop of an isothermal ideal gas and its low-Reynolds limit,
//      and the Blake-Kozeny law (MomentumModel::BlakeKozeny: Darcy with the
//      constant 150, Ergun without its inertial term);
//   B. the explicit isothermal mode (EnergyBalanceMode::ISOTHERMAL), which
//      drops the temperature unknowns from the DAE, checked against the
//      "fake isothermal" (zero heat-transfer coefficient, zero heat of
//      adsorption) breakthrough it replaces;
//   C. inlet-temperature ramps (IBoundaryStream / PressureBoundaryCondition
//      temperature_profile_), unit and integration level;
//   D. serialisation round trips and legacy-archive defaults for all three.
//
// The analytic references are written from first principles (validation
// helpers rule), not from the solver's discretisation.

#include "breakthrough_case.h"

#include <chrono>
#include <functional>
#include <filesystem>
#include <fstream>
#include <sstream>

namespace adsim {

  using namespace adsim::validation;

  namespace {

    // ---- A. steady single-layer flow deck -----------------------------------
    constexpr double kErgunLength = 1.0;      // [m]
    constexpr double kErgunDiameter = 0.1;    // [m]
    constexpr double kErgunVoidage = 0.4;     // [-]
    constexpr double kErgunDp = 2.0e-3;       // [m]
    constexpr double kErgunPout = 1.0e5;      // [Pa]
    constexpr double kErgunT = 298.15;        // [K]

    struct SteadyFlowResult {
      double dP = 0.0;   // inlet-face minus outlet-face pressure [Pa]
      double mu = 0.0;   // logged gas viscosity [Pa s]
      double M = 0.0;    // molar mass [kg/mol]
      double u_mid = 0.0; // superficial velocity at the middle face [m/s]
    };

    void BuildSteadyFlowDeck(UserData& data, double molar_flow, MomentumModel model, int n_cells) {
      data.fluid.AddComponent("Nitrogen");

      auto ads = data.adsorbent_directory.Add("Inert packing");
      auto lid = data.bed.layers_.Add("Packing");
      auto& layer = data.bed.layers_.Get(lid);
      layer.domain_.length = kErgunLength;
      layer.domain_.n_cells = n_cells;
      layer.eb_ = kErgunVoidage;
      layer.particle_diameter_ = kErgunDp;
      layer.composition_[ads] = 1.0;
      layer.momentum_model_ = model;
      auto disp = std::make_unique<DispersionModelConstant>();
      disp->dispersion_coefficient_ = 0.0;
      layer.equations_.dispersion_model_ = std::move(disp);
      std::get<CylinderBed>(data.bed.shape_).diameter = kErgunDiameter;

      auto& adsorbent = data.adsorbent_directory.Get(ads);
      adsorbent.ep_ = 0.0;

      // Isothermal: mu and T are uniform, so the steady Ergun law integrates in
      // closed form (below).
      data.bed.equations_.energy_balance_mode_ = EnergyBalanceMode::ISOTHERMAL;

      data.ic.pressure = kErgunPout;
      data.ic.temperature = kErgunT;
      data.ic.wall_temperature = kErgunT;
      data.ic.velocity = 0.0;
      data.ic.mole_fraction["Nitrogen"] = 1.0;

      data.cycle.AddStep("Flow");
      Step& step = data.cycle.GetStep("Flow");
      step.step_time_ = 20.0;

      auto inlet = std::make_unique<BoundaryStream>();
      inlet->name_ = "Feed";
      auto flow = std::make_unique<MolarFlowBoundaryValue>();
      static_cast<ConstantBoundaryValue*>(flow->profile_.get())->constant_value_ = molar_flow;
      inlet->boundary_flow_ = std::move(flow);
      inlet->mole_fraction["Nitrogen"] = 1.0;
      inlet->temperature_ = kErgunT;
      auto inlet_bc = std::make_unique<FlowBoundaryCondition>();
      inlet_bc->direction_ = BoundaryDirection::INLET;
      inlet_bc->location_ = BoundaryLocation::BOTTOM;
      inlet_bc->streams_.push_back(std::move(inlet));
      step.bottom_boundary_condition_ = std::move(inlet_bc);

      auto outlet_bc = std::make_unique<PressureBoundaryCondition>();
      outlet_bc->direction_ = BoundaryDirection::OUTLET;
      outlet_bc->location_ = BoundaryLocation::TOP;
      outlet_bc->pressure_ = std::make_unique<ConstantBoundaryValue>(kErgunPout);
      outlet_bc->mole_fraction_["Nitrogen"] = 1.0;
      step.top_boundary_condition_ = std::move(outlet_bc);

      data.cycle.number_of_cycles_ = 1;
      data.settings.reltol = 1e-10;
      data.settings.max_step = 1.0;
    }

    SteadyFlowResult RunSteadyFlow(double molar_flow, MomentumModel model, int n_cells = 100) {
      UserData data;
      BuildSteadyFlowDeck(data, molar_flow, model, n_cells);
      std::string err;
      EXPECT_TRUE(data.PerformChecks(err)) << err;
      CoupledSolver solver{data};
      EXPECT_GE(solver.IntegrateCycle(data.cycle, false), 0) << "Solver did not complete";

      const auto& log = solver.GetLog();
      const auto P = log.GetVariable(0, 0, VariableType::PRESSURE).values.back();
      const auto mu = log.GetVariable(0, 0, VariableType::VISCOSITY).values.back();
      const auto u = log.GetVariable(0, 0, VariableType::VELOCITY).values.back();
      const size_t N = P.size();
      SteadyFlowResult r;
      // Face pressures at z = 0 and z = L: the ghost / first-cell and
      // last-cell / ghost averages (the faces the boundary conditions act on).
      const double P_in = 0.5 * (P[0] + P[1]);
      const double P_out = 0.5 * (P[N - 2] + P[N - 1]);
      r.dP = P_in - P_out;
      r.mu = mu[N / 2];
      r.M = data.fluid.molar_mass_.at("Nitrogen");
      r.u_mid = u[u.size() / 2];
      EXPECT_NEAR(P_out, kErgunPout, 1e-6 * kErgunPout);
      return r;
    }

    // Closed-form steady pressure drop of an isothermal ideal gas through a
    // packed bed obeying  -dP/dz = a mu u + b rho u^2  with the mass flux
    // G = rho u constant and rho = P M / (R T):
    //   -P dP/dz = (R T / M) (a mu G + b G^2)
    //   =>  P_in^2 - P_out^2 = 2 L (R T / M) (a mu G + b G^2).
    double AnalyticPressureDrop(double a, double b, double mu, double M, double molar_flow) {
      const double area = 0.25 * M_PI * kErgunDiameter * kErgunDiameter;
      const double G = molar_flow * M / area;  // [kg/(m2 s)]
      const double rhs = 2.0 * kErgunLength * (kR * kErgunT / M) * (a * mu * G + b * G * G);
      return std::sqrt(kErgunPout * kErgunPout + rhs) - kErgunPout;
    }

    // Ergun (1952) coefficients, written out independently of the solver.
    double ErgunA(double constant) {
      const double e = kErgunVoidage;
      return constant * (1 - e) * (1 - e) / (e * e * e * kErgunDp * kErgunDp);
    }
    double ErgunB() {
      const double e = kErgunVoidage;
      return 1.75 * (1 - e) / (e * e * e * kErgunDp);
    }

    double Seconds(std::chrono::steady_clock::time_point a, std::chrono::steady_clock::time_point b) {
      return std::chrono::duration<double>(b - a).count();
    }

    BoundaryStream* FirstStream(Step& step, bool top) {
      auto* bc = dynamic_cast<FlowBoundaryCondition*>(
        top ? step.top_boundary_condition_.get() : step.bottom_boundary_condition_.get());
      if (!bc || bc->streams_.empty()) return nullptr;
      return dynamic_cast<BoundaryStream*>(bc->streams_.front().get());
    }

  } // namespace

  // ===========================================================================
  // A. Ergun pressure drop
  // ===========================================================================

  TEST(ErgunVelocity, ClosedFormRootIsExactAndSmooth) {
    const double A = 3.7e3, B = 8.2e3;
    for (double g : {-5.0e4, -12.0, -1e-9, 0.0, 1e-9, 3.0, 750.0, 2.0e5}) {
      const double u = ErgunVelocity(g, A, B);
      // Root of g = A u + B |u| u, with the sign of g.
      EXPECT_NEAR(A * u + B * std::abs(u) * u, g, 1e-12 * std::max(1.0, std::abs(g))) << g;
      EXPECT_EQ(ErgunVelocity(-g, A, B), -u) << g;  // odd
    }
    // B = 0 is Darcy exactly (no rounding difference).
    for (double g : {-4.0, 0.0, 2.5, 1.0e6}) EXPECT_EQ(ErgunVelocity(g, A, 0.0), g / A);
    // Smooth through g = 0: slope 1/A from both sides.
    const double h = 1e-8;
    EXPECT_NEAR(ErgunVelocity(h, A, B) / h, 1.0 / A, 1e-9 / A);
    EXPECT_NEAR(ErgunVelocity(-h, A, B) / -h, 1.0 / A, 1e-9 / A);
    // Degenerate input does not produce NaN.
    EXPECT_EQ(ErgunVelocity(1.0, 0.0, 0.0), 0.0);
    // The implicit residual the solver uses vanishes exactly at that root, and
    // is the Darcy residual when B = 0.
    for (double g : {-5.0e4, -12.0, 0.0, 3.0, 2.0e5}) {
      const double u = ErgunVelocity(g, A, B);
      EXPECT_NEAR(ErgunMomentumResidual(u, g, A, B), 0.0, 1e-12 * std::max(1.0, std::abs(u))) << g;
      EXPECT_EQ(ErgunMomentumResidual(0.3, g, A, 0.0), 0.3 - g / A);
    }
  }

  TEST(ErgunPressureDrop, LayerCoefficientsMatchTheLiterature) {
    const double e = 0.37, dp = 1.5e-3;
    EXPECT_DOUBLE_EQ(GetLayerPermeability(e, dp, MomentumModel::KarmenCozeny),
                     e * e * e * dp * dp / (180 * (1 - e) * (1 - e)));
    EXPECT_DOUBLE_EQ(GetLayerPermeability(e, dp, MomentumModel::Ergun),
                     e * e * e * dp * dp / (150 * (1 - e) * (1 - e)));
    EXPECT_EQ(GetLayerInertialCoefficient(e, dp, MomentumModel::KarmenCozeny), 0.0);
    EXPECT_DOUBLE_EQ(GetLayerInertialCoefficient(e, dp, MomentumModel::Ergun),
                     1.75 * (1 - e) / (e * e * e * dp));
    EXPECT_DOUBLE_EQ(GetLayerPermeability(e, dp, MomentumModel::BlakeKozeny),
                     e * e * e * dp * dp / (150 * (1 - e) * (1 - e)));
    EXPECT_EQ(GetLayerInertialCoefficient(e, dp, MomentumModel::BlakeKozeny), 0.0);
  }

  // Blake-Kozeny as Estupinan Perez et al. (2019, SI Eq. S5) write it, in the
  // interstitial velocity v and particle radius r_p:
  //   -dP/dz = (150/4) (1/r_p^2) ((1 - e)/e)^2 mu v.
  // With the engine's superficial velocity u = e v and d_p = 2 r_p this is
  // Darcy's law -dP/dz = (mu/K) u with the layer's Blake-Kozeny permeability.
  TEST(BlakeKozenyPressureDrop, PermeabilityIsThePapersInterstitialLaw) {
    const double mu = 1.8e-5, u = 0.07;  // [Pa s], superficial [m/s]
    for (double e : {0.3, 0.37, 0.45}) {
      for (double dp : {1.0e-3, 2.0e-3}) {
        const double rp = 0.5 * dp, v = u / e;
        const double grad_paper = (150.0 / 4.0) / (rp * rp) * ((1 - e) / e) * ((1 - e) / e) * mu * v;
        const double K = GetLayerPermeability(e, dp, MomentumModel::BlakeKozeny);
        EXPECT_NEAR(mu * u / K, grad_paper, 1e-12 * grad_paper) << "e " << e << ", dp " << dp;
        // Ergun's viscous term exactly, and 150/180 of the Kozeny-Carman one.
        EXPECT_DOUBLE_EQ(K, GetLayerPermeability(e, dp, MomentumModel::Ergun));
        EXPECT_NEAR(GetLayerPermeability(e, dp, MomentumModel::KarmenCozeny) / K, 150.0 / 180.0, 1e-14);
      }
    }
  }

  // Steady flow at a Reynolds number where the inertial term is about as large
  // as the viscous one: the simulated bed pressure drop must match the
  // closed-form Ergun result.
  TEST(ErgunPressureDrop, SteadyFlowMatchesAnalyticErgun) {
    const double molar_flow = 0.128;  // [mol/s] -> G ~ 0.46 kg/(m2 s), u ~ 0.4 m/s
    const auto r = RunSteadyFlow(molar_flow, MomentumModel::Ergun);
    const double dP_ergun = AnalyticPressureDrop(ErgunA(150.0), ErgunB(), r.mu, r.M, molar_flow);
    const double dP_darcy150 = AnalyticPressureDrop(ErgunA(150.0), 0.0, r.mu, r.M, molar_flow);

    // GUARD: the inertial term must matter here, or this is just Darcy again.
    ASSERT_GT(dP_ergun, 1.5 * dP_darcy150)
      << "Inertial term too small to test Ergun (dP_ergun " << dP_ergun
      << " Pa vs viscous-only " << dP_darcy150 << " Pa)";
    // Particle Reynolds number rho u dp / mu, for the record.
    const double area = 0.25 * M_PI * kErgunDiameter * kErgunDiameter;
    const double Re_p = (molar_flow * r.M / area) * kErgunDp / r.mu;
    RecordProperty("Re_p", std::to_string(Re_p));
    RecordProperty("dP_sim_Pa", std::to_string(r.dP));
    RecordProperty("dP_analytic_Pa", std::to_string(dP_ergun));

    EXPECT_NEAR(r.dP, dP_ergun, 1e-5 * dP_ergun)
      << "Re_p = " << Re_p << "; simulated dP " << r.dP << " Pa, analytic Ergun " << dP_ergun << " Pa";
  }

  // At low Reynolds number the inertial term vanishes: Ergun -> Darcy with the
  // constant 150, and the legacy Kozeny-Carman law (constant 180) gives exactly
  // 180/150 = 1.2 times that pressure drop.
  TEST(ErgunPressureDrop, ReducesToDarcy150AtLowReynolds) {
    const double molar_flow = 1.28e-4;  // [mol/s], 1000x below the case above
    const auto ergun = RunSteadyFlow(molar_flow, MomentumModel::Ergun);
    const auto kc = RunSteadyFlow(molar_flow, MomentumModel::KarmenCozeny);
    const double dP_darcy150 = AnalyticPressureDrop(ErgunA(150.0), 0.0, ergun.mu, ergun.M, molar_flow);
    const double dP_darcy180 = AnalyticPressureDrop(ErgunA(180.0), 0.0, kc.mu, kc.M, molar_flow);

    ASSERT_GT(ergun.dP, 0.0);
    EXPECT_NEAR(ergun.dP, dP_darcy150, 3e-3 * dP_darcy150);
    EXPECT_NEAR(kc.dP, dP_darcy180, 2e-3 * dP_darcy180);
    EXPECT_NEAR(kc.dP / ergun.dP, 1.2, 4e-3);
  }

  // Blake-Kozeny at the flow rate where Ergun's inertial term doubles the
  // pressure drop: no inertial term, so the simulated drop is the closed-form
  // Darcy result with the constant 150, well below Ergun's. For an isothermal
  // ideal gas P_in^2 - P_out^2 is proportional to the viscous constant, so the
  // Kozeny-Carman bed's P_in^2 - P_out^2 is exactly 180/150 times this one.
  TEST(BlakeKozenyPressureDrop, SteadyFlowMatchesAnalyticDarcy150) {
    const double molar_flow = 0.128;  // [mol/s], as in SteadyFlowMatchesAnalyticErgun
    const auto bk = RunSteadyFlow(molar_flow, MomentumModel::BlakeKozeny);
    const auto kc = RunSteadyFlow(molar_flow, MomentumModel::KarmenCozeny);
    const double dP_darcy150 = AnalyticPressureDrop(ErgunA(150.0), 0.0, bk.mu, bk.M, molar_flow);
    const double dP_ergun = AnalyticPressureDrop(ErgunA(150.0), ErgunB(), bk.mu, bk.M, molar_flow);

    ASSERT_GT(bk.dP, 0.0);
    EXPECT_NEAR(bk.dP, dP_darcy150, 1e-5 * dP_darcy150)
      << "simulated dP " << bk.dP << " Pa, analytic Blake-Kozeny " << dP_darcy150 << " Pa";
    EXPECT_LT(bk.dP, 0.75 * dP_ergun) << "Blake-Kozeny must not carry Ergun's inertial term";

    // Superficial velocity: u = F R T / (P A) at the middle face.
    const double area = 0.25 * M_PI * kErgunDiameter * kErgunDiameter;
    EXPECT_GT(bk.u_mid, 0.0);
    EXPECT_LT(bk.u_mid, molar_flow * kR * kErgunT / (kErgunPout * area));

    const auto sq = [](const SteadyFlowResult& r) {
      const double P_in = kErgunPout + r.dP;
      return P_in * P_in - kErgunPout * kErgunPout;
    };
    ASSERT_NEAR(kc.mu, bk.mu, 1e-12 * bk.mu);
    EXPECT_NEAR(sq(kc) / sq(bk), 180.0 / 150.0, 2e-5);
  }

  // ===========================================================================
  // B. Isothermal mode
  // ===========================================================================

  namespace {
    struct BreakthroughRun {
      std::vector<double> t, y;
      double seconds = 0.0;
      int system_size = 0;
    };

    BreakthroughRun RunBreakthrough(UserData& data,
                                    const std::function<void(const EquationLog&)>& inspect = {}) {
      std::string err;
      EXPECT_TRUE(data.PerformChecks(err)) << err;
      CoupledSolver solver{data};
      const auto t0 = std::chrono::steady_clock::now();
      EXPECT_GE(solver.IntegrateCycle(data.cycle, false), 0) << "Solver did not complete";
      const auto t1 = std::chrono::steady_clock::now();
      BreakthroughRun r;
      std::tie(r.t, r.y) = OutletSeries(solver.GetLog(), false);
      r.seconds = Seconds(t0, t1);
      r.system_size = data.state.GetSystemSize();
      if (inspect) inspect(solver.GetLog());
      return r;
    }

    double MaxAbsDifference(const BreakthroughRun& a, const BreakthroughRun& b, double t_end) {
      double worst = 0.0;
      for (int k = 0; k <= 400; ++k) {
        const double t = t_end * k / 400.0;
        worst = std::max(worst, std::abs(interp_linear(a.t, a.y, t) - interp_linear(b.t, b.y, t)));
      }
      return worst;
    }
  } // namespace

  // The isothermal switch reproduces the "fake isothermal" breakthrough (zero
  // HTCs and zero heat of adsorption, which leaves only the small compression /
  // adsorbed-phase heat effects) with three fewer unknowns per cell.
  TEST(IsothermalMode, MatchesZeroHeatTransferBreakthrough) {
    BreakthroughCase c;
    c.adsorbing = true;
    c.n_cells = 60;
    c.run_time = 400.0;
    c.reltol = 1e-8;

    UserData fake;
    Build(fake, c);
    const auto ref = RunBreakthrough(fake);

    UserData iso;
    Build(iso, c);
    iso.bed.equations_.energy_balance_mode_ = EnergyBalanceMode::ISOTHERMAL;
    // Every temperature is pinned at the initial-condition temperature, and the
    // (logged) energy sources are inactive.
    size_t pinned = 0;
    bool sources_zero = true;
    const auto run = RunBreakthrough(iso, [&](const EquationLog& log) {
      for (auto type : {VariableType::TEMPERATURE, VariableType::WALL_TEMPERATURE, VariableType::SOLID_TEMPERATURE}) {
        const auto v = log.GetVariable(0, 0, type);
        for (const auto& snap : v.values)
          for (size_t i = 1; i + 1 < snap.size(); ++i) {
            EXPECT_EQ(snap[i], c.bed_temperature);
            ++pinned;
          }
      }
      for (const auto& snap : log.GetVariable(0, 0, VariableType::ENERGY_SOURCE).values)
        for (double q : snap) sources_zero = sources_zero && (q == 0.0);
    });
    EXPECT_GT(pinned, 0u);
    EXPECT_TRUE(sources_zero);

    // The temperature unknowns are gone: fluid + wall at every cell incl.
    // ghosts, and the solid at every packed cell incl. ghosts unless the LTE
    // reference already aliases it to the fluid temperature
    // (IDASettings::alias_lte_solid_temperature_), in which case it was never
    // an unknown there either.
    const int n_total = c.n_cells + 2;
    const int removed_per_cell = fake.state.SolidTemperatureAliased() ? 2 : 3;
    EXPECT_EQ(ref.system_size - run.system_size, removed_per_cell * n_total);

    // The breakthrough curves agree closely.
    const double diff = MaxAbsDifference(ref, run, c.run_time);
    RecordProperty("max_abs_dy_over_yfeed", std::to_string(diff / c.y_feed));
    RecordProperty("seconds_fake_isothermal", std::to_string(ref.seconds));
    RecordProperty("seconds_isothermal", std::to_string(run.seconds));
    std::cout << "[IsothermalMode] max |dy|/y_feed = " << diff / c.y_feed
              << "; run time fake-isothermal " << ref.seconds << " s, isothermal " << run.seconds
              << " s; unknowns " << ref.system_size << " -> " << run.system_size << "\n";
    EXPECT_LT(diff, 5e-3 * c.y_feed);
  }

  // A user-set T_iso pins every logged temperature, and the inlet temperature,
  // the heat of adsorption and an electric energy source have no effect.
  TEST(IsothermalMode, HoldsUserTemperatureAndIgnoresHeatInputs) {
    auto make = [](double feed_T, bool heater, double hoa) {
      BreakthroughCase c;
      c.adsorbing = true;
      c.n_cells = 30;
      c.run_time = 150.0;
      c.reltol = 1e-8;
      c.feed_temperature = feed_T;
      auto data = std::make_unique<UserData>();
      Build(*data, c);
      auto& eqs = data->bed.equations_;
      eqs.energy_balance_mode_ = EnergyBalanceMode::ISOTHERMAL;
      eqs.isothermal_use_initial_temperature_ = false;
      eqs.isothermal_temperature_ = 320.0;
      auto& ads = data->adsorbent_directory.Get(data->bed.layers_.layers_.front().composition_.begin()->first);
      ads.adsorption_source_.heat_of_adsorption_["Oxygen"] = hoa;
      if (heater) {
        auto& src = data->cycle.GetStep("Feed").energy_sources_[data->bed.layers_.layers_.front().id_];
        src.enabled_ = true;
        src.power_ = 5000.0;
      }
      return data;
    };

    auto base = make(298.0, false, 0.0);
    CoupledSolver s_base{*base};
    ASSERT_GE(s_base.IntegrateCycle(base->cycle, false), 0);
    const auto [t0, y0] = OutletSeries(s_base.GetLog(), false);
    for (auto type : {VariableType::TEMPERATURE, VariableType::WALL_TEMPERATURE, VariableType::SOLID_TEMPERATURE})
      for (const auto& snap : s_base.GetLog().GetVariable(0, 0, type).values)
        for (size_t i = 1; i + 1 < snap.size(); ++i) ASSERT_EQ(snap[i], 320.0);
    // Stream-side temperatures are T_iso too.
    auto* feed = FirstStream(base->cycle.GetStep("Feed"), false);
    ASSERT_NE(feed, nullptr);
    for (double T : feed->stream_data_.temperature_boundary_) ASSERT_EQ(T, 320.0);

    auto hot = make(450.0, true, 40000.0);
    CoupledSolver s_hot{*hot};
    ASSERT_GE(s_hot.IntegrateCycle(hot->cycle, false), 0);
    const auto [t1, y1] = OutletSeries(s_hot.GetLog(), false);
    ASSERT_EQ(t0.size(), t1.size());
    for (size_t k = 0; k < y0.size(); ++k) {
      ASSERT_EQ(t0[k], t1[k]);
      ASSERT_EQ(y0[k], y1[k]);
    }
  }

  // Compressor power is still computed in an isothermal run, with the gas at
  // T_iso on the suction side.
  TEST(IsothermalMode, PumpPowerUsesIsothermalTemperature) {
    BreakthroughCase c;
    c.n_cells = 20;
    c.run_time = 20.0;
    c.reltol = 1e-8;
    c.feed_temperature = 400.0;  // must not matter
    UserData data;
    Build(data, c);
    data.bed.equations_.energy_balance_mode_ = EnergyBalanceMode::ISOTHERMAL;
    auto* feed = FirstStream(data.cycle.GetStep("Feed"), false);
    ASSERT_NE(feed, nullptr);
    feed->pump_->use_pump_ = true;
    std::string err;
    ASSERT_TRUE(data.PerformChecks(err)) << err;
    CoupledSolver solver{data};
    ASSERT_GE(solver.IntegrateCycle(data.cycle, false), 0);

    const auto& sd = feed->stream_data_;
    ASSERT_FALSE(sd.electrical_power_.empty());
    const size_t k = sd.electrical_power_.size() - 1;
    const double expected = feed->pump_->CalculatePower(sd.pressure_boundary_[k], sd.molar_flow_[k], c.bed_temperature);
    ASSERT_TRUE(std::isfinite(expected));
    EXPECT_GT(sd.electrical_power_[k], 0.0);
    EXPECT_DOUBLE_EQ(sd.electrical_power_[k], expected);
    EXPECT_EQ(sd.temperature_boundary_[k], c.bed_temperature);  // no compression heating
  }

  // ===========================================================================
  // C. Inlet temperature ramps
  // ===========================================================================

  TEST(InletTemperatureProfile, FollowsLinearAndExponentialRampsWithoutExtrapolating) {
    BoundaryStream s;
    s.SetDirection(BoundaryDirection::INLET);
    s.temperature_ = 300.0;
    // Constant (nullptr) -> temperature_.
    EXPECT_EQ(s.InletTemperature(0.0), 300.0);
    EXPECT_EQ(s.InletTemperature(1e4), 300.0);

    s.temperature_profile_ = CreateInletTemperatureProfile(BoundaryValueType::LINEAR, 300.0, 100.0);
    // Switching the type changes nothing until the values are edited.
    EXPECT_EQ(s.InletTemperature(50.0), 300.0);
    static_cast<LinearBoundaryValue*>(s.temperature_profile_.get())->end_value_ = 400.0;
    EXPECT_DOUBLE_EQ(s.InletTemperature(0.0), 300.0);
    EXPECT_DOUBLE_EQ(s.InletTemperature(25.0), 325.0);
    EXPECT_DOUBLE_EQ(s.InletTemperature(50.0), 350.0);
    EXPECT_DOUBLE_EQ(s.InletTemperature(100.0), 400.0);
    // Held after the ramp, clamped before it (no extrapolation; cf. bug E-05).
    EXPECT_DOUBLE_EQ(s.InletTemperature(250.0), 400.0);
    EXPECT_DOUBLE_EQ(s.InletTemperature(1e6), 400.0);
    EXPECT_DOUBLE_EQ(s.InletTemperature(-5.0), 300.0);
    // The temperature_ leaf is not what a ramp uses.
    s.temperature_ = 123.0;
    EXPECT_DOUBLE_EQ(s.InletTemperature(50.0), 350.0);

    auto exp_profile = CreateInletTemperatureProfile(BoundaryValueType::EXPONENTIAL, 300.0, 60.0);
    static_cast<ExponentialBoundaryValue*>(exp_profile.get())->end_value_ = 400.0;
    s.temperature_profile_ = std::move(exp_profile);
    EXPECT_DOUBLE_EQ(s.InletTemperature(0.0), 300.0);
    EXPECT_NEAR(s.InletTemperature(60.0), 399.0, 1e-9);           // 99 % at the ramp time
    EXPECT_NEAR(s.InletTemperature(30.0), 300.0 + 90.0, 1e-9);    // 1 - 100^-0.5 = 0.9
    EXPECT_LT(s.InletTemperature(1e3), 400.0 + 1e-9);             // saturates, never overshoots
    EXPECT_GT(s.InletTemperature(1e3), 399.999);

    // A ConstantBoundaryValue profile means "constant temperature_".
    s.temperature_ = 310.0;
    s.temperature_profile_ = std::make_unique<ConstantBoundaryValue>(999.0);
    EXPECT_EQ(s.InletTemperature(10.0), 310.0);
    EXPECT_EQ(CreateInletTemperatureProfile(BoundaryValueType::CONSTANT, 310.0, 5.0), nullptr);
  }

  TEST(InletTemperatureProfile, BufferTankAndIsothermalKeepPrecedence) {
    BoundaryStream s;
    s.SetDirection(BoundaryDirection::INLET);
    s.temperature_ = 300.0;
    s.temperature_profile_ = CreateInletTemperatureProfile(BoundaryValueType::LINEAR, 300.0, 100.0);
    static_cast<LinearBoundaryValue*>(s.temperature_profile_.get())->end_value_ = 400.0;

    CycleEndpointSummary summary;
    summary.is_valid = false;
    summary.avg_temperature = 280.0;
    s.buffer_tank_.Activate("Src", "Stream 1");
    s.buffer_tank_.summary_ = &summary;
    // Tank not ready: the profile is the first-cycle assumption.
    EXPECT_DOUBLE_EQ(s.InletTemperature(50.0), 350.0);
    // Tank ready: its temperature wins over the profile, as before.
    summary.is_valid = true;
    EXPECT_DOUBLE_EQ(s.InletTemperature(50.0), 280.0);
    // Isothermal run: every stream-side gas is at T_iso.
    s.isothermal_temperature_ = 305.0;
    EXPECT_DOUBLE_EQ(s.InletTemperature(50.0), 305.0);
    EXPECT_DOUBLE_EQ(s.AssumedInletTemperature(50.0), 305.0);

    // Same rules on a pressure-boundary inlet.
    PressureBoundaryCondition pbc;
    pbc.direction_ = BoundaryDirection::INLET;
    pbc.temperature_ = 300.0;
    pbc.temperature_profile_ = CreateInletTemperatureProfile(BoundaryValueType::LINEAR, 300.0, 10.0);
    static_cast<LinearBoundaryValue*>(pbc.temperature_profile_.get())->end_value_ = 350.0;
    EXPECT_DOUBLE_EQ(pbc.InletTemperature(5.0), 325.0);
    EXPECT_DOUBLE_EQ(pbc.InletTemperature(50.0), 350.0);
    pbc.buffer_tank_.Activate("Src", "Stream 1");
    pbc.buffer_tank_.summary_ = &summary;
    EXPECT_DOUBLE_EQ(pbc.InletTemperature(5.0), 280.0);
  }

  // In a run, the logged inlet stream temperature follows the ramp at several
  // times and holds the final value after it; the bed responds.
  TEST(InletTemperatureProfile, RunFollowsLinearRamp) {
    BreakthroughCase c;
    c.n_cells = 30;
    c.run_time = 100.0;
    c.reltol = 1e-8;
    UserData data;
    Build(data, c);
    auto* feed = FirstStream(data.cycle.GetStep("Feed"), false);
    ASSERT_NE(feed, nullptr);
    feed->temperature_profile_ = CreateInletTemperatureProfile(BoundaryValueType::LINEAR, 298.0, 40.0);
    static_cast<LinearBoundaryValue*>(feed->temperature_profile_.get())->end_value_ = 348.0;

    std::string err;
    ASSERT_TRUE(data.PerformChecks(err)) << err;
    CoupledSolver solver{data};
    ASSERT_GE(solver.IntegrateCycle(data.cycle, false), 0);

    const auto& sd = feed->stream_data_;
    ASSERT_GT(sd.time_series_.size(), 10u);
    // The stream log holds one point per accepted step and is read back by
    // linear interpolation, so the ramp's kink at t = 40 is rounded off by up
    // to the step length there (0.1 K with the error test's default settings);
    // sample the ramp and the plateau, not the kink itself.
    for (double t : {0.0, 10.0, 20.0, 30.0, 39.0, 60.0, 100.0}) {
      const double expected = 298.0 + 50.0 * std::min(t, 40.0) / 40.0;
      EXPECT_NEAR(sd.GetTemperature(t), expected, 0.05) << "t = " << t;
    }
    for (size_t k = 0; k < sd.time_series_.size(); ++k)
      EXPECT_LE(sd.temperature_boundary_[k], 348.0 + 1e-9) << "extrapolated past the ramp end";

    // The hot feed starts heating the first cell (adiabatic bed, large solid
    // heat capacity: a few kelvin within the run).
    const auto T = solver.GetLog().GetVariable(0, 0, VariableType::TEMPERATURE).values.back();
    EXPECT_GT(T[1], 298.5);
  }

  // ===========================================================================
  // D. Serialisation
  // ===========================================================================

  TEST(PhysicsParitySerialisation, RoundTripsAllThreeAndRecordsThemExplicitly) {
    BreakthroughCase c;
    UserData data;
    Build(data, c);
    // A layer added to the deck defaults to Ergun.
    const LayerID fresh = data.bed.layers_.Add("Fresh");
    EXPECT_EQ(data.bed.layers_.Get(fresh).momentum_model_, MomentumModel::Ergun);
    data.bed.layers_.Remove(fresh);
    // Non-default setting of each feature (BreakthroughCase pins Kozeny-Carman).
    data.bed.layers_.layers_.front().momentum_model_ = MomentumModel::Ergun;
    data.bed.equations_.energy_balance_mode_ = EnergyBalanceMode::ISOTHERMAL;
    data.bed.equations_.isothermal_use_initial_temperature_ = false;
    data.bed.equations_.isothermal_temperature_.SetExpr("310 [K]");
    auto* feed = FirstStream(data.cycle.GetStep("Feed"), false);
    ASSERT_NE(feed, nullptr);
    feed->temperature_profile_ = CreateInletTemperatureProfile(BoundaryValueType::EXPONENTIAL, 300.0, 42.0);
    static_cast<ExponentialBoundaryValue*>(feed->temperature_profile_.get())->end_value_ = 390.0;
    auto* outlet = dynamic_cast<PressureBoundaryCondition*>(data.cycle.GetStep("Feed").top_boundary_condition_.get());
    ASSERT_NE(outlet, nullptr);
    outlet->temperature_profile_ = CreateInletTemperatureProfile(BoundaryValueType::LINEAR, 300.0, 7.0);

    std::string json;
    SaveToString(data, json);
    EXPECT_NE(json.find("\"Momentum_Model\": 1"), std::string::npos);
    EXPECT_NE(json.find("\"Energy_Balance_Mode\": 1"), std::string::npos);
    EXPECT_NE(json.find("\"Temperature_Profile\""), std::string::npos);
    EXPECT_NE(json.find("\"Inlet_Temperature_Profile\""), std::string::npos);
    EXPECT_EQ(json.find("\"Permeability_Model\""), std::string::npos);

    UserData back;
    LoadFromString(back, json);
    EXPECT_EQ(back.bed.layers_.layers_.front().momentum_model_, MomentumModel::Ergun);
    EXPECT_EQ(back.bed.equations_.energy_balance_mode_, EnergyBalanceMode::ISOTHERMAL);
    EXPECT_FALSE(back.bed.equations_.isothermal_use_initial_temperature_);
    EXPECT_EQ(back.bed.equations_.isothermal_temperature_.GetExpr(), "310 [K]");
    auto* feed_back = FirstStream(back.cycle.GetStep("Feed"), false);
    ASSERT_NE(feed_back, nullptr);
    auto* e = dynamic_cast<ExponentialBoundaryValue*>(feed_back->temperature_profile_.get());
    ASSERT_NE(e, nullptr);
    EXPECT_EQ(e->end_value_.Get(), 390.0);
    EXPECT_EQ(e->step_time_.Get(), 42.0);
    auto* outlet_back = dynamic_cast<PressureBoundaryCondition*>(back.cycle.GetStep("Feed").top_boundary_condition_.get());
    ASSERT_NE(outlet_back, nullptr);
    EXPECT_NE(dynamic_cast<LinearBoundaryValue*>(outlet_back->temperature_profile_.get()), nullptr);

    // A Kozeny-Carman layer is recorded explicitly too (value 0).
    back.bed.layers_.layers_.front().momentum_model_ = MomentumModel::KarmenCozeny;
    std::string json2;
    SaveToString(back, json2);
    EXPECT_NE(json2.find("\"Momentum_Model\": 0"), std::string::npos);
    UserData back2;
    LoadFromString(back2, json2);
    EXPECT_EQ(back2.bed.layers_.layers_.front().momentum_model_, MomentumModel::KarmenCozeny);
  }

  // A deck saved before this change (a verbatim copy of an example deck from
  // that version) loads as Darcy / Kozeny-Carman, non-isothermal, with constant
  // inlet temperatures — even though a default-constructed layer is Ergun.
  TEST(PhysicsParitySerialisation, LegacyArchiveLoadsWithLegacyBehaviour) {
    const std::filesystem::path deck = "test/fixtures/legacy/Ward_2022_Q=50[ccm]_pre_physics_parity.json";
    std::ifstream in(deck);
    ASSERT_TRUE(in.good()) << "missing fixture " << deck;
    std::stringstream ss;
    ss << in.rdbuf();
    const std::string legacy = ss.str();
    ASSERT_NE(legacy.find("\"Permeability_Model\""), std::string::npos);
    ASSERT_EQ(legacy.find("\"Momentum_Model\""), std::string::npos);
    ASSERT_EQ(legacy.find("\"Energy_Balance_Mode\""), std::string::npos);

    EXPECT_EQ(MixedAdsorbentLayer{}.momentum_model_, MomentumModel::Ergun);

    UserData data;
    LoadFromJson(data, deck.string());
    ASSERT_FALSE(data.bed.layers_.layers_.empty());
    for (const auto& layer : data.bed.layers_.layers_)
      EXPECT_EQ(layer.momentum_model_, MomentumModel::KarmenCozeny) << layer.name_;
    EXPECT_EQ(data.bed.equations_.energy_balance_mode_, EnergyBalanceMode::NON_ISOTHERMAL);
    EXPECT_TRUE(data.bed.equations_.isothermal_use_initial_temperature_);

    size_t inlets = 0;
    for (auto& step : data.cycle.step_list_) {
      for (auto* bc : {step.top_boundary_condition_.get(), step.bottom_boundary_condition_.get()}) {
        if (auto* fbc = dynamic_cast<FlowBoundaryCondition*>(bc)) {
          for (auto& s : fbc->streams_) {
            ++inlets;
            EXPECT_EQ(s->temperature_profile_, nullptr);
            EXPECT_EQ(s->AssumedInletTemperature(1e3), s->temperature_.Get());
          }
        } else if (auto* pbc = dynamic_cast<PressureBoundaryCondition*>(bc)) {
          EXPECT_EQ(pbc->temperature_profile_, nullptr);
        }
      }
    }
    EXPECT_GT(inlets, 0u);

    // Re-saving records the legacy choices explicitly in the new format.
    std::string json;
    SaveToString(data, json);
    EXPECT_NE(json.find("\"Momentum_Model\": 0"), std::string::npos);
    EXPECT_NE(json.find("\"Energy_Balance_Mode\": 0"), std::string::npos);
  }

} // namespace adsim
