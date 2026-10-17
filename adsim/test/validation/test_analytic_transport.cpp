/*
 * Skarstrom - pressure and temperature swing adsorption simulator
 * Copyright (c) 2024-2026 Sam Affleck and contributors
 *
 * Licensed under the terms in the LICENSE file at the repository root.
 */

// Analytical validation of the transport equations.
//
// These tests target the parts of the model that the existing species mass
// balance cannot see. A mass balance confirms that whatever entered the bed is
// accounted for; it says nothing about *when* it arrives, or about the
// relationship between flow and pressure gradient. Both are pinned here against
// closed-form results.

#include "validation_helpers.h"

namespace adsim {

  using namespace adsim::validation;

  class AnalyticTransportTest : public UserDataTest {
  public:
    LayerID layer_id{};
    AdsorbentID ads_id{};

    static constexpr double kLength = 1.0;       // [m]
    static constexpr double kDiameter = 0.1;     // [m]
    static constexpr double kVoidFraction = 0.4; // [-] bed void; particle void set to zero
    static constexpr double kTemperature = 298.0;// [K]

    // A bed of inert, non-adsorbing packing: no isotherm and no mass transfer
    // model are registered, so the loading stays identically zero and the gas
    // sees nothing but void volume and flow resistance.
    void BuildInertBed(int n_cells, double particle_diameter) {
      AddComponents({"Nitrogen", "Oxygen"});

      auto ids = AddLayer("Inert Packing");
      layer_id = ids.first;
      ads_id = ids.second;

      auto& layer = GetLayer(layer_id);
      layer.domain_.length = kLength;
      layer.domain_.n_cells = n_cells;
      layer.eb_ = kVoidFraction;
      layer.particle_diameter_ = particle_diameter;
      std::get<CylinderBed>(data.bed.shape_).diameter = kDiameter;

      auto& adsorbent = GetAdsorbent(ads_id);
      adsorbent.ep_ = 0.0;  // total void fraction == bed void fraction
      adsorbent.solid_density_ = 1000.0;
      adsorbent.heat_capacity_ = 900.0;

      // Isothermal: no heat exchange with the wall or the surroundings, and the
      // feed enters at the bed temperature, so T stays uniform and constant.
      AddConstantInternalHeatTransferCoefficientModel(0.0, layer_id);
      AddConstantExternalHeatTransferCoefficientModel(0.0, layer_id);
    }

    double CrossSectionalArea() const { return 0.25 * M_PI * kDiameter * kDiameter; }
    double BedVolume() const { return CrossSectionalArea() * kLength; }

    // Bottom inlet at a prescribed molar flow, top outlet at a fixed pressure.
    void AddFlowThroughStep(const std::string& step_name, double step_time,
                            double molar_flow, double outlet_pressure,
                            const std::unordered_map<std::string, double>& feed) {
      AddStep(step_name, step_time);
      Step& step = data.cycle.GetStep(step_name);

      auto inlet = CreateNoValveInletStream(CreateMolarFlowBC(molar_flow), feed, kTemperature);
      auto inlet_bc = std::make_unique<FlowBoundaryCondition>();
      inlet_bc->direction_ = BoundaryDirection::INLET;
      inlet_bc->location_ = BoundaryLocation::BOTTOM;
      inlet_bc->streams_.push_back(std::move(inlet));
      step.bottom_boundary_condition_ = std::move(inlet_bc);

      auto outlet_bc = std::make_unique<PressureBoundaryCondition>();
      outlet_bc->direction_ = BoundaryDirection::OUTLET;
      outlet_bc->location_ = BoundaryLocation::TOP;
      outlet_bc->pressure_ = CreateConstantBoundaryValue(outlet_pressure);
      for (const auto& [component, y] : feed) outlet_bc->mole_fraction_[component] = y;
      step.top_boundary_condition_ = std::move(outlet_bc);
    }
  };

  // ===========================================================================
  // 1. Mean residence time of a non-adsorbing tracer
  // ===========================================================================
  //
  // Step-change a dilute, non-adsorbing tracer into the feed of an isothermal,
  // essentially isobaric bed. A tracer balance integrated over the whole
  // transient gives, with no approximation beyond a constant molar throughput,
  //
  //     int_0^inf ( 1 - y_out(t)/y_feed ) dt  =  N_void / n_dot
  //
  // where N_void = et * V * P / (R T) is the moles of gas the column holds. The
  // integral is the classic first moment of the breakthrough curve, and it is
  // *invariant to dispersion*: physical axial dispersion, numerical smearing
  // from the upwind scheme, and cell count all change the shape of the curve but
  // none of them move its centre of mass. That makes it a sharp test of the
  // quantities that do set the timing -- void volume, cross-sectional area, and
  // the velocity/flow coupling -- while being immune to everything that a
  // grid-converged shape comparison would be sensitive to.

  TEST_F(AnalyticTransportTest, InertTracerMeanResidenceTimeMatchesVoidVolume) {
    data.info.title = "Validation: inert tracer residence time";

    BuildInertBed(/*n_cells=*/100, /*particle_diameter=*/0.002);

    const double P_outlet = 1.0e5;   // [Pa]
    const double molar_flow = 0.005; // [mol/s]
    const double y_feed = 0.05;      // [-] dilute tracer

    AddInitialConditions(P_outlet, kTemperature, kTemperature,
                         {{"Nitrogen", 1.0}, {"Oxygen", 0.0}}, 0.0);

    const double run_time = 200.0;  // [s], roughly eight residence times
    AddFlowThroughStep("Feed", run_time, molar_flow, P_outlet,
                       {{"Nitrogen", 1.0 - y_feed}, {"Oxygen", y_feed}});
    ConfigureCycle(1);
    data.settings.reltol = 1e-9;
    data.settings.max_step = 0.5;

    CoupledSolver solver{data};
    std::string error;
    ASSERT_TRUE(data.PerformChecks(error)) << error;
    ASSERT_GE(solver.IntegrateCycle(data.cycle, false), 0) << "Solver did not complete";

    const auto& log = solver.GetLog();

    // --- analytical residence time ------------------------------------------
    const double N_void = kVoidFraction * BedVolume() * P_outlet / (kR * kTemperature);
    const double tau_expected = N_void / molar_flow;

    // The derivation assumes the column is essentially isobaric; confirm that.
    const auto P = log.GetVariable(0, 0, VariableType::PRESSURE);
    const double P_inlet = P.values.back()[1];
    ASSERT_LT(std::abs(P_inlet - P_outlet) / P_outlet, 0.01)
      << "Pressure drop is too large for the isobaric residence-time derivation";

    // --- measured first moment ----------------------------------------------
    auto [times, y_out] = TopFaceMoleFraction(log, 0, 0, "Oxygen");
    ASSERT_FALSE(times.empty());

    // The tail must have died away, otherwise the truncated integral is biased.
    ASSERT_NEAR(y_out.back(), y_feed, 1e-3 * y_feed)
      << "Tracer has not fully broken through; the first moment would be truncated";

    std::vector<double> unbroken;
    unbroken.reserve(y_out.size());
    for (double y : y_out) unbroken.push_back(1.0 - y / y_feed);
    const double tau_measured = TimeIntegral(times, unbroken);

    EXPECT_NEAR(tau_measured, tau_expected, 0.01 * tau_expected)
      << "Tracer mean residence time disagrees with the column void inventory "
      << "(measured " << tau_measured << " s, analytic " << tau_expected << " s)";
  }

  // ===========================================================================
  // 2. Steady-state pressure drop
  // ===========================================================================
  //
  // At steady state the momentum equation is Darcy's law with the Carman-Kozeny
  // permeability K = eb^3 dp^2 / (180 (1-eb)^2). For an isothermal ideal gas the
  // molar flux N = u P / (R T) is constant along the column, so
  //
  //     N = -( K / (mu R T) ) P dP/dz   =>   P_in^2 - P_out^2 = 2 N L mu R T / K
  //
  // which is the compressible Darcy (Forchheimer-free) result. This is the only
  // test in the suite that exercises the momentum equation quantitatively; every
  // other check is satisfied by any velocity field that conserves mass.

  TEST_F(AnalyticTransportTest, SteadyFlowPressureDropMatchesCarmanKozeny) {
    data.info.title = "Validation: Carman-Kozeny pressure drop";

    const double particle_diameter = 5.0e-4;  // [m], small enough for a real dP
    const int n_cells = 100;
    BuildInertBed(n_cells, particle_diameter);

    const double P_outlet = 1.0e5;  // [Pa]
    const double molar_flow = 0.02; // [mol/s]

    AddInitialConditions(P_outlet, kTemperature, kTemperature,
                         {{"Nitrogen", 1.0}, {"Oxygen", 0.0}}, 0.0);

    AddFlowThroughStep("Feed", 400.0, molar_flow, P_outlet,
                       {{"Nitrogen", 1.0}, {"Oxygen", 0.0}});
    ConfigureCycle(1);
    data.settings.reltol = 1e-9;
    data.settings.max_step = 2.0;

    CoupledSolver solver{data};
    std::string error;
    ASSERT_TRUE(data.PerformChecks(error)) << error;
    ASSERT_GE(solver.IntegrateCycle(data.cycle, false), 0) << "Solver did not complete";

    const auto& log = solver.GetLog();
    const auto& dom = log.GetDomain();
    const auto P = log.GetVariable(0, 0, VariableType::PRESSURE);
    const auto mu = log.GetVariable(0, 0, VariableType::VISCOSITY);
    ASSERT_FALSE(P.values.empty());

    // Compare between the first and last interior cell centres. Summing the
    // face-by-face discrete Darcy law over that span is what the analytical
    // integral represents, so the two are directly comparable.
    const auto& P_final = P.values.back();
    const size_t first = 1;
    const size_t last = P_final.size() - 2;
    const double span = static_cast<double>(last - first) * dom.cell_width[1];

    const double viscosity = mu.values.back()[1];  // [Pa-s], a fluid property input
    const double permeability =
      std::pow(kVoidFraction, 3) * particle_diameter * particle_diameter /
      (180.0 * std::pow(1.0 - kVoidFraction, 2));
    const double molar_flux = molar_flow / CrossSectionalArea();  // [mol/(m2-s)]

    const double expected_dP2 =
      2.0 * molar_flux * span * viscosity * kR * kTemperature / permeability;
    const double measured_dP2 = P_final[first] * P_final[first] - P_final[last] * P_final[last];

    // The drop must be big enough to actually resolve.
    ASSERT_GT(P_final[first] - P_final[last], 0.02 * P_outlet)
      << "Pressure drop too small to be a meaningful test";

    EXPECT_NEAR(measured_dP2, expected_dP2, 0.01 * expected_dP2)
      << "Steady-state pressure drop disagrees with compressible Darcy / "
         "Carman-Kozeny (measured P_in - P_out = "
      << P_final[first] - P_final[last] << " Pa, analytic = "
      << std::sqrt(P_final[last] * P_final[last] + expected_dP2) - P_final[last] << " Pa)";
  }

  // ===========================================================================
  // 3. Stoichiometric breakthrough time for a linear isotherm
  // ===========================================================================
  //
  // The mass-balance tests already confirm that everything fed in is accounted
  // for. What they cannot see is *when* it arrives. For a dilute, isothermal
  // tracer on a Henry isotherm the stoichiometric time -- the first moment of
  // the breakthrough curve, weighted by the outlet molar flux -- is
  //
  //     t_st = V P [ et/(R T) + rho_b k exp(Q/RT) ] / n_dot
  //
  // i.e. the column's total capacity for the tracer divided by the rate it is
  // fed. Like the inert-tracer moment this is invariant to dispersion and to
  // mass-transfer resistance: a slower LDF coefficient smears the front but
  // cannot move its centre of mass. The isotherm's exp(Q/RT) factor is what
  // separates this from the inert case, so a sign error in the van 't Hoff
  // exponent would show up here and nowhere else.
  //
  // The heat of adsorption is set to zero so the column stays isothermal and the
  // capacity is a single number rather than a trajectory.

  TEST_F(AnalyticTransportTest, LinearIsothermBreakthroughTimeMatchesCapacity) {
    data.info.title = "Validation: stoichiometric breakthrough time";

    BuildInertBed(/*n_cells=*/120, /*particle_diameter=*/0.002);

    // Henry isotherm on the tracer only.
    const double henry_k = 4.15e-10;  // [mol/(kg-Pa)]
    const double Qiso = 20000.0;      // [J/mol]
    auto isotherm = std::make_unique<HenryIsotherm>();
    isotherm->params_["Oxygen"].k = henry_k;
    isotherm->params_["Oxygen"].Q = Qiso;
    AddIsotherm({"Oxygen"}, ads_id, std::move(isotherm));
    AddLinearDrivingForceModel(ads_id, "Oxygen", 0.5);
    SetHOA(ads_id, "Oxygen", 0.0);  // isothermal by construction

    auto& adsorbent = GetAdsorbent(ads_id);
    adsorbent.solid_density_ = 1000.0;
    const double rho_b = (1.0 - kVoidFraction) * 1000.0;

    const double P_outlet = 1.0e5;
    const double molar_flow = 0.005;  // [mol/s]
    const double y_feed = 0.02;       // dilute

    AddInitialConditions(P_outlet, kTemperature, kTemperature,
                         {{"Nitrogen", 1.0}, {"Oxygen", 0.0}}, 0.0);

    const double run_time = 900.0;
    AddFlowThroughStep("Feed", run_time, molar_flow, P_outlet,
                       {{"Nitrogen", 1.0 - y_feed}, {"Oxygen", y_feed}});
    ConfigureCycle(1);
    data.settings.reltol = 1e-9;
    data.settings.max_step = 1.0;

    CoupledSolver solver{data};
    std::string error;
    ASSERT_TRUE(data.PerformChecks(error)) << error;
    ASSERT_GE(solver.IntegrateCycle(data.cycle, false), 0) << "Solver did not complete";

    // --- analytical stoichiometric time -------------------------------------
    const double gas_capacity = kVoidFraction / (kR * kTemperature);          // [mol/(m3-Pa)]
    const double solid_capacity = rho_b * henry_k * std::exp(Qiso / (kR * kTemperature));
    const double t_st_expected =
      BedVolume() * P_outlet * (gas_capacity + solid_capacity) / molar_flow;

    // The adsorbed capacity must dominate, otherwise this degenerates into the
    // inert-tracer test and says nothing about the isotherm.
    ASSERT_GT(solid_capacity, 2.0 * gas_capacity)
      << "Adsorbed capacity is too small for this test to exercise the isotherm";
    ASSERT_LT(t_st_expected, 0.5 * run_time) << "Run is too short to see full breakthrough";

    // --- measured first moment ----------------------------------------------
    const auto& log = solver.GetLog();
    auto [times, y_out] = TopFaceMoleFraction(log, 0, 0, "Oxygen");
    const auto n_dot = log.GetFaceVariable(0, 0, VariableType::MOLAR_FLOW_RATE);
    ASSERT_FALSE(times.empty());
    ASSERT_NEAR(y_out.back(), y_feed, 1e-3 * y_feed)
      << "Tracer has not fully broken through; the first moment would be truncated";

    // Weight by the outlet molar flux, not just the mole fraction: adsorption
    // removes moles from the gas so the outlet flow is not exactly the feed
    // flow during the transient.
    std::vector<double> unbroken(times.size(), 0.0);
    for (size_t n = 0; n < times.size(); ++n) {
      const double tracer_out = n_dot.values[n].back() * y_out[n];
      unbroken[n] = 1.0 - tracer_out / (molar_flow * y_feed);
    }
    const double t_st_measured = TimeIntegral(times, unbroken);

    EXPECT_NEAR(t_st_measured, t_st_expected, 0.02 * t_st_expected)
      << "Stoichiometric breakthrough time disagrees with the analytic capacity "
      << "(measured " << t_st_measured << " s, analytic " << t_st_expected << " s)";
  }

  // ===========================================================================
  // 4. Valve blowdown
  // ===========================================================================
  //
  // The valve law is linear in the pressure difference,
  // n_dot = cv (P_up - P_down) / 1e5. Discharging a rigid, isothermal vessel
  // through it therefore gives a first-order decay with an analytically known
  // time constant:
  //
  //     (V_void / R T) dP/dt = -(cv / 1e5) (P - P_down)
  //     P(t) = P_down + (P0 - P_down) exp(-t / tau),   tau = V_void * 1e5 / (R T cv)
  //
  // This exercises three things no other test touches: the valve boundary
  // condition, the coupling between a boundary molar flow and the interior
  // pressure field, and the accumulation term under a large pressure swing.
  // A heavy, high-heat-capacity packing pins the temperature so the isothermal
  // derivation holds; the test asserts that it did.

  TEST_F(AnalyticTransportTest, ValveBlowdownFollowsAnalyticExponentialDecay) {
    data.info.title = "Validation: valve blowdown";

    BuildInertBed(/*n_cells=*/40, /*particle_diameter=*/0.004);

    // Heavy thermal ballast: the gas holds ~0.3 % of the bed's heat capacity, so
    // expansion cooling cannot move the temperature appreciably.
    auto& adsorbent = GetAdsorbent(ads_id);
    adsorbent.solid_density_ = 1200.0;
    adsorbent.heat_capacity_ = 1000.0;

    const double P0 = 5.0e5;        // [Pa]
    const double P_downstream = 1.0e5;
    const double cv = 1.5e-3;       // valve coefficient

    AddInitialConditions(P0, kTemperature, kTemperature,
                         {{"Nitrogen", 1.0}, {"Oxygen", 0.0}}, 0.0);

    const double run_time = 400.0;
    AddStep("Blowdown", run_time);
    Step& step = data.cycle.GetStep("Blowdown");

    auto closed = std::make_unique<ClosedBoundaryCondition>();
    closed->location_ = BoundaryLocation::BOTTOM;
    step.bottom_boundary_condition_ = std::move(closed);

    auto valve_stream = CreateValveOutletStream(CreateConstantBoundaryValue(P_downstream), cv);
    auto outlet_bc = std::make_unique<FlowBoundaryCondition>();
    outlet_bc->direction_ = BoundaryDirection::OUTLET;
    outlet_bc->location_ = BoundaryLocation::TOP;
    outlet_bc->streams_.push_back(std::move(valve_stream));
    step.top_boundary_condition_ = std::move(outlet_bc);

    ConfigureCycle(1);
    data.settings.reltol = 1e-10;
    data.settings.max_step = 1.0;

    CoupledSolver solver{data};
    std::string error;
    ASSERT_TRUE(data.PerformChecks(error)) << error;
    ASSERT_GE(solver.IntegrateCycle(data.cycle, false), 0) << "Solver did not complete";

    const auto& log = solver.GetLog();
    const auto& dom = log.GetDomain();
    const auto P = log.GetVariable(0, 0, VariableType::PRESSURE);
    const auto T = log.GetVariable(0, 0, VariableType::TEMPERATURE);
    ASSERT_FALSE(P.values.empty());

    // --- analytical decay ----------------------------------------------------
    const double void_volume = kVoidFraction * BedVolume();
    const double tau = void_volume * 1.0e5 / (kR * kTemperature * cv);
    ASSERT_LT(tau, 0.5 * run_time) << "Run is too short to resolve several time constants";

    std::vector<double> P_mean, T_mean;
    for (size_t n = 0; n < P.values.size(); ++n) {
      P_mean.push_back(InteriorMean(P.values[n], dom.volume));
      T_mean.push_back(InteriorMean(T.values[n], dom.volume));
    }

    // The isothermal derivation must actually hold.
    for (double t : {50.0, 200.0, 400.0}) {
      EXPECT_NEAR(InterpolateSeries(T.times, T_mean, t), kTemperature, 1.0)
        << "Blowdown is not isothermal at t = " << t << " s; the closed form does not apply";
    }

    for (double t : {20.0, 60.0, 150.0, 300.0}) {
      const double expected = P_downstream + (P0 - P_downstream) * std::exp(-t / tau);
      const double actual = InterpolateSeries(P.times, P_mean, t);
      EXPECT_NEAR(actual, expected, 0.02 * (expected - P_downstream))
        << "Blowdown pressure deviates from P_down + (P0 - P_down) exp(-t/tau) at t = "
        << t << " s (simulated " << actual << " Pa, analytic " << expected << " Pa, tau = "
        << tau << " s)";
    }
  }

} // namespace adsim
