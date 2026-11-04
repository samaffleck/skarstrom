/*
 * Skarstrom - pressure and temperature swing adsorption simulator
 * Copyright (c) 2024-2026 Sam Affleck and contributors
 *
 * Licensed under the terms in the LICENSE file at the repository root.
 */

// Analytical verification of the dynamic tank (adsim/objects/tank.h): a
// well-mixed constant-volume vessel whose species moles n_i and temperature T
// are DAE unknowns next to the bed's, P = N R T / V. Its balances are
//
//     dn_i/dt = sum_p F_p,i,
//     d/dt( sum_i n_i c_v,i T ) = sum_p sum_i F_p,i c_p,i T_p + Q,
//
// with F_p,i the flow of component i into the tank through port or stream p
// and T_p the temperature of the gas that crosses p.
//
// 1. Valve exchange with a rigid, isothermal, inert bed (closed otherwise),
//    linear valve n_dot = k (P_bed - P_tank), k = Cv / 1e5:
//        d(P_bed - P_tank)/dt = -k R T (1/V_bed + 1/V_tank) (P_bed - P_tank),
//    an exponential with 1/tau = k R T (1/V_bed + 1/V_tank), ending at the
//    volume-weighted mean pressure. Guard: tau is half the single-vessel
//    V_bed / (k R T) for equal volumes. And the tank venting through its own
//    valve: P - P_d = (P_0 - P_d) exp(-t k R T / V).
// 2. Mixing: a tank fed at a constant flow F of pure A, drawn off by a
//    pressure hold (N constant to the hold band): 1 - y_A = exp(-F t / N), and
//    the mean residence time int_0^inf (1 - y_A) dt = N / F.
// 3. Mole balance through bed + tank: a closed loop (a Langmuir bed whose two
//    ends draw from and send to the same tank, with a step that feeds and
//    drains the tank and a step in which the outlet valve flows backwards)
//    conserves the total and the component inventories of bed + tank.
// 4. Energy: an adiabatic bed and an adiabatic tank exchanging gas conserve
//    U_bed + U_tank; a closed tank losing heat (UA to ambient) and a tank with
//    a wall node follow the one- and two-body exponentials.
// 5. Adiabatic filling: a tank filled at T_in from N_0 at T_0 obeys
//    N c_v T = N_0 c_v T_0 + c_p T_in (N - N_0), so T -> gamma T_in as
//    N / N_0 -> infinity.
// 7. CSS: with a large tank the bed reaches its own cyclic steady state long
//    before the tank does; the CSS test, which includes the tank's state,
//    must not stop there.
//
// Tolerances: see each test (from the neglected terms named there, or the
// integrator tolerance for the exact balances).

#include "validation_helpers.h"

#include <algorithm>
#include <cmath>

namespace adsim {

  using namespace adsim::validation;

  namespace {

    constexpr double kT = 298.15;
    constexpr double kEps = 0.4;
    constexpr double kRhoS = 1000.0;
    constexpr double kL = 1.0, kD = 0.1;
    constexpr int kCells = 10;

    double BedVolume() { return 0.25 * M_PI * kD * kD * kL; }
    double BedGasVolume() { return kEps * BedVolume(); }

    // A rigid inert bed (no macropores: void fraction exactly kEps) of the
    // given components, coarse packing (negligible pressure drop).
    struct BedOptions {
      std::vector<std::string> components{"Nitrogen"};
      bool isothermal = true;
      bool langmuir = false;
      double P0 = 1e5;
      std::unordered_map<std::string, double> y0{{"Nitrogen", 1.0}};
      double solid_heat_capacity = 900.0;
    };

    void BuildBed(UserData& data, const BedOptions& o) {
      data.bed.wall_.thickness_ = 0.005;
      for (const auto& c : o.components) data.fluid.AddComponent(c);
      auto ads = data.adsorbent_directory.Add("Packing");
      auto lid = data.bed.layers_.Add("Packing");
      auto& layer = data.bed.layers_.Get(lid);
      layer.momentum_model_ = MomentumModel::KarmenCozeny;
      layer.domain_.length = kL;
      layer.domain_.n_cells = kCells;
      layer.eb_ = kEps;
      layer.particle_diameter_ = 0.01;
      layer.composition_[ads] = 1.0;
      std::get<CylinderBed>(data.bed.shape_).diameter = kD;
      auto& adsorbent = data.adsorbent_directory.Get(ads);
      adsorbent.ep_ = 0.0;
      adsorbent.solid_density_ = kRhoS;
      adsorbent.heat_capacity_ = o.solid_heat_capacity;
      adsorbent.thermal_conductivity_ = 0.0;
      if (o.langmuir) {
        // Different affinities: the bed separates the two gases.
        auto iso = std::make_unique<ExtendedLangmuirIsotherm>();
        iso->params_[o.components[0]].m1 = 1.0;
        iso->params_[o.components[0]].b1 = 1e-5;
        iso->params_[o.components[0]].q1 = 0.0;
        iso->params_[o.components[1]].m1 = 1.0;
        iso->params_[o.components[1]].b1 = 2e-6;
        iso->params_[o.components[1]].q1 = 0.0;
        adsorbent.adsorption_source_.AddIsotherm(o.components, std::move(iso));
        for (const auto& c : o.components) {
          auto ldf = std::make_unique<LinearDrivingForce>();
          auto mtc = std::make_unique<MassTransferCoefficientConstant>();
          mtc->constant_mtc_ = 1.0;
          ldf->mtc_model_ = std::move(mtc);
          adsorbent.adsorption_source_.mass_transfer_models_[c] = std::move(ldf);
          adsorbent.adsorption_source_.heat_of_adsorption_[c] = 0.0;
        }
      }
      if (!o.isothermal) {
        // Adiabatic bed: no bed-wall or wall-ambient exchange.
        auto hi = std::make_unique<ConstantInternalHeatTransferModel>();
        hi->heat_transfer_coefficient_ = 0.0;
        layer.equations_.internal_heat_transfer_model_ = std::move(hi);
        auto ho = std::make_unique<ConstantExternalHeatTransferModel>();
        ho->heat_transfer_coefficient_ = 0.0;
        layer.equations_.external_heat_transfer_model_ = std::move(ho);
      }
      data.bed.equations_.energy_balance_mode_ =
          o.isothermal ? EnergyBalanceMode::ISOTHERMAL : EnergyBalanceMode::NON_ISOTHERMAL;
      data.ic.pressure = o.P0;
      data.ic.temperature = kT;
      data.ic.wall_temperature = kT;
      data.ic.velocity = 0.0;
      for (const auto& [c, y] : o.y0) data.ic.mole_fraction[c] = y;
      data.settings.reltol = 1e-9;
      data.settings.max_step = 1.0;
      data.settings.max_log_points_per_step_ = 400;
    }

    std::unique_ptr<ConstantBoundaryValue> Constant(double v) {
      auto bv = std::make_unique<ConstantBoundaryValue>();
      bv->constant_value_ = v;
      return bv;
    }

    void Close(Step& step, BoundaryLocation loc) {
      auto closed = std::make_unique<ClosedBoundaryCondition>();
      closed->location_ = loc;
      (loc == BoundaryLocation::TOP ? step.top_boundary_condition_ : step.bottom_boundary_condition_) =
          std::move(closed);
    }

    // A flow boundary with one stream on `loc` of `step`.
    template <class S>
    S* AddStream(Step& step, BoundaryLocation loc, BoundaryDirection dir, std::unique_ptr<S> stream) {
      S* raw = stream.get();
      auto bc = std::make_unique<FlowBoundaryCondition>();
      bc->direction_ = dir;
      bc->location_ = loc;
      bc->streams_.push_back(std::move(stream));
      (loc == BoundaryLocation::TOP ? step.top_boundary_condition_ : step.bottom_boundary_condition_) =
          std::move(bc);
      return raw;
    }

    // A linear valve stream (Cv in mol/(s bar)) ported to `tank`.
    std::unique_ptr<ValveBoundaryStream> PortValve(double cv, size_t tank, const std::string& name = "Valve") {
      auto v = std::make_unique<ValveBoundaryStream>();
      v->name_ = name;
      v->valve_.cv_ = cv;
      v->temperature_ = kT;
      v->tank_port_.tank_id_ = tank;
      return v;
    }

    // A set molar flow stream ported to `tank`.
    std::unique_ptr<BoundaryStream> PortFlow(double F, size_t tank, const std::string& name) {
      auto s = std::make_unique<BoundaryStream>();
      s->name_ = name;
      auto flow = std::make_unique<MolarFlowBoundaryValue>();
      static_cast<ConstantBoundaryValue*>(flow->profile_.get())->constant_value_ = F;
      s->boundary_flow_ = std::move(flow);
      s->temperature_ = kT;
      s->tank_port_.tank_id_ = tank;
      return s;
    }

    Tank& NewTank(UserData& data, double V, double P0, double T0,
                  const std::unordered_map<std::string, double>& y0) {
      const size_t id = data.tanks.Add("Tank");
      Tank& t = *data.tanks.Find(id);
      t.volume_ = V;
      t.initial_pressure_ = P0;
      t.initial_temperature_ = T0;
      for (const auto& c : data.fluid.components_) t.initial_mole_fraction_[c] = 0.0;
      for (const auto& [c, y] : y0) t.initial_mole_fraction_[c] = y;
      return t;
    }

    void Run(UserData& data, CoupledSolver& solver) {
      std::string error;
      ASSERT_TRUE(data.PerformChecks(error)) << error;
      ASSERT_GE(solver.IntegrateCycle(data.cycle, false), 0) << "Solver did not complete";
    }

    // Volume-weighted interior mean of a logged cell-centred row.
    double BedMean(const std::vector<double>& row) {
      double s = 0.0;
      for (size_t i = 1; i + 1 < row.size(); ++i) s += row[i];
      return s / static_cast<double>(row.size() - 2);
    }

    // Linear interpolation of a series at t.
    double At(const std::vector<double>& t, const std::vector<double>& v, double x) {
      return interp_linear(t, v, x);
    }

    // Total moles of bed + tank (gas + adsorbed) at logged row k of a step,
    // per component; uses only logged state.
    double BedInventory(const EquationLog& log, int cycle, int step, const std::string& comp, size_t k) {
      const auto& dom = log.GetDomain();
      const auto c = log.GetComponentVariable(cycle, step, ComponentVariableType::COMPONENT_MOLAR_CONCENTRATION, comp);
      const auto q = log.GetComponentVariable(cycle, step, ComponentVariableType::ADSORBED_CONCENTRATION, comp);
      double n = 0.0;
      for (size_t i = 1; i + 1 < c.values[k].size(); ++i) {
        n += c.values[k][i] * dom.volume[i] * dom.et[i];
        if (!q.values.empty()) n += q.values[k][i] * dom.volume[i] * dom.bulk_density[i];
      }
      return n;
    }

  }  // namespace

  // ---------------------------------------------------------------------------
  // 1. Valve exchange with the bed; tank vent
  // ---------------------------------------------------------------------------
  //
  // Neglected: the pressure drop inside the bed. The valve reads the bed-end
  // face, whose pressure differs from the bed's mean by about half the Darcy
  // drop of the flow, Delta p_D / 2 ~ 150 mu (1-e)^2 u L / (2 e^3 d^2) ~ 2 Pa
  // at the largest flow (u ~ 0.06 m/s, d = 1 cm), proportional to the flow and
  // so to P_bed - P_tank: a relative error of the rate of ~2 Pa / 1e5 Pa =
  // 2e-5. The bound, 1e-3, leaves 50x for the bed's acoustic transient at the
  // start (the bed's internal pressure relaxes in ~L / c ~ 3 ms).
  TEST(DynamicTank, ValveExchangeWithTheBedDecaysWithTheTwoVesselRate) {
    UserData data;
    BedOptions o;
    o.P0 = 2e5;
    BuildBed(data, o);
    const double V_bed = BedGasVolume();
    const double V_tank = V_bed;  // equal volumes: tau is half the single-vessel value
    Tank& tank = NewTank(data, V_tank, 1e5, kT, {{"Nitrogen", 1.0}});
    const size_t tank_id = tank.id_;
    const double tau_target = 20.0;
    const double k = 1.0 / (tau_target * kR * kT * (1.0 / V_bed + 1.0 / V_tank));  // [mol/(s Pa)]
    const double tau = 1.0 / (k * kR * kT * (1.0 / V_bed + 1.0 / V_tank));

    data.cycle.AddStep("Exchange");
    Step& step = data.cycle.GetStep("Exchange");
    step.step_time_ = 10.0 * tau;
    Close(step, BoundaryLocation::BOTTOM);
    AddStream(step, BoundaryLocation::TOP, BoundaryDirection::OUTLET, PortValve(k * 1e5, tank_id));
    data.cycle.number_of_cycles_ = 1;
    data.settings.max_step = 0.2;

    CoupledSolver solver{data};
    Run(data, solver);
    const auto& log = solver.GetLog();
    const TankLogData* t = log.GetTankLog(0, 0, tank_id);
    ASSERT_NE(t, nullptr);
    const auto P = log.GetVariable(0, 0, VariableType::PRESSURE);
    const auto& times = P.times;
    ASSERT_EQ(times.size(), t->pressure.size());
    std::vector<double> dP(times.size());
    for (size_t i = 0; i < times.size(); ++i) dP[i] = BedMean(P.values[i]) - t->pressure[i];

    const double t1 = 2.0 * tau, t2 = 6.0 * tau;
    const double tau_measured = (t2 - t1) / std::log(At(times, dP, t1) / At(times, dP, t2));
    const double P_end = (V_bed * 2e5 + V_tank * 1e5) / (V_bed + V_tank);
    std::cout << "[Tank exchange] tau " << tau_measured << " s (expected " << tau << "), P_tank,end "
              << t->pressure.back() << " (expected " << P_end << ")\n";
    EXPECT_NEAR(tau_measured / tau, 1.0, 1e-3);
    // Moles are conserved: the final pressure is the volume-weighted mean (to
    // exp(-10) of the initial difference, 4.5 Pa, plus the rate error above).
    EXPECT_NEAR(t->pressure.back(), P_end, 6.0);
    EXPECT_NEAR(BedMean(P.values.back()), P_end, 6.0);
    // GUARD: the single-vessel time constant (tank at fixed pressure) is twice tau.
    const double tau_single = V_bed / (k * kR * kT);
    EXPECT_GT(std::abs(tau_measured - tau_single) / tau_single, 0.4);
  }

  // Tank alone (bed closed), venting through its own linear valve with a
  // check valve. Neglected: nothing in the balance; the check valve's
  // smoothing acts only within 1 Pa of P_d, never reached (P - P_d >= 2e5
  // e^-8 = 67 Pa). Bound: 1e-4 on the rate (integrator tolerance 1e-9).
  TEST(DynamicTank, TankVentDecaysExponentiallyToTheDownstreamPressure) {
    UserData data;
    BuildBed(data, BedOptions{});
    const double V = 0.05, P0 = 3e5, P_d = 1e5;
    Tank& tank = NewTank(data, V, P0, kT, {{"Nitrogen", 1.0}});
    const size_t tank_id = tank.id_;
    auto valve = std::make_unique<TankValve>();
    valve->name_ = "Vent";
    valve->cv_ = 0.5;  // [mol/(s bar)]
    valve->downstream_pressure_ = P_d;
    tank.streams_.push_back(std::move(valve));
    const double k = 0.5 / 1e5;
    const double tau = V / (k * kR * kT);

    data.cycle.AddStep("Vent");
    Step& step = data.cycle.GetStep("Vent");
    step.step_time_ = 8.0 * tau;
    Close(step, BoundaryLocation::BOTTOM);
    Close(step, BoundaryLocation::TOP);
    data.cycle.number_of_cycles_ = 1;

    CoupledSolver solver{data};
    Run(data, solver);
    const auto& log = solver.GetLog();
    const TankLogData* t = log.GetTankLog(0, 0, tank_id);
    ASSERT_NE(t, nullptr);
    const auto& times = log.GetLog(0, 0).step_data.times;
    std::vector<double> dP(times.size());
    for (size_t i = 0; i < times.size(); ++i) dP[i] = t->pressure[i] - P_d;
    double worst = 0.0;
    for (size_t i = 0; i < times.size(); ++i)
      worst = std::max(worst, std::abs(dP[i] - (P0 - P_d) * std::exp(-times[i] / tau)) / (P0 - P_d));
    const double tau_measured = (6.0 * tau) / std::log(At(times, dP, tau) / At(times, dP, 7.0 * tau));
    std::cout << "[Tank vent] tau " << tau_measured << " s (expected " << tau << "), worst deviation " << worst << "\n";
    EXPECT_NEAR(tau_measured / tau, 1.0, 1e-4);
    EXPECT_LT(worst, 1e-5);
    // The vented moles are logged on the valve: its outflow integrates to the
    // inventory the tank lost.
    const StreamData* vent = solver.GetBoundaryLog().GetStreamData(0, 0, TankStreamLogName("Tank", "Vent"));
    ASSERT_NE(vent, nullptr);
    const double vented = IntegrateVariableInTime(vent->molar_flow_, vent->time_series_);
    const double lost = t->total_moles.front() - t->total_moles.back();
    EXPECT_NEAR(vented / lost, 1.0, 1e-4);
    // GUARD: the pressure fell by most of the initial difference.
    EXPECT_LT(dP.back(), 1e-3 * (P0 - P_d));
  }

  // ---------------------------------------------------------------------------
  // 2. Mixing: CSTR washout under a pressure hold
  // ---------------------------------------------------------------------------
  //
  // Neglected: the variation of N within the hold band (P within band P_set
  // above the set point: |dN| / N <= 1e-4), which perturbs the exponent
  // F t / N by at most 1e-4 of itself: |d y| <= 1e-4 max(x e^-x) = 3.7e-5.
  // The initial rise of the pressure into the band takes band N / F ~ 1e-2 s.
  // Bound: 1e-4 on y and on the mean residence time.
  TEST(DynamicTank, HeldTankWashesOutAsAStirredTank) {
    UserData data;
    BedOptions o;
    o.components = {"Nitrogen", "Oxygen"};
    o.y0 = {{"Nitrogen", 1.0}, {"Oxygen", 0.0}};
    BuildBed(data, o);
    const double V = 0.01, P_set = 2e5, F = 0.01;
    Tank& tank = NewTank(data, V, P_set, kT, {{"Nitrogen", 1.0}});
    const size_t tank_id = tank.id_;
    auto feed = std::make_unique<TankMakeUp>();
    feed->name_ = "Feed";
    feed->molar_flow_ = F;
    feed->temperature_ = kT;
    feed->mole_fraction_["Nitrogen"] = 0.0;
    feed->mole_fraction_["Oxygen"] = 1.0;
    tank.streams_.push_back(std::move(feed));
    auto hold = std::make_unique<TankHold>();
    hold->name_ = "Product";
    hold->set_pressure_ = P_set;
    hold->cv_ = 10.0;
    hold->downstream_pressure_ = 1e5;
    tank.streams_.push_back(std::move(hold));
    const double N = P_set * V / (kR * kT);
    const double tau = N / F;

    data.cycle.AddStep("Wash");
    Step& step = data.cycle.GetStep("Wash");
    step.step_time_ = 12.0 * tau;
    Close(step, BoundaryLocation::BOTTOM);
    Close(step, BoundaryLocation::TOP);
    data.cycle.number_of_cycles_ = 1;
    data.settings.max_log_points_per_step_ = 4000;
    data.settings.max_step = tau / 400.0;

    CoupledSolver solver{data};
    Run(data, solver);
    const auto& log = solver.GetLog();
    const TankLogData* t = log.GetTankLog(0, 0, tank_id);
    ASSERT_NE(t, nullptr);
    const auto& times = log.GetLog(0, 0).step_data.times;
    const auto& yO2 = t->mole_fraction.at("Oxygen");
    double worst = 0.0, P_max = 0.0;
    std::vector<double> remaining(times.size());
    for (size_t i = 0; i < times.size(); ++i) {
      worst = std::max(worst, std::abs(yO2[i] - (1.0 - std::exp(-times[i] / tau))));
      remaining[i] = 1.0 - yO2[i];
      P_max = std::max(P_max, t->pressure[i]);
    }
    const double T_end = times.back();
    const double moment = TimeIntegral(times, remaining);
    const double moment_expected = tau * (1.0 - std::exp(-T_end / tau));
    const StreamData* product = solver.GetBoundaryLog().GetStreamData(0, 0, TankStreamLogName("Tank", "Product"));
    ASSERT_NE(product, nullptr);
    const double drawn = IntegrateVariableInTime(product->molar_flow_, product->time_series_);
    std::cout << "[Tank washout] worst |y - y_exact| " << worst << ", mean residence time " << moment
              << " s (expected " << moment_expected << "), P_max - P_set " << P_max - P_set
              << " Pa, drawn / fed " << drawn / (F * T_end) << "\n";
    EXPECT_LT(worst, 1e-4);
    EXPECT_NEAR(moment / moment_expected, 1.0, 1e-4);
    // The hold kept the pressure within its band and drew off what was fed.
    EXPECT_LE(P_max - P_set, 1e-4 * P_set);
    EXPECT_NEAR(drawn / (F * T_end), 1.0, 1e-4);
    // GUARD: the composition changed completely, and without the draw-off
    // the feed would have raised the pressure by F t R T / V = 1.2e6 Pa.
    EXPECT_GT(yO2.back(), 0.9999);
    EXPECT_GT(F * T_end * kR * kT / V, 1e4 * (P_max - P_set));
  }

  // ---------------------------------------------------------------------------
  // 3. Mole balance through bed + tank
  // ---------------------------------------------------------------------------
  //
  // A closed loop: nothing enters or leaves bed + tank, so the inventory of
  // each component (bed gas + adsorbed + tank) is a constant of the motion.
  //   Circulate: the bed's bottom draws a set flow from the tank, its top
  //              sends its outflow back through a valve (the tank is fed and
  //              drained in the same step);
  //   Draw:      the bottom sends a set flow into the tank, the top is closed
  //              (the bed falls below the tank);
  //   Equalise:  the top valve, an OUTLET, now flows backwards (tank -> bed).
  // The tank is used at both bed ends, so the solver's band spans the system.
  // The ports take the bed's own boundary flows (the streams' flow functions),
  // so the inventory is conserved to the integrator: its nonlinear invariant
  // (P y / T and the loadings) is kept to the local error, reltol 1e-9 of the
  // inventory per step, over ~10^3 steps: bound 1e-6 relative.
  TEST(DynamicTank, BedAndTankConserveTheInventoryInAClosedLoop) {
    UserData data;
    BedOptions o;
    o.components = {"CarbonDioxide", "Nitrogen"};
    o.y0 = {{"CarbonDioxide", 0.3}, {"Nitrogen", 0.7}};
    o.langmuir = true;
    o.P0 = 1.5e5;
    BuildBed(data, o);
    Tank& tank = NewTank(data, 2.0 * BedGasVolume(), 1.5e5, kT, {{"CarbonDioxide", 0.3}, {"Nitrogen", 0.7}});
    const size_t tank_id = tank.id_;

    data.cycle.AddStep("Circulate");
    {
      Step& s = data.cycle.GetStep("Circulate");
      s.step_time_ = 60.0;
      AddStream(s, BoundaryLocation::BOTTOM, BoundaryDirection::INLET, PortFlow(0.01, tank_id, "Loop in"));
      AddStream(s, BoundaryLocation::TOP, BoundaryDirection::OUTLET, PortValve(0.05, tank_id, "Loop out"));
    }
    data.cycle.AddStep("Draw");
    {
      Step& s = data.cycle.GetStep("Draw");
      s.step_time_ = 40.0;
      AddStream(s, BoundaryLocation::BOTTOM, BoundaryDirection::OUTLET, PortFlow(0.01, tank_id, "Draw"));
      Close(s, BoundaryLocation::TOP);
    }
    data.cycle.AddStep("Equalise");
    ValveBoundaryStream* back = nullptr;
    {
      Step& s = data.cycle.GetStep("Equalise");
      s.step_time_ = 60.0;
      Close(s, BoundaryLocation::BOTTOM);
      back = AddStream(s, BoundaryLocation::TOP, BoundaryDirection::OUTLET, PortValve(0.05, tank_id, "Back"));
    }
    data.cycle.number_of_cycles_ = 4;
    data.settings.max_step = 0.5;

    CoupledSolver solver{data};
    Run(data, solver);
    const auto& log = solver.GetLog();

    double worst = 0.0, worst_component = 0.0, throughput = 0.0, min_back_flow = 0.0, y_swing = 0.0;
    double total0 = -1.0;
    std::unordered_map<std::string, double> component0;
    for (int cycle : log.sys_data.cycle_index_set_) {
      for (int step : log.sys_data.step_index_set_) {
        const TankLogData* t = log.GetTankLog(cycle, step, tank_id);
        ASSERT_NE(t, nullptr);
        const auto& times = log.GetLog(cycle, step).step_data.times;
        for (size_t k = 0; k < times.size(); ++k) {
          double total = 0.0;
          for (const auto& c : o.components) {
            const double n = BedInventory(log, cycle, step, c, k) + t->moles.at(c)[k];
            total += n;
            if (!component0.count(c)) component0[c] = n;
            worst_component = std::max(worst_component, std::abs(n - component0[c]) / component0[c]);
          }
          if (total0 < 0.0) total0 = total;
          worst = std::max(worst, std::abs(total - total0) / total0);
        }
        for (const auto& c : o.components) {
          throughput += IntegrateVariableInTime(t->flow_in.at(c), times);
          const auto& y = t->mole_fraction.at(c);
          y_swing = std::max(y_swing, *std::max_element(y.begin(), y.end()) - *std::min_element(y.begin(), y.end()));
        }
      }
    }
    // The equalisation valve (an outlet) flowed backwards, into the bed.
    // (back->stream_data_ holds the last cycle's run.)
    for (double f : back->stream_data_.molar_flow_) min_back_flow = std::min(min_back_flow, f);
    std::cout << "[Tank loop] worst inventory drift " << worst << " (component " << worst_component
              << "), tank throughput " << throughput << " mol vs inventory " << total0
              << " mol, most negative outlet flow " << min_back_flow << " mol/s, tank y swing " << y_swing << "\n";
    EXPECT_LT(worst, 1e-6);
    EXPECT_LT(worst_component, 1e-6);
    // GUARDS: the tank moved several inventories of gas, the outlet valve
    // flowed backwards at a sizeable rate, and the bed separated the gases.
    EXPECT_GT(throughput, 2.0 * total0);
    EXPECT_LT(min_back_flow, -1e-3);
    EXPECT_GT(y_swing, 1e-2);
  }

  // ---------------------------------------------------------------------------
  // 4. Energy
  // ---------------------------------------------------------------------------
  //
  // An adiabatic, inert bed (light packing, so the gas carries a large share
  // of the heat capacity) and an adiabatic tank exchange gas through a valve:
  // the gas leaving the bed expands (the bed cools), the gas entering the tank
  // is compressed (the tank warms), and U_bed + U_tank stays constant:
  //   U_bed = sum_i dV_i (e C_i c_v + rho_b c_s) T_i,  U_tank = N c_v T.
  // Neglected: the bed's temperature equation is written in the advective
  // (non-conservative) form, rho c_p u dT/dz, which conserves energy only to
  // the discretisation error of the temperature gradient along the bed; with
  // uniform expansion the bed cools uniformly, so the gradient is small. The
  // bound, 1e-2 of the energy the tank gained, is set from the measured
  // closure and its halving under grid refinement (first order: printed).
  TEST(DynamicTank, AdiabaticBedAndTankConserveTheInternalEnergy) {
    auto run = [](int cells, double& closure, double& tank_gain, double& bed_cooling) {
      UserData data;
      BedOptions o;
      o.isothermal = false;
      o.P0 = 3e5;
      o.solid_heat_capacity = 1.0;  // [J/(kg K)]: the gas dominates the heat capacity
      BuildBed(data, o);
      data.bed.layers_.layers_.front().domain_.n_cells = cells;
      Tank& tank = NewTank(data, BedGasVolume(), 1e5, kT, {{"Nitrogen", 1.0}});
      tank.thermal_ = TankThermalMode::ADIABATIC;
      const size_t tank_id = tank.id_;
      data.cycle.AddStep("Exchange");
      Step& s = data.cycle.GetStep("Exchange");
      s.step_time_ = 100.0;
      Close(s, BoundaryLocation::BOTTOM);
      AddStream(s, BoundaryLocation::TOP, BoundaryDirection::OUTLET, PortValve(0.02, tank_id));
      data.cycle.number_of_cycles_ = 1;
      data.settings.max_step = 0.5;
      CoupledSolver solver{data};
      Run(data, solver);
      const auto& log = solver.GetLog();
      const TankLogData* t = log.GetTankLog(0, 0, tank_id);
      const double cv = data.fluid.MolarHeatCapacity("Nitrogen") - kR;
      const auto T = log.GetVariable(0, 0, VariableType::TEMPERATURE);
      const auto C = log.GetVariable(0, 0, VariableType::MOLAR_CONCENTRATION);
      const auto& dom = log.GetDomain();
      auto U_bed = [&](size_t k) {
        double U = 0.0;
        for (size_t i = 1; i + 1 < T.values[k].size(); ++i)
          U += dom.volume[i] * (dom.et[i] * C.values[k][i] * cv + dom.bulk_density[i] * dom.Cps[i]) * T.values[k][i];
        return U;
      };
      auto U_tank = [&](size_t k) { return t->total_moles[k] * cv * t->temperature[k]; };
      const size_t last = T.values.size() - 1;
      const double dU_bed = U_bed(last) - U_bed(0);
      const double dU_tank = U_tank(last) - U_tank(0);
      closure = std::abs(dU_bed + dU_tank) / std::abs(dU_tank);
      tank_gain = t->temperature.back() - t->temperature.front();
      bed_cooling = kT - BedMean(T.values.back());
    };
    double c10, c20, g10, g20, b10, b20;
    run(10, c10, g10, b10);
    run(20, c20, g20, b20);
    std::cout << "[Tank energy] closure |dU_bed + dU_tank| / |dU_tank|: 10 cells " << c10 << ", 20 cells " << c20
              << "; tank warmed " << g20 << " K, bed cooled " << b20 << " K\n";
    EXPECT_LT(c20, 1e-2);
    EXPECT_LT(c20, c10 + 1e-12);  // refinement does not worsen it
    // GUARDS: both vessels changed temperature by tens of kelvin.
    EXPECT_GT(g20, 10.0);
    EXPECT_GT(b20, 10.0);
  }

  // A closed tank (no ports in use, no streams) with heat loss to ambient:
  // N c_v dT/dt = UA (T_a - T), so T - T_a = (T_0 - T_a) exp(-UA t / (N c_v)).
  // And with a wall node, the two-body system
  //   C_g dT/dt  = h_i A (T_w - T),
  //   C_w dT_w/dt = h_i A (T - T_w) - h_o A (T_w - T_a),
  // solved in closed form by its two eigenmodes. Neglected: nothing (linear
  // ODEs); bound 1e-6 K (integrator tolerance 1e-9 on ~50 K).
  TEST(DynamicTank, ClosedTankCoolsAsTheOneAndTwoBodyExponentials) {
    for (bool wall : {false, true}) {
      UserData data;
      BedOptions o;
      o.isothermal = false;
      BuildBed(data, o);
      const double V = 0.05, P0 = 5e5, T0 = 350.0, Ta = 300.0, Tw0 = 310.0;
      Tank& tank = NewTank(data, V, P0, T0, {{"Nitrogen", 1.0}});
      const size_t tank_id = tank.id_;
      tank.ambient_temperature_ = Ta;
      double ua = 2.0, hiA = 5.0, hoA = 1.0, Cw = 2000.0;
      if (wall) {
        tank.thermal_ = TankThermalMode::WALL;
        tank.wall_area_ = 1.0;
        tank.h_in_ = hiA;
        tank.h_out_ = hoA;
        tank.wall_heat_capacity_ = Cw;
        tank.initial_wall_temperature_ = Tw0;
      } else {
        tank.thermal_ = TankThermalMode::HEAT_LOSS;
        tank.ua_ = ua;
      }
      data.cycle.AddStep("Cool");
      Step& s = data.cycle.GetStep("Cool");
      s.step_time_ = 2000.0;
      Close(s, BoundaryLocation::BOTTOM);
      Close(s, BoundaryLocation::TOP);
      data.cycle.number_of_cycles_ = 1;
      data.settings.max_step = 5.0;
      CoupledSolver solver{data};
      Run(data, solver);
      const auto& log = solver.GetLog();
      const TankLogData* t = log.GetTankLog(0, 0, tank_id);
      ASSERT_NE(t, nullptr);
      const auto& times = log.GetLog(0, 0).step_data.times;
      const double N = P0 * V / (kR * T0);
      const double Cg = N * (data.fluid.MolarHeatCapacity("Nitrogen") - kR);
      double worst = 0.0;
      for (size_t k = 0; k < times.size(); ++k) {
        double T_exact = 0.0, Tw_exact = 0.0;
        const double tt = times[k];
        if (!wall) {
          T_exact = Ta + (T0 - Ta) * std::exp(-ua * tt / Cg);
        } else {
          // x = (T - Ta, Tw - Ta): x' = M x.
          const double m11 = -hiA / Cg, m12 = hiA / Cg;
          const double m21 = hiA / Cw, m22 = -(hiA + hoA) / Cw;
          const double tr = m11 + m22, det = m11 * m22 - m12 * m21;
          const double disc = std::sqrt(tr * tr / 4.0 - det);
          const double l1 = tr / 2.0 + disc, l2 = tr / 2.0 - disc;
          // Eigenvectors (m12, l - m11).
          const double x0 = T0 - Ta, w0 = Tw0 - Ta;
          // Solve c1 v1 + c2 v2 = (x0, w0), v = (m12, l - m11).
          const double a11 = m12, a12 = m12, a21 = l1 - m11, a22 = l2 - m11;
          const double dd = a11 * a22 - a12 * a21;
          const double c1 = (x0 * a22 - a12 * w0) / dd;
          const double c2 = (a11 * w0 - a21 * x0) / dd;
          T_exact = Ta + c1 * a11 * std::exp(l1 * tt) + c2 * a12 * std::exp(l2 * tt);
          Tw_exact = Ta + c1 * a21 * std::exp(l1 * tt) + c2 * a22 * std::exp(l2 * tt);
          worst = std::max(worst, std::abs(t->wall_temperature[k] - Tw_exact));
        }
        worst = std::max(worst, std::abs(t->temperature[k] - T_exact));
      }
      std::cout << "[Tank " << (wall ? "wall" : "heat loss") << "] worst |T - T_exact| " << worst
                << " K, T_end " << t->temperature.back() << " K\n";
      EXPECT_LT(worst, 1e-6);
      // GUARD: the gas cooled by most of its excess over the ambient.
      EXPECT_LT(t->temperature.back() - Ta, 0.5 * (T0 - Ta));
      // The pressure follows the temperature at constant moles.
      EXPECT_NEAR(t->pressure.back(), N * kR * t->temperature.back() / V, 1e-6 * P0);
    }
  }

  // ---------------------------------------------------------------------------
  // 5. Adiabatic filling
  // ---------------------------------------------------------------------------
  //
  // From the energy balance with a single inflow at T_in:
  //   d(N c_v T)/dt = F c_p T_in  =>  N c_v T = N_0 c_v T_0 + c_p T_in (N - N_0),
  // i.e. T = gamma T_in + (N_0 / N)(T_0 - gamma T_in), gamma = c_p / c_v: the
  // classic gamma T_in limit of filling an evacuated vessel. Exact; bound
  // 1e-6 relative (integrator tolerance 1e-9).
  TEST(DynamicTank, AdiabaticFillingApproachesGammaTimesTheInletTemperature) {
    UserData data;
    BedOptions o;
    o.isothermal = false;
    BuildBed(data, o);
    const double V = 0.01, P0 = 1e3, T_in = kT;
    Tank& tank = NewTank(data, V, P0, kT, {{"Nitrogen", 1.0}});
    tank.thermal_ = TankThermalMode::ADIABATIC;
    const size_t tank_id = tank.id_;
    auto feed = std::make_unique<TankMakeUp>();
    feed->name_ = "Fill";
    feed->molar_flow_ = 0.01;
    feed->temperature_ = T_in;
    feed->mole_fraction_["Nitrogen"] = 1.0;
    tank.streams_.push_back(std::move(feed));
    data.cycle.AddStep("Fill");
    Step& s = data.cycle.GetStep("Fill");
    s.step_time_ = 800.0;
    Close(s, BoundaryLocation::BOTTOM);
    Close(s, BoundaryLocation::TOP);
    data.cycle.number_of_cycles_ = 1;
    CoupledSolver solver{data};
    Run(data, solver);
    const auto& log = solver.GetLog();
    const TankLogData* t = log.GetTankLog(0, 0, tank_id);
    ASSERT_NE(t, nullptr);
    const double cp = data.fluid.MolarHeatCapacity("Nitrogen");
    const double gamma = cp / (cp - kR);
    const double N0 = P0 * V / (kR * kT);
    double worst = 0.0;
    for (size_t k = 0; k < t->temperature.size(); ++k) {
      const double N = t->total_moles[k];
      const double T_exact = gamma * T_in + (N0 / N) * (kT - gamma * T_in);
      worst = std::max(worst, std::abs(t->temperature[k] - T_exact) / T_exact);
    }
    std::cout << "[Tank fill] worst relative |T - T_exact| " << worst << ", T_end " << t->temperature.back()
              << " K (gamma T_in = " << gamma * T_in << " K), N / N0 " << t->total_moles.back() / N0 << "\n";
    EXPECT_LT(worst, 1e-6);
    // GUARD: the tank heated by over 100 K, to within 0.1 % of gamma T_in.
    EXPECT_GT(t->temperature.back() - kT, 100.0);
    EXPECT_NEAR(t->temperature.back() / (gamma * T_in), 1.0, 1e-3);
  }

  // ---------------------------------------------------------------------------
  // 7. CSS with a tank
  // ---------------------------------------------------------------------------
  //
  // An inert, isothermal bed is fed at the bottom and sends its outflow into a
  // large tank (time constant N / F ~ 100 cycles), which vents through a valve;
  // a blowdown step empties the bed. The bed has no memory beyond a cycle
  // except through the tank's pressure, which moves by ~1 % per cycle, so the
  // bed-only tests (its balance and its end-of-cycle state) pass within a few
  // cycles. The tank's composition and inventory still drift; the CSS rule,
  // which includes them, must not stop the run.
  TEST(DynamicTank, CSSWaitsForTheTank) {
    UserData data;
    BedOptions o;
    o.components = {"Nitrogen", "Oxygen"};
    o.y0 = {{"Nitrogen", 0.79}, {"Oxygen", 0.21}};
    BuildBed(data, o);
    Tank& tank = NewTank(data, 1.0, 1e5, kT, {{"Nitrogen", 1.0}});
    const size_t tank_id = tank.id_;
    auto vent = std::make_unique<TankValve>();
    vent->name_ = "Product";
    vent->cv_ = 0.5;
    vent->downstream_pressure_ = 1e5;
    tank.streams_.push_back(std::move(vent));

    data.cycle.AddStep("Feed");
    {
      Step& s = data.cycle.GetStep("Feed");
      s.step_time_ = 20.0;
      auto feed = PortFlow(0.02, 0, "Feed");
      feed->mole_fraction["Nitrogen"] = 0.79;
      feed->mole_fraction["Oxygen"] = 0.21;
      AddStream(s, BoundaryLocation::BOTTOM, BoundaryDirection::INLET, std::move(feed));
      AddStream(s, BoundaryLocation::TOP, BoundaryDirection::OUTLET, PortValve(0.2, tank_id, "To tank"));
    }
    data.cycle.AddStep("Blowdown");
    {
      Step& s = data.cycle.GetStep("Blowdown");
      s.step_time_ = 20.0;
      auto bc = std::make_unique<PressureBoundaryCondition>();
      bc->direction_ = BoundaryDirection::OUTLET;
      bc->location_ = BoundaryLocation::BOTTOM;
      bc->pressure_ = Constant(1e5);
      s.bottom_boundary_condition_ = std::move(bc);
      Close(s, BoundaryLocation::TOP);
    }
    data.cycle.number_of_cycles_ = 12;
    auto& css = data.settings.css;
    css.method_ = CSSMethod::MASS_BALANCE_COMPONENT;
    css.min_cycles_ = 3;
    css.consecutive_cycles_ = 2;
    data.settings.max_step = 0.5;

    CoupledSolver solver{data};
    Run(data, solver);
    const auto& log = solver.GetLog();
    EXPECT_FALSE(solver.ReachedCSS()) << "CSS declared while the tank still drifts";

    // At the last cycles the bed alone passes, the tank does not.
    const int last = *log.sys_data.cycle_index_set_.rbegin();
    const double tol_state = css.state_tolerance_;
    const double tol_balance = css.mass_balance_tolerance_;
    int bed_passing = 0;
    for (int c = last - 1; c <= last; ++c) {
      const double bed_state = CyclicSteadyState::BedStateChange(log, c);
      double bed_balance = 0.0;
      for (const auto& comp : o.components) {
        double in = 0.0, out = 0.0;
        for (int st : log.sys_data.step_index_set_)
          for (BoundaryLocation loc : {BoundaryLocation::BOTTOM, BoundaryLocation::TOP}) {
            const double n = GetTotalComponentMolesIn(log, c, st, loc, comp);
            (n > 0.0 ? in : out) += std::abs(n);
          }
        bed_balance = std::max(bed_balance, std::abs(in - out) / in);
      }
      const double tank_state = CyclicSteadyState::TankStateChange(log, c);
      std::cout << "[Tank CSS] cycle " << c + 1 << ": bed state change " << bed_state << ", bed balance "
                << bed_balance << ", tank state change " << tank_state << " (tolerances " << tol_state << ", "
                << tol_balance << ")\n";
      if (bed_state >= 0.0 && bed_state <= tol_state && bed_balance <= tol_balance) ++bed_passing;
      // The tank's drift is what keeps the rule from passing.
      EXPECT_GT(tank_state, 10.0 * tol_state);
    }
    // GUARD: the bed-only rule would have stopped the run.
    EXPECT_EQ(bed_passing, 2);
  }

}  // namespace adsim
