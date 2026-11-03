/*
 * Skarstrom - pressure and temperature swing adsorption simulator
 * Copyright (c) 2024-2026 Sam Affleck and contributors
 *
 * Licensed under the terms in the LICENSE file at the repository root.
 */

// Analytical validation of the UNIT OPERATIONS between and around beds:
// constant-speed evacuation, two-bed pressure equalisation and a buffer-tank
// fill.
//
// Both reduce to lumped ODEs because the bed is isothermal (isothermal mode),
// the pressure drop is negligible and, with linear isotherms and fast LDF, the
// inventory is linear in the pressure. Per unit bed volume a bed then holds
//
//     n(P) = a P,   a = eps / (R T) + rho_b K   [mol/(m3 Pa)],
//
// for any composition when every component has the same Henry slope K (the
// tests use the same K for both gases so that composition cannot enter).
//
// 1. Evacuation at a constant volumetric speed S (a volumetric-flow outlet,
//    the face concentration P/(R T)): V a dP/dt = -S P / (R T), so
//
//        P = P0 exp(-t / tau_S),   tau_S = V a R T / S = (eps V / S) (1 + rho_b K R T / eps).
//
//    Adsorption stretches the time constant by exactly the Henry retardation
//    factor; the test measures both the inert and the adsorbing decay rate.
//
// 2. Two-bed pressure equalisation through a linear valve,
//    n_dot = cv (P_hi - P_lo) / 1e5. The single-bed model represents the two
//    beds as two connected steps of the same bed (ED: depressurising, EU:
//    repressurising), each driven by its partner's pressure from the previous
//    cycle; at the cyclic steady state the two histories are each other's
//    partner exactly. Then V a dP_ED/dt = -n_dot = -V a dP_EU/dt, so the
//    total moles are conserved and the difference decays exponentially:
//
//        P_ED + P_EU = P_H + P_L                      (mole balance),
//        P_ED - P_EU = (P_H - P_L) exp(-t / tau_E),   tau_E = 1e5 V a / (2 cv),
//
//    and the final (equalised) pressure is (P_H + P_L)/2. The 2 in tau_E is the
//    signature of two beds exchanging gas, as against one bed venting to a
//    fixed pressure.
//    With a Langmuir isotherm the pressures no longer add up; the conserved
//    quantity is the inventory, n(P_ED) + n(P_EU) = n(P_H) + n(P_L).
//
// 3. Buffer-tank fill. A "Vent" step discharges the bed through a valve and a
//    later "Fill" step (other end closed) takes p % of the vented moles from a
//    buffer tank, spread evenly over the step. By the mole balance of the two
//    rigid isothermal beds the fill raises the pressure by
//
//        P_fill,end - P_vent,end = (p / 100) (P_0 - P_vent,end),
//
//    whatever the isotherm slope (the same capacity a appears on both sides),
//    and the fill pressure rises LINEARLY in time (constant molar inflow).
//
// Neglected terms: LDF lag O(1/(k tau)) (k = 20 1/s, tau >= 25 s: < 0.2 %),
// and the pressure drop across the bed (a few Pa on >= 5 kPa).

#include "validation_helpers.h"

#include <algorithm>

namespace adsim {

  using namespace adsim::validation;

  namespace {

    constexpr double kT = 298.0;
    constexpr double kEps = 0.4;
    constexpr double kRhoS = 1000.0;
    constexpr double kL = 1.0, kD = 0.1;
    constexpr double kLdf = 20.0;      // [1/s]
    constexpr double kHenry = 1.0e-6;  // [mol/(kg Pa)], both gases

    double Volume() { return 0.25 * M_PI * kD * kD * kL; }
    double RhoB() { return (1.0 - kEps) * kRhoS; }
    double Capacity(bool adsorbing) { return kEps / (kR * kT) + (adsorbing ? RhoB() * kHenry : 0.0); }

    // A small isothermal bed of N2/O2 (y_O2 = 0.21 everywhere), optionally with
    // the same Henry isotherm on both gases.
    // Langmuir variant: the same single-site Langmuir on both gases, so the
    // total loading m K P / (1 + K P) is again composition-independent.
    constexpr double kLangmuirM = 0.5;    // [mol/kg]
    constexpr double kLangmuirK = 5e-6;   // [1/Pa]: K P = 0.5 .. 1.5
    // Inventory per bed volume [mol/m3] with the Langmuir isotherm.
    double LangmuirInventory(double P) {
      return kEps * P / (kR * kT) + RhoB() * kLangmuirM * kLangmuirK * P / (1.0 + kLangmuirK * P);
    }

    void BuildBed(UserData& data, bool adsorbing, double P0, bool langmuir = false) {
      data.bed.wall_.thickness_ = 0.005;
      data.fluid.AddComponent("Nitrogen");
      data.fluid.AddComponent("Oxygen");
      auto ads = data.adsorbent_directory.Add("Packing");
      auto lid = data.bed.layers_.Add("Packing");
      auto& layer = data.bed.layers_.Get(lid);
      layer.momentum_model_ = MomentumModel::KarmenCozeny;
      layer.domain_.length = kL;
      layer.domain_.n_cells = 10;
      layer.eb_ = kEps;
      layer.particle_diameter_ = 0.005;  // coarse: negligible pressure drop
      layer.composition_[ads] = 1.0;
      std::get<CylinderBed>(data.bed.shape_).diameter = kD;
      auto& adsorbent = data.adsorbent_directory.Get(ads);
      adsorbent.ep_ = 0.0;
      adsorbent.solid_density_ = kRhoS;
      adsorbent.heat_capacity_ = 900.0;
      if (adsorbing && langmuir) {
        auto iso = std::make_unique<ExtendedLangmuirIsotherm>();
        for (const char* c : {"Nitrogen", "Oxygen"}) {
          iso->params_[c].m1 = kLangmuirM;
          iso->params_[c].b1 = kLangmuirK;
          iso->params_[c].q1 = 0.0;
        }
        adsorbent.adsorption_source_.AddIsotherm({"Nitrogen", "Oxygen"}, std::move(iso));
      } else if (adsorbing) {
        auto iso = std::make_unique<HenryIsotherm>();
        for (const char* c : {"Nitrogen", "Oxygen"}) {
          iso->params_[c].k = kHenry;
          iso->params_[c].Q = 0.0;
        }
        adsorbent.adsorption_source_.AddIsotherm({"Nitrogen", "Oxygen"}, std::move(iso));
      }
      if (adsorbing) {
        for (const char* c : {"Nitrogen", "Oxygen"}) {
          auto ldf = std::make_unique<LinearDrivingForce>();
          auto mtc = std::make_unique<MassTransferCoefficientConstant>();
          mtc->constant_mtc_ = kLdf;
          ldf->mtc_model_ = std::move(mtc);
          adsorbent.adsorption_source_.mass_transfer_models_[c] = std::move(ldf);
          adsorbent.adsorption_source_.heat_of_adsorption_[c] = 0.0;
        }
      }
      data.bed.equations_.energy_balance_mode_ = EnergyBalanceMode::ISOTHERMAL;
      data.ic.pressure = P0;
      data.ic.temperature = kT;
      data.ic.wall_temperature = kT;
      data.ic.velocity = 0.0;
      data.ic.mole_fraction["Nitrogen"] = 0.79;
      data.ic.mole_fraction["Oxygen"] = 0.21;
      data.settings.reltol = 1e-9;
      data.settings.max_step = 1.0;
    }

    template <class Map>
    void SetAir(Map& y) {
      y["Nitrogen"] = 0.79;
      y["Oxygen"] = 0.21;
    }

    std::unique_ptr<ConstantBoundaryValue> Constant(double v) {
      auto bv = std::make_unique<ConstantBoundaryValue>();
      bv->constant_value_ = v;
      return bv;
    }

    std::unique_ptr<LinearBoundaryValue> Ramp(double from, double to, double ramp_time) {
      auto bv = std::make_unique<LinearBoundaryValue>();
      bv->start_value_ = from;
      bv->end_value_ = to;
      bv->ramp_time_ = ramp_time;
      return bv;
    }

    void CloseBottom(Step& step) {
      auto closed = std::make_unique<ClosedBoundaryCondition>();
      closed->location_ = BoundaryLocation::BOTTOM;
      step.bottom_boundary_condition_ = std::move(closed);
    }

    // Volume-mean pressure per logged snapshot of one step.
    std::pair<std::vector<double>, std::vector<double>> MeanPressure(const EquationLog& log,
                                                                     int cycle, int step) {
      const auto P = log.GetVariable(cycle, step, VariableType::PRESSURE);
      std::vector<double> mean;
      for (const auto& snap : P.values) mean.push_back(InteriorMean(snap, log.GetDomain().volume));
      return {P.times, mean};
    }

    // Least-squares slope of ln(f) against t over the samples in [t_a, t_b].
    double LogSlope(const std::vector<double>& t, const std::vector<double>& f,
                    double t_a, double t_b) {
      double n = 0, st = 0, sy = 0, stt = 0, sty = 0;
      for (size_t i = 0; i < t.size(); ++i) {
        if (t[i] < t_a || t[i] > t_b || f[i] <= 0.0) continue;
        const double y = std::log(f[i]);
        n += 1; st += t[i]; sy += y; stt += t[i] * t[i]; sty += t[i] * y;
      }
      return (n * sty - st * sy) / (n * stt - st * st);
    }

    // -------------------------------------------------------------------------
    // 1. Constant-speed evacuation
    // -------------------------------------------------------------------------

    struct EvacuationResult { double tau_model, tau_fit, worst_rel; };

    EvacuationResult RunEvacuation(bool adsorbing) {
      const double P0 = 1.0e5;
      const double S = 6.0e-5;  // [m3/s]
      const double tau = Volume() * Capacity(adsorbing) * kR * kT / S;
      const double run_time = 3.0 * tau;

      UserData data;
      BuildBed(data, adsorbing, P0);
      data.settings.max_log_points_per_step_ = 300;
      data.cycle.AddStep("Evacuate");
      Step& step = data.cycle.GetStep("Evacuate");
      step.step_time_ = run_time;
      CloseBottom(step);
      auto stream = std::make_unique<BoundaryStream>();
      auto flow = std::make_unique<VolumetricFlowBoundaryValue>();
      static_cast<ConstantBoundaryValue*>(flow->profile_.get())->constant_value_ = S;
      stream->boundary_flow_ = std::move(flow);
      auto outlet = std::make_unique<FlowBoundaryCondition>();
      outlet->direction_ = BoundaryDirection::OUTLET;
      outlet->location_ = BoundaryLocation::TOP;
      outlet->streams_.push_back(std::move(stream));
      step.top_boundary_condition_ = std::move(outlet);
      data.cycle.number_of_cycles_ = 1;

      CoupledSolver solver{data};
      std::string error;
      EXPECT_TRUE(data.PerformChecks(error)) << error;
      EXPECT_GE(solver.IntegrateCycle(data.cycle, false), 0) << "Solver did not complete";
      const auto [t, P] = MeanPressure(solver.GetLog(), 0, 0);

      EvacuationResult r{tau, -1.0 / LogSlope(t, P, 0.1 * tau, run_time), 0.0};
      for (size_t i = 0; i < t.size(); ++i)
        r.worst_rel = std::max(r.worst_rel, std::abs(P[i] / (P0 * std::exp(-t[i] / tau)) - 1.0));
      return r;
    }

    // -------------------------------------------------------------------------
    // 2. Two-bed pressure equalisation (single-bed model, connected steps)
    // -------------------------------------------------------------------------

    constexpr double kPH = 3.0e5, kPL = 1.0e5;
    constexpr double kCvEq = 2.0e-3;   // [mol/s per bar]

    struct EqualisationRun {
      std::vector<double> t, P_ed, P_eu;
      double tau_model = 0.0;
      int cycles = 0;
    };

    EqualisationRun RunEqualisation(bool adsorbing, int n_cycles, bool langmuir = false) {
      // For the Langmuir bed, the linearised capacity at P_eq sets the step length.
      const double P_mid = 0.5 * (kPH + kPL);
      const double a = langmuir ? (LangmuirInventory(P_mid * 1.001) - LangmuirInventory(P_mid * 0.999))
                                    / (0.002 * P_mid)
                                : Capacity(adsorbing);
      const double tau_E = 1e5 * Volume() * a / (2.0 * kCvEq);
      const double t_eq = 5.0 * tau_E;       // ED and EU (equal, as connected steps must be)
      const double t_pb = 60.0;              // pressurise / blow down to a fixed pressure

      UserData data;
      BuildBed(data, adsorbing, 0.5 * (kPH + kPL), langmuir);
      data.settings.max_log_points_per_step_ = 200;

      // Pressurise to P_H through the top (a pressure boundary).
      data.cycle.AddStep("Pressurise");
      {
        Step& s = data.cycle.GetStep("Pressurise");
        s.step_time_ = t_pb;
        CloseBottom(s);
        auto bc = std::make_unique<PressureBoundaryCondition>();
        bc->direction_ = BoundaryDirection::INLET;
        bc->location_ = BoundaryLocation::TOP;
        bc->pressure_ = Ramp(0.5 * (kPH + kPL), kPH, 20.0);  // from the equalised pressure
        bc->temperature_ = kT;
        SetAir(bc->mole_fraction_);
        s.top_boundary_condition_ = std::move(bc);
      }
      // ED: vent through a valve into the partner bed.
      data.cycle.AddStep("ED");
      ValveBoundaryStream* ed_valve = nullptr;
      {
        Step& s = data.cycle.GetStep("ED");
        s.step_time_ = t_eq;
        CloseBottom(s);
        auto valve = std::make_unique<ValveBoundaryStream>();
        valve->boundary_value_ = Constant(kPL);   // first-cycle assumption
        valve->valve_.cv_ = kCvEq;
        ed_valve = valve.get();
        auto bc = std::make_unique<FlowBoundaryCondition>();
        bc->direction_ = BoundaryDirection::OUTLET;
        bc->location_ = BoundaryLocation::TOP;
        bc->streams_.push_back(std::move(valve));
        s.top_boundary_condition_ = std::move(bc);
      }
      // Blow down to P_L.
      data.cycle.AddStep("Blowdown");
      {
        Step& s = data.cycle.GetStep("Blowdown");
        s.step_time_ = t_pb;
        CloseBottom(s);
        auto bc = std::make_unique<PressureBoundaryCondition>();
        bc->direction_ = BoundaryDirection::OUTLET;
        bc->location_ = BoundaryLocation::TOP;
        bc->pressure_ = Ramp(0.5 * (kPH + kPL), kPL, 20.0);
        SetAir(bc->mole_fraction_);
        s.top_boundary_condition_ = std::move(bc);
      }
      // EU: fill through a valve from the partner bed.
      data.cycle.AddStep("EU");
      ValveBoundaryStream* eu_valve = nullptr;
      {
        Step& s = data.cycle.GetStep("EU");
        s.step_time_ = t_eq;
        CloseBottom(s);
        auto valve = std::make_unique<ValveBoundaryStream>();
        valve->boundary_value_ = Constant(kPH);   // first-cycle assumption
        valve->valve_.cv_ = kCvEq;
        valve->temperature_ = kT;
        SetAir(valve->mole_fraction);
        eu_valve = valve.get();
        auto bc = std::make_unique<FlowBoundaryCondition>();
        bc->direction_ = BoundaryDirection::INLET;
        bc->location_ = BoundaryLocation::TOP;
        bc->streams_.push_back(std::move(valve));
        s.top_boundary_condition_ = std::move(bc);
      }
      {
        auto a = MakeEndpoint(data.cycle.GetStep("ED"), *ed_valve);
        auto b = MakeEndpoint(data.cycle.GetStep("EU"), *eu_valve);
        EXPECT_TRUE(MakeMutualConnection(a, b));
      }
      data.cycle.number_of_cycles_ = n_cycles;

      CoupledSolver solver{data};
      std::string error;
      EXPECT_TRUE(data.PerformChecks(error)) << error;
      EXPECT_GE(solver.IntegrateCycle(data.cycle, false), 0) << "Solver did not complete";
      const auto& log = solver.GetLog();
      const int last = *log.sys_data.cycle_index_set_.rbegin();
      EqualisationRun r;
      // GUARD: the lagged coupling has converged (cyclic steady state): the
      // previous cycle ended the ED step at the same pressure.
      {
        const auto prev = MeanPressure(log, last - 1, 1).second;
        const auto curr = MeanPressure(log, last, 1).second;
        EXPECT_NEAR(prev.back(), curr.back(), 1e-6 * kPH) << "Not at the cyclic steady state";
      }
      r.tau_model = tau_E;
      r.cycles = last + 1;
      // The two steps are logged at their own instants: put EU on ED's clock
      // (time into the step).
      auto [t_ed, P_ed] = MeanPressure(log, last, 1);
      auto [t_eu, P_eu] = MeanPressure(log, last, 3);
      std::vector<double> tau_eu;
      for (double t : t_eu) tau_eu.push_back(t - t_eu.front());
      for (size_t i = 0; i < t_ed.size(); ++i) {
        r.t.push_back(t_ed[i] - t_ed.front());
        r.P_ed.push_back(P_ed[i]);
        r.P_eu.push_back(InterpolateSeries(tau_eu, P_eu, r.t.back()));
      }
      return r;
    }

  } // namespace

  TEST(UnitOperations, ConstantSpeedEvacuationDecaysExponentially) {
    const auto inert = RunEvacuation(false);
    const auto ads = RunEvacuation(true);
    const double R_H = Capacity(true) / Capacity(false);  // Henry retardation factor
    // GUARD: adsorption changes the time constant a lot.
    ASSERT_GT(R_H, 2.0);
    std::cout << "[Evacuation] inert tau " << inert.tau_fit << " s (model " << inert.tau_model
              << "), max rel dev " << inert.worst_rel << "; Henry tau " << ads.tau_fit
              << " s (model " << ads.tau_model << "), max rel dev " << ads.worst_rel << "\n";

    // Over three time constants (P falls 20x) the pressure stays on the
    // exponential, and the fitted time constants are the model's.
    EXPECT_LT(inert.worst_rel, 2e-3);
    EXPECT_LT(ads.worst_rel, 5e-3);
    EXPECT_NEAR(inert.tau_fit, inert.tau_model, 1e-3 * inert.tau_model);
    EXPECT_NEAR(ads.tau_fit, ads.tau_model, 2e-3 * ads.tau_model);
    // Adsorption stretches the decay by exactly the retardation factor.
    EXPECT_NEAR(ads.tau_fit / inert.tau_fit, R_H, 3e-3 * R_H);
  }

  TEST(UnitOperations, TwoBedEqualisationConservesMolesAndDecaysAtTwiceTheRate) {
    for (bool adsorbing : {false, true}) {
      const auto r = RunEqualisation(adsorbing, 30);
      ASSERT_GT(r.t.size(), 10u);
      const double P_eq = 0.5 * (kPH + kPL);

      // At the start of the step the beds hold P_H and P_L.
      EXPECT_NEAR(r.P_ed.front(), kPH, 1e-3 * kPH);
      EXPECT_NEAR(r.P_eu.front(), kPL, 1e-3 * kPL);

      double worst_sum = 0.0, worst_diff = 0.0;
      for (size_t i = 0; i < r.t.size(); ++i) {
        worst_sum = std::max(worst_sum, std::abs(r.P_ed[i] + r.P_eu[i] - (kPH + kPL)) / P_eq);
        const double diff_model = (kPH - kPL) * std::exp(-r.t[i] / r.tau_model);
        worst_diff = std::max(worst_diff, std::abs(r.P_ed[i] - r.P_eu[i] - diff_model) / (kPH - kPL));
      }
      const double tau_fit = -1.0 / LogSlope(r.t, [&] {
        std::vector<double> d;
        for (size_t i = 0; i < r.t.size(); ++i) d.push_back(r.P_ed[i] - r.P_eu[i]);
        return d;
      }(), 0.0, 3.0 * r.tau_model);
      std::cout << "[Equalisation] adsorbing=" << adsorbing << " cycles=" << r.cycles
                << " P_end ED/EU " << r.P_ed.back() << " / " << r.P_eu.back() << " (P_eq " << P_eq
                << ") tau fit " << tau_fit << " s (model " << r.tau_model << ")"
                << " max|sum|/P_eq " << worst_sum << " max|diff err|/dP0 " << worst_diff << "\n";

      // The step lasts 5 tau_E, so the beds end (P_H - P_L) exp(-5) / 2 either
      // side of P_eq.
      const double half_gap = 0.5 * (kPH - kPL) * std::exp(-r.t.back() / r.tau_model);
      EXPECT_NEAR(r.P_ed.back(), P_eq + half_gap, 1e-3 * P_eq) << "adsorbing=" << adsorbing;
      EXPECT_NEAR(r.P_eu.back(), P_eq - half_gap, 1e-3 * P_eq) << "adsorbing=" << adsorbing;
      EXPECT_LT(worst_sum, 2e-3) << "Total moles of the two beds not conserved (adsorbing=" << adsorbing << ")";
      EXPECT_LT(worst_diff, 3e-3) << "Pressure difference does not decay as exp(-t/tau_E)";
      // The two beds exchange gas: tau_E, not the single-bed 2 tau_E.
      EXPECT_NEAR(tau_fit, r.tau_model, 5e-3 * r.tau_model) << "adsorbing=" << adsorbing;
    }
  }

  // With a NONLINEAR (Langmuir) isotherm the pressures no longer add up: what
  // the two beds conserve is their inventory n(P) = eps P/(R T) + rho_b q*(P),
  // so at the end of the step n(P_ED) + n(P_EU) = n(P_H) + n(P_L), while
  // P_ED + P_EU falls short of P_H + P_L (more gas is held per pascal at low
  // pressure). The test asserts the first and, as a guard, that the second is
  // violated by far more than the first's tolerance.
  TEST(UnitOperations, TwoBedEqualisationWithLangmuirConservesTheInventory) {
    const auto r = RunEqualisation(true, 30, /*langmuir=*/true);
    ASSERT_GT(r.t.size(), 10u);
    EXPECT_NEAR(r.P_ed.front(), kPH, 1e-3 * kPH);
    EXPECT_NEAR(r.P_eu.front(), kPL, 1e-3 * kPL);
    const double n0 = LangmuirInventory(kPH) + LangmuirInventory(kPL);
    double worst_inventory = 0.0;
    for (size_t i = 0; i < r.t.size(); ++i)
      worst_inventory = std::max(worst_inventory,
        std::abs(LangmuirInventory(r.P_ed[i]) + LangmuirInventory(r.P_eu[i]) - n0) / n0);
    const double pressure_sum_defect = std::abs(r.P_ed.back() + r.P_eu.back() - (kPH + kPL)) / (kPH + kPL);
    std::cout << "[Equalisation] Langmuir P_end ED/EU " << r.P_ed.back() << " / " << r.P_eu.back()
              << " max inventory defect " << worst_inventory
              << " pressure-sum defect " << pressure_sum_defect << "\n";
    // Tolerance: the LDF lag (k = 20 1/s against tau ~ 100 s), as for Henry.
    EXPECT_LT(worst_inventory, 2e-3) << "The two beds' inventory is not conserved";
    // GUARD: the isotherm's curvature is visible in the pressures.
    EXPECT_GT(pressure_sum_defect, 10.0 * 2e-3);
  }

  TEST(UnitOperations, BufferTankFillDeliversTheVentedMoles) {
    const double P0 = 3.0e5, P_down = 1.0e5, cv = 2.0e-3, pct = 60.0;
    for (bool adsorbing : {false, true}) {
      const double tau = 1e5 * Volume() * Capacity(adsorbing) / cv;  // vent time constant
      const double t_step = 2.0 * tau;

      UserData data;
      BuildBed(data, adsorbing, P0);
      data.settings.max_log_points_per_step_ = 400;
      data.cycle.AddStep("Vent");
      {
        Step& s = data.cycle.GetStep("Vent");
        s.step_time_ = t_step;
        CloseBottom(s);
        auto valve = std::make_unique<ValveBoundaryStream>();
        valve->name_ = "Vent valve";
        valve->boundary_value_ = Constant(P_down);
        valve->valve_.cv_ = cv;
        auto bc = std::make_unique<FlowBoundaryCondition>();
        bc->direction_ = BoundaryDirection::OUTLET;
        bc->location_ = BoundaryLocation::TOP;
        bc->streams_.push_back(std::move(valve));
        s.top_boundary_condition_ = std::move(bc);
      }
      data.cycle.AddStep("Fill");
      {
        Step& s = data.cycle.GetStep("Fill");
        s.step_time_ = t_step;
        CloseBottom(s);
        auto in = std::make_unique<BoundaryStream>();
        in->name_ = "Tank";
        auto flow = std::make_unique<MolarFlowBoundaryValue>();
        static_cast<ConstantBoundaryValue*>(flow->profile_.get())->constant_value_ = 1e-4;  // unused: tank is ready
        in->boundary_flow_ = std::move(flow);
        SetAir(in->mole_fraction);
        in->temperature_ = kT;
        in->buffer_tank_.Activate("Vent", "Vent valve");
        in->buffer_tank_.percentage_ = pct;
        auto bc = std::make_unique<FlowBoundaryCondition>();
        bc->direction_ = BoundaryDirection::INLET;
        bc->location_ = BoundaryLocation::TOP;
        bc->streams_.push_back(std::move(in));
        s.top_boundary_condition_ = std::move(bc);
      }
      data.cycle.number_of_cycles_ = 1;

      CoupledSolver solver{data};
      std::string error;
      ASSERT_TRUE(data.PerformChecks(error)) << error;
      ASSERT_GE(solver.IntegrateCycle(data.cycle, false), 0) << "Solver did not complete";
      const auto& log = solver.GetLog();
      const auto [t_v, P_v] = MeanPressure(log, 0, 0);
      const auto [t_f, P_f] = MeanPressure(log, 0, 1);

      // GUARD: a substantial vent.
      const double vented = P0 - P_v.back();
      ASSERT_GT(vented, 0.5 * (P0 - P_down));
      const double rise_model = pct / 100.0 * vented;
      const double rise = P_f.back() - P_f.front();
      // Linear in time: compare with the straight line at mid-step.
      const double mid = InterpolateSeries(t_f, P_f, t_f.front() + 0.5 * t_step);
      std::cout << "[BufferTank] adsorbing=" << adsorbing << " vented " << vented << " Pa, fill rise "
                << rise << " Pa (model " << rise_model << "), mid-step " << mid - P_f.front()
                << " Pa (linear " << 0.5 * rise_model << ")\n";
      EXPECT_NEAR(P_f.front(), P_v.back(), 1e-6 * P0);
      EXPECT_NEAR(rise, rise_model, 3e-3 * rise_model)
        << "The fill did not deliver " << pct << " % of the vented moles (adsorbing=" << adsorbing << ")";
      EXPECT_NEAR(mid - P_f.front(), 0.5 * rise_model, 3e-3 * rise_model)
        << "The fill pressure does not rise linearly (constant tank flow)";
    }
  }

} // namespace adsim
