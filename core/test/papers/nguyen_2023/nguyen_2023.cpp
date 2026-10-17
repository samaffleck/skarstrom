// STL includes
#define _USE_MATH_DEFINES
#include <cmath>
#include <filesystem>

// AdSim includes
#include "../../../src/pch.h"
#include "adsim/solver/optimisation.h"
#include "adsim/solver/coupled_solver.h"
#include "adsim/solver/kpis.h"
#include "adsim/math/constants.h"
#include "adsim/equations/isotherm.h"
#include "adsim/boundary_conditions/boundary_conditions.h"
#include "adsim/boundary_conditions/boundary_stream.h"
#include "adsim/boundary_conditions/connection.h"
#include "adsim/boundary_conditions/step.h"
#include "adsim/equations/dispersion_model.h"
#include "adsim/objects/compressor.h"
#include "adsim/equations/heat_transfer_model.h"

// GTest includes
#include <gtest/gtest.h>

#include "nguyen_2023.h"
#include "../../fixtures/validation_metrics.h"
#include <limits>
#include <map>

namespace adsim {

  // ==========================================================================
  //  Quantitative comparison with Table 2 of Nguyen, Shimizu & Rajendran,
  //  Chem. Eng. J. 452 (2023) 139550, doi:10.1016/j.cej.2022.139550:
  //    * experimental CO2 purity / recovery / productivity (validation), and
  //    * the paper's own detailed-model purity / recovery / productivity
  //      (code-to-code check).
  //  Productivity follows Sect. 2.3: mol CO2 recovered in the evacuation step
  //  per m3 of ADSORBENT per second; the adsorbent volume is m_ads / rho_p
  //  (= (1 - eps) A L of the deck). Energy (Table 2 "En") is not compared: the
  //  deck uses pressure boundary conditions, not a vacuum-pump model.
  //  Tolerances (abs. percentage points for Pu/Re, relative for Pr) are set from
  //  the measured error, rounded up; see kNguyen2023Tol below.
  // ==========================================================================
  namespace {
    struct Nguyen2023Tol {
      double pu_exp, re_exp, pr_exp, pu_model, re_model, pr_model;
    };
    const std::map<std::string, Nguyen2023Tol> kNguyen2023Tol = {
      // measured: 2.0005 pp, 2.43014 pp, 0.0655223 rel | 1.0005 pp, 0.76986 pp, 0.0358563 rel
      { "ExperimentI",    { 2.4, 2.8, 0.076, 1.4, 1.1, 0.046 } },
      // measured: 1.85656 pp, 1.23426 pp, 0.0939452 rel | 0.956562 pp, 0.634264 pp, 0.127503 rel
      { "ExperimentII",   { 2.2, 1.6, 0.11, 1.3, 0.94, 0.15 } },
      // measured: 1.70017 pp, 2.49223 pp, 0.0847345 rel | 1.80017 pp, 0.892233 pp, 0.138574 rel
      { "ExperimentIII",  { 2.1, 2.9, 0.098, 2.2, 1.2, 0.16 } },
      // measured: 3.1767 pp, 2.6597 pp, 0.070217 rel | 2.1767 pp, 0.859704 pp, 0.1446 rel
      { "ExperimentIV",   { 3.7, 3.1, 0.081, 2.6, 1.2, 0.17 } },
      // measured: 4.8734 pp, 2.57905 pp, 0.0492303 rel | 2.0734 pp, 1.37905 pp, 0.135664 rel
      { "ExperimentV",    { 5.7, 3, 0.06, 2.4, 1.7, 0.16 } },
      // measured: 0.934884 pp, 6.75405 pp, 0.0659414 rel | 0.434884 pp, 5.05405 pp, 0.000777056 rel
      { "ExperimentVI",   { 1.3, 7.8, 0.076, 0.74, 5.9, 0.011 } },
      // measured: 0.70128 pp, 7.74829 pp, 0.139047 rel | 1.90128 pp, 8.04829 pp, 0.139047 rel
      { "ExperimentVII",  { 1.1, 9, 0.16, 2.3, 9.3, 0.16 } },
      // measured: 2.6684 pp, 9.21778 pp, 0.0790468 rel | 3.5684 pp, 7.11778 pp, 0.0214872 rel
      { "ExperimentVIII", { 3.1, 11, 0.091, 4.2, 8.2, 0.032 } },
      // measured: 1.32558 pp, 9.75953 pp, 0.115661 rel | 3.95558 pp, 8.65953 pp, 0.0359711 rel
      { "ExperimentIX",   { 1.7, 12, 0.14, 4.6, 10, 0.046 } },
      // measured: 0.0701674 pp, 8.26968 pp, 0.0171134 rel | 0.0701674 pp, 7.06968 pp, 0.173592 rel
      { "ExperimentX",    { 0.38, 9.6, 0.028, 0.38, 8.2, 0.2 } },
    };

    void ValidateNguyen2023(const std::string& experiment_name, const UserData& data,
                            const KPIParams& kpi, double purity, double recovery) {
      using namespace validation;
      // "Nguyen2023_Exp_VII" -> "ExperimentVII"
      const std::string id = "Experiment" + experiment_name.substr(experiment_name.rfind('_') + 1);
      const auto& tol = kNguyen2023Tol.at(id);
      const auto& log = kpi.eq_log_;
      const double t_cycle = log.sys_data.total_cycle_time;
      const double V_ads = 36.0e-3 / 570.0;  // m_ads / rho_p (Table 1)
      const double pr = CalculateComponentMolesInProduct(kpi) / (V_ads * t_cycle);

      CaseReport report("nguyen_2023", id);
      report.CheckAbs("CO2 purity vs experiment", purity, Reference("nguyen_2023", id, "CO2_purity_exp"), tol.pu_exp);
      report.CheckAbs("CO2 recovery vs experiment", recovery, Reference("nguyen_2023", id, "CO2_recovery_exp"), tol.re_exp);
      report.CheckRel("CO2 productivity vs experiment", pr, Reference("nguyen_2023", id, "CO2_productivity_exp"), tol.pr_exp);
      report.CheckAbs("CO2 purity vs paper model", purity, Reference("nguyen_2023", id, "CO2_purity_model"), tol.pu_model);
      report.CheckAbs("CO2 recovery vs paper model", recovery, Reference("nguyen_2023", id, "CO2_recovery_model"), tol.re_model);
      report.CheckRel("CO2 productivity vs paper model", pr, Reference("nguyen_2023", id, "CO2_productivity_model"), tol.pr_model);
      if (id == "ExperimentVII") {
        report.Info("CO2 purity vs repeat experiment VII.b", purity,
                    Reference("nguyen_2023", id, "CO2_purity_exp_b"), "abs_error", "repeat run of VII");
        report.Info("CO2 recovery vs repeat experiment VII.b", recovery,
                    Reference("nguyen_2023", id, "CO2_recovery_exp_b"), "abs_error", "repeat run of VII");
      }
      (void)data;
    }
  } // namespace

  // Helper function to set up and run Nguyen 2023 experiment cycle
  void RunNguyen2023Experiment(
    Nguyen2023Test& test_fixture,
    const std::string& experiment_name, 
    const std::string& title,
    const std::string& description,
    double feed_time,
    double interstitial_feed_vel,
    double bd_pressure
  ) {
    test_fixture.data.info.title = title;
    test_fixture.data.info.description = description;
    test_fixture.data.info.author = "Sam Affleck";
    
    const std::string CO2 = "CarbonDioxide";
    const std::string N2 = "Nitrogen";
    
    test_fixture.SetUpNguyen2023();

    // Cycle times
    // Paper terminates blowdown with a pressure trigger at P_I, using an exponential
    // decay P = P_H * exp(-lambda_BD * t) with lambda_BD = 0.153 /s (Table 1). The
    // trigger duration therefore depends on P_I: t = ln(P_H / P_I) / lambda_BD, rather
    // than a fixed 25 s. Match that so the blowdown is P_I-dependent like the paper's.
    const double lambda_bd = 0.153;                                  // [1/s], Table 1
    const double P_high = 1.0e5;                                     // [Pa]
    const double blowdown_time = std::log(P_high / bd_pressure) / lambda_bd; // [s]
    const double evac_time = 70; // [s]
    const double press_time = 29; // [s]

    const std::string feed = "Feed";
    const std::string blowdown = "Blowdown";
    const std::string evacuation = "Evacuation";
    const std::string press = "Pressurisation";

    test_fixture.AddStep(feed, feed_time);
    test_fixture.AddStep(blowdown, blowdown_time);
    test_fixture.AddStep(evacuation, evac_time);
    test_fixture.AddStep(press, press_time);

    const double feed_vel = interstitial_feed_vel * 0.4;
    // Feed temperature: not stated in the paper or SI text; 293 K assumed (same
    // as the assumed ambient / initial temperature, see nguyen_2023.h).
    const double feed_temp = 293;
    const double ambient_pressure = 1.0e5;
    const double evac_pressure = 0.031e5;
    const std::unordered_map<std::string, double> inlet_mf = {
      {CO2, 0.15 },
      {N2, 0.85 }
    };

    // FEED STEP
    Step& step_feed = test_fixture.data.cycle.GetStep(feed);
    auto feed_inlet_stream = test_fixture.CreateNoValveInletStream(
      test_fixture.CreateVelocityBC(feed_vel),
      inlet_mf,
      feed_temp
    );
    feed_inlet_stream->name_ = "Feed inlet";
    
    auto feed_flow_bc = std::make_unique<FlowBoundaryCondition>();
    feed_flow_bc->direction_ = BoundaryDirection::INLET;
    feed_flow_bc->location_ = BoundaryLocation::BOTTOM;
    feed_flow_bc->streams_.push_back(std::move(feed_inlet_stream));
    step_feed.bottom_boundary_condition_ = std::move(feed_flow_bc);

    auto feed_pressure_out = std::make_unique<PressureBoundaryCondition>();
    feed_pressure_out->direction_ = BoundaryDirection::OUTLET;
    feed_pressure_out->location_ = BoundaryLocation::TOP;
    feed_pressure_out->pressure_ = test_fixture.CreateConstantBoundaryValue(ambient_pressure);
    feed_pressure_out->mole_fraction_[CO2] = 0.0;
    feed_pressure_out->mole_fraction_[N2] = 1.0;
    step_feed.top_boundary_condition_ = std::move(feed_pressure_out);

    // PRESSURISATION
    auto rep_stream = test_fixture.CreatePressureInletBoundary(
      test_fixture.CreateLinearBoundaryValue(evac_pressure, ambient_pressure, press_time),
      inlet_mf,
      feed_temp
    );
    rep_stream->location_ = BoundaryLocation::BOTTOM;
    rep_stream->direction_ = BoundaryDirection::INLET;

    // BLOWDOWN
    Step& step_blowdown = test_fixture.data.cycle.GetStep(blowdown);
    auto blowdown_pressure = std::make_unique<PressureBoundaryCondition>();
    blowdown_pressure->direction_ = BoundaryDirection::OUTLET;
    blowdown_pressure->location_ = BoundaryLocation::TOP;
    blowdown_pressure->pressure_ = test_fixture.CreateExponentialBoundaryValue(bd_pressure, blowdown_time);
    blowdown_pressure->mole_fraction_[CO2] = 0.0;
    blowdown_pressure->mole_fraction_[N2] = 1.0;
    step_blowdown.top_boundary_condition_ = std::move(blowdown_pressure);
    
    auto blowdown_closed = std::make_unique<ClosedBoundaryCondition>();
    blowdown_closed->location_ = BoundaryLocation::BOTTOM;
    step_blowdown.bottom_boundary_condition_ = std::move(blowdown_closed);

    // EVACUATION
    Step& step_evac = test_fixture.data.cycle.GetStep(evacuation);
    auto evac_closed_top = std::make_unique<ClosedBoundaryCondition>();
    evac_closed_top->location_ = BoundaryLocation::TOP;
    step_evac.top_boundary_condition_ = std::move(evac_closed_top);
    
    auto evac_pressure_bc = std::make_unique<PressureBoundaryCondition>();
    evac_pressure_bc->direction_ = BoundaryDirection::OUTLET;
    evac_pressure_bc->location_ = BoundaryLocation::BOTTOM;
    evac_pressure_bc->pressure_ = test_fixture.CreateExponentialBoundaryValue(evac_pressure, evac_time);
    evac_pressure_bc->mole_fraction_[CO2] = 0.0;
    evac_pressure_bc->mole_fraction_[N2] = 1.0;
    step_evac.bottom_boundary_condition_ = std::move(evac_pressure_bc);
    
    // PRESSURISATION
    Step& step_press = test_fixture.data.cycle.GetStep(press);
    step_press.bottom_boundary_condition_ = std::move(rep_stream);

    auto press_closed = std::make_unique<ClosedBoundaryCondition>();
    press_closed->location_ = BoundaryLocation::TOP;
    step_press.top_boundary_condition_ = std::move(press_closed);

    // CONFIGURE CYCLE
    // Run to cyclic steady state (paper's mass-balance criterion), minimum 100 cycles.
    test_fixture.ConfigureCycle(250); // hard cap on cycle count
    test_fixture.data.settings.css.method_ = CSSMethod::MASS_BALANCE_COMPONENT;
    test_fixture.data.settings.css.min_cycles_ = 100;

    // KPI: Feed(0)@Bottom, Product: Evacuation(2)@Bottom, Product component: CO2(0)
    test_fixture.ConfigureKPI(0, {0, 3}, {BoundaryLocation::BOTTOM, BoundaryLocation::BOTTOM}, {2}, {BoundaryLocation::BOTTOM});

    test_fixture.data.settings.max_step = 2.0;

    CoupledSolver solver{test_fixture.data};

    std::string error{};
    EXPECT_TRUE(test_fixture.data.PerformChecks(error)) << "System checks failed: " << error;

    bool verbose = false;
    int retval = solver.IntegrateCycle(test_fixture.data.cycle, verbose);
    EXPECT_GE(retval, 0) << "Solver did not complete";

    EXPECT_TRUE(test_fixture.IsMassConserved(solver.GetLog())) << "Mass is not conserved.";
    EXPECT_TRUE(test_fixture.IsMolesConserved(solver.GetLog())) << "Mass is not conserved.";
    EXPECT_TRUE(test_fixture.IsComponentMolesConserved(solver.GetLog())) << "Component moles not conserved.";

    // --- CO2 purity/recovery at cyclic steady state (vs Nguyen 2023 Table 2) ---
    // Basic four-step cycle: Feed(0), Blowdown(1), Evacuation(2), Pressurisation(3).
    // CO2 product leaves at the BOTTOM during Evacuation; CO2 is fed at the BOTTOM
    // during both Feed and Pressurisation. Evaluate at the last (CSS) cycle.
    {
      const EquationLog& eq_log = solver.GetLog();
      const BoundaryLog& b_log = solver.GetBoundaryLog();
      const auto& cyc = eq_log.sys_data.cycle_index_set_;
      int last_cycle = cyc.empty() ? 0 : *cyc.rbegin();
      KPIParams kpi(eq_log, b_log);
      kpi.product_comp_ = "CarbonDioxide";
      kpi.cycle_index_ = last_cycle;
      kpi.feed_endpoints = KpiEndpointsForFaces(test_fixture.data.cycle,
        { { 0, BoundaryLocation::BOTTOM }, { 3, BoundaryLocation::BOTTOM } });
      kpi.product_endpoints = KpiEndpointsForFaces(test_fixture.data.cycle,
        { { 2, BoundaryLocation::BOTTOM } });
      double purity = CalculatePurity(kpi) * 100.0;
      double recovery = CalculateRecovery(kpi) * 100.0;
      std::cout << "\n[KPI] " << title << " | cycles run = " << (last_cycle + 1)
                << " | CO2 purity = " << purity << " % | CO2 recovery = " << recovery << " %\n" << std::endl;
      ValidateNguyen2023(experiment_name, test_fixture.data, kpi, purity, recovery);
    }

    test_fixture.SaveTest(experiment_name, solver);
    test_fixture.SaveResultsAsExample(experiment_name, solver);
  }

  // Helper function to set up and run Nguyen 2023 4 step LPP
  void RunNguyen2023LPP(
    Nguyen2023Test& test_fixture,
    const std::string& experiment_name, 
    const std::string& title,
    const std::string& description,
    double feed_time,
    double interstitial_feed_vel,
    double bd_p,
    double evac_p
  ) {
    test_fixture.data.info.title = title;
    test_fixture.data.info.description = description;
    test_fixture.data.info.author = "Sam Affleck";
    
    const std::string CO2 = "CarbonDioxide";
    const std::string N2 = "Nitrogen";
    
    test_fixture.SetUpNguyen2023();

    const double press_time = 50; // [s]
    const double feed1_time = press_time; // [s]
    const double feed2_time = feed_time - feed1_time; // [s]
    const double blowdown_time = 60; // [s]
    const double evac_time = 40; // [s]

    const std::string feed1 = "Feed 1";
    const std::string feed2 = "Feed 2";
    const std::string blowdown = "Blowdown";
    const std::string evacuation = "Evacuation";
    const std::string press = "Pressurisation";

    test_fixture.AddStep(feed1, feed1_time);
    test_fixture.AddStep(feed2, feed2_time);
    test_fixture.AddStep(blowdown, blowdown_time);
    test_fixture.AddStep(evacuation, evac_time);
    test_fixture.AddStep(press, press_time);

    // FEED1 STEP
    const double feed_vel = interstitial_feed_vel * 0.40;
    const double feed_temp = 293;
    const double ambient_pressure = 1.0e5;
    const double ambient_valve_cv = 10.0;
    const std::unordered_map<std::string, double> inlet_mf = {
      {CO2, 0.15 },
      {N2, 0.85 }
    };

    // Step 1: Feed 1
    Step& step_feed1 = test_fixture.data.cycle.GetStep(feed1);
    auto feed1_inlet_stream = test_fixture.CreateNoValveInletStream(
      test_fixture.CreateVelocityBC(feed_vel),
      inlet_mf,
      feed_temp
    );
    feed1_inlet_stream->name_ = "Feed 1 inlet";
    
    auto feed1_flow_bc = std::make_unique<FlowBoundaryCondition>();
    feed1_flow_bc->direction_ = BoundaryDirection::INLET;
    feed1_flow_bc->location_ = BoundaryLocation::BOTTOM;
    feed1_flow_bc->streams_.push_back(std::move(feed1_inlet_stream));
    step_feed1.bottom_boundary_condition_ = std::move(feed1_flow_bc);
    
    auto feed1_ambient_stream = test_fixture.CreateValveOutletStream(
      test_fixture.CreateConstantBoundaryValue(ambient_pressure),
      ambient_valve_cv
    );
    feed1_ambient_stream->name_ = "Feed 1 ambient Stream";
    
    auto feed1_LPP_stream = test_fixture.CreateNoValveOutletStream(
      test_fixture.CreateVelocityBC(0.0)
    );
    feed1_LPP_stream->name_ = "LPP Stream";

    // FEED 2 STEP
    Step& step_feed2 = test_fixture.data.cycle.GetStep(feed2);
    auto feed2_inlet_stream = test_fixture.CreateNoValveInletStream(
      test_fixture.CreateVelocityBC(feed_vel),
      inlet_mf,
      feed_temp
    );
    feed2_inlet_stream->name_ = "Feed 2 inlet";
    
    auto feed2_flow_bc_in = std::make_unique<FlowBoundaryCondition>();
    feed2_flow_bc_in->direction_ = BoundaryDirection::INLET;
    feed2_flow_bc_in->location_ = BoundaryLocation::BOTTOM;
    feed2_flow_bc_in->streams_.push_back(std::move(feed2_inlet_stream));
    step_feed2.bottom_boundary_condition_ = std::move(feed2_flow_bc_in);

    auto feed2_ambient_stream = test_fixture.CreateValveOutletStream(
      test_fixture.CreateConstantBoundaryValue(ambient_pressure),
      ambient_valve_cv
    );
    feed2_ambient_stream->name_ = "Feed 2 ambient Stream";
    
    auto feed2_flow_bc_out = std::make_unique<FlowBoundaryCondition>();
    feed2_flow_bc_out->direction_ = BoundaryDirection::OUTLET;
    feed2_flow_bc_out->location_ = BoundaryLocation::TOP;
    feed2_flow_bc_out->streams_.push_back(std::move(feed2_ambient_stream));
    step_feed2.top_boundary_condition_ = std::move(feed2_flow_bc_out);

    // ==== PRESSURISATION ====
    auto rep_stream = test_fixture.CreatePressureInletBoundary(
      test_fixture.CreateLinearBoundaryValue(evac_p, ambient_pressure, press_time),
      { 
        { CO2, 0.0 },
        { N2, 1.0 } 
      },
      293.0
    );
    rep_stream->location_ = BoundaryLocation::TOP;
    rep_stream->direction_ = BoundaryDirection::INLET;

    // Make connecting streams interact
    // TODO: Need to check how connections work in new architecture
    auto endA = MakeEndpoint(test_fixture.data.cycle.GetStep(press), *rep_stream);
    auto endB = MakeEndpoint(test_fixture.data.cycle.GetStep(feed1), *feed1_LPP_stream);
    MakeMutualConnection(endA, endB);

    auto feed1_flow_bc_out = std::make_unique<FlowBoundaryCondition>();
    feed1_flow_bc_out->direction_ = BoundaryDirection::OUTLET;
    feed1_flow_bc_out->location_ = BoundaryLocation::TOP;
    feed1_flow_bc_out->streams_.push_back(std::move(feed1_ambient_stream));
    feed1_flow_bc_out->streams_.push_back(std::move(feed1_LPP_stream));
    step_feed1.top_boundary_condition_ = std::move(feed1_flow_bc_out);

    // BLOWDOWN
    Step& step_blowdown = test_fixture.data.cycle.GetStep(blowdown);
    auto blowdown_pressure = std::make_unique<PressureBoundaryCondition>();
    blowdown_pressure->direction_ = BoundaryDirection::OUTLET;
    blowdown_pressure->location_ = BoundaryLocation::TOP;
    blowdown_pressure->pressure_ = test_fixture.CreateExponentialBoundaryValue(bd_p, blowdown_time);
    blowdown_pressure->mole_fraction_[CO2] = 0.0;
    blowdown_pressure->mole_fraction_[N2] = 1.0;
    step_blowdown.top_boundary_condition_ = std::move(blowdown_pressure);
    
    auto blowdown_closed = std::make_unique<ClosedBoundaryCondition>();
    blowdown_closed->location_ = BoundaryLocation::BOTTOM;
    step_blowdown.bottom_boundary_condition_ = std::move(blowdown_closed);

    // EVACUATION
    Step& step_evac = test_fixture.data.cycle.GetStep(evacuation);
    auto evac_closed_top = std::make_unique<ClosedBoundaryCondition>();
    evac_closed_top->location_ = BoundaryLocation::TOP;
    step_evac.top_boundary_condition_ = std::move(evac_closed_top);
    
    auto evac_pressure = std::make_unique<PressureBoundaryCondition>();
    evac_pressure->direction_ = BoundaryDirection::OUTLET;
    evac_pressure->location_ = BoundaryLocation::BOTTOM;
    evac_pressure->pressure_ = test_fixture.CreateExponentialBoundaryValue(evac_p, evac_time);
    evac_pressure->mole_fraction_[CO2] = 0.0;
    evac_pressure->mole_fraction_[N2] = 1.0;
    step_evac.bottom_boundary_condition_ = std::move(evac_pressure);
    
    // PRESSURISATION
    Step& step_press = test_fixture.data.cycle.GetStep(press);
    step_press.top_boundary_condition_ = std::move(rep_stream);

    auto press_closed = std::make_unique<ClosedBoundaryCondition>();
    press_closed->location_ = BoundaryLocation::BOTTOM;
    step_press.bottom_boundary_condition_ = std::move(press_closed);

    // CONFIGURE CYCLE
    // Run to cyclic steady state (paper's mass-balance criterion), minimum 100 cycles.
    test_fixture.ConfigureCycle(250); // hard cap on cycle count
    test_fixture.data.settings.css.method_ = CSSMethod::MASS_BALANCE_COMPONENT;
    test_fixture.data.settings.css.min_cycles_ = 100;

    // KPI: Feed1(0),Feed2(1)@Bottom; Product: Evacuation(3)@Bottom; Recycled product:
    // light-product pressurisation (4)@Top (CO2 re-entering the bed); Product component: CO2(0)
    test_fixture.ConfigureKPI(0, {0, 1}, {BoundaryLocation::BOTTOM, BoundaryLocation::BOTTOM}, {3}, {BoundaryLocation::BOTTOM},
                              {4}, {BoundaryLocation::TOP});

    test_fixture.data.settings.max_step = 2.0;

    CoupledSolver solver{test_fixture.data};

    std::string error{};
    EXPECT_TRUE(test_fixture.data.PerformChecks(error)) << "System checks failed: " << error;

    bool verbose = false;
    int retval = solver.IntegrateCycle(test_fixture.data.cycle, verbose);
    EXPECT_GE(retval, 0) << "Solver did not complete";

    EXPECT_TRUE(test_fixture.IsMassConserved(solver.GetLog())) << "Mass is not conserved.";
    EXPECT_TRUE(test_fixture.IsMolesConserved(solver.GetLog())) << "Mass is not conserved.";
    EXPECT_TRUE(test_fixture.IsComponentMolesConserved(solver.GetLog())) << "Component moles not conserved.";

    // --- CO2 purity/recovery at cyclic steady state (vs Nguyen 2023 Table 2) ---
    // LPP cycle: Feed1(0), Feed2(1), Blowdown(2), Evacuation(3), Pressurisation(4).
    // CO2 product leaves at the BOTTOM during Evacuation; CO2 is fed at the BOTTOM
    // during both feed steps. LPP repressurisation uses light product, not feed.
    {
      const EquationLog& eq_log = solver.GetLog();
      const BoundaryLog& b_log = solver.GetBoundaryLog();
      const auto& cyc = eq_log.sys_data.cycle_index_set_;
      int last_cycle = cyc.empty() ? 0 : *cyc.rbegin();
      KPIParams kpi(eq_log, b_log);
      kpi.product_comp_ = "CarbonDioxide";
      kpi.cycle_index_ = last_cycle;
      kpi.feed_endpoints = KpiEndpointsForFaces(test_fixture.data.cycle,
        { { 0, BoundaryLocation::BOTTOM }, { 1, BoundaryLocation::BOTTOM } });
      kpi.product_endpoints = KpiEndpointsForFaces(test_fixture.data.cycle,
        { { 3, BoundaryLocation::BOTTOM } });
      // Validation KPIs follow the papers: purity of the evacuation product, recovery
      // = CO2 in that product / CO2 fed in the feed steps. (Folding the light-product
      // pressurisation stream into the product, as the deck's KPI config does,
      // makes the product-total moles nearly cancel and the purity meaningless.)
      double purity = CalculatePurity(kpi) * 100.0;
      double recovery = CalculateRecovery(kpi) * 100.0;
      std::cout << "\n[KPI] " << title << " | cycles run = " << (last_cycle + 1)
                << " | CO2 purity = " << purity << " % | CO2 recovery = " << recovery << " %\n" << std::endl;
      ValidateNguyen2023(experiment_name, test_fixture.data, kpi, purity, recovery);
    }

    test_fixture.SaveTest(experiment_name, solver);
    test_fixture.SaveResultsAsExample(experiment_name, solver);
  }
  TEST_F(Nguyen2023Test, ExperimentI) {
    RunNguyen2023Experiment(*this,
      "Nguyen2023_Exp_I",
      "Nguyen 2023 - Experiment I.",
      "Basic four-step cycle for a binary CO2/N2 VPSA cycle using CALF-20 taken from https://doi.org/10.1016/j.cej.2022.139550. This configuration corresponds to experiment I in the paper.",
      176,    // feed_time [s]
      0.093,  // interstitial_feed_vel
      0.119e5 // bd_pressure [Pa]
    );
  }

  TEST_F(Nguyen2023Test, ExperimentII) {
    RunNguyen2023Experiment(*this,
      "Nguyen2023_Exp_II",
      "Nguyen 2023 - Experiment II.",
      "Basic four-step cycle for a binary CO2/N2 VPSA cycle using CALF-20 taken from https://doi.org/10.1016/j.cej.2022.139550. This configuration corresponds to experiment II in the paper.",
      174,    // feed_time [s]
      0.065,  // interstitial_feed_vel
      0.148e5 // bd_pressure [Pa]
    );
  }

  TEST_F(Nguyen2023Test, ExperimentIII) {
    RunNguyen2023Experiment(*this,
      "Nguyen2023_Exp_III",
      "Nguyen 2023 - Experiment III.",
      "Basic four-step cycle for a binary CO2/N2 VPSA cycle using CALF-20 taken from https://doi.org/10.1016/j.cej.2022.139550. This configuration corresponds to experiment III in the paper.",
      174,    // feed_time [s]
      0.058,  // interstitial_feed_vel
      0.189e5 // bd_pressure [Pa]
    );
  }
  
  TEST_F(Nguyen2023Test, ExperimentIV) {
    RunNguyen2023Experiment(*this,
      "Nguyen2023_Exp_IV",
      "Nguyen 2023 - Experiment IV.",
      "Basic four-step cycle for a binary CO2/N2 VPSA cycle using CALF-20 taken from https://doi.org/10.1016/j.cej.2022.139550. This configuration corresponds to experiment IV in the paper.",
      174,    // feed_time [s]
      0.055,  // interstitial_feed_vel
      0.233e5 // bd_pressure [Pa]
    );
  }

  TEST_F(Nguyen2023Test, ExperimentV) {
    RunNguyen2023Experiment(*this,
      "Nguyen2023_Exp_V",
      "Nguyen 2023 - Experiment V.",
      "Basic four-step cycle for a binary CO2/N2 VPSA cycle using CALF-20 taken from https://doi.org/10.1016/j.cej.2022.139550. This configuration corresponds to experiment V in the paper.",
      174,    // feed_time [s]
      0.046,  // interstitial_feed_vel
      0.263e5 // bd_pressure [Pa]
    );
  }

  TEST_F(Nguyen2023Test, ExperimentVI) {
    RunNguyen2023LPP(*this,
      "Nguyen2023_Exp_VI",
      "Nguyen 2023 - Experiment VI.",
      "Four-step LPP cycle for a binary CO2/N2 VPSA cycle using CALF-20 taken from https://doi.org/10.1016/j.cej.2022.139550. This configuration corresponds to experiment VI in the paper.",
      193,    // feed_time [s]
      0.041,  // interstitial_feed_vel
      0.063e5, // bd_pressure [Pa]
      0.031e5 // evac_pressure [Pa]
    );
  }

  TEST_F(Nguyen2023Test, ExperimentVII) {
    RunNguyen2023LPP(*this,
      "Nguyen2023_Exp_VII",
      "Nguyen 2023 - Experiment VII.",
      "Four-step LPP cycle for a binary CO2/N2 VPSA cycle using CALF-20 taken from https://doi.org/10.1016/j.cej.2022.139550. This configuration corresponds to experiment VII in the paper.",
      157,    // feed_time [s]
      0.024,  // interstitial_feed_vel
      0.046e5, // bd_pressure [Pa]
      0.031e5 // evac_pressure [Pa]
    );
  }

  TEST_F(Nguyen2023Test, ExperimentVIII) {
    RunNguyen2023LPP(*this,
      "Nguyen2023_Exp_VIII",
      "Nguyen 2023 - Experiment VIII.",
      "Four-step LPP cycle for a binary CO2/N2 VPSA cycle using CALF-20 taken from https://doi.org/10.1016/j.cej.2022.139550. This configuration corresponds to experiment VIII in the paper.",
      162,    // feed_time [s]
      0.023,  // interstitial_feed_vel
      0.060e5, // bd_pressure [Pa]
      0.034e5 // evac_pressure [Pa]
    );
  }

  TEST_F(Nguyen2023Test, ExperimentIX) {
    RunNguyen2023LPP(*this,
      "Nguyen2023_Exp_IX",
      "Nguyen 2023 - Experiment IX.",
      "Four-step LPP cycle for a binary CO2/N2 VPSA cycle using CALF-20 taken from https://doi.org/10.1016/j.cej.2022.139550. This configuration corresponds to experiment IX in the paper.",
      156,    // feed_time [s]
      0.021,  // interstitial_feed_vel
      0.066e5, // bd_pressure [Pa]
      0.031e5 // evac_pressure [Pa]
    );
  }

  TEST_F(Nguyen2023Test, ExperimentX) {
    RunNguyen2023LPP(*this,
      "Nguyen2023_Exp_X",
      "Nguyen 2023 - Experiment X.",
      "Four-step LPP cycle for a binary CO2/N2 VPSA cycle using CALF-20 taken from https://doi.org/10.1016/j.cej.2022.139550. This configuration corresponds to experiment X in the paper.",
      156,    // feed_time [s]
      0.022,  // interstitial_feed_vel
      0.080e5, // bd_pressure [Pa]
      0.031e5 // evac_pressure [Pa]
    );
  }

} // namespace adsim
