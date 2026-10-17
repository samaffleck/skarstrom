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

// Paper-specific fixture + setup (Perez et al. 2019)
#include "perez_2019.h"
#include "../../fixtures/validation_metrics.h"
#include <limits>
#include <map>

namespace adsim {

  // ==========================================================================
  //  Quantitative comparison with Table 3 of Estupinan Perez, Sarkar &
  //  Rajendran, Sep. Purif. Technol. 224 (2019) 553-563,
  //  doi:10.1016/j.seppur.2019.05.039: CO2 purity / recovery at CSS against
  //  the experiment (validation) and the paper's own model (code-to-code).
  //  Absolute error in percentage points; tolerances from the measured error,
  //  rounded up (kPerezTol).
  // ==========================================================================
  namespace {
    struct PerezTol { double pu_exp, re_exp, pu_model, re_model; };
    const std::map<std::string, PerezTol> kPerezTol = {
      { "ConfigA",  { 1.1, 12, 0.67, 9.6 } },  // measured 0.760915, 9.62585, 0.360915, 8.32585 pp
      { "ConfigB",  { 0.41, 9.9, 1.1, 7.2 } },  // measured 0.103064, 8.52423, 0.796936, 6.22423 pp
      { "ConfigC1", { 1.6, 9.5, 4.1, 5.3 } },  // measured 1.21589, 8.18045, 3.48411, 4.58045 pp
      { "ConfigD1", { 3.9, 4.8, 2.8, 4 } },  // measured 3.37943, 4.12065, 2.37943, 3.42065 pp
      { "ConfigE",  { 5.8, 5, 8, 5.7 } },  // measured 4.99608, 4.276, 6.89608, 4.876 pp
      { "ConfigF",  { 11, 3.9, 13, 5 } },  // measured 8.85068, 3.31355, 11.0507, 4.31355 pp
    };

    void ValidatePerez(const std::string& id, double purity, double recovery) {
      using namespace validation;
      const auto& tol = kPerezTol.at(id);
      CaseReport report("perez_2019", id);
      report.CheckAbs("CO2 purity vs experiment", purity, Reference("perez_2019", id, "CO2_purity_exp"), tol.pu_exp);
      report.CheckAbs("CO2 recovery vs experiment", recovery, Reference("perez_2019", id, "CO2_recovery_exp"), tol.re_exp);
      report.CheckAbs("CO2 purity vs paper model", purity, Reference("perez_2019", id, "CO2_purity_model"), tol.pu_model);
      report.CheckAbs("CO2 recovery vs paper model", recovery, Reference("perez_2019", id, "CO2_recovery_model"), tol.re_model);
      if (id == "ConfigC1" || id == "ConfigD1") {
        const std::string rep = id == "ConfigC1" ? "C2" : "D2";
        report.Info("CO2 purity vs repeat experiment " + rep, purity,
                    Reference("perez_2019", id, "CO2_purity_exp_" + rep), "abs_error", "repeat run");
        report.Info("CO2 recovery vs repeat experiment " + rep, recovery,
                    Reference("perez_2019", id, "CO2_recovery_exp_" + rep), "abs_error", "repeat run");
      }
    }
  } // namespace

  // ---- Basic 4-step VSA cycle (Perez 2019 Table 3, experiments A/B/C) ----
  void Perez2019Test::RunPerez4Step(const std::string& save_name, const std::string& title,
      double feed_time, double interstitial_feed_vel, double bd_pressure, double evac_pressure,
      double blowdown_time, double evac_time, double press_time,
      double target_purity, double target_recovery) {
    data.info.title = title;
    data.info.description = title + " - Perez et al. (2019) basic four-step VSA cycle, CO2/N2 on 13X.";
    data.info.author = "Sam Affleck";

    const std::string CO2 = "CarbonDioxide";
    const std::string N2 = "Nitrogen";
    SetUpPerez2019();

    const std::string feed = "Feed";
    const std::string blowdown = "Blowdown";
    const std::string evacuation = "Evacuation";
    const std::string press = "Pressurisation";
    AddStep(feed, feed_time);
    AddStep(blowdown, blowdown_time);
    AddStep(evacuation, evac_time);
    AddStep(press, press_time);

    const double feed_vel = interstitial_feed_vel * 0.41;
    const double feed_temp = 298;
    const double ambient_pressure = 1.0e5;
    const std::unordered_map<std::string, double> inlet_mf = { {CO2, 0.15}, {N2, 0.85} };

    // FEED
    Step& step_feed = data.cycle.GetStep(feed);
    auto feed_inlet_stream = CreateNoValveInletStream(CreateVelocityBC(feed_vel), inlet_mf, feed_temp);
    feed_inlet_stream->name_ = "Feed inlet";
    auto feed_flow_bc = std::make_unique<FlowBoundaryCondition>();
    feed_flow_bc->direction_ = BoundaryDirection::INLET;
    feed_flow_bc->location_ = BoundaryLocation::BOTTOM;
    feed_flow_bc->streams_.push_back(std::move(feed_inlet_stream));
    step_feed.bottom_boundary_condition_ = std::move(feed_flow_bc);
    auto feed_pressure_out = std::make_unique<PressureBoundaryCondition>();
    feed_pressure_out->direction_ = BoundaryDirection::OUTLET;
    feed_pressure_out->location_ = BoundaryLocation::TOP;
    feed_pressure_out->pressure_ = CreateConstantBoundaryValue(ambient_pressure);
    feed_pressure_out->mole_fraction_[CO2] = 0.0;
    feed_pressure_out->mole_fraction_[N2] = 1.0;
    step_feed.top_boundary_condition_ = std::move(feed_pressure_out);

    // PRESSURISATION (feed gas, bottom)
    auto rep_stream = CreatePressureInletBoundary(
      CreateLinearBoundaryValue(evac_pressure, ambient_pressure, press_time), inlet_mf, feed_temp);
    rep_stream->location_ = BoundaryLocation::BOTTOM;
    rep_stream->direction_ = BoundaryDirection::INLET;

    // BLOWDOWN (co-current, top)
    Step& step_blowdown = data.cycle.GetStep(blowdown);
    auto blowdown_pressure = std::make_unique<PressureBoundaryCondition>();
    blowdown_pressure->direction_ = BoundaryDirection::OUTLET;
    blowdown_pressure->location_ = BoundaryLocation::TOP;
    blowdown_pressure->pressure_ = CreateExponentialBoundaryValue(bd_pressure, blowdown_time);
    blowdown_pressure->mole_fraction_[CO2] = 0.0;
    blowdown_pressure->mole_fraction_[N2] = 1.0;
    step_blowdown.top_boundary_condition_ = std::move(blowdown_pressure);
    auto blowdown_closed = std::make_unique<ClosedBoundaryCondition>();
    blowdown_closed->location_ = BoundaryLocation::BOTTOM;
    step_blowdown.bottom_boundary_condition_ = std::move(blowdown_closed);

    // EVACUATION (counter-current, bottom = CO2 product)
    Step& step_evac = data.cycle.GetStep(evacuation);
    auto evac_closed_top = std::make_unique<ClosedBoundaryCondition>();
    evac_closed_top->location_ = BoundaryLocation::TOP;
    step_evac.top_boundary_condition_ = std::move(evac_closed_top);
    auto evac_pressure_bc = std::make_unique<PressureBoundaryCondition>();
    evac_pressure_bc->direction_ = BoundaryDirection::OUTLET;
    evac_pressure_bc->location_ = BoundaryLocation::BOTTOM;
    evac_pressure_bc->pressure_ = CreateExponentialBoundaryValue(evac_pressure, evac_time);
    evac_pressure_bc->mole_fraction_[CO2] = 0.0;
    evac_pressure_bc->mole_fraction_[N2] = 1.0;
    step_evac.bottom_boundary_condition_ = std::move(evac_pressure_bc);

    // PRESSURISATION assignment
    Step& step_press = data.cycle.GetStep(press);
    step_press.bottom_boundary_condition_ = std::move(rep_stream);
    auto press_closed = std::make_unique<ClosedBoundaryCondition>();
    press_closed->location_ = BoundaryLocation::TOP;
    step_press.top_boundary_condition_ = std::move(press_closed);

    // Run to CSS, min 100 cycles
    ConfigureCycle(250);
    data.settings.css.method_ = CSSMethod::MASS_BALANCE_COMPONENT;
    data.settings.css.min_cycles_ = 100;
    ConfigureKPI(0, {0, 3}, {BoundaryLocation::BOTTOM, BoundaryLocation::BOTTOM}, {2}, {BoundaryLocation::BOTTOM});
    data.settings.max_step = 2.0;

    CoupledSolver solver{data};
    std::string error{};
    EXPECT_TRUE(data.PerformChecks(error)) << "System checks failed: " << error;
    EXPECT_GE(solver.IntegrateCycle(data.cycle, false), 0) << "Solver did not complete";
    EXPECT_TRUE(IsMassConserved(solver.GetLog())) << "Mass is not conserved.";
    EXPECT_TRUE(IsMolesConserved(solver.GetLog())) << "Moles not conserved.";
    EXPECT_TRUE(IsComponentMolesConserved(solver.GetLog())) << "Component moles not conserved.";

    const EquationLog& eq_log = solver.GetLog();
    const BoundaryLog& b_log = solver.GetBoundaryLog();
    const auto& cyc = eq_log.sys_data.cycle_index_set_;
    int last_cycle = cyc.empty() ? 0 : *cyc.rbegin();
    KPIParams kpi(eq_log, b_log);
    kpi.product_comp_ = "CarbonDioxide";
    kpi.cycle_index_ = last_cycle;
    kpi.feed_endpoints = KpiEndpointsForFaces(data.cycle,
      { { 0, BoundaryLocation::BOTTOM }, { 3, BoundaryLocation::BOTTOM } });
    kpi.product_endpoints = KpiEndpointsForFaces(data.cycle,
      { { 2, BoundaryLocation::BOTTOM } });
    double purity = CalculatePurity(kpi) * 100.0;
    double recovery = CalculateRecovery(kpi) * 100.0;
    std::cout << "\n[KPI] Perez 2019 " << title << " | cycles run = " << (last_cycle + 1)
              << " | CO2 purity = " << purity << " % (paper " << target_purity << ") | CO2 recovery = "
              << recovery << " % (paper " << target_recovery << ")\n" << std::endl;
    ValidatePerez("Config" + save_name.substr(save_name.rfind('_') + 1), purity, recovery);

    SaveTest(save_name, solver);
    SaveResultsAsExample(save_name, solver);
  }

  // ---- 4-step LPP cycle (Perez 2019 Table 3, experiments D/E/F) ----
  void Perez2019Test::RunPerezLPP(const std::string& save_name, const std::string& title,
      double feed_time, double interstitial_feed_vel, double bd_pressure, double evac_pressure,
      double blowdown_time, double evac_time, double press_time,
      double target_purity, double target_recovery) {
    data.info.title = title;
    data.info.description = title + " - Perez et al. (2019) four-step VSA cycle with light-product pressurisation, CO2/N2 on 13X.";
    data.info.author = "Sam Affleck";

    const std::string CO2 = "CarbonDioxide";
    const std::string N2 = "Nitrogen";
    SetUpPerez2019();

    const double feed1_time = press_time;                 // ADS during the other column's LPP
    const double feed2_time = feed_time - press_time;     // remainder of the adsorption step
    const std::string feed1 = "Feed 1";
    const std::string feed2 = "Feed 2";
    const std::string blowdown = "Blowdown";
    const std::string evacuation = "Evacuation";
    const std::string press = "Pressurisation";
    AddStep(feed1, feed1_time);
    AddStep(feed2, feed2_time);
    AddStep(blowdown, blowdown_time);
    AddStep(evacuation, evac_time);
    AddStep(press, press_time);

    const double feed_vel = interstitial_feed_vel * 0.41;
    const double feed_temp = 298;
    const double ambient_pressure = 1.0e5;
    const double ambient_valve_cv = 10.0;
    const std::unordered_map<std::string, double> inlet_mf = { {CO2, 0.15}, {N2, 0.85} };

    // FEED 1 (dual outlet: ambient vent + recycled LPP stream)
    Step& step_feed1 = data.cycle.GetStep(feed1);
    auto feed1_inlet_stream = CreateNoValveInletStream(CreateVelocityBC(feed_vel), inlet_mf, feed_temp);
    feed1_inlet_stream->name_ = "Feed 1 inlet";
    auto feed1_flow_bc = std::make_unique<FlowBoundaryCondition>();
    feed1_flow_bc->direction_ = BoundaryDirection::INLET;
    feed1_flow_bc->location_ = BoundaryLocation::BOTTOM;
    feed1_flow_bc->streams_.push_back(std::move(feed1_inlet_stream));
    step_feed1.bottom_boundary_condition_ = std::move(feed1_flow_bc);
    auto feed1_ambient_stream = CreateValveOutletStream(CreateConstantBoundaryValue(ambient_pressure), ambient_valve_cv);
    feed1_ambient_stream->name_ = "Feed 1 ambient Stream";
    auto feed1_LPP_stream = CreateNoValveOutletStream(CreateVelocityBC(0.0));
    feed1_LPP_stream->name_ = "LPP Stream";

    // FEED 2 (single ambient outlet)
    Step& step_feed2 = data.cycle.GetStep(feed2);
    auto feed2_inlet_stream = CreateNoValveInletStream(CreateVelocityBC(feed_vel), inlet_mf, feed_temp);
    feed2_inlet_stream->name_ = "Feed 2 inlet";
    auto feed2_flow_bc_in = std::make_unique<FlowBoundaryCondition>();
    feed2_flow_bc_in->direction_ = BoundaryDirection::INLET;
    feed2_flow_bc_in->location_ = BoundaryLocation::BOTTOM;
    feed2_flow_bc_in->streams_.push_back(std::move(feed2_inlet_stream));
    step_feed2.bottom_boundary_condition_ = std::move(feed2_flow_bc_in);
    auto feed2_ambient_stream = CreateValveOutletStream(CreateConstantBoundaryValue(ambient_pressure), ambient_valve_cv);
    feed2_ambient_stream->name_ = "Feed 2 ambient Stream";
    auto feed2_flow_bc_out = std::make_unique<FlowBoundaryCondition>();
    feed2_flow_bc_out->direction_ = BoundaryDirection::OUTLET;
    feed2_flow_bc_out->location_ = BoundaryLocation::TOP;
    feed2_flow_bc_out->streams_.push_back(std::move(feed2_ambient_stream));
    step_feed2.top_boundary_condition_ = std::move(feed2_flow_bc_out);

    // PRESSURISATION (light product, top), connected to the Feed1 LPP stream
    auto rep_stream = CreatePressureInletBoundary(
      CreateLinearBoundaryValue(evac_pressure, ambient_pressure, press_time), { {CO2, 0.0}, {N2, 1.0} }, feed_temp);
    rep_stream->location_ = BoundaryLocation::TOP;
    rep_stream->direction_ = BoundaryDirection::INLET;
    auto endA = MakeEndpoint(data.cycle.GetStep(press), *rep_stream);
    auto endB = MakeEndpoint(data.cycle.GetStep(feed1), *feed1_LPP_stream);
    MakeMutualConnection(endA, endB);
    auto feed1_flow_bc_out = std::make_unique<FlowBoundaryCondition>();
    feed1_flow_bc_out->direction_ = BoundaryDirection::OUTLET;
    feed1_flow_bc_out->location_ = BoundaryLocation::TOP;
    feed1_flow_bc_out->streams_.push_back(std::move(feed1_ambient_stream));
    feed1_flow_bc_out->streams_.push_back(std::move(feed1_LPP_stream));
    step_feed1.top_boundary_condition_ = std::move(feed1_flow_bc_out);

    // BLOWDOWN (co-current, top)
    Step& step_blowdown = data.cycle.GetStep(blowdown);
    auto blowdown_pressure = std::make_unique<PressureBoundaryCondition>();
    blowdown_pressure->direction_ = BoundaryDirection::OUTLET;
    blowdown_pressure->location_ = BoundaryLocation::TOP;
    blowdown_pressure->pressure_ = CreateExponentialBoundaryValue(bd_pressure, blowdown_time);
    blowdown_pressure->mole_fraction_[CO2] = 0.0;
    blowdown_pressure->mole_fraction_[N2] = 1.0;
    step_blowdown.top_boundary_condition_ = std::move(blowdown_pressure);
    auto blowdown_closed = std::make_unique<ClosedBoundaryCondition>();
    blowdown_closed->location_ = BoundaryLocation::BOTTOM;
    step_blowdown.bottom_boundary_condition_ = std::move(blowdown_closed);

    // EVACUATION (counter-current, bottom = CO2 product)
    Step& step_evac = data.cycle.GetStep(evacuation);
    auto evac_closed_top = std::make_unique<ClosedBoundaryCondition>();
    evac_closed_top->location_ = BoundaryLocation::TOP;
    step_evac.top_boundary_condition_ = std::move(evac_closed_top);
    auto evac_pressure_bc = std::make_unique<PressureBoundaryCondition>();
    evac_pressure_bc->direction_ = BoundaryDirection::OUTLET;
    evac_pressure_bc->location_ = BoundaryLocation::BOTTOM;
    evac_pressure_bc->pressure_ = CreateExponentialBoundaryValue(evac_pressure, evac_time);
    evac_pressure_bc->mole_fraction_[CO2] = 0.0;
    evac_pressure_bc->mole_fraction_[N2] = 1.0;
    step_evac.bottom_boundary_condition_ = std::move(evac_pressure_bc);

    // PRESSURISATION assignment
    Step& step_press = data.cycle.GetStep(press);
    step_press.top_boundary_condition_ = std::move(rep_stream);
    auto press_closed = std::make_unique<ClosedBoundaryCondition>();
    press_closed->location_ = BoundaryLocation::BOTTOM;
    step_press.bottom_boundary_condition_ = std::move(press_closed);

    // Run to CSS, min 100 cycles
    ConfigureCycle(250);
    data.settings.css.method_ = CSSMethod::MASS_BALANCE_COMPONENT;
    data.settings.css.min_cycles_ = 100;
    ConfigureKPI(0, {0, 1}, {BoundaryLocation::BOTTOM, BoundaryLocation::BOTTOM}, {3}, {BoundaryLocation::BOTTOM}, {4}, {BoundaryLocation::TOP});
    data.settings.max_step = 2.0;

    CoupledSolver solver{data};
    std::string error{};
    EXPECT_TRUE(data.PerformChecks(error)) << "System checks failed: " << error;
    EXPECT_GE(solver.IntegrateCycle(data.cycle, false), 0) << "Solver did not complete";
    EXPECT_TRUE(IsMassConserved(solver.GetLog())) << "Mass is not conserved.";
    EXPECT_TRUE(IsMolesConserved(solver.GetLog())) << "Moles not conserved.";
    EXPECT_TRUE(IsComponentMolesConserved(solver.GetLog())) << "Component moles not conserved.";

    const EquationLog& eq_log = solver.GetLog();
    const BoundaryLog& b_log = solver.GetBoundaryLog();
    const auto& cyc = eq_log.sys_data.cycle_index_set_;
    int last_cycle = cyc.empty() ? 0 : *cyc.rbegin();
    KPIParams kpi(eq_log, b_log);
    kpi.product_comp_ = "CarbonDioxide";
    kpi.cycle_index_ = last_cycle;
    kpi.feed_endpoints = KpiEndpointsForFaces(data.cycle,
      { { 0, BoundaryLocation::BOTTOM }, { 1, BoundaryLocation::BOTTOM } });
    kpi.product_endpoints = KpiEndpointsForFaces(data.cycle,
      { { 3, BoundaryLocation::BOTTOM } });
    // Validation KPIs follow the papers: purity of the evacuation product, recovery
    // = CO2 in that product / CO2 fed in the feed steps. (Folding the light-product
    // pressurisation stream into the product, as the deck's KPI config does,
    // makes the product-total moles nearly cancel and the purity meaningless.)
    double purity = CalculatePurity(kpi) * 100.0;
    double recovery = CalculateRecovery(kpi) * 100.0;
    std::cout << "\n[KPI] Perez 2019 " << title << " | cycles run = " << (last_cycle + 1)
              << " | CO2 purity = " << purity << " % (paper " << target_purity << ") | CO2 recovery = "
              << recovery << " % (paper " << target_recovery << ")\n" << std::endl;
    ValidatePerez("Config" + save_name.substr(save_name.rfind('_') + 1), purity, recovery);

    SaveTest(save_name, solver);
    SaveResultsAsExample(save_name, solver);
  }

  TEST_F(Perez2019Test, ConfigD1_ConnectedBC) {
    const std::string CO2 = "CarbonDioxide";
    const std::string N2  = "Nitrogen";

    SetUpPerez2019();

    const double press_time    = 150.0;
    const double feed1_time    = press_time;
    const double feed2_time    = 117.0;
    const double blowdown_time = 21.0;
    const double evac_time     = 167.0;

    const std::string feed1      = "Feed 1";
    const std::string feed2      = "Feed 2";
    const std::string blowdown   = "Blowdown";
    const std::string evacuation = "Evacuation";
    const std::string press      = "Pressurisation";

    AddStep(feed1, feed1_time);
    AddStep(feed2, feed2_time);
    AddStep(blowdown, blowdown_time);
    AddStep(evacuation, evac_time);
    AddStep(press, press_time);

    const double feed_vel        = 0.2;
    const double feed_temp       = 298.0;
    const double ambient_pressure = 1.0e5;
    const double ambient_valve_cv = 10.0;
    const std::unordered_map<std::string, double> inlet_mf = { {CO2, 0.15}, {N2, 0.85} };

    // Feed 1 – bottom inlet
    Step& step_feed1 = data.cycle.GetStep(feed1);
    auto feed1_inlet = CreateNoValveInletStream(CreateVelocityBC(feed_vel), inlet_mf, feed_temp);
    feed1_inlet->name_ = "Feed 1 inlet";
    auto feed1_in_bc = std::make_unique<FlowBoundaryCondition>();
    feed1_in_bc->direction_ = BoundaryDirection::INLET;
    feed1_in_bc->location_  = BoundaryLocation::BOTTOM;
    feed1_in_bc->streams_.push_back(std::move(feed1_inlet));
    step_feed1.bottom_boundary_condition_ = std::move(feed1_in_bc);

    // Feed 1 – top outlets (ambient + LPP connected stream)
    auto feed1_ambient = CreateValveOutletStream(CreateConstantBoundaryValue(ambient_pressure), ambient_valve_cv);
    feed1_ambient->name_ = "Feed 1 ambient Stream";

    auto feed1_lpp = CreateNoValveOutletStream(CreateVelocityBC(0.0));
    feed1_lpp->name_ = "LPP Stream";

    // Feed 2
    Step& step_feed2 = data.cycle.GetStep(feed2);
    auto feed2_inlet = CreateNoValveInletStream(CreateVelocityBC(feed_vel), inlet_mf, feed_temp);
    feed2_inlet->name_ = "Feed 2 inlet";
    auto feed2_in_bc = std::make_unique<FlowBoundaryCondition>();
    feed2_in_bc->direction_ = BoundaryDirection::INLET;
    feed2_in_bc->location_  = BoundaryLocation::BOTTOM;
    feed2_in_bc->streams_.push_back(std::move(feed2_inlet));
    step_feed2.bottom_boundary_condition_ = std::move(feed2_in_bc);

    auto feed2_ambient = CreateValveOutletStream(CreateConstantBoundaryValue(ambient_pressure), ambient_valve_cv);
    feed2_ambient->name_ = "Feed 2 ambient Stream";
    auto feed2_out_bc = std::make_unique<FlowBoundaryCondition>();
    feed2_out_bc->direction_ = BoundaryDirection::OUTLET;
    feed2_out_bc->location_  = BoundaryLocation::TOP;
    feed2_out_bc->streams_.push_back(std::move(feed2_ambient));
    step_feed2.top_boundary_condition_ = std::move(feed2_out_bc);

    // Pressurisation top inlet – pressure BC connected to LPP stream
    auto rep_stream = CreatePressureInletBoundary(
        CreateLinearBoundaryValue(0.03e5, 101325.0, press_time),
        { {CO2, 0.0}, {N2, 1.0} },
        298.0
    );
    rep_stream->location_  = BoundaryLocation::TOP;
    rep_stream->direction_ = BoundaryDirection::INLET;

    // Wire up the connection
    auto endA = MakeEndpoint(data.cycle.GetStep(press), *rep_stream);
    auto endB = MakeEndpoint(data.cycle.GetStep(feed1), *feed1_lpp);
    MakeMutualConnection(endA, endB);

    // Commit Feed 1 top outlets (ambient + LPP)
    auto feed1_out_bc = std::make_unique<FlowBoundaryCondition>();
    feed1_out_bc->direction_ = BoundaryDirection::OUTLET;
    feed1_out_bc->location_  = BoundaryLocation::TOP;
    feed1_out_bc->streams_.push_back(std::move(feed1_ambient));
    feed1_out_bc->streams_.push_back(std::move(feed1_lpp));
    step_feed1.top_boundary_condition_ = std::move(feed1_out_bc);

    // Blowdown
    Step& step_blowdown = data.cycle.GetStep(blowdown);
    auto bd_pressure = std::make_unique<PressureBoundaryCondition>();
    bd_pressure->direction_ = BoundaryDirection::OUTLET;
    bd_pressure->location_  = BoundaryLocation::TOP;
    bd_pressure->pressure_  = CreateExponentialBoundaryValue(0.072e5, blowdown_time);
    bd_pressure->mole_fraction_[CO2] = 0.0;
    bd_pressure->mole_fraction_[N2]  = 1.0;
    step_blowdown.top_boundary_condition_ = std::move(bd_pressure);
    auto bd_closed = std::make_unique<ClosedBoundaryCondition>();
    bd_closed->location_ = BoundaryLocation::BOTTOM;
    step_blowdown.bottom_boundary_condition_ = std::move(bd_closed);

    // Evacuation
    Step& step_evac = data.cycle.GetStep(evacuation);
    auto evac_closed = std::make_unique<ClosedBoundaryCondition>();
    evac_closed->location_ = BoundaryLocation::TOP;
    step_evac.top_boundary_condition_ = std::move(evac_closed);
    auto evac_pressure = std::make_unique<PressureBoundaryCondition>();
    evac_pressure->direction_ = BoundaryDirection::OUTLET;
    evac_pressure->location_  = BoundaryLocation::BOTTOM;
    evac_pressure->pressure_  = CreateExponentialBoundaryValue(0.03e5, evac_time);
    evac_pressure->mole_fraction_[CO2] = 0.0;
    evac_pressure->mole_fraction_[N2]  = 1.0;
    step_evac.bottom_boundary_condition_ = std::move(evac_pressure);

    // Pressurisation
    Step& step_press = data.cycle.GetStep(press);
    step_press.top_boundary_condition_ = std::move(rep_stream);
    auto press_closed = std::make_unique<ClosedBoundaryCondition>();
    press_closed->location_ = BoundaryLocation::BOTTOM;
    step_press.bottom_boundary_condition_ = std::move(press_closed);

    ConfigureCycle(5);
    ConfigureKPI(0, {0, 1}, {BoundaryLocation::BOTTOM, BoundaryLocation::BOTTOM}, {3}, {BoundaryLocation::BOTTOM});
    data.settings.max_step = 2.0;

    CoupledSolver solver{data};
    std::string error{};
    ASSERT_TRUE(data.PerformChecks(error)) << "System checks failed: " << error;
    int retval = solver.IntegrateCycle(data.cycle, false);
    ASSERT_GE(retval, 0) << "Solver did not complete";

    // ---- validate connected-stream consistency on the last cycle ----
    const auto& b_log = solver.GetBoundaryLog();
    const int last_cycle = data.cycle.number_of_cycles_ - 1;

    // Step indices (insertion order)
    constexpr int feed1_idx = 0;
    constexpr int press_idx = 4;

    // Names under which the two sides are logged.
    // Flow streams  : stream->name_
    // Pressure BCs  : "<step_name> Pressure Boundary <Top|Bottom>"
    const std::string lpp_name      = "LPP Stream";
    const std::string press_bc_name = press + " Pressure Boundary Top";

    // Helper: normalise a time vector so it starts at 0
    auto normalize = [](std::vector<double> t) {
      if (!t.empty()) {
        double t0 = t.front();
        for (auto& v : t) v -= t0;
      }
      return t;
    };

    // Helper: linear interpolation
    auto interp = [](const std::vector<double>& xp,
                     const std::vector<double>& yp,
                     double x) -> double {
      if (x <= xp.front()) return yp.front();
      if (x >= xp.back())  return yp.back();
      auto it = std::lower_bound(xp.begin(), xp.end(), x);
      size_t i = std::distance(xp.begin(), it);
      double alpha = (x - xp[i-1]) / (xp[i] - xp[i-1]);
      return yp[i-1] + alpha * (yp[i] - yp[i-1]);
    };

    // ---- CO2 mole fraction ----
    auto lpp_mf   = b_log.GetComponentVariableData(last_cycle, feed1_idx, lpp_name,
                        BoundaryComponentVariableType::MOLE_FRACTION, CO2);
    auto press_mf = b_log.GetComponentVariableData(last_cycle, press_idx, press_bc_name,
                        BoundaryComponentVariableType::MOLE_FRACTION, CO2);

    ASSERT_FALSE(lpp_mf.variable_timeseries.empty())   << "LPP stream CO2 mole fraction not logged";
    ASSERT_FALSE(press_mf.variable_timeseries.empty()) << "Pressurisation BC CO2 mole fraction not logged";

    {
      auto t_lpp   = normalize(lpp_mf.time);
      auto t_press = normalize(press_mf.time);

      // Sanity: the LPP stream should carry some N2-rich gas (y_CO2 < 0.15)
      double mean_lpp_co2 = 0.0;
      for (double v : lpp_mf.variable_timeseries) mean_lpp_co2 += v;
      mean_lpp_co2 /= static_cast<double>(lpp_mf.variable_timeseries.size());
      EXPECT_LT(mean_lpp_co2, 0.15) << "LPP stream should be N2-enriched (light product)";

      // The pressurisation inlet must receive the same composition as the LPP outlet.
      size_t n_compared = 0;
      for (size_t i = 0; i < t_lpp.size(); ++i) {
        // t=0 is logged before the step is advanced, so compare only t>0.
        if (t_lpp[i] <= 1e-12) continue;
        double y_lpp   = lpp_mf.variable_timeseries[i];
        double y_press = interp(t_press, press_mf.variable_timeseries, t_lpp[i]);
        EXPECT_NEAR(y_lpp, y_press, 5e-3)
            << "CO2 mole fraction mismatch at step-local t=" << t_lpp[i]
            << "s  (LPP=" << y_lpp << ", Pressurisation BC=" << y_press << ")";
        n_compared++;
      }
      EXPECT_GT(n_compared, static_cast<size_t>(0)) << "No t>0 samples were compared for mole fraction";
    }

    // ---- Fluid temperature ----
    auto lpp_temp   = b_log.GetVariableData(last_cycle, feed1_idx, lpp_name,
                          BoundaryVariableType::TEMPERATURE);
    auto press_temp = b_log.GetVariableData(last_cycle, press_idx, press_bc_name,
                          BoundaryVariableType::TEMPERATURE);

    ASSERT_FALSE(lpp_temp.variable_timeseries.empty())   << "LPP stream temperature not logged";
    ASSERT_FALSE(press_temp.variable_timeseries.empty()) << "Pressurisation BC temperature not logged";

    {
      auto t_lpp   = normalize(lpp_temp.time);
      auto t_press = normalize(press_temp.time);

      size_t n_compared = 0;
      for (size_t i = 0; i < t_lpp.size(); ++i) {
        if (t_lpp[i] <= 1e-12) continue;
        double T_lpp   = lpp_temp.variable_timeseries[i];
        double T_press = interp(t_press, press_temp.variable_timeseries, t_lpp[i]);
        EXPECT_NEAR(T_lpp, T_press, 1.0)
            << "Temperature mismatch at step-local t=" << t_lpp[i]
            << "s  (LPP=" << T_lpp << "K, Pressurisation BC=" << T_press << "K)";
        n_compared++;
      }
      EXPECT_GT(n_compared, static_cast<size_t>(0)) << "No t>0 samples were compared for temperature";
    }
  }

  TEST_F(Perez2019Test, ConfigD1) {
    data.info.title = "Perez 2019 - Configuration D1";
    data.info.description = "Five-step pressure swing adsorption cycle for CO2/N2 separation based on Perez et al. 2019, configuration D1. "
                            "This cycle features feed splitting with two separate feed steps, allowing for product stream recycling. "
                            "The cycle uses 13X zeolite adsorbent and includes inter-step connections where the low-pressure product stream "
                            "from the first feed step is recycled and used for bed pressurisation, improving process efficiency. "
                            "The cycle consists of two feed steps, blowdown, evacuation, and pressurisation steps. "
                            "Feed 1 has dual outlets allowing product stream splitting, while Feed 2 operates with a single outlet. "
                            "The simulation demonstrates a more complex PSA cycle configuration with stream connections between steps.";
    data.info.author = "Sam Affleck";
    
    const std::string CO2 = "CarbonDioxide";
    const std::string N2 = "Nitrogen";
    const std::string z13X = "13X";
    
    SetUpPerez2019();

    // Cycle
    const double press_time = 27.5; // [s]
    const double feed1_time = press_time; // [s]
    const double feed2_time = 117; // [s]
    const double blowdown_time = 21; // [s]
    const double evac_time = 167; // [s]

    const std::string feed1 = "Feed 1";
    const std::string feed2 = "Feed 2";
    const std::string blowdown = "Blowdown";
    const std::string evacuation = "Evacuation";
    const std::string press = "Pressurisation";

    AddStep(feed1, feed1_time);
    AddStep(feed2, feed2_time);
    AddStep(blowdown, blowdown_time);
    AddStep(evacuation, evac_time);
    AddStep(press, press_time);

    // FEED1 STEP
    const double interstitial_feed_vel = 0.08;
    const double feed_vel = interstitial_feed_vel * 0.41;
    const double feed_temp = 298;
    const double ambient_pressure = 1.0e5;
    const double ambient_valve_cv = 10.0;
    const std::unordered_map<std::string, double> inlet_mf = {
      {CO2, 0.15 },
      {N2, 0.85 }
    };

    // Step 1: Feed 1
    Step& step_feed1 = data.cycle.GetStep(feed1);
    auto feed1_inlet_stream = CreateNoValveInletStream(
      CreateVelocityBC(feed_vel),
      inlet_mf,
      feed_temp
    );
    feed1_inlet_stream->name_ = "Feed 1 inlet";
    
    auto feed1_flow_bc = std::make_unique<FlowBoundaryCondition>();
    feed1_flow_bc->direction_ = BoundaryDirection::INLET;
    feed1_flow_bc->location_ = BoundaryLocation::BOTTOM;
    feed1_flow_bc->streams_.push_back(std::move(feed1_inlet_stream));
    step_feed1.bottom_boundary_condition_ = std::move(feed1_flow_bc);
    
    auto feed1_ambient_stream = CreateValveOutletStream(
      CreateConstantBoundaryValue(ambient_pressure),
      ambient_valve_cv
    );
    feed1_ambient_stream->name_ = "Feed 1 ambient Stream";
    
    auto feed1_LPP_stream = CreateNoValveOutletStream(
      CreateVelocityBC(0.0)
    );
    feed1_LPP_stream->name_ = "LPP Stream";

    // FEED 2 STEP
    Step& step_feed2 = data.cycle.GetStep(feed2);
    auto feed2_inlet_stream = CreateNoValveInletStream(
      CreateVelocityBC(feed_vel),
      inlet_mf,
      feed_temp
    );
    feed2_inlet_stream->name_ = "Feed 2 inlet";
    
    auto feed2_flow_bc_in = std::make_unique<FlowBoundaryCondition>();
    feed2_flow_bc_in->direction_ = BoundaryDirection::INLET;
    feed2_flow_bc_in->location_ = BoundaryLocation::BOTTOM;
    feed2_flow_bc_in->streams_.push_back(std::move(feed2_inlet_stream));
    step_feed2.bottom_boundary_condition_ = std::move(feed2_flow_bc_in);

    auto feed2_ambient_stream = CreateValveOutletStream(
      CreateConstantBoundaryValue(ambient_pressure),
      ambient_valve_cv
    );
    feed2_ambient_stream->name_ = "Feed 2 ambient Stream";
    
    auto feed2_flow_bc_out = std::make_unique<FlowBoundaryCondition>();
    feed2_flow_bc_out->direction_ = BoundaryDirection::OUTLET;
    feed2_flow_bc_out->location_ = BoundaryLocation::TOP;
    feed2_flow_bc_out->streams_.push_back(std::move(feed2_ambient_stream));
    step_feed2.top_boundary_condition_ = std::move(feed2_flow_bc_out);

    // ==== PRESSURISATION ====
    auto rep_stream = CreatePressureInletBoundary(
      CreateLinearBoundaryValue(0.03e5, 101325, press_time),
      { 
        { CO2, 0.0 },
        { N2, 1.0 } 
      },
      298.0
    );
    rep_stream->location_ = BoundaryLocation::TOP;
    rep_stream->direction_ = BoundaryDirection::INLET;

    // Make connecting streams interact
    // TODO: Need to check how connections work in new architecture
    auto endA = MakeEndpoint(data.cycle.GetStep(press), *rep_stream);
    auto endB = MakeEndpoint(data.cycle.GetStep(feed1), *feed1_LPP_stream);
    MakeMutualConnection(endA, endB);

    auto feed1_flow_bc_out = std::make_unique<FlowBoundaryCondition>();
    feed1_flow_bc_out->direction_ = BoundaryDirection::OUTLET;
    feed1_flow_bc_out->location_ = BoundaryLocation::TOP;
    feed1_flow_bc_out->streams_.push_back(std::move(feed1_ambient_stream));
    feed1_flow_bc_out->streams_.push_back(std::move(feed1_LPP_stream));
    step_feed1.top_boundary_condition_ = std::move(feed1_flow_bc_out);

    // BLOWDOWN
    Step& step_blowdown = data.cycle.GetStep(blowdown);
    auto blowdown_pressure = std::make_unique<PressureBoundaryCondition>();
    blowdown_pressure->direction_ = BoundaryDirection::OUTLET;
    blowdown_pressure->location_ = BoundaryLocation::TOP;
    blowdown_pressure->pressure_ = CreateExponentialBoundaryValue(0.072e5, blowdown_time);
    blowdown_pressure->mole_fraction_[CO2] = 0.0;
    blowdown_pressure->mole_fraction_[N2] = 1.0;
    step_blowdown.top_boundary_condition_ = std::move(blowdown_pressure);
    
    auto blowdown_closed = std::make_unique<ClosedBoundaryCondition>();
    blowdown_closed->location_ = BoundaryLocation::BOTTOM;
    step_blowdown.bottom_boundary_condition_ = std::move(blowdown_closed);

    // EVACUATION
    Step& step_evac = data.cycle.GetStep(evacuation);
    auto evac_closed_top = std::make_unique<ClosedBoundaryCondition>();
    evac_closed_top->location_ = BoundaryLocation::TOP;
    step_evac.top_boundary_condition_ = std::move(evac_closed_top);
    
    auto evac_pressure = std::make_unique<PressureBoundaryCondition>();
    evac_pressure->direction_ = BoundaryDirection::OUTLET;
    evac_pressure->location_ = BoundaryLocation::BOTTOM;
    evac_pressure->pressure_ = CreateExponentialBoundaryValue(0.03e5, evac_time);
    evac_pressure->mole_fraction_[CO2] = 0.0;
    evac_pressure->mole_fraction_[N2] = 1.0;
    step_evac.bottom_boundary_condition_ = std::move(evac_pressure);
    
    // PRESSURISATION
    Step& step_press = data.cycle.GetStep(press);
    step_press.top_boundary_condition_ = std::move(rep_stream);
    
    auto press_closed = std::make_unique<ClosedBoundaryCondition>();
    press_closed->location_ = BoundaryLocation::BOTTOM;
    step_press.bottom_boundary_condition_ = std::move(press_closed);

    // CONFIGURE CYCLE — run to cyclic steady state (Perez 2019: <1% component mass-balance
    // error for 5 consecutive cycles, ~80-150 cycles to CSS). Minimum 100 cycles, cap 250.
    ConfigureCycle(250);
    data.settings.css.method_ = CSSMethod::MASS_BALANCE_COMPONENT;
    data.settings.css.min_cycles_ = 100;

    // KPI: Feed1(0),Feed2(1)@Bottom, Product: Evacuation(3)@Bottom, Product component: CO2(0)
    ConfigureKPI(0, {0, 1}, {BoundaryLocation::BOTTOM, BoundaryLocation::BOTTOM}, {3}, {BoundaryLocation::BOTTOM});

    data.settings.max_step = 2.0;

    CoupledSolver solver{data};

    std::string error{};
    std::string msg{};
    EXPECT_TRUE(data.PerformChecks(error)) << "System checks failed: " << error;

    bool verbose = false;
    int retval = solver.IntegrateCycle(data.cycle, verbose);
    EXPECT_GE(retval, 0) << "Solver did not complete";
    
    // solver.PrintStatistics();

    EXPECT_TRUE(IsMassConserved(solver.GetLog())) << "Mass is not conserved.";
    EXPECT_TRUE(IsMolesConserved(solver.GetLog())) << "Mass is not conserved.";
    EXPECT_TRUE(IsComponentMolesConserved(solver.GetLog())) << "Component moles not conserved.";

    // --- CO2 purity/recovery at CSS vs Perez (2019) Table 3, Exp D1 (Model 96.1% / 86.8%) ---
    // LPP cycle: Feed1(0), Feed2(1), Blowdown(2), Evacuation(3), Pressurisation(4).
    // CO2 product leaves at BOTTOM during Evacuation; CO2 fed at BOTTOM during both feed
    // steps; LPP pressurisation (step 4, top) re-uses light product (internal, not product).
    {
      const EquationLog& eq_log = solver.GetLog();
      const BoundaryLog& b_log = solver.GetBoundaryLog();
      const auto& cyc = eq_log.sys_data.cycle_index_set_;
      int last_cycle = cyc.empty() ? 0 : *cyc.rbegin();
      KPIParams kpi(eq_log, b_log);
      kpi.product_comp_ = "CarbonDioxide";
      kpi.cycle_index_ = last_cycle;
      kpi.feed_endpoints = KpiEndpointsForFaces(data.cycle,
        { { 0, BoundaryLocation::BOTTOM }, { 1, BoundaryLocation::BOTTOM } });
      kpi.product_endpoints = KpiEndpointsForFaces(data.cycle,
        { { 3, BoundaryLocation::BOTTOM } });
      // Validation KPIs follow the papers: purity of the evacuation product, recovery
      // = CO2 in that product / CO2 fed in the feed steps. (Folding the light-product
      // pressurisation stream into the product, as the deck's KPI config does,
      // makes the product-total moles nearly cancel and the purity meaningless.)
      double purity = CalculatePurity(kpi) * 100.0;
      double recovery = CalculateRecovery(kpi) * 100.0;
      std::cout << "\n[KPI] Perez 2019 ConfigD1 | cycles run = " << (last_cycle + 1)
                << " | CO2 purity = " << purity << " % (paper 96.1) | CO2 recovery = "
                << recovery << " % (paper 86.8)\n" << std::endl;
      ValidatePerez("ConfigD1", purity, recovery);
    }

    SaveTest("Perez2019_D1", solver);
    SaveResultsAsExample("Perez2019_D1", solver);
  }

  TEST_F(Perez2019Test, ConfigC1) {
    data.info.title = "Perez 2019 - Configuration C1";
    data.info.description = "Four-step pressure swing adsorption cycle for CO2/N2 separation based on Perez et al. 2019, configuration C1. "
                            "This is a standard PSA cycle configuration using 13X zeolite adsorbent. The cycle consists of feed, blowdown, "
                            "evacuation, and pressurisation steps. During the feed step, a CO2/N2 mixture enters at the bottom of the bed "
                            "and exits at the top. The blowdown step reduces pressure at the top outlet, followed by evacuation at the bottom "
                            "to achieve deep vacuum conditions. Pressurisation occurs at the bottom inlet using the feed gas composition. "
                            "This configuration represents a typical four-step PSA process for gas separation applications.";
    data.info.author = "Sam Affleck";
    
    const std::string CO2 = "CarbonDioxide";
    const std::string N2 = "Nitrogen";
    const std::string z13X = "13X";
    
    SetUpPerez2019();

    // Cycle
    const double feed_time = 114; // [s]
    const double blowdown_time = 19; // [s]
    const double evac_time = 57; // [s]
    const double press_time = 30; // [s]

    const std::string feed = "Feed";
    const std::string blowdown = "Blowdown";
    const std::string evacuation = "Evacuation";
    const std::string press = "Pressurisation";

    AddStep(feed, feed_time);
    AddStep(blowdown, blowdown_time);
    AddStep(evacuation, evac_time);
    AddStep(press, press_time);

    const double interstitial_feed_vel = 0.078;
    const double feed_vel = interstitial_feed_vel * 0.41;
    const double feed_temp = 298;
    const double ambient_pressure = 1.0e5;
    const double bd_pressure = 0.107e5;
    const double evac_pressure = 0.032e5;
    const std::unordered_map<std::string, double> inlet_mf = {
      {CO2, 0.15 },
      {N2, 0.85 }
    };

    // FEED STEP
    Step& step_feed = data.cycle.GetStep(feed);
    auto feed_inlet_stream = CreateNoValveInletStream(
      CreateVelocityBC(feed_vel),
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
    feed_pressure_out->pressure_ = CreateConstantBoundaryValue(ambient_pressure);
    feed_pressure_out->mole_fraction_[CO2] = 0.0;
    feed_pressure_out->mole_fraction_[N2] = 1.0;
    step_feed.top_boundary_condition_ = std::move(feed_pressure_out);

    // PRESSURISATION
    auto rep_stream = CreatePressureInletBoundary(
      CreateLinearBoundaryValue(evac_pressure, ambient_pressure, press_time),
      inlet_mf,
      feed_temp
    );
    rep_stream->location_ = BoundaryLocation::BOTTOM;
    rep_stream->direction_ = BoundaryDirection::INLET;

    // BLOWDOWN
    Step& step_blowdown = data.cycle.GetStep(blowdown);
    auto blowdown_pressure = std::make_unique<PressureBoundaryCondition>();
    blowdown_pressure->direction_ = BoundaryDirection::OUTLET;
    blowdown_pressure->location_ = BoundaryLocation::TOP;
    blowdown_pressure->pressure_ = CreateExponentialBoundaryValue(bd_pressure, blowdown_time);
    blowdown_pressure->mole_fraction_[CO2] = 0.0;
    blowdown_pressure->mole_fraction_[N2] = 1.0;
    step_blowdown.top_boundary_condition_ = std::move(blowdown_pressure);
    
    auto blowdown_closed = std::make_unique<ClosedBoundaryCondition>();
    blowdown_closed->location_ = BoundaryLocation::BOTTOM;
    step_blowdown.bottom_boundary_condition_ = std::move(blowdown_closed);

    // EVACUATION
    Step& step_evac = data.cycle.GetStep(evacuation);
    auto evac_closed_top = std::make_unique<ClosedBoundaryCondition>();
    evac_closed_top->location_ = BoundaryLocation::TOP;
    step_evac.top_boundary_condition_ = std::move(evac_closed_top);
    
    auto evac_pressure_bc = std::make_unique<PressureBoundaryCondition>();
    evac_pressure_bc->direction_ = BoundaryDirection::OUTLET;
    evac_pressure_bc->location_ = BoundaryLocation::BOTTOM;
    evac_pressure_bc->pressure_ = CreateExponentialBoundaryValue(evac_pressure, evac_time);
    evac_pressure_bc->mole_fraction_[CO2] = 0.0;
    evac_pressure_bc->mole_fraction_[N2] = 1.0;
    step_evac.bottom_boundary_condition_ = std::move(evac_pressure_bc);
    
    // PRESSURISATION
    Step& step_press = data.cycle.GetStep(press);
    step_press.bottom_boundary_condition_ = std::move(rep_stream);
    
    auto press_closed = std::make_unique<ClosedBoundaryCondition>();
    press_closed->location_ = BoundaryLocation::TOP;
    step_press.top_boundary_condition_ = std::move(press_closed);

    // CONFIGURE CYCLE — run to cyclic steady state (Perez 2019: <1% component mass-balance
    // error for 5 consecutive cycles, ~80-150 cycles to CSS). Minimum 100 cycles, cap 250.
    ConfigureCycle(250);
    data.settings.css.method_ = CSSMethod::MASS_BALANCE_COMPONENT;
    data.settings.css.min_cycles_ = 100;

    // KPI: Feed(0)@Bottom, Product: Evacuation(2)@Bottom, Product component: CO2(0)
    ConfigureKPI(0, {0, 3}, {BoundaryLocation::BOTTOM, BoundaryLocation::BOTTOM}, {2}, {BoundaryLocation::BOTTOM});

    data.settings.max_step = 2.0;

    CoupledSolver solver{data};

    std::string error{};
    EXPECT_TRUE(data.PerformChecks(error)) << "System checks failed: " << error;

    bool verbose = false;
    int retval = solver.IntegrateCycle(data.cycle, verbose);
    EXPECT_GE(retval, 0) << "Solver did not complete";

    // solver.PrintStatistics();

    EXPECT_TRUE(IsMassConserved(solver.GetLog())) << "Mass is not conserved.";
    EXPECT_TRUE(IsMolesConserved(solver.GetLog())) << "Mass is not conserved.";
    EXPECT_TRUE(IsComponentMolesConserved(solver.GetLog())) << "Component moles not conserved.";

    // --- CO2 purity/recovery at CSS vs Perez (2019) Table 3, Exp C1 (Model 83.4% / 68.2%) ---
    // Basic 4-step cycle: Feed(0), Blowdown(1), Evacuation(2), Pressurisation(3).
    // CO2 product leaves at the BOTTOM during Evacuation; CO2 is fed at the BOTTOM during
    // both Feed and (feed-)Pressurisation. Evaluated at the last (CSS) cycle.
    {
      const EquationLog& eq_log = solver.GetLog();
      const BoundaryLog& b_log = solver.GetBoundaryLog();
      const auto& cyc = eq_log.sys_data.cycle_index_set_;
      int last_cycle = cyc.empty() ? 0 : *cyc.rbegin();
      KPIParams kpi(eq_log, b_log);
      kpi.product_comp_ = "CarbonDioxide";
      kpi.cycle_index_ = last_cycle;
      kpi.feed_endpoints = KpiEndpointsForFaces(data.cycle,
        { { 0, BoundaryLocation::BOTTOM }, { 3, BoundaryLocation::BOTTOM } });
      kpi.product_endpoints = KpiEndpointsForFaces(data.cycle,
        { { 2, BoundaryLocation::BOTTOM } });
      double purity = CalculatePurity(kpi) * 100.0;
      double recovery = CalculateRecovery(kpi) * 100.0;
      std::cout << "\n[KPI] Perez 2019 ConfigC1 | cycles run = " << (last_cycle + 1)
                << " | CO2 purity = " << purity << " % (paper 83.4) | CO2 recovery = "
                << recovery << " % (paper 68.2)\n" << std::endl;
      ValidatePerez("ConfigC1", purity, recovery);
    }

    SaveTest("Perez2019_C1", solver);
    SaveResultsAsExample("Perez2019_C1", solver);
  }

  // Remaining Table 3 experiments via the shared runners.
  // 4-step: A, B (C is ConfigC1 above). LPP: E, F (D is ConfigD1 above).
  TEST_F(Perez2019Test, ConfigA) {
    RunPerez4Step("Perez2019_A", "Config A", 99, 0.207, 0.060e5, 0.032e5, 19, 57, 30, 99.3, 60.0);
  }
  TEST_F(Perez2019Test, ConfigB) {
    RunPerez4Step("Perez2019_B", "Config B", 99, 0.146, 0.072e5, 0.032e5, 19, 57, 30, 95.2, 64.1);
  }
  TEST_F(Perez2019Test, ConfigE) {
    RunPerezLPP("Perez2019_E", "Config E", 140, 0.059, 0.091e5, 0.031e5, 21, 167, 27.5, 90.4, 91.1);
  }
  TEST_F(Perez2019Test, ConfigF) {
    RunPerezLPP("Perez2019_F", "Config F", 117, 0.061, 0.134e5, 0.032e5, 21, 167, 27.5, 79.9, 94.6);
  }

} // namespace adsim
