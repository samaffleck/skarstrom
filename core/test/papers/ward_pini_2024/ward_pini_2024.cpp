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

// Paper fixture (Ward & Pini 2024) — see README.md in this folder.
#include "../../fixtures/user_data_test.h"
#include "../../fixtures/validation_metrics.h"
#include <limits>

namespace adsim {

  class WardPini2024Test : public UserDataTest {};

  // Validation tolerances vs Table 4 (abs. percentage points; productivity
  // relative), set from the measured error.
  namespace {
    constexpr double kTolPurity = 2.7;  // measured 2.30363 pp
    constexpr double kTolRecovery = 2.1;  // measured 1.79488 pp
    constexpr double kTolProductivity = 0.018;  // measured 0.0071846 (relative)
  }

  TEST_F(WardPini2024Test, FourStepPSA) {
    data.info.title = "Ward 2024";
    data.info.description = "Four-step pressure swing adsorption cycle for CO2/N2 separation based on Ward & Pini 2024. ";
    data.info.author = "Sam Affleck";
    
    const std::string CO2 = "CarbonDioxide";
    const std::string N2 = "Nitrogen";
    const std::string mof_name = "Mg-MOF-74";
    const std::string ac_name = "Activated Carbon";
    
    AddComponents({ CO2, N2 });

    auto mof_id = data.adsorbent_directory.Add(mof_name);
    auto ac_id = data.adsorbent_directory.Add(ac_name);

    auto layer_id = data.bed.layers_.Add("Adsorbent");
    auto& mof = data.adsorbent_directory.Get(mof_id);
    auto& layer = data.bed.layers_.Get(layer_id);
    layer.momentum_model_ = MomentumModel::KarmenCozeny;  // fixture predates Ergun: keep its Darcy / Kozeny-Carman pressure drop

    double length = 1.0;
    double eb = 0.52;
    layer.domain_.length = length;
    layer.particle_diameter_ = 0.002;
    layer.eb_ = eb;
    layer.domain_.n_cells = 30;
    layer.composition_[mof_id] = 1.0;

    double Di = 2 * 0.145;
    double Do = 2 * 0.162;
    double tw = (Do - Di) / 2.0;
    double bulk_density = 437.0;
    double solid_density = bulk_density / (1 - eb);
    mof.solid_density_ = solid_density;
    mof.ep_ = 0.35;
    mof.heat_capacity_ = 900.0;
    mof.thermal_conductivity_ = 0.5;

    std::get<CylinderBed>(data.bed.shape_).diameter = Di;

    data.bed.wall_.heat_capacity_ = 502;
    data.bed.wall_.density_ = 7800;
    data.bed.wall_.thermal_conductivity_ = 16;
    data.bed.wall_.thickness_ = tw;
    data.bed.ambient_temperature_ = 293;

    // Extended (competitive) DSL. Ward & Pini (2024), Sect. 2.1: "The competitive
    // adsorption equilibrium of CO2/N2 on the surface of each adsorbent has been
    // described using the extended dual-site Langmuir (DSL) isotherm model." This
    // deck previously used the non-competitive DualSiteLangmuirIsotherm_Concentration.
    auto dsl = std::make_unique<ExtendedDualSiteLangmuirIsotherm_Concentration>(CO2);
    dsl->params_[CO2].m1 = 6.80;
    dsl->params_[CO2].b1 = 1.81e-7;
    dsl->params_[CO2].q1 = 39300;
    dsl->params_[CO2].m2 = 9.90;
    dsl->params_[CO2].b2 = 1.06e-6;
    dsl->params_[CO2].q2 = 21200;

    dsl->params_[N2].m1 = 14.0;
    dsl->params_[N2].b1 = 3.45e-6;
    dsl->params_[N2].q1 = 15500.0;
    dsl->params_[N2].m2 = 0.0;
    dsl->params_[N2].b2 = 0.0;
    dsl->params_[N2].q2 = 0.0;
    
    AddIsotherm({ CO2, N2 }, mof_id, std::move(dsl));

    AddLinearDrivingForceModel(mof_id, CO2, 1.0);
    AddLinearDrivingForceModel(mof_id, N2, 1.0);

    SetHOA(mof_id, CO2, 30.85e3);
    SetHOA(mof_id, N2, 15.93e3);

    // Adsorbed-phase heat capacity (Table 2: Cp,a = Cp,g = 30.7 J/mol/K). Without this
    // the field defaults to 0 (CONSTANT mode), dropping the adsorbed-phase thermal
    // ballast so the bed overheats and CO2 capacity/recovery collapse.
    data.fluid.adsorbed_heat_capacity_[CO2] = 30.7;
    data.fluid.adsorbed_heat_capacity_[N2] = 30.7;

    // Ward & Pini use the Haghpanah 2013 column energy balance, whose heat of
    // adsorption is the internal-energy convention (the balance carries the extra
    // -Cpa*T*dq/dt adsorbed-phase sensible term). The SetHOA values above are that
    // (larger) quantity, so select the matching convention.
    data.fluid.heat_of_adsorption_convention_ = HeatOfAdsorptionConvention::INTERNAL_ENERGY;

    AddConstantInternalHeatTransferCoefficientModel(8.6, layer_id);
    AddConstantExternalHeatTransferCoefficientModel(2.5, layer_id);

    layer.equations_.dispersion_model_ = CreateAxialDispersionModel(AxialDispersionType::RUTHVEN);

    AddInitialConditions(
      1.0e5, 298, 298,
      {
        {"CarbonDioxide", 0.0 },
        {"Nitrogen", 1.0 }
      }
    );

    // Cycle
    const double feed_time = 57.6; // [s]
    const double blowdown_time = 60.7; // [s]
    const double evac_time = 200; // [s]
    const double press_time = 30; // [s]

    const std::string feed = "Feed";
    const std::string blowdown = "Blowdown";
    const std::string evacuation = "Evacuation";
    const std::string press = "Pressurisation";

    AddStep(feed, feed_time);
    AddStep(blowdown, blowdown_time);
    AddStep(evacuation, evac_time);
    AddStep(press, press_time);

    const double interstitial_feed_vel = 1.53;
    const double feed_vel = interstitial_feed_vel * eb;
    const double feed_temp = 298.15;
    const double ambient_pressure = 1.0e5;
    const double bd_pressure = 0.0905e5;
    const double evac_pressure = 0.01e5;
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

    // CONFIGURE CYCLE
    // Run to cyclic steady state (paper's per-component mass-balance criterion,
    // Wilkins/Rajendran 2022), minimum 5 cycles, hard cap 250.
    ConfigureCycle(250);
    data.settings.css.method_ = CSSMethod::MASS_BALANCE_COMPONENT;
    data.settings.css.min_cycles_ = 5;

    // KPI: Feed(0)@Bottom, Product: Evacuation(2)@Bottom, Product component: CO2(0)
    ConfigureKPI(0, {0, 3}, {BoundaryLocation::BOTTOM, BoundaryLocation::BOTTOM}, {2}, {BoundaryLocation::BOTTOM});

    data.settings.max_step = 2.0; // [s]

    CoupledSolver solver{data};

    std::string error{};
    EXPECT_TRUE(data.PerformChecks(error)) << "System checks failed: " << error;

    bool verbose = false;
    int retval = solver.IntegrateCycle(data.cycle, verbose);
    EXPECT_GE(retval, 0) << "Solver did not complete";

    EXPECT_TRUE(IsMassConserved(solver.GetLog())) << "Mass is not conserved.";
    EXPECT_TRUE(IsMolesConserved(solver.GetLog())) << "Mass is not conserved.";
    EXPECT_TRUE(IsComponentMolesConserved(solver.GetLog())) << "Component moles not conserved.";

    // --- CO2 purity/recovery at cyclic steady state ---
    // Ward & Pini (2024) Table 4, Mg-MOF-74 single-adsorbent, "maximum separation
    // effectiveness": PuCO2 = 98.5 %, ReCO2 = 89.0 %.
    //   purity   = CO2_out(evac) / (CO2+N2)_out(evac)          [product = Evacuation@BOTTOM]
    //   recovery = CO2_out(evac) / [CO2_in(feed) + CO2_in(press)]  [feed = Feed(0)+Press(3)@BOTTOM]
    {
      const EquationLog& eq_log = solver.GetLog();
      const BoundaryLog& b_log = solver.GetBoundaryLog();
      const auto& cyc = eq_log.sys_data.cycle_index_set_;
      int last_cycle = cyc.empty() ? 0 : *cyc.rbegin();
      KPIParams kpi(eq_log, b_log);
      kpi.product_comp_ = CO2;
      kpi.cycle_index_ = last_cycle;
      kpi.feed_endpoints = KpiEndpointsForFaces(data.cycle,
        { { 0, BoundaryLocation::BOTTOM }, { 3, BoundaryLocation::BOTTOM } });
      kpi.product_endpoints = KpiEndpointsForFaces(data.cycle,
        { { 2, BoundaryLocation::BOTTOM } });
      double purity = CalculatePurity(kpi) * 100.0;
      double recovery = CalculateRecovery(kpi) * 100.0;
      std::cout << "\n[KPI] Ward & Pini 2024 Mg-MOF-74 | cycles run = " << (last_cycle + 1)
                << " | CO2 purity = " << purity << " % (paper 98.5)"
                << " | CO2 recovery = " << recovery << " % (paper 89.0)\n" << std::endl;

      // Quantitative comparison with Table 4 (paper's detailed-model optimum).
      // Productivity per Sect. 4.2: CO2 from the evacuation step per bed volume
      // V_bed = L A per total cycle time. Energy (904 kWh/t) is not compared: the
      // deck uses pressure boundary conditions, not vacuum-pump models.
      using namespace validation;
      const double V_bed = M_PI * 0.145 * 0.145 * length;
      const double pr = CalculateComponentMolesInProduct(kpi) / (V_bed * eq_log.sys_data.total_cycle_time);
      CaseReport report("ward_pini_2024", "FourStepPSA");
      report.CheckAbs("CO2 purity", purity, Reference("ward_pini_2024", "FourStepPSA", "CO2_purity"), kTolPurity);
      report.CheckAbs("CO2 recovery", recovery, Reference("ward_pini_2024", "FourStepPSA", "CO2_recovery"), kTolRecovery);
      report.CheckRel("CO2 productivity", pr, Reference("ward_pini_2024", "FourStepPSA", "CO2_productivity"), kTolProductivity);
    }

    SaveTest("Ward2024", solver);
    SaveResultsAsExample("Ward2024", solver);
  }

} // namespace adsim
