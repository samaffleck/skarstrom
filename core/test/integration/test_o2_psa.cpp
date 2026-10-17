// STL includes
#define _USE_MATH_DEFINES
#include <cmath>
#include <filesystem>

// AdSim includes
#include "../src/pch.h"
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

// Test helpers
#include "../fixtures/user_data_test.h"


namespace adsim {

  TEST_F(UserDataTest, PressurisationTest1) {
    AddComponents({ "Oxygen", "Nitrogen" });

    auto [layer_id, ads_id] = AddLayer("13X");

    auto ssl = std::make_unique<LangmuirIsotherm>();
    ssl->params_["Oxygen"].m1 = 0.0;
    ssl->params_["Oxygen"].b1 = 3e-11;
    ssl->params_["Oxygen"].q1 = 35000;
    AddIsotherm({ "Oxygen" }, ads_id, std::move(ssl));

    AddLinearDrivingForceModel(ads_id, "Oxygen", 0.5);
    AddLinearDrivingForceModel(ads_id, "Nitrogen", 0.5);

    AddConstantInternalHeatTransferCoefficientModel(0, layer_id);
    AddConstantExternalHeatTransferCoefficientModel(0, layer_id);

    data.bed.ambient_temperature_ = 298.0;

    auto& layer = GetLayer(layer_id);
    layer.domain_.n_cells = 20;
    layer.eb_ = 0.1;
    layer.domain_.length = 1.0;
    std::get<CylinderBed>(data.bed.shape_).diameter = 0.1;
    GetAdsorbent(ads_id).ep_ = 0.0;

    auto const_disp = std::make_unique<DispersionModelConstant>();
    const_disp->dispersion_coefficient_ = 0.0;
    layer.equations_.dispersion_model_ = std::move(const_disp);

    auto const_heat_disp = std::make_unique<HeatDispersionModelConstant>();
    const_heat_disp->dispersion_coefficient_ = 0.0;
    layer.equations_.heat_dispersion_model_ = std::move(const_heat_disp);

    double initial_pressure = 1.0e5; 

    AddInitialConditions(initial_pressure, 298, 298,
      { 
        {"Oxygen", 0.5 },
        {"Nitrogen", 0.5 } 
      }
    );

    AddStep("Press", 300);
   
    Step& step_press = data.cycle.GetStep("Press");
    auto press_inlet = CreatePressureInletBoundary(
      CreateLinearBoundaryValue(initial_pressure, 10e5, 300),
      {
        {"Oxygen", 0.0 },
        {"Nitrogen", 1.0 }
      },
      298
    );
    press_inlet->location_ = BoundaryLocation::BOTTOM;
    press_inlet->direction_ = BoundaryDirection::INLET;
    step_press.bottom_boundary_condition_ = std::move(press_inlet);

    auto press_closed = std::make_unique<ClosedBoundaryCondition>();
    press_closed->location_ = BoundaryLocation::TOP;
    step_press.top_boundary_condition_ = std::move(press_closed);

    ConfigureCycle(1);

    CoupledSolver solver{data};

    std::string error{};
    EXPECT_TRUE(data.cycle.PerformChecks(data.fluid, error)) << "System checks failed: " << error;
    
    bool verbose = false;
    int retval = solver.IntegrateCycle(data.cycle, verbose);
    EXPECT_GE(retval, 0) << "Solver did not complete";

    SaveTest("PressurisationTest1", solver);

    AvailableParameters available_parameters;
    available_parameters.UpdateParameterList(data);
    available_parameters.Print();

    OptimizationParameters optimisation_params;
    optimisation_params.Add(available_parameters[0]);
    // First available parameter is the Press step time; keep it away from the
    // degenerate 0 s lower bound so every trial produces a solvable cycle.
    optimisation_params[0].lower_bound = 60.0;
    optimisation_params[0].upper_bound = 600.0;

    auto& cycle = data.cycle;
    auto cost_function = [&solver, &cycle] () -> double {
      solver.IntegrateCycle(cycle, false);
      KPIParams kpi_params(solver.ExtractLog(), solver.ExtractBoundaryLog());
      double purity = CalculatePurity(kpi_params);
      double recovery = CalculateRecovery(kpi_params);
      return -(0.5 * purity + 0.5 * recovery);
    };

    // Run the optimiser against the real solver objective with a small
    // evaluation budget so optimiser-solver integration stays covered.
    OptimizerSettings opt_settings;
    opt_settings.n_init_samples = 2;
    opt_settings.n_iterations = 4;

    GlobalOptimizer bo(optimisation_params, cost_function, opt_settings);
    std::vector<double> best = bo.OptimizeParameters();

    ASSERT_EQ(best.size(), optimisation_params.Size());
    EXPECT_EQ(bo.GetIterationCount(),
              opt_settings.n_init_samples + opt_settings.n_iterations);
    EXPECT_GE(best[0], optimisation_params[0].lower_bound);
    EXPECT_LE(best[0], optimisation_params[0].upper_bound);

    // Re-run at the optimiser's best point: cost must be finite and
    // non-positive since purity and recovery are non-negative.
    optimisation_params[0].SetValue(best[0]);
    double best_cost = cost_function();
    EXPECT_TRUE(std::isfinite(best_cost));
    EXPECT_LE(best_cost, 1e-12);

    EXPECT_TRUE(IsMassConserved(solver.GetLog())) << "Total mole fraction is less than 1.";
    EXPECT_TRUE(IsMolesConserved(solver.GetLog())) << "Mass is not conserved.";
    EXPECT_TRUE(IsComponentMolesConserved(solver.GetLog())) << "Component moles not conserved.";
  }

  TEST_F(UserDataTest, IsothermalMassSource) {
    data.info.title = "Isothermal Mass Source";
    data.info.description = "Isothermal adsorption test with mass source term";
    data.info.author = "Sam Affleck";
    
    AddComponents({ "Oxygen", "Nitrogen" });

    const double temp = 298;
    data.bed.ambient_temperature_ = temp;
    const double internal_diameter = 0.2; //[m]
    const double external_diameter = 0.21;
    data.bed.wall_.thickness_ = (external_diameter - internal_diameter) / 2.0;

    auto [layer_id, ads_id] = AddLayer("13X");
    auto& layer = GetLayer(layer_id);
    layer.domain_.length = 1.0;
    std::get<CylinderBed>(data.bed.shape_).diameter = internal_diameter;
    layer.particle_diameter_ = 0.01;
    layer.eb_ = 0.4;
    layer.domain_.n_cells = 100;

    auto& ads_13X = GetAdsorbent(ads_id);
    ads_13X.solid_density_ = 1000;

    const double area = data.bed.GetArea();

    const auto volume = area * layer.domain_.length;
    const auto mass = ads_13X.solid_density_ * (1 - layer.eb_) * volume;

    // const double k = 1e-4;
    const double k = 0.0;
    const double Q = 0.0;

    auto ssl = std::make_unique<HenryIsotherm>();
    ssl->params_["Oxygen"].k = k;
    ssl->params_["Oxygen"].Q = Q;
    AddIsotherm({ "Oxygen" }, ads_id, std::move(ssl));

    const double pressure = 1e5;
    const auto sat_cap = pressure * k; // [mol/kg]
    const auto sat_moles = sat_cap * mass; // [mol]

    AddLinearDrivingForceModel(ads_id, "Oxygen", 1.5);

    layer.equations_.dispersion_model_ = CreateAxialDispersionModel(AxialDispersionType::CONSTANT);
    auto* disp_model = dynamic_cast<DispersionModelConstant*>(layer.equations_.dispersion_model_.get());
    if (disp_model) {
      disp_model->dispersion_coefficient_ = 0.0;
    }

    AddConstantInternalHeatTransferCoefficientModel(0, layer_id);
    AddConstantExternalHeatTransferCoefficientModel(0, layer_id);
    
    AddInitialConditions(pressure, temp, temp,
    { {"Oxygen", 1.0 },
      {"Nitrogen", 0.0 }
    }, 0);

    double feed_time = 220;
    AddStep("Feed", feed_time);

    double flow_rate = 1.0; // [mol/s]
    Step& step_feed = data.cycle.GetStep("Feed");
    auto feed_inlet_stream = CreateNoValveInletStream(
      CreateMolarFlowBC(flow_rate),
      {
        {"Oxygen", 1.0},
        {"Nitrogen", 0.0}
      },
      temp
    );
    
    auto feed_flow_bc = std::make_unique<FlowBoundaryCondition>();
    feed_flow_bc->direction_ = BoundaryDirection::INLET;
    feed_flow_bc->location_ = BoundaryLocation::BOTTOM;
    feed_flow_bc->streams_.push_back(std::move(feed_inlet_stream));
    step_feed.bottom_boundary_condition_ = std::move(feed_flow_bc);

    auto feed_pressure_out = std::make_unique<PressureBoundaryCondition>();
    feed_pressure_out->direction_ = BoundaryDirection::OUTLET;
    feed_pressure_out->location_ = BoundaryLocation::TOP;
    feed_pressure_out->pressure_ = CreateConstantBoundaryValue(pressure);
    feed_pressure_out->mole_fraction_["Oxygen"] = 1.0;
    feed_pressure_out->mole_fraction_["Nitrogen"] = 0.0;
    step_feed.top_boundary_condition_ = std::move(feed_pressure_out);

    ConfigureCycle(1);
    data.settings.max_step = 1.0;

    CoupledSolver solver{data};

    std::string error{};
    EXPECT_TRUE(data.PerformChecks(error)) << "System checks failed: " << error;
    
    bool verbose = false;
    int retval = solver.IntegrateCycle(data.cycle, verbose);
    EXPECT_GE(retval, 0) << "Solver did not complete";

    SaveTest("IsothermalMassSource", solver);

    const double expected_breakthrough_time = sat_moles / flow_rate;

    EXPECT_TRUE(IsMassConserved(solver.GetLog())) << "Total mole fraction is less than 1.";
    EXPECT_TRUE(IsMolesConserved(solver.GetLog())) << "Mass is not conserved.";
    EXPECT_TRUE(IsComponentMolesConserved(solver.GetLog())) << "Component moles not conserved.";
  }

  TEST_F(UserDataTest, NoMassOrHeatTransfer) {
    data.info.title = "No Mass or Heat Transfer";
    data.info.description = "Simple test case with no mass or heat transfer mechanisms";
    data.info.author = "Sam Affleck";
    
    AddComponents({ "Oxygen", "Nitrogen" });

    const std::string z13X = "13X";
    auto [layer_id, ads_id] = AddLayer(z13X);
    auto& _13x = GetLayer(layer_id);
    std::get<CylinderBed>(data.bed.shape_).diameter = 0.1;
    _13x.domain_.length = 1.0;
    _13x.eb_ = 0.4;

    AddConstantInternalHeatTransferCoefficientModel(0, layer_id);
    AddConstantExternalHeatTransferCoefficientModel(0, layer_id);

    AddInitialConditions(101325, 298, 298,
      { {"Oxygen", 0.0 },
      {"Nitrogen", 1.0 }
      }, 0.378);

    double feed_time = 100;
    AddStep("Feed", feed_time);

    double flow_rate = 1.0; // mol/s
    Step& step_feed = data.cycle.GetStep("Feed");
    auto feed_inlet_stream = CreateNoValveInletStream(
      CreateMolarFlowBC(flow_rate),
      {
        {"Oxygen", 1.0},
        {"Nitrogen", 0.0}
      },
      298
    );
    
    auto feed_flow_bc = std::make_unique<FlowBoundaryCondition>();
    feed_flow_bc->direction_ = BoundaryDirection::INLET;
    feed_flow_bc->location_ = BoundaryLocation::BOTTOM;
    feed_flow_bc->streams_.push_back(std::move(feed_inlet_stream));
    step_feed.bottom_boundary_condition_ = std::move(feed_flow_bc);

    auto feed_pressure_out = std::make_unique<PressureBoundaryCondition>();
    feed_pressure_out->direction_ = BoundaryDirection::OUTLET;
    feed_pressure_out->location_ = BoundaryLocation::TOP;
    feed_pressure_out->pressure_ = CreateConstantBoundaryValue(101325);
    feed_pressure_out->mole_fraction_["Oxygen"] = 0.0;
    feed_pressure_out->mole_fraction_["Nitrogen"] = 1.0;
    step_feed.top_boundary_condition_ = std::move(feed_pressure_out);

    ConfigureCycle(1);
    data.settings.max_step = 1.0;
    
    {
      CoupledSolver solver{data};
  
      std::string error{};
      EXPECT_TRUE(data.PerformChecks(error)) << "System checks failed: " << error;
      
      bool verbose = false;
      int retval = solver.IntegrateCycle(data.cycle, verbose);
      EXPECT_GE(retval, 0) << "Solver did not complete";  
      SaveTest("NoMassOrHeatTransfer", solver);
      EXPECT_TRUE(IsMassConserved(solver.GetLog())) << "Total mole fraction is less than 1.";
      EXPECT_TRUE(IsMolesConserved(solver.GetLog())) << "Mass is not conserved.";
      EXPECT_TRUE(IsComponentMolesConserved(solver.GetLog())) << "Component moles not conserved.";
    }

    {
      GetLayer(layer_id).domain_.n_cells = 70;
      CoupledSolver solver{data};
  
      std::string error{};
      EXPECT_TRUE(data.PerformChecks(error)) << "System checks failed: " << error;
      
      bool verbose = false;
      int retval = solver.IntegrateCycle(data.cycle, verbose);
      EXPECT_GE(retval, 0) << "Repeat solver did not complete";   
    }

  }

  TEST_F(UserDataTest, TwoLayerTest1) {
    data.info.title = "Two Layer Test";
    data.info.description = "Multi-layer adsorption bed with 13X and 5A zeolites";
    data.info.author = "Sam Affleck";
    
    AddComponents({ "Oxygen", "Nitrogen" });

    auto [layer_13x, ads_13x] = AddLayer("13X");
    auto [layer_5a, ads_5a] = AddLayer("5A");
    auto& z13X = GetLayer(layer_13x);
    auto& z5A = GetLayer(layer_5a);

    const double Di = 0.1;

    z13X.domain_.length = 1.0;
    z13X.eb_ = 0.37;
    std::get<CylinderBed>(data.bed.shape_).diameter = Di;

    z5A.domain_.length = 0.3;
    z5A.eb_ = 0.42;

    auto& ads_13X = GetAdsorbent(ads_13x);
    ads_13X.solid_density_ = 1000.0;

    auto& ads_5A = GetAdsorbent(ads_5a);
    ads_5A.solid_density_ = 920.0;

    auto ssl = std::make_unique<LangmuirIsotherm>("Oxygen");
    ssl->params_["Oxygen"].m1 = 1.5;
    ssl->params_["Oxygen"].b1 = 3e-11;
    ssl->params_["Oxygen"].q1 = 35000;
    AddIsotherm({ "Oxygen" }, ads_13x, std::move(ssl));

    auto ssl_5A = std::make_unique<LangmuirIsotherm>("Oxygen");
    ssl_5A->params_["Oxygen"].m1 = 0.5;
    ssl_5A->params_["Oxygen"].b1 = 3e-11;
    ssl_5A->params_["Oxygen"].q1 = 35000;
    AddIsotherm({ "Oxygen" }, ads_5a, std::move(ssl_5A));

    AddLinearDrivingForceModelArrhenius(ads_13x, "Oxygen", 120, 1e4);
    AddLinearDrivingForceModel(ads_13x, "Nitrogen", 0.5);
    AddLinearDrivingForceModel(ads_5a, "Oxygen", 0.5);
    AddLinearDrivingForceModel(ads_5a, "Nitrogen", 0.5);

    AddConstantInternalHeatTransferCoefficientModel(0, layer_13x);
    AddConstantExternalHeatTransferCoefficientModel(0, layer_13x);
    AddConstantInternalHeatTransferCoefficientModel(0, layer_5a);
    AddConstantExternalHeatTransferCoefficientModel(0, layer_5a);

    AddInitialConditions(101325, 298, 298,
      { {"Oxygen", 0.01 }, {"Nitrogen", 0.99 } });

    AddStep("Feed", 500.0);

    Step& step_feed = data.cycle.GetStep("Feed");
    auto pump_feed = std::make_unique<PumpBoundaryStream>();
    pump_feed->name_ = "Pump Feed";
    
    auto pump_power_model = std::make_unique<ConstantPumpPowerModel>();
    pump_power_model->power_ = 10.0;
    pump_feed->power_ = std::move(pump_power_model);
    pump_feed->mole_fraction = {{"Oxygen", 0.15}, {"Nitrogen", 0.85}};
    pump_feed->temperature_ = 318.0;
    
    auto feed_flow_bc = std::make_unique<FlowBoundaryCondition>();
    feed_flow_bc->direction_ = BoundaryDirection::INLET;
    feed_flow_bc->location_ = BoundaryLocation::BOTTOM;
    feed_flow_bc->streams_.push_back(std::move(pump_feed));
    step_feed.bottom_boundary_condition_ = std::move(feed_flow_bc);

    auto feed_pressure_out = std::make_unique<PressureBoundaryCondition>();
    feed_pressure_out->direction_ = BoundaryDirection::OUTLET;
    feed_pressure_out->location_ = BoundaryLocation::TOP;
    feed_pressure_out->pressure_ = CreateConstantBoundaryValue(101325);
    feed_pressure_out->mole_fraction_["Oxygen"] = 0.01;
    feed_pressure_out->mole_fraction_["Nitrogen"] = 0.99;
    step_feed.top_boundary_condition_ = std::move(feed_pressure_out);

    ConfigureCycle(1);

    CoupledSolver solver{data};

    std::string error{};
    EXPECT_TRUE(data.PerformChecks(error)) << "System checks failed: " << error;
    
    bool verbose = false;
    int retval = solver.IntegrateCycle(data.cycle, verbose);
    EXPECT_GE(retval, 0) << "Solver did not complete";

    // solver.SaveLogs();
    // solver.PrintStatistics();

    EXPECT_TRUE(IsMassConserved(solver.GetLog())) << "Mass is not conserved.";
    EXPECT_TRUE(IsMolesConserved(solver.GetLog())) << "Mass is not conserved.";
    EXPECT_TRUE(IsComponentMolesConserved(solver.GetLog())) << "Component moles not conserved.";

    SaveTest("TwoLayerTest1", solver);
  }

  TEST_F(UserDataTest, PressurisationTest_FractionMode) {
    AddComponents({ "Oxygen", "Nitrogen" });

    auto [layer_id, ads_id] = AddLayer("13X");

    auto ssl = std::make_unique<LangmuirIsotherm>();
    ssl->params_["Oxygen"].m1 = 0.0;
    ssl->params_["Oxygen"].b1 = 3e-11;
    ssl->params_["Oxygen"].q1 = 35000;
    AddIsotherm({ "Oxygen" }, ads_id, std::move(ssl));

    AddLinearDrivingForceModel(ads_id, "Oxygen", 0.5);
    AddLinearDrivingForceModel(ads_id, "Nitrogen", 0.5);

    AddConstantInternalHeatTransferCoefficientModel(0, layer_id);
    AddConstantExternalHeatTransferCoefficientModel(0, layer_id);

    data.bed.ambient_temperature_ = 298.0;

    auto& layer = GetLayer(layer_id);
    layer.domain_.n_cells = 20;
    layer.eb_ = 0.1;
    layer.domain_.length = 1.0;
    std::get<CylinderBed>(data.bed.shape_).diameter = 0.1;
    GetAdsorbent(ads_id).ep_ = 0.0;

    auto const_disp = std::make_unique<DispersionModelConstant>();
    const_disp->dispersion_coefficient_ = 0.0;
    layer.equations_.dispersion_model_ = std::move(const_disp);

    auto const_heat_disp = std::make_unique<HeatDispersionModelConstant>();
    const_heat_disp->dispersion_coefficient_ = 0.0;
    layer.equations_.heat_dispersion_model_ = std::move(const_heat_disp);

    double initial_pressure = 1.0e5;

    AddInitialConditions(initial_pressure, 298, 298,
      {
        {"Oxygen", 0.5 },
        {"Nitrogen", 0.5 }
      }
    );

    // Use FRACTION mode: total_cycle_time=300s, step_fraction=1.0
    data.cycle.step_time_mode_ = StepTimeMode::FRACTION;
    data.cycle.total_cycle_time_ = 300.0;
    data.cycle.AddStep("Press");
    Step& step_press = data.cycle.GetStep("Press");
    step_press.step_fraction_ = 1.0;

    auto press_inlet = CreatePressureInletBoundary(
      CreateLinearBoundaryValue(initial_pressure, 10e5, 300),
      {
        {"Oxygen", 0.0 },
        {"Nitrogen", 1.0 }
      },
      298
    );
    press_inlet->location_ = BoundaryLocation::BOTTOM;
    press_inlet->direction_ = BoundaryDirection::INLET;
    step_press.bottom_boundary_condition_ = std::move(press_inlet);

    auto press_closed = std::make_unique<ClosedBoundaryCondition>();
    press_closed->location_ = BoundaryLocation::TOP;
    step_press.top_boundary_condition_ = std::move(press_closed);

    ConfigureCycle(1);

    // Verify ResolveStepTimes works correctly
    std::string error{};
    EXPECT_TRUE(data.cycle.PerformChecks(data.fluid, error)) << "System checks failed: " << error;
    EXPECT_DOUBLE_EQ(step_press.step_time_, 300.0) << "Fraction mode did not resolve step time correctly";
    EXPECT_DOUBLE_EQ(data.cycle.GetCycleTime(), 300.0);

    CoupledSolver solver{data};

    bool verbose = false;
    int retval = solver.IntegrateCycle(data.cycle, verbose);
    EXPECT_GE(retval, 0) << "Solver did not complete";

    // Verify parameter registry registers fraction parameters
    AvailableParameters available_parameters;
    available_parameters.UpdateParameterList(data);

    bool found_cycle_time = false;
    bool found_step_fraction = false;
    for (const auto& param : available_parameters) {
      if (param.name == "total cycle time") found_cycle_time = true;
      if (param.name == "step fraction") found_step_fraction = true;
    }
    EXPECT_TRUE(found_cycle_time) << "total_cycle_time_ not registered in parameter registry";
    EXPECT_TRUE(found_step_fraction) << "step_fraction_ not registered in parameter registry";

    EXPECT_TRUE(IsMassConserved(solver.GetLog())) << "Total mole fraction is less than 1.";
    EXPECT_TRUE(IsMolesConserved(solver.GetLog())) << "Mass is not conserved.";
    EXPECT_TRUE(IsComponentMolesConserved(solver.GetLog())) << "Component moles not conserved.";
  }


} // end adsim namespace
