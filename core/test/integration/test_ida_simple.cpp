// Adsim includes
#include "../src/pch.h"
#include "adsim/math/constants.h"
#include "adsim/objects/domain.h"
#include "adsim/objects/porous_media.h"
#include "adsim/solver/ida_settings.h"
#include "adsim/equations/equation.h"
#include "adsim/objects/user_data.h"
#include "adsim/solver/log.h"
#include "adsim/solver/coupled_solver.h"

// GTest includes
#include <gtest/gtest.h>

// STL includes
#include <algorithm>
#include <string>


namespace adsim {

  // Test fixture base class
  class IDASimpleTest : public ::testing::Test {
  protected:
    UserData data;
    
  };
/*
  TEST_F(IDASimpleTest, ConvectionDiffusionSolves) {
    auto& ic = data.ic;
    auto& fluid = data.fluid;
    auto& cycle = data.cycle;
    auto& bed = data.bed;
    auto& dir = data.adsorbent_directory;

    // Create adsorbent material in the directory
    auto ads_id = dir.Add("Default Solid");
    auto& adsorbent = dir.Get(ads_id);

    // Create layer and add adsorbent to its composition
    auto layer_id = bed.layers_.Add("Layer");
    auto& media = bed.layers_.Get(layer_id);
    media.momentum_model_ = MomentumModel::KarmenCozeny;  // fixture predates Ergun: keep its Darcy / Kozeny-Carman pressure drop
    media.AddAdsorbent(ads_id);
    auto& domain = media.domain_;

    domain.n_cells = 30;

    const std::string CO2 = "CarbonDioxide";
    const std::string N2 = "Nitrogen";
    fluid.AddComponent(CO2);
    fluid.AddComponent(N2);

    auto ssl = CreateIsotherm(IsothermType::SINGLE_SITE_LANGMUIR_P, CO2); // Use default params
    adsorbent.adsorption_source_.AddIsotherm({CO2}, std::move(ssl));

    media.equations_.dispersion_model_ = CreateAxialDispersionModel(AxialDispersionType::RUTHVEN);
    media.equations_.heat_dispersion_model_ = CreateAxialHeatDispersionModel(AxialHeatDispersionType::RUTHVEN);

    ic.mole_fraction[CO2] = 0.0;
    ic.mole_fraction[N2] = 1.0;
    ic.pressure = 1e5;
    
    Step step{}; 
    // auto feed_boundary = std::make_unique<PressureBoundaryCondition>();
    // feed_boundary->direction_ = BoundaryDirection::INLET;
    // feed_boundary->location_ = BoundaryLocation::BOTTOM;
    // auto inlet_pressure = std::make_unique<ConstantBoundaryValue>();
    // inlet_pressure->constant_value_ = 1.1e5;
    // feed_boundary->pressure_ = std::move(inlet_pressure);
    // feed_boundary->temperature_ = 293.0;
    // feed_boundary->mole_fraction_ = {
    //   { CO2, 0.15 },
    //   { N2, 0.85 }
    // };
    // step.bottom_boundary_condition_ = std::move(feed_boundary);

    auto feed_boundary = std::make_unique<FlowBoundaryCondition>();
    feed_boundary->direction_ = BoundaryDirection::INLET;
    feed_boundary->location_ = BoundaryLocation::BOTTOM;
    auto inlet_velocity = std::make_unique<VelocityBoundaryValue>();
    inlet_velocity->velocity_ = 0.5;
    auto inlet_boundary = std::make_unique<BoundaryStream>();
    inlet_boundary->direction_ = BoundaryDirection::INLET;
    inlet_boundary->boundary_flow_ = std::move(inlet_velocity);
    inlet_boundary->temperature_ = 293.0;
    inlet_boundary->mole_fraction = {
      { CO2, 0.15 },
      { N2, 0.85 }
    };
    feed_boundary->streams_.push_back(std::move(inlet_boundary));
    step.bottom_boundary_condition_ = std::move(feed_boundary);
    
    auto outlet_boundary = std::make_unique<PressureBoundaryCondition>();
    outlet_boundary->direction_ = BoundaryDirection::OUTLET;
    outlet_boundary->location_ = BoundaryLocation::TOP;
    auto outlet_pressure = std::make_unique<ConstantBoundaryValue>();
    outlet_pressure->constant_value_ = 1.0e5;
    outlet_boundary->pressure_ = std::move(outlet_pressure);
    step.top_boundary_condition_ = std::move(outlet_boundary);

    // auto closed_bottom = std::make_unique<ClosedBoundaryCondition>();
    // closed_bottom->location_ = BoundaryLocation::BOTTOM;
    // step.bottom_boundary_condition_ = std::move(closed_bottom);

    // auto closed_top = std::make_unique<ClosedBoundaryCondition>();
    // closed_top->location_ = BoundaryLocation::TOP;
    // step.top_boundary_condition_ = std::move(closed_top);

    CoupledSolver solver{data};
    std::string error_message{};
    bool checks_successful = solver.Initialise(step, error_message);
    if (!checks_successful) {
      std::cout << error_message << std::endl;
      FAIL() << "Solver failed initialisation!";
    } 

    solver.EnableLogging(true);

    int retval = solver.IntegrateStep(step, true);
    
    EXPECT_GE(retval, 0) << "IDASolve should succeed";
    
    solver.SaveLogs();
    solver.PrintStatistics();
  }

  TEST_F(IDASimpleTest, CycleIntegrationSolves) {
    auto& ic = data.ic;
    auto& fluid = data.fluid;
    auto& cycle = data.cycle;
    auto& bed = data.bed;
    auto& dir = data.adsorbent_directory;

    // Create adsorbent material in the directory
    auto ads_id = dir.Add("Default Solid");
    auto& adsorbent = dir.Get(ads_id);

    // Create layer and add adsorbent to its composition
    auto layer_id = bed.layers_.Add("Layer");
    auto& media = bed.layers_.Get(layer_id);
    media.momentum_model_ = MomentumModel::KarmenCozeny;  // fixture predates Ergun: keep its Darcy / Kozeny-Carman pressure drop
    media.AddAdsorbent(ads_id);
    auto& domain = media.domain_;

    domain.n_cells = 30;

    const std::string CO2 = "CarbonDioxide";
    const std::string N2 = "Nitrogen";
    fluid.AddComponent(CO2);
    fluid.AddComponent(N2);

    auto ssl = CreateIsotherm(IsothermType::SINGLE_SITE_LANGMUIR_P, CO2);
    adsorbent.adsorption_source_.AddIsotherm({CO2}, std::move(ssl));

    media.equations_.dispersion_model_ = CreateAxialDispersionModel(AxialDispersionType::RUTHVEN);
    media.equations_.heat_dispersion_model_ = CreateAxialHeatDispersionModel(AxialHeatDispersionType::RUTHVEN);

    ic.mole_fraction[CO2] = 0.0;
    ic.mole_fraction[N2] = 1.0;
    ic.pressure = 1e5;

    // Configure cycle
    cycle.number_of_cycles_ = 2;

    // Step 1: Flow inlet at bottom (15% CO2, 85% N2) with pressure outlet at top
    Step step1;
    step1.step_name_ = "Step1_FlowUpward";
    step1.step_time_ = 50.0;
    
    auto feed_boundary1 = std::make_unique<FlowBoundaryCondition>();
    feed_boundary1->direction_ = BoundaryDirection::INLET;
    feed_boundary1->location_ = BoundaryLocation::BOTTOM;
    auto inlet_velocity1 = std::make_unique<VelocityBoundaryValue>();
    inlet_velocity1->velocity_ = 0.5;
    auto inlet_stream1 = std::make_unique<BoundaryStream>();
    inlet_stream1->direction_ = BoundaryDirection::INLET;
    inlet_stream1->boundary_flow_ = std::move(inlet_velocity1);
    inlet_stream1->temperature_ = 293.0;
    inlet_stream1->mole_fraction = {
      { CO2, 0.15 },
      { N2, 0.85 }
    };
    feed_boundary1->streams_.push_back(std::move(inlet_stream1));
    step1.bottom_boundary_condition_ = std::move(feed_boundary1);
    
    auto outlet_boundary1 = std::make_unique<PressureBoundaryCondition>();
    outlet_boundary1->direction_ = BoundaryDirection::OUTLET;
    outlet_boundary1->location_ = BoundaryLocation::TOP;
    auto outlet_pressure1 = std::make_unique<ConstantBoundaryValue>();
    outlet_pressure1->constant_value_ = 1.0e5;
    outlet_boundary1->pressure_ = std::move(outlet_pressure1);
    step1.top_boundary_condition_ = std::move(outlet_boundary1);

    // Step 2: Flow inlet at top (0% CO2, 100% N2) with pressure outlet at bottom (opposite direction)
    Step step2;
    step2.step_name_ = "Step2_FlowDownward";
    step2.step_time_ = 50.0;
    
    auto feed_boundary2 = std::make_unique<FlowBoundaryCondition>();
    feed_boundary2->direction_ = BoundaryDirection::INLET;
    feed_boundary2->location_ = BoundaryLocation::TOP;
    auto inlet_velocity2 = std::make_unique<VelocityBoundaryValue>();
    inlet_velocity2->velocity_ = 0.5;
    auto inlet_stream2 = std::make_unique<BoundaryStream>();
    inlet_stream2->direction_ = BoundaryDirection::INLET;
    inlet_stream2->boundary_flow_ = std::move(inlet_velocity2);
    inlet_stream2->temperature_ = 293.0;
    inlet_stream2->mole_fraction = {
      { CO2, 0.0 },
      { N2, 1.0 }
    };
    feed_boundary2->streams_.push_back(std::move(inlet_stream2));
    step2.top_boundary_condition_ = std::move(feed_boundary2);
    
    auto outlet_boundary2 = std::make_unique<PressureBoundaryCondition>();
    outlet_boundary2->direction_ = BoundaryDirection::OUTLET;
    outlet_boundary2->location_ = BoundaryLocation::BOTTOM;
    auto outlet_pressure2 = std::make_unique<ConstantBoundaryValue>();
    outlet_pressure2->constant_value_ = 1.0e5;
    outlet_boundary2->pressure_ = std::move(outlet_pressure2);
    step2.bottom_boundary_condition_ = std::move(outlet_boundary2);

    // Add steps to cycle
    cycle.step_list_.push_back(step1);
    cycle.step_list_.push_back(step2);

    CoupledSolver solver{data};
    solver.EnableLogging(true);

    int retval = solver.IntegrateCycle(cycle, true);
    
    EXPECT_GE(retval, 0) << "IDACycle should succeed";
    
    solver.SaveLogs();
    solver.PrintStatistics();
  }
  
  int main(int argc, char **argv) {
    ::testing::InitGoogleTest(&argc, argv);
    return RUN_ALL_TESTS();
  }
  */
} // End adsim namespace
