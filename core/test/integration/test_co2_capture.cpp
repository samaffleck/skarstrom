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

  TEST_F(UserDataTest, DAC) {
    data.info.title = "DAC example";
    data.info.description = "Direct air capture CO2 battery.";
    data.info.author = "Sam Affleck";
    
    AddComponents({ "CarbonDioxide", "Nitrogen" });

    const double temp = 293;
    data.bed.ambient_temperature_ = temp;
    const double internal_diameter = 0.16; //[m]
    const double external_diameter = 0.18;
    data.bed.wall_.thickness_ = (external_diameter - internal_diameter) / 2.0;

    double bulk_density = 450;
    double eb = 0.4;
    double solid_density = bulk_density / (1 - eb);

    std::string adsorbent_name = "Lewatit VP OC 1065";
    auto [layer_id, ads_id] = AddLayer(adsorbent_name);
    auto& layer = GetLayer(layer_id);
    layer.domain_.length = 1.65;
    std::get<CylinderBed>(data.bed.shape_).diameter = internal_diameter;
    layer.particle_diameter_ = 0.002;
    layer.eb_ = eb;
    layer.domain_.n_cells = 30;

    auto& ads_13X = GetAdsorbent(ads_id);
    ads_13X.solid_density_ = solid_density;

    const double area = data.bed.GetArea();

    const auto volume = area * layer.domain_.length;
    const auto mass = ads_13X.solid_density_ * (1 - layer.eb_) * volume;
    
    auto toth = std::make_unique<Toth2Isotherm>();
    toth->params_["CarbonDioxide"].q_sat = 3.715;
    toth->params_["CarbonDioxide"].x = 0.0;
    toth->params_["CarbonDioxide"].T_0 = 324.0;
    toth->params_["CarbonDioxide"].b0 = 0.027;
    toth->params_["CarbonDioxide"].dH0 = 108836;
    toth->params_["CarbonDioxide"].w0 = 0.28;
    toth->params_["CarbonDioxide"].alpha = 0.446;
    AddIsotherm({ "CarbonDioxide" }, ads_id, std::move(toth));

    const double pressure = 1.01325e5;

    AddLinearDrivingForceModel(ads_id, "CarbonDioxide", 1.0);
     
    AddInitialConditions(pressure, temp, temp,
    { 
      {"CarbonDioxide", 0.0 },
      {"Nitrogen", 1.0 }
    });

    double feed_time = 3600; // [1h]
    AddStep("Charge", feed_time);

    Step& step_feed = data.cycle.GetStep("Charge");
    auto feed_inlet_stream = CreateNoValveInletStream(
      CreateVelocityBC(0.5),
      {
        {"CarbonDioxide", 400.0e-6},
        {"Nitrogen", (1 - 400.0e-6)}
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
    feed_pressure_out->mole_fraction_["CarbonDioxide"] = 0.0;
    feed_pressure_out->mole_fraction_["Nitrogen"] = 1.0;
    step_feed.top_boundary_condition_ = std::move(feed_pressure_out);

    ConfigureCycle(1);
    CoupledSolver solver{data};

    std::string error{};
    EXPECT_TRUE(data.PerformChecks(error)) << "System checks failed: " << error;
    
    bool verbose = false;
    int retval = solver.IntegrateCycle(data.cycle, verbose);
    EXPECT_GE(retval, 0) << "Solver did not complete";

    SaveTest("DAC", solver);

    EXPECT_TRUE(IsMassConserved(solver.GetLog())) << "Total mole fraction is less than 1.";
    EXPECT_TRUE(IsMolesConserved(solver.GetLog())) << "Mass is not conserved.";
    EXPECT_TRUE(IsComponentMolesConserved(solver.GetLog())) << "Component moles not conserved.";
  }

} // end adsim namespace
