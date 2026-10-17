// STL includes
#include <filesystem>

// ads includes
#include "adsim/equations/dispersion_model.h"
#include "adsim/equations/heat_transfer_model.h"

// Gtest includes
#include <gtest/gtest.h>

// Test helpers
#include "../fixtures/user_data_test.h"

namespace adsim {

  TEST_F(UserDataTest, BreakthroughTest2) {
    AddComponents({ "Oxygen", "Nitrogen" });

    auto [layer_id, ads_id] = AddLayer("13X");

    auto ssl = std::make_unique<LangmuirIsotherm>();
    ssl->params_["Oxygen"].m1 = 1.5;
    ssl->params_["Oxygen"].b1 = 3e-11;
    ssl->params_["Oxygen"].q1 = 35000;
    AddIsotherm({ "Oxygen" }, ads_id, std::move(ssl));

    // AddLinearDrivingForceModel(ads_id, "Oxygen", 0.5);
    AddLinearDrivingForceModel_MacroporeKnudsen(ads_id, "Oxygen", 1e-9, 3, 0.2);
    AddLinearDrivingForceModelSecondOrder(ads_id, "Nitrogen", 0.5);

    AddConstantInternalHeatTransferCoefficientModel(0, layer_id);
    AddConstantExternalHeatTransferCoefficientModel(0, layer_id);

    AddInitialConditions(101325, 298, 298,
      {
        {"Oxygen", 0.01 },
        {"Nitrogen", 0.99 }
      }
    );

    AddStep("Feed", 500);

    Step& step_feed = data.cycle.GetStep("Feed");
    auto feed_pressure_in = CreatePressureInletBoundary(
      CreateConstantBoundaryValue(111325),
      { {"Oxygen", 0.15},
      {"Nitrogen", 0.85}},
      273
    );
    feed_pressure_in->location_ = BoundaryLocation::TOP;
    feed_pressure_in->direction_ = BoundaryDirection::INLET;
    step_feed.top_boundary_condition_ = std::move(feed_pressure_in);

    auto feed_pressure_out = std::make_unique<PressureBoundaryCondition>();
    feed_pressure_out->direction_ = BoundaryDirection::OUTLET;
    feed_pressure_out->location_ = BoundaryLocation::BOTTOM;
    feed_pressure_out->pressure_ = CreateConstantBoundaryValue(101325);
    feed_pressure_out->mole_fraction_["Oxygen"] = 0.01;
    feed_pressure_out->mole_fraction_["Nitrogen"] = 0.99;
    step_feed.bottom_boundary_condition_ = std::move(feed_pressure_out);

    ConfigureCycle(1);

    CoupledSolver solver{data};

    std::string error{};
    EXPECT_TRUE(data.cycle.PerformChecks(data.fluid, error)) << "System checks failed: " << error;

    bool verbose = false;
    int retval = solver.IntegrateCycle(data.cycle, verbose);
    EXPECT_GE(retval, 0) << "Solver did not complete";

    SaveTest("BreakthroughTest2", solver);

    EXPECT_TRUE(IsMassConserved(solver.GetLog())) << "Total mole fraction is less than 1.";
  }

  TEST_F(UserDataTest, BreakthroughTest3) {
    AddComponents({ "Oxygen", "Nitrogen" });

    auto [layer_id, ads_id] = AddLayer("13X");

    auto& layer = GetLayer(layer_id);
    layer.domain_.length = 1.0;
    layer.eb_ = 0.6;

    auto& ads = GetAdsorbent(ads_id);
    ads.solid_density_ = 1000;

    const double internal_diameter = 1.1284;
    const double external_diameter = 1.1384;
    std::get<CylinderBed>(data.bed.shape_).diameter = internal_diameter;
    data.bed.wall_.thickness_ = (external_diameter - internal_diameter) / 2.0;

    auto iso = std::make_unique<HenryIsotherm>();
    iso->params_["Oxygen"].k = 1e-5;
    iso->params_["Oxygen"].Q = 0;
    AddIsotherm({ "Oxygen" }, ads_id, std::move(iso));

    AddLinearDrivingForceModel(ads_id, "Oxygen", 0.5);

    AddConstantInternalHeatTransferCoefficientModel(0, layer_id);
    AddConstantExternalHeatTransferCoefficientModel(0, layer_id);

    AddInitialConditions(1e5, 298, 298,
      {
        {"Oxygen", 0.0 },
        {"Nitrogen", 1.0 }
      }
    );

    AddStep("Feed", 1000);

    Step& step_feed = data.cycle.GetStep("Feed");
    auto feed_inlet_stream = CreateNoValveInletStream(
      CreateMolarFlowBC(1.0), // 1 mol/s
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
    feed_pressure_out->pressure_ = CreateConstantBoundaryValue(1e5);
    feed_pressure_out->mole_fraction_["Oxygen"] = 0.0;
    feed_pressure_out->mole_fraction_["Nitrogen"] = 1.0;
    step_feed.top_boundary_condition_ = std::move(feed_pressure_out);

    ConfigureCycle(1);

    CoupledSolver solver{data};

    std::string error{};
    EXPECT_TRUE(data.cycle.PerformChecks(data.fluid, error)) << "System checks failed: " << error;

    bool verbose = false;
    int retval = solver.IntegrateCycle(data.cycle, verbose);
    EXPECT_GE(retval, 0) << "Solver did not complete";

    SaveTest("BreakthroughTest3", solver);

    EXPECT_TRUE(IsMassConserved(solver.GetLog())) << "Total mole fraction is less than 1.";
  }

  TEST_F(UserDataTest, ValvePressurisationTest1) {
    AddComponents({ "Oxygen", "Nitrogen" });

    auto [layer_id, ads_id] = AddLayer("13X");

    auto ssl = std::make_unique<LangmuirIsotherm>();
    ssl->params_["Oxygen"].m1 = 1.5;
    ssl->params_["Oxygen"].b1 = 3e-11;
    ssl->params_["Oxygen"].q1 = 35000;
    AddIsotherm({ "Oxygen" }, ads_id, std::move(ssl));

    AddLinearDrivingForceModel(ads_id, "Oxygen", 0.5);
    AddLinearDrivingForceModel(ads_id, "Nitrogen", 0.5);

    AddConstantInternalHeatTransferCoefficientModel(0, layer_id);
    AddConstantExternalHeatTransferCoefficientModel(0, layer_id);

    AddInitialConditions(
      101325, 298, 298,
      {
        {"Oxygen", 0.15 },
        {"Nitrogen", 0.85 }
      }
    );

    AddStep("Press", 300);

    Step& step_press = data.cycle.GetStep("Press");
    auto valve_inlet_stream = CreateValveInletStream(
      CreateConstantBoundaryValue(111325),
      0.03,
      {
        {"Oxygen", 0.15 },
        {"Nitrogen", 0.85 }
      },
      298.0
      // 318.0
    );

    auto press_flow_bc = std::make_unique<FlowBoundaryCondition>();
    press_flow_bc->direction_ = BoundaryDirection::INLET;
    press_flow_bc->location_ = BoundaryLocation::BOTTOM;
    press_flow_bc->streams_.push_back(std::move(valve_inlet_stream));
    step_press.bottom_boundary_condition_ = std::move(press_flow_bc);

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

    SaveTest("ValvePressurisationTest1", solver);

    EXPECT_TRUE(IsMassConserved(solver.GetLog())) << "Total mole fraction is less than 1.";
  }

  TEST_F(UserDataTest, ValveDepressurisationTest1) {
    AddComponents({ "Oxygen", "Nitrogen" });

    auto [layer_id, ads_id] = AddLayer("13X");

    auto ssl = std::make_unique<LangmuirIsotherm>();
    ssl->params_["Oxygen"].m1 = 1.5;
    ssl->params_["Oxygen"].b1 = 3e-11;
    ssl->params_["Oxygen"].q1 = 35000;
    AddIsotherm({ "Oxygen" }, ads_id, std::move(ssl));

    AddLinearDrivingForceModel(ads_id, "Oxygen", 0.5);
    AddLinearDrivingForceModel(ads_id, "Nitrogen", 0.5);

    AddConstantInternalHeatTransferCoefficientModel(0, layer_id);
    AddConstantExternalHeatTransferCoefficientModel(0, layer_id);

    AddInitialConditions(101325, 298, 298,
      { {"Oxygen", 0.15 },
      {"Nitrogen", 0.85 } });

    AddStep("Depress", 300);

    Step& step_depress = data.cycle.GetStep("Depress");
    auto depress_closed = std::make_unique<ClosedBoundaryCondition>();
    depress_closed->location_ = BoundaryLocation::BOTTOM;
    step_depress.bottom_boundary_condition_ = std::move(depress_closed);

    auto valve_outlet_stream = CreateValveOutletStream(
      CreateConstantBoundaryValue(10132.5),
      0.03
    );

    auto depress_flow_bc = std::make_unique<FlowBoundaryCondition>();
    depress_flow_bc->direction_ = BoundaryDirection::OUTLET;
    depress_flow_bc->location_ = BoundaryLocation::TOP;
    depress_flow_bc->streams_.push_back(std::move(valve_outlet_stream));
    step_depress.top_boundary_condition_ = std::move(depress_flow_bc);

    ConfigureCycle(1);

    CoupledSolver solver{data};

    std::string error{};
    EXPECT_TRUE(data.cycle.PerformChecks(data.fluid, error)) << "System checks failed: " << error;

    bool verbose = false;
    int retval = solver.IntegrateCycle(data.cycle, verbose);
    EXPECT_GE(retval, 0) << "Solver did not complete";

    SaveTest("ValveDepressurisationTest1", solver);

    EXPECT_TRUE(IsMassConserved(solver.GetLog())) << "Total mole fraction is less than 1.";
  }

  TEST_F(UserDataTest, DepressurisationTest1) {
    AddComponents({ "Oxygen", "Nitrogen" });

    auto [layer_id, ads_id] = AddLayer("13X");

    auto ssl = std::make_unique<LangmuirIsotherm>();
    ssl->params_["Oxygen"].m1 = 1.5;
    ssl->params_["Oxygen"].b1 = 3e-11;
    ssl->params_["Oxygen"].q1 = 35000;
    AddIsotherm({ "Oxygen" }, ads_id, std::move(ssl));

    AddLinearDrivingForceModel(ads_id, "Oxygen", 0.5);
    AddLinearDrivingForceModel(ads_id, "Nitrogen", 0.5);

    AddConstantInternalHeatTransferCoefficientModel(0, layer_id);
    AddConstantExternalHeatTransferCoefficientModel(0, layer_id);

    AddInitialConditions(101325, 298, 298,
      { {"Oxygen", 0.15 },
      {"Nitrogen", 0.85 } });

    AddStep("Depress", 300);

    Step& step_depress = data.cycle.GetStep("Depress");
    auto depress_closed = std::make_unique<ClosedBoundaryCondition>();
    depress_closed->location_ = BoundaryLocation::BOTTOM;
    step_depress.bottom_boundary_condition_ = std::move(depress_closed);

    auto no_valve_outlet_stream = CreateNoValveOutletStream(
      CreateVelocityBC(0.001)
    );

    auto depress_flow_bc = std::make_unique<FlowBoundaryCondition>();
    depress_flow_bc->direction_ = BoundaryDirection::OUTLET;
    depress_flow_bc->location_ = BoundaryLocation::TOP;
    depress_flow_bc->streams_.push_back(std::move(no_valve_outlet_stream));
    step_depress.top_boundary_condition_ = std::move(depress_flow_bc);

    ConfigureCycle(1);

    CoupledSolver solver{data};

    std::string error{};
    EXPECT_TRUE(data.cycle.PerformChecks(data.fluid, error)) << "System checks failed: " << error;

    bool verbose = false;
    int retval = solver.IntegrateCycle(data.cycle, verbose);
    EXPECT_GE(retval, 0) << "Solver did not complete";

    SaveTest("DepressurisationTest1", solver);

    EXPECT_TRUE(IsMassConserved(solver.GetLog())) << "Total mole fraction is less than 1.";
  }


  TEST_F(UserDataTest, DepressurisationTest2) {
    AddComponents({ "Oxygen", "Nitrogen" });

    auto [layer_id, ads_id] = AddLayer("13X");

    auto ssl = std::make_unique<LangmuirIsotherm>();
    ssl->params_["Oxygen"].m1 = 1.5;
    ssl->params_["Oxygen"].b1 = 3e-11;
    ssl->params_["Oxygen"].q1 = 35000;
    AddIsotherm({ "Oxygen" }, ads_id, std::move(ssl));

    AddLinearDrivingForceModel(ads_id, "Oxygen", 0.5);
    AddLinearDrivingForceModel(ads_id, "Nitrogen", 0.5);

    AddConstantInternalHeatTransferCoefficientModel(0, layer_id);
    AddConstantExternalHeatTransferCoefficientModel(0, layer_id);

    AddInitialConditions(101325, 298, 298,
      { {"Oxygen", 0.15 },
      {"Nitrogen", 0.85 } });

    AddStep("Depress", 300);

    Step& step_depress = data.cycle.GetStep("Depress");
    auto depress_closed = std::make_unique<ClosedBoundaryCondition>();
    depress_closed->location_ = BoundaryLocation::BOTTOM;
    step_depress.bottom_boundary_condition_ = std::move(depress_closed);

    auto depress_pressure_out = std::make_unique<PressureBoundaryCondition>();
    depress_pressure_out->direction_ = BoundaryDirection::OUTLET;
    depress_pressure_out->location_ = BoundaryLocation::TOP;
    depress_pressure_out->mole_fraction_ = {{"Oxygen", 0.15}, {"Nitrogen", 0.85}};
    depress_pressure_out->pressure_ = CreateLinearBoundaryValue(101325, 10132.5, 300);
    step_depress.top_boundary_condition_ = std::move(depress_pressure_out);

    ConfigureCycle(1);

    CoupledSolver solver{data};

    std::string error{};
    EXPECT_TRUE(data.cycle.PerformChecks(data.fluid, error)) << "System checks failed: " << error;

    bool verbose = false;
    int retval = solver.IntegrateCycle(data.cycle, verbose);
    EXPECT_GE(retval, 0) << "Solver did not complete";

    SaveTest("DepressurisationTest2", solver);

    EXPECT_TRUE(IsMassConserved(solver.GetLog())) << "Total mole fraction is less than 1.";
  }




  TEST_F(UserDataTest, FiveStepProductPurge) {
    AddComponents({ "CarbonDioxide", "Nitrogen" });

    auto [layer_id, ads_id] = AddLayer("13X");

    auto dsl = std::make_unique<DualSiteLangmuirIsotherm>();
    dsl->params_["CarbonDioxide"].m1 = 3.16;
    dsl->params_["CarbonDioxide"].b1 = 9.23e-13;
    dsl->params_["CarbonDioxide"].q1 = 47612.86;
    dsl->params_["CarbonDioxide"].m2 = 3.16;
    dsl->params_["CarbonDioxide"].b2 = 4.88e-13;
    dsl->params_["CarbonDioxide"].q2 = 38545.63;
    dsl->params_["Nitrogen"].m1 = 1.2;
    dsl->params_["Nitrogen"].b1 = 5.85e-8;
    dsl->params_["Nitrogen"].q1 = 16059.73;
    dsl->params_["Nitrogen"].m2 = 1.2;
    dsl->params_["Nitrogen"].b2 = 5.85e-8;
    dsl->params_["Nitrogen"].q2 = 16059.73;
    AddIsotherm({ "CarbonDioxide", "Nitrogen" }, ads_id, std::move(dsl));

    AddLinearDrivingForceModel(ads_id, "CarbonDioxide", 0.1);
    AddLinearDrivingForceModel(ads_id, "Nitrogen", 0.1);

    SetHOA(ads_id, "CarbonDioxide", 40000.0);
    SetHOA(ads_id, "Nitrogen", 0.0);

    AddConstantInternalHeatTransferCoefficientModel(59.1, layer_id);
    AddConstantExternalHeatTransferCoefficientModel(0, layer_id);

    auto& layer = GetLayer(layer_id);
    layer.equations_.dispersion_model_ = std::make_unique<DispersionModelRuthven>();
    layer.equations_.heat_dispersion_model_ = std::make_unique<HeatDispersionModelConstant>();

    AddInitialConditions
    (
      101325, 298, 298,
      {
        {"CarbonDioxide", 0.0 },
        {"Nitrogen", 1.0 }
      }
    );

    // Cycle
    AddStep("Feed1", 100);
    AddStep("Feed2", 100);
    AddStep("Blowdown", 100);
    AddStep("Evacuation", 100);
    AddStep("Press", 100);
    AddStep("Purge", 100);

    const double feed_velocity = 0.5;

    // FEED1 STEP
    std::vector<std::unique_ptr<IBoundaryStream>> feed1_streams_in;
    feed1_streams_in.push_back
    (
      CreateNoValveInletStream
      (
        CreateVelocityBC(feed_velocity),
        {
          {"CarbonDioxide", 0.15 },
          {"Nitrogen", 0.85 }
        },
        298.0
      )
    );
    feed1_streams_in[0]->name_ = "Feed 1 inlet";


    std::vector<std::unique_ptr<IBoundaryStream>> feed1_out_streams;
    feed1_out_streams.push_back
    (
      CreateValveOutletStream
      (
        CreateConstantBoundaryValue(101325),
        100.0
      )
    );
    feed1_out_streams[0]->name_ = "Ambient Stream";

    auto feed_press_stream = CreateNoValveOutletStream(
      CreateVelocityBC(0)
    );
    feed_press_stream->name_ = "Feed 1 Press Stream";

    // FEED2 STEP
    std::vector<std::unique_ptr<IBoundaryStream>> feed2_streams_in;
    feed2_streams_in.push_back
    (
      CreateNoValveInletStream
      (
        CreateVelocityBC(feed_velocity),
        {
          {"CarbonDioxide", 0.15 },
          {"Nitrogen", 0.85 }
        },
        298.0
      )
    );
    feed2_streams_in[0]->name_ = "Feed 2 inlet";

    std::vector<std::unique_ptr<IBoundaryStream>> feed2_out_streams;
    feed2_out_streams.push_back
    (
      CreateValveOutletStream
      (
        CreateConstantBoundaryValue(101325),
        100.0
      )
    );
    feed2_out_streams[0]->name_ = "Ambient Stream";

    auto feed_purge_stream = CreateNoValveOutletStream(
      CreateVelocityBC(0)
    );
    feed_purge_stream->name_ = "Feed 2 Purge Stream";

    // ==== BLOWDOWN ====
    std::vector<std::unique_ptr<IBoundaryStream>> blow_streams_out;
    auto blow_out = CreateValveOutletStream(
      CreateConstantBoundaryValue(4e4),
      1
    );
    blow_out->name_ = "Blowdown waste";
    blow_streams_out.push_back(std::move(blow_out));

    // EVACUATION STEP
    std::vector<std::unique_ptr<IBoundaryStream>> evac_streams_out;
    auto evac_open = CreateValveOutletStream(
      CreateConstantBoundaryValue(1e4),
      1
    );
    evac_open->name_ = "Product stream";
    evac_streams_out.push_back(std::move(evac_open));

    // ==== Pressurisation ====
    auto rep_stream = CreatePressureInletBoundary(
      CreateExponentialBoundaryValue(101325, 100),
      { {"CarbonDioxide", 0.0 },
      {"Nitrogen", 1.0 } },
      298.0
    );
    rep_stream->location_ = BoundaryLocation::TOP;
    rep_stream->direction_ = BoundaryDirection::INLET;

    // ==== PURGE ====
    std::vector<std::unique_ptr<IBoundaryStream>> purge_streams_in;
    auto purge_stream = CreateNoValveInletStream(
      CreateMolarFlowBC(0),
      {
        {"CarbonDioxide", 0.0 },
        {"Nitrogen", 1.0 }
      },
      298.0
    );
    purge_stream->name_ = "Purge stream";

    // Make connecting streams interact
    auto endA = MakeEndpoint(data.cycle.GetStep("Press"), *rep_stream);
    auto endB = MakeEndpoint(data.cycle.GetStep("Feed1"), *feed_press_stream);
    MakeMutualConnection(endA, endB);

    ASSERT_TRUE(endA.conn->IsConnected());
    ASSERT_TRUE(endB.conn->IsConnected());

    auto endC = MakeEndpoint(data.cycle.GetStep("Purge"), *purge_stream);
    auto endD = MakeEndpoint(data.cycle.GetStep("Feed2"), *feed_purge_stream);
    MakeMutualConnection(endC, endD);

    ASSERT_TRUE(endC.conn->IsConnected());
    ASSERT_TRUE(endD.conn->IsConnected());

    feed2_out_streams.push_back(
      std::move(feed_purge_stream)
    );

    purge_streams_in.push_back(
      std::move(purge_stream)
    );

    feed1_out_streams.push_back(
      std::move(feed_press_stream)
    );

    // Add all boundary condition streams
    AddFlowBC(std::move(feed1_out_streams), BoundaryLocation::TOP, BoundaryDirection::OUTLET, "Feed1");
    AddFlowBC(std::move(feed1_streams_in), BoundaryLocation::BOTTOM, BoundaryDirection::INLET, "Feed1");

    AddFlowBC(std::move(feed2_out_streams), BoundaryLocation::TOP, BoundaryDirection::OUTLET, "Feed2");
    AddFlowBC(std::move(feed2_streams_in), BoundaryLocation::BOTTOM, BoundaryDirection::INLET, "Feed2");

    AddFlowBC(std::move(blow_streams_out), BoundaryLocation::TOP, BoundaryDirection::OUTLET, "Blowdown");
    AddClosedBC(BoundaryLocation::BOTTOM, "Blowdown");

    AddClosedBC(BoundaryLocation::TOP, "Evacuation");
    AddFlowBC(std::move(evac_streams_out), BoundaryLocation::BOTTOM, BoundaryDirection::OUTLET, "Evacuation");

    AddBC(std::move(rep_stream), BoundaryLocation::TOP, "Press");
    AddClosedBC(BoundaryLocation::BOTTOM, "Press");

    AddFlowBC(std::move(purge_streams_in), BoundaryLocation::TOP, BoundaryDirection::INLET, "Purge");
    AddPressureOutletBC(
      CreateConstantBoundaryValue(
        101325.0
      ),
      BoundaryLocation::BOTTOM,
      "Purge",
      {
        {"CarbonDioxide", 0.0 },
        {"Nitrogen", 1.0 }
      }
    );

    ConfigureCycle(5);

    std::string error{};
    std::string msg{};
    CoupledSolver solver{data};

    EXPECT_TRUE(data.cycle.PerformChecks(data.fluid, error)) << "System checks failed: " << error;

    bool verbose = false;
    int retval = solver.IntegrateCycle(data.cycle, verbose);
    EXPECT_GE(retval, 0) << "Solver did not complete";

    SaveTest("FiveStepProductPurge", solver);

    EXPECT_TRUE(IsMassConserved(solver.GetLog())) << "Total mole fraction is less than 1.";

    // TODO: Below needs adapting for interacting steps.
    // EXPECT_TRUE(CheckMassAccumulation(solver.GetLog())) << "Mass is not conserved.";
  }

  TEST_F(UserDataTest, ExtendedLangmuirIsothermBreakthrough2) {
    AddComponents({ "CarbonDioxide", "Nitrogen" });

    auto [layer_id, ads_id] = AddLayer("13X");

    auto ssl = std::make_unique<ExtendedLangmuirIsotherm>();
    ssl->params_["CarbonDioxide"].m1 = 3.16;
    ssl->params_["CarbonDioxide"].b1 = 9.23e-13;
    ssl->params_["CarbonDioxide"].q1 = 47612.86;
    ssl->params_["Nitrogen"].m1 = 1.2;
    ssl->params_["Nitrogen"].b1 = 5.85e-8;
    ssl->params_["Nitrogen"].q1 = 16059.73;
    AddIsotherm({ "CarbonDioxide", "Nitrogen" }, ads_id, std::move(ssl));

    AddLinearDrivingForceModel(ads_id, "CarbonDioxide", 0.1);
    AddLinearDrivingForceModel(ads_id, "Nitrogen", 0.1);

    SetHOA(ads_id, "CarbonDioxide", 40000.0);
    SetHOA(ads_id, "Nitrogen", 0.0);

    AddConstantInternalHeatTransferCoefficientModel(59.1, layer_id);
    AddConstantExternalHeatTransferCoefficientModel(0, layer_id);

    auto& layer = GetLayer(layer_id);
    layer.equations_.dispersion_model_ = std::make_unique<DispersionModelRuthven>();
    layer.equations_.heat_dispersion_model_ = std::make_unique<HeatDispersionModelConstant>();

    AddInitialConditions
    (
      101325, 298, 298,
      {
        {"CarbonDioxide", 0.0 },
        {"Nitrogen", 1.0 }
      }
    );

    // ==== FEED STEP ====
    AddStep("Feed", 1000);

    std::vector<std::unique_ptr<IBoundaryStream>> feed_streams_in;
    feed_streams_in.push_back
    (
      CreateNoValveInletStream
      (
        CreateVelocityBC(0.5),
        {
          {"CarbonDioxide", 0.15 },
          {"Nitrogen", 0.85 }
        },
        298.0
      )
    );
    AddFlowBC(std::move(feed_streams_in), BoundaryLocation::BOTTOM, BoundaryDirection::INLET, "Feed");

    AddPressureOutletBC
    (
      CreateConstantBoundaryValue(101325),
      BoundaryLocation::TOP,
      "Feed",
      {
        {"CarbonDioxide", 0.0 },
        {"Nitrogen", 1.0 }
      }
    );

    ConfigureCycle(1);

    std::string error{};
    std::string msg{};
    CoupledSolver solver{data};

    EXPECT_TRUE(data.cycle.PerformChecks(data.fluid, error)) << "System checks failed: " << error;

    bool verbose = false;
    int retval = solver.IntegrateCycle(data.cycle, verbose);
    EXPECT_GE(retval, 0) << "Solver did not complete";

    SaveTest("ExtendedLangmuirIsothermBreakthrough2", solver);

    EXPECT_TRUE(IsMassConserved(solver.GetLog())) << "Total mole fraction is less than 1.";
  }

} // end adsim namespace
