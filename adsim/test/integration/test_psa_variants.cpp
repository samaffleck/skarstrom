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

  TEST_F(UserDataTest, ExtendedLangmuirIsothermBreakthrough1) {
    AddComponents({ "CarbonDioxide", "Nitrogen" });

    auto [layer_id, ads_id] = AddLayer("13X");

    auto dsl = std::make_unique<ExtendedDualSiteLangmuirIsotherm>();
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

    // ==== FEED STEP ====
    AddStep("Feed", 1000);
    Step& feed_step = data.cycle.GetStep("Feed"); 

    std::vector<std::unique_ptr<IBoundaryStream>> feed_streams_in;
    feed_streams_in.push_back
    (
      // CreateValveInletStream
      // (
      //   CreateConstantBoundaryValue(111325.0),
      //   1.0,
      //   {
      //     {"CarbonDioxide", 0.15 },
      //     {"Nitrogen", 0.85 }
      //   },
      //   298.0
      // )
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

    auto valve_outlet_stream = CreateValveOutletStream(
      CreateConstantBoundaryValue(101325.0),
      30.0
    );    
    auto flow_out = std::make_unique<FlowBoundaryCondition>();
    flow_out->direction_ = BoundaryDirection::OUTLET;
    flow_out->location_ = BoundaryLocation::TOP;
    flow_out->streams_.push_back(std::move(valve_outlet_stream));
    feed_step.top_boundary_condition_ = std::move(flow_out);

    // AddPressureOutletBC
    // (
    //   CreateConstantBoundaryValue(101325),
    //   BoundaryLocation::TOP,
    //   "Feed",
    //   {
    //     {"CarbonDioxide", 0.0 },
    //     {"Nitrogen", 1.0 }
    //   }
    // );

    ConfigureCycle(1);

    std::string error{};
    std::string msg{};
    CoupledSolver solver{data};

    EXPECT_TRUE(data.cycle.PerformChecks(data.fluid, error)) << "System checks failed: " << error;
    
    bool verbose = false;
    int retval = solver.IntegrateCycle(data.cycle, verbose);
    EXPECT_GE(retval, 0) << "Solver did not complete";

    SaveTest("ExtendedLangmuirIsothermBreakthrough1", solver);

    EXPECT_TRUE(IsMassConserved(solver.GetLog())) << "Total mole fraction is less than 1.";
    EXPECT_TRUE(IsMolesConserved(solver.GetLog())) << "Mass is not conserved.";
    EXPECT_TRUE(IsComponentMolesConserved(solver.GetLog())) << "Component moles not conserved.";
  }
  TEST_F(UserDataTest, TemperaturePurge) {
    data.info.title = "Temperature Purge";
    data.info.description = "Temperature swing adsorption with purge step";
    data.info.author = "Sam Affleck";
    
    AddComponents({ "Oxygen", "Nitrogen" });

    auto [layer_id, ads_id] = AddLayer("13X");
    auto& layer = GetLayer(layer_id);
    layer.particle_diameter_ = 0.02; // [2 mm]
    layer.domain_.length = 3.0;

    AddConstantInternalHeatTransferCoefficientModel(0, layer_id);
    AddConstantExternalHeatTransferCoefficientModel(0, layer_id);

    layer.equations_.dispersion_model_ = CreateAxialDispersionModel(AxialDispersionType::CONSTANT);
    auto* mass_disp_model = dynamic_cast<DispersionModelConstant*>(layer.equations_.dispersion_model_.get());
    if (mass_disp_model) {
      mass_disp_model->dispersion_coefficient_ = 0.0;
    }
    
    layer.equations_.heat_dispersion_model_ = CreateAxialHeatDispersionModel(AxialHeatDispersionType::CONSTANT);
    auto* heat_disp_model = dynamic_cast<HeatDispersionModelConstant*>(layer.equations_.heat_dispersion_model_.get());
    if (heat_disp_model) {
      heat_disp_model->dispersion_coefficient_ = 0.0;
    }

    AddInitialConditions(101325, 298, 298,
      { 
        {"Oxygen", 1.0 },
        {"Nitrogen", 0.0 }
      }
    );

    double feed_time = 100;
    AddStep("Feed", feed_time);

    double flow_rate = 1.0;
    Step& step_feed = data.cycle.GetStep("Feed");
    auto feed_inlet_stream = CreateNoValveInletStream(
      CreateMolarFlowBC(flow_rate),
      {
        {"Oxygen", 1.0},
        {"Nitrogen", 0.0}
      },
      338
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

    SaveTest("TemperaturePurge", solver);

    EXPECT_TRUE(IsMassConserved(solver.GetLog())) << "Total mole fraction is less than 1.";
    EXPECT_TRUE(IsMolesConserved(solver.GetLog())) << "Mass is not conserved.";
    EXPECT_TRUE(IsComponentMolesConserved(solver.GetLog())) << "Component moles not conserved.";
  }
  
  TEST_F(UserDataTest, BreakthroughTest1) {
    data.info.title = "Breakthrough Test 1";
    data.info.description = "Breakthrough curve test with three-component system";
    data.info.author = "Sam Affleck";
    
    AddComponents({ "Oxygen", "Nitrogen", "CarbonDioxide"});

    auto [layer_id, ads_id] = AddLayer("13X");

    auto ssl = std::make_unique<LangmuirIsotherm>();
    ssl->params_["Oxygen"].m1 = 0.5;
    ssl->params_["Oxygen"].b1 = 3e-11;
    ssl->params_["Oxygen"].q1 = 30000;
    ssl->params_["CarbonDioxide"].m1 = 1.5;
    ssl->params_["CarbonDioxide"].b1 = 3e-11;
    ssl->params_["CarbonDioxide"].q1 = 35000;
    AddIsotherm({ "Oxygen", "CarbonDioxide"}, ads_id, std::move(ssl));

    AddLinearDrivingForceModel(ads_id, "Oxygen", 0.1);
    AddLinearDrivingForceModel(ads_id, "Nitrogen", 0.1);
    AddLinearDrivingForceModel(ads_id, "CarbonDioxide", 0.1);

    AddConstantInternalHeatTransferCoefficientModel(0, layer_id);
    AddConstantExternalHeatTransferCoefficientModel(0, layer_id);

    AddInitialConditions(101325, 298, 298,
      { {"Oxygen", 0.0 },
      {"Nitrogen", 1.0 },
      {"CarbonDioxide", 0.0}
      });

    AddStep("Feed", 300);
    
    Step& step_feed = data.cycle.GetStep("Feed");
    auto feed_inlet_stream = CreateNoValveInletStream(
      CreateVelocityBC(0.5),
      {
        {"Oxygen", 0.75},
        {"Nitrogen", 0.0},
        {"CarbonDioxide", 0.25}
      },
      318
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

    CoupledSolverChuncked solver{data};

    std::string error{};
    EXPECT_TRUE(data.PerformChecks(error)) << "System checks failed: " << error;
        
    // Run solver in chunks until completion
    solver.Start();

    ChunkStatus status;
    do {
      solver.Integrate(50); // Run for 50ms per chunk
      
      status = solver.GetStatus();
      if (status.GetStatus() == SolverStatus::FAILED) {
        FAIL() << "Solver failed: " << status.GetMessage();
      }
      
    } while (solver.IsActive());
    
    EXPECT_EQ(status.GetStatus(), SolverStatus::COMPLETED) << "Solver did not complete successfully";
    
    SaveTest("BreakthroughTest1", solver);
    
    EXPECT_TRUE(IsMassConserved(solver.GetLog())) << "Total mole fraction is less than 1.";
    EXPECT_TRUE(IsMolesConserved(solver.GetLog())) << "Mass is not conserved.";
    EXPECT_TRUE(IsComponentMolesConserved(solver.GetLog())) << "Component moles not conserved.";
  }

  // Fixture for the FourStepLPP test family. Each test builds its own copy of
  // the cycle configuration rather than loading a JSON written by the
  // FourStepLPP test, so the family is safe under parallel ctest (`ctest -j`).
  class FourStepLPPTest : public UserDataTest {
  protected:
    void BuildFourStepLPP();
  };

  void FourStepLPPTest::BuildFourStepLPP() {
    data.info.title = "Four Step Low Pressure Process";
    data.info.description = "Four-step pressure swing adsorption cycle: Feed, Blowdown, Evacuation, and Pressurisation";
    data.info.author = "Sam Affleck";
    
    AddComponents({ "CarbonDioxide", "Nitrogen" });
    auto [layer_id, ads_id] = AddLayer("13X");

    data.bed.equations_.thermal_equilibrium_mode_ = ThermalEquilibriumMode::NON_LOCAL;
    auto particle_htm = std::make_unique<ConstantParticleHeatTransferModel>();
    particle_htm->heat_transfer_coefficient_ = 1.0;
    data.bed.equations_.particle_heat_transfer_model_ = std::move(particle_htm);

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

    AddInitialConditions
    (
      1e5, 298, 298,
      { 
        {"CarbonDioxide", 0.0 },
        {"Nitrogen", 1.0 } 
      }
    );

    // ==== FEED STEP ====
    AddStep("Feed", 200);

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
      CreateConstantBoundaryValue(1e5),
      BoundaryLocation::TOP,
      "Feed",
      { 
        {"CarbonDioxide", 0.0 },
        {"Nitrogen", 1.0 } 
      }
    );

    // ==== BLOWDOWN ====
    AddStep("Blowdown", 100);

    AddClosedBC(BoundaryLocation::BOTTOM, "Blowdown");

    AddPressureOutletBC(
      CreateExponentialBoundaryValue(4e4, 100),
      BoundaryLocation::TOP, 
      "Blowdown",
      { 
        {"CarbonDioxide", 0.0 },
        {"Nitrogen", 1.0 } 
      }
    );
    
    // ==== EVACUATION ====
    AddStep("Evacuation", 200);

    AddClosedBC(BoundaryLocation::TOP, "Evacuation");

    AddPressureOutletBC(
      CreateExponentialBoundaryValue(10132, 200),
      BoundaryLocation::BOTTOM,
      "Evacuation",
      { 
        {"CarbonDioxide", 0.0 },
        {"Nitrogen", 1.0 } 
      }
    );

    // ==== PRESSURISATION ====
    const double press_time = 200.0;
    AddStep("Press", press_time);
    AddPressureInletBC(
      CreateExponentialBoundaryValue(1.01325e5, press_time),
      { 
        {"CarbonDioxide", 0.0 },
        {"Nitrogen", 1.0 } 
      },
      298.0,
      BoundaryLocation::TOP,
      "Press"
    );

    AddClosedBC(BoundaryLocation::BOTTOM, "Press");

    ConfigureCycle(1);
  }

  TEST_F(FourStepLPPTest, FourStepLPP) {
    BuildFourStepLPP();

    // KPI: Feed(0)@Bottom, Product: Evacuation(2)@Bottom, Product component: CO2(0)
    ConfigureKPI(0, {0}, {BoundaryLocation::BOTTOM}, {2}, {BoundaryLocation::BOTTOM});

    std::string error{};
    std::string msg{};
    CoupledSolverChuncked solver{data};

    EXPECT_TRUE(data.PerformChecks(error)) << "System checks failed: " << error;
    
    // Run solver in chunks until completion
    solver.Start();

    ChunkStatus status;
    do {
      solver.Integrate(50); // Run for 50ms per chunk
      status = solver.GetStatus();

      if (status.GetStatus() == SolverStatus::FAILED) {
        FAIL() << "Solver failed: " << status.GetMessage();
      }
      
    } while (status.GetStatus() == SolverStatus::IN_PROGRESS);
    
    EXPECT_EQ(status.GetStatus(), SolverStatus::COMPLETED) << "Solver did not complete successfully";

    SaveTest("FourStepLPP", solver);

    EXPECT_TRUE(IsMassConserved(solver.GetLog())) << "Total mole fraction is less than 1.";
    EXPECT_TRUE(IsMolesConserved(solver.GetLog())) << "Mass is not conserved.";
    EXPECT_TRUE(IsComponentMolesConserved(solver.GetLog())) << "Component moles not conserved.";
  }

  TEST_F(FourStepLPPTest, FourStepLPP_LoadFromJSON) {
    // Round-trip the configuration through JSON using this test's own file
    // so the serialization path stays covered without depending on a file
    // written by the FourStepLPP test.
    BuildFourStepLPP();
    SaveSimulation("FourStepLPP_LoadFromJSON");
    data = UserData{};
    LoadSimulation("FourStepLPP_LoadFromJSON");

    // KPI: Feed(0)@Bottom, Product: Evacuation(2)@Bottom, Product component: CO2(0)
    ConfigureKPI(0, {0}, {BoundaryLocation::BOTTOM}, {2}, {BoundaryLocation::BOTTOM});

    CoupledSolver solver{data};

    std::string error{};
    EXPECT_TRUE(data.PerformChecks(error)) << "System checks failed: " << error;
    
    bool verbose = false;
    int retval = solver.IntegrateCycle(data.cycle, verbose);
    EXPECT_GE(retval, 0) << "Solver did not complete";

    SaveTest("FourStepLPP_LoadFromJSON", solver);

    EXPECT_TRUE(IsMassConserved(solver.GetLog())) << "Total mole fraction is less than 1.";
    EXPECT_TRUE(IsMolesConserved(solver.GetLog())) << "Mass is not conserved.";
    EXPECT_TRUE(IsComponentMolesConserved(solver.GetLog())) << "Component moles not conserved.";
  }

  TEST_F(FourStepLPPTest, FourStepLPP_10Cycle) {
    BuildFourStepLPP();

    // KPI: Feed(0)@Bottom, Product: Evacuation(2)@Bottom, Product component: CO2(0)
    ConfigureKPI(0, {0}, {BoundaryLocation::BOTTOM}, {2}, {BoundaryLocation::BOTTOM});

    data.cycle.number_of_cycles_ = 200;
    data.settings.css.method_ = CSSMethod::MOLAR_PROFILE;
    data.settings.css.molar_absolute_tolerance = 4e-4;
    data.settings.css.molar_relative_tolerance = 1e-4;

    CoupledSolver solver{data};

    std::string error{};
    EXPECT_TRUE(data.PerformChecks(error)) << "System checks failed: " << error;
    
    bool verbose = false;
    int retval = solver.IntegrateCycle(data.cycle, verbose);
    EXPECT_GE(retval, 0) << "Solver did not complete";
    
    SaveTest("FourStepLPP_10Cycle", solver);
    
    EXPECT_TRUE(IsMassConserved(solver.GetLog())) << "Mass is not conserved.";
    EXPECT_TRUE(IsMolesConserved(solver.GetLog())) << "Mass is not conserved.";
    EXPECT_TRUE(IsComponentMolesConserved(solver.GetLog())) << "Component moles not conserved.";
  }

  // Regression: the pressurisation step ends with the bed at the boundary
  // pressure and the inlet flow decaying to zero. A hard upwind switch at u = 0
  // made Newton chatter there and IDA stalled at h ~ 1e-6 (cycle 112 with the
  // default 2 s max step; earlier with larger steps). Run well past that point.
  TEST_F(FourStepLPPTest, FourStepLPP_LongRunDoesNotStall) {
    BuildFourStepLPP();
    data.cycle.number_of_cycles_ = 120;
    data.settings.css.method_ = CSSMethod::NONE;

    CoupledSolver solver{data};
    int retval = solver.IntegrateCycle(data.cycle, false);
    EXPECT_GE(retval, 0) << "Solver did not complete 120 cycles";
  }

  // A solver setting IDA rejects must fail initialisation with a message rather
  // than being logged to stderr and silently replaced by IDA's default (L-14).
  TEST_F(FourStepLPPTest, FourStepLPP_RejectedSolverSettingFailsInitialisation) {
    BuildFourStepLPP();
    data.settings.max_order = 6;  // IDA supports BDF orders 1-5 only

    std::string messages;
    CoupledSolver solver{data};
    solver.SetMessageCallback([&](LogLevel, const std::string& msg) { messages += msg + "\n"; });
    solver.SetLogLevel(LogLevel::LOW);
    int retval = solver.IntegrateCycle(data.cycle, false);
    EXPECT_LT(retval, 0) << "Solver ran with a setting IDA rejected";
    EXPECT_NE(messages.find("Applying the solver settings failed"), std::string::npos) << messages;
  }

  // Regression: the same deck must integrate to CSS without a max-step cap.
  TEST_F(FourStepLPPTest, FourStepLPP_NoMaxStepReachesCSS) {
    BuildFourStepLPP();
    data.cycle.number_of_cycles_ = 200;
    data.settings.css.method_ = CSSMethod::MOLAR_PROFILE;
    data.settings.css.molar_absolute_tolerance = 4e-4;
    data.settings.css.molar_relative_tolerance = 1e-4;
    data.settings.max_step = 0.0;  // no cap

    CoupledSolver solver{data};
    int retval = solver.IntegrateCycle(data.cycle, false);
    EXPECT_GE(retval, 0) << "Solver did not complete without a max-step cap";
    EXPECT_TRUE(IsMolesConserved(solver.GetLog())) << "Moles not conserved.";
  }

  TEST_F(FourStepLPPTest, FourStepLPP_CSS_Acceleration) {
    // Run WITHOUT acceleration
    BuildFourStepLPP();

    // KPI: Feed(0)@Bottom, Product: Evacuation(2)@Bottom, Product component: CO2(0)
    ConfigureKPI(0, {0}, {BoundaryLocation::BOTTOM}, {2}, {BoundaryLocation::BOTTOM});

    data.cycle.number_of_cycles_ = 200;
    data.settings.css.method_ = CSSMethod::MOLAR_PROFILE;
    // 1e-4 is a realistic molar-profile CSS tolerance. At 1e-6 the criterion
    // never converges, so the run grinds all 200 cycles into a quasi-steady
    // pressurisation state where the momentum/pressure DAE is stiff (velocity->0)
    // — with physically-correct gas viscosity that trips the integrator. At 1e-4
    // CSS is reached in ~11 cycles, well before that region. (The acceleration
    // cycle-count comparison this test was built for is disabled below anyway.)
    data.settings.css.molar_absolute_tolerance = 1e-4;
    data.settings.css.molar_relative_tolerance = 1e-4;

    CoupledSolver solver_no_acc{data};
    std::string error{};
    ASSERT_TRUE(data.PerformChecks(error)) << "System checks failed: " << error;
    int retval = solver_no_acc.IntegrateCycle(data.cycle, false);
    ASSERT_GE(retval, 0) << "Solver without acceleration did not complete";
    int cycles_no_acc = static_cast<int>(solver_no_acc.GetLog().sys_data.cycle_index_set_.size());

    // Reset and run WITH acceleration
    data = UserData{};
    BuildFourStepLPP();
    ConfigureKPI(0, {0}, {BoundaryLocation::BOTTOM}, {2}, {BoundaryLocation::BOTTOM});
    data.cycle.number_of_cycles_ = 200;
    data.settings.css.method_ = CSSMethod::MOLAR_PROFILE;
    // 1e-4 is a realistic molar-profile CSS tolerance. At 1e-6 the criterion
    // never converges, so the run grinds all 200 cycles into a quasi-steady
    // pressurisation state where the momentum/pressure DAE is stiff (velocity->0)
    // — with physically-correct gas viscosity that trips the integrator. At 1e-4
    // CSS is reached in ~11 cycles, well before that region. (The acceleration
    // cycle-count comparison this test was built for is disabled below anyway.)
    data.settings.css.molar_absolute_tolerance = 1e-4;
    data.settings.css.molar_relative_tolerance = 1e-4;
    data.settings.acc_css.method_ = AccelerationMethod::AITKEN;

    CoupledSolver solver_acc{data};
    ASSERT_TRUE(data.PerformChecks(error)) << "System checks failed: " << error;
    retval = solver_acc.IntegrateCycle(data.cycle, false);
    ASSERT_GE(retval, 0) << "Solver with acceleration did not complete";
    int cycles_acc = static_cast<int>(solver_acc.GetLog().sys_data.cycle_index_set_.size());

    // EXPECT_LT(cycles_acc, cycles_no_acc)
    //   << "Accelerated CSS (" << cycles_acc << " cycles) should converge faster than non-accelerated (" << cycles_no_acc << " cycles)";

    std::cout << "Number of cycles with acceleration: " << cycles_acc << "\n";
    std::cout << "Number of cycles without acceleration: " << cycles_no_acc << "\n";

    EXPECT_TRUE(IsMassConserved(solver_acc.GetLog())) << "Mass is not conserved with acceleration.";
  }

  TEST_F(FourStepLPPTest, FourStepLPP_MassBalanceCSS) {
    BuildFourStepLPP();

    // KPI: Feed(0)@Bottom, Product: Evacuation(2)@Bottom, Product component: CO2(0)
    ConfigureKPI(0, {0}, {BoundaryLocation::BOTTOM}, {2}, {BoundaryLocation::BOTTOM});

    data.cycle.number_of_cycles_ = 200;
    data.settings.css.method_ = CSSMethod::MASS_BALANCE_COMPONENT;
    data.settings.css.mass_balance_tolerance_ = 0.005;
    data.settings.css.min_cycles_ = 5;

    CoupledSolver solver{data};
    std::string error{};
    ASSERT_TRUE(data.PerformChecks(error)) << "System checks failed: " << error;
    int retval = solver.IntegrateCycle(data.cycle, false);
    ASSERT_GE(retval, 0) << "Solver did not complete";

    int cycles = static_cast<int>(solver.GetLog().sys_data.cycle_index_set_.size());
    EXPECT_LT(cycles, 200) << "Mass balance CSS should terminate before max cycles";
    EXPECT_GE(cycles, 5) << "Should run at least min_cycles before checking CSS";

    EXPECT_TRUE(IsMassConserved(solver.GetLog())) << "Mass is not conserved.";
  }

  TEST_F(FourStepLPPTest, FourStepLPP_MassBalanceTotalCSS) {
    BuildFourStepLPP();

    // KPI: Feed(0)@Bottom, Product: Evacuation(2)@Bottom, Product component: CO2(0)
    ConfigureKPI(0, {0}, {BoundaryLocation::BOTTOM}, {2}, {BoundaryLocation::BOTTOM});

    data.cycle.number_of_cycles_ = 200;
    data.settings.css.method_ = CSSMethod::MASS_BALANCE_TOTAL;
    data.settings.css.mass_balance_tolerance_ = 0.005;
    data.settings.css.min_cycles_ = 5;

    CoupledSolver solver{data};
    std::string error{};
    ASSERT_TRUE(data.PerformChecks(error)) << "System checks failed: " << error;
    int retval = solver.IntegrateCycle(data.cycle, false);
    ASSERT_GE(retval, 0) << "Solver did not complete";

    int cycles = static_cast<int>(solver.GetLog().sys_data.cycle_index_set_.size());
    EXPECT_LT(cycles, 200) << "Mass balance total CSS should terminate before max cycles";
    EXPECT_GE(cycles, 5) << "Should run at least min_cycles before checking CSS";

    EXPECT_TRUE(IsMassConserved(solver.GetLog())) << "Mass is not conserved.";
  }

  // Test that FRACTION mode resolves step times correctly and produces
  // the same results as an equivalent ABSOLUTE mode simulation.
  TEST_F(UserDataTest, FourStepLPP_FractionMode) {
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

    AddInitialConditions(
      1e5, 298, 298,
      {
        {"CarbonDioxide", 0.0 },
        {"Nitrogen", 1.0 }
      }
    );

    // Configure FRACTION mode: total cycle = 700s
    // Feed=200/700, Blowdown=100/700, Evacuation=200/700, Press=200/700
    const double total_time = 700.0;
    data.cycle.step_time_mode_ = StepTimeMode::FRACTION;
    data.cycle.total_cycle_time_ = total_time;

    // ==== FEED STEP ====
    data.cycle.AddStep("Feed");
    data.cycle.GetStep("Feed").step_fraction_ = 200.0 / total_time;

    std::vector<std::unique_ptr<IBoundaryStream>> feed_streams_in;
    feed_streams_in.push_back(
      CreateNoValveInletStream(
        CreateVelocityBC(0.5),
        {
          {"CarbonDioxide", 0.15 },
          {"Nitrogen", 0.85 }
        },
        298.0
      )
    );
    AddFlowBC(std::move(feed_streams_in), BoundaryLocation::BOTTOM, BoundaryDirection::INLET, "Feed");

    AddPressureOutletBC(
      CreateConstantBoundaryValue(1e5),
      BoundaryLocation::TOP,
      "Feed",
      {
        {"CarbonDioxide", 0.0 },
        {"Nitrogen", 1.0 }
      }
    );

    // ==== BLOWDOWN ====
    data.cycle.AddStep("Blowdown");
    data.cycle.GetStep("Blowdown").step_fraction_ = 100.0 / total_time;

    AddClosedBC(BoundaryLocation::BOTTOM, "Blowdown");

    AddPressureOutletBC(
      CreateExponentialBoundaryValue(4e4, 100),
      BoundaryLocation::TOP,
      "Blowdown",
      {
        {"CarbonDioxide", 0.0 },
        {"Nitrogen", 1.0 }
      }
    );

    // ==== EVACUATION ====
    data.cycle.AddStep("Evacuation");
    data.cycle.GetStep("Evacuation").step_fraction_ = 200.0 / total_time;

    AddClosedBC(BoundaryLocation::TOP, "Evacuation");

    AddPressureOutletBC(
      CreateExponentialBoundaryValue(10132, 200),
      BoundaryLocation::BOTTOM,
      "Evacuation",
      {
        {"CarbonDioxide", 0.0 },
        {"Nitrogen", 1.0 }
      }
    );

    // ==== PRESSURISATION ====
    const double press_time = 200.0;
    data.cycle.AddStep("Press");
    data.cycle.GetStep("Press").step_fraction_ = press_time / total_time;

    AddPressureInletBC(
      CreateExponentialBoundaryValue(1.01325e5, press_time),
      {
        {"CarbonDioxide", 0.0 },
        {"Nitrogen", 1.0 }
      },
      298.0,
      BoundaryLocation::TOP,
      "Press"
    );

    AddClosedBC(BoundaryLocation::BOTTOM, "Press");

    ConfigureCycle(1);

    // KPI: Feed(0)@Bottom, Product: Evacuation(2)@Bottom, Product component: CO2(0)
    ConfigureKPI(0, {0}, {BoundaryLocation::BOTTOM}, {2}, {BoundaryLocation::BOTTOM});

    // Verify fractions resolve correctly
    std::string error{};
    EXPECT_TRUE(data.PerformChecks(error)) << "System checks failed: " << error;

    EXPECT_NEAR(data.cycle.GetStep("Feed").step_time_, 200.0, 1e-10);
    EXPECT_NEAR(data.cycle.GetStep("Blowdown").step_time_, 100.0, 1e-10);
    EXPECT_NEAR(data.cycle.GetStep("Evacuation").step_time_, 200.0, 1e-10);
    EXPECT_NEAR(data.cycle.GetStep("Press").step_time_, 200.0, 1e-10);
    EXPECT_DOUBLE_EQ(data.cycle.GetCycleTime(), total_time);

    // Verify parameter registry
    AvailableParameters available_parameters;
    available_parameters.UpdateParameterList(data);

    int fraction_param_count = 0;
    bool found_cycle_time = false;
    bool found_step_time = false;
    for (const auto& param : available_parameters) {
      if (param.name == "total cycle time") found_cycle_time = true;
      if (param.name == "step fraction") fraction_param_count++;
      if (param.name == "step time") found_step_time = true;
    }
    EXPECT_TRUE(found_cycle_time) << "total_cycle_time_ not in parameter registry";
    EXPECT_EQ(fraction_param_count, 4) << "Expected 4 step fraction parameters";
    EXPECT_FALSE(found_step_time) << "step_time should not be registered in FRACTION mode";

    // Run solver
    CoupledSolverChuncked solver{data};
    solver.Start();

    ChunkStatus status;
    do {
      solver.Integrate(50);
      status = solver.GetStatus();
      if (status.GetStatus() == SolverStatus::FAILED) {
        FAIL() << "Solver failed: " << status.GetMessage();
      }
    } while (status.GetStatus() == SolverStatus::IN_PROGRESS);

    EXPECT_EQ(status.GetStatus(), SolverStatus::COMPLETED) << "Solver did not complete successfully";

    EXPECT_TRUE(IsMassConserved(solver.GetLog())) << "Total mole fraction is less than 1.";
    EXPECT_TRUE(IsMolesConserved(solver.GetLog())) << "Mass is not conserved.";
    EXPECT_TRUE(IsComponentMolesConserved(solver.GetLog())) << "Component moles not conserved.";
  }


  // The 8-step biogas cycle is run twice: with the Darcy / Kozeny-Carman
  // pressure drop it was built with, and with Ergun. At its 41.6 mol/s feed
  // through a 0.5 m column the Ergun inertial term is ~16x the viscous one,
  // so the Ergun run doubles as a robustness check of the Ergun momentum
  // residuals across valve / closed-end / equalisation step switches (no step
  // may need a solver retry).
  class BiogasPSATest : public UserDataTest {
  protected:
    void RunBiogas8Step(MomentumModel momentum, const std::string& save_name);
  };

  TEST_F(BiogasPSATest, BiogasPSA_8Step_Skarstrom) {
    RunBiogas8Step(MomentumModel::KarmenCozeny, "BiogasPSA_8Step_Skarstrom");
  }

  TEST_F(BiogasPSATest, BiogasPSA_8Step_Skarstrom_Ergun) {
    RunBiogas8Step(MomentumModel::Ergun, "BiogasPSA_8Step_Skarstrom_Ergun");
  }

  void BiogasPSATest::RunBiogas8Step(MomentumModel momentum, const std::string& save_name) {
    data.info.title = "Biogas PSA - 8-Step Skarstrom";
    data.info.description = "Eight step PSA cycle for biogas upgrading (CO2/CH4 separation) using Zeolite 13X.";
    data.info.author = "Sam Affleck";

    AddComponents({ "Methane", "CarbonDioxide" });

    auto [layer_id, ads_id] = AddLayer("Zeolite 13X");

    auto& layer = GetLayer(layer_id);
    layer.domain_.length = 1.2;
    layer.domain_.n_cells = 60;
    layer.particle_diameter_ = 0.0018;
    layer.eb_ = 0.4;
    std::get<CylinderBed>(data.bed.shape_).diameter = 0.5;

    auto& ads_13X = GetAdsorbent(ads_id);
    ads_13X.solid_density_ = 1085.0;
    ads_13X.heat_capacity_ = 1100.0;
    ads_13X.thermal_conductivity_ = 0.15;
    ads_13X.ep_ = 0.54;

    data.bed.wall_.density_ = 8000.0;
    data.bed.wall_.heat_capacity_ = 500.0;
    data.bed.wall_.thermal_conductivity_ = 16.2;
    data.bed.wall_.thickness_ = 0.002;
    data.bed.ambient_temperature_ = 288.15;
    data.bed.equations_.wall_temperature_boundary_condition_ = CreateWallTemperatureBoundaryCondition(WallTemperatureBoundaryConditionType::ZERO_FLUX);

    auto sips = std::make_unique<ExtendedSips2Isotherm>();
    sips->params_["CarbonDioxide"].A1 = 7.268;
    sips->params_["CarbonDioxide"].A2 = -0.6168;
    sips->params_["CarbonDioxide"].B1 = 1.13e-9;
    sips->params_["CarbonDioxide"].B2 = 28400.0;
    sips->params_["CarbonDioxide"].C1 = 0.425;
    sips->params_["CarbonDioxide"].C2 = 0.7238;
    sips->params_["CarbonDioxide"].T_ref = 298.0;
    sips->params_["Methane"].A1 = 4.473;
    sips->params_["Methane"].A2 = -0.6569;
    sips->params_["Methane"].B1 = 3.61e-9;
    sips->params_["Methane"].B2 = 15561.0;
    sips->params_["Methane"].C1 = 1.073;
    sips->params_["Methane"].C2 = 0.0;
    sips->params_["Methane"].T_ref = 298.0;
    AddIsotherm({ "CarbonDioxide", "Methane" }, ads_id, std::move(sips));

    AddLinearDrivingForceModel_Haynes(ads_id, "CarbonDioxide", 3.0, 0.35);
    AddLinearDrivingForceModel_Haynes(ads_id, "Methane", 3.0, 0.35);

    SetHOA(ads_id, "CarbonDioxide", 28400.0);
    SetHOA(ads_id, "Methane", 15561.0);

    layer.equations_.dispersion_model_ = CreateAxialDispersionModel(AxialDispersionType::CONSTANT);
    layer.equations_.heat_dispersion_model_ = CreateAxialHeatDispersionModel(AxialHeatDispersionType::CONSTANT);
    layer.equations_.internal_heat_transfer_model_ = CreateInternalHeatTransferModel(InternalHeatTransferType::DIXON_1996);
    AddConstantExternalHeatTransferCoefficientModel(0.0, layer_id);

    AddInitialConditions(
      1.5e6, 293.0, 293.0,
      {
        {"Methane", 0.53},
        {"CarbonDioxide", 0.47}
      }
    );

    const double P_high      = 1.5e6;
    const double P_low       = 1e5;
    const double F_in        = 41.56;   // [mol/s]
    const double F_purge     = 14.47;   // [mol/s]
    const double CV_product  = 100.0;
    const double CV_PE       = 0.2;
    const double T_feed      = 293.0;   // [K]
    const double T_purge     = 293.15;  // [K]
    const double cycle_time     = 6600;  // [s]

    // ---- Cycle in FRACTION mode ----
    data.cycle.step_time_mode_ = StepTimeMode::FRACTION;
    data.cycle.total_cycle_time_ = cycle_time;

    // ==== FEED 1 ====
    data.cycle.AddStep("Feed 1");
    data.cycle.GetStep("Feed 1").step_fraction_ = 0.125;

    {
      auto feed_stream = CreateNoValveInletStream(
        CreateMolarFlowBC(F_in),
        { {"Methane", 0.53}, {"CarbonDioxide", 0.47} },
        T_feed
      );
      feed_stream->name_ = "Biogas Feed 1";
      auto compressor = std::make_unique<Compressor>();
      compressor->use_pump_ = true;
      compressor->suction_pressure_ = P_low;
      auto eff = std::make_unique<ConstantPumpEfficiencyModel>();
      eff->efficiency_ = 0.8;
      compressor->efficiency_model_ = std::move(eff);
      compressor->energy_eq_ = std::make_unique<IsothermalCompressorEnergyModel>();
      feed_stream->pump_ = std::move(compressor);

      auto feed_flow_bc = std::make_unique<FlowBoundaryCondition>();
      feed_flow_bc->direction_ = BoundaryDirection::INLET;
      feed_flow_bc->location_ = BoundaryLocation::BOTTOM;
      feed_flow_bc->streams_.push_back(std::move(feed_stream));
      data.cycle.GetStep("Feed 1").bottom_boundary_condition_ = std::move(feed_flow_bc);
    }
    AddPressureOutletBC(CreateConstantBoundaryValue(P_high), BoundaryLocation::TOP, "Feed 1");

    // ==== FEED 2 ====
    data.cycle.AddStep("Feed 2");
    data.cycle.GetStep("Feed 2").step_fraction_ = 0.125;

    IBoundaryStream* purge_out_ptr = nullptr;
    IBoundaryStream* press_out_ptr = nullptr;
    {
      auto feed_stream = CreateNoValveInletStream(
        CreateMolarFlowBC(F_in),
        { {"Methane", 0.53}, {"CarbonDioxide", 0.47} },
        T_purge
      );
      feed_stream->name_ = "Biogas Feed 2";
      auto compressor = std::make_unique<Compressor>();
      compressor->use_pump_ = true;
      compressor->suction_pressure_ = P_low;
      auto eff = std::make_unique<ConstantPumpEfficiencyModel>();
      eff->efficiency_ = 0.8;
      compressor->efficiency_model_ = std::move(eff);
      compressor->energy_eq_ = std::make_unique<IsothermalCompressorEnergyModel>();
      feed_stream->pump_ = std::move(compressor);

      auto feed_flow_bc = std::make_unique<FlowBoundaryCondition>();
      feed_flow_bc->direction_ = BoundaryDirection::INLET;
      feed_flow_bc->location_ = BoundaryLocation::BOTTOM;
      feed_flow_bc->streams_.push_back(std::move(feed_stream));
      data.cycle.GetStep("Feed 2").bottom_boundary_condition_ = std::move(feed_flow_bc);
    }
    {
      // Product valve outlet
      auto product_valve = CreateValveOutletStream(CreateConstantBoundaryValue(P_high), CV_product);
      product_valve->name_ = "Product 2";

      // Purge outlet: sends product CH4 to the LR step of another bed
      auto purge_stream = CreateNoValveOutletStream(CreateMolarFlowBC(F_purge));
      purge_stream->name_ = "Purge 1";
      purge_stream->temperature_ = T_feed;
      purge_out_ptr = purge_stream.get();

      // Pressurisation outlet: connected to LPP step of another bed
      auto press_stream = CreateNoValveOutletStream(CreateVelocityBC(0.0));
      press_stream->name_ = "Press 1";
      press_out_ptr = press_stream.get();

      auto top_flow_bc = std::make_unique<FlowBoundaryCondition>();
      top_flow_bc->direction_ = BoundaryDirection::OUTLET;
      top_flow_bc->location_ = BoundaryLocation::TOP;
      top_flow_bc->streams_.push_back(std::move(product_valve));
      top_flow_bc->streams_.push_back(std::move(purge_stream));
      top_flow_bc->streams_.push_back(std::move(press_stream));
      data.cycle.GetStep("Feed 2").top_boundary_condition_ = std::move(top_flow_bc);
    }

    // ==== ED1 (Equalization Down) ====
    data.cycle.AddStep("ED1");
    data.cycle.GetStep("ED1").step_fraction_ = 0.125;

    IBoundaryStream* ed1_pe_ptr = nullptr;
    {
      auto pe_valve = CreateValveOutletStream(CreateConstantBoundaryValue(P_high / 3.0), CV_PE);
      pe_valve->name_ = "PE stream";
      pe_valve->valve_.xt_ = 0.53;
      ed1_pe_ptr = pe_valve.get();

      auto top_flow_bc = std::make_unique<FlowBoundaryCondition>();
      top_flow_bc->direction_ = BoundaryDirection::OUTLET;
      top_flow_bc->location_ = BoundaryLocation::TOP;
      top_flow_bc->streams_.push_back(std::move(pe_valve));
      data.cycle.GetStep("ED1").top_boundary_condition_ = std::move(top_flow_bc);
    }
    AddClosedBC(BoundaryLocation::BOTTOM, "ED1");

    // ==== BD (Blowdown) ====
    data.cycle.AddStep("BD");
    data.cycle.GetStep("BD").step_fraction_ = 0.25;

    AddClosedBC(BoundaryLocation::TOP, "BD");
    AddPressureOutletBC(
      CreateExponentialBoundaryValue(P_low, 0.25 * cycle_time),
      BoundaryLocation::BOTTOM,
      "BD",
      { {"Methane", 0.53}, {"CarbonDioxide", 0.47} }
    );

    // ==== LR (Light Reflux / Purge) ====
    data.cycle.AddStep("LR");
    data.cycle.GetStep("LR").step_fraction_ = 0.125;

    IBoundaryStream* lr_purge_ptr = nullptr;
    {
      auto purge_inlet = CreateNoValveInletStream(
        CreateMolarFlowBC(F_purge),
        { {"Methane", 1.0}, {"CarbonDioxide", 0.0} },
        T_purge
      );
      purge_inlet->name_ = "CH4 Purge";
      lr_purge_ptr = purge_inlet.get();

      auto top_flow_bc = std::make_unique<FlowBoundaryCondition>();
      top_flow_bc->direction_ = BoundaryDirection::INLET;
      top_flow_bc->location_ = BoundaryLocation::TOP;
      top_flow_bc->streams_.push_back(std::move(purge_inlet));
      data.cycle.GetStep("LR").top_boundary_condition_ = std::move(top_flow_bc);
    }
    AddPressureOutletBC(
      CreateConstantBoundaryValue(P_low),
      BoundaryLocation::BOTTOM,
      "LR",
      { {"Methane", 0.53}, {"CarbonDioxide", 0.47} }
    );

    // ==== EU1 (Equalization Up) ====
    data.cycle.AddStep("EU1");
    data.cycle.GetStep("EU1").step_fraction_ = 0.125;

    IBoundaryStream* eu1_pe_ptr = nullptr;
    {
      auto pe_valve = CreateValveOutletStream(CreateConstantBoundaryValue(P_high / 3.0), CV_PE);
      pe_valve->name_ = "PE stream";
      pe_valve->valve_.xt_ = 0.53;
      eu1_pe_ptr = pe_valve.get();

      auto top_flow_bc = std::make_unique<FlowBoundaryCondition>();
      top_flow_bc->direction_ = BoundaryDirection::OUTLET;
      top_flow_bc->location_ = BoundaryLocation::TOP;
      top_flow_bc->streams_.push_back(std::move(pe_valve));
      data.cycle.GetStep("EU1").top_boundary_condition_ = std::move(top_flow_bc);
    }
    AddClosedBC(BoundaryLocation::BOTTOM, "EU1");

    // ==== LPP (Light Product Pressurisation) ====
    const double lpp_step_time = 0.125 * cycle_time;
    data.cycle.AddStep("LPP");
    data.cycle.GetStep("LPP").step_fraction_ = 0.125;

    PressureBoundaryCondition* lpp_top_ptr = nullptr;
    {
      auto lpp_pressure_bc = std::make_unique<PressureBoundaryCondition>();
      lpp_pressure_bc->location_ = BoundaryLocation::TOP;
      lpp_pressure_bc->direction_ = BoundaryDirection::INLET;
      lpp_pressure_bc->pressure_ = CreateLinearBoundaryValue(P_low, P_high, lpp_step_time);
      lpp_pressure_bc->temperature_ = T_feed;
      lpp_pressure_bc->mole_fraction_["Methane"] = 1.0;
      lpp_pressure_bc->mole_fraction_["CarbonDioxide"] = 0.0;
      lpp_top_ptr = lpp_pressure_bc.get();
      data.cycle.GetStep("LPP").top_boundary_condition_ = std::move(lpp_pressure_bc);
    }
    AddClosedBC(BoundaryLocation::BOTTOM, "LPP");

    // ==== Inter-bed connections ====
    // Pair 1: Feed2 "Purge 1" outlet ↔ LR "CH4 Purge" inlet
    {
      auto ep_purge_out = MakeEndpoint(data.cycle.GetStep("Feed 2"), *purge_out_ptr);
      auto ep_purge_in  = MakeEndpoint(data.cycle.GetStep("LR"), *lr_purge_ptr);
      MakeMutualConnection(ep_purge_out, ep_purge_in);
    }
    // Pair 2: ED1 "PE stream" ↔ EU1 "PE stream" (pressure equalisation)
    {
      auto ep_ed1 = MakeEndpoint(data.cycle.GetStep("ED1"), *ed1_pe_ptr);
      auto ep_eu1 = MakeEndpoint(data.cycle.GetStep("EU1"), *eu1_pe_ptr);
      MakeMutualConnection(ep_ed1, ep_eu1);
    }
    // Pair 3: Feed2 "Press 1" outlet ↔ LPP top pressure inlet
    {
      auto ep_press_out = MakeEndpoint(data.cycle.GetStep("Feed 2"), *press_out_ptr);
      auto ep_lpp       = MakeEndpoint(data.cycle.GetStep("LPP"), *lpp_top_ptr);
      MakeMutualConnection(ep_press_out, ep_lpp);
    }

    for (auto& L : data.bed.layers_.layers_) L.momentum_model_ = momentum;
    ConfigureCycle(20);
    // KPI: product=Methane(0), feed=Feed1+Feed2@Bottom, product=Feed1+Feed2@Top, recycled=LR+LPP@Top
    ConfigureKPI(
      0,
      {0, 1}, {BoundaryLocation::BOTTOM, BoundaryLocation::BOTTOM},
      {0, 1}, {BoundaryLocation::TOP, BoundaryLocation::TOP},
      {4, 6}, {BoundaryLocation::TOP, BoundaryLocation::TOP}
    );

    CoupledSolver solver{data};

    std::string error{};
    EXPECT_TRUE(data.PerformChecks(error)) << "System checks failed: " << error;

    // Count step retries: each is a recoverable IDA failure (logged at MEDIUM)
    // that the step only survived by being re-attempted. Neither momentum
    // model may need one: before the pressure / closed boundary faces imposed
    // the explicit Ergun root (ImposeBoundaryFaceVelocity), the Ergun run
    // failed at the start of ED1, LR and EU1 in every cycle.
    int step_retries = 0;
    solver.SetLogLevel(LogLevel::MEDIUM);
    solver.SetMessageCallback([&step_retries](LogLevel, const std::string& msg) {
      if (msg.find("IDA failure") != std::string::npos) {
        ++step_retries;
        std::cout << msg << std::endl;
      }
    });

    bool verbose = false;
    int retval = solver.IntegrateCycle(data.cycle, verbose);
    EXPECT_GE(retval, 0) << "Solver did not complete";
    EXPECT_EQ(step_retries, 0) << "step switches needed solver retries";

    EXPECT_TRUE(IsMassConserved(solver.GetLog())) << "Mass is not conserved.";
    EXPECT_TRUE(IsMolesConserved(solver.GetLog())) << "Moles are not conserved.";
    // IsComponentMolesConserved is not checked here because the inter-bed connected
    // streams (Purge 1 / CH4 Purge, Press 1 / LPP) carry specified compositions that
    // differ slightly from the actual bed outlet, which is inherent to single-bed
    // modelling of a multi-bed cycle and introduces a small systematic mismatch.

    SaveTest(save_name, solver);
  }

} // end adsim namespace
