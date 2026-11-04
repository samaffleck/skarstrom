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

  // Moles of a component the bed accumulated in a step (gas + solid) minus
  // the net moles that crossed its two ends (GetTotalComponentMolesIn, the
  // stream integrals the KPIs and the cycle balance use): zero for a closed
  // balance. And the moles of it that crossed the two ends, either way.
  static double ComponentBalanceError(const UserDataTest& t, const EquationLog& log, const std::string& comp) {
    return t.CalculateTotalAccumulated(log, 0, 0, comp)
         - (GetTotalComponentMolesIn(log, 0, 0, BoundaryLocation::BOTTOM, comp)
            + GetTotalComponentMolesIn(log, 0, 0, BoundaryLocation::TOP, comp));
  }
  static double ComponentThroughput(const EquationLog& log, const std::string& comp) {
    return std::abs(GetTotalComponentMolesIn(log, 0, 0, BoundaryLocation::BOTTOM, comp))
         + std::abs(GetTotalComponentMolesIn(log, 0, 0, BoundaryLocation::TOP, comp));
  }

  // The direct-air-capture charge, run with the given advection scheme.
  // Besides the generic conservation checks it bounds the logged balances at
  // the level upwinding closes them. While the high-order schemes clamped
  // their bed-end faces to [0, 1] (WP20-1), the solver sent the clamped face
  // value out of the bed while the stream logs integrate the end cell's
  // value, and the pure-N2 outlet cell's mole-fraction sum drifted (4.0e-5
  // with WENO): the total balance was off by 3.45e-2 mol with WENO (1.1e-5
  // of the 3.1e3 mol throughput) and the CO2 balance by 1.2e-6 mol with WENO
  // and van Leer (1.9e-6 of the CO2 fed).
  static void RunDAC(UserDataTest& t, AdvectionScheme scheme, const std::string& save_name) {
    auto& data = t.data;
    data.info.title = "DAC example";
    data.info.description = "Direct air capture CO2 battery.";
    data.info.author = "Sam Affleck";
    
    t.AddComponents({ "CarbonDioxide", "Nitrogen" });

    const double temp = 293;
    data.bed.ambient_temperature_ = temp;
    const double internal_diameter = 0.16; //[m]
    const double external_diameter = 0.18;
    data.bed.wall_.thickness_ = (external_diameter - internal_diameter) / 2.0;

    double bulk_density = 450;
    double eb = 0.4;
    double solid_density = bulk_density / (1 - eb);

    std::string adsorbent_name = "Lewatit VP OC 1065";
    auto [layer_id, ads_id] = t.AddLayer(adsorbent_name);
    auto& layer = t.GetLayer(layer_id);
    layer.domain_.length = 1.65;
    std::get<CylinderBed>(data.bed.shape_).diameter = internal_diameter;
    layer.particle_diameter_ = 0.002;
    layer.eb_ = eb;
    layer.domain_.n_cells = 30;

    auto& ads_13X = t.GetAdsorbent(ads_id);
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
    t.AddIsotherm({ "CarbonDioxide" }, ads_id, std::move(toth));

    const double pressure = 1.01325e5;

    t.AddLinearDrivingForceModel(ads_id, "CarbonDioxide", 1.0);
     
    t.AddInitialConditions(pressure, temp, temp,
    { 
      {"CarbonDioxide", 0.0 },
      {"Nitrogen", 1.0 }
    });

    double feed_time = 3600; // [1h]
    t.AddStep("Charge", feed_time);

    Step& step_feed = data.cycle.GetStep("Charge");
    auto feed_inlet_stream = t.CreateNoValveInletStream(
      t.CreateVelocityBC(0.5),
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
    feed_pressure_out->pressure_ = t.CreateConstantBoundaryValue(pressure);
    feed_pressure_out->mole_fraction_["CarbonDioxide"] = 0.0;
    feed_pressure_out->mole_fraction_["Nitrogen"] = 1.0;
    step_feed.top_boundary_condition_ = std::move(feed_pressure_out);

    t.ConfigureCycle(1);
    data.settings.advection_scheme_ = scheme;
    CoupledSolver solver{data};

    std::string error{};
    EXPECT_TRUE(data.PerformChecks(error)) << "System checks failed: " << error;
    
    bool verbose = false;
    int retval = solver.IntegrateCycle(data.cycle, verbose);
    EXPECT_GE(retval, 0) << "Solver did not complete";

    if (!save_name.empty()) t.SaveTest(save_name, solver);

    EXPECT_TRUE(t.IsMassConserved(solver.GetLog())) << "Total mole fraction is less than 1.";
    EXPECT_TRUE(t.IsMolesConserved(solver.GetLog())) << "Mass is not conserved.";
    EXPECT_TRUE(t.IsComponentMolesConserved(solver.GetLog())) << "Component moles not conserved.";

    // Logged balances at the upwind level. Upwinding closes the total to
    // 2.7e-5 mol of a 3.1e3 mol throughput (8.5e-9) and CO2 to 2.2e-10 mol of
    // the 0.62 mol fed (3.5e-10); the high-order schemes to 2.7e-5 mol and
    // 2.7e-9 mol. Over 15 integration paths (reltol 5e-7 to 2e-6, initial
    // step 1e-10 to 1e-9) every scheme stays at or below 3.0e-5 mol (total)
    // and 1.1e-8 mol (CO2). The bounds, 1e-7 of each throughput (a tenth of
    // reltol), are a factor of 10 (total) and 6 (CO2) above those, and a
    // factor of 110 (total, WENO) and 18 (CO2) below the errors of the
    // clamped faces.
    const auto& log = solver.GetLog();
    const double co2_error = ComponentBalanceError(t, log, "CarbonDioxide");
    const double total_error = co2_error + ComponentBalanceError(t, log, "Nitrogen");
    const double co2_through = ComponentThroughput(log, "CarbonDioxide");
    const double total_through = co2_through + ComponentThroughput(log, "Nitrogen");
    std::cout << "[DAC] scheme " << static_cast<int>(scheme) << ": CO2 balance error " << co2_error
              << " mol of " << co2_through << ", total " << total_error << " mol of " << total_through << std::endl;
    EXPECT_LT(std::abs(co2_error), 1e-7 * co2_through) << "CO2 balance does not close";
    EXPECT_LT(std::abs(total_error), 1e-7 * total_through) << "total mole balance does not close";
  }

  TEST_F(UserDataTest, DAC) {
    RunDAC(*this, AdvectionScheme::Upwind, "DAC");
  }

  TEST_F(UserDataTest, DACWENO) {
    RunDAC(*this, AdvectionScheme::WENO, "");
  }

  TEST_F(UserDataTest, DACVanLeer) {
    RunDAC(*this, AdvectionScheme::VanLeer, "");
  }

} // end adsim namespace
