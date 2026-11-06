// test_heating_medium_solver.cpp
//
// A step with a hot heating medium on the wall's outer face heats a closed,
// non-adsorbing bed towards the medium temperature, and steps without a
// medium stay on the ambient exchange. Three steps:
//   Hold  (no medium, ambient = T0): nothing moves;
//   Heat  (medium at T_m > T0):      the wall and the bed rise towards T_m
//                                    and never pass it;
//   Rest  (no medium, ambient = T0): the bed cools back towards the ambient.
// Short and coarse (fast tier); the analytic heat-up check is separate.

#include <algorithm>
#include <cmath>
#include <string>

#include <gtest/gtest.h>

#include "../src/pch.h"
#include "adsim/boundary_conditions/step.h"
#include "adsim/solver/coupled_solver.h"

#include "../fixtures/user_data_test.h"

namespace adsim {

  namespace {

    double InteriorMean(const std::vector<double>& v) {
      double s = 0.0;
      for (size_t i = 1; i + 1 < v.size(); ++i) s += v[i];
      return s / static_cast<double>(v.size() - 2);
    }

  } // namespace

  class HeatingMediumSolverTest : public UserDataTest {};

  TEST_F(HeatingMediumSolverTest, HotMediumHeatsAClosedBedTowardsItsTemperature) {
    data.info.title = "Heating medium: closed bed heat-up";
    AddComponents({"Nitrogen", "Oxygen"});
    auto [layer_id, ads_id] = AddLayer("Packing");
    auto& layer = GetLayer(layer_id);
    layer.domain_.length = 0.2;
    layer.domain_.n_cells = 20;
    layer.eb_ = 0.4;
    layer.particle_diameter_ = 0.002;
    std::get<CylinderBed>(data.bed.shape_).diameter = 0.02;
    auto& ads = GetAdsorbent(ads_id);
    ads.solid_density_ = 1000.0;
    ads.heat_capacity_ = 900.0;
    ads.ep_ = 0.0;

    data.bed.wall_.thickness_ = 0.002;
    data.bed.wall_.density_ = 7900.0;
    data.bed.wall_.heat_capacity_ = 500.0;
    data.bed.wall_.thermal_conductivity_ = 16.2;
    data.bed.equations_.wall_temperature_boundary_condition_ =
        CreateWallTemperatureBoundaryCondition(WallTemperatureBoundaryConditionType::ZERO_FLUX);

    AddConstantInternalHeatTransferCoefficientModel(100.0, layer_id);
    AddConstantExternalHeatTransferCoefficientModel(10.0, layer_id);

    const double T0 = 300.0, T_medium = 380.0;
    data.bed.ambient_temperature_ = T0;
    AddInitialConditions(1e5, T0, T0, {{"Nitrogen", 1.0}, {"Oxygen", 0.0}});

    for (const char* name : {"Hold", "Heat", "Rest"}) {
      AddStep(name, name == std::string("Heat") ? 300.0 : 100.0);
      AddClosedBC(BoundaryLocation::BOTTOM, name);
      AddClosedBC(BoundaryLocation::TOP, name);
    }
    auto& medium = data.cycle.GetStep("Heat").heating_medium_;
    medium.enabled_ = true;
    medium.temperature_ = T_medium;
    medium.heat_transfer_coefficient_ = 500.0;
    ConfigureCycle(1);

    CoupledSolver solver{data};
    std::string error;
    ASSERT_TRUE(data.PerformChecks(error)) << error;
    ASSERT_GE(solver.IntegrateCycle(data.cycle, false), 0) << "solver did not complete";

    const auto& log = solver.GetLog();
    const auto hold_T = log.GetVariable(0, 0, VariableType::TEMPERATURE);
    const auto heat_T = log.GetVariable(0, 1, VariableType::TEMPERATURE);
    const auto heat_Tw = log.GetVariable(0, 1, VariableType::WALL_TEMPERATURE);
    const auto heat_ho = log.GetVariable(0, 1, VariableType::EXTERNAL_HEAT_TRANSFER_COEFFICIENT);
    const auto rest_T = log.GetVariable(0, 2, VariableType::TEMPERATURE);
    const auto rest_ho = log.GetVariable(0, 2, VariableType::EXTERNAL_HEAT_TRANSFER_COEFFICIENT);
    ASSERT_FALSE(hold_T.values.empty());
    ASSERT_FALSE(heat_T.values.empty());
    ASSERT_FALSE(rest_T.values.empty());

    // Hold: bed, wall and surroundings all at T0, so nothing moves.
    EXPECT_NEAR(InteriorMean(hold_T.values.back()), T0, 1e-6);

    // Heat: the logged outer coefficient is the medium's.
    EXPECT_DOUBLE_EQ(heat_ho.values.back()[5], 500.0);
    double previous = T0;
    for (const auto& profile : heat_T.values) {
      const double T = InteriorMean(profile);
      EXPECT_GE(T, previous - 1e-6);       // heats monotonically
      EXPECT_LE(T, T_medium + 1e-6);       // never passes the medium
      previous = T;
    }
    const double T_heated = InteriorMean(heat_T.values.back());
    const double Tw_heated = InteriorMean(heat_Tw.values.back());
    EXPECT_GT(T_heated, T0 + 0.8 * (T_medium - T0)) << "bed did not approach the medium temperature";
    EXPECT_GT(Tw_heated, T_heated);  // heat flows medium -> wall -> bed
    EXPECT_LE(Tw_heated, T_medium + 1e-6);

    // Rest: back on the layer's coefficient and the ambient, so it cools.
    EXPECT_DOUBLE_EQ(rest_ho.values.back()[5], 10.0);
    EXPECT_LT(InteriorMean(rest_T.values.back()), T_heated - 1e-3);
  }

} // namespace adsim
