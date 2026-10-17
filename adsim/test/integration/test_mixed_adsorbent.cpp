// Integration test for multi-adsorbent (mixed) layers.
//
// Core equivalence property: a layer packed with a SINGLE adsorbent at 100%
// must produce identical results to the same layer packed with TWO adsorbents
// that are physically identical, each at 50% mass fraction. The per-adsorbent
// source terms sum as sum_k rho_b,k * dq_k/dt = (0.5 + 0.5) * rho_b * dq/dt,
// so the two configurations are mathematically the same simulation.

#define _USE_MATH_DEFINES
#include <cmath>

#include "../src/pch.h"
#include "adsim/solver/coupled_solver.h"
#include "adsim/math/constants.h"
#include "adsim/equations/isotherm.h"
#include "adsim/equations/dispersion_model.h"
#include "adsim/equations/heat_transfer_model.h"
#include "adsim/boundary_conditions/boundary_conditions.h"
#include "adsim/boundary_conditions/boundary_stream.h"
#include "adsim/boundary_conditions/step.h"

#include <gtest/gtest.h>

#include "../fixtures/user_data_test.h"


namespace adsim {

  class MixedAdsorbentTest : public UserDataTest {
  public:
    const std::string CO2 = "CarbonDioxide";
    const std::string N2  = "Nitrogen";

    // Adds a CO2 DualSiteLangmuir (concentration form) isotherm + constant-MTC
    // LDF + heat of adsorption to the given adsorbent. Identical parameters for
    // every adsorbent so that a split packing is physically equivalent.
    void ConfigureAdsorbent(AdsorbentID ads_id) {
      auto& ads = GetAdsorbent(ads_id);
      ads.solid_density_        = 1100.0;
      ads.heat_capacity_        = 900.0;
      ads.thermal_conductivity_ = 0.5;
      ads.ep_                   = 0.35;

      auto dsl = std::make_unique<DualSiteLangmuirIsotherm_Concentration>(CO2);
      dsl->params_[CO2].m1 = 2.5;  dsl->params_[CO2].b1 = 4.0e-7; dsl->params_[CO2].q1 = 30000;
      dsl->params_[CO2].m2 = 3.0;  dsl->params_[CO2].b2 = 1.0e-7; dsl->params_[CO2].q2 = 25000;
      AddIsotherm({CO2}, ads_id, std::move(dsl));

      AddLinearDrivingForceModel(ads_id, CO2, 0.2);
      AddLinearDrivingForceModel(ads_id, N2, 0.0);
      SetHOA(ads_id, CO2, 30000.0);
      SetHOA(ads_id, N2, 0.0);
    }

    // Builds a one-layer breakthrough column. If `split` is true the layer holds
    // two identical adsorbents at 50% each; otherwise a single adsorbent at 100%.
    void BuildColumn(bool split) {
      AddComponents({CO2, N2});

      auto [layer_id, ads_id] = AddLayer("Sorbent A");
      auto& layer = GetLayer(layer_id);
      layer.domain_.length  = 0.3;
      layer.particle_diameter_ = 0.0015;
      layer.eb_             = 0.4;
      layer.domain_.n_cells = 20;
      std::get<CylinderBed>(data.bed.shape_).diameter = 0.02;
      ConfigureAdsorbent(ads_id);

      if (split) {
        // Replace the 100% composition with two identical adsorbents at 50%.
        auto ads_id2 = data.adsorbent_directory.Add("Sorbent A copy");
        ConfigureAdsorbent(ads_id2);
        layer.composition_.clear();
        layer.composition_[ads_id]  = 0.5;
        layer.composition_[ads_id2] = 0.5;
      }

      layer.equations_.dispersion_model_ = CreateAxialDispersionModel(AxialDispersionType::RUTHVEN);
      layer.equations_.heat_dispersion_model_ = CreateAxialHeatDispersionModel(AxialHeatDispersionType::CONSTANT);
      AddConstantInternalHeatTransferCoefficientModel(8.0, layer_id);
      AddConstantExternalHeatTransferCoefficientModel(10.0, layer_id);
      data.bed.ambient_temperature_ = 293.0;

      AddInitialConditions(1.0e5, 293.0, 293.0, {{CO2, 0.0}, {N2, 1.0}});

      AddStep("Feed", 200.0);
      auto inlet = CreateNoValveInletStream(CreateVelocityBC(0.05), {{CO2, 0.2}, {N2, 0.8}}, 298.0);
      inlet->name_ = "Feed inlet";
      auto in_bc = std::make_unique<FlowBoundaryCondition>();
      in_bc->direction_ = BoundaryDirection::INLET;
      in_bc->location_  = BoundaryLocation::BOTTOM;
      in_bc->streams_.push_back(std::move(inlet));
      data.cycle.GetStep("Feed").bottom_boundary_condition_ = std::move(in_bc);

      auto out_bc = std::make_unique<PressureBoundaryCondition>();
      out_bc->pressure_  = CreateConstantBoundaryValue(1.0e5);
      out_bc->direction_ = BoundaryDirection::OUTLET;
      out_bc->location_  = BoundaryLocation::TOP;
      out_bc->mole_fraction_[CO2] = 0.2;
      out_bc->mole_fraction_[N2]  = 0.8;
      data.cycle.GetStep("Feed").top_boundary_condition_ = std::move(out_bc);

      ConfigureCycle(1);
    }
  };

  TEST_F(MixedAdsorbentTest, SplitAdsorbentEqualsSingle) {
    // Run the single-adsorbent column.
    BuildColumn(/*split=*/false);
    CoupledSolver solver_single{data};
    ASSERT_GE(solver_single.IntegrateCycle(data.cycle, false), 0);
    EXPECT_EQ(data.state.GetAdsorbentIDs().size(), 1u);  // one adsorbent
    auto [t_single, co2_single] = GetOutletMoleFraction(solver_single, CO2, 0, 0, BoundaryLocation::TOP);

    // Reset and run the equivalent split (50/50) column.
    data = UserData{};
    UserDataTest::SetUp();
    BuildColumn(/*split=*/true);
    CoupledSolver solver_split{data};
    ASSERT_GE(solver_split.IntegrateCycle(data.cycle, false), 0);
    EXPECT_EQ(data.state.GetAdsorbentIDs().size(), 2u);  // two adsorbents in one layer
    auto [t_split, co2_split] = GetOutletMoleFraction(solver_split, CO2, 0, 0, BoundaryLocation::TOP);

    // The two outlet breakthrough curves must coincide. Compare the split result
    // against the single result interpolated onto the split time grid.
    ASSERT_FALSE(co2_single.empty());
    ASSERT_FALSE(co2_split.empty());
    ASSERT_TRUE(is_strictly_increasing(t_single));

    double max_abs_diff = 0.0;
    for (size_t i = 0; i < t_split.size(); ++i) {
      double ref = interp_linear(t_single, co2_single, t_split[i]);
      max_abs_diff = std::max(max_abs_diff, std::abs(co2_split[i] - ref));
    }
    EXPECT_LT(max_abs_diff, 1e-6)
      << "A 50/50 split of an identical adsorbent must reproduce the single-adsorbent breakthrough curve";
  }

} // end adsim namespace
