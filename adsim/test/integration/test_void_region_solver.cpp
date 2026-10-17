// Integration tests for the top/bottom void space feature.
//
// Verifies that the solver runs end-to-end with voids enabled and that
// physically expected effects appear in the output.

#define _USE_MATH_DEFINES
#include <cmath>

#include "../src/pch.h"
#include "adsim/equations/dispersion_model.h"
#include "adsim/equations/heat_transfer_model.h"
#include "adsim/equations/isotherm.h"
#include "adsim/solver/coupled_solver.h"

#include <gtest/gtest.h>

#include "../fixtures/user_data_test.h"


namespace adsim {

  // Build a simple inert-tracer feed problem on a small bed.
  class VoidRegionSolverTest : public UserDataTest {
  protected:
    void BuildInertTracerFeed() {
      data.info.title = "Void region inert tracer";

      AddComponents({"CarbonDioxide", "Nitrogen"});

      auto [layer_id, ads_id] = AddLayer("Inert Solid");
      auto& layer = GetLayer(layer_id);
      layer.domain_.length = 0.3;
      layer.domain_.n_cells = 20;
      layer.particle_diameter_ = 1e-3;
      layer.eb_ = 0.4;

      auto& ads = GetAdsorbent(ads_id);
      ads.solid_density_ = 1500.0;
      ads.heat_capacity_ = 900.0;
      ads.thermal_conductivity_ = 0.5;
      ads.ep_ = 0.35;

      // Henry isotherm with k=0 -> effectively inert.
      auto henry = std::make_unique<HenryIsotherm>();
      henry->params_["CarbonDioxide"].k = 0.0;
      henry->params_["CarbonDioxide"].Q = 0.0;
      henry->params_["Nitrogen"].k = 0.0;
      henry->params_["Nitrogen"].Q = 0.0;
      AddIsotherm({"CarbonDioxide", "Nitrogen"}, ads_id, std::move(henry));

      AddLinearDrivingForceModel(ads_id, "CarbonDioxide", 1.0);
      AddLinearDrivingForceModel(ads_id, "Nitrogen", 1.0);
      SetHOA(ads_id, "CarbonDioxide", 0.0);
      SetHOA(ads_id, "Nitrogen", 0.0);

      AddConstantInternalHeatTransferCoefficientModel(10.0, layer_id);
      AddConstantExternalHeatTransferCoefficientModel(10.0, layer_id);

      layer.equations_.dispersion_model_ = CreateAxialDispersionModel(AxialDispersionType::RUTHVEN);
      layer.equations_.heat_dispersion_model_ = std::make_unique<HeatDispersionModelConstant>();

      std::get<CylinderBed>(data.bed.shape_).diameter = 0.05;

      AddInitialConditions(
        1.0e5, 298.0, 298.0,
        {{"CarbonDioxide", 0.0}, {"Nitrogen", 1.0}});

      AddStep("Feed", 60.0);

      // Bottom: flow inlet with CO2 tracer
      Step& step_feed = data.cycle.GetStep("Feed");
      auto inlet = CreateNoValveInletStream(
        CreateVolumeFlowBC(5e-5),
        {{"CarbonDioxide", 1.0}, {"Nitrogen", 0.0}},
        298.0);

      auto bottom = std::make_unique<FlowBoundaryCondition>();
      bottom->direction_ = BoundaryDirection::INLET;
      bottom->location_ = BoundaryLocation::BOTTOM;
      bottom->streams_.push_back(std::move(inlet));
      step_feed.bottom_boundary_condition_ = std::move(bottom);

      auto top = std::make_unique<PressureBoundaryCondition>();
      top->direction_ = BoundaryDirection::OUTLET;
      top->location_ = BoundaryLocation::TOP;
      top->pressure_ = CreateConstantBoundaryValue(1e5);
      top->mole_fraction_["CarbonDioxide"] = 0.0;
      top->mole_fraction_["Nitrogen"] = 1.0;
      step_feed.top_boundary_condition_ = std::move(top);

      ConfigureCycle(1);
    }
  };

  TEST_F(VoidRegionSolverTest, SolverCompletesWithVoidsDisabled) {
    BuildInertTracerFeed();
    std::string err;
    ASSERT_TRUE(data.PerformChecks(err)) << err;

    CoupledSolver solver{data};
    int retval = solver.IntegrateCycle(data.cycle, false);
    EXPECT_GE(retval, 0) << "solver did not complete";
    EXPECT_TRUE(IsMassConserved(solver.GetLog()));
  }

  TEST_F(VoidRegionSolverTest, SolverCompletesWithBothVoids) {
    BuildInertTracerFeed();

    // Enable voids comparable to ~10% of the bed volume each.
    const double area = data.bed.GetArea();
    const double layer_volume = data.bed.layers_.layers_[0].domain_.length.Get() * area;
    data.bed.bottom_void_.enabled_ = true;
    data.bed.bottom_void_.volume_ = 0.1 * layer_volume;
    data.bed.top_void_.enabled_ = true;
    data.bed.top_void_.volume_ = 0.1 * layer_volume;

    std::string err;
    ASSERT_TRUE(data.PerformChecks(err)) << err;

    CoupledSolver solver{data};
    int retval = solver.IntegrateCycle(data.cycle, false);
    EXPECT_GE(retval, 0) << "solver did not complete with voids enabled";
    EXPECT_TRUE(IsMassConserved(solver.GetLog()))
      << "mole fractions should sum to 1 even with voids";
  }

  TEST_F(VoidRegionSolverTest, VoidCellSizeIsCorrect) {
    // Verify the lumped void physics: with the void enabled, the cell width
    // equals user volume / bed area, and the cell volume equals user volume.
    BuildInertTracerFeed();
    const double area = data.bed.GetArea();
    const double V_user = 5e-5;

    data.bed.bottom_void_.enabled_ = true;
    data.bed.bottom_void_.volume_ = V_user;

    std::string err;
    ASSERT_TRUE(data.PerformChecks(err)) << err;
    CoupledSolver solver{data};
    int retval = solver.IntegrateCycle(data.cycle, false);
    EXPECT_GE(retval, 0) << "solver should complete with bottom void";

    const auto& dv = data.state.domain_vectors;
    const int iv = data.bed.GetBottomVoidCellIndex();
    EXPECT_NEAR(dv.volume[iv], V_user, 1e-15);
    EXPECT_NEAR(dv.cell_width[iv], V_user / area, 1e-15);
  }

  TEST_F(VoidRegionSolverTest, UserResistanceRaisesInletPressure) {
    // A non-zero VoidRegion::resistance_coefficient_ adds an in-series
    // resistor at the void↔bed face. For a fixed feed flow rate Q, this
    // forces the inlet (void-side) pressure higher than the no-resistance
    // baseline by approximately u * R.
    auto inlet_pressure = [&](double R) {
      data = UserData{};
      SetUp();
      BuildInertTracerFeed();

      const double area = data.bed.GetArea();
      const double layer_volume = data.bed.layers_.layers_[0].domain_.length.Get() * area;
      data.bed.bottom_void_.enabled_ = true;
      data.bed.bottom_void_.volume_ = 0.1 * layer_volume;
      data.bed.bottom_void_.resistance_coefficient_ = R;

      std::string err;
      EXPECT_TRUE(data.PerformChecks(err)) << err;
      CoupledSolver solver{data};
      int retval = solver.IntegrateCycle(data.cycle, false);
      EXPECT_GE(retval, 0);

      const auto& dv = data.state.domain_vectors;
      const int iv = data.bed.GetBottomVoidCellIndex();
      return data.state.pressure.var[iv];
    };

    const double P_no_R   = inlet_pressure(0.0);
    const double P_with_R = inlet_pressure(1e6);  // sizeable resistance
    EXPECT_GT(P_with_R, P_no_R)
      << "user resistance should raise the void-side pressure "
      << "(no-R: " << P_no_R << " Pa, with-R: " << P_with_R << " Pa)";
  }

  TEST_F(VoidRegionSolverTest, VoidEnabledShiftsBreakthroughLater) {
    // With a sizable bottom void, the manifold smooths the step feed and
    // adds gas-phase residence time, so the outlet CO2 0.5-crossing time
    // must be later than the no-void baseline.
    auto run = [&](double bottom_void_volume) {
      data = UserData{};
      SetUp();
      BuildInertTracerFeed();

      if (bottom_void_volume > 0.0) {
        data.bed.bottom_void_.enabled_ = true;
        data.bed.bottom_void_.volume_ = bottom_void_volume;
      }

      std::string err;
      EXPECT_TRUE(data.PerformChecks(err)) << err;
      CoupledSolver solver{data};
      int retval = solver.IntegrateCycle(data.cycle, false);
      EXPECT_GE(retval, 0);

      auto [t, y_co2] = GetOutletMoleFraction(solver, "CarbonDioxide", 0, 0, BoundaryLocation::TOP);
      double t_half = -1.0;
      for (size_t i = 0; i < t.size(); ++i) {
        if (y_co2[i] >= 0.5) { t_half = t[i]; break; }
      }
      return t_half;
    };

    const double t_no_void = run(0.0);
    ASSERT_GT(t_no_void, 0.0) << "baseline breakthrough not detected";

    // 30% of bed volume.
    const double layer_volume =
        data.bed.layers_.layers_[0].domain_.length.Get() * data.bed.GetArea();
    const double V_void = 0.3 * layer_volume;
    const double Q = 5e-5;  // matches BuildInertTracerFeed
    const double tau = V_void / Q;

    const double t_with_void = run(V_void);
    ASSERT_GT(t_with_void, 0.0) << "void breakthrough not detected";
    EXPECT_GT(t_with_void, t_no_void + 0.3 * tau)
      << "void should delay outlet breakthrough by ~tau = V_void/Q "
      << "(no-void: " << t_no_void << "s, with-void: " << t_with_void
      << "s, tau: " << tau << "s)";
  }

} // end adsim namespace
