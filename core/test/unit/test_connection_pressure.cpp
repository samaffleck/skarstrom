// test_connection_pressure.cpp
//
// A connected set-flow stream (BoundaryStream, no device) replays its
// partner's logged flow. Gas passes through the connection only from the
// higher to the lower pressure: the replayed flow is scaled by
// ConnectionFlowWeight of the pressure difference at the same instant
// (closure smoothed over kConnectionPressureBand).

#include <algorithm>
#include <cmath>
#include <string>
#include <vector>

#include "../src/pch.h"
#include "adsim/objects/deck_builder.h"
#include "adsim/solver/coupled_solver.h"
#include "adsim/boundary_conditions/connection.h"
#include "adsim/boundary_conditions/boundary_conditions.h"
#include "adsim/boundary_conditions/boundary_stream.h"
#include "adsim/boundary_conditions/step.h"

#include <gtest/gtest.h>

#include "../fixtures/user_data_test.h"

namespace adsim {

  // ---- the smoothed closure ------------------------------------------------

  // Shut at or below zero, exactly open beyond the band, and C1 at both ends
  // of the band (a hard switch has a kink that stalls Newton as the flow
  // decays through zero).
  TEST(ConnectionFlowWeight, IsSmoothAndExactOutsideTheBand) {
    const double band = kConnectionPressureBand;
    EXPECT_EQ(ConnectionFlowWeight(-1e5), 0.0);
    EXPECT_EQ(ConnectionFlowWeight(0.0), 0.0);
    EXPECT_EQ(ConnectionFlowWeight(band), 1.0);
    EXPECT_EQ(ConnectionFlowWeight(2e5), 1.0);
    EXPECT_TRUE(std::isnan(ConnectionFlowWeight(std::nan(""))));

    const double h = 1e-6 * band;
    auto slope = [&](double x) { return (ConnectionFlowWeight(x + h) - ConnectionFlowWeight(x - h)) / (2.0 * h); };
    EXPECT_NEAR(slope(0.0), 0.0, 1e-6 / band);
    EXPECT_NEAR(slope(band), 0.0, 1e-6 / band);
    double previous = 0.0;
    for (int k = 1; k <= 100; ++k) {
      const double w = ConnectionFlowWeight(band * k / 100.0);
      EXPECT_GE(w, previous);
      EXPECT_LE(w, 1.0);
      previous = w;
    }
  }

  // ---- a light-product pressurisation pair ---------------------------------

  class ConnectionPressureTest : public UserDataTest {
  protected:
    static constexpr double kAtm = 101325.0;
    const std::unordered_map<std::string, double> feed{ {"CarbonDioxide", 0.15}, {"Nitrogen", 0.85} };

    void SetUp() override {
      UserDataTest::SetUp();
      AddComponents({ "CarbonDioxide", "Nitrogen" });
      AddInitialConditions(kAtm, 293, 293, { {"CarbonDioxide", 0.0}, {"Nitrogen", 1.0} });
      auto [layer_id, ads_id] = AddLayer("Adsorbent");
      GetLayer(layer_id).domain_.n_cells = 10;
      AddIsotherm({ "Nitrogen" }, ads_id, std::make_unique<HenryIsotherm>("Nitrogen"));
      AddIsotherm({ "CarbonDioxide" }, ads_id, std::make_unique<LangmuirIsotherm>("CarbonDioxide"));
      AddLinearDrivingForceModel(ads_id, "CarbonDioxide", 0.5);
      AddLinearDrivingForceModel(ads_id, "Nitrogen", 0.5);
      AddConstantInternalHeatTransferCoefficientModel(0.0, layer_id);
      AddConstantExternalHeatTransferCoefficientModel(0.0, layer_id);
    }

    FlowBoundaryCondition& Flow(const std::string& step, BoundaryLocation loc) {
      auto& s = data.cycle.GetStep(step);
      auto& bc = loc == BoundaryLocation::TOP ? s.top_boundary_condition_ : s.bottom_boundary_condition_;
      auto* flow = dynamic_cast<FlowBoundaryCondition*>(bc.get());
      EXPECT_NE(flow, nullptr);
      return *flow;
    }
    PressureBoundaryCondition& Pressure(const std::string& step, BoundaryLocation loc) {
      auto& s = data.cycle.GetStep(step);
      auto& bc = loc == BoundaryLocation::TOP ? s.top_boundary_condition_ : s.bottom_boundary_condition_;
      auto* p = dynamic_cast<PressureBoundaryCondition*>(bc.get());
      EXPECT_NE(p, nullptr);
      return *p;
    }

    // Feed (the donor): a small feed at the bottom; at the top a set-flow
    // light-product stream "LPP" connected to the receiver and a vent whose
    // check valve keeps ambient gas out. The draw is larger than the feed,
    // so the donor's pressure falls through the step.
    // Blowdown: the bottom to 0.3 atm.
    // Pressurisation (the receiver): the top ramps from the bed pressure to
    // 1 atm, through the connection; it overtakes the falling donor before
    // the end of the step.
    void LppCycle(int cycles) {
      AddStep("Feed", 20.0);
      {
        std::vector<std::unique_ptr<IBoundaryStream>> streams;
        streams.push_back(CreateNoValveInletStream(CreateVelocityBC(0.005), feed, 293));
        AddFlowBC(std::move(streams), BoundaryLocation::BOTTOM, BoundaryDirection::INLET, "Feed");
      }
      {
        std::vector<std::unique_ptr<IBoundaryStream>> streams;
        streams.push_back(CreateNoValveOutletStream(CreateVelocityBC(0.02)));
        streams.back()->name_ = "LPP";
        streams.push_back(CreateValveOutletStream(CreateConstantBoundaryValue(kAtm), 0.05));
        streams.back()->name_ = "Vent";
        AddFlowBC(std::move(streams), BoundaryLocation::TOP, BoundaryDirection::OUTLET, "Feed");
        dynamic_cast<ValveBoundaryStream&>(*Flow("Feed", BoundaryLocation::TOP).streams_.at(1)).check_valve_ = true;
      }
      AddStep("Blowdown", 10.0);
      AddPressureOutletBC(CreateConstantBoundaryValue(0.3 * kAtm), BoundaryLocation::BOTTOM, "Blowdown");
      AddClosedBC(BoundaryLocation::TOP, "Blowdown");
      AddStep("Pressurisation", 20.0);
      AddClosedBC(BoundaryLocation::BOTTOM, "Pressurisation");
      AddPressureInletBC(CreateLinearBoundaryValue(kAtm, 20.0), feed, 293, BoundaryLocation::TOP, "Pressurisation");
      ConfigureCycle(cycles);
      data.settings.css.method_ = CSSMethod::NONE;
      data.settings.warnings_.vacuum = false;
      data.cycle.EnsureIds();
      build::Connect(data, Lpp().endpoint_id_, Pressure("Pressurisation", BoundaryLocation::TOP).endpoint_id_);
    }

    IBoundaryStream& Lpp() { return *Flow("Feed", BoundaryLocation::TOP).streams_.at(0); }

    void Run(CoupledSolver& solver) {
      std::string error;
      ASSERT_TRUE(data.PerformChecks(error)) << error;
      ASSERT_GE(solver.IntegrateCycle(data.cycle, false), 0);
    }
  };

  // The donor falls below the receiver mid-step. From then on the LPP stream
  // carries nothing, apart from what the smooth closure lets through within
  // kConnectionPressureBand of the crossing; before it the draw replays the
  // receiver's intake.
  TEST_F(ConnectionPressureTest, NoDrawWhileTheDonorIsBelowTheReceiver) {
    LppCycle(3);
    CoupledSolver solver{data};
    Run(solver);

    const StreamData& donor = Lpp().stream_data_;
    const StreamData& receiver = Pressure("Pressurisation", BoundaryLocation::TOP).stream_data_;
    ASSERT_GT(donor.time_series_.size(), 10u);

    double below = 0.0, max_dp = 0.0, drawn = 0.0, drawn_below = 0.0, max_flow = 0.0;
    for (size_t k = 0; k < donor.time_series_.size(); ++k) {
      const double t = donor.time_series_[k];
      const double out = donor.GetOutflow(t);
      const double dp = receiver.GetPressure(t) - donor.pressure_boundary_[k];
      max_flow = std::max(max_flow, std::abs(out));
      if (dp > kConnectionPressureBand) {
        // Donor clearly below the receiver: nothing leaves it.
        EXPECT_LE(out, 1e-12) << "t = " << t << " s, dp = " << dp << " Pa";
      }
      if (k > 0) {
        const double dt = t - donor.time_series_[k - 1];
        const double o = 0.5 * (std::max(out, 0.0) + std::max(donor.GetOutflow(donor.time_series_[k - 1]), 0.0));
        drawn += o * dt;
        if (dp > 0.0) {
          below += dt;
          drawn_below += o * dt;
        }
      }
      max_dp = std::max(max_dp, dp);
    }
    // The case is the one the test is about: the donor ends well below the
    // receiver, and gas was drawn before the crossing.
    ASSERT_GT(below, 2.0) << "the donor never fell below the receiver";
    ASSERT_GT(max_dp, 1000.0);
    ASSERT_GT(drawn, 0.0);
    EXPECT_GT(max_flow, 0.0);
    EXPECT_LT(drawn_below, 1e-3 * drawn);
  }

} // namespace adsim
