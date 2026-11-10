// test_node_mixing.cpp
//
// Mixing at an outlet flow boundary with several streams
// (FlowBoundaryCondition::Mixing). A light-product draw (a connected set-flow
// outlet) that takes more gas than the bed supplies through the same end is
// made up by a vent that flows back in. The gas the draw carries is then the
// mix of the bed's net outflow and the vent's own gas, so every component is
// conserved at the node. Before, every stream carried the bed-end gas: the
// back-fed gas was a copy of it, and the pressurisation step put back about
// twice what left the bed, so a trace at the top grew every cycle.

#include <algorithm>
#include <cmath>
#include <string>
#include <vector>

#include "../src/pch.h"
#include "adsim/objects/deck_builder.h"
#include "adsim/solver/coupled_solver.h"
#include "adsim/boundary_conditions/boundary_conditions.h"
#include "adsim/boundary_conditions/boundary_stream.h"
#include "adsim/boundary_conditions/step.h"

#include <gtest/gtest.h>

#include "../fixtures/user_data_test.h"

namespace adsim {

  class NodeMixingTest : public UserDataTest {
  protected:
    static constexpr double kAtm = 101325.0;
    static constexpr double kTrace = 1e-4;  // initial CO2 mole fraction in the bed

    void SetUp() override {
      UserDataTest::SetUp();
      AddComponents({ "CarbonDioxide", "Nitrogen" });
      // A CO2 trace in an N2 bed; nothing feeds CO2.
      AddInitialConditions(kAtm, 293, 293, { {"CarbonDioxide", kTrace}, {"Nitrogen", 1.0 - kTrace} });
      auto [layer_id, ads_id] = AddLayer("Adsorbent");
      GetLayer(layer_id).domain_.n_cells = 10;
      AddIsotherm({ "Nitrogen" }, ads_id, std::make_unique<HenryIsotherm>("Nitrogen"));
      AddIsotherm({ "CarbonDioxide" }, ads_id, std::make_unique<HenryIsotherm>("CarbonDioxide"));
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
    IBoundaryStream& Vent() { return *Flow("Feed", BoundaryLocation::TOP).streams_.at(0); }
    IBoundaryStream& Draw() { return *Flow("Feed", BoundaryLocation::TOP).streams_.at(1); }

    // Feed: N2 in at the bottom at a small set flow; at the top a vent valve
    // to 1 atm and the light-product draw "LPP", connected to the
    // Pressurisation pressure inlet. Blowdown: the bottom to 0.3 atm.
    // Pressurisation: the top ramps back to 1 atm with the light product.
    // The draw is larger than the feed, so the vent back-feeds the rest.
    // `vent_n2`: the vent's composition is N2 (else all zero).
    void LppCycle(int cycles, bool vent_n2) {
      const std::unordered_map<std::string, double> n2{ {"CarbonDioxide", 0.0}, {"Nitrogen", 1.0} };
      AddStep("Feed", 20.0);
      {
        std::vector<std::unique_ptr<IBoundaryStream>> streams;
        streams.push_back(CreateNoValveInletStream(CreateVelocityBC(0.002), n2, 293));
        streams.back()->name_ = "Feed in";
        AddFlowBC(std::move(streams), BoundaryLocation::BOTTOM, BoundaryDirection::INLET, "Feed");
      }
      {
        std::vector<std::unique_ptr<IBoundaryStream>> streams;
        streams.push_back(CreateValveOutletStream(CreateConstantBoundaryValue(kAtm), 10.0));
        streams.back()->name_ = "Vent";
        streams.back()->mole_fraction["CarbonDioxide"] = 0.0;
        streams.back()->mole_fraction["Nitrogen"] = vent_n2 ? 1.0 : 0.0;
        streams.push_back(CreateNoValveOutletStream(CreateVelocityBC(0.0)));
        streams.back()->name_ = "LPP";
        AddFlowBC(std::move(streams), BoundaryLocation::TOP, BoundaryDirection::OUTLET, "Feed");
      }
      AddStep("Blowdown", 10.0);
      AddPressureOutletBC(CreateConstantBoundaryValue(0.3 * kAtm), BoundaryLocation::BOTTOM, "Blowdown");
      AddClosedBC(BoundaryLocation::TOP, "Blowdown");
      AddStep("Pressurisation", 20.0);
      AddClosedBC(BoundaryLocation::BOTTOM, "Pressurisation");
      AddPressureInletBC(CreateLinearBoundaryValue(kAtm, 20.0), n2, 293, BoundaryLocation::TOP, "Pressurisation");
      ConfigureCycle(cycles);
      data.settings.css.method_ = CSSMethod::NONE;
      data.settings.warnings_.vacuum = false;
      data.cycle.EnsureIds();
      auto* receiver = dynamic_cast<PressureBoundaryCondition*>(
          data.cycle.GetStep("Pressurisation").top_boundary_condition_.get());
      ASSERT_NE(receiver, nullptr);
      build::Connect(data, Draw().endpoint_id_, receiver->endpoint_id_);
    }

    void Run(CoupledSolver& solver) {
      std::string error;
      ASSERT_TRUE(data.PerformChecks(error)) << error;
      ASSERT_GE(solver.IntegrateCycle(data.cycle, false), 0);
    }

    // Trapezoid of the logged flow out of the bed (total, or one component).
    static double MolesOut(const StreamData& d, const std::string& component = "") {
      const auto& t = d.time_series_;
      const double sign = d.location_ == BoundaryLocation::TOP ? 1.0 : -1.0;
      auto f = [&](size_t k) {
        return component.empty() ? d.GetOutflow(t[k]) : sign * d.component_molar_flow_.at(component)[k];
      };
      double n = 0.0;
      for (size_t k = 1; k < t.size(); ++k) n += 0.5 * (f(k) + f(k - 1)) * (t[k] - t[k - 1]);
      return n;
    }
    static double MaxMoleFraction(const StreamData& d, const std::string& component) {
      const auto& y = d.mole_fraction_boundary_.at(component);
      return y.empty() ? 0.0 : *std::max_element(y.begin(), y.end());
    }
  };

  // With an N2 vent the back-fed gas is N2: the draw carries the bed's net
  // outflow mixed with it, the pressurisation puts back only the CO2 that
  // left the bed, and the trace is flushed out by the N2 feed instead of
  // growing.
  TEST_F(NodeMixingTest, BackFedVentGasIsTheVentsOwnAndATraceDoesNotGrow) {
    LppCycle(12, /*vent_n2=*/true);
    CoupledSolver solver{data};
    Run(solver);
    const StreamData& vent = Vent().stream_data_;
    const StreamData& draw = Draw().stream_data_;
    ASSERT_GT(draw.time_series_.size(), 10u);

    // The scenario: the vent back-feeds the draw (net into the bed through
    // the vent), and the draw takes more than the bed gives through the top.
    const double vent_out = MolesOut(vent);
    const double draw_out = MolesOut(draw);
    ASSERT_LT(vent_out, -0.1 * draw_out) << "the vent does not back-feed the draw";

    // The back-fed gas is the vent's N2, in the component flows and the log.
    EXPECT_LE(std::abs(MolesOut(vent, "CarbonDioxide")), 1e-12);
    EXPECT_LE(MaxMoleFraction(vent, "CarbonDioxide"), 1e-12);

    // The node balance of CO2: what the draw carries is what leaves the bed
    // through the top, which is at most the bed-end fraction of the net
    // outflow; the trace therefore falls cycle on cycle.
    EXPECT_LT(MaxMoleFraction(draw, "CarbonDioxide"), 0.1 * kTrace);
  }

  // A vent with an all-zero composition has no gas of its own: the back-fed
  // gas is the bed-end gas, as without mixing, so the draw and the vent log
  // the same composition, and the trace is not flushed.
  TEST_F(NodeMixingTest, VentWithoutCompositionBackFeedsTheBedEndGas) {
    LppCycle(4, /*vent_n2=*/false);
    CoupledSolver solver{data};
    Run(solver);
    const StreamData& vent = Vent().stream_data_;
    const StreamData& draw = Draw().stream_data_;
    ASSERT_EQ(vent.time_series_.size(), draw.time_series_.size());
    ASSERT_LT(MolesOut(vent), -0.1 * MolesOut(draw)) << "the vent does not back-feed the draw";
    for (const std::string c : {"CarbonDioxide", "Nitrogen"}) {
      const auto& yv = vent.mole_fraction_boundary_.at(c);
      const auto& yd = draw.mole_fraction_boundary_.at(c);
      for (size_t k = 0; k < yv.size(); ++k)
        EXPECT_NEAR(yd[k], yv[k], 1e-12) << c << " at t = " << draw.time_series_[k] << " s";
    }
    EXPECT_GT(MaxMoleFraction(draw, "CarbonDioxide"), 0.1 * kTrace);
  }

} // namespace adsim
