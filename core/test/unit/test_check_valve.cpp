// test_check_valve.cpp
//
// The check valve of a valve stream (ValveBoundaryStream::check_valve_: no
// flow against the stream's direction, closure smoothed over
// Valve::kCheckValveBand) and the end-of-run warning for gas that flows back
// into the bed through an unconnected outlet (OperatingLimit::OutletBackflow).

#include <cmath>
#include <filesystem>
#include <string>
#include <vector>

#include "../src/pch.h"
#include "adsim/objects/valve.h"
#include "adsim/objects/deck_builder.h"
#include "adsim/solver/operating_limits.h"
#include "adsim/solver/coupled_solver.h"
#include "adsim/simulation.h"
#include "adsim/boundary_conditions/boundary_conditions.h"
#include "adsim/boundary_conditions/boundary_stream.h"
#include "adsim/boundary_conditions/step.h"

#include <gtest/gtest.h>

#include "../fixtures/user_data_test.h"

namespace adsim {

  // ---- the smoothed closure ------------------------------------------------

  // Shut against the reverse pressure drop, exactly the valve law beyond the
  // band, and continuous with a continuous slope at both ends of the band
  // (a hard max(dp, 0) has a kink at 0 that stalls Newton).
  TEST(CheckValve, ClosureIsSmoothAndExactOutsideTheBand) {
    const double band = Valve::kCheckValveBand;
    EXPECT_EQ(Valve::CheckValvePressureDrop(-1e5), 0.0);
    EXPECT_EQ(Valve::CheckValvePressureDrop(-1e-9), 0.0);
    EXPECT_EQ(Valve::CheckValvePressureDrop(0.0), 0.0);
    EXPECT_EQ(Valve::CheckValvePressureDrop(band), band);
    EXPECT_EQ(Valve::CheckValvePressureDrop(3.0 * band), 3.0 * band);
    EXPECT_EQ(Valve::CheckValvePressureDrop(2e5), 2e5);
    EXPECT_TRUE(std::isnan(Valve::CheckValvePressureDrop(std::nan(""))));

    const double h = 1e-6 * band;
    auto slope = [&](double x) {
      return (Valve::CheckValvePressureDrop(x + h) - Valve::CheckValvePressureDrop(x - h)) / (2.0 * h);
    };
    EXPECT_NEAR(slope(0.0), 0.0, 1e-5);   // matches the closed side
    EXPECT_NEAR(slope(band), 1.0, 1e-5);  // matches the open side
    // Monotonic and never more than the plain pressure drop inside the band.
    double previous = 0.0;
    for (int k = 1; k <= 100; ++k) {
      const double x = band * k / 100.0;
      const double s = Valve::CheckValvePressureDrop(x);
      EXPECT_GE(s, previous);
      EXPECT_LE(s, x + 1e-12);
      previous = s;
    }
  }

  // ---- runs ----------------------------------------------------------------

  class CheckValveTest : public UserDataTest {
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
    ValveBoundaryStream& Valve(const std::string& step, BoundaryLocation loc, size_t k = 0) {
      auto* v = dynamic_cast<ValveBoundaryStream*>(Flow(step, loc).streams_.at(k).get());
      EXPECT_NE(v, nullptr);
      return *v;
    }

    // A feed from the bottom with the top vented through a valve to 1 atm,
    // then a blowdown of the bottom to 0.9 atm: every feed after the first
    // starts with the bed below the vent's downstream pressure, so without a
    // check valve gas flows back in through the vent.
    void VentedCycle() {
      AddStep("Feed", 60.0);
      {
        std::vector<std::unique_ptr<IBoundaryStream>> streams;
        streams.push_back(CreateNoValveInletStream(CreateVelocityBC(0.1), feed, 293));
        AddFlowBC(std::move(streams), BoundaryLocation::BOTTOM, BoundaryDirection::INLET, "Feed");
      }
      {
        std::vector<std::unique_ptr<IBoundaryStream>> streams;
        streams.push_back(CreateValveOutletStream(CreateConstantBoundaryValue(kAtm), 0.05));
        streams.back()->name_ = "Vent";
        AddFlowBC(std::move(streams), BoundaryLocation::TOP, BoundaryDirection::OUTLET, "Feed");
      }
      AddStep("Blowdown", 10.0);
      AddPressureOutletBC(CreateConstantBoundaryValue(0.9 * kAtm), BoundaryLocation::BOTTOM, "Blowdown");
      AddClosedBC(BoundaryLocation::TOP, "Blowdown");
      ConfigureCycle(2);
      data.settings.css.method_ = CSSMethod::NONE;
    }

    void Run(CoupledSolver& solver) {
      std::string error;
      ASSERT_TRUE(data.PerformChecks(error)) << error;
      ASSERT_GE(solver.IntegrateCycle(data.cycle, false), 0);
    }

    static std::vector<const OperatingWarning*> Backflow(const OperatingLimitsReport& r) {
      std::vector<const OperatingWarning*> out;
      for (const auto& w : r.warnings)
        if (w.limit == OperatingLimit::OutletBackflow) out.push_back(&w);
      return out;
    }

    // Gross moles in through the bottom feed of `step`.
    double FeedMoles(const std::string& step) {
      const auto& ints = Flow(step, BoundaryLocation::BOTTOM).streams_.at(0)->stream_data_.integrals_;
      EXPECT_TRUE(ints.valid);
      return ints.moles_in;
    }
  };

  // Without a check valve the vent back-feeds the bed at the start of the
  // feed: one warning, for the feed step's vent, naming the moles and the
  // share of the feed and suggesting the check valve.
  TEST_F(CheckValveTest, BackflowThroughAnUnconnectedValveOutletIsReported) {
    VentedCycle();
    CoupledSolver solver{data};
    Run(solver);
    const auto& vent = Valve("Feed", BoundaryLocation::TOP).stream_data_.integrals_;
    ASSERT_TRUE(vent.valid);
    const double feed_moles = FeedMoles("Feed");
    ASSERT_GT(vent.moles_in, 10.0 * kOutletBackflowFraction * feed_moles);

    const auto r = CheckOperatingLimits(data, solver.GetLog());
    const auto w = Backflow(r);
    ASSERT_EQ(w.size(), 1u);
    EXPECT_EQ(w[0]->step, 0);
    EXPECT_NEAR(w[0]->value, vent.moles_in / feed_moles, 1e-9 * w[0]->value);
    EXPECT_EQ(w[0]->limit_value, kOutletBackflowFraction);
    EXPECT_NE(w[0]->message.find("\"Vent\""), std::string::npos) << w[0]->message;
    EXPECT_NE(w[0]->message.find("\"Feed\""), std::string::npos) << w[0]->message;
    EXPECT_NE(w[0]->message.find("check_valve"), std::string::npos) << w[0]->message;

    // settings.warn_outlet_backflow switches it off.
    data.settings.warnings_.outlet_backflow = false;
    EXPECT_TRUE(Backflow(CheckOperatingLimits(data, solver.GetLog())).empty());
  }

  // With the check valve the vent stays shut until the feed has raised the
  // bed above the downstream pressure: nothing flows back, the run
  // integrates, and no warning is given.
  TEST_F(CheckValveTest, CheckValveStopsTheBackflow) {
    VentedCycle();
    Valve("Feed", BoundaryLocation::TOP).check_valve_ = true;
    CoupledSolver solver{data};
    Run(solver);
    const auto& vent = Valve("Feed", BoundaryLocation::TOP).stream_data_.integrals_;
    ASSERT_TRUE(vent.valid);
    const double feed_moles = FeedMoles("Feed");
    EXPECT_LT(vent.moles_in, 1e-6 * feed_moles);
    EXPECT_GT(vent.moles_out, 0.02 * feed_moles);  // it still vents once open
    for (double f : Valve("Feed", BoundaryLocation::TOP).stream_data_.molar_flow_)
      EXPECT_GE(f, -1e-9 * feed_moles);  // top outlet: never into the bed
    EXPECT_TRUE(Backflow(CheckOperatingLimits(data, solver.GetLog())).empty());
  }

  // A connection carries flow both ways by design: an outlet connected to
  // another step's inlet that takes gas back is neither warned about nor
  // closed by its check valve. Here the receiving step runs at the feed
  // pressure and the donor after a blowdown, so all the gas through the
  // connection flows the "wrong" way: out of the receiver through its
  // inlet, into the donor through its outlet.
  TEST_F(CheckValveTest, ConnectedEndpointsAreNotCheckedOrClosed) {
    AddStep("Feed", 20.0);
    {
      std::vector<std::unique_ptr<IBoundaryStream>> streams;
      streams.push_back(CreateNoValveInletStream(CreateVelocityBC(0.05), feed, 293));
      AddFlowBC(std::move(streams), BoundaryLocation::BOTTOM, BoundaryDirection::INLET, "Feed");
    }
    {
      // Vented through a check valve, so this unconnected outlet does not
      // back-feed the bed after the donor step.
      std::vector<std::unique_ptr<IBoundaryStream>> streams;
      streams.push_back(CreateValveOutletStream(CreateConstantBoundaryValue(kAtm), 0.05));
      streams.back()->name_ = "Vent";
      AddFlowBC(std::move(streams), BoundaryLocation::TOP, BoundaryDirection::OUTLET, "Feed");
    }
    AddStep("Receiver", 10.0);
    AddClosedBC(BoundaryLocation::BOTTOM, "Receiver");
    {
      std::vector<std::unique_ptr<IBoundaryStream>> streams;
      streams.push_back(CreateValveInletStream(CreateConstantBoundaryValue(0.5 * kAtm), 0.05, feed, 293));
      streams.back()->name_ = "PE in";
      AddFlowBC(std::move(streams), BoundaryLocation::TOP, BoundaryDirection::INLET, "Receiver");
    }
    AddStep("Blowdown", 10.0);
    AddPressureOutletBC(CreateConstantBoundaryValue(0.5 * kAtm), BoundaryLocation::BOTTOM, "Blowdown");
    AddClosedBC(BoundaryLocation::TOP, "Blowdown");
    AddStep("Donor", 10.0);
    AddClosedBC(BoundaryLocation::BOTTOM, "Donor");
    {
      std::vector<std::unique_ptr<IBoundaryStream>> streams;
      streams.push_back(CreateValveOutletStream(CreateConstantBoundaryValue(kAtm), 0.05));
      streams.back()->name_ = "PE out";
      AddFlowBC(std::move(streams), BoundaryLocation::TOP, BoundaryDirection::OUTLET, "Donor");
    }
    ConfigureCycle(3);
    data.settings.css.method_ = CSSMethod::NONE;
    data.settings.warnings_.vacuum = false;
    data.cycle.EnsureIds();
    auto& out = Valve("Donor", BoundaryLocation::TOP);
    auto& in = Valve("Receiver", BoundaryLocation::TOP);
    build::Connect(data, out.endpoint_id_, in.endpoint_id_);
    Valve("Feed", BoundaryLocation::TOP).check_valve_ = true;
    out.check_valve_ = true;  // ignored: the streams are connected
    in.check_valve_ = true;

    CoupledSolver solver{data};
    Run(solver);
    const double feed_moles = FeedMoles("Feed");
    const auto& donor = out.stream_data_.integrals_;
    ASSERT_TRUE(donor.valid);
    // The backflow the test is about: well above the warning's threshold,
    // and not stopped by the check valve.
    ASSERT_GT(donor.moles_in, 10.0 * kOutletBackflowFraction * feed_moles);
    const auto r = CheckOperatingLimits(data, solver.GetLog());
    for (const auto* w : Backflow(r)) ADD_FAILURE() << w->message;
  }

  // ---- the deck ------------------------------------------------------------

  // The flag survives a save and load (ValveBoundaryStream cereal version 2).
  TEST_F(CheckValveTest, FlagRoundTrips) {
    VentedCycle();
    Valve("Feed", BoundaryLocation::TOP).check_valve_ = true;
    data.settings.warnings_.outlet_backflow = false;
    std::string text;
    SaveToString(data, text);
    EXPECT_NE(text.find("\"Check_Valve\": true"), std::string::npos);
    UserData loaded;
    LoadFromJsonText(loaded, text);
    auto& s = loaded.cycle.GetStep("Feed");
    auto* flow = dynamic_cast<FlowBoundaryCondition*>(s.top_boundary_condition_.get());
    ASSERT_NE(flow, nullptr);
    auto* v = dynamic_cast<ValveBoundaryStream*>(flow->streams_.at(0).get());
    ASSERT_NE(v, nullptr);
    EXPECT_TRUE(v->check_valve_);
    EXPECT_FALSE(loaded.settings.warnings_.outlet_backflow);
  }

  // A deck saved before the flag (version 1 valve streams) loads with every
  // check valve off and the warning on.
  TEST(CheckValveDeck, OlderDecksLoadWithTheCheckValveOff) {
    const auto deck = std::filesystem::current_path() / ".." / "examples" / "decks" / "Perez2019_D1.json";
    UserData data;
    ASSERT_NO_THROW(LoadFromJson(data, deck.string()));
    int valves = 0;
    for (auto& step : data.cycle.step_list_)
      for (auto* bc : {step.bottom_boundary_condition_.get(), step.top_boundary_condition_.get()})
        if (auto* flow = dynamic_cast<FlowBoundaryCondition*>(bc))
          for (auto& st : flow->streams_)
            if (auto* v = dynamic_cast<ValveBoundaryStream*>(st.get())) {
              ++valves;
              EXPECT_FALSE(v->check_valve_);
            }
    EXPECT_GT(valves, 0);
    EXPECT_TRUE(data.settings.warnings_.outlet_backflow);
  }

} // namespace adsim
