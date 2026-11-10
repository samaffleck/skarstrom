// test_connection_pressure.cpp
//
// Connections and the pressures of the two beds. A connected endpoint
// replays its partner's logged run: a set-flow stream replays the partner's
// flow, both ways, and a connected pressure boundary follows its own
// profile. Neither enforces the order of the pressures; the end-of-run
// warning connection_pressure reports gas that passes from the lower to the
// higher pressure.

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
#include "adsim/solver/operating_limits.h"

#include <gtest/gtest.h>

#include "../fixtures/user_data_test.h"

namespace adsim {

  // ---- light-product pressurisation ---------------------------------------

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
    PressureBoundaryCondition& Receiver() {
      auto* p = dynamic_cast<PressureBoundaryCondition*>(
          data.cycle.GetStep("Pressurisation").top_boundary_condition_.get());
      EXPECT_NE(p, nullptr);
      return *p;
    }
    IBoundaryStream& Lpp() { return *Flow("Feed", BoundaryLocation::TOP).streams_.at(0); }

    // Feed (the donor): fed at the bottom from 1 atm; at the top a set-flow
    // light-product stream "LPP", so the donor's top sits a little below
    // 1 atm. Blowdown: the bottom to 0.3 atm. Pressurisation (the
    // receiver): the top ramps from the bed pressure to 1.1 atm, above the
    // donor, so late in the step the receiver takes in gas while it is the
    // higher bed.
    void LppCycle(int cycles, bool connect = true) {
      AddStep("Feed", 20.0);
      AddPressureInletBC(CreateConstantBoundaryValue(kAtm), feed, 293, BoundaryLocation::BOTTOM, "Feed");
      {
        std::vector<std::unique_ptr<IBoundaryStream>> streams;
        streams.push_back(CreateNoValveOutletStream(CreateVelocityBC(0.02)));
        streams.back()->name_ = "LPP";
        AddFlowBC(std::move(streams), BoundaryLocation::TOP, BoundaryDirection::OUTLET, "Feed");
      }
      AddStep("Blowdown", 10.0);
      AddPressureOutletBC(CreateConstantBoundaryValue(0.3 * kAtm), BoundaryLocation::BOTTOM, "Blowdown");
      AddClosedBC(BoundaryLocation::TOP, "Blowdown");
      AddStep("Pressurisation", 20.0);
      AddClosedBC(BoundaryLocation::BOTTOM, "Pressurisation");
      AddPressureInletBC(CreateLinearBoundaryValue(1.1 * kAtm, 20.0), feed, 293, BoundaryLocation::TOP, "Pressurisation");
      ConfigureCycle(cycles);
      data.settings.css.method_ = CSSMethod::NONE;
      data.settings.warnings_.vacuum = false;
      data.cycle.EnsureIds();
      if (!connect) return;
      build::Connect(data, Lpp().endpoint_id_, Receiver().endpoint_id_);
    }

    // Pressure equalisation through a valve: Feed fills the bed to 3 atm;
    // in "Equalise down" (the donor) gas leaves the top through a valve
    // outlet connected to the top of "Equalise up" (the receiver), a
    // set-flow inlet; Blowdown takes the bed to 1 atm between them. The
    // donor's valve reads the receiver's logged pressure: the beds equalise
    // early in the step, and the donor's bottom vent then takes it below
    // the receiver, so the valve draws gas back.
    void EqualisationCycle(int cycles) {
      AddStep("Feed", 20.0);
      AddPressureInletBC(CreateConstantBoundaryValue(3.0 * kAtm), feed, 293, BoundaryLocation::BOTTOM, "Feed");
      AddClosedBC(BoundaryLocation::TOP, "Feed");
      AddStep("Equalise down", 20.0);
      // A vent at the bottom keeps the donor falling after the beds
      // have equalised, so the valve then draws gas back from the receiver.
      {
        std::vector<std::unique_ptr<IBoundaryStream>> streams;
        streams.push_back(CreateValveOutletStream(CreateConstantBoundaryValue(kAtm), 0.05));
        streams.back()->name_ = "Vent";
        AddFlowBC(std::move(streams), BoundaryLocation::BOTTOM, BoundaryDirection::OUTLET, "Equalise down");
      }
      {
        std::vector<std::unique_ptr<IBoundaryStream>> streams;
        streams.push_back(CreateValveOutletStream(CreateConstantBoundaryValue(2.0 * kAtm), 10.0));
        streams.back()->name_ = "PE out";
        AddFlowBC(std::move(streams), BoundaryLocation::TOP, BoundaryDirection::OUTLET, "Equalise down");
      }
      AddStep("Blowdown", 10.0);
      AddPressureOutletBC(CreateConstantBoundaryValue(kAtm), BoundaryLocation::BOTTOM, "Blowdown");
      AddClosedBC(BoundaryLocation::TOP, "Blowdown");
      AddStep("Equalise up", 20.0);
      AddClosedBC(BoundaryLocation::BOTTOM, "Equalise up");
      {
        std::vector<std::unique_ptr<IBoundaryStream>> streams;
        streams.push_back(CreateNoValveInletStream(CreateVelocityBC(0.0), feed, 293));
        streams.back()->name_ = "PE in";
        AddFlowBC(std::move(streams), BoundaryLocation::TOP, BoundaryDirection::INLET, "Equalise up");
      }
      ConfigureCycle(cycles);
      data.settings.css.method_ = CSSMethod::NONE;
      data.settings.warnings_.vacuum = false;
      data.cycle.EnsureIds();
      build::Connect(data, Flow("Equalise down", BoundaryLocation::TOP).streams_.at(0)->endpoint_id_,
                     Flow("Equalise up", BoundaryLocation::TOP).streams_.at(0)->endpoint_id_);
    }

    void Run(CoupledSolver& solver) {
      std::string error;
      ASSERT_TRUE(data.PerformChecks(error)) << error;
      ASSERT_GE(solver.IntegrateCycle(data.cycle, false), 0);
    }

    // Trapezoids of the parts of the logged flow that leave / enter the bed.
    static double Moles(const StreamData& d, double sign) {
      double n = 0.0;
      for (size_t k = 1; k < d.time_series_.size(); ++k)
        n += 0.5 * (std::max(sign * d.GetOutflow(d.time_series_[k]), 0.0)
                    + std::max(sign * d.GetOutflow(d.time_series_[k - 1]), 0.0))
             * (d.time_series_[k] - d.time_series_[k - 1]);
      return n;
    }
    static double MolesOut(const StreamData& d) { return Moles(d, 1.0); }
    static double MolesIn(const StreamData& d) { return Moles(d, -1.0); }

    std::vector<const OperatingWarning*> ConnectionWarnings(const OperatingLimitsReport& r) {
      std::vector<const OperatingWarning*> out;
      for (const auto& w : r.warnings)
        if (w.limit == OperatingLimit::ConnectionPressure) out.push_back(&w);
      return out;
    }
  };

  // Unconnected, the receiver follows its ramp to 1.1 atm, and no
  // connection warning is given.
  TEST_F(ConnectionPressureTest, UnconnectedPressureBoundaryFollowsItsRamp) {
    LppCycle(2, /*connect=*/false);
    CoupledSolver solver{data};
    Run(solver);
    const StreamData& r = Receiver().stream_data_;
    ASSERT_FALSE(r.time_series_.empty());
    EXPECT_NEAR(r.pressure_boundary_.back(), 1.1 * kAtm, 1.0);
    EXPECT_TRUE(ConnectionWarnings(CheckOperatingLimits(data, solver.GetLog())).empty());
  }

  // Connected: the receiver follows its own ramp to 1.1 atm, above the
  // donor, whose light-product stream replays what the receiver took in, so
  // the moles through the connection balance. The receiver then takes in
  // gas while it is the higher bed, and the end-of-run warning reports it.
  TEST_F(ConnectionPressureTest, ConnectedRampFollowsItsProfileAndUphillFlowIsReported) {
    LppCycle(10);
    CoupledSolver solver{data};
    Run(solver);
    const StreamData& donor = Lpp().stream_data_;
    const StreamData& receiver = Receiver().stream_data_;
    ASSERT_GT(receiver.time_series_.size(), 10u);
    EXPECT_NEAR(receiver.pressure_boundary_.back(), 1.1 * kAtm, 1.0);

    const double given = MolesOut(donor) - MolesIn(donor);
    const double taken = MolesIn(receiver) - MolesOut(receiver);
    ASSERT_GT(given, 0.0);
    EXPECT_NEAR(taken, given, 1e-3 * given) << "given " << given << " mol, taken " << taken << " mol";

    const auto r = CheckOperatingLimits(data, solver.GetLog());
    const auto w = ConnectionWarnings(r);
    ASSERT_EQ(w.size(), 1u);
    EXPECT_EQ(w[0]->step, 2);
    EXPECT_GT(w[0]->value, kConnectionReverseFlowFraction);
    EXPECT_NE(w[0]->message.find("\"Pressurisation\""), std::string::npos) << w[0]->message;
    EXPECT_NE(w[0]->message.find("\"LPP\""), std::string::npos) << w[0]->message;
    // settings.warn_connection_pressure switches it off.
    data.settings.warnings_.connection_pressure = false;
    EXPECT_TRUE(ConnectionWarnings(CheckOperatingLimits(data, solver.GetLog())).empty());
  }

  // A set-flow stream replays a connected valve's flow as it is, both ways:
  // the valve already sets its flow from the pressure difference. After the
  // beds equalise, the donor's valve draws small amounts of gas back while
  // its bed is a few Pa below the receiver's log. A one-sided closure on the
  // replay gave none of it back, and the connection made gas (2 % of the
  // equalisation moles on Vysyaraju2024_PSA_exp1). The moles the donor
  // gives, net, are the moles the receiver takes in.
  TEST_F(ConnectionPressureTest, ValvePartnerFlowIsReplayedBothWays) {
    EqualisationCycle(20);
    CoupledSolver solver{data};
    Run(solver);
    const StreamData& donor = Flow("Equalise down", BoundaryLocation::TOP).streams_.at(0)->stream_data_;
    const StreamData& receiver = Flow("Equalise up", BoundaryLocation::TOP).streams_.at(0)->stream_data_;
    ASSERT_GT(donor.time_series_.size(), 10u);
    ASSERT_GT(receiver.time_series_.size(), 10u);

    // The scenario: the beds equalised and the valve then drew gas back.
    const double given = MolesOut(donor) - MolesIn(donor);
    const double taken = MolesIn(receiver) - MolesOut(receiver);
    ASSERT_GT(given, 0.0);
    ASSERT_GT(MolesIn(donor), 1e-3 * given) << "the donor's valve drew no gas back";
    EXPECT_NEAR(taken, given, 1e-4 * given) << "given " << given << " mol, taken " << taken << " mol";
  }

} // namespace adsim
