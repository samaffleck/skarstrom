// test_connection_pressure.cpp
//
// Gas passes through a connection only from the higher to the lower
// pressure. A connected set-flow stream scales the flow it replays by
// ConnectionFlowWeight of the pressure difference at the same instant; a
// connected pressure boundary caps its imposed pressure one band inside the
// partner's logged pressure (ConnectionCappedPressure), so the pair
// balances once the replay has settled. Both closures are smoothed over
// kConnectionPressureBand.

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
    EXPECT_NEAR(slope(0.0), 0.0, 1e-5 / band);
    EXPECT_NEAR(slope(band), 0.0, 1e-5 / band);
    double previous = 0.0;
    for (int k = 1; k <= 100; ++k) {
      const double w = ConnectionFlowWeight(band * k / 100.0);
      EXPECT_GE(w, previous);
      EXPECT_LE(w, 1.0);
      previous = w;
    }
  }

  // The target is kept exactly while it is a band inside the bound, the
  // bound at or past it, and the blend between never passes the bound and
  // is C1 at both ends.
  TEST(ConnectionCappedPressure, KeepsTheTargetOnTheRightSideOfThePartner) {
    const double band = kConnectionPressureBand;
    const double p = 1e5;
    // Inlet: at most the partner.
    EXPECT_EQ(ConnectionCappedPressure(p - band, p, true), p - band);
    EXPECT_EQ(ConnectionCappedPressure(0.5 * p, p, true), 0.5 * p);
    EXPECT_EQ(ConnectionCappedPressure(p, p, true), p);
    EXPECT_EQ(ConnectionCappedPressure(2 * p, p, true), p);
    // Outlet: at least the partner.
    EXPECT_EQ(ConnectionCappedPressure(p + band, p, false), p + band);
    EXPECT_EQ(ConnectionCappedPressure(0.5 * p, p, false), p);
    EXPECT_TRUE(std::isnan(ConnectionCappedPressure(std::nan(""), p, true)));
    const double h = 1e-4;
    auto slope = [&](double target) {
      return (ConnectionCappedPressure(target + h, p, true) - ConnectionCappedPressure(target - h, p, true)) / (2 * h);
    };
    EXPECT_NEAR(slope(p - band), 1.0, 1e-4);
    EXPECT_NEAR(slope(p), 0.0, 1e-4);
    for (int k = 0; k <= 100; ++k) {
      const double target = p - band + 2.0 * band * k / 100.0;
      const double c = ConnectionCappedPressure(target, p, true);
      EXPECT_LE(c, p);
      EXPECT_GE(c, std::min(target, p) - 1e-9);
    }
  }

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
    // donor, so the cap acts before the end of the step.
    // `receiver_stream`: the receiver's top is a set-flow inlet stream and
    // its bottom a pressure outlet ramped to 1.2 atm instead, so the receiver rises above
    // the donor mid-step (a stream-stream connection).
    void LppCycle(int cycles, bool connect = true, bool receiver_stream = false) {
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
      if (receiver_stream) {
        // A pressure outlet above the bed: gas flows in through it.
        AddPressureOutletBC(CreateLinearBoundaryValue(1.2 * kAtm, 20.0), BoundaryLocation::BOTTOM, "Pressurisation");
        std::vector<std::unique_ptr<IBoundaryStream>> streams;
        streams.push_back(CreateNoValveInletStream(CreateVelocityBC(0.02), feed, 293));
        streams.back()->name_ = "LPP in";
        AddFlowBC(std::move(streams), BoundaryLocation::TOP, BoundaryDirection::INLET, "Pressurisation");
      } else {
        AddClosedBC(BoundaryLocation::BOTTOM, "Pressurisation");
        AddPressureInletBC(CreateLinearBoundaryValue(1.1 * kAtm, 20.0), feed, 293, BoundaryLocation::TOP, "Pressurisation");
      }
      ConfigureCycle(cycles);
      data.settings.css.method_ = CSSMethod::NONE;
      data.settings.warnings_.vacuum = false;
      data.cycle.EnsureIds();
      if (!connect) return;
      const size_t receiver_id = receiver_stream
          ? Flow("Pressurisation", BoundaryLocation::TOP).streams_.at(0)->endpoint_id_
          : Receiver().endpoint_id_;
      build::Connect(data, Lpp().endpoint_id_, receiver_id);
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
    for (size_t k = 0; k < r.time_series_.size(); ++k)
      EXPECT_EQ(Receiver().ImposedPressure(r.time_series_[k]), Receiver().pressure_->GetValue(r.time_series_[k]));
    EXPECT_NEAR(r.pressure_boundary_.back(), 1.1 * kAtm, 1.0);
    EXPECT_TRUE(ConnectionWarnings(CheckOperatingLimits(data, solver.GetLog())).empty());
  }

  // Connected: the receiver's ramp is capped a band below the donor's
  // logged pressure, so it never rises above the bed that feeds it and does
  // not reach its target; the donor replays what it took in, and once the
  // replay has settled the moles the receiver takes in equal those the
  // donor gives. The end-of-run warning reports the capped ramp.
  TEST_F(ConnectionPressureTest, ConnectedRampIsCappedAndTheMolesBalance) {
    LppCycle(40);
    CoupledSolver solver{data};
    Run(solver);
    const StreamData& donor = Lpp().stream_data_;
    const StreamData& receiver = Receiver().stream_data_;
    ASSERT_GT(receiver.time_series_.size(), 10u);

    double max_above = -1e30;
    for (size_t k = 0; k < receiver.time_series_.size(); ++k)
      max_above = std::max(max_above, receiver.pressure_boundary_[k] - donor.GetPressure(receiver.time_series_[k]));
    // The face carries the solver tolerance around the imposed value.
    EXPECT_LE(max_above, -kConnectionPressureBand + 1.0);
    // The cap acted: the ramp ends at 1.1 atm, the receiver near the donor.
    EXPECT_LT(receiver.pressure_boundary_.back(), kAtm);

    const double given = MolesOut(donor);
    const double taken = MolesIn(receiver);
    ASSERT_GT(given, 0.0);
    EXPECT_NEAR(taken, given, 1e-3 * given) << "given " << given << " mol, taken " << taken << " mol";
    EXPECT_LT(MolesIn(donor), 1e-3 * given);

    const auto r = CheckOperatingLimits(data, solver.GetLog());
    const auto w = ConnectionWarnings(r);
    ASSERT_EQ(w.size(), 1u);
    EXPECT_EQ(w[0]->step, 2);
    EXPECT_GT(w[0]->value, 2000.0);
    EXPECT_NE(w[0]->message.find("\"Pressurisation\""), std::string::npos) << w[0]->message;
    EXPECT_NE(w[0]->message.find("\"LPP\""), std::string::npos) << w[0]->message;
    // settings.warn_connection_pressure switches it off.
    data.settings.warnings_.connection_pressure = false;
    EXPECT_TRUE(ConnectionWarnings(CheckOperatingLimits(data, solver.GetLog())).empty());
  }

  // Stream-stream: each side scales the other's flow by ConnectionFlowWeight,
  // so nothing leaves the donor while it is below the receiver (beyond the
  // band), and the receiver takes in nothing while it is above the donor.
  TEST_F(ConnectionPressureTest, SetFlowPartnerIsNotDrawnFromBelowTheReceiver) {
    LppCycle(4, /*connect=*/true, /*receiver_stream=*/true);
    CoupledSolver solver{data};
    Run(solver);
    const StreamData& donor = Lpp().stream_data_;
    const StreamData& receiver = Flow("Pressurisation", BoundaryLocation::TOP).streams_.at(0)->stream_data_;
    ASSERT_GT(donor.time_series_.size(), 10u);

    // The logs compared are of the same cycle, while each side read the
    // other's previous one, so the check allows a band of lag.
    double below = 0.0, drawn = 0.0, drawn_below = 0.0;
    for (size_t k = 1; k < donor.time_series_.size(); ++k) {
      const double t = donor.time_series_[k];
      const double dt = t - donor.time_series_[k - 1];
      const double dp = receiver.GetPressure(t) - donor.pressure_boundary_[k];
      const double out = 0.5 * (std::max(donor.GetOutflow(t), 0.0) + std::max(donor.GetOutflow(t - dt), 0.0));
      drawn += out * dt;
      if (dp >= kConnectionPressureBand) EXPECT_LE(donor.GetOutflow(t), 1e-12) << "t = " << t << " s, dp = " << dp << " Pa";
      if (dp > 0.0) { below += dt; drawn_below += out * dt; }
    }
    ASSERT_GT(below, 1.0) << "the donor never fell below the receiver";
    ASSERT_GT(drawn, 0.0);
    EXPECT_LT(drawn_below, 1e-3 * drawn);
    for (size_t k = 0; k < receiver.time_series_.size(); ++k) {
      const double t = receiver.time_series_[k];
      if (receiver.pressure_boundary_[k] >= donor.GetPressure(t) + kConnectionPressureBand)
        EXPECT_LE(-receiver.GetOutflow(t), 1e-12) << "t = " << t << " s";
    }
    data.settings.warnings_.outlet_backflow = false;  // the bottom outlet feeds the receiver on purpose
    EXPECT_TRUE(ConnectionWarnings(CheckOperatingLimits(data, solver.GetLog())).empty());
  }

} // namespace adsim
