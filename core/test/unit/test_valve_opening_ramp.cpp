// test_valve_opening_ramp.cpp
//
// The opening ramp of a valve (Valve::opening_start_, opening_ramp_time_):
// the valve law's flow is scaled by an opening fraction that rises, C1
// smooth, from the start opening to fully open over the ramp time. A
// pressurisation through a valve onto an evacuated bed lifts the bed less
// with a longer ramp.

#include <cmath>
#include <sstream>
#include <string>
#include <vector>

#include "../src/pch.h"
#include "adsim/objects/valve.h"
#include "adsim/solver/operating_limits.h"
#include "adsim/solver/coupled_solver.h"
#include "adsim/boundary_conditions/boundary_conditions.h"
#include "adsim/boundary_conditions/boundary_stream.h"
#include "adsim/boundary_conditions/step.h"

#include <cereal/archives/json.hpp>
#include <gtest/gtest.h>

#include "../fixtures/user_data_test.h"

namespace adsim {

  // ---- the multiplier ------------------------------------------------------

  // The defaults (start 1, ramp 0) give exactly 1 at every time, so a deck
  // without a ramp runs as before.
  TEST(ValveOpeningRamp, DefaultsAreFullyOpen) {
    Valve v;
    EXPECT_EQ(v.opening_start_.Get(), 1.0);
    EXPECT_EQ(v.opening_ramp_time_.Get(), 0.0);
    for (double t : {-1.0, 0.0, 1e-9, 1.0, 1e6}) EXPECT_EQ(v.OpeningFraction(t), 1.0);
    // A ramp with a start of 1, or a start with no ramp, is also fully open.
    EXPECT_EQ(Valve::OpeningFraction(0.5, 1.0, 10.0), 1.0);
    EXPECT_EQ(Valve::OpeningFraction(0.0, 0.1, 0.0), 1.0);
  }

  // The start value at the start of the step, fully open from the end of the
  // ramp, monotonic in between, and a continuous slope (zero) at both ends.
  TEST(ValveOpeningRamp, MultiplierIsSmoothAtBothEnds) {
    const double s = 0.05, T = 20.0;
    EXPECT_EQ(Valve::OpeningFraction(-1.0, s, T), s);
    EXPECT_EQ(Valve::OpeningFraction(0.0, s, T), s);
    EXPECT_NEAR(Valve::OpeningFraction(0.5 * T, s, T), 0.5 * (1.0 + s), 1e-15);
    EXPECT_EQ(Valve::OpeningFraction(T, s, T), 1.0);
    EXPECT_EQ(Valve::OpeningFraction(2.0 * T, s, T), 1.0);

    const double h = 1e-6 * T;
    auto slope = [&](double t) {
      return (Valve::OpeningFraction(t + h, s, T) - Valve::OpeningFraction(t - h, s, T)) / (2.0 * h);
    };
    EXPECT_NEAR(slope(0.0), 0.0, 1e-6);
    EXPECT_NEAR(slope(T), 0.0, 1e-6);
    EXPECT_NEAR(slope(0.5 * T), 1.5 * (1.0 - s) / T, 1e-6);  // peak slope of the cubic
    double previous = s;
    for (int k = 1; k <= 100; ++k) {
      const double f = Valve::OpeningFraction(T * k / 100.0, s, T);
      EXPECT_GE(f, previous);
      EXPECT_LE(f, 1.0);
      previous = f;
    }
  }

  // The ramp is saved with the valve and loaded back.
  TEST(ValveOpeningRamp, RoundTripsThroughTheDeck) {
    Valve v;
    v.opening_start_ = 0.2;
    v.opening_ramp_time_ = 15.0;
    std::stringstream ss;
    {
      cereal::JSONOutputArchive ar(ss);
      ar(cereal::make_nvp("valve", v));
    }
    Valve loaded;
    {
      cereal::JSONInputArchive ar(ss);
      ar(cereal::make_nvp("valve", loaded));
    }
    EXPECT_EQ(loaded.opening_start_.Get(), 0.2);
    EXPECT_EQ(loaded.opening_ramp_time_.Get(), 15.0);
  }

  // ---- a pressurisation run ------------------------------------------------

  class ValveOpeningRampTest : public UserDataTest {
  protected:
    static constexpr double kAtm = 101325.0;
    const std::unordered_map<std::string, double> feed{ {"CarbonDioxide", 0.15}, {"Nitrogen", 0.85} };

    void SetUp() override {
      UserDataTest::SetUp();
      AddComponents({ "CarbonDioxide", "Nitrogen" });
      // An evacuated bed.
      AddInitialConditions(0.1 * kAtm, 293, 293, { {"CarbonDioxide", 0.0}, {"Nitrogen", 1.0} });
      auto [layer_id, ads_id] = AddLayer("Adsorbent");
      GetLayer(layer_id).domain_.n_cells = 10;
      AddIsotherm({ "Nitrogen" }, ads_id, std::make_unique<HenryIsotherm>("Nitrogen"));
      AddIsotherm({ "CarbonDioxide" }, ads_id, std::make_unique<LangmuirIsotherm>("CarbonDioxide"));
      AddLinearDrivingForceModel(ads_id, "CarbonDioxide", 0.5);
      AddLinearDrivingForceModel(ads_id, "Nitrogen", 0.5);
      AddConstantInternalHeatTransferCoefficientModel(0.0, layer_id);
      AddConstantExternalHeatTransferCoefficientModel(0.0, layer_id);

      // One step: the bottom fed through a valve from a 1 atm header, the top closed.
      AddStep("Pressurisation", kStepTime);
      std::vector<std::unique_ptr<IBoundaryStream>> streams;
      streams.push_back(CreateValveInletStream(CreateConstantBoundaryValue(kAtm), kCv, feed, 293));
      streams.back()->name_ = "Header";
      AddFlowBC(std::move(streams), BoundaryLocation::BOTTOM, BoundaryDirection::INLET, "Pressurisation");
      AddClosedBC(BoundaryLocation::TOP, "Pressurisation");
      ConfigureCycle(1);
      data.settings.css.method_ = CSSMethod::NONE;
    }

    static constexpr double kStepTime = 60.0;
    static constexpr double kCv = 1.0;

    ValveBoundaryStream& Header() {
      auto& s = data.cycle.GetStep("Pressurisation");
      auto* flow = dynamic_cast<FlowBoundaryCondition*>(s.bottom_boundary_condition_.get());
      EXPECT_NE(flow, nullptr);
      auto* v = dynamic_cast<ValveBoundaryStream*>(flow->streams_.at(0).get());
      EXPECT_NE(v, nullptr);
      return *v;
    }

    // Peak approach to fluidisation and step-end pressure (max over the bed).
    std::pair<double, double> Run(double start, double ramp) {
      Header().valve_.opening_start_ = start;
      Header().valve_.opening_ramp_time_ = ramp;
      std::string error;
      EXPECT_TRUE(data.PerformChecks(error)) << error;
      CoupledSolver solver{data};
      EXPECT_GE(solver.IntegrateCycle(data.cycle, false), 0);
      const auto r = CheckOperatingLimits(data, solver.GetLog());
      EXPECT_EQ(r.steps.size(), 1u);
      return {r.steps.at(0).max_approach_to_fluidisation, r.steps.at(0).max_pressure};
    }
  };

  // A longer ramp from a small opening lowers the peak approach to
  // fluidisation; the bed still reaches the header pressure in the step.
  TEST_F(ValveOpeningRampTest, LongerRampLowersThePeakApproachToFluidisation) {
    const auto [open, p_open] = Run(1.0, 0.0);
    const auto [short_ramp, p_short] = Run(0.05, 5.0);
    const auto [long_ramp, p_long] = Run(0.05, 20.0);
    EXPECT_GT(open, 0.0);
    EXPECT_LT(short_ramp, 0.9 * open);
    EXPECT_LT(long_ramp, 0.9 * short_ramp);
    for (double p : {p_open, p_short, p_long}) EXPECT_GT(p, 0.95 * kAtm);
  }

  // An opening start outside (0, 1] or a negative ramp time fails the checks.
  TEST_F(ValveOpeningRampTest, ChecksRejectAnInvalidRamp) {
    std::string error;
    Header().valve_.opening_start_ = 0.0;
    EXPECT_FALSE(data.PerformChecks(error));
    EXPECT_NE(error.find("opening start"), std::string::npos) << error;
    Header().valve_.opening_start_ = 0.5;
    Header().valve_.opening_ramp_time_ = -1.0;
    error.clear();
    EXPECT_FALSE(data.PerformChecks(error));
    EXPECT_NE(error.find("ramp time"), std::string::npos) << error;
  }

} // end adsim namespace
