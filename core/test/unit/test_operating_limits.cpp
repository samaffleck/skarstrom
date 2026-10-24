// test_operating_limits.cpp
//
// CheckOperatingLimits (the end-of-run report: pressure, wall temperature and
// approach to fluidisation per step of the last cycle, against the wall's
// pressure design) and WallDesignWarnings (the start-of-run thickness check).

#include <algorithm>
#include <string>
#include <vector>

#include "../src/pch.h"
#include "adsim/solver/operating_limits.h"
#include "adsim/solver/deck_warnings.h"
#include "adsim/solver/coupled_solver.h"
#include "adsim/objects/wall_design.h"
#include "adsim/simulation.h"
#include "adsim/boundary_conditions/boundary_conditions.h"
#include "adsim/boundary_conditions/boundary_stream.h"
#include "adsim/boundary_conditions/step.h"

#include <gtest/gtest.h>

#include "../fixtures/user_data_test.h"

namespace adsim {

  class OperatingLimitsTest : public UserDataTest {
  protected:
    static constexpr double kAtm = 101325.0;
    LayerID layer_id_{};
    const std::unordered_map<std::string, double> feed{ {"CarbonDioxide", 0.15}, {"Nitrogen", 0.85} };

    void SetUp() override {
      UserDataTest::SetUp();
      AddComponents({ "CarbonDioxide", "Nitrogen" });
      AddInitialConditions(kAtm, 293, 293, { {"CarbonDioxide", 0.0}, {"Nitrogen", 1.0} });
      auto [layer_id, ads_id] = AddLayer("Adsorbent");
      layer_id_ = layer_id;
      GetLayer(layer_id).domain_.n_cells = 10;
      AddIsotherm({ "Nitrogen" }, ads_id, std::make_unique<HenryIsotherm>("Nitrogen"));
      AddIsotherm({ "CarbonDioxide" }, ads_id, std::make_unique<LangmuirIsotherm>("CarbonDioxide"));
      AddLinearDrivingForceModel(ads_id, "CarbonDioxide", 0.5);
      AddLinearDrivingForceModel(ads_id, "Nitrogen", 0.5);
      AddConstantInternalHeatTransferCoefficientModel(0.0, layer_id);
      AddConstantExternalHeatTransferCoefficientModel(0.0, layer_id);
    }

    // Upward feed at `velocity` with the top at `feed_outlet`, then a downward
    // purge (top inlet) with the bottom at 0.9 atm.
    void Cycle(double velocity, double feed_outlet = kAtm) {
      AddStep("Feed", 30.0);
      {
        std::vector<std::unique_ptr<IBoundaryStream>> streams;
        streams.push_back(CreateNoValveInletStream(CreateVelocityBC(velocity), feed, 293));
        AddFlowBC(std::move(streams), BoundaryLocation::BOTTOM, BoundaryDirection::INLET, "Feed");
      }
      AddPressureOutletBC(CreateConstantBoundaryValue(feed_outlet), BoundaryLocation::TOP, "Feed");
      AddStep("Purge", 10.0);
      {
        std::vector<std::unique_ptr<IBoundaryStream>> streams;
        streams.push_back(CreateNoValveInletStream(CreateVelocityBC(0.02), feed, 293));
        AddFlowBC(std::move(streams), BoundaryLocation::TOP, BoundaryDirection::INLET, "Purge");
      }
      AddPressureOutletBC(CreateConstantBoundaryValue(0.9 * kAtm), BoundaryLocation::BOTTOM, "Purge");
      ConfigureCycle(2);
      data.settings.css.method_ = CSSMethod::NONE;
    }

    // Runs the cycle; returns the solver's LOW-level messages.
    std::vector<std::string> Run(CoupledSolver& solver) {
      std::string error;
      EXPECT_TRUE(data.PerformChecks(error)) << error;
      std::vector<std::string> messages;
      solver.SetLogLevel(LogLevel::LOW);
      solver.SetMessageCallback([&messages](LogLevel, const std::string& msg) { messages.push_back(msg); });
      EXPECT_GE(solver.IntegrateCycle(data.cycle, false), 0);
      return messages;
    }

    static int Count(const OperatingLimitsReport& r, OperatingLimit limit) {
      int n = 0;
      for (const auto& w : r.warnings) n += (w.limit == limit);
      return n;
    }
  };

  // Without a wall design only fluidisation is checked; a gentle feed keeps
  // the approach well below 1, and the extremes are reported for every step.
  TEST_F(OperatingLimitsTest, GentleFeedWithoutADesignGivesNoWarnings) {
    Cycle(0.02);
    CoupledSolver solver{data};
    Run(solver);
    const auto r = CheckOperatingLimits(data, solver.GetLog());
    EXPECT_EQ(r.cycle, 1);
    EXPECT_FALSE(r.design_checked);
    ASSERT_EQ(r.steps.size(), 2u);
    EXPECT_EQ(r.steps[0].step_name, "Feed");
    EXPECT_GT(r.steps[0].max_pressure, kAtm);
    EXPECT_GT(r.steps[0].max_approach_to_fluidisation, 0.0);
    EXPECT_LT(r.steps[0].max_approach_to_fluidisation, 1.0);
    EXPECT_LT(r.steps[1].min_pressure, 0.95 * kAtm);
    EXPECT_TRUE(r.warnings.empty());
  }

  // A fast upward feed lifts the bed: one warning, naming the feed step of the
  // last cycle, with where and when it peaks.
  TEST_F(OperatingLimitsTest, FastUpwardFeedFluidises) {
    Cycle(2.0);
    CoupledSolver solver{data};
    Run(solver);
    const auto r = CheckOperatingLimits(data, solver.GetLog());
    ASSERT_EQ(Count(r, OperatingLimit::Fluidisation), 1);
    const auto& w = r.warnings.front();
    EXPECT_EQ(w.step, 0);
    EXPECT_GT(w.value, 1.0);
    EXPECT_NE(w.message.find("step \"Feed\" of cycle 2"), std::string::npos) << w.message;
    const auto& s = r.steps[0];
    EXPECT_GE(s.fluidisation_position, 0.0);
    EXPECT_LE(s.fluidisation_position, GetLayer(layer_id_).domain_.length.Get());
    EXPECT_GE(s.fluidisation_time, 0.0);
    EXPECT_LE(s.fluidisation_time, 30.0 + 1e-9);
  }

  // With a design in use: the feed at 1.2 atm passes a 0.1 bar(g) design
  // pressure, the purge goes below atmospheric without vacuum service, and
  // the wall passes a design temperature set below the feed temperature.
  TEST_F(OperatingLimitsTest, DesignLimitsAreReportedPerStep) {
    Cycle(0.02, 1.2 * kAtm);
    data.bed.wall_.design_pressure_ = 1e4;
    data.bed.wall_.design_temperature_ = 280.0;
    CoupledSolver solver{data};
    Run(solver);
    auto r = CheckOperatingLimits(data, solver.GetLog());
    EXPECT_TRUE(r.design_checked);
    ASSERT_EQ(Count(r, OperatingLimit::DesignPressure), 1);
    EXPECT_GE(Count(r, OperatingLimit::Vacuum), 1);
    EXPECT_EQ(Count(r, OperatingLimit::DesignTemperature), 2);
    bool purge_vacuum = false;
    for (const auto& w : r.warnings) {
      if (w.limit == OperatingLimit::DesignPressure) EXPECT_EQ(w.step, 0);
      if (w.limit == OperatingLimit::Vacuum) purge_vacuum |= (w.step == 1);
    }
    EXPECT_TRUE(purge_vacuum);

    data.bed.wall_.vacuum_service_ = true;
    data.bed.wall_.design_pressure_ = 10e5;
    data.bed.wall_.design_temperature_ = 350.0;
    r = CheckOperatingLimits(data, solver.GetLog());
    EXPECT_TRUE(r.warnings.empty());

    // A pressure within 2 kPa of the design pressure is not reported.
    data.bed.wall_.design_pressure_ = r.steps[0].max_pressure - kAtm - 1000.0;
    r = CheckOperatingLimits(data, solver.GetLog());
    EXPECT_TRUE(r.warnings.empty());
  }

  // The solver logs the warnings once, at the end of the run, at LOW level.
  TEST_F(OperatingLimitsTest, SolverLogsTheWarningsAtTheEndOfTheRun) {
    Cycle(2.0);
    CoupledSolver solver{data};
    const auto messages = Run(solver);
    int found = 0;
    size_t at = 0, completed = 0;
    for (size_t i = 0; i < messages.size(); ++i) {
      if (messages[i].find("approach to fluidisation reaches") != std::string::npos) { ++found; at = i; }
      if (messages[i].rfind("Solver completed", 0) == 0) completed = i;
    }
    EXPECT_EQ(found, 1);
    EXPECT_GT(at, completed);
  }

  // ---- WallDesignWarnings -----------------------------------------------

  TEST_F(OperatingLimitsTest, ThinWallIsReportedBeforeTheRun) {
    auto& wall = data.bed.wall_;
    wall.ApplyMaterial(WallMaterial::SS316L);
    wall.thickness_ = 0.0001;
    EXPECT_TRUE(WallDesignWarnings(data.bed).empty());  // no design in use

    wall.design_pressure_ = 10e5;
    const auto w = WallDesignWarnings(data.bed);
    ASSERT_EQ(w.size(), 1u);
    EXPECT_NE(w[0].find("thinner than"), std::string::npos) << w[0];

    wall.thickness_ = DesignBedWall(data.bed).required_thickness;
    EXPECT_TRUE(WallDesignWarnings(data.bed).empty());

    wall.ApplyMaterial(WallMaterial::Custom);
    ASSERT_EQ(WallDesignWarnings(data.bed).size(), 1u);
    EXPECT_NE(WallDesignWarnings(data.bed)[0].find("cannot be checked"), std::string::npos);
  }

  // The design fields survive a save and load (BedWall cereal version 4).
  TEST_F(OperatingLimitsTest, WallDesignFieldsRoundTrip) {
    auto& wall = data.bed.wall_;
    wall.design_code_ = WallDesignCode::ASME_B31_3;
    wall.design_pressure_ = 6e5;
    wall.corrosion_allowance_ = 0.0015;
    wall.vacuum_service_ = true;
    std::string text;
    SaveToString(data, text);
    UserData loaded;
    LoadFromJsonText(loaded, text);
    EXPECT_EQ(loaded.bed.wall_.design_code_, WallDesignCode::ASME_B31_3);
    EXPECT_DOUBLE_EQ(loaded.bed.wall_.design_pressure_.Get(), 6e5);
    EXPECT_DOUBLE_EQ(loaded.bed.wall_.corrosion_allowance_.Get(), 0.0015);
    EXPECT_TRUE(loaded.bed.wall_.vacuum_service_);
  }

  // ---- WarningSettings ---------------------------------------------------

  // A switched-off warning is neither listed nor logged; the summaries stay.
  TEST_F(OperatingLimitsTest, SwitchedOffWarningsAreNotReported) {
    Cycle(2.0, 1.2 * kAtm);
    data.bed.wall_.design_pressure_ = 1e4;
    data.settings.warnings_.fluidisation = false;
    data.settings.warnings_.design_pressure = false;
    CoupledSolver solver{data};
    const auto messages = Run(solver);
    for (const auto& m : messages) {
      EXPECT_EQ(m.find("approach to fluidisation reaches"), std::string::npos) << m;
      EXPECT_EQ(m.find("above the wall's design pressure"), std::string::npos) << m;
    }
    auto r = CheckOperatingLimits(data, solver.GetLog());
    EXPECT_EQ(Count(r, OperatingLimit::Fluidisation), 0);
    EXPECT_EQ(Count(r, OperatingLimit::DesignPressure), 0);
    EXPECT_GE(Count(r, OperatingLimit::Vacuum), 1);  // still on
    EXPECT_GT(r.steps[0].max_approach_to_fluidisation, 1.0);  // the summary is complete

    data.settings.warnings_.fluidisation = true;
    data.settings.warnings_.design_pressure = true;
    r = CheckOperatingLimits(data, solver.GetLog());
    EXPECT_EQ(Count(r, OperatingLimit::Fluidisation), 1);
    EXPECT_EQ(Count(r, OperatingLimit::DesignPressure), 1);
  }

  // The start-of-run wall-thickness warning has its own switch.
  TEST_F(OperatingLimitsTest, WallThicknessWarningCanBeSwitchedOff) {
    Cycle(0.02);
    data.bed.wall_.ApplyMaterial(WallMaterial::SS316L);
    data.bed.wall_.thickness_ = 0.0001;
    data.bed.wall_.design_pressure_ = 10e5;
    auto thin = [](const std::vector<std::string>& m) {
      return std::count_if(m.begin(), m.end(), [](const std::string& s) { return s.find("thinner than") != std::string::npos; });
    };
    {
      CoupledSolver solver{data};
      EXPECT_EQ(thin(Run(solver)), 1);
    }
    data.settings.warnings_.wall_thickness = false;
    CoupledSolver solver{data};
    EXPECT_EQ(thin(Run(solver)), 0);
  }

  TEST_F(OperatingLimitsTest, WarningSettingsRoundTrip) {
    data.settings.warnings_.fluidisation = false;
    data.settings.warnings_.pressure_step_change = false;
    std::string text;
    SaveToString(data, text);
    UserData loaded;
    LoadFromJsonText(loaded, text);
    EXPECT_FALSE(loaded.settings.warnings_.fluidisation);
    EXPECT_FALSE(loaded.settings.warnings_.pressure_step_change);
    EXPECT_TRUE(loaded.settings.warnings_.vacuum);
  }

} // namespace adsim
