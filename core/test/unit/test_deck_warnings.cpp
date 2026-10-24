// test_deck_warnings.cpp
//
// BoundaryPressureStepChangeWarnings: the once-per-run solver-log warning for a
// constant pressure boundary far from the pressure its step starts at (an
// instantaneous pressurisation or blowdown the integrator has to resolve).

#include <string>
#include <vector>

#include "../src/pch.h"
#include "adsim/solver/deck_warnings.h"
#include "adsim/solver/coupled_solver.h"
#include "adsim/boundary_conditions/boundary_conditions.h"
#include "adsim/boundary_conditions/boundary_stream.h"
#include "adsim/boundary_conditions/step.h"

#include <gtest/gtest.h>

#include "../fixtures/user_data_test.h"

namespace adsim {

  class DeckWarningsTest : public UserDataTest {
  protected:
    static constexpr double kAtm = 101325.0;
    const std::unordered_map<std::string, double> feed{ {"CarbonDioxide", 0.15}, {"Nitrogen", 0.85} };

    void SetUp() override {
      UserDataTest::SetUp();
      AddComponents({ "CarbonDioxide", "Nitrogen" });
      AddInitialConditions(kAtm, 293, 293, { {"CarbonDioxide", 0.0}, {"Nitrogen", 1.0} });
    }

    // The four-step cycle of the Ergun pressure-switch regression test: feed,
    // blowdown (top, exponential to 0.3 bar), evacuation (bottom, exponential
    // to 0.05 bar), then re-pressurisation through the bottom with `repressurise`.
    void FourStepCycle(std::unique_ptr<IBoundaryValue> repressurise) {
      AddStep("Feed", 25.0);
      {
        std::vector<std::unique_ptr<IBoundaryStream>> streams;
        streams.push_back(CreateNoValveInletStream(CreateVelocityBC(0.5), feed, 293));
        AddFlowBC(std::move(streams), BoundaryLocation::BOTTOM, BoundaryDirection::INLET, "Feed");
      }
      AddPressureOutletBC(CreateConstantBoundaryValue(kAtm), BoundaryLocation::TOP, "Feed");
      AddStep("Blowdown", 25.0);
      AddPressureOutletBC(CreateExponentialBoundaryValue(30000.0, 1.0), BoundaryLocation::TOP, "Blowdown");
      AddClosedBC(BoundaryLocation::BOTTOM, "Blowdown");
      AddStep("Evacuation", 25.0);
      AddClosedBC(BoundaryLocation::TOP, "Evacuation");
      AddPressureOutletBC(CreateExponentialBoundaryValue(5000.0, 1.0), BoundaryLocation::BOTTOM, "Evacuation");
      AddStep("Pressurisation", 25.0);
      AddClosedBC(BoundaryLocation::TOP, "Pressurisation");
      AddPressureInletBC(std::move(repressurise), feed, 293, BoundaryLocation::BOTTOM, "Pressurisation");
      ConfigureCycle(3);
    }

    std::vector<std::string> Warnings() {
      return BoundaryPressureStepChangeWarnings(data.cycle, data.ic.pressure.Get());
    }
  };

  TEST_F(DeckWarningsTest, ConstantPressureAfterEvacuationIsAStepChange) {
    FourStepCycle(CreateConstantBoundaryValue(kAtm));
    const auto w = Warnings();
    ASSERT_EQ(w.size(), 1u) << (w.empty() ? "" : w[0]);
    EXPECT_NE(w[0].find("\"Pressurisation\" (bottom)"), std::string::npos) << w[0];
    EXPECT_NE(w[0].find("1.01 bar"), std::string::npos) << w[0];
    EXPECT_NE(w[0].find("0.05 bar at the end of step \"Evacuation\""), std::string::npos) << w[0];
    EXPECT_NE(w[0].find("ramp"), std::string::npos) << w[0];
  }

  TEST_F(DeckWarningsTest, RampedPressurisationGivesNoWarning) {
    FourStepCycle(CreateExponentialBoundaryValue(kAtm, 1.0));
    EXPECT_TRUE(Warnings().empty());
  }

  TEST_F(DeckWarningsTest, FirstStepIsCheckedAgainstTheInitialPressure) {
    FourStepCycle(CreateExponentialBoundaryValue(kAtm, 1.0));
    EXPECT_TRUE(Warnings().empty());  // Feed's constant 1 atm outlet matches the 1 atm initial pressure
    data.ic.pressure = 0.5 * kAtm;
    const auto w = Warnings();
    ASSERT_EQ(w.size(), 1u);
    EXPECT_NE(w[0].find("\"Feed\" (top)"), std::string::npos) << w[0];
    EXPECT_NE(w[0].find("initial pressure 0.507 bar"), std::string::npos) << w[0];
  }

  TEST_F(DeckWarningsTest, LooksBackPastStepsWithoutAPressureBoundary) {
    // Vent (top, 1 atm) -> Hold (both ends closed) -> Vacuum (top, constant
    // 0.3 bar): the vacuum step starts near 1 atm, set two steps earlier.
    AddStep("Vent", 10.0);
    AddPressureOutletBC(CreateConstantBoundaryValue(kAtm), BoundaryLocation::TOP, "Vent");
    AddClosedBC(BoundaryLocation::BOTTOM, "Vent");
    AddStep("Hold", 10.0);
    AddClosedBC(BoundaryLocation::TOP, "Hold");
    AddClosedBC(BoundaryLocation::BOTTOM, "Hold");
    AddStep("Vacuum", 10.0);
    AddPressureOutletBC(CreateConstantBoundaryValue(30000.0), BoundaryLocation::TOP, "Vacuum");
    AddClosedBC(BoundaryLocation::BOTTOM, "Vacuum");
    ConfigureCycle(1);
    const auto w = Warnings();
    // Vacuum after Vent (1 -> 0.3 bar), and cyclically Vent after Vacuum (0.3 -> 1 bar).
    ASSERT_EQ(w.size(), 2u);
    EXPECT_NE(w[0].find("\"Vent\" (top)"), std::string::npos) << w[0];
    EXPECT_NE(w[0].find("at the end of step \"Vacuum\""), std::string::npos) << w[0];
    EXPECT_NE(w[1].find("\"Vacuum\" (top)"), std::string::npos) << w[1];
    EXPECT_NE(w[1].find("at the end of step \"Vent\""), std::string::npos) << w[1];
  }

  TEST_F(DeckWarningsTest, SmallDifferencesAreNotReported) {
    AddStep("A", 10.0);
    AddPressureOutletBC(CreateConstantBoundaryValue(kAtm), BoundaryLocation::TOP, "A");
    AddClosedBC(BoundaryLocation::BOTTOM, "A");
    AddStep("B", 10.0);
    AddPressureOutletBC(CreateConstantBoundaryValue(0.9 * kAtm), BoundaryLocation::TOP, "B");  // 10 %
    AddClosedBC(BoundaryLocation::BOTTOM, "B");
    ConfigureCycle(1);
    EXPECT_TRUE(Warnings().empty());
    EXPECT_EQ(BoundaryPressureStepChangeWarnings(data.cycle, data.ic.pressure.Get(), 0.05).size(), 2u);
  }

  // The solver logs the warning once per run, at LOW level, not once per cycle.
  TEST_F(DeckWarningsTest, SolverLogsTheWarningOncePerRun) {
    auto [layer_id, ads_id] = AddLayer("Adsorbent");
    GetLayer(layer_id).domain_.n_cells = 10;
    AddIsotherm({ "Nitrogen" }, ads_id, std::make_unique<HenryIsotherm>("Nitrogen"));
    AddIsotherm({ "CarbonDioxide" }, ads_id, std::make_unique<LangmuirIsotherm>("CarbonDioxide"));
    AddLinearDrivingForceModel(ads_id, "CarbonDioxide", 0.5);
    AddLinearDrivingForceModel(ads_id, "Nitrogen", 0.5);
    AddConstantInternalHeatTransferCoefficientModel(0.0, layer_id);
    AddConstantExternalHeatTransferCoefficientModel(0.0, layer_id);
    FourStepCycle(CreateConstantBoundaryValue(kAtm));
    for (auto& step : data.cycle.step_list_) step.step_time_ = 3.0;
    data.settings.css.method_ = CSSMethod::NONE;

    std::string error;
    ASSERT_TRUE(data.PerformChecks(error)) << error;
    CoupledSolver solver{data};
    int warnings = 0;
    solver.SetLogLevel(LogLevel::LOW);
    solver.SetMessageCallback([&warnings](LogLevel level, const std::string& msg) {
      if (msg.find("step change in boundary pressure") != std::string::npos) {
        EXPECT_EQ(level, LogLevel::LOW);
        ++warnings;
      }
    });
    ASSERT_GE(solver.IntegrateCycle(data.cycle, false), 0);
    EXPECT_EQ(warnings, 1);

    // Switched off in the deck's settings (warning_settings.h): not logged.
    data.settings.warnings_.pressure_step_change = false;
    warnings = 0;
    CoupledSolver quiet{data};
    quiet.SetLogLevel(LogLevel::LOW);
    quiet.SetMessageCallback([&warnings](LogLevel, const std::string& msg) {
      if (msg.find("step change in boundary pressure") != std::string::npos) ++warnings;
    });
    ASSERT_GE(quiet.IntegrateCycle(data.cycle, false), 0);
    EXPECT_EQ(warnings, 0);
  }

} // namespace adsim
