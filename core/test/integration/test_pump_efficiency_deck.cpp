// test_pump_efficiency_deck.cpp
//
// Regression test for F-04: a deck whose pump efficiency is 0 (e.g. the
// all-zero default polynomial) used to pass UserData::PerformChecks and then
// divided the pump power by zero. The checks must now reject any pump in use
// whose efficiency is not > 0 and <= 100 %, naming the stream and the step.

#include <string>

#include "../src/pch.h"
#include "adsim/boundary_conditions/boundary_conditions.h"
#include "adsim/boundary_conditions/boundary_stream.h"
#include "adsim/boundary_conditions/step.h"
#include "adsim/objects/compressor.h"

#include <gtest/gtest.h>

#include "../fixtures/user_data_test.h"

namespace adsim {

  // One-step deck: a compressor-driven pump inlet at the bottom and a pressure
  // outlet (with its optional vacuum-pump power calculation) at the top.
  class PumpEfficiencyDeckTest : public UserDataTest {
  protected:
    void SetUp() override {
      UserDataTest::SetUp();

      AddComponents({ "Oxygen", "Nitrogen" });

      auto [layer_id, ads_id] = AddLayer("13X");
      GetLayer(layer_id).domain_.n_cells = 10;
      AddConstantInternalHeatTransferCoefficientModel(0, layer_id);
      AddConstantExternalHeatTransferCoefficientModel(0, layer_id);

      AddInitialConditions(101325, 298, 298, { {"Oxygen", 0.21}, {"Nitrogen", 0.79} });

      AddStep("Feed", 5.0);
      auto stream = std::make_unique<PumpBoundaryStream>();
      stream->name_ = "Compressor";
      stream->temperature_ = 298;
      stream->mole_fraction = { {"Oxygen", 0.21}, {"Nitrogen", 0.79} };
      std::vector<std::unique_ptr<IBoundaryStream>> streams;
      streams.push_back(std::move(stream));
      AddFlowBC(std::move(streams), BoundaryLocation::BOTTOM, BoundaryDirection::INLET, "Feed");

      AddPressureOutletBC(CreateConstantBoundaryValue(101325), BoundaryLocation::TOP, "Feed");
      PressureOutlet().pump_ = CreatePump(PumpType::VACUUM_PUMP);
      PressureOutlet().pump_->use_pump_ = true;

      ConfigureCycle(1);
    }

    IPump& CompressorPump() {
      auto* fbc = dynamic_cast<FlowBoundaryCondition*>(
        data.cycle.GetStep("Feed").bottom_boundary_condition_.get());
      return *dynamic_cast<PumpBoundaryStream&>(*fbc->streams_.front()).pump_;
    }

    PressureBoundaryCondition& PressureOutlet() {
      return *dynamic_cast<PressureBoundaryCondition*>(
        data.cycle.GetStep("Feed").top_boundary_condition_.get());
    }

    static void SetConstantEfficiency(IPump& pump, double eta) {
      auto model = std::make_unique<ConstantPumpEfficiencyModel>();
      model->efficiency_ = eta;
      pump.efficiency_model_ = std::move(model);
    }
  };

  TEST_F(PumpEfficiencyDeckTest, ValidEfficienciesPassChecks) {
    std::string error;
    EXPECT_TRUE(data.PerformChecks(error)) << error;
  }

  TEST_F(PumpEfficiencyDeckTest, ZeroEfficiencyFailsChecksNamingStreamAndStep) {
    SetConstantEfficiency(CompressorPump(), 0.0);
    std::string error;
    EXPECT_FALSE(data.PerformChecks(error));
    EXPECT_NE(error.find("Pump efficiency must be > 0 and <= 100 %"), std::string::npos) << error;
    EXPECT_NE(error.find("stream 'Compressor'"), std::string::npos) << error;
    EXPECT_NE(error.find("step 'Feed'"), std::string::npos) << error;
  }

  TEST_F(PumpEfficiencyDeckTest, DefaultPolynomialEfficiencyFailsChecks) {
    CompressorPump().efficiency_model_ = CreatePumpEfficiencyModel(PumpEfficiencyType::POLYNOMIAL);
    std::string error;
    EXPECT_FALSE(data.PerformChecks(error));
    EXPECT_NE(error.find("Polynomial efficiency model gives 0 (0 %)"), std::string::npos) << error;
  }

  TEST_F(PumpEfficiencyDeckTest, EfficiencyAboveHundredPercentFailsChecks) {
    SetConstantEfficiency(*PressureOutlet().pump_, 1.2);
    std::string error;
    EXPECT_FALSE(data.PerformChecks(error));
    EXPECT_NE(error.find("Pump efficiency must be > 0 and <= 100 %"), std::string::npos) << error;
    EXPECT_NE(error.find("120 %"), std::string::npos) << error;
    EXPECT_NE(error.find("pressure boundary"), std::string::npos) << error;
    EXPECT_NE(error.find("step 'Feed'"), std::string::npos) << error;
  }

} // end adsim namespace
