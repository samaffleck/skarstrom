// test_ergun_pressure_switch.cpp
//
// Regression test for an Ergun layer at a step switch that moves the boundary
// pressure by a large factor: a four-step PSA cycle (feed, blowdown to
// 0.3 bar, evacuation to 0.05 bar, re-pressurisation with the feed at a
// constant 1 bar) on a bed with a small bottom void. The deck is the
// web-workspace defaults a new user gets (Ergun momentum, Henry / Langmuir
// isotherms, constant LDF), and it is not a demanding problem: with the Darcy
// (Kozeny-Carman) pressure drop it never needs a solver retry.
//
// With Ergun, the start of every re-pressurisation step failed IDA's corrector
// down to the minimum step. The boundary face's Ergun coefficients took the
// gas density and viscosity from the ghost cell, whose state is the boundary
// condition's mirror (P = 2 Pf - P1): at the switch from 0.05 bar outlet to
// 1 bar inlet the ghost density jumps about forty-fold, and the inertial
// coefficient's 1 / sqrt(rho) dependence on it, linearised at the pre-switch
// state, sent Newton's second iterate further away than its first on every
// attempt. The restart from a consistent initial condition rescued each cycle
// until its own Newton line search failed (cycle 10 natively; cycle 74 in the
// browser). The boundary faces now take their transport properties from the
// interior neighbour (ErgunFaceCoefficients in adsorbent_equations.cpp).
//
// Both momentum models must complete the cycles with no step retry.

#include <string>

#include "../src/pch.h"
#include "adsim/solver/coupled_solver.h"
#include "adsim/boundary_conditions/boundary_conditions.h"
#include "adsim/boundary_conditions/boundary_stream.h"
#include "adsim/boundary_conditions/step.h"
#include "adsim/equations/heat_transfer_model.h"
#include "adsim/equations/isotherm.h"

#include <gtest/gtest.h>

#include "../fixtures/user_data_test.h"

namespace adsim {

  class ErgunPressureSwitchTest : public UserDataTest {
  protected:
    static constexpr double kStepTime = 25.0;       // [s], every step
    static constexpr double kFeedPressure = 101325; // [Pa]
    static constexpr double kBlowdownPressure = 30000.0;   // [Pa]
    static constexpr double kEvacuationPressure = 5000.0;  // [Pa]
    static constexpr double kFeedVelocity = 0.5;    // [m/s]
    static constexpr double kTemperature = 293.0;   // [K]
    // Natively the unpatched boundary face survived nine cycles on the
    // restart before its consistent-IC solve failed; a dozen covers that with
    // margin, and the zero-retry check below catches the fault in cycle 1.
    static constexpr int kCycles = 12;

    void Run(MomentumModel momentum) {
      const std::string CO2 = "CarbonDioxide";
      const std::string N2 = "Nitrogen";
      AddComponents({ CO2, N2 });
      const std::unordered_map<std::string, double> feed{ {CO2, 0.15}, {N2, 0.85} };

      auto [layer_id, ads_id] = AddLayer("Adsorbent 1");
      auto& layer = GetLayer(layer_id);
      layer.momentum_model_ = momentum;
      layer.domain_.n_cells = 30;
      layer.domain_.length = 1.0;
      layer.eb_ = 0.37;
      layer.particle_diameter_ = 1.0e-3;
      std::get<CylinderBed>(data.bed.shape_).diameter = 0.1;

      auto& ads = GetAdsorbent(ads_id);
      ads.solid_density_ = 1200.0;
      ads.heat_capacity_ = 500.0;
      ads.thermal_conductivity_ = 0.5;
      ads.ep_ = 0.35;
      // The default Henry (N2) and Langmuir (CO2) parameters of a new adsorbent.
      AddIsotherm({ N2 }, ads_id, std::make_unique<HenryIsotherm>(N2));
      AddIsotherm({ CO2 }, ads_id, std::make_unique<LangmuirIsotherm>(CO2));
      AddLinearDrivingForceModel(ads_id, CO2, 0.5);
      AddLinearDrivingForceModel(ads_id, N2, 0.5);
      SetHOA(ads_id, CO2, 0.0);
      SetHOA(ads_id, N2, 0.0);
      AddConstantInternalHeatTransferCoefficientModel(0.0, layer_id);
      AddConstantExternalHeatTransferCoefficientModel(0.0, layer_id);

      // A small bottom header: the first cell the feed meets is a void cell,
      // which made the unpatched failure terminal sooner (the ghost cell then
      // supplies the whole boundary face's coefficients).
      data.bed.bottom_void_.enabled_ = true;
      data.bed.bottom_void_.volume_ = 1.0e-4;

      AddInitialConditions(kFeedPressure, kTemperature, kTemperature, { {CO2, 0.0}, {N2, 1.0} });

      // 1. Feed: velocity inlet at the bottom, product out of the top at 1 bar.
      AddStep("Feed", kStepTime);
      {
        std::vector<std::unique_ptr<IBoundaryStream>> streams;
        streams.push_back(CreateNoValveInletStream(CreateVelocityBC(kFeedVelocity), feed, kTemperature));
        AddFlowBC(std::move(streams), BoundaryLocation::BOTTOM, BoundaryDirection::INLET, "Feed");
      }
      AddPressureOutletBC(CreateConstantBoundaryValue(kFeedPressure), BoundaryLocation::TOP, "Feed");

      // 2. Blowdown through the top to 0.3 bar, bottom closed.
      AddStep("Blowdown", kStepTime);
      AddPressureOutletBC(CreateExponentialBoundaryValue(kBlowdownPressure, 1.0), BoundaryLocation::TOP, "Blowdown");
      AddClosedBC(BoundaryLocation::BOTTOM, "Blowdown");

      // 3. Evacuation through the bottom to 0.05 bar, top closed.
      AddStep("Evacuation", kStepTime);
      AddClosedBC(BoundaryLocation::TOP, "Evacuation");
      AddPressureOutletBC(CreateExponentialBoundaryValue(kEvacuationPressure, 1.0), BoundaryLocation::BOTTOM, "Evacuation");

      // 4. Re-pressurisation with the feed at a constant 1 bar: the boundary
      // pressure of the bottom face jumps twenty-fold at the switch.
      AddStep("Pressurisation", kStepTime);
      AddClosedBC(BoundaryLocation::TOP, "Pressurisation");
      AddPressureInletBC(CreateConstantBoundaryValue(kFeedPressure), feed, kTemperature,
                         BoundaryLocation::BOTTOM, "Pressurisation");

      ConfigureCycle(kCycles);
      data.settings.css.method_ = CSSMethod::NONE;

      std::string error;
      ASSERT_TRUE(data.PerformChecks(error)) << error;

      CoupledSolver solver{data};
      // Each step retry is a recoverable IDA failure (logged at MEDIUM) that
      // the step survived only by being re-attempted; a step switch of this
      // cycle must not need one.
      int step_retries = 0;
      solver.SetLogLevel(LogLevel::MEDIUM);
      solver.SetMessageCallback([&step_retries](LogLevel, const std::string& msg) {
        if (msg.find("IDA failure") != std::string::npos) {
          ++step_retries;
          std::cout << "[ErgunPressureSwitch] " << msg << std::endl;
        }
      });
      const int retval = solver.IntegrateCycle(data.cycle, false);
      ASSERT_GE(retval, 0) << "Solver did not complete ("
                           << MomentumModelList[static_cast<int>(momentum)] << ")";
      EXPECT_EQ(step_retries, 0) << "step switches needed solver retries ("
                                 << MomentumModelList[static_cast<int>(momentum)] << ")";

      const auto& cycles = solver.GetLog().sys_data.cycle_index_set_;
      ASSERT_FALSE(cycles.empty());
      EXPECT_EQ(*cycles.rbegin() + 1, kCycles);
      EXPECT_TRUE(IsMassConserved(solver.GetLog()));
      EXPECT_TRUE(IsMolesConserved(solver.GetLog()));
    }
  };

  TEST_F(ErgunPressureSwitchTest, ConstantPressureRepressurisationErgun) {
    Run(MomentumModel::Ergun);
  }

  TEST_F(ErgunPressureSwitchTest, ConstantPressureRepressurisationKozenyCarman) {
    Run(MomentumModel::KarmenCozeny);
  }

} // namespace adsim
