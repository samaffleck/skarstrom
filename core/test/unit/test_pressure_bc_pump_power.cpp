/*
 * Skarstrom - pressure and temperature swing adsorption simulator
 * Copyright (c) 2024-2026 Sam Affleck and contributors
 *
 * Licensed under the terms in the LICENSE file at the repository root.
 */

// Regression tests for the compressor power logged by PressureBoundaryCondition.
//
// Compressor::CalculatePower(P2, n_dot, T1) evaluates the suction volumetric
// flow Q1 = R T1 n_dot / P1, so T1 must be the pre-compression (suction)
// temperature. With an ADIABATIC energy model the pressure BC also logs the
// post-compression temperature T2 = T1 (P2/P1)^((n-1)/n) for the stream that
// enters the bed; the power must still be evaluated at T1, not T2.

// STL includes
#include <cmath>
#include <memory>
#include <string>

// adsim includes
#include "adsim/boundary_conditions/boundary_conditions.h"
#include "adsim/boundary_conditions/boundary_stream.h"
#include "adsim/equations/equation.h"
#include "adsim/objects/compressor.h"
#include "adsim/solver/log.h"

// Gtest includes
#include <gtest/gtest.h>


namespace adsim {

  namespace {

    constexpr int kCells = 6; // includes the two ghost cells

    // Minimal state for PressureBoundaryCondition::Log: uniform pressure,
    // temperature and molar flow, and a single-component (N2) gas.
    Equations MakeUniformState(double pressure, double temperature, double molar_flow) {
      Equations state;
      state.pressure.var.assign(kCells, pressure);
      state.fluid_temperature.var.assign(kCells, temperature);
      state.molar_flow_rate.var.assign(kCells + 1, molar_flow);
      Equation y{ "Mole Fraction", CellType::CellCentre, EquationType::Algebraic, nullptr };
      y.var.assign(kCells, 1.0);
      state.mole_fraction.emplace("N2", std::move(y));
      return state;
    }

    void ConfigurePump(IPump& pump, double polytropic_index, double efficiency) {
      pump.use_pump_ = true;
      pump.poly_index_ = polytropic_index;
      pump.energy_eq_ = CreateCompressorEnergyModel(CompressorEnergyType::ADIABATIC);
      auto eff = std::make_unique<ConstantPumpEfficiencyModel>();
      eff->efficiency_ = efficiency;
      pump.efficiency_model_ = std::move(eff);
    }

  } // namespace

  // Pressure-inlet BC with a compressor, adiabatic energy model, P2/P1 = 5.
  // The suction side is the specified inlet stream (temperature_), so the
  // logged power / COMPRESSOR_POWER KPI must equal CalculatePower at temperature_.
  TEST(PressureBCPumpPower, InletCompressorPowerUsesSuctionTemperature) {
    const double P1 = 1.0e5;      // compressor suction pressure [Pa]
    const double P2 = 5.0e5;      // bed boundary (discharge) pressure [Pa]
    const double T_suction = 300.0;
    const double T_bed = 350.0;   // distinct from both T1 and T2
    const double molar_flow = 0.1;
    const double n = 1.4;

    for (BoundaryLocation loc : { BoundaryLocation::BOTTOM, BoundaryLocation::TOP }) {
      PressureBoundaryCondition bc;
      bc.location_ = loc;
      bc.direction_ = BoundaryDirection::INLET;
      bc.temperature_ = T_suction;
      bc.mole_fraction_ = { { "N2", 1.0 } };
      bc.pump_ = CreatePump(PumpType::COMPRESSOR);
      ConfigurePump(*bc.pump_, n, 0.72);
      static_cast<Compressor&>(*bc.pump_).suction_pressure_ = P1;

      const Equations state = MakeUniformState(P2, T_bed, molar_flow);
      bc.Log(state, 0.0);

      const double expected_power = bc.pump_->CalculatePower(P2, molar_flow, T_suction);
      const double expected_T2 = T_suction * std::pow(P2 / P1, (n - 1.0) / n);
      ASSERT_GT(expected_power, 0.0);

      const auto& sd = bc.stream_data_;
      ASSERT_EQ(sd.electrical_power_.size(), 1u);
      EXPECT_NEAR(sd.electrical_power_[0], expected_power, 1e-9 * expected_power);

      const auto kpi = GetBoundaryVariableVector(BoundaryVariableType::COMPRESSOR_POWER, sd);
      ASSERT_EQ(kpi.size(), 1u);
      EXPECT_NEAR(kpi[0], expected_power, 1e-9 * expected_power);

      // The logged stream temperature is still the post-compression temperature
      // of the gas entering the bed (what Apply() imposes).
      ASSERT_EQ(sd.temperature_boundary_.size(), 1u);
      EXPECT_NEAR(sd.temperature_boundary_[0], expected_T2, 1e-9 * expected_T2);
    }
  }

  // Pressure-outlet BC with a vacuum pump: the suction side is the bed, so the
  // power is evaluated at the boundary-adjacent bed temperature.
  TEST(PressureBCPumpPower, OutletVacuumPumpPowerUsesBedTemperature) {
    const double P_bed = 2.0e4;   // vacuum pump suction pressure [Pa]
    const double P_discharge = 1.0e5;
    const double T_bed = 320.0;
    const double molar_flow = 0.05;
    const double n = 1.4;

    for (BoundaryLocation loc : { BoundaryLocation::BOTTOM, BoundaryLocation::TOP }) {
      PressureBoundaryCondition bc;
      bc.location_ = loc;
      bc.direction_ = BoundaryDirection::OUTLET;
      bc.temperature_ = 290.0; // inlet temperature: must not be used for an outlet
      bc.pump_ = CreatePump(PumpType::VACUUM_PUMP);
      ConfigurePump(*bc.pump_, n, 0.6);
      static_cast<VacuumPump&>(*bc.pump_).discharge_pressure_ = P_discharge;

      const Equations state = MakeUniformState(P_bed, T_bed, molar_flow);
      bc.Log(state, 0.0);

      const double expected_power = bc.pump_->CalculatePower(P_bed, molar_flow, T_bed);
      ASSERT_GT(expected_power, 0.0);

      const auto& sd = bc.stream_data_;
      ASSERT_EQ(sd.electrical_power_.size(), 1u);
      EXPECT_NEAR(sd.electrical_power_[0], expected_power, 1e-9 * expected_power);
      ASSERT_EQ(sd.temperature_boundary_.size(), 1u);
      EXPECT_NEAR(sd.temperature_boundary_[0], T_bed, 1e-9 * T_bed);
    }
  }

} // end adsim namespace
