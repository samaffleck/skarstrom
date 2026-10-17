/*
 * Skarstrom - pressure and temperature swing adsorption simulator
 * Copyright (c) 2024-2026 Sam Affleck and contributors
 *
 * Licensed under the terms in the LICENSE file at the repository root.
 */
#include "../pch.h"
#include "adsim/objects/compressor.h"
#include "adsim/math/constants.h"

// STL includes
#include <iomanip>
#include <iterator>
#include <sstream>


namespace adsim {

  double ConstantPumpPowerModel::CalculatePower(const PumpPowerModelParams &params) const {
    return power_;
  }

  double PolyPumpPowerModel::CalculatePower(const PumpPowerModelParams &params) const {
    return poly_.Evaluate(params.pressure);
  }

  double IsothermalCompressorEnergyModel::CalculateOutletTemperature(const EnergyModelParams& params) const {
    return params.T1;
  }

  double AdiabaticCompressorEnergyModel::CalculateOutletTemperature(const EnergyModelParams& params) const {
    const double P1 = params.P1;
    const double n = params.n;
    return params.T1 * std::pow((params.P2 / P1), ((n - 1) / n));
  }

  double ConstantPumpEfficiencyModel::CalculateEfficiency(double pressure) const {
    return efficiency_;
  }

  double PolyPumpEfficiencyModel::CalculateEfficiency(double pressure) const {
    double P2 = pressure * pressure;
    double P3 = P2 * pressure;
    return a_ + b_ * pressure + c_ * P3 + d_ * P3 * pressure + e_ * P3 * P2;
  }

  double VacuumPumpEfficiencyModel::CalculateEfficiency(double pressure) const {
    double p_ratio = pressure / 1e5;
    return (eta_max_ * alpha_ * p_ratio) / (1.0 + alpha_ * p_ratio);
  }

  namespace {

    // An efficiency (stored as a fraction) is valid when 0 < eta <= 1, i.e.
    // > 0 and <= 100 %. False for NaN.
    bool IsValidEfficiency(double eta) {
      return eta > 0.0 && eta <= 1.0;
    }

    std::string FormatNumber(double x) {
      std::ostringstream os;
      os << std::setprecision(6) << x;
      return os.str();
    }

    const char* EfficiencyModelName(PumpEfficiencyType type) {
      const auto i = static_cast<std::size_t>(type);
      return i < std::size(PumpEfficiencyTypeList) ? PumpEfficiencyTypeList[i] : "Unknown";
    }

  } // namespace

  bool IPumpEfficiencyModel::PerformChecks(double reference_pressure, std::string& error_msg) const {
    const double eta = CalculateEfficiency(reference_pressure);
    if (IsValidEfficiency(eta))
      return true;

    error_msg += "Pump efficiency must be > 0 and <= 100 % (entered as a fraction, 0 < efficiency <= 1), "
      "but the " + std::string(EfficiencyModelName(GetType())) + " efficiency model gives " +
      FormatNumber(eta) + " (" + FormatNumber(100.0 * eta) + " %)";
    if (GetType() != PumpEfficiencyType::CONSTANT)
      error_msg += " at a suction pressure of " + FormatNumber(reference_pressure) + " Pa";
    error_msg += ". ";
    return false;
  }

  bool VacuumPumpEfficiencyModel::PerformChecks(double /*reference_pressure*/, std::string& error_msg) const {
    const double eta_max = eta_max_;
    const double alpha = alpha_;
    if (IsValidEfficiency(eta_max) && alpha > 0.0 && std::isfinite(alpha))
      return true;

    error_msg += "Pump efficiency must be > 0 and <= 100 %: the Vacuum Pump efficiency model needs a "
      "maximum efficiency 0 < eta_max <= 1 (got " + FormatNumber(eta_max) + ", i.e. " +
      FormatNumber(100.0 * eta_max) + " %) and alpha > 0 (got " + FormatNumber(alpha) + "). ";
    return false;
  }

  IPump::IPump() {
    energy_eq_ = CreateCompressorEnergyModel(CompressorEnergyType::ISOTHERMAL);
    efficiency_model_ = CreatePumpEfficiencyModel(PumpEfficiencyType::CONSTANT);
  }

  bool IPump::PerformChecks(std::string& error_msg) const {
    if (!use_pump_) return true;

    const double n = poly_index_;
    if (!(n > 1.0) || !std::isfinite(n)) {
      error_msg += "Pump polytropic index must be greater than 1 (got " + std::to_string(n) +
        "); the polytropic work n/(n-1) is undefined at n = 1. Use a value slightly above 1 "
        "(e.g. 1.001) for near-isothermal compression. ";
      return false;
    }

    return CheckEfficiency(error_msg);
  }

  bool IPump::CheckEfficiency(std::string& error_msg) const {
    if (!efficiency_model_) {
      error_msg += "Pump efficiency model is null. ";
      return false;
    }
    // Blocking for both bounds: eta <= 0 gives an infinite/negative power, and
    // eta > 100 % is unphysical (the engine has no non-blocking warning channel).
    return efficiency_model_->PerformChecks(EfficiencyReferencePressure(), error_msg);
  }

  namespace {

    // Polytropic "head factor" x = (P2/P1)^((n-1)/n) - 1 of the work equation
    //   W = (1/eta) * (n/(n-1)) * P1 * Q1 * x.
    double PolytropicHeadFactor(double n, double P1, double P2) {
      return std::pow(P2 / P1, (n - 1) / n) - 1;
    }

    // Smallest pressure rise (P2/P1 - 1) above which the flow is the exact
    // inverse of the power. As x -> 0 a fixed power would drive an unbounded
    // flow (infinite at P2 == P1, e.g. at t = 0 of a pressurisation step that
    // starts at the suction pressure).
    constexpr double kMinPressureRatioRise = 1e-3;

    // C1-smooth floor on the head factor, used ONLY to invert power -> flow:
    //   g(x) = x                       for x >= x0
    //   g(x) = x0/2 + x^2 / (2 x0)     for 0 < x < x0
    //   g(x) = x0/2                    for x <= 0 (P2 <= P1)
    // where x0 is the head factor at P2/P1 = 1 + kMinPressureRatioRise. g and
    // dg/dx are continuous (IDA's Newton iteration needs a smooth residual)
    // and g >= x0/2 > 0, so the flow is finite everywhere, non-increasing in
    // P2/P1, and capped at its value for g = x0/2. When the discharge is at or
    // below the suction pressure (gas would pass freely) the pump delivers
    // that capped flow instead of an infinite flow (P2 == P1) or, as before,
    // no flow at all (P2 < P1, where the head factor turns negative).
    double RegularisedHeadFactor(double x, double n) {
      const double x0 = std::pow(1.0 + kMinPressureRatioRise, (n - 1) / n) - 1;
      if (x >= x0) return x;
      if (x <= 0.0) return 0.5 * x0;
      return 0.5 * x0 + 0.5 * x * x / x0;
    }

    // Floor on the efficiency used in the power/flow equations. An efficiency
    // <= 0 (or NaN) is rejected by IPump::PerformChecks; this only keeps the
    // maths finite and positive if one slips through (e.g. a pressure-
    // dependent model going non-positive at a run-time pressure the checks
    // could not sample). Far below any real efficiency, so valid models are
    // unaffected (the Vacuum Pump model only reaches it below ~0.01 Pa).
    constexpr double kMinPumpEfficiency = 1e-6;

    double SafeEfficiency(double eff) {
      return eff > kMinPumpEfficiency ? eff : kMinPumpEfficiency;
    }

  } // namespace

  inline double CompressorPower(double eff, double n, double P1, double P2, double Q1) {
    Q1 = std::abs(Q1);
    double power = (1 / SafeEfficiency(eff)) * (n / (n - 1)) * P1 * Q1 * PolytropicHeadFactor(n, P1, P2);
    return power;
  }

  // Inverse of CompressorPower for the suction volumetric flow Q1: exact for
  // P2/P1 >= 1 + kMinPressureRatioRise, regularised below (see above).
  inline double CompressorVolumetricFlowrate(double eff, double n, double P1, double P2, double power) {
    power = std::max(power, 0.0);
    const double x = RegularisedHeadFactor(PolytropicHeadFactor(n, P1, P2), n);
    return SafeEfficiency(eff) * ((n - 1) / n) * power / (P1 * x);
  }

  // The efficiency is always evaluated at the pump's SUCTION pressure P1 (the
  // suction or connected stream for a compressor, the bed for a vacuum pump):
  // the efficiency correlations are defined against the inlet pressure, and
  // using the same pressure in CalculatePower and CalculateMolarFlow keeps the
  // flow the exact inverse of the power for every efficiency model.
  double Compressor::CalculatePower(double P2, double mol_flow, double T1) const {
    if (!use_pump_) return 0.0;

    const double P1 = suction_pressure_;
    const double eff = efficiency_model_->CalculateEfficiency(P1);
    const double n = poly_index_;
    const double R = constants::kIdealGasConstant;
    const double Q1 = R * T1 * mol_flow / P1;

    double power = CompressorPower(eff, n, P1, P2, Q1);
    return std::max(power, 0.0); // No turbines
  }

  double Compressor::CalculateMolarFlow(double P2, double power, double T1) const {
    return CalculateMolarFlow(P2, power, T1, suction_pressure_);
  }

  double Compressor::CalculateMolarFlow(double P2, double power, double T1, double P1) const {
    if (!use_pump_) return 0.0;

    const double eff = efficiency_model_->CalculateEfficiency(P1);
    const double n = poly_index_;
    const double R = constants::kIdealGasConstant;
    const double Q1 = CompressorVolumetricFlowrate(eff, n, P1, P2, power);
    const double mol_flow = Q1 * P1 / (R * T1);
    return std::max(mol_flow, 0.0);
  }

  double Compressor::CalculateOutletTemperature(double outlet_pressure, double inlet_temperature) const {
    return CalculateOutletTemperature(outlet_pressure, inlet_temperature, suction_pressure_);
  }

  double Compressor::CalculateOutletTemperature(double outlet_pressure, double inlet_temperature, double suction_pressure) const {
    assert(energy_eq_);
    EnergyModelParams params;
    params.n = poly_index_;
    params.P1 = suction_pressure;
    params.P2 = outlet_pressure;
    params.T1 = inlet_temperature;
    return energy_eq_->CalculateOutletTemperature(params);
  }

  double VacuumPump::CalculatePower(double P1, double mol_flow, double T1) const {
    if (!use_pump_) return 0.0;

    const double P2 = discharge_pressure_;
    const double eff = efficiency_model_->CalculateEfficiency(P1);
    const double n = poly_index_;
    const double R = constants::kIdealGasConstant;
    const double Q1 = R * T1 * mol_flow / P1;

    double power = CompressorPower(eff, n, P1, P2, Q1);
    return std::max(power, 0.0); // No turbines
  }

  double VacuumPump::CalculateMolarFlow(double P1, double power, double T1) const {
    return CalculateMolarFlow(P1, power, T1, discharge_pressure_);
  }

  double VacuumPump::CalculateMolarFlow(double P1, double power, double T1, double P2) const {
    if (!use_pump_) return 0.0;

    const double eff = efficiency_model_->CalculateEfficiency(P1);
    const double n = poly_index_;
    const double R = constants::kIdealGasConstant;
    const double Q1 = CompressorVolumetricFlowrate(eff, n, P1, P2, power);
    const double mol_flow = Q1 * P1 / (R * T1);
    return std::max(mol_flow, 0.0);
  }

  double VacuumPump::CalculateOutletTemperature(double inlet_pressure, double inlet_temperature) const {
    assert(energy_eq_);
    EnergyModelParams params;
    params.n = poly_index_;
    params.P1 = inlet_pressure;
    params.P2 = discharge_pressure_;
    params.T1 = inlet_temperature;
    return energy_eq_->CalculateOutletTemperature(params);  
  }

} // end adsim namespace

CEREAL_REGISTER_TYPE(adsim::ConstantPumpPowerModel)
CEREAL_REGISTER_TYPE(adsim::PolyPumpPowerModel)
CEREAL_REGISTER_POLYMORPHIC_RELATION(adsim::IPumpPowerModel, adsim::ConstantPumpPowerModel)
CEREAL_REGISTER_POLYMORPHIC_RELATION(adsim::IPumpPowerModel, adsim::PolyPumpPowerModel)
CEREAL_CLASS_VERSION(adsim::ConstantPumpPowerModel, 2);
CEREAL_CLASS_VERSION(adsim::PolyPumpPowerModel, 2);

CEREAL_REGISTER_TYPE(adsim::IsothermalCompressorEnergyModel)
CEREAL_REGISTER_TYPE(adsim::AdiabaticCompressorEnergyModel)
CEREAL_REGISTER_POLYMORPHIC_RELATION(adsim::ICompressorEnergyModel, adsim::IsothermalCompressorEnergyModel)
CEREAL_REGISTER_POLYMORPHIC_RELATION(adsim::ICompressorEnergyModel, adsim::AdiabaticCompressorEnergyModel)
CEREAL_CLASS_VERSION(adsim::IsothermalCompressorEnergyModel, 1);
CEREAL_CLASS_VERSION(adsim::AdiabaticCompressorEnergyModel, 1);

CEREAL_REGISTER_TYPE(adsim::ConstantPumpEfficiencyModel)
CEREAL_REGISTER_TYPE(adsim::VacuumPumpEfficiencyModel)
CEREAL_REGISTER_TYPE(adsim::PolyPumpEfficiencyModel)
CEREAL_REGISTER_POLYMORPHIC_RELATION(adsim::IPumpEfficiencyModel, adsim::ConstantPumpEfficiencyModel)
CEREAL_REGISTER_POLYMORPHIC_RELATION(adsim::IPumpEfficiencyModel, adsim::VacuumPumpEfficiencyModel)
CEREAL_REGISTER_POLYMORPHIC_RELATION(adsim::IPumpEfficiencyModel, adsim::PolyPumpEfficiencyModel)
CEREAL_CLASS_VERSION(adsim::ConstantPumpEfficiencyModel, 2);
CEREAL_CLASS_VERSION(adsim::VacuumPumpEfficiencyModel, 2);
CEREAL_CLASS_VERSION(adsim::PolyPumpEfficiencyModel, 2);

CEREAL_REGISTER_TYPE(adsim::Compressor)
CEREAL_REGISTER_TYPE(adsim::VacuumPump)
CEREAL_REGISTER_POLYMORPHIC_RELATION(adsim::IPump, adsim::Compressor)
CEREAL_REGISTER_POLYMORPHIC_RELATION(adsim::IPump, adsim::VacuumPump)
CEREAL_CLASS_VERSION(adsim::IPump, 3);
CEREAL_CLASS_VERSION(adsim::Compressor, 2);
CEREAL_CLASS_VERSION(adsim::VacuumPump, 2);
