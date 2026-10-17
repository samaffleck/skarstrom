/*
 * Skarstrom - pressure and temperature swing adsorption simulator
 * Copyright (c) 2024-2026 Sam Affleck and contributors
 *
 * Licensed under the terms in the LICENSE file at the repository root.
 */
#include "../pch.h"
#include "adsim/boundary_conditions/boundary_stream.h"
#include "adsim/boundary_conditions/boundary_value.h"
#include "adsim/boundary_conditions/step.h"
#include "adsim/math/constants.h"
#include "adsim/math/general_functions.h"
#include "adsim/equations/equation.h"


namespace adsim {

  void IBoundaryStream::Log(const Equations& state, double time, BoundaryLocation location) {
    auto& component_mole_fractions = state.mole_fraction;
    auto& T = state.fluid_temperature.var;
    auto& P = state.pressure.var;

    int i = 0;
    int d = 1;
    int N = static_cast<int>(T.size()) - 1;
    if (location == BoundaryLocation::TOP) {
      i = N;
      d = -1;
    }

    const double boundary_temperature = 0.5 * (T[i] + T[i + d]);
    const double boundary_pressure = 0.5 * (P[i] + P[i + d]);

    double molar_flow = GetMolarFlow(state, time, location);
    double enthalpy_flow = GetEnthalpyFlow(state, time, location);
    constexpr double kTinyFlow = 1e-12;

    // Log the stream-side temperature. For inlets this should reflect the
    // connected/specified stream value (including pump heating via enthalpy flow).
    double temperature = boundary_temperature;
    if (std::abs(molar_flow) > kTinyFlow) {
      temperature = enthalpy_flow / molar_flow;
    } else if (direction_ == BoundaryDirection::INLET) {
      temperature = InletTemperature(time);
    }

    std::unordered_map<std::string, double> component_molar_flow{};
    std::unordered_map<std::string, double> logged_mole_fraction{};
    for (const auto& [component, yi] : component_mole_fractions) {
      const double n_i = GetComponentMolarFlow(state, time, location, component);
      component_molar_flow[component] = n_i;

      double y = 0.0;
      if (std::abs(molar_flow) > kTinyFlow) {
        y = n_i / molar_flow;
      } else if (direction_ == BoundaryDirection::INLET) {
        y = InletMoleFraction(time, component);
      } else {
        y = 0.5 * (yi.var[i] + yi.var[i + d]);
      }

      if (!std::isfinite(y)) y = 0.0;
      if (y < 0.0 && std::abs(y) < 1e-12) y = 0.0;
      if (y < 0.0) y = 0.0;
      logged_mole_fraction[component] = y;
    }

    double sum_y = 0.0;
    for (const auto& [component, y] : logged_mole_fraction) {
      (void)component;
      sum_y += y;
    }
    if (sum_y > 0.0) {
      for (auto& [component, y] : logged_mole_fraction) {
        (void)component;
        y /= sum_y;
      }
    }

    double power = GetPower(boundary_pressure, molar_flow, time);

    stream_data_.Log(
      molar_flow,
      enthalpy_flow,
      temperature,
      boundary_pressure,
      power,
      time,
      component_molar_flow,
      logged_mole_fraction
    );
  }

  bool IBoundaryStream::PerformChecks(std::string& error_msg, const Fluid& fluid) {
    bool connection_checks_passed = connection_.PerformChecks(error_msg);
    if (!connection_checks_passed)
      return false;

    if (!buffer_tank_.PerformChecks(error_msg))
      return false;

    bool specific_check_passed = PerformSpecificChecks(error_msg, fluid);
    if (!specific_check_passed)
      return false;

    // Key the inlet composition exactly by the fluid's components (back-fill
    // missing ones with 0, drop removed ones) so the residual never meets a
    // missing key. Unconditional: a component swapped for another leaves the
    // count unchanged but the keys wrong.
    ReconcileMoleFractionMap(mole_fraction, fluid.components_);

    if (direction_ == BoundaryDirection::OUTLET)
      return true;

    // Is an INLET
    double totalMoleFractions = 0;
    for (const auto& [component, yi] : mole_fraction) {
      (void)component;
      totalMoleFractions += yi;
    }

    if (std::abs(1 - totalMoleFractions) > 1e-10) {
      error_msg += "Sum of mole fractions in the Flow step does not equal 1. ";
      return false;
    }

    return true;  
  }

  ValveBoundaryStream::ValveBoundaryStream() {
    boundary_value_ = std::make_unique<ConstantBoundaryValue>(101325.0);
  }

  // void ValveBoundaryStream::Initialise(double step_time) {}

  void ValveBoundaryStream::Initialise(Equations& state, double step_time, BoundaryLocation location, BoundaryDirection dir) {
    int i = 0;
    int d = 1;
    auto& P = state.pressure.var;
    int N = P.size() - 1;
    if (location == BoundaryLocation::TOP) {
      i = N; d = -1;
    }

    double Pstart = 0.5 * (P[i] + P[i + d]);
    boundary_value_->InitialiseStartValue(Pstart);
    boundary_value_->InitialiseStepTime(step_time);
    direction_ = dir;
    step_time_cached_ = step_time;
  }

  double ValveBoundaryStream::GetMolarFlow(const Equations& state, double time, BoundaryLocation loc) {
    const auto& T = state.fluid_temperature.var;
    const auto& P = state.pressure.var;
    const auto& C = state.concentration.var;
    const auto& Cp = state.heat_capacity.var;
    const auto& Mr = state.molar_mass.var;
    const double R = constants::kIdealGasConstant;

    int N = (int)P.size() - 1;
    int i = 0;
    int d = 1;
    if (loc == BoundaryLocation::TOP) {
      i = N;
      d = -1;
    }

    // double C_boundary = 0.5 * (C[i] + C[i + d]);
    double T_boundary = 0.5 * (T[i] + T[i + d]);
    double P_boundary = 0.5 * (P[i] + P[i + d]);
    // double P_boundary = C_boundary * R * T_boundary;
    double Cp_boundary = 0.5 * (Cp[i] + Cp[i + d]);
    double Mr_boundary = 0.5 * (Mr[i] + Mr[i + d]);

    auto get_molar_flow = [&](){
      if (direction_ == BoundaryDirection::INLET) {
        // Upstream gas is the buffer tank's once it is ready (the tank's
        // pressure is still the user-specified upstream boundary value).
        double mol_flow = valve_.GetMolarFlowRate(boundary_value_->GetValue(time), P_boundary, Mr_boundary, InletTemperature(time), Cp_boundary);
        if (loc == BoundaryLocation::TOP)
          return -mol_flow;
        else
          return mol_flow;
      } else { // Outlet
        double mol_flow = valve_.GetMolarFlowRate(P_boundary, boundary_value_->GetValue(time), Mr_boundary, T_boundary, Cp_boundary);
        if (loc == BoundaryLocation::TOP)
          return mol_flow;
        else
          return -mol_flow;
      }
    };

    if (!connection_.IsConnected()) {
      return get_molar_flow();
    } else { // Connected, get flow rate from another stream
      assert(connection_.stream_ != nullptr);

      if (connection_.stream_->molar_flow_.size() <= 0) { 
        return get_molar_flow(); // If no data, use the user specified assumption for first cycle
      }

      if (direction_ == BoundaryDirection::INLET) {
        double mol_flow = valve_.GetMolarFlowRate(connection_.stream_->GetPressure(time), P_boundary, Mr_boundary, AssumedInletTemperature(time), Cp_boundary);
        if (loc == BoundaryLocation::TOP)
          return -mol_flow;
        else
          return mol_flow;
      } else { // Outlet
        double mol_flow = valve_.GetMolarFlowRate(P_boundary, connection_.stream_->GetPressure(time), Mr_boundary, T_boundary, Cp_boundary);
        if (loc == BoundaryLocation::TOP)
          return mol_flow;
        else
          return -mol_flow;
      }
    }
  }

  double ValveBoundaryStream::GetComponentMolarFlow(const Equations& state, double time, BoundaryLocation loc, const std::string& component) {
    double n = GetMolarFlow(state, time, loc);
    auto& yi = state.mole_fraction.at(component).var;

    double inlet_mole_fraction = 0.0;
    if (direction_ == BoundaryDirection::INLET) {
      inlet_mole_fraction = InletMoleFraction(time, component);
    }

    int N = (int)yi.size() - 1;
    if (direction_ == BoundaryDirection::INLET)
      return n * inlet_mole_fraction;
    else { // OUTLET
      if (loc == BoundaryLocation::TOP)
        return n * 0.5 * (yi[N] + yi[N - 1]); 
      else // BOTTOM
        return n * 0.5 * (yi[0] + yi[1]);
    }
  }

  double ValveBoundaryStream::GetEnthalpyFlow(const Equations& state, double time, BoundaryLocation loc) {
    double n = GetMolarFlow(state, time, loc);
    auto& T = state.fluid_temperature.var;

    double inlet_temperature = 0.0;
    if (direction_ == BoundaryDirection::INLET) {
      inlet_temperature = InletTemperature(time);
    }

    int N = (int)T.size() - 1;
    if (direction_ == BoundaryDirection::INLET)
      return n * inlet_temperature;
    else { // OUTLET
      if (loc == BoundaryLocation::TOP)
        return n * 0.5 * (T[N] + T[N - 1]); 
      else // BOTTOM
        return n * 0.5 * (T[0] + T[1]);
    }
  }

  bool ValveBoundaryStream::PerformSpecificChecks(std::string& error_msg, const Fluid& fluid) {
    return true;
  }

  BoundaryStream::BoundaryStream() {
    boundary_flow_ = std::make_unique<VelocityBoundaryValue>();
    pump_ = CreatePump(PumpType::COMPRESSOR);
  }

  double BoundaryStream::GetMolarFlow(const Equations& state, double time, BoundaryLocation loc) {
    double mol_flow = 0.0;

    if (buffer_tank_.IsActive() && direction_ == BoundaryDirection::INLET) {
      // Constant inlet flow rate: source's last-cycle total moles, scaled by
      // pct and spread over the consumer step's duration. Falls back to the
      // user-specified first-cycle assumption until the source has run.
      if (buffer_tank_.IsReady()) {
        mol_flow = buffer_tank_.GetMolarFlow(step_time_cached_);
      } else {
        mol_flow = boundary_flow_->GetMolarFlow(state, time, loc);
      }
    } else if (!connection_.IsConnected()) {
      mol_flow = boundary_flow_->GetMolarFlow(state, time, loc);
    } else {
      assert(connection_.stream_ != nullptr);
      if (connection_.stream_->molar_flow_.size() > 0) {
        mol_flow = connection_.stream_->GetMolarFlow(time);
      } else {
        mol_flow = boundary_flow_->GetMolarFlow(state, time, loc); // First cycle uses the specified values
      }
    }

    mol_flow = std::abs(mol_flow); // We give the right sign with the logic below:
    if (
      direction_ == BoundaryDirection::OUTLET && loc == BoundaryLocation::TOP ||
      direction_ == BoundaryDirection::INLET && loc == BoundaryLocation::BOTTOM
    ) {
      return mol_flow;
    }
    else
      return -mol_flow;
  }

  double BoundaryStream::GetComponentMolarFlow(const Equations& state, double time, BoundaryLocation loc, const std::string& component) {
    double n = GetMolarFlow(state, time, loc);
    auto& yi = state.mole_fraction.at(component).var;

    double inlet_mole_fraction = 0.0;
    if (direction_ == BoundaryDirection::INLET) {
      inlet_mole_fraction = InletMoleFraction(time, component);
    }

    int N = (int)yi.size() - 1;

    if (direction_ == BoundaryDirection::INLET)
      return n * inlet_mole_fraction;
    else { // OUTLET
      if (loc == BoundaryLocation::TOP)
        return n * 0.5 * (yi[N] + yi[N - 1]);
      else // BOTTOM
        return n * 0.5 * (yi[0] + yi[1]);
    }
  }

  double BoundaryStream::GetEnthalpyFlow(const Equations& state, double time, BoundaryLocation loc) {
    double n = GetMolarFlow(state, time, loc);
    auto& T = state.fluid_temperature.var;
    auto& P = state.pressure.var;

    double inlet_temperature = 0.0;
    double outlet_pressure = 0.0;
    const int N = (int)P.size() - 1;

    if (loc == BoundaryLocation::TOP) {
      outlet_pressure = 0.5 * (P[N] + P[N - 1]);
    } else {
      outlet_pressure = 0.5 * (P[0] + P[1]);
    }

    if (direction_ == BoundaryDirection::INLET) {
      inlet_temperature = InletTemperature(time);

      if (pump_->use_pump_ && !isothermal_temperature_) { // no compression heating in an isothermal run
        // Increases the inlet temperature due to compression
        double outlet_temperature = pump_->CalculateOutletTemperature(outlet_pressure, inlet_temperature);
        inlet_temperature = outlet_temperature; // The outlet temp from the compressor is used as the inlet temp to the bed.
      }
    }

    if (direction_ == BoundaryDirection::INLET)
      return n * inlet_temperature;
    else { // OUTLET
      if (loc == BoundaryLocation::TOP)
        return n * 0.5 * (T[N] + T[N - 1]);
      else // BOTTOM
        return n * 0.5 * (T[0] + T[1]);
    }
  }

  double BoundaryStream::GetPower(double pressure, double molar_flow, double time) {
    // Suction side is the user-specified inlet gas (as before; isothermal runs
    // use the isothermal temperature, see AssumedInletTemperature).
    return pump_->CalculatePower(pressure, molar_flow, AssumedInletTemperature(time));
  }

  bool BoundaryStream::PerformSpecificChecks(std::string& error_msg, const Fluid& fluid) {
    // The optional pump only computes the logged power here, but an invalid
    // polytropic index or efficiency would still make that power non-finite.
    if (pump_ && !pump_->PerformChecks(error_msg)) {
      error_msg += "(pump of stream '" + name_ + "') ";
      return false;
    }
    return true;
  }

  // ========================================
  // Pump boundary stream
  // ========================================
  PumpBoundaryStream::PumpBoundaryStream() {
    pump_ = CreatePump(PumpType::COMPRESSOR);
    pump_->use_pump_ = true;
    power_ = CreatePumpPowerModel(PumpPowerType::CONSTANT);
  }

  double PumpBoundaryStream::GetMolarFlow(const Equations& state, double time, BoundaryLocation loc) {    
    const auto& T = state.fluid_temperature.var;
    const auto& P = state.pressure.var;
    const double R = constants::kIdealGasConstant;

    int N = (int)P.size() - 1;

    double P_boundary = 0.0;
    double T_boundary = 0.0;

    if (loc == BoundaryLocation::TOP) {
      P_boundary = 0.5 * (P[N] + P[N - 1]);
      T_boundary = 0.5 * (T[N] + T[N - 1]);
    } else {
      P_boundary = 0.5 * (P[0] + P[1]);
      T_boundary = 0.5 * (T[0] + T[1]);
    }

    auto get_molar_flow = [&](){
      if (direction_ == BoundaryDirection::INLET) { // Compressor
        // Suction-side temperature: the buffer tank's once it is ready, the
        // same value GetEnthalpyFlow compresses from.
        double mol_flow = pump_->CalculateMolarFlow(P_boundary, power_->CalculatePower(PumpPowerModelParams{P_boundary}), InletTemperature(time));
        if (loc == BoundaryLocation::TOP)
          return -mol_flow;
        else
          return mol_flow;
      } else { // Outlet -> Vacuum pump
        double mol_flow = pump_->CalculateMolarFlow(P_boundary, power_->CalculatePower(PumpPowerModelParams{P_boundary}), T_boundary);
        if (loc == BoundaryLocation::TOP)
          return mol_flow;
        else
          return -mol_flow;
      }
    };

    if (!connection_.IsConnected()) {
      return get_molar_flow();
    } else { // Connected, get flow rate from another stream
      assert(connection_.stream_ != nullptr);

      if (connection_.stream_->molar_flow_.size() <= 0) { 
        return get_molar_flow(); // If no data, use the user specified assumption for first cycle
      }

      // The connected stream's pressure replaces the user's suction/discharge
      // pressure for this evaluation only. It is passed as an argument rather
      // than assigned to the ExpressionParam, which would overwrite the user's
      // expression string (and the saved deck) on every residual evaluation.
      const double P_connected = connection_.stream_->GetPressure(time);

      if (direction_ == BoundaryDirection::INLET) {
        if (auto* comp = dynamic_cast<Compressor*>(pump_.get())){
          // Suction side of the compressor is the connected stream, so use its
          // temperature (the same one GetEnthalpyFlow uses for this case).
          const double T_connected = isothermal_temperature_ ? *isothermal_temperature_ : connection_.stream_->GetTemperature(time);
          double mol_flow = comp->CalculateMolarFlow(P_boundary, power_->CalculatePower(PumpPowerModelParams{P_boundary}), T_connected, P_connected);
          if (loc == BoundaryLocation::TOP)
            return -mol_flow;
          else
            return mol_flow;
        }
      } else { // Outlet
        if (auto* vp = dynamic_cast<VacuumPump*>(pump_.get())){
          // Suction side of the vacuum pump is the bed, so use the bed-side
          // temperature as the unconnected branch does.
          double mol_flow = vp->CalculateMolarFlow(P_boundary, power_->CalculatePower(PumpPowerModelParams{P_boundary}), T_boundary, P_connected);
          if (loc == BoundaryLocation::TOP)
            return mol_flow;
          else
            return -mol_flow;
        }
      }
      // Fallback if dynamic_cast fails
      return get_molar_flow();
    } 
  }

  bool PumpBoundaryStream::PerformSpecificChecks(std::string& error_msg, const Fluid& /*fluid*/) {
    if (!pump_) {
      error_msg += "Pump is null. ";
      return false;
    }

    if (!power_) {
      error_msg += "Pump power model is null. ";
      return false;
    }

    switch (direction_) {
    case BoundaryDirection::OUTLET:
      if (pump_->GetType() != PumpType::VACUUM_PUMP) {
        error_msg += "Direction = OUTLET but Pump Type is not Vacuum Pump. ";
        return false;
      }
      break;

    case BoundaryDirection::INLET:
      if (pump_->GetType() != PumpType::COMPRESSOR) {
        error_msg += "Direction = INLET but Pump Type is not Compressor. ";
        return false;
      }
      break;

    default:
      error_msg += "Unknown boundary direction. ";
      return false;
    }

    if (!pump_->PerformChecks(error_msg)) {
      error_msg += "(pump of stream '" + name_ + "') ";
      return false;
    }

    return true; // all checks passed
  }

  double PumpBoundaryStream::GetComponentMolarFlow(const Equations& state, double time, BoundaryLocation loc, const std::string& component) {
    double n = GetMolarFlow(state, time, loc);
    auto& yi = state.mole_fraction.at(component).var;

    double inlet_mole_fraction = 0.0;
    if (direction_ == BoundaryDirection::INLET) {
      inlet_mole_fraction = InletMoleFraction(time, component);
    }

    int N = (int)yi.size() - 1;

    if (direction_ == BoundaryDirection::INLET)
      return n * inlet_mole_fraction;
    else { // OUTLET
      if (loc == BoundaryLocation::TOP)
        return n * 0.5 * (yi[N] + yi[N - 1]);
      else // BOTTOM
        return n * 0.5 * (yi[0] + yi[1]);
    }
  }
  
  double PumpBoundaryStream::GetEnthalpyFlow(const Equations& state, double time, BoundaryLocation loc) {
    double n = GetMolarFlow(state, time, loc);
    auto& T = state.fluid_temperature.var;
    auto& P = state.pressure.var;

    double inlet_temperature = 0.0;
    double outlet_pressure = 0.0;
    const int N = (int)P.size() - 1;

    if (loc == BoundaryLocation::TOP) {
      outlet_pressure = 0.5 * (P[N] + P[N - 1]);
    } else {
      outlet_pressure = 0.5 * (P[0] + P[1]);
    }

    if (direction_ == BoundaryDirection::INLET) {
      inlet_temperature = InletTemperature(time);

      if (pump_->use_pump_ && !isothermal_temperature_) { // no compression heating in an isothermal run
        // Increases the inlet temperature due to compression
        double outlet_temperature = pump_->CalculateOutletTemperature(outlet_pressure, inlet_temperature);
        // When connected (with logged data) the suction pressure is the
        // connected stream's pressure, consistent with GetMolarFlow.
        if (connection_.IsConnected() && connection_.stream_->molar_flow_.size() > 0) {
          if (auto* comp = dynamic_cast<Compressor*>(pump_.get())) {
            outlet_temperature = comp->CalculateOutletTemperature(outlet_pressure, inlet_temperature, connection_.stream_->GetPressure(time));
          }
        }
        inlet_temperature = outlet_temperature; // The outlet temp from the compressor is used as the inlet temp to the bed.
      }
    }

    if (direction_ == BoundaryDirection::INLET)
      return n * inlet_temperature;
    else { // OUTLET
      if (loc == BoundaryLocation::TOP)
        return n * 0.5 * (T[N] + T[N - 1]);
      else // BOTTOM
        return n * 0.5 * (T[0] + T[1]);
    }
  }
  
  double PumpBoundaryStream::GetPower(double pressure, double molar_flow, double /*time*/) {
    return power_->CalculatePower(PumpPowerModelParams{pressure});
  }

  void StreamData::Clear() {
    molar_flow_.clear();
    enthalpy_flow_.clear();
    temperature_boundary_.clear();
    pressure_boundary_.clear();
    electrical_power_.clear();
    time_series_.clear();
    component_molar_flow_.clear();
    mole_fraction_boundary_.clear();
  }

  void StreamData::Log(
    double molar_flow, 
    double enthalpy_flow,
    double temperature_boundary, 
    double pressure_boundary, 
    double power,
    double time,
    std::unordered_map<std::string, double> component_molar_flow,
    std::unordered_map<std::string, double> mole_fraction_boundary
  ) {
    molar_flow_.emplace_back(molar_flow);
    enthalpy_flow_.emplace_back(enthalpy_flow);
    temperature_boundary_.emplace_back(temperature_boundary);
    pressure_boundary_.emplace_back(pressure_boundary);
    electrical_power_.emplace_back(power);
    time_series_.emplace_back(time);

    for (const auto& [comp, mf] : component_molar_flow) {
      component_molar_flow_[comp].emplace_back(mf);
    }

    if (mole_fraction_boundary.empty()) return;  // begin() below would dereference end()
    double sum_indep = 0.0;
    const std::string last_comp = mole_fraction_boundary.begin()->first;
    for (const auto& [comp, x] : mole_fraction_boundary) {
      if (comp == last_comp) continue; // ignore if accidentally passed
      mole_fraction_boundary_[comp].emplace_back(x);
      sum_indep += x;
    }
    double x_dep = 1.0 - sum_indep;
    mole_fraction_boundary_[last_comp].emplace_back(x_dep);
  }

  // InterpolateTimeSeries is defined once in outlet_conditions.cpp.

  double StreamData::GetMolarFlow(double time) const {
    return InterpolateTimeSeries(time_series_, molar_flow_, time);
  }

  double StreamData::GetEnthalpyFlow(double time) const {
    return InterpolateTimeSeries(time_series_, enthalpy_flow_, time);
  }

  double StreamData::GetTemperature(double time) const {
    return InterpolateTimeSeries(time_series_, temperature_boundary_, time);
  }

  double StreamData::GetPressure(double time) const {
    return InterpolateTimeSeries(time_series_, pressure_boundary_, time);
  }

  double StreamData::GetComponentFlow(double time, const std::string& component) const {
    return InterpolateTimeSeries(time_series_, component_molar_flow_.at(component), time);
  }

  double StreamData::GetMoleFraction(double time, const std::string& component) const {
    // Read by connected inlets inside the IDA residual: a component the source
    // stream never logged contributes 0 rather than throwing.
    const auto it = mole_fraction_boundary_.find(component);
    if (it == mole_fraction_boundary_.end()) return 0.0;
    return InterpolateTimeSeries(time_series_, it->second, time);
  }

  double StreamData::GetAveragePurity(const std::string& component) const {
    double comp_molar_flow = std::abs(IntegrateVariableInTime(component_molar_flow_.at(component), time_series_));
    double total_molar_flow = std::abs(IntegrateVariableInTime(molar_flow_, time_series_));
    
    double purity = 0.0;
    if (total_molar_flow > 0) {
      purity = comp_molar_flow / total_molar_flow; 
    }

    return purity;
  }

} // end adsim namespace

CEREAL_REGISTER_TYPE(adsim::BoundaryStream)
CEREAL_REGISTER_TYPE(adsim::ValveBoundaryStream)
CEREAL_REGISTER_TYPE(adsim::PumpBoundaryStream)
CEREAL_REGISTER_POLYMORPHIC_RELATION(adsim::IBoundaryStream, adsim::BoundaryStream)
CEREAL_REGISTER_POLYMORPHIC_RELATION(adsim::IBoundaryStream, adsim::ValveBoundaryStream)
CEREAL_REGISTER_POLYMORPHIC_RELATION(adsim::IBoundaryStream, adsim::PumpBoundaryStream)
CEREAL_CLASS_VERSION(adsim::IBoundaryStream, 4);
CEREAL_CLASS_VERSION(adsim::BoundaryStream, 1);
CEREAL_CLASS_VERSION(adsim::ValveBoundaryStream, 1);
CEREAL_CLASS_VERSION(adsim::PumpBoundaryStream, 1);
