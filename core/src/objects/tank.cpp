/*
 * Skarstrom - pressure and temperature swing adsorption simulator
 * Copyright (c) 2024-2026 Sam Affleck and contributors
 *
 * Licensed under the terms in the LICENSE file at the repository root.
 */
#include "../pch.h"
#include "adsim/objects/tank.h"
#include "adsim/objects/fluid.h"
#include "adsim/math/constants.h"

#include <algorithm>
#include <cmath>

namespace adsim {

  namespace detail {
    double SmoothPositive(double x, double width) {
      if (!(width > 0.0) || x >= width) return x > 0.0 ? x : 0.0;
      if (x <= -width) return 0.0;
      // C1 quadratic blend on (-w, w): never negative, so a smoothed check
      // valve or regulator never passes flow backwards.
      return (x + width) * (x + width) / (4.0 * width);
    }
  }

  namespace {
    bool CheckComposition(ExpressionMap& y, const Fluid& fluid, const std::string& what,
                          std::string& error_msg) {
      ReconcileMoleFractionMap(y, fluid.components_);
      double sum = 0.0;
      for (const auto& [component, yi] : y) {
        if (yi < 0.0) {
          error_msg += what + ": negative mole fraction of " + component + ". ";
          return false;
        }
        sum += yi;
      }
      if (std::abs(sum - 1.0) > 1e-10) {
        error_msg += what + ": mole fractions sum to " + std::to_string(sum) + ", not 1. ";
        return false;
      }
      return true;
    }

    bool Positive(double v, const std::string& what, std::string& error_msg) {
      if (v > 0.0 && std::isfinite(v)) return true;
      error_msg += what + " must be positive (is " + std::to_string(v) + "). ";
      return false;
    }

    bool NonNegative(double v, const std::string& what, std::string& error_msg) {
      if (v >= 0.0 && std::isfinite(v)) return true;
      error_msg += what + " must not be negative (is " + std::to_string(v) + "). ";
      return false;
    }
  }

  bool ITankStream::PerformChecks(std::string& error_msg, const Fluid& fluid) {
    if (!CanFlowIn()) {
      // The composition is never used; keep its keys in step with the fluid.
      ReconcileMoleFractionMap(mole_fraction_, fluid.components_);
      return true;
    }
    if (!CheckComposition(mole_fraction_, fluid, "Tank stream '" + name_ + "'", error_msg)) return false;
    return Positive(temperature_, "Temperature of tank stream '" + name_ + "'", error_msg);
  }

  double TankMakeUp::FlowIn(const TankLiveState& s) const {
    if (mode_ == TankMakeUpMode::SET_FLOW) return molar_flow_;
    return cv_ * detail::SmoothPositive(header_pressure_ - s.pressure, kTankValveSmoothingPressure) / 1e5;
  }

  bool TankMakeUp::PerformChecks(std::string& error_msg, const Fluid& fluid) {
    if (!ITankStream::PerformChecks(error_msg, fluid)) return false;
    const std::string w = "Make-up '" + name_ + "'";
    if (mode_ == TankMakeUpMode::SET_FLOW) return NonNegative(molar_flow_, w + " molar flow", error_msg);
    return NonNegative(cv_, w + " Cv", error_msg) && Positive(header_pressure_, w + " header pressure", error_msg);
  }

  double TankValve::FlowIn(const TankLiveState& s) const {
    const double dP = s.pressure - downstream_pressure_;
    if (check_valve_) return -cv_ * detail::SmoothPositive(dP, kTankValveSmoothingPressure) / 1e5;
    return -cv_ * dP / 1e5;
  }

  bool TankValve::PerformChecks(std::string& error_msg, const Fluid& fluid) {
    if (!ITankStream::PerformChecks(error_msg, fluid)) return false;
    const std::string w = "Tank valve '" + name_ + "'";
    return NonNegative(cv_, w + " Cv", error_msg)
        && NonNegative(downstream_pressure_, w + " downstream pressure", error_msg);
  }

  double TankRegulator::Outflow(double P, double cv, double P_set, double P_full, double P_down) {
    const double x = (P - P_set) / (P_full - P_set);
    const double opening = detail::SmoothPositive(x, kTankRegulatorSmoothing)
                         - detail::SmoothPositive(x - 1.0, kTankRegulatorSmoothing);
    return cv * opening * detail::SmoothPositive(P - P_down, kTankValveSmoothingPressure) / 1e5;
  }

  double TankRegulator::FlowIn(const TankLiveState& s) const {
    return -Outflow(s.pressure, cv_, set_pressure_, full_open_pressure_, downstream_pressure_);
  }

  bool TankRegulator::PerformChecks(std::string& error_msg, const Fluid& fluid) {
    if (!ITankStream::PerformChecks(error_msg, fluid)) return false;
    const std::string w = "Regulator '" + name_ + "'";
    if (!NonNegative(cv_, w + " Cv", error_msg)) return false;
    if (!Positive(set_pressure_, w + " set pressure", error_msg)) return false;
    if (!(full_open_pressure_ > set_pressure_)) {
      error_msg += w + ": the full-open pressure must exceed the set pressure. ";
      return false;
    }
    return NonNegative(downstream_pressure_, w + " downstream pressure", error_msg);
  }

  double TankHold::FlowIn(const TankLiveState& s) const {
    const double P_set = set_pressure_;
    return -TankRegulator::Outflow(s.pressure, cv_, P_set, P_set * (1.0 + band_), downstream_pressure_);
  }

  bool TankHold::PerformChecks(std::string& error_msg, const Fluid& fluid) {
    if (!ITankStream::PerformChecks(error_msg, fluid)) return false;
    const std::string w = "Pressure hold '" + name_ + "'";
    return Positive(set_pressure_, w + " set pressure", error_msg)
        && Positive(band_, w + " band", error_msg)
        && Positive(cv_, w + " Cv", error_msg)
        && NonNegative(downstream_pressure_, w + " downstream pressure", error_msg);
  }

  std::unique_ptr<ITankStream> CreateTankStream(TankStreamKind kind) {
    switch (kind) {
      case TankStreamKind::MAKE_UP: return std::make_unique<TankMakeUp>();
      case TankStreamKind::VALVE: return std::make_unique<TankValve>();
      case TankStreamKind::REGULATOR: return std::make_unique<TankRegulator>();
      case TankStreamKind::HOLD: return std::make_unique<TankHold>();
    }
    return nullptr;
  }

  bool Tank::PerformChecks(std::string& error_msg, const Fluid& fluid) {
    const std::string w = "Tank '" + name_ + "'";
    if (beds_ != 1) {
      error_msg += w + ": a tank shared by " + std::to_string(beds_)
                 + " beds is not supported yet (beds must be 1). ";
      return false;
    }
    if (!Positive(volume_, w + " volume", error_msg)) return false;
    if (!Positive(initial_pressure_, w + " initial pressure", error_msg)) return false;
    if (!Positive(initial_temperature_, w + " initial temperature", error_msg)) return false;
    if (!CheckComposition(initial_mole_fraction_, fluid, w + " initial composition", error_msg)) return false;
    switch (thermal_) {
      case TankThermalMode::ISOTHERMAL:
      case TankThermalMode::ADIABATIC:
        break;
      case TankThermalMode::HEAT_LOSS:
        if (!NonNegative(ua_, w + " UA", error_msg)) return false;
        if (!Positive(ambient_temperature_, w + " ambient temperature", error_msg)) return false;
        break;
      case TankThermalMode::WALL:
        if (!NonNegative(wall_area_, w + " wall area", error_msg)) return false;
        if (!NonNegative(h_in_, w + " inner heat-transfer coefficient", error_msg)) return false;
        if (!NonNegative(h_out_, w + " outer heat-transfer coefficient", error_msg)) return false;
        if (!Positive(wall_heat_capacity_, w + " wall heat capacity", error_msg)) return false;
        if (!Positive(ambient_temperature_, w + " ambient temperature", error_msg)) return false;
        if (!Positive(initial_wall_temperature_, w + " initial wall temperature", error_msg)) return false;
        break;
      default:
        error_msg += w + ": unknown thermal mode. ";
        return false;
    }
    for (auto& stream : streams_) {
      if (!stream) {
        error_msg += w + ": empty stream. ";
        return false;
      }
      if (!stream->PerformChecks(error_msg, fluid)) {
        error_msg += "(" + w + ") ";
        return false;
      }
    }
    return true;
  }

  void Tank::EvaluateExpressions(const EvalContext& ctx) {
    volume_.Evaluate(ctx);
    initial_pressure_.Evaluate(ctx);
    initial_temperature_.Evaluate(ctx);
    EvaluateExpressionMap(initial_mole_fraction_, ctx);
    ambient_temperature_.Evaluate(ctx);
    ua_.Evaluate(ctx);
    wall_area_.Evaluate(ctx);
    h_in_.Evaluate(ctx);
    h_out_.Evaluate(ctx);
    wall_heat_capacity_.Evaluate(ctx);
    initial_wall_temperature_.Evaluate(ctx);
    for (auto& stream : streams_)
      if (stream) stream->EvaluateExpressions(ctx);
  }

  double Tank::InitialTotalMoles() const {
    return initial_pressure_ * volume_ / (constants::kIdealGasConstant * initial_temperature_);
  }

  ITankStream* Tank::FindStream(size_t endpoint_id) {
    if (endpoint_id == 0) return nullptr;
    for (auto& s : streams_)
      if (s && s->endpoint_id_ == endpoint_id) return s.get();
    return nullptr;
  }

  size_t TankDirectory::Add(const std::string& name) {
    for (const auto& t : tanks_) next_id_ = std::max(next_id_, t.id_);
    Tank tank(name);
    tank.id_ = ++next_id_;
    tanks_.push_back(std::move(tank));
    return tanks_.back().id_;
  }

  void TankDirectory::Remove(size_t id) {
    tanks_.erase(std::remove_if(tanks_.begin(), tanks_.end(),
                                [id](const Tank& t) { return t.id_ == id; }),
                 tanks_.end());
  }

  Tank* TankDirectory::Find(size_t id) {
    if (id == 0) return nullptr;
    for (auto& t : tanks_)
      if (t.id_ == id) return &t;
    return nullptr;
  }

  const Tank* TankDirectory::Find(size_t id) const {
    if (id == 0) return nullptr;
    for (const auto& t : tanks_)
      if (t.id_ == id) return &t;
    return nullptr;
  }

  Tank* TankDirectory::FindByName(const std::string& name) {
    for (auto& t : tanks_)
      if (t.name_ == name) return &t;
    return nullptr;
  }

} // end adsim namespace

CEREAL_REGISTER_TYPE(adsim::TankMakeUp)
CEREAL_REGISTER_TYPE(adsim::TankValve)
CEREAL_REGISTER_TYPE(adsim::TankRegulator)
CEREAL_REGISTER_TYPE(adsim::TankHold)
CEREAL_REGISTER_POLYMORPHIC_RELATION(adsim::ITankStream, adsim::TankMakeUp)
CEREAL_REGISTER_POLYMORPHIC_RELATION(adsim::ITankStream, adsim::TankValve)
CEREAL_REGISTER_POLYMORPHIC_RELATION(adsim::ITankStream, adsim::TankRegulator)
CEREAL_REGISTER_POLYMORPHIC_RELATION(adsim::ITankStream, adsim::TankHold)
