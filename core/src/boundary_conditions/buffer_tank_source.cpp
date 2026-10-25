/*
 * Skarstrom - pressure and temperature swing adsorption simulator
 * Copyright (c) 2024-2026 Sam Affleck and contributors
 *
 * Licensed under the terms in the LICENSE file at the repository root.
 */
#include "../pch.h"
#include "adsim/boundary_conditions/buffer_tank_source.h"
#include "adsim/boundary_conditions/boundary_stream.h"
#include "adsim/boundary_conditions/connection.h"


namespace adsim {

  void BufferTankSource::Activate(const std::string& source_step_name,
                                  const std::string& source_endpoint_name) {
    active_ = true;
    source_step_name_ = source_step_name;
    source_endpoint_name_ = source_endpoint_name;
    // Named source: Cycle::ResolveBufferTankReferences fills the IDs in.
    source_step_id_ = 0;
    source_endpoint_id_ = 0;
  }

  void BufferTankSource::Deactivate() {
    active_ = false;
    summary_ = nullptr;
  }

  bool BufferTankSource::PerformChecks(std::string& error_msg) const {
    if (!active_) return true;

    const double pct = percentage_.Get();
    if (pct < 0.0 || pct > 100.0) {
      error_msg += "Buffer tank percentage must be between 0 and 100 (got "
                + std::to_string(pct) + "). ";
      return false;
    }
    if (source_step_name_.empty()) {
      error_msg += "Buffer tank inlet has no source step selected. ";
      return false;
    }
    if (source_endpoint_name_.empty()) {
      error_msg += "Buffer tank inlet has no source endpoint selected. ";
      return false;
    }
    return true;
  }

  double BufferTankSource::GetMolarFlow(double consumer_step_time) const {
    if (!IsReady() || consumer_step_time <= 0.0) return 0.0;
    return (percentage_.Get() / 100.0) * summary_->total_moles / consumer_step_time;
  }

  double BufferTankSource::GetTemperature() const {
    if (!IsReady()) return 0.0;
    return summary_->avg_temperature;
  }

  double BufferTankSource::GetMoleFraction(const std::string& component) const {
    if (!IsReady()) return 0.0;
    const auto it = summary_->avg_mole_fraction.find(component);
    if (it == summary_->avg_mole_fraction.end()) return 0.0;
    return it->second;
  }

  namespace {
    // True once an interacting-step connection has logged data to read.
    bool ConnectionHasData(const Connection& connection) {
      return connection.IsConnected() &&
             connection.stream_ != nullptr &&
             !connection.stream_->molar_flow_.empty();
    }
  } // anonymous namespace

  double SelectInletTemperature(const BufferTankSource& buffer_tank,
                                const Connection& connection,
                                double assumed_temperature,
                                double time) {
    if (buffer_tank.IsActive()) {
      if (buffer_tank.IsReady()) return buffer_tank.GetTemperature();
    } else if (ConnectionHasData(connection)) {
      return connection.stream_->GetTemperature(time);
    }
    return assumed_temperature; // First-cycle assumption / unconnected value
  }

  double SelectInletMoleFraction(const BufferTankSource& buffer_tank,
                                 const Connection& connection,
                                 const ExpressionMap& assumed_mole_fraction,
                                 const std::string& component,
                                 double time) {
    if (buffer_tank.IsActive()) {
      if (buffer_tank.IsReady()) return buffer_tank.GetMoleFraction(component);
    } else if (ConnectionHasData(connection)) {
      return connection.stream_->GetMoleFraction(time, component);
    }
    // First-cycle assumption / unconnected value
    const auto it = assumed_mole_fraction.find(component);
    return (it != assumed_mole_fraction.end()) ? it->second.Get() : 0.0;
  }

} // end adsim namespace
