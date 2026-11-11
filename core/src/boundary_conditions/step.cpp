/*
 * Skarstrom - pressure and temperature swing adsorption simulator
 * Copyright (c) 2024-2026 Sam Affleck and contributors
 *
 * Licensed under the terms in the LICENSE file at the repository root.
 */
#include "../pch.h"
#include "adsim/boundary_conditions/step.h"
#include "adsim/equations/equation.h"


namespace adsim {

  bool HeatingMedium::PerformChecks(const std::string& step_name, std::string& error_message) const {
    if (!enabled_) return true;
    const std::string where = "the heating medium of step '" + step_name + "'";
    const double h = heat_transfer_coefficient_.Get();
    if (!std::isfinite(h) || h < 0.0) {
      error_message += "The heat transfer coefficient of " + where + " is " + std::to_string(h)
                       + " W/(m2 K); it must be zero or positive. ";
      return false;
    }
    const std::pair<const char*, double> exchanger[] = {
        {"internal exchanger area per volume", exchanger_area_per_volume_.Get()},
        {"internal exchanger heat transfer coefficient", exchanger_heat_transfer_coefficient_.Get()}};
    for (const auto& item : exchanger) {
      if (!std::isfinite(item.second) || item.second < 0.0) {
        error_message += std::string("The ") + item.first + " of " + where + " is " + std::to_string(item.second)
                         + "; it must be zero or positive. ";
        return false;
      }
    }
    if (!CheckProfileFactor(temperature_profile_.get(), "temperature profile of " + where, error_message))
      return false;
    // A ramp's start and end bound every temperature it takes (all the
    // profile kinds are monotonic between them).
    std::vector<double> temperatures;
    if (auto* l = dynamic_cast<const LinearBoundaryValue*>(temperature_profile_.get()))
      temperatures = {l->start_value_.Get(), l->end_value_.Get()};
    else if (auto* e = dynamic_cast<const ExponentialBoundaryValue*>(temperature_profile_.get()))
      temperatures = {e->start_value_.Get(), e->end_value_.Get()};
    else if (auto* g = dynamic_cast<const LogLinearBoundaryValue*>(temperature_profile_.get()))
      temperatures = {g->start_value_.Get(), g->end_value_.Get()};
    else if (auto* d = dynamic_cast<const DoubleExponentialBoundaryValue*>(temperature_profile_.get()))
      temperatures = {d->start_value_.Get(), d->end_value_.Get()};
    else
      temperatures = {temperature_.Get()};
    for (double T : temperatures) {
      if (!std::isfinite(T) || T <= 0.0) {
        error_message += "The temperature of " + where + " is " + std::to_string(T)
                         + " K; it must be above absolute zero. ";
        return false;
      }
    }
    return true;
  }

  Step::Step() {
    auto top_boundary = std::make_shared<PressureBoundaryCondition>();
    top_boundary->location_ = BoundaryLocation::TOP;
    top_boundary->direction_ = BoundaryDirection::OUTLET;
    top_boundary_condition_ = std::move(top_boundary);

    auto bottom_boundary = std::make_shared<FlowBoundaryCondition>();
    bottom_boundary->location_ = BoundaryLocation::BOTTOM;
    auto stream = std::make_unique<BoundaryStream>();
    bottom_boundary->direction_ = BoundaryDirection::INLET;
    stream->name_ = "Stream 1";
    bottom_boundary->streams_.emplace_back(std::move(stream));
    bottom_boundary_condition_ = std::move(bottom_boundary);
  }

  void Step::Initialise(Equations& state) {
    bottom_boundary_condition_->ClearLog();
    top_boundary_condition_->ClearLog();
  
    bottom_boundary_condition_->location_ = BoundaryLocation::BOTTOM;
    bottom_boundary_condition_->Initialise(state, step_time_);

    top_boundary_condition_->location_ = BoundaryLocation::TOP;
    top_boundary_condition_->Initialise(state, step_time_);
    
    state.SetCurrentStep(this);
  }

  void Step::LogOutletConditions(Equations& state, double time) {
    top_boundary_condition_->Log(state, time);
    bottom_boundary_condition_->Log(state, time);
  }

  void Step::ClearOutletConditions() {
    top_boundary_condition_->ClearLog();
    bottom_boundary_condition_->ClearLog();
  }

  namespace {
    void UpdateSummariesOnBC(IBoundaryCondition* bc, int cycle_index) {
      if (!bc) return;
      if (auto* fbc = dynamic_cast<FlowBoundaryCondition*>(bc)) {
        for (auto& stream : fbc->streams_) {
          if (!stream) continue;
          stream->cycle_summary_ = ComputeSummary(stream->stream_data_, cycle_index);
        }
      } else if (auto* pbc = dynamic_cast<PressureBoundaryCondition*>(bc)) {
        pbc->cycle_summary_ = ComputeSummary(pbc->stream_data_, cycle_index);
      }
    }

    void ClearSummariesOnBC(IBoundaryCondition* bc) {
      if (!bc) return;
      if (auto* fbc = dynamic_cast<FlowBoundaryCondition*>(bc)) {
        for (auto& stream : fbc->streams_) {
          if (stream) stream->cycle_summary_.Clear();
        }
      } else if (auto* pbc = dynamic_cast<PressureBoundaryCondition*>(bc)) {
        pbc->cycle_summary_.Clear();
      }
    }

    CycleEndpointSummary* FindSummaryOnBC(IBoundaryCondition* bc,
                                          const std::string& endpoint_name) {
      if (!bc) return nullptr;
      if (auto* fbc = dynamic_cast<FlowBoundaryCondition*>(bc)) {
        for (auto& stream : fbc->streams_) {
          if (stream && stream->name_ == endpoint_name) {
            return &stream->cycle_summary_;
          }
        }
      } else if (auto* pbc = dynamic_cast<PressureBoundaryCondition*>(bc)) {
        if (endpoint_name == kPressureBoundaryEndpointName) {
          return &pbc->cycle_summary_;
        }
      }
      return nullptr;
    }

    // KPI endpoints: a flow stream is named by its stream name; a pressure
    // boundary is named to MATCH the BoundaryLog key produced by
    // GetAllBoundaryStreams ("<step> Pressure Boundary <Top|Bottom>") so that
    // KPI lookups (GetVariableData) resolve the logged per-endpoint flow.
    void CollectEndpoints(const std::string& step_name, const IBoundaryCondition* bc,
                          std::vector<Step::Endpoint>& out) {
      if (!bc) return;
      if (auto* fbc = dynamic_cast<const FlowBoundaryCondition*>(bc)) {
        for (const auto& stream : fbc->streams_) {
          if (stream) out.push_back({stream->name_, bc->location_, stream->endpoint_id_});
        }
      } else if (auto* pbc = dynamic_cast<const PressureBoundaryCondition*>(bc)) {
        out.push_back({Step::PressureEndpointName(step_name, bc->location_), bc->location_, pbc->endpoint_id_});
      }
    }

    // Endpoint names for BufferTankSource / cycle-summary resolution: a pressure
    // boundary uses the canonical kPressureBoundaryEndpointName (what
    // FindSummaryForEndpoint matches against), NOT the log-key form above.
    void CollectEndpointNames(const IBoundaryCondition* bc, std::vector<std::string>& out) {
      if (!bc) return;
      if (auto* fbc = dynamic_cast<const FlowBoundaryCondition*>(bc)) {
        for (const auto& stream : fbc->streams_) {
          if (stream) out.push_back(stream->name_);
        }
      } else if (dynamic_cast<const PressureBoundaryCondition*>(bc)) {
        out.emplace_back(kPressureBoundaryEndpointName);
      }
    }
  } // anonymous namespace

  void Step::UpdateCycleSummary(int cycle_index) {
    UpdateSummariesOnBC(top_boundary_condition_.get(), cycle_index);
    UpdateSummariesOnBC(bottom_boundary_condition_.get(), cycle_index);
  }

  void Step::ClearCycleSummary() {
    ClearSummariesOnBC(top_boundary_condition_.get());
    ClearSummariesOnBC(bottom_boundary_condition_.get());
  }

  CycleEndpointSummary* Step::FindSummaryForEndpoint(const std::string& endpoint_name) {
    if (auto* s = FindSummaryOnBC(top_boundary_condition_.get(), endpoint_name)) return s;
    return FindSummaryOnBC(bottom_boundary_condition_.get(), endpoint_name);
  }

  std::vector<Step::Endpoint> Step::GetEndpoints() const {
    std::vector<Step::Endpoint> endpoints;
    CollectEndpoints(step_name_, top_boundary_condition_.get(), endpoints);
    CollectEndpoints(step_name_, bottom_boundary_condition_.get(), endpoints);
    return endpoints;
  }

  std::vector<Step::LoggedEndpoint> Step::LoggedEndpoints() {
    std::vector<LoggedEndpoint> out;
    auto add = [&](IBoundaryCondition* bc) {
      if (!bc) return;
      if (auto* fbc = dynamic_cast<FlowBoundaryCondition*>(bc)) {
        for (auto& stream : fbc->streams_) {
          if (stream) out.push_back({&stream->stream_data_, bc->location_, stream.get(), nullptr});
        }
      } else if (auto* pbc = dynamic_cast<PressureBoundaryCondition*>(bc)) {
        out.push_back({&pbc->stream_data_, bc->location_, nullptr, pbc});
      }
    };
    add(bottom_boundary_condition_.get());
    add(top_boundary_condition_.get());
    return out;
  }

  void Step::SampleEndpoints(const Equations& state, double time, std::vector<EndpointSample>& out) {
    // The outlet node mixing of each end (FlowBoundaryCondition::Mixing), as
    // FlowBoundaryCondition::Log applies it.
    auto mixing = [&](IBoundaryCondition* bc) {
      const auto* fbc = dynamic_cast<const FlowBoundaryCondition*>(bc);
      return fbc ? fbc->Mixing(state, time) : NodeMixing{};
    };
    const NodeMixing bottom = mixing(bottom_boundary_condition_.get());
    const NodeMixing top = mixing(top_boundary_condition_.get());
    for (const auto& ep : LoggedEndpoints()) {
      if (ep.stream) out.push_back(ep.stream->Sample(state, time, ep.location,
                                                     ep.location == BoundaryLocation::TOP ? &top : &bottom));
      else if (ep.pressure) out.push_back(ep.pressure->Sample(state, time));
    }
  }

  int Step::NumberOfOpenEnds() const {
    auto open = [](const IBoundaryCondition* bc) {
      return bc && !dynamic_cast<const ClosedBoundaryCondition*>(bc);
    };
    return (open(bottom_boundary_condition_.get()) ? 1 : 0) + (open(top_boundary_condition_.get()) ? 1 : 0);
  }

  std::vector<std::string> Step::GetEndpointNames() const {
    std::vector<std::string> names;
    CollectEndpointNames(top_boundary_condition_.get(), names);
    CollectEndpointNames(bottom_boundary_condition_.get(), names);
    return names;
  }

} // end adsim namespace
