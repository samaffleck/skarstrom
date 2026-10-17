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

  std::vector<std::string> Step::GetEndpointNames() const {
    std::vector<std::string> names;
    CollectEndpointNames(top_boundary_condition_.get(), names);
    CollectEndpointNames(bottom_boundary_condition_.get(), names);
    return names;
  }

} // end adsim namespace

CEREAL_CLASS_VERSION(adsim::Step, 5);
