/*
 * Skarstrom - pressure and temperature swing adsorption simulator
 * Copyright (c) 2024-2026 Sam Affleck and contributors
 *
 * Licensed under the terms in the LICENSE file at the repository root.
 */
#include "../pch.h"
#include "adsim/solver/cycle.h"
#include "adsim/boundary_conditions/boundary_conditions.h"
#include "adsim/boundary_conditions/boundary_stream.h"

#include <algorithm>
#include <map>
#include <unordered_set>
#include <sstream>


namespace adsim {

  namespace {
    // Walks every stream / pressure boundary on a step and invokes
    // visit(buffer_tank, owner_step, endpoint_name, conn) on each one. The
    // Connection pointer is provided so callers can enforce mutual-exclusion
    // checks against the existing interacting-step Connection.
    template <typename Visitor>
    void ForEachBufferTank(Step& step, Visitor visit) {
      auto handle = [&](IBoundaryCondition* bc) {
        if (!bc) return;
        if (auto* fbc = dynamic_cast<FlowBoundaryCondition*>(bc)) {
          for (auto& stream : fbc->streams_) {
            if (!stream) continue;
            visit(stream->buffer_tank_, step, stream->name_, &stream->connection_);
          }
        } else if (auto* pbc = dynamic_cast<PressureBoundaryCondition*>(bc)) {
          visit(pbc->buffer_tank_, step,
                std::string(kPressureBoundaryEndpointName), &pbc->connection_);
        }
      };
      handle(step.top_boundary_condition_.get());
      handle(step.bottom_boundary_condition_.get());
    }

    // Invokes visit(conn, stream_data, endpoint_name) for every connectable
    // endpoint carried by `bc` (its pressure endpoint, or each flow stream).
    template <typename Visitor>
    void ForEachConnection(IBoundaryCondition* bc, Visitor visit) {
      if (auto* fbc = dynamic_cast<FlowBoundaryCondition*>(bc)) {
        for (auto& stream : fbc->streams_) {
          if (stream) visit(stream->connection_, stream->stream_data_, stream->name_);
        }
      } else if (auto* pbc = dynamic_cast<PressureBoundaryCondition*>(bc)) {
        visit(pbc->connection_, pbc->stream_data_, std::string(kPressureBoundaryEndpointName));
      }
    }

    // Every CONNECTED endpoint must have exactly one partner (same pair id) and
    // point at that partner's StreamData. An endpoint whose partner was
    // destroyed keeps a dangling StreamData* the solver would dereference.
    // Read-only: GetConnectedEndpoints would silently disconnect such orphans.
    bool CheckConnectionPartners(std::vector<Step>& step_list, std::string& error_message) {
      struct Seen { const Connection* conn; const StreamData* data; std::string label; };
      std::map<size_t, std::vector<Seen>> by_id;
      for (Step& step : step_list) {
        auto collect = [&](Connection& conn, StreamData& data, const std::string& name) {
          if (conn.IsConnected())
            by_id[conn.GetID()].push_back({&conn, &data, "step '" + step.step_name_ + "' endpoint '" + name + "'"});
        };
        ForEachConnection(step.top_boundary_condition_.get(), collect);
        ForEachConnection(step.bottom_boundary_condition_.get(), collect);
      }
      for (const auto& [id, eps] : by_id) {
        const bool paired = eps.size() == 2
          && eps[0].conn->stream_ == eps[1].data
          && eps[1].conn->stream_ == eps[0].data;
        if (!paired) {
          error_message += "Interacting-step connection on " + eps[0].label
                        + " has no valid partner endpoint (it was removed or replaced). "
                          "Disconnect and reconnect it. ";
          return false;
        }
      }
      return true;
    }
  } // anonymous namespace

  bool Cycle::PerformChecks(const Fluid& fluid, std::string& error_message) {
    ResolveBufferTankReferences();  // also EnsureIds; refreshes source names from IDs
    if (step_time_mode_ == StepTimeMode::FRACTION) {
      double sum = 0.0;
      for (const auto& step : step_list_) {
        sum += step.step_fraction_.Get();
      }
      if (std::abs(sum - 1.0) > 1e-6) {
        error_message += "Step fractions must sum to 1.0, but currently sum to " + std::to_string(sum) + ". ";
        return false;
      }
    }

    ResolveStepTimes();

    if (step_list_.empty()) {
      error_message += "Cycle steps are empty. You need to add at least one step to the cycle. ";
      return false;
    }

    // Every step must have a strictly positive duration. A zero-length step makes
    // the chunked solver spin forever (its integration loop never runs, so the
    // step never completes) and the non-chunked solver fail with IDA_ILL_INPUT.
    // Checked after ResolveStepTimes() so FRACTION mode is covered too.
    for (const Step& s : step_list_) {
      const double t = s.step_time_.Get();
      if (!(t > 0.0)) { // also rejects NaN
        error_message += "Step '" + s.step_name_ + "' has a step time of " + std::to_string(t)
                        + " s. Every step must have a step time greater than zero. ";
        return false;
      }
    }

    std::unordered_set<std::string> step_names;
    step_names.reserve(step_list_.size());
    for (const Step& s : step_list_) {
      if (!step_names.insert(s.step_name_).second) {
        error_message += "Duplicate step name '" + s.step_name_ + "'. Step names must be unique. ";
        return false;
      }
    }

    // TODO: Check the first step. If one of the bc's supplies a pressure, change the IC 


    // For each step in the cycle, make sure both steps aren't INLETS
    for (const auto& step : step_list_) {
      auto& top_bc = step.top_boundary_condition_;
      auto& bot_bc = step.bottom_boundary_condition_;

      // Report a missing boundary condition instead of dereferencing it (here
      // and in the location check below).
      if (!top_bc || !bot_bc) {
        error_message += std::string(!top_bc ? "Top" : "Bottom") + " boundary condition of step '"
                        + step.step_name_ + "' is null. Every step needs a top and a bottom boundary condition. ";
        return false;
      }

      if (top_bc->direction_ == BoundaryDirection::INLET && bot_bc->direction_ == BoundaryDirection::INLET) {
        error_message += "Both boundary conditions are of type 'INLET' in " + step.step_name_ + ". You must have one OUTLET. ";
        return false;
      }
    }

    auto validate_boundary = [&](IBoundaryCondition* bc_ptr,
                                std::string_view which,
                                const Step& step) -> bool {
      if (bc_ptr == nullptr) {
        error_message += std::string(which) + " boundary condition of step '" + step.step_name_ + "' is null. ";
        return false;
      }

      if (!bc_ptr->PerformChecks(error_message, fluid)) {
        error_message += std::string(which) + " boundary condition of step '" + step.step_name_ + "'. ";
        return false;
      }

      if (auto fbc = dynamic_cast<const FlowBoundaryCondition*>(bc_ptr)) {
        for (const auto& stream_ptr : fbc->streams_) {
          if (auto valve = dynamic_cast<const ValveBoundaryStream*>(stream_ptr.get())) {
            if (valve->connection_.IsConnected()) {
              // Check that the connection has valid stream data
              if (valve->connection_.stream_ == nullptr) {
                error_message += std::string(which) + " valve stream in step '" + step.step_name_
                                + "' is CONNECTED but has no stream data. ";
                return false;
              }
            }
          }
        }
      } else if (auto pbc = dynamic_cast<const PressureBoundaryCondition*>(bc_ptr)) {
        if (pbc->connection_.IsConnected()) {
          // Check that the connection has valid stream data
          if (pbc->connection_.stream_ == nullptr) {
            error_message += std::string(which) + " pressure boundary in step '" + step.step_name_
                            + "' is CONNECTED but has no stream data. ";
            return false;
          }
        }
      }

      return true;
    };

    for (const Step& step : step_list_) {
      if (step.bottom_boundary_condition_->location_ == step.top_boundary_condition_->location_) {
        error_message += "In step '" + step.step_name_
                        + "', top and bottom boundary conditions have the same location; they must be different. ";
        return false;
      }

      if (!validate_boundary(step.top_boundary_condition_.get(), "top", step))    return false;
      if (!validate_boundary(step.bottom_boundary_condition_.get(), "bottom", step)) return false;
    }

    if (!CheckConnectionPartners(step_list_, error_message)) return false;

    auto pairs = GetConnectedEndpoints(step_list_);
    for (const auto& ep : pairs) {
      auto& ep1 = ep.first;
      auto& ep2 = ep.second;
      if (std::abs(ep1.step->step_time_.Get() - ep2.step->step_time_.Get()) > 0.01) {
        error_message += "Interacting steps don't have an equal step time. ";
        error_message += " Time of " + ep1.step->step_name_ + " = " + std::to_string(ep1.step->step_time_.Get()) + ". ";
        error_message += " Time of " + ep2.step->step_name_ + " = " + std::to_string(ep2.step->step_time_.Get()) + ". ";
        return false;
      }
    }

    // Bind every active BufferTankSource to its source endpoint before we
    // validate references. Resolution silently leaves summary_ == nullptr
    // for unresolved references; the loop below turns those into errors.
    ResolveBufferTankSources();

    for (Step& step : step_list_) {
      bool step_ok = true;
      ForEachBufferTank(step,
        [&](BufferTankSource& bt, Step& owner, const std::string& owner_endpoint,
            Connection* conn) {
          if (!bt.IsActive()) return;

          if (!bt.PerformChecks(error_message)) {
            error_message += "(step '" + owner.step_name_ + "', endpoint '" + owner_endpoint + "'). ";
            step_ok = false;
            return;
          }
          if (conn && conn->IsConnected()) {
            error_message += "Step '" + owner.step_name_ + "' endpoint '" + owner_endpoint
                          + "' has both an interacting-step connection and a buffer tank source. "
                            "Use only one. ";
            step_ok = false;
            return;
          }
          if (bt.source_step_id_ != 0 ? bt.source_step_id_ == owner.id_
                                      : bt.source_step_name_ == owner.step_name_) {
            error_message += "Buffer tank on step '" + owner.step_name_
                          + "' cannot reference itself. ";
            step_ok = false;
            return;
          }
          Step* src = FindStep(bt.source_step_name_);
          if (!src) {
            error_message += "Buffer tank on step '" + owner.step_name_
                          + "' references unknown source step '" + bt.source_step_name_ + "'. ";
            step_ok = false;
            return;
          }
          if (!src->FindSummaryForEndpoint(bt.source_endpoint_name_)) {
            error_message += "Buffer tank on step '" + owner.step_name_
                          + "' references unknown endpoint '" + bt.source_endpoint_name_
                          + "' on step '" + bt.source_step_name_ + "'. ";
            step_ok = false;
            return;
          }
        });
      if (!step_ok) return false;
    }

    return true;
  }

  void Cycle::AddStep(const std::string& step_name) {
    auto it = std::find_if(step_list_.begin(), step_list_.end(),
      [&step_name](const Step& step) {
      return step.step_name_ == step_name;
    });

    if (it == step_list_.end()) {
      Step newStep;
      newStep.step_name_ = step_name;
      step_list_.emplace_back(newStep);
      EnsureIds();
    }
  }


  void Cycle::AddStep(const std::string& step_name, size_t pos) {
    for (const auto& step : step_list_) {
      if (step.step_name_ == step_name)
        return; // Already exists
    }

    Step newStep;
    newStep.step_name_ = step_name;
    step_list_.insert(step_list_.begin() + pos, newStep);
    EnsureIds();
  }


  void Cycle::RemoveStep(const std::string& step_name) {
    auto it = std::find_if(step_list_.begin(), step_list_.end(),
      [&step_name](const Step& step) {
      return step.step_name_ == step_name;
    });

    if (it != step_list_.end()) {
      DisconnectBoundaryCondition(step_list_, it->top_boundary_condition_.get());
      DisconnectBoundaryCondition(step_list_, it->bottom_boundary_condition_.get());
      step_list_.erase(it);
    }
  }

  void Cycle::ReplaceBoundaryCondition(std::shared_ptr<IBoundaryCondition>& slot,
                                       std::shared_ptr<IBoundaryCondition> replacement) {
    DisconnectBoundaryCondition(step_list_, slot.get());
    slot = std::move(replacement);
    EnsureIds();  // a different boundary type exposes new endpoints
  }

  void Cycle::ReplaceStream(std::unique_ptr<IBoundaryStream>& slot,
                            std::unique_ptr<IBoundaryStream> replacement) {
    // Changing a stream's device type keeps the endpoint (the GUI carries
    // its name, temperature and composition over), so it keeps its ID and
    // every reference to it.
    const size_t keep_id = slot ? slot->endpoint_id_ : 0;
    if (slot) DisconnectWithPartner(step_list_, slot->connection_);
    slot = std::move(replacement);
    if (slot && slot->endpoint_id_ == 0) slot->endpoint_id_ = keep_id;
    EnsureIds();
  }

  // ---- Identity ------------------------------------------------------------

  namespace {
    // Visits every endpoint ID slot in assignment order: steps in order, each
    // step's bottom boundary before its top, a flow boundary's streams in
    // order, a pressure boundary as one endpoint.
    template <class Fn>
    void ForEachEndpointId(std::vector<Step>& steps, Fn fn) {
      for (Step& step : steps) {
        for (IBoundaryCondition* bc : {step.bottom_boundary_condition_.get(), step.top_boundary_condition_.get()}) {
          if (auto* fbc = dynamic_cast<FlowBoundaryCondition*>(bc)) {
            for (auto& s : fbc->streams_)
              if (s) fn(s->endpoint_id_);
          } else if (auto* pbc = dynamic_cast<PressureBoundaryCondition*>(bc)) {
            fn(pbc->endpoint_id_);
          }
        }
      }
    }

    // Gives every 0 / repeated ID in `slots` (in order) the next free ID.
    void AssignIds(const std::vector<size_t*>& slots) {
      size_t next = 0;
      for (size_t* id : slots) next = std::max(next, *id);
      std::unordered_set<size_t> seen;
      for (size_t* id : slots) {
        if (*id == 0 || !seen.insert(*id).second) {
          *id = ++next;
          seen.insert(*id);
        }
      }
    }
  }  // namespace

  void Cycle::EnsureIds() {
    std::vector<size_t*> step_ids;
    for (Step& s : step_list_) step_ids.push_back(&s.id_);
    AssignIds(step_ids);
    std::vector<size_t*> endpoint_ids;
    ForEachEndpointId(step_list_, [&](size_t& id) { endpoint_ids.push_back(&id); });
    AssignIds(endpoint_ids);
  }

  std::vector<Cycle::EndpointHandle> Cycle::Endpoints() {
    std::vector<EndpointHandle> out;
    for (size_t i = 0; i < step_list_.size(); ++i) {
      Step& step = step_list_[i];
      for (IBoundaryCondition* bc : {step.bottom_boundary_condition_.get(), step.top_boundary_condition_.get()}) {
        if (auto* fbc = dynamic_cast<FlowBoundaryCondition*>(bc)) {
          for (auto& s : fbc->streams_) {
            if (!s) continue;
            EndpointHandle h;
            h.id = s->endpoint_id_;
            h.step_index = i;
            h.step = &step;
            h.boundary = bc;
            h.stream = s.get();
            h.location = bc->location_;
            h.log_name = s->name_;
            h.summary_name = s->name_;
            out.push_back(std::move(h));
          }
        } else if (auto* pbc = dynamic_cast<PressureBoundaryCondition*>(bc)) {
          EndpointHandle h;
          h.id = pbc->endpoint_id_;
          h.step_index = i;
          h.step = &step;
          h.boundary = bc;
          h.pressure = pbc;
          h.location = bc->location_;
          h.log_name = Step::PressureEndpointName(step.step_name_, bc->location_);
          h.summary_name = kPressureBoundaryEndpointName;
          out.push_back(std::move(h));
        }
      }
    }
    return out;
  }

  std::optional<Cycle::EndpointHandle> Cycle::FindEndpoint(size_t endpoint_id) {
    if (endpoint_id == 0) return std::nullopt;
    for (auto& h : Endpoints())
      if (h.id == endpoint_id) return h;
    return std::nullopt;
  }

  Step* Cycle::FindStepById(size_t step_id) {
    const int i = StepIndexById(step_id);
    return i < 0 ? nullptr : &step_list_[static_cast<size_t>(i)];
  }

  int Cycle::StepIndexById(size_t step_id) const {
    if (step_id == 0) return -1;
    for (size_t i = 0; i < step_list_.size(); ++i)
      if (step_list_[i].id_ == step_id) return static_cast<int>(i);
    return -1;
  }

  void Cycle::RemoveStream(FlowBoundaryCondition& bc, size_t index) {
    if (index >= bc.streams_.size()) return;
    if (bc.streams_[index]) DisconnectWithPartner(step_list_, bc.streams_[index]->connection_);
    bc.streams_.erase(bc.streams_.begin() + static_cast<std::ptrdiff_t>(index));
  }

  namespace {
    // Clears every Connection on a step so it starts detached from the rest
    // of the cycle. We use Disconnect() (not DisconnectKeepID) to reset the
    // pair id — a duplicate must not collide with the original's pairing.
    void DisconnectAllOnStep(Step& step) {
      auto clear_bc = [](IBoundaryCondition* bc) {
        if (auto* flow = dynamic_cast<FlowBoundaryCondition*>(bc)) {
          for (auto& stream : flow->streams_) {
            if (stream) stream->connection_.Disconnect();
          }
        } else if (auto* p = dynamic_cast<PressureBoundaryCondition*>(bc)) {
          p->connection_.Disconnect();
        }
      };
      clear_bc(step.bottom_boundary_condition_.get());
      clear_bc(step.top_boundary_condition_.get());
    }

    std::string MakeUniqueStepName(const std::vector<Step>& steps,
                                   const std::string& base) {
      auto taken = [&](const std::string& candidate) {
        return std::any_of(steps.begin(), steps.end(),
          [&](const Step& s) { return s.step_name_ == candidate; });
      };
      std::string candidate = base + " (copy)";
      if (!taken(candidate)) return candidate;
      for (int n = 2; ; ++n) {
        candidate = base + " (copy " + std::to_string(n) + ")";
        if (!taken(candidate)) return candidate;
      }
    }
  } // anonymous namespace

  size_t Cycle::DuplicateStep(size_t index) {
    if (index >= step_list_.size()) return static_cast<size_t>(-1);

    // Deep-copy via cereal JSON round-trip — handles every polymorphic BC
    // variant (and any future ones) without needing a hand-written Clone().
    // JSON is what the rest of the codebase already registers types with.
    Step copy;
    {
      std::stringstream buf;
      {
        cereal::JSONOutputArchive out(buf);
        out(step_list_[index]);
      }
      cereal::JSONInputArchive in(buf);
      in(copy);
    }

    DisconnectAllOnStep(copy);
    copy.step_name_ = MakeUniqueStepName(step_list_, step_list_[index].step_name_);

    const size_t insert_at = index + 1;
    step_list_.insert(step_list_.begin() + insert_at, std::move(copy));
    // The copy repeats the original's IDs; it comes later, so it is the one
    // that gets new ones.
    EnsureIds();
    return insert_at;
  }


  Step& Cycle::GetStep(const std::string& step_name) {
    for (auto& step : step_list_) {
      if (step.step_name_ == step_name)
        return step;
    }

    throw "Step does not exist";
  }

  Step* Cycle::FindStep(const std::string& step_name) {
    for (auto& step : step_list_) {
      if (step.step_name_ == step_name) return &step;
    }
    return nullptr;
  }

  void Cycle::ResolveBufferTankReferences() {
    EnsureIds();
    const auto endpoints = Endpoints();
    for (Step& step : step_list_) {
      ForEachBufferTank(step,
        [&](BufferTankSource& bt, Step& /*owner*/, const std::string& /*owner_endpoint*/,
            Connection* /*conn*/) {
          if (!bt.IsActive()) return;
          const EndpointHandle* hit = nullptr;
          if (bt.source_endpoint_id_ != 0) {
            for (const auto& h : endpoints)
              if (h.id == bt.source_endpoint_id_) { hit = &h; break; }
          }
          if (!hit) {
            // By name, with FindSummaryForEndpoint's precedence: the first
            // step of that name, its top boundary before its bottom.
            Step* src = FindStep(bt.source_step_name_);
            for (const auto& h : endpoints) {
              if (h.step != src || h.summary_name != bt.source_endpoint_name_) continue;
              if (!hit || (h.location == BoundaryLocation::TOP && hit->location != BoundaryLocation::TOP)) hit = &h;
            }
          }
          if (!hit) return;  // PerformChecks reports the unknown source
          bt.source_step_id_ = hit->step->id_;
          bt.source_endpoint_id_ = hit->id;
          bt.source_step_name_ = hit->step->step_name_;
          bt.source_endpoint_name_ = hit->summary_name;
        });
    }
  }

  void Cycle::ResolveBufferTankSources() {
    ResolveBufferTankReferences();
    for (Step& step : step_list_) {
      ForEachBufferTank(step,
        [&](BufferTankSource& bt, Step& /*owner*/, const std::string& /*owner_endpoint*/,
            Connection* /*conn*/) {
          if (!bt.IsActive()) { bt.summary_ = nullptr; return; }
          bt.summary_ = nullptr;
          if (auto h = FindEndpoint(bt.source_endpoint_id_)) {
            bt.summary_ = h->stream ? &h->stream->cycle_summary_ : &h->pressure->cycle_summary_;
          }
        });
    }
  }

  void Cycle::RenameBufferTankSourceStep(const std::string& old_name, const std::string& new_name) {
    for (Step& step : step_list_) {
      ForEachBufferTank(step,
        [&](BufferTankSource& bt, Step& /*owner*/, const std::string& /*owner_endpoint*/,
            Connection* /*conn*/) {
          if (bt.source_step_name_ == old_name) bt.source_step_name_ = new_name;
        });
    }
  }

  void Cycle::RenameBufferTankSourceEndpoint(const std::string& step_name,
                                             const std::string& old_endpoint,
                                             const std::string& new_endpoint) {
    for (Step& step : step_list_) {
      ForEachBufferTank(step,
        [&](BufferTankSource& bt, Step& /*owner*/, const std::string& /*owner_endpoint*/,
            Connection* /*conn*/) {
          if (bt.source_step_name_ == step_name && bt.source_endpoint_name_ == old_endpoint)
            bt.source_endpoint_name_ = new_endpoint;
        });
    }
  }

  void Cycle::ResolveStepTimes() {
    if (step_time_mode_ == StepTimeMode::FRACTION) {
      for (auto& step : step_list_) {
        step.step_time_ = step.step_fraction_.Get() * total_cycle_time_.Get();
      }
    }
  }

  double Cycle::GetCycleTime() const {
    if (step_time_mode_ == StepTimeMode::FRACTION) {
      return total_cycle_time_.Get();
    }
    double time = 0.0;
    for (const auto& step : step_list_) {
      time += step.step_time_.Get();
    }
    return time;
  }

  void Cycle::EvaluateExpressions(const EvalContext& ctx) {
    total_cycle_time_.Evaluate(ctx);
    for (auto& step : step_list_) {
      step.EvaluateExpressions(ctx);
    }
  }

  double Cycle::GetTotalSimulationTime() const {
    return GetCycleTime() * number_of_cycles_;
  }

  EndpointPairList GetConnectedEndpoints(std::vector<Step>& step_list) {
    EndpointPairList endpoint_list;
    endpoint_list.reserve(step_list.size());

    std::unordered_map<size_t, ConnectionEndpoint> id_to_endpoint;
    id_to_endpoint.reserve(step_list.size());

    auto pair_up = [&](const Connection& connection, const ConnectionEndpoint& ep) {
      const size_t id = connection.GetID();
      auto it = id_to_endpoint.find(id);
      if (it == id_to_endpoint.end()) {
        // first sighting: remember this endpoint
        id_to_endpoint.emplace(id, ep);
      } else {
        // second sighting: emit the pair and erase the remembered one
        endpoint_list.emplace_back(ep, it->second);
        id_to_endpoint.erase(it);
      }
    };

    auto visit_bc = [&](Step& step, IBoundaryCondition* bc) {
      if (auto* flow = dynamic_cast<FlowBoundaryCondition*>(bc)) {
        for (auto& stream : flow->streams_) {
          if (!stream->connection_.IsConnected()) continue;
          ConnectionEndpoint ep{};
          ep.conn = &stream->connection_;
          ep.step = &step;
          ep.stream = stream.get();
          pair_up(stream->connection_, ep);
        }
      } else if (auto* pressure = dynamic_cast<PressureBoundaryCondition*>(bc)) {
        if (pressure->connection_.IsConnected()) {
          ConnectionEndpoint ep{};
          ep.conn = &pressure->connection_;
          ep.step = &step;
          ep.p_boundary = pressure;
          pair_up(pressure->connection_, ep);
        }
      }
    };

    for (Step& step : step_list) {
      visit_bc(step, step.bottom_boundary_condition_.get());
      visit_bc(step, step.top_boundary_condition_.get());
    }

    // Note: any unpaired endpoints remain in id_to_endpoint need to be disconnected
    for (auto& [id, ep] : id_to_endpoint) {
      ep.conn->Disconnect();
    }

    return endpoint_list;
  }

  void DisconnectWithPartner(std::vector<Step>& step_list, Connection& conn) {
    if (!conn.IsConnected()) return;
    const size_t id = conn.GetID();
    conn.Disconnect();
    auto disconnect_match = [id](Connection& other, StreamData&, const std::string&) {
      if (other.IsConnected() && other.GetID() == id) other.Disconnect();
    };
    for (Step& step : step_list) {
      ForEachConnection(step.top_boundary_condition_.get(), disconnect_match);
      ForEachConnection(step.bottom_boundary_condition_.get(), disconnect_match);
    }
  }

  void DisconnectBoundaryCondition(std::vector<Step>& step_list, IBoundaryCondition* bc) {
    ForEachConnection(bc, [&](Connection& conn, StreamData&, const std::string&) {
      DisconnectWithPartner(step_list, conn);
    });
  }

} // end adsim namespace

CEREAL_CLASS_VERSION(adsim::Cycle, 2);
