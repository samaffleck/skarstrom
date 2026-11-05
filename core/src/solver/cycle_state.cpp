/*
 * Skarstrom - pressure and temperature swing adsorption simulator
 * Copyright (c) 2024-2026 Sam Affleck and contributors
 *
 * Licensed under the terms in the LICENSE file at the repository root.
 */
#include "../pch.h"
#include "adsim/solver/cycle_state.h"

#include <algorithm>
#include <cmath>
#include <sstream>
#include <stdexcept>
#include <unordered_map>

#include "adsim/objects/user_data.h"

namespace adsim {

  namespace {

    // The DAE layout of an initialised state: one CycleStateField per
    // equation, with the entries in the solver's banded order (cell-major:
    // for every cell, every equation active there). This is the order
    // Equations::PackBanded uses; CaptureCycleState checks the total against
    // the system size.
    std::vector<CycleStateField> BuildFields(const Equations& state, std::vector<size_t>* differential) {
      const auto& equations = state.GetEquations();
      // Component / adsorbent of the per-component equations.
      std::unordered_map<const Equation*, std::pair<std::string, size_t>> owner;
      for (const auto& [component, eq] : state.mole_fraction) owner[&eq] = {component, 0};
      for (const auto& [ads, comp_map] : state.ads_concentration)
        for (const auto& [component, eq] : comp_map) owner[&eq] = {component, ads};

      std::vector<CycleStateField> fields(equations.size());
      std::unordered_map<const Equation*, size_t> slot;
      size_t max_size = 0;
      for (size_t k = 0; k < equations.size(); ++k) {
        const Equation* eq = equations[k];
        CycleStateField& f = fields[k];
        f.name = eq->name;
        f.differential = eq->GetEquationType() == EquationType::Differential;
        if (auto it = owner.find(eq); it != owner.end()) {
          f.component = it->second.first;
          f.adsorbent = it->second.second;
        }
        slot[eq] = k;
        max_size = std::max(max_size, eq->Size());
      }

      const auto* void_mask = &state.domain_vectors.is_void_cell;
      size_t g = 0;
      for (size_t cell = 0; cell < max_size; ++cell) {
        for (const Equation* eq : equations) {
          if (!eq->IsActiveAt(cell, void_mask)) continue;
          CycleStateField& f = fields[slot[eq]];
          f.index.push_back(g);
          f.cell.push_back(cell);
          if (differential && eq->GetIDValue(cell) != 0.0) differential->push_back(g);
          ++g;
        }
      }
      return fields;
    }

    void BuildMetadata(UserData& data, CycleState& s) {
      s.components = data.fluid.components_;
      s.layers.clear();
      for (const auto& layer : data.bed.layers_.layers_)
        s.layers.push_back({layer.id_, layer.name_, layer.domain_.n_cells});
      s.steps.clear();
      for (const auto& step : data.cycle.step_list_)
        s.steps.push_back({step.id_, step.step_name_, step.step_time_.Get(), {}});
      for (const auto& h : data.cycle.Endpoints())
        s.steps[h.step_index].endpoint_ids.push_back(h.id);
    }

    // True when the deck has a connection or a buffer tank, i.e. something that
    // reads another endpoint's history or summary during the cycle.
    bool DeckReadsHistory(UserData& data) {
      for (const auto& h : data.cycle.Endpoints()) {
        if (h.stream && (h.stream->connection_.IsConnected() || h.stream->buffer_tank_.IsActive())) return true;
        if (h.pressure && (h.pressure->connection_.IsConnected() || h.pressure->buffer_tank_.IsActive())) return true;
      }
      return false;
    }

    std::string Join(const std::vector<std::string>& v) {
      std::string out;
      for (size_t i = 0; i < v.size(); ++i) out += (i ? ", " : "") + v[i];
      return out;
    }

  }  // namespace

  bool CycleState::HasHistory() const {
    for (const auto& e : endpoints)
      if (!e.history.time_series_.empty() || e.summary.is_valid) return true;
    return false;
  }

  std::vector<double> CycleState::Differential() const {
    std::vector<double> out;
    out.reserve(differential_index.size());
    for (size_t g : differential_index) out.push_back(y.at(g));
    return out;
  }

  void CycleState::SetDifferential(const std::vector<double>& values) {
    if (values.size() != differential_index.size())
      throw std::invalid_argument("CycleState::SetDifferential: expected " + std::to_string(differential_index.size())
                                  + " values, got " + std::to_string(values.size()));
    for (size_t k = 0; k < values.size(); ++k) y.at(differential_index[k]) = values[k];
  }

  const CycleStateField* CycleState::FindField(const std::string& name, const std::string& component,
                                               size_t adsorbent) const {
    for (const auto& f : fields) {
      // "Mole Fraction" with a component matches "<component> Mole Fraction".
      const bool name_ok = f.name == name || (!component.empty() && f.name == component + " " + name);
      if (!name_ok) continue;
      if (!component.empty() && f.component != component) continue;
      if (adsorbent != 0 && f.adsorbent != adsorbent) continue;
      return &f;
    }
    return nullptr;
  }

  std::vector<double> CycleState::FieldValues(const std::string& name, const std::string& component,
                                              size_t adsorbent) const {
    const CycleStateField* f = FindField(name, component, adsorbent);
    if (!f) throw std::out_of_range("CycleState: no field \"" + name + "\"" +
                                    (component.empty() ? std::string{} : " for " + component));
    std::vector<double> out;
    out.reserve(f->index.size());
    for (size_t g : f->index) out.push_back(y.at(g));
    return out;
  }

  CycleState CaptureCycleState(UserData& data, int cycles) {
    CycleState s;
    BuildMetadata(data, s);
    s.cycles = cycles;

    const int n = data.state.GetSystemSize();
    s.fields = BuildFields(data.state, &s.differential_index);
    size_t total = 0;
    for (const auto& f : s.fields) total += f.index.size();
    if (n <= 0 || total != static_cast<size_t>(n))
      throw std::logic_error("CaptureCycleState: the state layout (" + std::to_string(total)
                             + " entries) does not match the system size (" + std::to_string(n) + ")");
    s.y.resize(static_cast<size_t>(n));
    s.yp.resize(static_cast<size_t>(n));
    data.state.PackBanded(s.y.data(), s.yp.data());

    for (const auto& h : data.cycle.Endpoints()) {
      CycleStateEndpoint e;
      e.endpoint_id = h.id;
      e.step_id = h.step->id_;
      e.name = h.log_name;
      if (h.stream) {
        e.history = h.stream->stream_data_;
        e.summary = h.stream->cycle_summary_;
      } else if (h.pressure) {
        e.history = h.pressure->stream_data_;
        e.summary = h.pressure->cycle_summary_;
      }
      s.endpoints.push_back(std::move(e));
    }
    return s;
  }

  std::string CycleStateIncompatibility(const CycleState& state, UserData& data) {
    if (state.Empty()) return "the state is empty";
    if (state.components != data.fluid.components_)
      return "different components (state: " + Join(state.components) + "; deck: "
             + Join(data.fluid.components_) + ")";

    const auto& layers = data.bed.layers_.layers_;
    if (state.layers.size() != layers.size())
      return "different layers (state: " + std::to_string(state.layers.size()) + ", deck: "
             + std::to_string(layers.size()) + ")";
    for (size_t i = 0; i < layers.size(); ++i) {
      if (state.layers[i].id != layers[i].id_)
        return "different layers (layer " + std::to_string(i + 1) + " is layer[" + std::to_string(state.layers[i].id)
               + "] in the state, layer[" + std::to_string(layers[i].id_) + "] in the deck)";
      if (state.layers[i].cells != layers[i].domain_.n_cells)
        return "different cell counts in layer \"" + layers[i].name_ + "\" (state: "
               + std::to_string(state.layers[i].cells) + ", deck: " + std::to_string(layers[i].domain_.n_cells)
               + "); remapping a state onto another grid is not supported yet";
    }

    CycleState deck;
    BuildMetadata(data, deck);
    if (state.steps.size() != deck.steps.size())
      return "different step structure (state: " + std::to_string(state.steps.size()) + " steps, deck: "
             + std::to_string(deck.steps.size()) + ")";
    for (size_t i = 0; i < deck.steps.size(); ++i) {
      if (state.steps[i].id != deck.steps[i].id || state.steps[i].endpoint_ids != deck.steps[i].endpoint_ids)
        return "different step structure at step " + std::to_string(i + 1) + " (\"" + deck.steps[i].name + "\")";
    }

    if (state.HasHistory() && DeckReadsHistory(data)) {
      for (size_t i = 0; i < deck.steps.size(); ++i) {
        const double a = state.steps[i].time, b = deck.steps[i].time;
        if (std::abs(a - b) > 1e-12 * std::max(1.0, std::abs(b))) {
          std::ostringstream ss;
          ss << "step \"" << deck.steps[i].name << "\" changed its time (state: " << a << " s, deck: " << b
             << " s) and the deck's connections or buffer tanks read the state's step histories; "
             << "resampling histories onto new step times is not supported yet";
          return ss.str();
        }
      }
    }
    return {};
  }

  std::string CycleStateLayoutMismatch(const CycleState& state, const UserData& data) {
    std::vector<CycleStateField> fields = BuildFields(data.state, nullptr);
    if (fields.size() != state.fields.size())
      return "different state layout (" + std::to_string(state.fields.size()) + " fields in the state, "
             + std::to_string(fields.size()) + " in the deck: e.g. a changed energy balance or adsorbent)";
    for (size_t k = 0; k < fields.size(); ++k) {
      const auto& a = state.fields[k];
      const auto& b = fields[k];
      if (a.name != b.name || a.adsorbent != b.adsorbent || a.index != b.index)
        return "different state layout at field \"" + b.name + "\" (state: \"" + a.name + "\", "
               + std::to_string(a.index.size()) + " entries; deck: " + std::to_string(b.index.size()) + " entries)";
    }
    return {};
  }

  void ApplyCycleStateDAE(const CycleState& state, UserData& data) {
    if (static_cast<int>(state.y.size()) != data.state.GetSystemSize() || state.yp.size() != state.y.size())
      throw std::invalid_argument("ApplyCycleStateDAE: state size does not match the system");
    data.state.UnpackBanded(state.y.data(), state.yp.data());
  }

  void ApplyCycleStateHistories(const CycleState& state, UserData& data) {
    std::unordered_map<size_t, const CycleStateEndpoint*> by_id;
    for (const auto& e : state.endpoints) by_id[e.endpoint_id] = &e;
    for (const auto& h : data.cycle.Endpoints()) {
      StreamData* history = h.stream ? &h.stream->stream_data_ : h.pressure ? &h.pressure->stream_data_ : nullptr;
      CycleEndpointSummary* summary =
          h.stream ? &h.stream->cycle_summary_ : h.pressure ? &h.pressure->cycle_summary_ : nullptr;
      if (!history || !summary) continue;
      auto it = by_id.find(h.id);
      // Assign in place: connections and buffer tanks hold pointers to these
      // objects, which must stay where they are.
      if (it == by_id.end()) {
        history->Clear();
        summary->Clear();
      } else {
        *history = it->second->history;
        *summary = it->second->summary;
      }
    }
  }

}  // namespace adsim
