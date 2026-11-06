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
#include <limits>
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
      // Cell fields have one entry per cell, face fields one fewer.
      for (size_t k = 0; k < equations.size(); ++k) fields[k].face = equations[k]->Size() < max_size;

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
      s.bottom_void = data.bed.bottom_void_.enabled_;
      s.top_void = data.bed.top_void_.enabled_;
      s.steps.clear();
      for (const auto& step : data.cycle.step_list_)
        s.steps.push_back({step.id_, step.step_name_, step.step_time_.Get(), {}});
      for (const auto& h : data.cycle.Endpoints())
        s.steps[h.step_index].endpoint_ids.push_back(h.id);
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
    }
    if (state.bottom_void != data.bed.bottom_void_.enabled_ || state.top_void != data.bed.top_void_.enabled_)
      return "different void regions (a void region was added or removed)";

    CycleState deck;
    BuildMetadata(data, deck);
    if (state.steps.size() != deck.steps.size())
      return "different step structure (state: " + std::to_string(state.steps.size()) + " steps, deck: "
             + std::to_string(deck.steps.size()) + ")";
    for (size_t i = 0; i < deck.steps.size(); ++i) {
      if (state.steps[i].id != deck.steps[i].id || state.steps[i].endpoint_ids != deck.steps[i].endpoint_ids)
        return "different step structure at step " + std::to_string(i + 1) + " (\"" + deck.steps[i].name + "\")";
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


  // ---- grid transfer ---------------------------------------------------------

  namespace {

    // The column split into segments along the cell array: the ghost cells,
    // the void cells (one cell each) and the layers (their cells). A segment
    // has the same index on both grids; only layer segments change size.
    struct Segment {
      size_t first = 0;  // first cell
      size_t cells = 0;
      bool layer = false;
    };

    std::vector<Segment> Segments(const std::vector<int>& layer_cells, bool bottom_void, bool top_void) {
      std::vector<Segment> segs;
      size_t c = 0;
      auto add = [&](size_t n, bool layer) { segs.push_back({c, n, layer}); c += n; };
      add(1, false);
      if (bottom_void) add(1, false);
      for (int n : layer_cells) add(static_cast<size_t>(n), true);
      if (top_void) add(1, false);
      add(1, false);
      return segs;
    }

    size_t TotalCells(const std::vector<Segment>& segs) { return segs.back().first + segs.back().cells; }

    // Lower edge of cell c in the segment coordinate (segment index + fraction).
    double CellEdge(const std::vector<Segment>& segs, size_t c) {
      for (size_t s = 0; s < segs.size(); ++s)
        if (c < segs[s].first + segs[s].cells)
          return static_cast<double>(s) + static_cast<double>(c - segs[s].first) / static_cast<double>(segs[s].cells);
      return static_cast<double>(segs.size());
    }

    // Overlap-weighted piecewise-constant remap of one segment: keeps the sum
    // of value x width (the cells of a layer are uniform). Exact copy when the
    // counts agree.
    void RemapSegment(const double* src, size_t n, double* dst, size_t m) {
      if (n == m) {
        std::copy(src, src + n, dst);
        return;
      }
      // On [0, n*m]: source cell i is [i m, (i+1) m), target cell j is [j n, (j+1) n).
      for (size_t j = 0; j < m; ++j) {
        const size_t lo = j * n, hi = (j + 1) * n;
        double acc = 0.0;
        for (size_t i = lo / m; i < n && i * m < hi; ++i) {
          const size_t a = std::max(lo, i * m), b = std::min(hi, (i + 1) * m);
          if (b > a) acc += src[i] * static_cast<double>(b - a);
        }
        dst[j] = acc / static_cast<double>(n);
      }
    }

    // A cell field as one value per cell (NaN where the field has no entry).
    std::vector<double> CellArray(const CycleStateField& f, const std::vector<double>& v, size_t cells) {
      std::vector<double> out(cells, std::numeric_limits<double>::quiet_NaN());
      for (size_t k = 0; k < f.index.size(); ++k)
        if (f.cell[k] < cells) out[f.cell[k]] = v[f.index[k]];
      return out;
    }

    // Cell arrays remapped segment by segment; a segment with any missing
    // source entry stays missing (NaN) on the target.
    std::vector<double> RemapCells(const std::vector<double>& src, const std::vector<Segment>& from,
                                   const std::vector<Segment>& to) {
      std::vector<double> out(TotalCells(to), std::numeric_limits<double>::quiet_NaN());
      for (size_t s = 0; s < to.size(); ++s) {
        const double* a = src.data() + from[s].first;
        if (std::any_of(a, a + from[s].cells, [](double x) { return std::isnan(x); })) continue;
        RemapSegment(a, from[s].cells, out.data() + to[s].first, to[s].cells);
      }
      return out;
    }

    // A face field (faces between cells c and c + 1) interpolated linearly in
    // the segment coordinate.
    std::vector<double> RemapFaces(const CycleStateField& f, const std::vector<double>& v,
                                   const std::vector<Segment>& from, const std::vector<Segment>& to) {
      std::vector<std::pair<double, double>> pts;
      for (size_t k = 0; k < f.index.size(); ++k) pts.emplace_back(CellEdge(from, f.cell[k] + 1), v[f.index[k]]);
      std::sort(pts.begin(), pts.end());
      const size_t faces = TotalCells(to) - 1;
      std::vector<double> out(faces, std::numeric_limits<double>::quiet_NaN());
      if (pts.empty()) return out;
      for (size_t j = 0; j < faces; ++j) {
        const double u = CellEdge(to, j + 1);
        auto it = std::lower_bound(pts.begin(), pts.end(), std::make_pair(u, -std::numeric_limits<double>::infinity()));
        if (it == pts.begin()) out[j] = it->second;
        else if (it == pts.end()) out[j] = pts.back().second;
        else {
          const auto& [u1, v1] = *it;
          const auto& [u0, v0] = *(it - 1);
          out[j] = u1 > u0 ? v0 + (v1 - v0) * (u - u0) / (u1 - u0) : v1;
        }
      }
      return out;
    }

    const CycleStateField* Match(const std::vector<CycleStateField>& fields, const CycleStateField& f) {
      for (const auto& g : fields)
        if (g.name == f.name && g.adsorbent == f.adsorbent) return &g;
      return nullptr;
    }

    void Write(const CycleStateField& f, const std::vector<double>& values, std::vector<double>& v) {
      for (size_t k = 0; k < f.index.size(); ++k) {
        const size_t c = f.cell[k];
        if (c < values.size() && !std::isnan(values[c])) v[f.index[k]] = values[c];
      }
    }

  }  // namespace

  bool CycleStateNeedsRemap(const CycleState& state, const UserData& data) {
    const auto& layers = data.bed.layers_.layers_;
    if (state.layers.size() != layers.size()) return false;
    for (size_t i = 0; i < layers.size(); ++i)
      if (state.layers[i].cells != layers[i].domain_.n_cells) return true;
    return false;
  }

  CycleState RemapCycleState(const CycleState& state, UserData& data) {
    const std::string why = CycleStateIncompatibility(state, data);
    if (!why.empty()) throw std::invalid_argument("RemapCycleState: the state does not fit the deck: " + why);

    // The deck's layout and initial conditions (the fallback for entries the
    // state cannot supply).
    data.Initialise();
    CycleState out;
    BuildMetadata(data, out);
    out.cycles = state.cycles;
    out.endpoints = state.endpoints;
    out.fields = BuildFields(data.state, &out.differential_index);
    const int n = data.state.GetSystemSize();
    out.y.assign(static_cast<size_t>(n), 0.0);
    out.yp.assign(static_cast<size_t>(n), 0.0);
    data.state.PackBanded(out.y.data(), out.yp.data());

    std::vector<int> from_cells, to_cells;
    for (const auto& l : state.layers) from_cells.push_back(l.cells);
    for (const auto& l : out.layers) to_cells.push_back(l.cells);
    const auto from = Segments(from_cells, state.bottom_void, state.top_void);
    const auto to = Segments(to_cells, out.bottom_void, out.top_void);
    const size_t from_n = TotalCells(from), to_n = TotalCells(to);

    // Every field (value and derivative) by the plain remap.
    std::vector<std::vector<double>> cells_y(out.fields.size());  // remapped cell arrays (y) by target field
    for (size_t k = 0; k < out.fields.size(); ++k) {
      const CycleStateField& t = out.fields[k];
      const CycleStateField* f = Match(state.fields, t);
      if (!f) continue;
      const bool face = t.face && f->face;
      if (face) {
        Write(t, RemapFaces(*f, state.y, from, to), out.y);
        Write(t, RemapFaces(*f, state.yp, from, to), out.yp);
      } else {
        cells_y[k] = RemapCells(CellArray(*f, state.y, from_n), from, to);
        Write(t, cells_y[k], out.y);
        Write(t, RemapCells(CellArray(*f, state.yp, from_n), from, to), out.yp);
      }
    }

    // The gas, on the layers whose cell count changed: remap the component
    // concentrations c_i = C y_i (conserved) instead of C (or P) and y_i.
    const CycleStateField* p_from = state.FindField("Pressure");
    const CycleStateField* p_to = out.FindField("Pressure");
    const CycleStateField* t_from = state.FindField("Fluid Temperature");
    const CycleStateField* c_from = state.FindField("Concentration");
    const CycleStateField* c_to = out.FindField("Concentration");
    if (p_from && p_to) {
      const auto P = CellArray(*p_from, state.y, from_n);
      const auto T = t_from ? CellArray(*t_from, state.y, from_n) : std::vector<double>(from_n, 1.0);
      const auto C = c_from ? CellArray(*c_from, state.y, from_n) : std::vector<double>{};
      std::vector<double> total(from_n), heat(from_n);
      for (size_t c = 0; c < from_n; ++c) {
        total[c] = c_from ? C[c] : P[c] / (R * T[c]);
        heat[c] = total[c] * T[c];  // C T (= P / R when C is derived)
      }
      // The fluid temperature, weighted by the gas concentration: keeps the
      // gas's sensible heat (C T) and makes P = C R T the plain cell average
      // of P, which is monotone. Remapping C and T separately would not: their
      // product overshoots where the temperature varies.
      const auto total2 = RemapCells(total, from, to);
      const auto heat2 = RemapCells(heat, from, to);
      std::vector<double> T2(to_n, std::numeric_limits<double>::quiet_NaN());
      for (size_t c = 0; c < to_n; ++c)
        if (total2[c] > 0.0) T2[c] = heat2[c] / total2[c];
      const CycleStateField* t_to = out.FindField("Fluid Temperature");

      std::vector<std::vector<double>> conc;  // remapped c_i per component
      std::vector<const CycleStateField*> y_to;
      for (const auto& comp : state.components) {
        const CycleStateField* yf = state.FindField("Mole Fraction", comp);
        const CycleStateField* yt = out.FindField("Mole Fraction", comp);
        if (!yf || !yt) continue;
        auto y = CellArray(*yf, state.y, from_n);
        for (size_t c = 0; c < from_n; ++c) y[c] *= total[c];
        conc.push_back(RemapCells(y, from, to));
        y_to.push_back(yt);
      }
      for (size_t s = 0; s < to.size(); ++s) {
        if (!to[s].layer || to[s].cells == from[s].cells) continue;
        for (size_t c = to[s].first; c < to[s].first + to[s].cells; ++c) {
          double sum = 0.0;
          for (const auto& ci : conc) sum += ci[c];
          if (!(sum > 0.0) || std::isnan(T2[c])) continue;
          std::vector<double> cell_values(to_n, std::numeric_limits<double>::quiet_NaN());
          for (size_t i = 0; i < conc.size(); ++i) {
            cell_values[c] = conc[i][c] / sum;
            Write(*y_to[i], cell_values, out.y);
          }
          if (t_from && t_to) {
            cell_values[c] = T2[c];
            Write(*t_to, cell_values, out.y);
          }
          cell_values[c] = sum * R * T2[c];
          Write(*p_to, cell_values, out.y);
          if (c_to) {
            cell_values[c] = sum;
            Write(*c_to, cell_values, out.y);
          }
        }
      }
    }
    return out;
  }


  // ---- step-time changes and the combined preparation -------------------------

  namespace {
    bool TimesDiffer(double a, double b) { return std::abs(a - b) > 1e-12 * std::max(1.0, std::abs(b)); }

    std::string FormatSeconds(double t) {
      std::ostringstream ss;
      ss << t << " s";
      return ss.str();
    }
  }  // namespace

  bool CycleStateNeedsResample(const CycleState& state, const UserData& data) {
    const auto& steps = data.cycle.step_list_;
    if (state.steps.size() != steps.size()) return false;
    for (size_t i = 0; i < steps.size(); ++i)
      if (TimesDiffer(state.steps[i].time, steps[i].step_time_.Get())) return true;
    return false;
  }

  CycleState ResampleCycleStateHistories(const CycleState& state, const UserData& data, std::string* note) {
    CycleState out = state;
    const auto& steps = data.cycle.step_list_;
    if (state.steps.size() != steps.size())
      throw std::invalid_argument("ResampleCycleStateHistories: different step structure");
    std::vector<std::string> changed;
    for (size_t i = 0; i < steps.size(); ++i) {
      const double t_old = state.steps[i].time, t_new = steps[i].step_time_.Get();
      if (!TimesDiffer(t_old, t_new)) continue;
      changed.push_back("\"" + steps[i].step_name_ + "\" " + FormatSeconds(t_old) + " -> " + FormatSeconds(t_new));
      out.steps[i].time = t_new;
      const bool scalable = t_old > 0.0 && t_new > 0.0;
      const double ratio = scalable ? t_new / t_old : 0.0;
      for (auto& e : out.endpoints) {
        if (e.step_id != state.steps[i].id) continue;
        if (!scalable) {
          // A zero-length step on either side: nothing to stretch; the
          // consumers fall back to their first-cycle values.
          e.history.Clear();
          e.summary.Clear();
          continue;
        }
        // The step's logged time series, stretched onto the new duration
        // (same shape, step-relative time), and the moles of its summary
        // scaled with the duration (same average flow).
        for (double& t : e.history.time_series_) t *= ratio;
        e.summary.total_moles *= ratio;
      }
    }
    if (note && !changed.empty()) {
      std::string list;
      for (size_t k = 0; k < changed.size(); ++k) list += (k ? ", " : "") + changed[k];
      *note = "histories stretched onto the new step times (" + list + ")";
    }
    return out;
  }

  CycleState PrepareCycleState(const CycleState& state, UserData& data, std::string* note) {
    const std::string why = CycleStateIncompatibility(state, data);
    if (!why.empty()) throw std::invalid_argument("the state does not fit the deck: " + why);
    std::vector<std::string> notes;
    CycleState out = state;
    if (CycleStateNeedsResample(out, data)) {
      std::string n;
      out = ResampleCycleStateHistories(out, data, &n);
      notes.push_back(n);
    }
    if (CycleStateNeedsRemap(out, data)) {
      std::string grid;
      for (size_t i = 0; i < out.layers.size(); ++i)
        grid += (i ? ", " : "") + std::to_string(out.layers[i].cells) + " -> "
                + std::to_string(data.bed.layers_.layers_[i].domain_.n_cells);
      out = RemapCycleState(out, data);
      notes.push_back("remapped onto the deck's grid (cells per layer " + grid + ")");
    }
    if (note) {
      note->clear();
      for (size_t k = 0; k < notes.size(); ++k) *note += (k ? "; " : "") + notes[k];
    }
    return out;
  }

}  // namespace adsim
