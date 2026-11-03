/*
 * Skarstrom - pressure and temperature swing adsorption simulator
 * Copyright (c) 2024-2026 Sam Affleck and contributors
 *
 * Licensed under the terms in the LICENSE file at the repository root.
 */
#include "../pch.h"
#include "adsim/solver/bed_series.h"

#include <algorithm>
#include <cmath>
#include <limits>
#include <optional>

#include "adsim/boundary_conditions/step.h"

namespace adsim {

  namespace {

    constexpr double kNaN = std::numeric_limits<double>::quiet_NaN();
    // Below this total molar flow [mol/s] the streams at one end are averaged
    // instead of mixed by their flows (the stream logs use the same floor).
    constexpr double kTinyFlow = 1e-12;

    bool IsPorousOnly(const BedVariable& v) {
      return v.is_component ? IsPorousOnlyComponentVariable(v.comp) : IsPorousOnlyVariable(v.bulk);
    }

    // The variables whose bed-end value is logged per endpoint.
    bool HasEndpointLog(const BedVariable& v) {
      if (v.is_component) return v.comp == ComponentVariableType::MOLE_FRACTION;
      return v.bulk == VariableType::PRESSURE || v.bulk == VariableType::TEMPERATURE;
    }

    // Linear interpolation of (x, y) at xq, x ascending, clamped to the ends.
    double Interp(const double* x, const double* y, size_t n, double xq) {
      if (n == 0) return kNaN;
      if (xq <= x[0]) return y[0];
      if (xq >= x[n - 1]) return y[n - 1];
      const size_t k = static_cast<size_t>(std::upper_bound(x, x + n, xq) - x);  // 1..n-1
      const double dx = x[k] - x[k - 1];
      if (!(dx > 0.0)) return y[k];
      const double w = (xq - x[k - 1]) / dx;
      return y[k - 1] + w * (y[k] - y[k - 1]);
    }

    double Lerp(double x0, double y0, double x1, double y1, double xq) {
      if (!(x1 != x0)) return y1;
      return y0 + (xq - x0) / (x1 - x0) * (y1 - y0);
    }

    // The logged endpoints at one bed end over one step.
    struct EndStreams {
      std::vector<const StreamData*> streams;
      bool by_index = false;  // the stream rows are the equation-log rows

      // `v` (a vector of one of the streams) at equation-log row t, whose time
      // from the step start is t_rel.
      double At(const StreamData& sd, const std::vector<double>& v, size_t t, double t_rel) const {
        const auto& ts = sd.time_series_;
        if (v.empty() || v.size() != ts.size()) return kNaN;
        if (by_index) return t < v.size() ? v[t] : kNaN;
        return Interp(ts.data(), v.data(), ts.size(), t_rel);
      }
    };

    // Every endpoint the step has at `location`, with its log for this cycle
    // and step. Empty when there is none (a closed end), when the sources are
    // not given, or when one of them was not logged.
    EndStreams FindEndStreams(const BedBoundarySources& sources, int cycle, int step, BoundaryLocation location,
                              const std::vector<double>& eq_times) {
      EndStreams out;
      if (!sources.boundary_log || !sources.steps || step < 0 ||
          step >= static_cast<int>(sources.steps->size()))
        return out;
      for (const auto& ep : (*sources.steps)[static_cast<size_t>(step)].GetEndpoints()) {
        if (ep.location != location) continue;
        const StreamData* sd = sources.boundary_log->GetStreamData(cycle, step, ep.name);
        if (!sd || sd->time_series_.empty()) return {};
        if (!out.streams.empty() && sd->time_series_.size() != out.streams.front()->time_series_.size()) return {};
        out.streams.push_back(sd);
      }
      if (out.streams.empty()) return out;
      // The solver logs a stream row with every equation-log row (time from
      // the step start); read them by index when they line up.
      const auto& ts = out.streams.front()->time_series_;
      if (ts.size() == eq_times.size()) {
        auto close = [](double a, double b) { return std::abs(a - b) <= 1e-9 * std::max(1.0, std::abs(b)); };
        out.by_index = close(ts.front(), eq_times.front() - eq_times.front()) &&
                       close(ts.back(), eq_times.back() - eq_times.front());
      }
      return out;
    }

    // Bed-end value of `var` from the endpoint logs at row t, or nullopt.
    std::optional<double> EndpointValue(const BedVariable& var, const EndStreams& e, size_t t, double t_rel) {
      if (e.streams.empty()) return std::nullopt;
      auto values_of = [&](const StreamData& sd) -> const std::vector<double>* {
        if (!var.is_component) {
          return var.bulk == VariableType::PRESSURE ? &sd.pressure_boundary_ : &sd.temperature_boundary_;
        }
        auto it = sd.mole_fraction_boundary_.find(var.component);
        return it == sd.mole_fraction_boundary_.end() ? nullptr : &it->second;
      };
      // Every stream at one end sees the same face pressure.
      if (!var.is_component && var.bulk == VariableType::PRESSURE) {
        const double p = e.At(*e.streams.front(), e.streams.front()->pressure_boundary_, t, t_rel);
        return std::isfinite(p) ? std::optional<double>(p) : std::nullopt;
      }
      double weighted = 0.0, flow = 0.0, plain = 0.0;
      for (const StreamData* sd : e.streams) {
        const std::vector<double>* v = values_of(*sd);
        if (!v) return std::nullopt;
        const double x = e.At(*sd, *v, t, t_rel);
        if (!std::isfinite(x)) return std::nullopt;
        const double n = std::abs(e.At(*sd, sd->molar_flow_, t, t_rel));
        plain += x;
        if (std::isfinite(n)) {
          weighted += n * x;
          flow += n;
        }
      }
      if (e.streams.size() == 1) return plain;
      if (flow > kTinyFlow) return weighted / flow;
      return plain / static_cast<double>(e.streams.size());
    }

    // Bed-end face value of a cell row (ghost cells included).
    double FaceValue(const std::vector<double>& row, ParamSpatialLoc where, const BedVariable& var,
                     const DomainVectors& domain) {
      const size_t N = row.size();
      if (N < 2) return N == 1 ? row[0] : kNaN;
      // A mole-fraction ghost holds the face value itself: the inlet
      // composition, or a copy of the end cell at an outlet or closed end.
      if (var.is_component && var.comp == ComponentVariableType::MOLE_FRACTION)
        return where == ParamSpatialLoc::Bottom ? row[0] : row[N - 1];
      if (IsPorousOnly(var) && !domain.bed_cell_indices.empty()) {
        const int idx = where == ParamSpatialLoc::Bottom ? domain.bed_cell_indices.front()
                                                         : domain.bed_cell_indices.back();
        if (idx >= 0 && static_cast<size_t>(idx) < N) return row[static_cast<size_t>(idx)];
      }
      return where == ParamSpatialLoc::Bottom ? 0.5 * (row[0] + row[1]) : 0.5 * (row[N - 1] + row[N - 2]);
    }

  }  // namespace

  std::vector<double> CellCentrePositions(const DomainVectors& domain) {
    const auto& z = domain.positions;
    std::vector<double> out;
    if (z.size() < 2) return out;
    out.reserve(z.size() - 1);
    for (size_t i = 1; i < z.size(); ++i) out.push_back(0.5 * (z[i - 1] + z[i]));
    return out;
  }

  bool IsFaceRow(const DomainVectors& domain, size_t row_size) {
    return !domain.positions.empty() && row_size == domain.positions.size();
  }

  VariableLogEntry GetBedVariable(const EquationLog& eq_log, int cycle_index, int step_index,
                                  const BedVariable& var) {
    return var.is_component ? eq_log.GetComponentVariable(cycle_index, step_index, var.comp, var.component)
                            : eq_log.GetVariable(cycle_index, step_index, var.bulk);
  }

  bool ExtractBedTimeSeries(const EquationLog& eq_log, const BedVariable& var, ParamSpatialLoc where,
                            double z_pos, std::vector<double>& out_times, std::vector<double>& out_values,
                            int cycle_index, const BedBoundarySources& sources) {
    out_times.clear();
    out_values.clear();
    const auto& sim_log = eq_log.GetSimulationLog();
    if (sim_log.empty()) return false;
    if (cycle_index >= 0 && sim_log.find(cycle_index) == sim_log.end()) return false;
    const int cycle = cycle_index >= 0 ? cycle_index : sim_log.rbegin()->first;

    const DomainVectors& domain = eq_log.GetDomain();
    const std::vector<double>& faces = domain.positions;
    if (faces.empty()) return false;
    const std::vector<double> centres = CellCentrePositions(domain);
    const double z_bottom = faces.front();
    const double z_top = faces.back();
    const double zq = std::clamp(z_pos, z_bottom, z_top);

    const auto& cycle_log = sim_log.at(cycle);
    // Logged times are absolute solver times: re-base them on the cycle start
    // (the first log of the cycle's first step).
    double cycle_start = 0.0;
    for (const auto& [step_idx, step_log] : cycle_log) {
      (void)step_idx;
      if (!step_log.step_data.times.empty()) {
        cycle_start = step_log.step_data.times.front();
        break;
      }
    }

    for (const auto& [step_idx, step_log] : cycle_log) {
      (void)step_log;
      const VariableLogEntry vle = GetBedVariable(eq_log, cycle, step_idx, var);
      if (vle.times.empty()) continue;

      // The endpoint logs at each end, looked up once per step and only when
      // a row needs them.
      std::optional<EndStreams> ends[2];
      auto end_streams = [&](ParamSpatialLoc end) -> const EndStreams& {
        const int k = end == ParamSpatialLoc::Bottom ? 0 : 1;
        if (!ends[k])
          ends[k] = FindEndStreams(sources, cycle, step_idx,
                                   end == ParamSpatialLoc::Bottom ? BoundaryLocation::BOTTOM : BoundaryLocation::TOP,
                                   vle.times);
        return *ends[k];
      };
      auto boundary_value = [&](const std::vector<double>& row, size_t t, ParamSpatialLoc end) {
        if (HasEndpointLog(var)) {
          if (auto v = EndpointValue(var, end_streams(end), t, vle.times[t] - vle.times.front())) return *v;
        }
        return FaceValue(row, end, var, domain);
      };

      for (size_t t = 0; t < vle.times.size(); ++t) {
        out_times.push_back(vle.times[t] - cycle_start);
        if (t >= vle.values.size()) {
          out_values.push_back(kNaN);
          continue;
        }
        const std::vector<double>& row = vle.values[t];
        double value = kNaN;
        if (IsFaceRow(domain, row.size())) {
          switch (where) {
            case ParamSpatialLoc::Bottom: value = row.front(); break;
            case ParamSpatialLoc::Top: value = row.back(); break;
            case ParamSpatialLoc::Position: value = Interp(faces.data(), row.data(), row.size(), zq); break;
          }
        } else if (row.size() == centres.size() + 2) {
          const size_t n = centres.size();  // interior cells, row[1..n]
          if (n == 0) {  // ghost cells only (no bed cells): nothing in between
            out_values.push_back(FaceValue(row, where == ParamSpatialLoc::Top ? where : ParamSpatialLoc::Bottom,
                                           var, domain));
            continue;
          }
          switch (where) {
            case ParamSpatialLoc::Bottom: value = boundary_value(row, t, ParamSpatialLoc::Bottom); break;
            case ParamSpatialLoc::Top: value = boundary_value(row, t, ParamSpatialLoc::Top); break;
            case ParamSpatialLoc::Position:
              if (zq <= centres.front())
                value = Lerp(z_bottom, boundary_value(row, t, ParamSpatialLoc::Bottom), centres.front(), row[1], zq);
              else if (zq >= centres.back())
                value = Lerp(centres.back(), row[n], z_top, boundary_value(row, t, ParamSpatialLoc::Top), zq);
              else
                value = Interp(centres.data(), row.data() + 1, n, zq);
              break;
          }
        }
        out_values.push_back(value);
      }
    }
    return !out_times.empty();
  }

}  // namespace adsim
