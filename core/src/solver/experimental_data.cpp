/*
 * Skarstrom - pressure and temperature swing adsorption simulator
 * Copyright (c) 2024-2026 Sam Affleck and contributors
 *
 * Licensed under the terms in the LICENSE file at the repository root.
 */
#include "../pch.h"
#include "adsim/solver/experimental_data.h"
#include "adsim/solver/parameterisation.h"
#include "adsim/math/unit_parser.h"
#include "adsim/simulation.h"  // BulkVariableNames, ComponentVariableNames (in enum order)

#include <charconv>
#include <algorithm>
#include <filesystem>
#include <cstdlib>
#include <cctype>
#include <string_view>
#include <cmath>
#include <limits>
#include <optional>
#include <stdexcept>

namespace adsim {

  static inline void chomp_cr(std::string& s) {
    if (!s.empty() && s.back() == '\r') s.pop_back();
  }

  static inline std::string trim(std::string s) {
    auto is_space = [](unsigned char c) { return std::isspace(c); };
    s.erase(s.begin(), std::find_if(s.begin(), s.end(), [&](unsigned char c) { return !is_space(c); }));
    s.erase(std::find_if(s.rbegin(), s.rend(), [&](unsigned char c) { return !is_space(c); }).base(), s.end());
    return s;
  }

static bool parse_double(std::string_view sv, double& out) {
    while (!sv.empty() && std::isspace(static_cast<unsigned char>(sv.front()))) sv.remove_prefix(1);
    while (!sv.empty() && std::isspace(static_cast<unsigned char>(sv.back())))  sv.remove_suffix(1);
    if (sv.empty()) return false;

    const char* b = sv.data();
    char* end;
    out = std::strtod(b, &end);

    return end == b + sv.size();
}

  static void split_row_2cols(const std::string& line, std::string& c0, std::string& c1) {
    // Allow ',', ';', or '\t'
    size_t pos = std::string::npos;
    for (char d : {',', ';', '\t'}) {
      size_t p = line.find(d);
      if (p != std::string::npos) { pos = p; break; }
    }
    if (pos == std::string::npos) {
      c0 = trim(line);
      c1.clear();
    }
    else {
      c0 = trim(line.substr(0, pos));
      c1 = trim(line.substr(pos + 1));
    }
  }

  bool ExperimentalData::LoadData(const std::filesystem::path& path, bool headers_included) {
    std::ifstream in(path);
    if (!in) return false;

    XYData d;
    d.name = path.stem().string();

    std::string line;
    bool first_line = true;

    while (std::getline(in, line)) {
      chomp_cr(line);
      if (line.empty()) continue;

      if (first_line && headers_included) {
        first_line = false;
        std::string h0, h1;
        split_row_2cols(line, h0, h1);
        if (!h0.empty()) d.x_label = h0;
        if (!h1.empty()) d.y_label = h1;
        continue;
      }
      first_line = false;

      std::string s0, s1;
      split_row_2cols(line, s0, s1);
      if (s0.empty() || s1.empty()) continue;

      double xv, yv;
      if (!parse_double(s0, xv)) continue;
      if (!parse_double(s1, yv)) continue;

      d.x.push_back(xv);
      d.y.push_back(yv);
    }

    if (d.x.empty() || d.y.empty() || d.x.size() != d.y.size()) 
      return false;

    data_.push_back(std::move(d));
    return true;
  }
  
  bool ExperimentalData::AddData(const XYData& data) {
    data_.push_back(data);
    data_.back().name = "Experiment " + std::to_string(data_.size());
    return true;
  }

  bool ExperimentalData::DeleteDataByIndex(size_t idx) {
    if (idx >= data_.size()) return false;
    data_.erase(data_.begin() + static_cast<std::ptrdiff_t>(idx));
    return true;
  }

  // -------------------------------------------------------------------------
  // Observable bindings
  // -------------------------------------------------------------------------

  namespace {

    struct KpiField { const char* name; double KPIRecord::*field; const char* unit; };
    // sensitivity::KPINames() with the SI unit of each.
    constexpr KpiField kKpiFields[] = {
      {"purity", &KPIRecord::purity, "mol/mol"},
      {"recovery", &KPIRecord::recovery, "mol/mol"},
      {"productivity_mass", &KPIRecord::productivity_mass, "kg/(kg*s)"},
      {"productivity_volume", &KPIRecord::productivity_volume, "kg/(m3*s)"},
      {"energy_consumption", &KPIRecord::energy_consumption, "J/kg"},
      {"production_rate", &KPIRecord::production_rate, "mol/s"},
      {"electrical_energy", &KPIRecord::electrical_energy, "J/kg"},
      {"thermal_energy", &KPIRecord::thermal_energy, "J/kg"},
      {"cooling_duty", &KPIRecord::cooling_duty, "J/kg"},
    };

    const KpiField* FindKpi(const std::string& name) {
      for (const auto& k : kKpiFields)
        if (name == k.name) return &k;
      return nullptr;
    }

    std::optional<size_t> IndexOf(const std::vector<std::string>& names, const std::string& name) {
      for (size_t i = 0; i < names.size(); ++i)
        if (names[i] == name) return i;
      return std::nullopt;
    }

    std::string Join(const std::vector<std::string>& v) {
      std::string out;
      for (const auto& s : v) out += (out.empty() ? "" : ", ") + s;
      return out;
    }

    std::vector<std::string> KpiNameList() {
      std::vector<std::string> out;
      for (const auto& k : kKpiFields) out.emplace_back(k.name);
      return out;
    }

    // "[Pa]" -> "Pa"
    std::string StripBrackets(std::string unit) {
      if (unit.size() >= 2 && unit.front() == '[' && unit.back() == ']') unit = unit.substr(1, unit.size() - 2);
      return unit;
    }

    // SI value = data value * factor + offset, for the dataset's y unit.
    UnitInfo YUnit(const ObservableBinding& b) {
      UnitInfo u;
      if (b.y_unit.empty()) return u;
      std::string err;
      if (!ParseUnitExpression(b.y_unit, u, &err))
        throw std::invalid_argument("dataset y unit '" + b.y_unit + "' does not parse: " + err);
      return u;
    }

  }  // namespace

  std::string ObservableSIUnit(const ObservableBinding& b) {
    switch (b.kind) {
      case ObservableKind::Bulk:
        if (auto i = IndexOf(BulkVariableNames(), b.variable))
          return StripBrackets(GetVariableUnits(static_cast<VariableType>(*i)));
        return {};
      case ObservableKind::Component:
        if (auto i = IndexOf(ComponentVariableNames(), b.variable))
          return StripBrackets(GetComponentVariableUnits(static_cast<ComponentVariableType>(*i)));
        return {};
      case ObservableKind::Kpi:
        if (const KpiField* k = FindKpi(b.variable)) return k->unit;
        return {};
      case ObservableKind::None:
        break;
    }
    return {};
  }

  bool NormaliseObservableBinding(ObservableBinding& b, const XYData& data,
                                  const std::vector<std::string>& components, std::string& error) {
    auto fail = [&](const std::string& message) { error = message; return false; };
    const int kind = static_cast<int>(b.kind);
    if (kind < 0 || kind > static_cast<int>(ObservableKind::Kpi)) return fail("unknown observable kind");
    const int loc = static_cast<int>(b.location);
    if (loc < 0 || loc > static_cast<int>(ObservableLocation::Position)) return fail("unknown observable location");

    if (!b.uncertainty.empty()) {
      if (b.uncertainty.size() != data.y.size())
        return fail("uncertainty has " + std::to_string(b.uncertainty.size()) + " values but the dataset has " +
                    std::to_string(data.y.size()) + " points");
      for (double s : b.uncertainty)
        if (!(std::isfinite(s) && s > 0.0)) return fail("every uncertainty must be a finite number > 0");
    }
    if (!(std::isfinite(b.weight) && b.weight > 0.0)) return fail("the weight must be a finite number > 0");
    if (!b.Bound()) return true;

    const bool has_component = std::find(components.begin(), components.end(), b.component) != components.end();
    switch (b.kind) {
      case ObservableKind::Bulk:
        if (!IndexOf(BulkVariableNames(), b.variable))
          return fail("unknown bulk variable '" + b.variable + "'; one of: " + Join(BulkVariableNames()));
        if (!b.component.empty()) return fail("'" + b.variable + "' is a bulk variable; it takes no component");
        break;
      case ObservableKind::Component:
        if (!IndexOf(ComponentVariableNames(), b.variable))
          return fail("unknown component variable '" + b.variable + "'; one of: " + Join(ComponentVariableNames()));
        if (!has_component) return fail("the deck has no component '" + b.component + "'");
        break;
      case ObservableKind::Kpi:
        if (!FindKpi(b.variable)) return fail("unknown KPI '" + b.variable + "'; one of: " + Join(KpiNameList()));
        if (!has_component) return fail("the deck has no component '" + b.component + "'");
        break;
      case ObservableKind::None:
        break;
    }

    if (b.kind == ObservableKind::Kpi) {
      // "s" is the default of a new binding: it means "not set" here.
      if (b.x_unit != "cycle" && !b.x_unit.empty() && b.x_unit != "-" && b.x_unit != "s")
        return fail("a KPI dataset's x is the cycle number; x unit '" + b.x_unit + "' is not allowed");
      b.x_unit = "cycle";
      b.x_scale = 1.0;
    } else {
      if (b.location == ObservableLocation::Position && !(std::isfinite(b.z) && b.z >= 0.0))
        return fail("the position z must be a finite number >= 0 [m]");
      if (b.cycle < -1) return fail("the cycle must be a 0-based index or -1 (the last logged cycle)");
      if (b.x_unit.empty()) b.x_unit = "s";
      UnitInfo xu;
      std::string err;
      if (!ParseUnitExpression(b.x_unit, xu, &err)) return fail("x unit '" + b.x_unit + "' does not parse: " + err);
      if (!(xu.dimension == Dimension{0, 0, 1, 0, 0}) || xu.offset != 0.0)
        return fail("x unit '" + b.x_unit + "' is not a unit of time");
      b.x_scale = xu.factor;
    }

    if (!b.y_unit.empty()) {
      UnitInfo yu, si;
      std::string err;
      if (!ParseUnitExpression(b.y_unit, yu, &err)) return fail("y unit '" + b.y_unit + "' does not parse: " + err);
      const std::string si_text = ObservableSIUnit(b);
      if (ParseUnitExpression(si_text, si)) {
        if (!(yu.dimension == si.dimension))
          return fail("y unit '" + b.y_unit + "' is a " + DimensionName(yu.dimension) + "; '" + b.variable +
                      "' is in " + si_text + " (" + DimensionName(si.dimension) + ")");
      } else if (b.y_unit != si_text) {
        return fail("'" + b.variable + "' is in " + si_text + "; give the y unit as that or leave it empty");
      }
    }
    return true;
  }

  SimulatedCounterpart SimulatedSeries(const XYData& data, const EquationLog& eq_log,
                                       const std::vector<KPIRecord>& kpis, const BedBoundarySources& sources) {
    const ObservableBinding& b = data.observable;
    if (!b.Bound()) throw std::invalid_argument("dataset '" + data.name + "' is not bound to an observable");
    SimulatedCounterpart out;
    out.si_unit = ObservableSIUnit(b);
    if (out.si_unit.empty())
      throw std::invalid_argument("dataset '" + data.name + "': unknown variable '" + b.variable + "'");
    const UnitInfo yu = YUnit(b);
    auto to_data = [&](double si) { return (si - yu.offset) / yu.factor; };

    if (b.kind == ObservableKind::Kpi) {
      const KpiField* k = FindKpi(b.variable);
      bool component_seen = false;
      for (const auto& r : kpis) {
        if (r.component != b.component) continue;
        component_seen = true;
        const double v = r.*(k->field);
        if (!std::isfinite(v)) continue;
        out.x.push_back(static_cast<double>(r.cycle - b.cycle_offset));
        out.y.push_back(to_data(v));
      }
      if (!component_seen)
        throw std::out_of_range("dataset '" + data.name + "': the run has no KPIs for '" + b.component +
                                "' (is it tracked?)");
      return out;
    }

    const ParamSpatialLoc where = b.location == ObservableLocation::Top      ? ParamSpatialLoc::Top
                                : b.location == ObservableLocation::Bottom   ? ParamSpatialLoc::Bottom
                                                                             : ParamSpatialLoc::Position;
    std::vector<double> t, v;
    bool ok = false;
    if (b.kind == ObservableKind::Bulk) {
      const auto i = IndexOf(BulkVariableNames(), b.variable);
      ok = ExtractBulkTimeSeries(eq_log, static_cast<VariableType>(*i), where, b.z, t, v, b.cycle, sources);
    } else {
      const auto i = IndexOf(ComponentVariableNames(), b.variable);
      ok = ExtractComponentTimeSeries(eq_log, static_cast<ComponentVariableType>(*i), b.component, where, b.z, t, v,
                                      b.cycle, sources);
    }
    if (!ok || t.empty())
      throw std::out_of_range("dataset '" + data.name + "': the run logged no '" + b.variable + "'" +
                              (b.cycle >= 0 ? " at cycle " + std::to_string(b.cycle) : std::string()));
    out.x.reserve(t.size());
    out.y.reserve(v.size());
    for (size_t i = 0; i < t.size() && i < v.size(); ++i) {
      out.x.push_back(t[i] / b.x_scale);
      out.y.push_back(to_data(v[i]));
    }
    return out;
  }

  DatasetComparison CompareWithSimulation(const XYData& data, const SimulatedCounterpart& sim) {
    DatasetComparison c;
    const ObservableBinding& b = data.observable;
    const bool by_cycle = b.kind == ObservableKind::Kpi;
    double sse = 0.0, wsse = 0.0, max_abs = 0.0;
    size_t n = 0;
    if (!sim.x.empty()) {
      const auto [lo_it, hi_it] = std::minmax_element(sim.x.begin(), sim.x.end());
      const double lo = *lo_it, hi = *hi_it;
      for (size_t i = 0; i < data.x.size() && i < data.y.size(); ++i) {
        double ys = std::numeric_limits<double>::quiet_NaN();
        if (by_cycle) {
          const double cyc = std::round(data.x[i]);
          for (size_t j = 0; j < sim.x.size(); ++j)
            if (sim.x[j] == cyc) { ys = sim.y[j]; break; }
        } else if (data.x[i] >= lo && data.x[i] <= hi) {
          ys = InterpLinear(sim.x, sim.y, data.x[i]);
        }
        if (!std::isfinite(ys) || !std::isfinite(data.y[i])) continue;
        const double r = ys - data.y[i];
        const double sigma = b.uncertainty.size() == data.y.size() ? b.uncertainty[i] : 1.0;
        sse += r * r;
        wsse += (r / sigma) * (r / sigma);
        max_abs = std::max(max_abs, std::abs(r));
        ++n;
      }
    }
    c.points = n;
    if (n == 0) return c;
    c.rmse = std::sqrt(sse / static_cast<double>(n));
    c.max_abs_error = max_abs;
    c.weighted_sse = b.weight * wsse;
    double ymin = std::numeric_limits<double>::infinity(), ymax = -ymin;
    for (double y : data.y)
      if (std::isfinite(y)) { ymin = std::min(ymin, y); ymax = std::max(ymax, y); }
    if (ymax > ymin) c.nrmse = c.rmse / (ymax - ymin);
    return c;
  }

} // end adsim namespace
