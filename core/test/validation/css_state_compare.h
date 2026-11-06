/*
 * Skarstrom - pressure and temperature swing adsorption simulator
 * Copyright (c) 2024-2026 Sam Affleck and contributors
 *
 * Licensed under the terms in the LICENSE file at the repository root.
 */
#pragma once

// Distances between two cyclic-steady-state results of the same deck: their
// end-of-cycle bed states, endpoint histories and KPIs. Used by the CSS
// verification tests (test_css_reference.cpp, test_css_start_independence.cpp)
// to compare a run that stopped at CSS with a reference.
//
// The bed-state distance uses the norm of the CSS state test
// (CyclicSteadyState::EndOfCycleStateChange), so a distance can be read
// against the stop rule's state tolerance directly: per field, the RMS over
// the interior cells of the difference, divided by a scale, and the largest
// over the fields. The scales are:
//   * adsorbed loadings: the RMS of the adsorbent's total loading in the
//     reference, floored at 1e-8 mol/kg (a trace component's round-off does
//     not dominate, nor an inert packing's);
//   * gas mole fractions: 1 (absolute);
//   * pressure and the fluid, solid and wall temperatures: their RMS in the
//     reference.
// The state test leaves pressure out; it is included here, as are all the
// differential fields of the state. Velocity (algebraic, recomputed from the
// state) is not compared.
//
// Histories: each endpoint's logged run is compared through its summary, the
// moles it passed and their average composition. Moles are relative to the
// largest |moles| of any endpoint in the reference; the average mole
// fractions are absolute and compared only where an endpoint carries at
// least 1e-3 of that largest flow (a closed or idle endpoint has no
// meaningful composition).
//
// KPIs: the largest absolute difference of purity and recovery [mol/mol]
// over the tracked components.

#include <algorithm>
#include <cmath>
#include <map>
#include <string>
#include <vector>

#include "adsim/simulation.h"
#include "adsim/solver/cycle_state.h"

namespace adsim::validation {

  struct CssDistance {
    double state = 0.0;    // bed state, in the CSS state-test norm
    double history = 0.0;  // endpoint histories
    double kpi = 0.0;      // purity / recovery, absolute
    std::string state_field;  // the field that set `state`
  };

  namespace detail {
    inline bool IsTemperatureOrPressure(const std::string& name) {
      return name == "Pressure" || name.find("Temperature") != std::string::npos;
    }
    inline bool IsLoading(const std::string& name) {
      return name.find("Adsorbed Concentration") != std::string::npos;
    }
    inline bool IsMoleFraction(const std::string& name) {
      return name.find("Mole Fraction") != std::string::npos;
    }
    // Interior entries of a field: every entry except the first and last cell
    // (the ghost cells).
    inline std::vector<double> Interior(const CycleState& s, const CycleStateField& f) {
      std::vector<double> out;
      if (f.index.empty()) return out;
      size_t lo = *std::min_element(f.cell.begin(), f.cell.end());
      size_t hi = *std::max_element(f.cell.begin(), f.cell.end());
      for (size_t k = 0; k < f.index.size(); ++k)
        if (f.cell[k] != lo && f.cell[k] != hi) out.push_back(s.y[f.index[k]]);
      return out;
    }
    inline double Rms(const std::vector<double>& v) {
      if (v.empty()) return 0.0;
      double s = 0.0;
      for (double x : v) s += x * x;
      return std::sqrt(s / static_cast<double>(v.size()));
    }
  }  // namespace detail

  // Bed-state distance of `a` from the reference `ref` (same deck, same grid).
  inline double StateDistance(const CycleState& a, const CycleState& ref, std::string* worst = nullptr) {
    using namespace detail;
    // Total loading per adsorbent in the reference, for the loading scale.
    std::map<size_t, std::vector<double>> total_loading;
    for (const auto& f : ref.fields) {
      if (!IsLoading(f.name)) continue;
      const auto v = Interior(ref, f);
      auto& t = total_loading[f.adsorbent];
      if (t.empty()) t.assign(v.size(), 0.0);
      for (size_t i = 0; i < v.size() && i < t.size(); ++i) t[i] += v[i];
    }
    double d = 0.0;
    for (const auto& f : ref.fields) {
      if (!f.differential || f.face) continue;
      double scale = 0.0;
      if (IsLoading(f.name)) scale = std::max(Rms(total_loading[f.adsorbent]), 1e-8);
      else if (IsMoleFraction(f.name)) scale = 1.0;
      else if (IsTemperatureOrPressure(f.name)) scale = Rms(Interior(ref, f));
      else continue;
      const CycleStateField* g = a.FindField(f.name, f.component, f.adsorbent);
      if (!g || !(scale > 0.0)) continue;
      const auto x = Interior(a, *g), r = Interior(ref, f);
      if (x.size() != r.size()) return std::numeric_limits<double>::infinity();
      std::vector<double> diff(x.size());
      for (size_t i = 0; i < x.size(); ++i) diff[i] = x[i] - r[i];
      const double e = Rms(diff) / scale;
      if (e > d) {
        d = e;
        if (worst) *worst = f.name;
      }
    }
    return d;
  }

  // Endpoint-history distance of `a` from `ref`.
  inline double HistoryDistance(const CycleState& a, const CycleState& ref) {
    double n_scale = 0.0;
    for (const auto& e : ref.endpoints)
      if (e.summary.is_valid) n_scale = std::max(n_scale, std::abs(e.summary.total_moles));
    if (!(n_scale > 0.0)) return 0.0;
    double d = 0.0;
    for (const auto& r : ref.endpoints) {
      if (!r.summary.is_valid) continue;
      const CycleStateEndpoint* x = nullptr;
      for (const auto& e : a.endpoints)
        if (e.endpoint_id == r.endpoint_id) x = &e;
      if (!x || !x->summary.is_valid) return std::numeric_limits<double>::infinity();
      d = std::max(d, std::abs(x->summary.total_moles - r.summary.total_moles) / n_scale);
      if (std::abs(r.summary.total_moles) < 1e-3 * n_scale) continue;
      for (const auto& [component, y] : r.summary.avg_mole_fraction) {
        const auto it = x->summary.avg_mole_fraction.find(component);
        if (it != x->summary.avg_mole_fraction.end()) d = std::max(d, std::abs(it->second - y));
      }
    }
    return d;
  }

  // KPI distance (purity and recovery of every tracked component).
  inline double KpiDistance(const std::vector<KPIRecord>& a, const std::vector<KPIRecord>& ref) {
    double d = 0.0;
    for (const auto& r : ref) {
      if (!r.valid) continue;
      for (const auto& k : a) {
        if (k.component != r.component) continue;
        if (!k.valid) return std::numeric_limits<double>::infinity();
        d = std::max({d, std::abs(k.purity - r.purity), std::abs(k.recovery - r.recovery)});
      }
    }
    return d;
  }

  inline CssDistance Distance(const CycleState& a, const std::vector<KPIRecord>& a_kpis, const CycleState& ref,
                              const std::vector<KPIRecord>& ref_kpis) {
    CssDistance d;
    d.state = StateDistance(a, ref, &d.state_field);
    d.history = HistoryDistance(a, ref);
    d.kpi = KpiDistance(a_kpis, ref_kpis);
    return d;
  }

}  // namespace adsim::validation
