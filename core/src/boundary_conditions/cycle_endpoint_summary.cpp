/*
 * Skarstrom - pressure and temperature swing adsorption simulator
 * Copyright (c) 2024-2026 Sam Affleck and contributors
 *
 * Licensed under the terms in the LICENSE file at the repository root.
 */
#include "../pch.h"
#include "adsim/boundary_conditions/cycle_endpoint_summary.h"
#include "adsim/boundary_conditions/boundary_stream.h"
#include "adsim/math/general_functions.h"

#include <cmath>


namespace adsim {

  void CycleEndpointSummary::Clear() {
    is_valid = false;
    total_moles = 0.0;
    avg_temperature = 0.0;
    avg_mole_fraction.clear();
    source_cycle_index = -1;
  }

  CycleEndpointSummary ComputeSummary(const StreamData& s, int cycle_index) {
    CycleEndpointSummary summary;
    summary.source_cycle_index = cycle_index;

    const auto& t = s.time_series_;
    const auto& n_dot = s.molar_flow_;
    const auto& T = s.temperature_boundary_;

    if (t.size() < 2 || n_dot.size() < 2 || T.size() < 2) {
      return summary; // is_valid stays false; consumer falls back to first-cycle assumption
    }

    summary.total_moles = IntegrateVariableInTime(n_dot, t);

    // Temperature average uses |n_dot| as weight so back-flow can't blow up the mean.
    std::vector<double> abs_n_dot(n_dot.size());
    std::vector<double> T_weighted(n_dot.size());
    for (size_t i = 0; i < n_dot.size(); ++i) {
      abs_n_dot[i] = std::abs(n_dot[i]);
      T_weighted[i] = T[i] * abs_n_dot[i];
    }
    const double abs_total = IntegrateVariableInTime(abs_n_dot, t);
    if (abs_total > 1e-30) {
      summary.avg_temperature = IntegrateVariableInTime(T_weighted, t) / abs_total;
    } else if (!T.empty()) {
      summary.avg_temperature = T.back();
    }

    // y_i_avg = (∫ n_i dt) / (∫ n_dot dt)
    for (const auto& [comp, comp_flow] : s.component_molar_flow_) {
      double y_avg = 0.0;
      if (std::abs(summary.total_moles) > 1e-30) {
        y_avg = IntegrateVariableInTime(comp_flow, t) / summary.total_moles;
      } else {
        const auto it_y = s.mole_fraction_boundary_.find(comp);
        if (it_y != s.mole_fraction_boundary_.end() && !it_y->second.empty()) {
          y_avg = it_y->second.back();
        }
      }
      if (!std::isfinite(y_avg) || y_avg < 0.0) y_avg = 0.0;
      summary.avg_mole_fraction[comp] = y_avg;
    }

    double sum_y = 0.0;
    for (const auto& [c, y] : summary.avg_mole_fraction) { (void)c; sum_y += y; }
    if (sum_y > 0.0) {
      for (auto& [c, y] : summary.avg_mole_fraction) { (void)c; y /= sum_y; }
    }

    summary.is_valid = true;
    return summary;
  }

} // end adsim namespace
