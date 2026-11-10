/*
 * Skarstrom - pressure and temperature swing adsorption simulator
 * Copyright (c) 2024-2026 Sam Affleck and contributors
 *
 * Licensed under the terms in the LICENSE file at the repository root.
 */
#include "../pch.h"
#include "adsim/solver/step_end_pressure.h"

#include <algorithm>
#include <cctype>
#include <cmath>
#include <stdexcept>

#include "adsim/boundary_conditions/step.h"
#include "adsim/solver/bed_series.h"
#include "adsim/solver/log.h"

namespace adsim {

  std::vector<StepEndPressure> ComputeStepEndPressures(const std::vector<Step>& steps, const EquationLog& eq_log,
                                                       const BoundaryLog& b_log, int cycle) {
    std::vector<StepEndPressure> out;
    const auto& sim_log = eq_log.GetSimulationLog();
    const auto cycle_it = sim_log.find(cycle);
    if (cycle_it == sim_log.end()) return out;

    // The bed-end series of the whole cycle, as Simulation::Series gives
    // them: the step boundaries are found below from the rows of each step.
    const BedVariable pressure = BedVariable::Bulk(VariableType::PRESSURE);
    const BedBoundarySources sources{&b_log, &steps};
    std::vector<double> t_bottom, p_bottom, t_top, p_top;
    ExtractBedTimeSeries(eq_log, pressure, ParamSpatialLoc::Bottom, 0.0, t_bottom, p_bottom, cycle, sources);
    ExtractBedTimeSeries(eq_log, pressure, ParamSpatialLoc::Top, 0.0, t_top, p_top, cycle, sources);

    const DomainVectors& domain = eq_log.GetDomain();
    out.reserve(steps.size());
    for (size_t s = 0; s < steps.size(); ++s) {
      StepEndPressure r;
      r.cycle = cycle;
      r.step = static_cast<int>(s);
      r.step_id = steps[s].id_;
      r.step_name = steps[s].step_name_;
      out.push_back(std::move(r));
    }

    // ExtractBedTimeSeries walks the logged steps in order and appends every
    // row of each; `row` is the running row count.
    size_t row = 0;
    for (const auto& [step_idx, step_log] : cycle_it->second) {
      (void)step_log;
      const VariableLogEntry vle = GetBedVariable(eq_log, cycle, step_idx, pressure);
      if (vle.times.empty()) continue;
      row += vle.times.size();
      if (step_idx < 0 || static_cast<size_t>(step_idx) >= out.size()) continue;
      StepEndPressure& r = out[static_cast<size_t>(step_idx)];
      const size_t last = row - 1;
      if (last < t_bottom.size()) {
        r.time = t_bottom[last];
        r.bottom = p_bottom[last];
      }
      if (last < p_top.size()) r.top = p_top[last];

      if (vle.values.size() < vle.times.size()) continue;
      const std::vector<double>& cells = vle.values.back();
      std::vector<int> packed = domain.bed_cell_indices;
      if (packed.empty())
        for (int i = 1; i + 1 < static_cast<int>(cells.size()); ++i) packed.push_back(i);
      double lo = std::numeric_limits<double>::infinity(), hi = -lo, sum = 0.0, width = 0.0;
      for (int i : packed) {
        if (i < 0 || static_cast<size_t>(i) >= cells.size()) continue;
        const double p = cells[static_cast<size_t>(i)];
        const double w = static_cast<size_t>(i) < domain.cell_width.size() ? domain.cell_width[static_cast<size_t>(i)]
                                                                            : 1.0;
        lo = std::min(lo, p);
        hi = std::max(hi, p);
        sum += w * p;
        width += w;
      }
      if (width > 0.0) {
        r.bed_min = lo;
        r.bed_max = hi;
        r.bed_mean = sum / width;
      }
    }
    return out;
  }

  const std::vector<std::string>& StepOutputNames() {
    static const std::vector<std::string> kNames = {"end_pressure_bottom", "end_pressure_top", "end_pressure_min",
                                                    "end_pressure_mean", "end_pressure_max"};
    return kNames;
  }

  double StepOutputValue(const StepEndPressure& r, std::string_view name) {
    if (name == "end_pressure_bottom") return r.bottom;
    if (name == "end_pressure_top") return r.top;
    if (name == "end_pressure_min") return r.bed_min;
    if (name == "end_pressure_mean") return r.bed_mean;
    if (name == "end_pressure_max") return r.bed_max;
    throw std::invalid_argument("unknown step output '" + std::string(name) + "'");
  }

  bool IsStepReference(std::string_view subject) {
    return subject.size() > 6 && subject.substr(0, 5) == "step[" && subject.back() == ']';
  }

  int ResolveStepReference(const std::vector<Step>& steps, std::string_view subject) {
    if (!IsStepReference(subject))
      throw std::invalid_argument("'" + std::string(subject) + "' is not a step reference step[<ID or name>]");
    const std::string key(subject.substr(5, subject.size() - 6));
    const bool number = !key.empty() && std::all_of(key.begin(), key.end(), [](unsigned char c) {
      return std::isdigit(c) != 0;
    });
    if (number && key.size() <= 18) {
      const unsigned long long id = std::stoull(key);
      for (size_t s = 0; s < steps.size(); ++s)
        if (steps[s].id_ == id) return static_cast<int>(s);
    }
    for (size_t s = 0; s < steps.size(); ++s)
      if (steps[s].step_name_ == key) return static_cast<int>(s);
    throw std::invalid_argument("'" + std::string(subject) + "' names no step of the cycle");
  }

}  // namespace adsim
