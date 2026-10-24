/*
 * Skarstrom - pressure and temperature swing adsorption simulator
 * Copyright (c) 2024-2026 Sam Affleck and contributors
 *
 * Licensed under the terms in the LICENSE file at the repository root.
 */
#include "../pch.h"
#include "adsim/solver/operating_limits.h"
#include "adsim/objects/user_data.h"
#include "adsim/objects/wall_design.h"
#include "adsim/solver/log.h"

#include <cmath>
#include <iomanip>
#include <sstream>

namespace adsim {

  namespace {

    // Below these margins a limit is not counted as passed. Logged values
    // carry solver noise around a boundary value, and decks often take 1 bar
    // as atmospheric (1.3 kPa below kAtmosphericPressure).
    constexpr double kPressureMargin = 0.02 * kAtmosphericPressure;  // [Pa], 2 kPa
    constexpr double kTemperatureMargin = 0.01;    // [K]

    std::string Num(double v, int precision = 3) {
      std::ostringstream ss;
      ss << std::setprecision(precision) << v;
      return ss.str();
    }

    std::string BarGauge(double pascal_gauge) { return Num(pascal_gauge / 1e5) + " bar(g)"; }
    std::string BarAbs(double pascal) { return Num(pascal / 1e5) + " bar(a)"; }
    std::string Celsius(double kelvin) { return Num(kelvin - 273.15, 4) + " °C"; }

    std::string InStep(const StepOperatingSummary& s, int cycle) {
      return "in step \"" + s.step_name + "\" of cycle " + std::to_string(cycle + 1);
    }

  } // namespace

  OperatingLimitsReport CheckOperatingLimits(const UserData& data, const EquationLog& log, int cycle) {
    OperatingLimitsReport report;
    const auto& sim_log = log.GetSimulationLog();
    if (sim_log.empty()) return report;
    const auto cycle_it = cycle < 0 ? std::prev(sim_log.end()) : sim_log.find(cycle);
    if (cycle_it == sim_log.end()) return report;
    report.cycle = cycle_it->first;

    const auto& dom = log.GetDomain();
    const auto& bed_cells = dom.bed_cell_indices;
    // The column's cells: packed and void, no ghosts.
    std::vector<int> column_cells = bed_cells;
    if (dom.bottom_void_cell_index > 0) column_cells.push_back(dom.bottom_void_cell_index);
    if (dom.top_void_cell_index > 0) column_cells.push_back(dom.top_void_cell_index);
    // Faces between two packed cells (face n lies between cells n and n + 1).
    std::vector<int> bed_faces;
    for (size_t k = 0; k + 1 < bed_cells.size(); ++k)
      if (bed_cells[k + 1] == bed_cells[k] + 1) bed_faces.push_back(bed_cells[k]);
    const double z0 = (dom.bed_face_index_begin < static_cast<int>(dom.positions.size()))
                          ? dom.positions[dom.bed_face_index_begin] : 0.0;

    struct Inherited {
      bool max_pressure = false, min_pressure = false, max_wall_temperature = false, max_approach = false;
    };
    std::vector<Inherited> inherited;

    const auto& steps = data.cycle.step_list_;
    for (const auto& [step_index, entry] : cycle_it->second) {
      StepOperatingSummary s;
      s.step = step_index;
      if (step_index >= 0 && step_index < static_cast<int>(steps.size()))
        s.step_name = steps[step_index].step_name_;
      const auto& times = entry.step_data.times;

      auto rows = [&](VariableType type) -> const std::vector<Var1D>* {
        const auto it = entry.variables.find(type);
        return it == entry.variables.end() ? nullptr : &it->second.values;
      };
      auto scan = [](const std::vector<Var1D>* r, const std::vector<int>& idx, auto&& visit) {
        if (!r) return;
        for (size_t t = 0; t < r->size(); ++t)
          for (int i : idx)
            if (i >= 0 && i < static_cast<int>((*r)[t].size())) visit(t, i, (*r)[t][i]);
      };

      // A step's first row is the state its predecessor ended in, so an
      // extreme reached only there is inherited: the step before reports it.
      Inherited from_start;
      // Keeps the extreme (the larger value when `larger`) and whether it is
      // reached only in the first row.
      auto track = [](double& extreme, bool& at_start, double v, size_t t, bool larger) {
        const bool beyond = std::isnan(extreme) || (larger ? v > extreme : v < extreme);
        if (beyond) { extreme = v; at_start = t == 0; }
        else if (t > 0 && std::abs(v - extreme) <= 1e-9 * std::max(std::abs(extreme), 1.0)) at_start = false;
      };
      scan(rows(VariableType::PRESSURE), column_cells, [&](size_t t, int, double p) {
        if (!std::isfinite(p)) return;
        track(s.max_pressure, from_start.max_pressure, p, t, true);
        track(s.min_pressure, from_start.min_pressure, p, t, false);
      });
      scan(rows(VariableType::WALL_TEMPERATURE), bed_cells, [&](size_t t, int, double T) {
        if (std::isfinite(T)) track(s.max_wall_temperature, from_start.max_wall_temperature, T, t, true);
      });
      scan(rows(VariableType::APPROACH_TO_FLUIDISATION), bed_faces, [&](size_t t, int n, double a) {
        if (!std::isfinite(a)) return;
        const double before = s.max_approach_to_fluidisation;
        track(s.max_approach_to_fluidisation, from_start.max_approach, a, t, true);
        if (s.max_approach_to_fluidisation != before) {
          s.fluidisation_position = n < static_cast<int>(dom.positions.size()) ? dom.positions[n] - z0 : 0.0;
          s.fluidisation_time = (t < times.size() && !times.empty()) ? times[t] - times.front() : 0.0;
        }
      });
      // A step logged at one instant only has nothing to inherit.
      if (times.size() <= 1) from_start = {};
      inherited.push_back(from_start);
      report.steps.push_back(std::move(s));
    }

    // A warning switched off in the deck's settings is not reported.
    const WarningSettings& on = data.settings.warnings_;
    auto enabled = [&](OperatingLimit limit) {
      switch (limit) {
        case OperatingLimit::DesignPressure:    return on.design_pressure;
        case OperatingLimit::Vacuum:            return on.vacuum;
        case OperatingLimit::DesignTemperature: return on.design_temperature;
        case OperatingLimit::Fluidisation:      return on.fluidisation;
      }
      return true;
    };
    auto warn = [&](OperatingLimit limit, const StepOperatingSummary& s, double value, double limit_value,
                    std::string message) {
      if (!enabled(limit)) return;
      report.warnings.push_back({ limit, s.step, value, limit_value, "Warning: " + std::move(message) });
    };

    const auto& wall = data.bed.wall_;
    report.design_checked = wall.DesignInUse();
    if (report.design_checked) {
      const double P_design = wall.design_pressure_.Get();
      const double T_design = wall.design_temperature_.Get();
      for (size_t k = 0; k < report.steps.size(); ++k) {
        const auto& s = report.steps[k];
        const double p_gauge = s.max_pressure - kAtmosphericPressure;
        if (!inherited[k].max_pressure && std::isfinite(p_gauge)
            && p_gauge > std::max(P_design, 0.0) + kPressureMargin) {
          warn(OperatingLimit::DesignPressure, s, p_gauge, std::max(P_design, 0.0),
               "the pressure reaches " + BarGauge(p_gauge) + " " + InStep(s, report.cycle)
               + ", above the wall's design pressure of " + BarGauge(std::max(P_design, 0.0))
               + ". Raise the design pressure and recalculate the wall thickness, or lower the operating pressure.");
        }
        if (!wall.vacuum_service_ && !inherited[k].min_pressure && std::isfinite(s.min_pressure)
            && s.min_pressure < kAtmosphericPressure - kPressureMargin) {
          warn(OperatingLimit::Vacuum, s, s.min_pressure, kAtmosphericPressure,
               "the pressure falls to " + BarAbs(s.min_pressure) + " " + InStep(s, report.cycle)
               + ", below atmospheric, but the wall is not designed for vacuum service (external pressure).");
        }
        if (!inherited[k].max_wall_temperature && std::isfinite(s.max_wall_temperature)
            && s.max_wall_temperature > T_design + kTemperatureMargin) {
          warn(OperatingLimit::DesignTemperature, s, s.max_wall_temperature, T_design,
               "the wall temperature reaches " + Celsius(s.max_wall_temperature) + " " + InStep(s, report.cycle)
               + ", above the wall's design temperature of " + Celsius(T_design)
               + " (the allowable stress falls with temperature).");
        }
      }
    }

    for (size_t k = 0; k < report.steps.size(); ++k) {
      const auto& s = report.steps[k];
      if (!inherited[k].max_approach && std::isfinite(s.max_approach_to_fluidisation)
          && s.max_approach_to_fluidisation > 1.0) {
        std::ostringstream where;
        where << " (" << Num(s.fluidisation_position) << " m from the bottom of the bed, "
              << Num(s.fluidisation_time) << " s into the step)";
        warn(OperatingLimit::Fluidisation, s, s.max_approach_to_fluidisation, 1.0,
             "the approach to fluidisation reaches " + Num(s.max_approach_to_fluidisation) + " "
             + InStep(s, report.cycle) + where.str()
             + ": the upward pressure drop exceeds the weight of the bed, which would lift unless it is held down.");
      }
    }
    return report;
  }

} // namespace adsim
