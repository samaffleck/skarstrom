#pragma once

// STL includes
#include <algorithm>
#include <cmath>
#include <map>
#include <string>
#include <unordered_map>
#include <vector>

// adsim includes
#include "adsim/solver/log.h"
#include "adsim/solver/kpis.h"
#include "adsim/solver/experimental_data.h"


namespace adsim {

  // ---------------------------------------------------------------------------
  // Spatial location for time-series extraction
  // ---------------------------------------------------------------------------

  enum class ParamSpatialLoc { Top = 0, Bottom = 1, Position = 2 };

  // ---------------------------------------------------------------------------
  // Linear interpolation
  // Clamps to boundary values outside the input range.
  // ---------------------------------------------------------------------------

  inline double InterpLinear(const std::vector<double>& x,
                              const std::vector<double>& y,
                              double x_query)
  {
    if (x.empty()) return 0.0;
    auto it  = std::upper_bound(x.begin(), x.end(), x_query);
    size_t idx = static_cast<size_t>(std::distance(x.begin(), it));
    if (idx == 0)         return y.front();
    if (idx >= x.size())  return y.back();
    size_t i0 = idx - 1;
    double w = (x_query - x[i0]) / (x[idx] - x[i0]);
    return y[i0] + w * (y[idx] - y[i0]);
  }

  // ---------------------------------------------------------------------------
  // CycleStartTime
  //
  // The solver logs ABSOLUTE IDA times (IDA time runs continuously across all
  // steps and cycles), so a cycle's time axis must be re-based on the first
  // time logged in that cycle. The first log of every step is taken at the
  // step's start time, so this is the cycle's start time. Uses the shared
  // per-step time vector, so it does not depend on which variable is queried.
  // Returns 0 if the cycle has no logged times.
  // ---------------------------------------------------------------------------

  inline double CycleStartTime(const std::map<int, Log>& cycle_log)
  {
    for (const auto& [step_idx, step_log] : cycle_log) {
      if (!step_log.step_data.times.empty())
        return step_log.step_data.times.front();
    }
    return 0.0;
  }

  // ---------------------------------------------------------------------------
  // ExtractBulkTimeSeries
  //
  // Extracts a time series of a bulk VariableType from one cycle of eq_log
  // (the FINAL logged cycle by default, or `cycle_index` when it is >= 0),
  // concatenating across all steps. Times are relative to the start of that
  // cycle (first point at t = 0, last point at ~cycle duration), independent
  // of how many cycles were run before it.
  //
  // For each timestep the spatial profile is interpolated at query_z using
  // InterpLinear and the cell-centre positions from eq_log.GetDomain().
  //
  // Returns false (and empty out vectors) if no data is available (including
  // a `cycle_index` that is not in the log).
  // ---------------------------------------------------------------------------

  inline bool ExtractBulkTimeSeries(
    const EquationLog&   eq_log,
    VariableType         var_type,
    ParamSpatialLoc      spatial_loc,
    double               z_pos,         // [m], used when spatial_loc == Position
    std::vector<double>& out_times,     // [s] from start of the selected cycle
    std::vector<double>& out_values,
    int                  cycle_index = -1)  // < 0: the final logged cycle
  {
    out_times.clear();
    out_values.clear();
    const auto& sim_log  = eq_log.GetSimulationLog();
    if (sim_log.empty()) return false;
    if (cycle_index >= 0 && sim_log.find(cycle_index) == sim_log.end()) return false;

    const int last_cycle = cycle_index >= 0 ? cycle_index : sim_log.rbegin()->first;
    const auto& positions = eq_log.GetDomain().positions;

    double query_z = 0.0;
    switch (spatial_loc) {
      case ParamSpatialLoc::Top:
        query_z = positions.empty() ? 0.0 : positions.back();
        break;
      case ParamSpatialLoc::Bottom:
        query_z = positions.empty() ? 0.0 : positions.front();
        break;
      case ParamSpatialLoc::Position:
        query_z = z_pos;
        break;
    }

    const auto& cycle_log = sim_log.at(last_cycle);
    const double cycle_start = CycleStartTime(cycle_log);

    for (const auto& [step_idx, step_log] : cycle_log) {
      VariableLogEntry vle = eq_log.GetVariable(last_cycle, step_idx, var_type);
      if (vle.times.empty()) continue;

      for (size_t t = 0; t < vle.times.size(); ++t) {
        out_times.push_back(vle.times[t] - cycle_start);
        double val = 0.0;
        if (t < vle.values.size() && !positions.empty()) {
          val = InterpLinear(positions, vle.values[t], query_z);
        }
        out_values.push_back(val);
      }
    }

    return !out_times.empty();
  }

  // ---------------------------------------------------------------------------
  // ExtractComponentTimeSeries
  //
  // Same as ExtractBulkTimeSeries but for a ComponentVariableType.
  // ---------------------------------------------------------------------------

  inline bool ExtractComponentTimeSeries(
    const EquationLog&    eq_log,
    ComponentVariableType comp_var_type,
    const std::string&    component_name,
    ParamSpatialLoc       spatial_loc,
    double                z_pos,
    std::vector<double>&  out_times,
    std::vector<double>&  out_values,
    int                   cycle_index = -1)  // < 0: the final logged cycle
  {
    out_times.clear();
    out_values.clear();
    const auto& sim_log = eq_log.GetSimulationLog();
    if (sim_log.empty()) return false;
    if (cycle_index >= 0 && sim_log.find(cycle_index) == sim_log.end()) return false;

    const int last_cycle = cycle_index >= 0 ? cycle_index : sim_log.rbegin()->first;
    const auto& positions = eq_log.GetDomain().positions;

    double query_z = 0.0;
    switch (spatial_loc) {
      case ParamSpatialLoc::Top:    query_z = positions.empty() ? 0.0 : positions.back();  break;
      case ParamSpatialLoc::Bottom: query_z = positions.empty() ? 0.0 : positions.front(); break;
      case ParamSpatialLoc::Position: query_z = z_pos; break;
    }

    const auto& cycle_log = sim_log.at(last_cycle);
    const double cycle_start = CycleStartTime(cycle_log);

    for (const auto& [step_idx, step_log] : cycle_log) {
      VariableLogEntry vle =
        eq_log.GetComponentVariable(last_cycle, step_idx, comp_var_type, component_name);
      if (vle.times.empty()) continue;

      for (size_t t = 0; t < vle.times.size(); ++t) {
        out_times.push_back(vle.times[t] - cycle_start);
        double val = 0.0;
        if (t < vle.values.size() && !positions.empty()) {
          val = InterpLinear(positions, vle.values[t], query_z);
        }
        out_values.push_back(val);
      }
    }

    return !out_times.empty();
  }

  // ---------------------------------------------------------------------------
  // ComputeRMSE
  //
  // Computes the root-mean-squared error between the simulated time series
  // and the experimental data points, using linear interpolation of the
  // simulated series at each experimental x value.
  //
  // time_unit_mul converts experimental x-axis values from user units to
  // the same unit as sim_times (seconds). E.g. 60.0 for minutes.
  //
  // Returns 1e10 if either series is empty or interpolation cannot proceed.
  // ---------------------------------------------------------------------------

  inline double ComputeRMSE(
    const std::vector<double>& sim_times,    // sorted, in base time unit [s]
    const std::vector<double>& sim_values,
    const XYData&              exp_data,
    double                     time_unit_mul = 1.0)
  {
    if (sim_times.empty() || exp_data.x.empty()) return 1e10;
    double sse = 0.0;
    int    n   = 0;
    for (size_t i = 0; i < exp_data.x.size(); ++i) {
      double x_s   = exp_data.x[i] * time_unit_mul;
      double y_sim = InterpLinear(sim_times, sim_values, x_s);
      double res   = y_sim - exp_data.y[i];
      sse += res * res;
      ++n;
    }
    return (n > 0) ? std::sqrt(sse / n) : 1e10;
  }

  // ---------------------------------------------------------------------------
  // ComputeKPIRMSE
  //
  // Computes RMSE between computed KPI values and experimental data where
  // exp_data.x contains cycle numbers and exp_data.y the observed KPI values.
  //
  // cycle_offset is added to each CSV cycle number to map it to the solver's
  // cycle index (e.g. -1 converts 1-based CSV numbering to 0-based solver).
  //
  // Only experimental points whose cycle index exists in the simulation log
  // contribute to the RMSE; unmatched points are silently skipped.
  // ---------------------------------------------------------------------------

  inline double ComputeKPIRMSE(
    const EquationLog&     eq_log,
    const BoundaryLog&     b_log,
    const XYData&          exp_data,
    int                    cycle_offset,
    KPIVariableType        kpi_type,
    const std::string&     product_comp,
    const std::vector<KPIEndpointRef>& feed_endpoints,
    const std::vector<KPIEndpointRef>& product_endpoints,
    const std::vector<KPIEndpointRef>& recycle_endpoints = {})
  {
    const auto& sim_log = eq_log.GetSimulationLog();
    double sse = 0.0;
    int    n   = 0;

    for (size_t i = 0; i < exp_data.x.size(); ++i) {
      int cycle_idx = static_cast<int>(std::round(exp_data.x[i])) + cycle_offset;
      if (sim_log.find(cycle_idx) == sim_log.end()) continue;

      KPIParams kpi_p(eq_log, b_log);
      kpi_p.cycle_index_                = cycle_idx;
      kpi_p.product_comp_               = product_comp;
      kpi_p.product_endpoints           = product_endpoints;
      kpi_p.feed_endpoints              = feed_endpoints;
      kpi_p.recycle_endpoints           = recycle_endpoints;

      double kpi_val = 0.0;
      switch (kpi_type) {
        case KPIVariableType::PURITY:              kpi_val = CalculatePurity(kpi_p);                   break;
        case KPIVariableType::RECOVERY:            kpi_val = CalculateRecovery(kpi_p);                 break;
        case KPIVariableType::PRODUCTIVITY_MASS:   kpi_val = CalculateProductivityMass(kpi_p);         break;
        case KPIVariableType::PRODUCTIVITY_VOLUME: kpi_val = CalculateProductivityVolume(kpi_p);       break;
        case KPIVariableType::ENERGY_CONSUMPTION:  kpi_val = CalculateSpecificEnergyConsumption(kpi_p);break;
        case KPIVariableType::PRODUCTION_RATE:     kpi_val = CalculateProductionRate(kpi_p);           break;
      }

      double res = kpi_val - exp_data.y[i];
      sse += res * res;
      ++n;
    }

    return (n > 0) ? std::sqrt(sse / n) : 1e10;
  }

} // end adsim namespace
