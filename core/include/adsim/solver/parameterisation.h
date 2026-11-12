#pragma once

// STL includes
#include <algorithm>
#include <cmath>
#include <limits>
#include <map>
#include <string>
#include <unordered_map>
#include <vector>

// adsim includes
#include "adsim/solver/bed_series.h"
#include "adsim/solver/log.h"
#include "adsim/solver/kpis.h"
#include "adsim/solver/experimental_data.h"
#include "adsim/solver/kpi_history.h"


namespace adsim {

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
  // Time series of a bulk VariableType over one cycle of eq_log (the FINAL
  // logged cycle by default, or `cycle_index` when it is >= 0), concatenated
  // across all steps. Times are relative to the start of that cycle (first
  // point at t = 0, last point at ~cycle duration), independent of how many
  // cycles were run before it.
  //
  // Each row is read at the bed end or at z_pos with the positions and
  // bed-end values of ExtractBedTimeSeries (bed_series.h): cell centres
  // (never the ghost cells) or faces, and at Top / Bottom the physical
  // boundary value. Pass `sources` (the run's BoundaryLog and steps) so the
  // bed-end pressure, temperature and mole fraction are the endpoint logs,
  // as Simulation::Series does; without them they are the face values
  // reconstructed from the ghost cells.
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
    int                  cycle_index = -1,  // < 0: the final logged cycle
    const BedBoundarySources& sources = {})
  {
    return ExtractBedTimeSeries(eq_log, BedVariable::Bulk(var_type), spatial_loc, z_pos, out_times, out_values,
                                cycle_index, sources);
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
    int                   cycle_index = -1,  // < 0: the final logged cycle
    const BedBoundarySources& sources = {})
  {
    return ExtractBedTimeSeries(eq_log, BedVariable::Component(comp_var_type, component_name), spatial_loc, z_pos,
                                out_times, out_values, cycle_index, sources);
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
        case KPIVariableType::THERMAL_ENERGY:      kpi_val = CalculateSpecificThermalEnergy(kpi_p);    break;
        case KPIVariableType::COOLING_DUTY:        kpi_val = CalculateSpecificCoolingDuty(kpi_p);      break;
      }

      double res = kpi_val - exp_data.y[i];
      sse += res * res;
      ++n;
    }

    return (n > 0) ? std::sqrt(sse / n) : 1e10;
  }

  // ---------------------------------------------------------------------------
  // The simulated counterpart of a bound dataset, and the comparison
  //
  // SimulatedSeries reads what a dataset's observable binding names
  // (experimental_data.h) from a run: a bulk or component variable at the top,
  // the bottom or a position over one logged cycle (ExtractBedTimeSeries, with
  // the bed-end values of Simulation::Series when `sources` is given), or a
  // per-cycle KPI from `kpis` (Simulation::KPIs(): every completed cycle, not
  // only the logged ones). The result is in the DATASET's units: x is time
  // from the start of the cycle divided by x_scale, or the cycle number
  // (solver index - cycle_offset); y is converted from SI to y_unit. NaN KPIs
  // (undefined for the run) are left out.
  //
  // Throws std::invalid_argument for an unbound or invalid binding, and
  // std::out_of_range when the run did not log the cycle or has no data for
  // the observable.
  // ---------------------------------------------------------------------------

  struct SimulatedCounterpart {
    std::vector<double> x;  // in the dataset's x unit (time / x_scale, or cycle number)
    std::vector<double> y;  // in the dataset's y unit
    std::string si_unit;    // the observable's SI unit (before conversion)
  };

  SimulatedCounterpart SimulatedSeries(const XYData& data, const EquationLog& eq_log,
                                       const std::vector<KPIRecord>& kpis,
                                       const BedBoundarySources& sources = {});

  // The comparison of a dataset with its simulated counterpart, in the
  // dataset's y unit. A time series is interpolated linearly at the measured x
  // values it covers (points outside the simulated time range are skipped,
  // never extrapolated); a KPI is matched by cycle number. Points where the
  // simulation is not finite are skipped.
  //   rmse          sqrt(mean(r^2)), r = simulated - measured
  //   nrmse         rmse / (max - min of ALL the measured y); NaN for a flat dataset
  //   max_abs_error max |r|
  //   points        how many measured points were compared (0: the others are NaN)
  //   weighted_sse  weight * sum((r / sigma)^2), sigma the point's uncertainty
  //                 (1 when the dataset has none): a least-squares objective term
  struct DatasetComparison {
    double rmse = std::numeric_limits<double>::quiet_NaN();
    double nrmse = std::numeric_limits<double>::quiet_NaN();
    double max_abs_error = std::numeric_limits<double>::quiet_NaN();
    size_t points = 0;
    double weighted_sse = std::numeric_limits<double>::quiet_NaN();
  };

  DatasetComparison CompareWithSimulation(const XYData& data, const SimulatedCounterpart& simulated);

} // end adsim namespace
