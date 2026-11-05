/*
 * Skarstrom - pressure and temperature swing adsorption simulator
 * Copyright (c) 2024-2026 Sam Affleck and contributors
 *
 * Licensed under the terms in the LICENSE file at the repository root.
 */
#pragma once

// KPIs of every cycle of a run.
//
// The equation and boundary logs keep only the last settings.max_cycles_to_keep
// cycles (their spatial profiles are large), so KPIs computed from them on
// demand cover only those cycles. The solvers therefore evaluate the KPIs of
// every tracked component when a cycle completes and keep the few numbers in a
// KPIHistory, which is never evicted: a run's KPI history spans all its cycles.

#include <limits>
#include <map>
#include <string>
#include <vector>

#include "adsim/solver/kpi_config.h"

namespace adsim {

  class EquationLog;
  struct BoundaryLog;

  // KPIs of one tracked component at one cycle. A KPI that is not defined for
  // the run (no product drawn, zero denominator) is NaN, and so is every KPI
  // when the component's KPI configuration is incomplete (`valid` false: no
  // feed or no product endpoint).
  struct KPIRecord {
    int cycle = 0;
    std::string component;
    bool valid = false;
    double purity = std::numeric_limits<double>::quiet_NaN();               // [mol/mol]
    double recovery = std::numeric_limits<double>::quiet_NaN();             // [mol/mol]
    double productivity_mass = std::numeric_limits<double>::quiet_NaN();    // [kg/(kg s)]
    double productivity_volume = std::numeric_limits<double>::quiet_NaN();  // [kg/(m3 s)]
    double energy_consumption = std::numeric_limits<double>::quiet_NaN();   // [J/kg]
    double production_rate = std::numeric_limits<double>::quiet_NaN();      // [mol/s]
  };

  // KPIs of every tracked component of `config` at one logged cycle, computed
  // from the logs. The endpoint references are re-bound by ID against `steps`
  // first (on a copy), so a step or stream renamed since they were resolved
  // still counts. Empty when no component is tracked.
  std::vector<KPIRecord> ComputeCycleKPIs(const KPIConfig& config, const std::vector<Step>& steps,
                                          const EquationLog& eq_log, const BoundaryLog& b_log,
                                          int cycle);

  // The KPIs of each completed cycle, recorded when the cycle completed.
  class KPIHistory {
  public:
    void Clear() { records_.clear(); }
    // Compute and keep the KPIs of `cycle` (replaces an earlier record of it).
    void Record(const KPIConfig& config, const std::vector<Step>& steps,
                const EquationLog& eq_log, const BoundaryLog& b_log, int cycle);

    bool Has(int cycle) const { return records_.find(cycle) != records_.end(); }
    // The records of one cycle (empty when it was not recorded).
    const std::vector<KPIRecord>& Cycle(int cycle) const;
    // Recorded cycle indices, ascending.
    std::vector<int> Cycles() const;
    const std::map<int, std::vector<KPIRecord>>& Records() const { return records_; }

  private:
    std::map<int, std::vector<KPIRecord>> records_;
  };

}  // namespace adsim
