/*
 * Skarstrom - pressure and temperature swing adsorption simulator
 * Copyright (c) 2024-2026 Sam Affleck and contributors
 *
 * Licensed under the terms in the LICENSE file at the repository root.
 */
#include "../pch.h"
#include "adsim/solver/kpi_history.h"

#include <cmath>

#include "adsim/solver/kpis.h"
#include "adsim/solver/log.h"

namespace adsim {

  std::vector<KPIRecord> ComputeCycleKPIs(const KPIConfig& config_in, const std::vector<Step>& steps,
                                          const EquationLog& eq_log, const BoundaryLog& b_log,
                                          int cycle) {
    std::vector<KPIRecord> out;
    KPIConfig config = config_in;
    config.ResolveIds(steps);
    for (const auto& cfg : config.per_component) {
      if (!cfg.tracked) continue;
      KPIParams p(eq_log, b_log);
      p.cycle_index_ = cycle;
      p.product_comp_ = cfg.product_component;
      p.feed_endpoints = cfg.feed_endpoints;
      p.product_endpoints = cfg.product_endpoints;
      p.recycle_endpoints = cfg.recycle_endpoints;
      KPIRecord r;
      r.cycle = cycle;
      r.component = cfg.product_component;
      r.valid = p.Check();
      if (r.valid) {
        constexpr double nan = std::numeric_limits<double>::quiet_NaN();
        bool ok = true;
        auto value = [&](double v) { return (ok && std::isfinite(v)) ? v : nan; };
        double v = CalculatePurity(p, &ok);                     r.purity = value(v);
        v = CalculateRecovery(p, &ok);                          r.recovery = value(v);
        v = CalculateProductivityMass(p, &ok);                  r.productivity_mass = value(v);
        v = CalculateProductivityVolume(p, &ok);                r.productivity_volume = value(v);
        v = CalculateSpecificEnergyConsumption(p, &ok);         r.energy_consumption = value(v);
        v = CalculateProductionRate(p, &ok);                    r.production_rate = value(v);
      }
      out.push_back(std::move(r));
    }
    return out;
  }

  void KPIHistory::Record(const KPIConfig& config, const std::vector<Step>& steps,
                          const EquationLog& eq_log, const BoundaryLog& b_log, int cycle) {
    records_[cycle] = ComputeCycleKPIs(config, steps, eq_log, b_log, cycle);
  }

  const std::vector<KPIRecord>& KPIHistory::Cycle(int cycle) const {
    static const std::vector<KPIRecord> kEmpty;
    const auto it = records_.find(cycle);
    return it == records_.end() ? kEmpty : it->second;
  }

  std::vector<int> KPIHistory::Cycles() const {
    std::vector<int> out;
    out.reserve(records_.size());
    for (const auto& [c, rows] : records_) {
      (void)rows;
      out.push_back(c);
    }
    return out;
  }

}  // namespace adsim
