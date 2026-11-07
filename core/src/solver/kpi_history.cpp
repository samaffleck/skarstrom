/*
 * Skarstrom - pressure and temperature swing adsorption simulator
 * Copyright (c) 2024-2026 Sam Affleck and contributors
 *
 * Licensed under the terms in the LICENSE file at the repository root.
 */
#include "../pch.h"
#include "adsim/solver/kpi_history.h"

#include <algorithm>
#include <cmath>

#include "adsim/solver/kpis.h"
#include "adsim/solver/log.h"

namespace adsim {

  double CycleMassBalanceError(const EquationLog& eq_log, const BoundaryLog& b_log, int cycle,
                               const std::string& component, KPIIntegration integration) {
    constexpr double nan = std::numeric_limits<double>::quiet_NaN();
    const auto& steps = eq_log.sys_data.step_list_;
    const int n_steps = static_cast<int>(steps.size());
    if (n_steps == 0) return nan;
    const StepBalance* first = b_log.GetStepBalance(cycle, 0);
    const StepBalance* last = b_log.GetStepBalance(cycle, n_steps - 1);
    if (!first || !last) return nan;
    const auto& comps = first->components;
    const auto it = std::find(comps.begin(), comps.end(), component);
    if (it == comps.end()) return nan;
    const size_t c = static_cast<size_t>(it - comps.begin());
    if (c >= first->inventory_start.size() || c >= last->inventory_end.size()) return nan;
    const double accumulation = last->inventory_end[c] - first->inventory_start[c];

    double net_in = 0.0, entered = 0.0;
    for (int s = 0; s < n_steps; ++s) {
      const auto names = b_log.stream_names.find(steps[s]);
      if (names == b_log.stream_names.end()) continue;
      for (const auto& name : names->second) {
        const StreamData* d = b_log.GetStreamData(cycle, s, name);
        if (!d) return nan;
        KPIEndpointRef ref;
        ref.step_index = s;
        ref.endpoint_name = name;
        ref.location = d->location_;
        const double n = EndpointComponentMolesIn(b_log, cycle, ref, component, integration);
        net_in += n;
        entered += std::max(n, 0.0);
      }
    }
    if (!std::isfinite(net_in) || !(entered > 0.0)) return nan;
    return (net_in - accumulation) / entered;
  }

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
        r.electrical_energy = r.energy_consumption;
        v = CalculateSpecificThermalEnergy(p, &ok);             r.thermal_energy = value(v);
        v = CalculateSpecificCoolingDuty(p, &ok);               r.cooling_duty = value(v);
        // Per step: the step's heat over the same net product mass.
        ok = true;
        const double mass = NetProductMass(p, ok);
        const int n_steps = static_cast<int>(eq_log.sys_data.step_list_.size());
        for (int s = 0; s < n_steps; ++s) {
          const StepThermalDuty d = CalculateStepThermalDuty(p, s);
          bool ok_h = ok, ok_c = ok;
          const double h = KPIRatio(d.heating, mass, ok_h);
          const double c = KPIRatio(d.cooling, mass, ok_c);
          r.step_thermal_energy.push_back(ok_h ? h : nan);
          r.step_cooling_duty.push_back(ok_c ? c : nan);
        }

        // The same KPIs from the quadratures.
        p.integration_ = KPIIntegration::QUADRATURE;
        ok = true;
        v = CalculatePurity(p, &ok);                    r.quadrature_purity = value(v);
        v = CalculateRecovery(p, &ok);                  r.quadrature_recovery = value(v);
        v = CalculateProductivityMass(p, &ok);          r.quadrature_productivity_mass = value(v);
        v = CalculateProductivityVolume(p, &ok);        r.quadrature_productivity_volume = value(v);
        v = CalculateSpecificEnergyConsumption(p, &ok); r.quadrature_energy_consumption = value(v);
        v = CalculateProductionRate(p, &ok);            r.quadrature_production_rate = value(v);
      }
      r.mass_balance_error =
          CycleMassBalanceError(eq_log, b_log, cycle, cfg.product_component, KPIIntegration::TRAPEZOID);
      r.quadrature_mass_balance_error =
          CycleMassBalanceError(eq_log, b_log, cycle, cfg.product_component, KPIIntegration::QUADRATURE);
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
