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
#include <set>

#include "adsim/solver/bed_series.h"
#include "adsim/solver/kpis.h"
#include "adsim/solver/log.h"

namespace adsim {

  namespace {
    constexpr double kNaN = std::numeric_limits<double>::quiet_NaN();

    // Moles of `component` that entered the bed over the cycle through every
    // endpoint of every step (the inflow only), from the KPI quadratures.
    double AllEndpointInflow(const EquationLog& eq_log, const BoundaryLog& b_log, int cycle,
                             const std::string& component) {
      const auto& steps = eq_log.sys_data.step_list_;
      double entered = 0.0;
      for (int s = 0; s < static_cast<int>(steps.size()); ++s) {
        const auto names = b_log.stream_names.find(steps[s]);
        if (names == b_log.stream_names.end()) continue;
        for (const auto& name : names->second) {
          const StreamData* d = b_log.GetStreamData(cycle, s, name);
          if (!d) return kNaN;
          KPIEndpointRef ref;
          ref.step_index = s;
          ref.endpoint_name = name;
          ref.location = d->location_;
          entered += std::max(EndpointComponentMolesIn(b_log, cycle, ref, component), 0.0);
        }
      }
      return entered;
    }

    // Mass [kg] of the solids that carry the logged reactions.
    double CatalystMass(const EquationLog& eq_log) {
      const DomainVectors& domain = eq_log.GetDomain();
      std::set<AdsorbentID> catalysts;
      for (const auto& r : eq_log.Reactions()) catalysts.insert(r.adsorbent);
      double mass = 0.0;
      for (size_t i = 0; i < domain.cell_adsorbents.size() && i < domain.volume.size(); ++i)
        for (const auto& ca : domain.cell_adsorbents[i])
          if (catalysts.count(ca.id)) mass += ca.bulk_density * domain.volume[i];
      return mass;
    }

    double RatioOrNaN(double num, double den) {
      bool ok = true;
      const double r = KPIRatio(num, den, ok);
      return ok ? r : kNaN;
    }
  }  // namespace

  PeakTemperature CyclePeakTemperature(const EquationLog& eq_log, int cycle) {
    PeakTemperature peak;
    const auto& sim_log = eq_log.GetSimulationLog();
    const auto it = sim_log.find(cycle);
    if (it == sim_log.end()) return peak;
    const DomainVectors& domain = eq_log.GetDomain();
    const std::vector<double> centres = CellCentrePositions(domain);
    const VariableType type = eq_log.NonLocal() ? VariableType::SOLID_TEMPERATURE : VariableType::TEMPERATURE;
    // The packed cells; every interior cell when the domain lists none.
    std::vector<size_t> cells;
    for (int i : domain.bed_cell_indices)
      if (i > 0) cells.push_back(static_cast<size_t>(i));
    if (cells.empty())
      for (size_t i = 1; i <= centres.size(); ++i) cells.push_back(i);
    double cycle_start = kNaN;
    for (const auto& [step, log] : it->second) {
      const VariableLogEntry e = eq_log.GetVariable(cycle, step, type);
      if (e.times.empty()) continue;
      if (!std::isfinite(cycle_start)) cycle_start = e.times.front();
      for (size_t t = 0; t < e.values.size() && t < e.times.size(); ++t) {
        const auto& row = e.values[t];
        for (size_t i : cells) {
          // Row index i is cell i (ghost cell 0 first); its centre is centres[i - 1].
          if (i >= row.size() || i - 1 >= centres.size()) continue;
          if (!std::isfinite(row[i])) continue;
          if (!(row[i] <= peak.temperature)) {  // also takes the first (NaN) comparison
            peak.temperature = row[i];
            peak.position = centres[i - 1];
            peak.time = e.times[t] - cycle_start;
          }
        }
      }
    }
    return peak;
  }

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
    double generated = 0.0;  // by the reactions (0 without)
    for (int s = 0; s < n_steps; ++s) {
      generated += eq_log.ReactionGeneration(cycle, s, component);
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
    // A component that (to round-off) never entered has no balance: its
    // round-off over its round-off is meaningless.
    double inventory = 0.0;
    for (double v : first->inventory_start) inventory += std::abs(v);
    if (!std::isfinite(net_in) || !(entered > 1e-12 * inventory)) return nan;
    return (net_in + generated - accumulation) / entered;
  }

  double CycleMassBalanceError(const EquationLog& eq_log, const BoundaryLog& b_log, int cycle,
                               const std::string& component) {
    return CycleMassBalanceError(eq_log, b_log, cycle, component, KPIIntegration::QUADRATURE);
  }

  std::vector<KPIRecord> ComputeCycleKPIs(const KPIConfig& config_in, const std::vector<Step>& steps,
                                          const EquationLog& eq_log, const BoundaryLog& b_log,
                                          int cycle) {
    std::vector<KPIRecord> out;
    KPIConfig config = config_in;
    config.ResolveIds(steps);
    // Reaction KPIs: only for a deck with reactions and a key reactant.
    const ReactionKPIConfig& rx = config.reaction;
    const bool reactive = rx.IsSet() && !eq_log.Reactions().empty();
    double g_key = kNaN, f_key = kNaN, sty_denominator = kNaN;
    PeakTemperature peak;
    if (reactive) {
      g_key = 0.0;
      for (int s = 0; s < static_cast<int>(eq_log.sys_data.step_list_.size()); ++s)
        g_key += eq_log.ReactionGeneration(cycle, s, rx.key_reactant);
      // F_A,in: the key reactant's KPI feed endpoints, else every endpoint's inflow of it.
      const ComponentKPIConfig* key_cfg = config.FindByComponent(rx.key_reactant);
      if (key_cfg && key_cfg->HasAnyFeed()) {
        f_key = 0.0;
        for (const auto& ep : key_cfg->feed_endpoints)
          f_key += EndpointComponentMolesIn(b_log, cycle, ep, rx.key_reactant);
      } else {
        f_key = AllEndpointInflow(eq_log, b_log, cycle, rx.key_reactant);
      }
      const double extent = rx.basis == SpaceTimeYieldBasis::BED_VOLUME ? PackedBedVolume(eq_log.GetDomain()) : CatalystMass(eq_log);
      sty_denominator = eq_log.sys_data.total_cycle_time * extent;
      peak = CyclePeakTemperature(eq_log, cycle);
    }

    for (const auto& cfg : config.per_component) {
      const bool is_key = reactive && cfg.product_component == rx.key_reactant;
      const bool is_product = reactive && rx.IsProduct(cfg.product_component);
      if (!cfg.tracked && !is_key && !is_product) continue;
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

        // The same KPIs from the trapezoid over the logged samples.
        p.integration_ = KPIIntegration::TRAPEZOID;
        ok = true;
        v = CalculatePurity(p, &ok);                    r.trapezoid_purity = value(v);
        v = CalculateRecovery(p, &ok);                  r.trapezoid_recovery = value(v);
        v = CalculateProductivityMass(p, &ok);          r.trapezoid_productivity_mass = value(v);
        v = CalculateProductivityVolume(p, &ok);        r.trapezoid_productivity_volume = value(v);
        v = CalculateSpecificEnergyConsumption(p, &ok); r.trapezoid_energy_consumption = value(v);
        v = CalculateProductionRate(p, &ok);            r.trapezoid_production_rate = value(v);
      }
      r.mass_balance_error =
          CycleMassBalanceError(eq_log, b_log, cycle, cfg.product_component, KPIIntegration::QUADRATURE);
      r.trapezoid_mass_balance_error =
          CycleMassBalanceError(eq_log, b_log, cycle, cfg.product_component, KPIIntegration::TRAPEZOID);
      if (is_key || is_product) {
        r.reaction_kpis = true;
        double g = 0.0;
        for (int s = 0; s < static_cast<int>(eq_log.sys_data.step_list_.size()); ++s)
          g += eq_log.ReactionGeneration(cycle, s, cfg.product_component);
        r.reaction_generation = g;
        if (is_key) r.conversion = RatioOrNaN(-g_key, f_key);
        if (is_product) {
          r.selectivity = RatioOrNaN(g, -g_key);
          r.reaction_yield = RatioOrNaN(g, f_key);
          const auto& mr = eq_log.sys_data.molar_mass;
          const auto m = mr.find(cfg.product_component);
          if (m != mr.end()) r.space_time_yield = RatioOrNaN(g * m->second, sty_denominator);
          r.space_time_yield_unit = rx.basis == SpaceTimeYieldBasis::BED_VOLUME ? "kg/(m3 s)" : "kg/(kg s)";
        }
        r.peak_temperature = peak.temperature;
        r.peak_temperature_position = peak.position;
        r.peak_temperature_time = peak.time;
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
