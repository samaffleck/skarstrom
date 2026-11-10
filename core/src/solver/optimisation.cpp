/*
 * Skarstrom - pressure and temperature swing adsorption simulator
 * Copyright (c) 2024-2026 Sam Affleck and contributors
 *
 * Licensed under the terms in the LICENSE file at the repository root.
 */
#include "adsim/solver/optimisation.h"

#include <cmath>
#include <limits>
#include <memory>
#include <stdexcept>

namespace adsim {

  namespace {

    constexpr double kNaN = std::numeric_limits<double>::quiet_NaN();

    // The KPI of a record for an optimiser KPI type; NaN when the run does not
    // report it (LCOC needs cost data the run does not use).
    double KPIValue(const KPIRecord& k, KPIType type) {
      switch (type) {
        case KPIType::PURITY: return k.purity;
        case KPIType::RECOVERY: return k.recovery;
        case KPIType::PRODUCTIVITY_MASS: return k.productivity_mass;
        case KPIType::PRODUCTIVITY_VOLUME: return k.productivity_volume;
        case KPIType::ENERGY_CONSUMPTION: return k.energy_consumption;
        case KPIType::PRODUCTION_RATE: return k.production_rate;
        case KPIType::THERMAL_ENERGY: return k.thermal_energy;
        case KPIType::COOLING_DUTY: return k.cooling_duty;
        case KPIType::LCOC: return kNaN;
        default: return kNaN;  // the step KPI types: StepKPIValue
      }
      return kNaN;
    }

    // The step output of a step-end pressure record for a step KPI type.
    double StepKPIValue(const StepEndPressure& r, KPIType type) {
      switch (type) {
        case KPIType::END_PRESSURE_BOTTOM: return r.bottom;
        case KPIType::END_PRESSURE_TOP: return r.top;
        case KPIType::END_PRESSURE_MIN: return r.bed_min;
        case KPIType::END_PRESSURE_MEAN: return r.bed_mean;
        case KPIType::END_PRESSURE_MAX: return r.bed_max;
        default: return kNaN;
      }
    }

    constexpr KPIType kStepKPIs[] = {KPIType::END_PRESSURE_BOTTOM, KPIType::END_PRESSURE_TOP,
                                     KPIType::END_PRESSURE_MIN, KPIType::END_PRESSURE_MEAN,
                                     KPIType::END_PRESSURE_MAX};

    constexpr KPIType kRunKPIs[] = {KPIType::PURITY,           KPIType::RECOVERY,
                                    KPIType::PRODUCTIVITY_MASS, KPIType::PRODUCTIVITY_VOLUME,
                                    KPIType::ENERGY_CONSUMPTION, KPIType::PRODUCTION_RATE,
                                    KPIType::THERMAL_ENERGY,    KPIType::COOLING_DUTY};

    void LoadDeck(Simulation& sim, const std::string& deck_json,
                  const std::vector<std::pair<std::string, std::string>>& settings) {
      try {
        sim.LoadJson(deck_json);
      } catch (const std::exception& e) {
        throw std::invalid_argument(std::string("the deck does not load: ") + e.what());
      }
      for (const auto& [path, text] : settings) {
        try {
          sim.SetParameterText(path, text);
        } catch (const std::exception& e) {
          throw std::invalid_argument("setting '" + path + "': " + e.what());
        }
      }
    }

    OptimizationEvaluation Evaluate(const std::string& deck_json, const std::vector<OptimizationVariable>& variables,
                                    const std::vector<double>& x, const OptimizationCostFunction& cost,
                                    const SimulationOptimizerSettings& settings, CycleStateCache* cache, int index) {
      OptimizationEvaluation ev;
      ev.index = index;
      ev.x = x;
      Simulation sim;
      try {
        LoadDeck(sim, deck_json, settings.settings);
        for (size_t i = 0; i < variables.size(); ++i) sim.SetParameter(variables[i].path, x[i]);
      } catch (const std::exception& e) {
        ev.status = RunStatus::Failed;
        ev.message = e.what();
        return ev;
      }

      RunOptions options;
      if (settings.cycles > 0) options.cycles = settings.cycles;
      options.cancel = settings.cancel;
      CachedRunInfo info;
      const RunResult r = RunWithStateCache(sim, cache, x, options, &info);
      ev.status = r.status;
      ev.message = r.message;
      ev.warm = info.warm;
      ev.distance = info.distance;
      ev.note = info.note;
      ev.css_reached = r.css_reached;
      ev.cycles_run = r.cycles_run;
      ev.wall_time_s = r.wall_time_s;
      ev.accepted = r.Ok() && (!settings.require_css || r.css_reached);
      if (r.Ok() && settings.require_css && !r.css_reached) {
        ev.message = "did not reach cyclic steady state in " + std::to_string(r.cycles_run) + " cycles";
      }
      if (!ev.accepted) return ev;

      const std::vector<KPIRecord> records = sim.KPIs(-1);
      std::vector<StepEndPressure> step_pressures;
      for (const std::string& component : cost.GetReferencedComponents()) {
        auto& values = ev.kpis[component];
        if (IsStepReference(component)) {
          for (KPIType t : kStepKPIs) values[t] = kNaN;
          if (step_pressures.empty()) step_pressures = sim.StepEndPressures(-1);
          const int s = ResolveStepReference(sim.Data().cycle.step_list_, component);
          if (s >= 0 && static_cast<size_t>(s) < step_pressures.size())
            for (KPIType t : kStepKPIs) values[t] = StepKPIValue(step_pressures[static_cast<size_t>(s)], t);
          continue;
        }
        for (KPIType t : kRunKPIs) values[t] = kNaN;
        for (const auto& k : records) {
          if (k.component != component) continue;
          for (KPIType t : kRunKPIs) values[t] = KPIValue(k, t);
          break;
        }
      }
      ev.cost = cost.ComputeCost(ev.kpis);
      return ev;
    }

  }  // namespace

  void CheckOptimization(const std::string& deck_json, const std::vector<OptimizationVariable>& variables,
                         const OptimizationCostFunction& cost, const SimulationOptimizerSettings& settings) {
    if (variables.empty()) throw std::invalid_argument("an optimisation needs at least one variable");
    if (cost.GetObjectives().empty()) throw std::invalid_argument("an optimisation needs at least one objective");
    if (settings.n_init_samples + settings.n_iterations < 1)
      throw std::invalid_argument("an optimisation needs at least one evaluation");
    Simulation sim;
    LoadDeck(sim, deck_json, settings.settings);
    ParameterRegistry& registry = sim.Parameters();
    std::vector<std::string> seen;
    for (const auto& v : variables) {
      Parameter* p = registry.Find(v.path);
      if (!p) throw std::invalid_argument("variable '" + v.path + "': no such parameter in the deck");
      if (!p->Info().tunable) throw std::invalid_argument("variable '" + v.path + "' is not a tunable parameter");
      const std::string& path = p->Info().path;
      for (const auto& s : seen)
        if (s == path) throw std::invalid_argument("variable '" + v.path + "' is listed twice");
      seen.push_back(path);
      if (!std::isfinite(v.lower) || !std::isfinite(v.upper) || !(v.lower < v.upper)) {
        throw std::invalid_argument("variable '" + v.path + "' needs finite bounds with lower < upper");
      }
    }
    auto CheckTerm = [&](const std::string& component, KPIType type) {
      if (type == KPIType::LCOC) {
        throw std::invalid_argument("the levelised cost (LCOC) is not available to OptimizeSimulation");
      }
      if (IsStepReference(component) != IsStepKPIType(type)) {
        throw std::invalid_argument(IsStepReference(component)
                                        ? "'" + component + "' is a step: its terms take a step-end pressure KPI"
                                        : "a step-end pressure KPI needs a step reference step[<ID or name>], not '" +
                                              component + "'");
      }
      if (IsStepReference(component)) {
        try {
          ResolveStepReference(sim.Data().cycle.step_list_, component);
        } catch (const std::exception& e) {
          throw std::invalid_argument(e.what());
        }
        return;
      }
      bool tracked = false;
      for (const auto& c : sim.Data().kpi_config.per_component)
        if (c.tracked && c.product_component == component) tracked = true;
      if (!tracked) {
        throw std::invalid_argument("component '" + component +
                                    "' is not tracked by the deck's KPI configuration");
      }
    };
    for (const auto& [id, term] : cost.GetObjectives()) {
      (void)id;
      CheckTerm(term.product_component, term.type);
    }
    for (const auto& [id, term] : cost.GetConstraints()) {
      (void)id;
      CheckTerm(term.product_component, term.type);
    }
  }

  SimulationOptimizationResult OptimizeSimulation(const std::string& deck_json,
                                                  const std::vector<OptimizationVariable>& variables,
                                                  const OptimizationCostFunction& cost,
                                                  const SimulationOptimizerSettings& settings) {
    CheckOptimization(deck_json, variables, cost, settings);
    SimulationOptimizationResult result;

    std::vector<double> lower, upper;
    for (const auto& v : variables) {
      lower.push_back(v.lower);
      upper.push_back(v.upper);
    }
    std::unique_ptr<CycleStateCache> cache;
    if (settings.warm_start) cache = std::make_unique<CycleStateCache>(lower, upper, settings.cache_capacity);

    OptimizerSettings search;
    search.n_init_samples = settings.n_init_samples;
    search.n_iterations = settings.n_iterations;
    auto objective = [&](const std::vector<double>& x) -> double {
      OptimizationEvaluation ev;
      if (result.cancelled) {
        // dlib asks for its whole budget: answer the rest without running.
        ev.index = static_cast<int>(result.evaluations.size());
        ev.x = x;
        ev.status = RunStatus::Cancelled;
        ev.message = "optimisation cancelled";
      } else {
        ev = Evaluate(deck_json, variables, x, cost, settings, cache.get(),
                      static_cast<int>(result.evaluations.size()));
        if (ev.status == RunStatus::Cancelled) result.cancelled = true;
        if (settings.on_evaluation) settings.on_evaluation(ev);
      }
      result.evaluations.push_back(ev);
      return ev.cost;
    };
    GlobalOptimizer optimizer(lower, upper, objective, search);
    result.x = optimizer.OptimizeParameters();
    result.search_cost = optimizer.GetBestValue();
    for (const auto& ev : result.evaluations) {
      if (ev.x == result.x) {
        result.best_index = ev.index;
        break;
      }
    }
    result.cost = result.search_cost;

    if (settings.rerun_cold && !result.cancelled) {
      result.cold = Evaluate(deck_json, variables, result.x, cost, settings, nullptr, -1);
      result.rerun_cold = true;
      result.cost = result.cold.cost;
      if (result.cold.status == RunStatus::Cancelled) result.cancelled = true;
      if (settings.on_evaluation) settings.on_evaluation(result.cold);
    }
    return result;
  }

}  // namespace adsim
