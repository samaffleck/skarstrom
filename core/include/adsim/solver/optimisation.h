/*
 * Skarstrom - pressure and temperature swing adsorption simulator
 * Copyright (c) 2024-2026 Sam Affleck and contributors
 *
 * Licensed under the terms in the LICENSE file at the repository root.
 */
#pragma once

// STL includes
#include <algorithm>
#include <cmath>
#include <functional>
#include <iostream>
#include <limits>
#include <stdexcept>
#include <string>
#include <string_view>
#include <unordered_map>
#include <unordered_set>
#include <utility>
#include <vector>

// dlib (global optimisation)
#include <dlib/global_optimization.h>

// AdSim includes
#include "adsim/solver/kpis.h"
#include "adsim/objects/user_data.h"
#include "adsim/parameters/parameter_registry.h"
#include "adsim/simulation.h"            // RunStatus, CancelPoll
#include "adsim/solver/state_cache.h"    // CycleStateCache


namespace adsim {

  // ENERGY_CONSUMPTION is the specific electrical energy (compressor and
  // vacuum-pump work per kg of net product); ELECTRICAL_ENERGY is the same
  // enumerator under its explicit name, so objectives and constraints that
  // use either mean the same KPI. New KPIs are appended at the end.
  enum class KPIType {
    PURITY,
    RECOVERY,
    PRODUCTIVITY_MASS,
    PRODUCTIVITY_VOLUME,
    ENERGY_CONSUMPTION,
    PRODUCTION_RATE,
    LCOC,
    THERMAL_ENERGY,
    COOLING_DUTY,
    // Aliases last: an enumerator after them would take the alias's value + 1.
    ELECTRICAL_ENERGY = ENERGY_CONSUMPTION
  };

  inline constexpr std::string_view kKPITypeList[] = {
    "Purity",
    "Recovery",
    "Productivity (mass)",
    "Productivity (volume)",
    "Electrical Energy",
    "Production Rate",
    "LCOC",
    "Thermal Energy",
    "Cooling Duty"
  };

  static std::string GetKPITypeUnits(KPIType type) {
    switch (type) {
      case KPIType::PURITY: return "[mol/mol]";
      case KPIType::RECOVERY: return "[mol/mol]";
      case KPIType::PRODUCTIVITY_MASS: return "[kg_gas/(kg_solid·s)]";
      case KPIType::PRODUCTIVITY_VOLUME: return "[kg_gas/(m3·s)]";
      case KPIType::ENERGY_CONSUMPTION: return "[J/kg]";
      case KPIType::PRODUCTION_RATE: return "[mol/s]";
      case KPIType::LCOC: return "[£/kg]";
      case KPIType::THERMAL_ENERGY: return "[J/kg]";
      case KPIType::COOLING_DUTY: return "[J/kg]";
      default: return "[-]";
    }
  }

  enum class ObjectiveDirection { MINIMIZE, MAXIMIZE };

  enum class ConstraintSense {
    GREATER_OR_EQUAL,  // actual >= target; penalty when actual < target
    LESS_OR_EQUAL      // actual <= target; penalty when actual > target
  };

  inline constexpr std::string_view kConstraintSenseList[] = { "≥", "≤" };

  struct ObjectiveTerm {
    std::string product_component{};  // which component's KPI this term refers to
    KPIType type{};
    double weight{1.0};
    ObjectiveDirection direction{ObjectiveDirection::MINIMIZE};
  };

  struct ConstraintTerm {
    std::string product_component{};  // which component's KPI this term refers to
    KPIType type{};
    double target{0.0};
    double penalty_scale{1.0};
    ConstraintSense sense{ConstraintSense::GREATER_OR_EQUAL};
  };

  using ObjectiveId = size_t;
  using ConstraintId = size_t;

  // KPI values nested by component then KPI type. For the cost function a
  // non-finite value (NaN/Inf) means "undefined" -- see ComputeCost.
  using PerComponentKPIValues =
    std::unordered_map<std::string, std::unordered_map<KPIType, double>>;

  // Cost reported for a candidate that cannot be scored: an objective or
  // constraint KPI is undefined (e.g. no product drawn), or the solve failed
  // (the desktop optimisation loop uses the same 1e10). Large but FINITE, so
  // dlib's MaxLIPO model stays well-conditioned, and positive so the candidate
  // is never attractive whichever way the objective is directed (the cost is
  // always minimised; a maximise term contributes -w*v). The web optimiser
  // settles the same case as double max (adsim_opt_tell).
  inline constexpr double kOptimisationPenaltyCost = 1e10;

  class OptimizationCostFunction {
  public:
    ObjectiveId AddObjective(const std::string& product_component, KPIType type,
                             double weight, ObjectiveDirection direction) {
      ObjectiveId id = next_objective_id_++;
      objectives_[id] = ObjectiveTerm{product_component, type, weight, direction};
      return id;
    }

    void RemoveObjective(ObjectiveId id) {
      objectives_.erase(id);
    }

    ConstraintId AddConstraint(const std::string& product_component, KPIType type,
                               double target, double penalty_scale, ConstraintSense sense) {
      ConstraintId id = next_constraint_id_++;
      constraints_[id] = ConstraintTerm{product_component, type, target, penalty_scale, sense};
      return id;
    }

    void RemoveConstraint(ConstraintId id) {
      constraints_.erase(id);
    }

    // Scalar cost to MINIMISE: +w*v for a minimise objective, -w*v for a
    // maximise one, plus penalty_scale*violation per constraint. A referenced
    // KPI whose value is non-finite is undefined (kpis.h convention: the KPI
    // returns 0 with defined=false; the caller stores NaN here), and the
    // candidate gets kOptimisationPenaltyCost instead of being scored -- a
    // placeholder 0 would make e.g. "minimise energy" chase the degenerate
    // no-product region. The result is always finite.
    double ComputeCost(const PerComponentKPIValues& kpi_by_component) const {
      auto LookupKPI = [&](const std::string& component, KPIType type) -> const double* {
        auto cit = kpi_by_component.find(component);
        if (cit == kpi_by_component.end()) return nullptr;
        auto kit = cit->second.find(type);
        if (kit == cit->second.end()) return nullptr;
        return &kit->second;
      };

      double cost = 0.0;
      for (const auto& [id, term] : objectives_) {
        (void)id;
        const double* v = LookupKPI(term.product_component, term.type);
        if (!v) continue;
        if (!std::isfinite(*v)) return kOptimisationPenaltyCost;
        if (term.direction == ObjectiveDirection::MINIMIZE) {
          cost += term.weight * (*v);
        } else {
          cost += term.weight * (-(*v));
        }
      }
      for (const auto& [id, term] : constraints_) {
        (void)id;
        const double* v = LookupKPI(term.product_component, term.type);
        if (!v) continue;
        if (!std::isfinite(*v)) return kOptimisationPenaltyCost;
        double violation = (term.sense == ConstraintSense::GREATER_OR_EQUAL)
            ? std::max(0.0, term.target - (*v))
            : std::max(0.0, (*v) - term.target);
        cost += term.penalty_scale * violation;
      }
      return std::isfinite(cost) ? cost : kOptimisationPenaltyCost;
    }

    // Returns the set of distinct components referenced by any enabled term.
    std::vector<std::string> GetReferencedComponents() const {
      std::vector<std::string> out;
      auto Add = [&](const std::string& name) {
        if (std::find(out.begin(), out.end(), name) == out.end())
          out.push_back(name);
      };
      for (const auto& [id, term] : objectives_)  { (void)id; Add(term.product_component); }
      for (const auto& [id, term] : constraints_) { (void)id; Add(term.product_component); }
      return out;
    }

    const std::unordered_map<ObjectiveId, ObjectiveTerm>& GetObjectives() const {
      return objectives_;
    }

    std::unordered_map<ObjectiveId, ObjectiveTerm>& GetObjectives() {
      return objectives_;
    }

    const std::unordered_map<ConstraintId, ConstraintTerm>& GetConstraints() const {
      return constraints_;
    }

    std::unordered_map<ConstraintId, ConstraintTerm>& GetConstraints() {
      return constraints_;
    }

  private:
    std::unordered_map<ObjectiveId, ObjectiveTerm> objectives_{};
    std::unordered_map<ConstraintId, ConstraintTerm> constraints_{};
    ObjectiveId next_objective_id_{1};
    ConstraintId next_constraint_id_{1};
  };

  // A tunable model value. value_ref / expr_param point straight into the
  // model, so an OptimizationParameter is only valid while the object it was
  // enumerated from is alive; see AvailableParameters for the lifetime rules.
  struct OptimizationParameter {
    std::string name;
    std::string category;
    std::string subcategory;
    std::string units;
    // Stable identity assigned by AvailableParameters (empty for parameters
    // built by hand). Unlike value_ref it survives model edits, so callers
    // that hold a selection across edits store keys and re-resolve them with
    // AvailableParameters::Find after re-enumerating.
    std::string key;
    double& value_ref;
    ExpressionParam* expr_param = nullptr;
    double lower_bound;
    double upper_bound;
    double default_value;  // value_ref's value when this entry was built

    // For plain doubles (user variables, mole fractions, etc.)
    OptimizationParameter(
      const std::string& name,
      double& value,
      double lower,
      double upper)
        : name(name),
          category(""),
          subcategory(""),
          units(""),
          value_ref(value),
          expr_param(nullptr),
          lower_bound(lower),
          upper_bound(upper),
          default_value(value) {}

    OptimizationParameter(
      const std::string& name,
      const std::string& category,
      const std::string& subcategory,
      const std::string& units,
      double& value,
      double lower,
      double upper)
        : name(name),
          category(category),
          subcategory(subcategory),
          units(units),
          value_ref(value),
          expr_param(nullptr),
          lower_bound(lower),
          upper_bound(upper),
          default_value(value) {}

    // For ExpressionParam-backed parameters — SetValue keeps expr_str_ in sync
    // so EvaluateAllExpressions can safely re-run without clobbering optimizer changes.
    OptimizationParameter(
      const std::string& name,
      const std::string& category,
      const std::string& subcategory,
      const std::string& units,
      ExpressionParam& param,
      double lower,
      double upper)
        : name(name),
          category(category),
          subcategory(subcategory),
          units(units),
          value_ref(param.Get()),
          expr_param(&param),
          lower_bound(lower),
          upper_bound(upper),
          default_value(param.Get()) {}

    void SetValue(double v) {
      if (expr_param) {
        *expr_param = v;  // updates both value_ and expr_str_
      } else {
        value_ref = v;
      }
    }
  };

  class OptimizationParameters {
  public:
    void Add(const std::string& name, double& value, double lower, double upper) {
      if (Contains(value)) return;
      parameters_.emplace_back(name, value, lower, upper);
    }

    void Add(
      const std::string& name,
      const std::string& category,
      const std::string& subcategory,
      const std::string& units,
      double& value,
      double lower,
      double upper
    ) {
      if (Contains(value)) return;
      parameters_.emplace_back(name, category, subcategory, units, value, lower, upper);
    }

    void Add(const OptimizationParameter& param) {
      if (Contains(param.value_ref)) return;
      parameters_.push_back(param);
    }

    bool Contains(const double& value) const {
      for (const auto& p : parameters_) {
        if (&p.value_ref == &value) return true;
      }
      return false;
    }

    size_t Size() const { return parameters_.size(); }

    OptimizationParameter& operator[](size_t i) { return parameters_[i]; }
    const OptimizationParameter& operator[](size_t i) const {
      return parameters_[i];
    }

    auto begin() { return parameters_.begin(); }
    auto end() { return parameters_.end(); }
    auto begin() const { return parameters_.begin(); }
    auto end() const { return parameters_.end(); }

  private:
    std::vector<OptimizationParameter> parameters_;
  };

  // Settings for GlobalOptimizer. n_init_samples is kept alongside
  // n_iterations so existing call sites can express the same evaluation
  // budget as before: the optimizer performs
  // (n_init_samples + n_iterations) objective evaluations in total.
  struct OptimizerSettings {
    int n_iterations = 30;
    int n_init_samples = 5;
  };

  // Derivative-free global optimizer over the selected parameters.
  // Wraps dlib::find_min_global (MaxLIPO + trust region), which targets
  // expensive black-box objectives with few evaluations. Deterministic for
  // a given objective and bounds.
  //
  // Two ways to give the objective:
  //   * over OptimizationParameters: each candidate is written into the model
  //     (OptimizationParameter::SetValue) and objective_fn() is called;
  //     CurrentPoint() holds the candidate during the call, e.g. to key a
  //     CycleStateCache (adsim/solver/state_cache.h);
  //   * over plain bounds: objective_fn(x) receives the candidate.
  class GlobalOptimizer {
  public:
    GlobalOptimizer(
      OptimizationParameters& params,
      std::function<double()> objective_fn,
      OptimizerSettings settings = {}
    ) : params_ref_(&params),
        objective_function_(std::move(objective_fn)),
        settings_(settings) {
      for (const auto& p : params) {
        lower_.push_back(p.lower_bound);
        upper_.push_back(p.upper_bound);
      }
    }

    // Throws std::invalid_argument when the bounds differ in size.
    GlobalOptimizer(
      std::vector<double> lower,
      std::vector<double> upper,
      std::function<double(const std::vector<double>&)> objective_fn,
      OptimizerSettings settings = {}
    ) : vector_objective_(std::move(objective_fn)),
        settings_(settings),
        lower_(std::move(lower)),
        upper_(std::move(upper)) {
      if (lower_.size() != upper_.size())
        throw std::invalid_argument("GlobalOptimizer: lower and upper bounds differ in size");
    }

    std::vector<double> OptimizeParameters() {
      if (params_ref_) {
        // The parameters and their bounds may have changed since construction.
        lower_.resize(params_ref_->Size());
        upper_.resize(params_ref_->Size());
        for (size_t i = 0; i < params_ref_->Size(); ++i) {
          lower_[i] = (*params_ref_)[i].lower_bound;
          upper_[i] = (*params_ref_)[i].upper_bound;
        }
      }
      const long n = static_cast<long>(lower_.size());
      dlib::matrix<double, 0, 1> lower_bound(n);
      dlib::matrix<double, 0, 1> upper_bound(n);

      for (long i = 0; i < n; ++i) {
        lower_bound(i) = lower_[static_cast<size_t>(i)];
        upper_bound(i) = upper_[static_cast<size_t>(i)];
      }

      const int max_calls =
        std::max(1, settings_.n_init_samples + settings_.n_iterations);

      auto result = dlib::find_min_global(
        [this](const dlib::matrix<double, 0, 1>& x) {
          return EvaluateSample(x);
        },
        lower_bound,
        upper_bound,
        dlib::max_function_calls(static_cast<size_t>(max_calls))
      );

      std::vector<double> optimal_values(static_cast<size_t>(n));
      for (long i = 0; i < n; ++i) {
        optimal_values[static_cast<size_t>(i)] = result.x(i);
      }
      best_value_ = result.y;
      return optimal_values;
    }

    int GetIterationCount() const { return iteration_count_; }
    // The objective at the returned optimum, as the search evaluated it.
    double GetBestValue() const { return best_value_; }
    // The candidate being evaluated (valid during an objective call).
    const std::vector<double>& CurrentPoint() const { return current_; }
    const std::vector<double>& Lower() const { return lower_; }
    const std::vector<double>& Upper() const { return upper_; }

  private:
    double EvaluateSample(const dlib::matrix<double, 0, 1>& x) {
      iteration_count_++;

      current_.resize(static_cast<size_t>(x.size()));
      for (long i = 0; i < x.size(); ++i) current_[static_cast<size_t>(i)] = x(i);
      double cost = 0.0;
      if (params_ref_) {
        for (size_t i = 0; i < params_ref_->Size(); ++i) {
          (*params_ref_)[i].SetValue(current_[i]);
        }
        cost = objective_function_();
      } else {
        cost = vector_objective_(current_);
      }

      // dlib must only ever see finite values: settle a non-finite objective
      // (whatever its source) as the penalty cost.
      return std::isfinite(cost) ? cost : kOptimisationPenaltyCost;
    }

    OptimizationParameters* params_ref_ = nullptr;
    std::function<double()> objective_function_;
    std::function<double(const std::vector<double>&)> vector_objective_;
    OptimizerSettings settings_;
    std::vector<double> lower_;
    std::vector<double> upper_;
    std::vector<double> current_;
    double best_value_ = std::numeric_limits<double>::quiet_NaN();
    int iteration_count_ = 0;
  };

  // ---- optimising a deck -----------------------------------------------------
  //
  // OptimizeSimulation runs the whole loop on a deck: GlobalOptimizer over
  // registry paths, one run per candidate, the cost from the candidate's KPIs
  // at its last cycle (OptimizationCostFunction). Each candidate is run on a
  // fresh copy of the deck, warm-started from the cyclic steady state of the
  // nearest candidate evaluated so far (CycleStateCache, scaled by the
  // variables' bounds); the optimum is then re-run from the deck's initial
  // conditions and that cold value is the reported one. The implementation is
  // in core/src/solver/optimisation.cpp.

  // One decision variable: a tunable registry path and its bounds (SI units).
  struct OptimizationVariable {
    std::string path;
    double lower = 0.0;
    double upper = 0.0;
  };

  struct OptimizationEvaluation {
    int index = 0;                // 0-based evaluation number (-1: the cold re-run)
    std::vector<double> x;        // the candidate
    double cost = kOptimisationPenaltyCost;
    PerComponentKPIValues kpis;   // the referenced components' KPIs (NaN: undefined)
    RunStatus status = RunStatus::NotRun;
    std::string message;
    bool accepted = false;        // the run finished (and reached CSS when required)
    bool warm = false;            // started from a cached state
    double distance = std::numeric_limits<double>::quiet_NaN();  // scaled distance to that state
    std::string note;             // why a warm start was skipped or repeated
    bool css_reached = false;
    int cycles_run = 0;
    double wall_time_s = 0.0;
  };

  struct SimulationOptimizerSettings {
    // GlobalOptimizer's budget: n_init_samples + n_iterations evaluations.
    int n_iterations = 30;
    int n_init_samples = 5;
    // Warm-start each candidate from the nearest cached cyclic steady state.
    bool warm_start = true;
    // Re-run the optimum from the deck's initial conditions (the reported
    // value). Off: the search's own value is reported.
    bool rerun_cold = true;
    size_t cache_capacity = CycleStateCache::kDefaultCapacity;
    // A candidate that stops at the cycle limit without reaching CSS gets the
    // penalty cost.
    bool require_css = false;
    // Overrides the deck's number of cycles for every run (0: the deck's).
    int cycles = 0;
    // Registry text sets applied to the deck before any run (path, text).
    std::vector<std::pair<std::string, std::string>> settings;
    CancelPoll cancel;
    // Called after every evaluation, and after the cold re-run (index -1).
    std::function<void(const OptimizationEvaluation&)> on_evaluation;
  };

  struct SimulationOptimizationResult {
    std::vector<double> x;          // the optimum
    // The optimum's cost: the cold re-run's when rerun_cold, else the search's.
    double cost = kOptimisationPenaltyCost;
    double search_cost = kOptimisationPenaltyCost;  // the optimum's cost as the search evaluated it
    int best_index = -1;            // the evaluation that found the optimum
    bool rerun_cold = false;        // `cold` holds the re-run
    OptimizationEvaluation cold;    // the cold re-run of the optimum
    std::vector<OptimizationEvaluation> evaluations;
    bool cancelled = false;         // stopped early by settings.cancel
  };

  // Checks the set-up without running: the deck loads, every variable is a
  // tunable parameter with finite bounds lower < upper (no path twice), every
  // cost term names a component the deck's KPI configuration tracks and a KPI
  // a run reports (not LCOC), and there is at least one objective. Throws
  // std::invalid_argument with the reason.
  void CheckOptimization(const std::string& deck_json, const std::vector<OptimizationVariable>& variables,
                         const OptimizationCostFunction& cost, const SimulationOptimizerSettings& settings = {});

  // Runs the optimisation (after CheckOptimization). Never throws for a failed
  // run: the candidate gets kOptimisationPenaltyCost. The deck text is not
  // changed; apply result.x to it to keep the optimum.
  SimulationOptimizationResult OptimizeSimulation(const std::string& deck_json,
                                                  const std::vector<OptimizationVariable>& variables,
                                                  const OptimizationCostFunction& cost,
                                                  const SimulationOptimizerSettings& settings = {});

  // The tunable values of a UserData, as optimiser decision variables: the
  // tunable subset of the ParameterRegistry (adsim/parameters/
  // parameter_registry.h), in registry order. Each entry's key is the
  // registry path, so a key survives renames of the step, stream, layer or
  // adsorbent it belongs to.
  //
  // Lifetime: each entry's value_ref / expr_param points into `data` as it was
  // at the last UpdateParameterList call. Structural edits (changing a BC,
  // stream, profile, isotherm or rate-model type; adding or removing steps,
  // streams, layers or user variables; loading or resetting the model) destroy
  // the referenced objects, so entries must not be dereferenced after such an
  // edit. Anything that outlives an edit (a GUI selection, bounds the user
  // typed) should be keyed by OptimizationParameter::key and re-resolved with
  // Find after calling UpdateParameterList again: Find then returns the entry
  // backed by the live object, or nullptr when the parameter no longer exists.
  class AvailableParameters {
  public:
    void UpdateParameterList(UserData& data) {
      parameters_.clear();
      ParameterRegistry registry(data);
      for (Parameter* p : registry.Tunable()) {
        const ParameterInfo& i = p->Info();
        if (ExpressionParam* e = p->ExpressionPointer()) {
          parameters_.emplace_back(i.name, i.category, i.subcategory, i.unit, *e, i.lower, i.upper);
        } else {
          parameters_.emplace_back(i.name, i.category, i.subcategory, i.unit, *p->ValuePointer(),
                                   i.lower, i.upper);
        }
        parameters_.back().key = i.path;
      }
    }

    // Look up an entry by OptimizationParameter::key; nullptr if absent.
    // The pointer is valid until the next UpdateParameterList call.
    OptimizationParameter* Find(std::string_view key) {
      for (auto& p : parameters_)
        if (p.key == key) return &p;
      return nullptr;
    }

    const OptimizationParameter* Find(std::string_view key) const {
      return const_cast<AvailableParameters*>(this)->Find(key);
    }

    // Builds an OptimizationParameters from keys, in the order given. Keys
    // that no longer resolve are skipped. Like Find, the result references the
    // model as of the last UpdateParameterList call.
    OptimizationParameters SelectByKeys(const std::vector<std::string>& keys) {
      OptimizationParameters result;
      for (const auto& key : keys) {
        if (const auto* p = Find(key)) result.Add(*p);
      }
      return result;
    }

    void Print() const {
      std::cout << "==== Available Parameters ====\n";
      for (const auto& param : parameters_) {
        std::cout << param.category << " > " << param.subcategory << " > "
                  << param.name << " [" << param.units << "]\n";
      }
    }

    size_t Size() const { return parameters_.size(); }

    const std::vector<OptimizationParameter>& GetParameters() const { return parameters_; }
    std::vector<OptimizationParameter>& GetParameters() { return parameters_; }

    const OptimizationParameter& operator[](size_t i) const {
      return parameters_[i];
    }

    OptimizationParameters SelectForOptimization(
      const std::vector<size_t>& selected_indices
    ) {
      OptimizationParameters result;
      for (size_t idx : selected_indices) {
        result.Add(parameters_[idx]);
      }
      return result;
    }

    std::vector<OptimizationParameter> FilterByCategory(
      const std::string& category
    ) const {
      std::vector<OptimizationParameter> result;
      for (const auto& p : parameters_) {
        if (p.category == category) {
          result.push_back(p);
        }
      }
      return result;
    }

    auto begin() { return parameters_.begin(); }
    auto end() { return parameters_.end(); }
    auto begin() const { return parameters_.begin(); }
    auto end() const { return parameters_.end(); }

  private:
    std::vector<OptimizationParameter> parameters_{};
  };

} // end adsim namespace
