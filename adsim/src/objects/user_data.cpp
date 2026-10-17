/*
 * Skarstrom - pressure and temperature swing adsorption simulator
 * Copyright (c) 2024-2026 Sam Affleck and contributors
 *
 * Licensed under the terms in the LICENSE file at the repository root.
 */
#include "../pch.h"
#include "adsim/objects/user_data.h"
#include "adsim/math/expression_param.h"

#include <cmath>


namespace adsim {

  namespace {
    // Resolve user-variable expressions (each variable may reference others)
    // via fixed-point iteration. Returns true if all converged; false (with
    // err_msg populated) if a cycle is detected or some expression failed to
    // parse after the iteration budget is exhausted.
    bool ResolveVariableExpressions(UserVariableRegistry& reg,
                                    const EvalContext&    base_ctx,
                                    std::string&          err_msg) {
      auto& vars = reg.GetAll();
      const size_t budget = vars.size() + 1; // linear chain converges in N rounds

      bool any_unresolved = true;
      std::string last_err;
      for (size_t round = 0; round < budget && any_unresolved; ++round) {
        any_unresolved = false;
        for (auto& v : vars) {
          if (v.value_expr.empty()) continue; // plain numeric var: nothing to do
          ExpressionParam tmp;
          tmp.SetExpr(v.value_expr);
          if (tmp.Evaluate(base_ctx)) {
            if (tmp.HasInlineUnit() &&
                v.unit_category != UnitCategory::Custom &&
                v.unit_category != UnitCategory::Dimensionless &&
                tmp.InlineUnitCategory() != v.unit_category) {
              // The variable has an explicit, non-default category that
              // contradicts its inline unit — a genuine mismatch, reject it.
              err_msg = "Variable '" + v.name + "': inline unit is "
                      + UnitCategoryName(tmp.InlineUnitCategory())
                      + " but variable category is "
                      + UnitCategoryName(v.unit_category);
              return false;
            }
            // Dimensionless is the "unspecified" default; the inline `[unit]`
            // already sets the value's dimension via tmp.Get(), so we do NOT
            // persist the adopted category onto the variable. Writing it back
            // would overwrite the user's own category (the web exposes a picker)
            // and make a later edit to a different dimension look like a mismatch
            // against the auto-adopted one — matching the frontend rule that
            // Dimensionless adopts the inline unit without changing category.
            v.value = tmp.Get();
          } else {
            any_unresolved = true;
            last_err = "Variable '" + v.name + "': " + tmp.LastError();
          }
        }
      }

      if (any_unresolved) {
        err_msg = last_err.empty()
          ? "Variable dependency cycle or unresolvable expression"
          : last_err;
        return false;
      }
      return true;
    }

    std::string FormatNonFinite(double v) {
      if (std::isnan(v)) return "NaN";
      return v > 0 ? "inf" : "-inf";
    }
  } // anonymous namespace

  std::vector<std::string> UserData::FindNonFiniteExpressions() {
    std::vector<std::string> out;

    // User variables are resolved once with placeholder state (they are never
    // re-evaluated during the run), so a non-finite variable is a non-finite
    // constant even when its expression mentions T, P, c, y or v.
    for (const auto& v : variable_registry.GetAll()) {
      if (std::isfinite(v.value)) continue;
      const std::string expr = v.value_expr.empty() ? FormatNonFinite(v.value) : v.value_expr;
      out.push_back("Variable '" + v.name + "': \"" + expr + "\" = " + FormatNonFinite(v.value));
    }

    // Audit walk: the EvaluateExpressions cascade of EvaluateAllExpressions, in
    // read-only audit mode. Named sub-trees are visited first so their
    // parameters get a specific scope; each container's full cascade is then
    // visited under a coarser scope, which records only what was not yet
    // recorded (Record de-duplicates), so nothing is missed. The exception is
    // a step's own fields, walked one by one (see the step-time note below).
    NonFiniteExpressionAudit audit;
    EvalContext ctx;
    ctx.registry = &variable_registry;
    ctx.audit = &audit;
    auto visit = [&](std::string scope, auto&& fn) {
      audit.scope = std::move(scope);
      fn();
    };

    // Only the step-time inputs of the active mode. In FRACTION mode step_time_
    // is derived (Cycle::ResolveStepTimes, run by Cycle::PerformChecks after
    // this check), so a stale value must not block the check that refreshes it.
    const bool fraction_mode = cycle.step_time_mode_ == StepTimeMode::FRACTION;
    if (fraction_mode)
      visit("Cycle total time", [&] { cycle.total_cycle_time_.Evaluate(ctx); });
    for (auto& step : cycle.step_list_) {
      const std::string s = "Step '" + step.step_name_ + "'";
      if (fraction_mode)
        visit(s + " step fraction", [&] { step.step_fraction_.Evaluate(ctx); });
      else
        visit(s + " step time", [&] { step.step_time_.Evaluate(ctx); });
      // The lambdas passed to visit() below only ever capture named locals, never
      // a structured binding: Apple clang before 16 rejects capturing a binding.
      for (auto& entry : step.energy_sources_) {
        auto& src = entry.second;
        visit(s + " energy source", [&] { src.EvaluateExpressions(ctx); });
      }
      const std::pair<const char*, IBoundaryCondition*> sides[] = {
        {"top", step.top_boundary_condition_.get()},
        {"bottom", step.bottom_boundary_condition_.get()}};
      for (const auto& side_bc : sides) {
        const char* side = side_bc.first;
        IBoundaryCondition* bc = side_bc.second;
        if (!bc) continue;
        if (auto* fbc = dynamic_cast<FlowBoundaryCondition*>(bc)) {
          for (auto& stream : fbc->streams_) {
            if (!stream) continue;
            visit(s + " " + side + " stream '" + stream->name_ + "'",
                  [&] { stream->EvaluateExpressions(ctx); });
          }
        } else if (auto* pbc = dynamic_cast<PressureBoundaryCondition*>(bc)) {
          if (pbc->pressure_)
            visit(s + " " + side + " pressure", [&] { pbc->pressure_->EvaluateExpressions(ctx); });
          visit(s + " " + side + " inlet temperature", [&] { pbc->temperature_.Evaluate(ctx); });
        }
        visit(s + " " + side + " boundary", [&] { bc->EvaluateExpressions(ctx); });
      }
    }

    visit("Bed ambient temperature", [&] { bed.ambient_temperature_.Evaluate(ctx); });
    visit("Bed geometry", [&] { EvaluateBedShapeExpressions(bed.shape_, ctx); });
    visit("Bed wall", [&] { bed.wall_.EvaluateExpressions(ctx); });
    visit("Bed bottom void", [&] { bed.bottom_void_.EvaluateExpressions(ctx); });
    visit("Bed top void", [&] { bed.top_void_.EvaluateExpressions(ctx); });
    for (auto& layer : bed.layers_.layers_) {
      const std::string l = "Layer '" + layer.name_ + "'";
      visit(l + " bed void fraction", [&] { layer.eb_.Evaluate(ctx); });
      visit(l + " particle diameter", [&] { layer.particle_diameter_.Evaluate(ctx); });
      visit(l + " length", [&] { layer.domain_.EvaluateExpressions(ctx); });
      visit(l + " transport model", [&] { layer.EvaluateExpressions(ctx); });
    }
    visit("Bed heat-transfer model", [&] { bed.EvaluateExpressions(ctx); });

    for (auto& entry : const_cast<std::unordered_map<AdsorbentID, Adsorbent>&>(
             adsorbent_directory.GetAll())) {
      Adsorbent& ads = entry.second;
      const std::string a = "Adsorbent '" + ads.name_ + "'";
      visit(a + " solid density", [&] { ads.solid_density_.Evaluate(ctx); });
      visit(a + " heat capacity", [&] { ads.heat_capacity_.Evaluate(ctx); });
      visit(a + " thermal conductivity", [&] { ads.thermal_conductivity_.Evaluate(ctx); });
      visit(a + " particle void fraction", [&] { ads.ep_.Evaluate(ctx); });
      auto& src = ads.adsorption_source_;
      for (auto& [component, dH_entry] : src.heat_of_adsorption_) {
        auto& dH = dH_entry;
        visit(a + " " + component + " heat of adsorption", [&] { dH.Evaluate(ctx); });
      }
      for (auto& [component, mtm_entry] : src.mass_transfer_models_) {
        auto& mtm = mtm_entry;
        if (mtm) visit(a + " " + component + " mass transfer", [&] { mtm->EvaluateExpressions(ctx); });
      }
      visit(a + " isotherm", [&] { ads.EvaluateExpressions(ctx); });
    }

    visit("Initial temperature", [&] { ic.temperature.Evaluate(ctx); });
    visit("Initial wall temperature", [&] { ic.wall_temperature.Evaluate(ctx); });
    visit("Initial pressure", [&] { ic.pressure.Evaluate(ctx); });
    visit("Initial velocity", [&] { ic.velocity.Evaluate(ctx); });
    visit("Initial conditions", [&] { ic.EvaluateExpressions(ctx); });
    visit("Solver settings", [&] { settings.EvaluateExpressions(ctx); });
    for (auto& [component, cp_entry] : fluid.adsorbed_heat_capacity_) {
      auto& cp = cp_entry;
      visit(component + " adsorbed-phase heat capacity", [&] { cp.Evaluate(ctx); });
    }
    visit("Fluid", [&] { fluid.EvaluateExpressions(ctx); });

    for (const auto& e : audit.entries)
      out.push_back(e.scope + ": \"" + e.expr + "\" = " + FormatNonFinite(e.value));
    return out;
  }

  bool UserData::PerformChecks(std::string& err_msg) {
    // First, so a deck with an inf/NaN parameter is reported as such rather
    // than tripping a downstream check in some less obvious way (NaN fails
    // every comparison).
    const auto non_finite = FindNonFiniteExpressions();
    if (!non_finite.empty()) {
      err_msg = "Non-finite parameter value (inf/NaN, e.g. a division by zero): ";
      for (size_t i = 0; i < non_finite.size(); ++i)
        err_msg += (i ? "; " : "") + non_finite[i];
      return false;
    }

    bool checks_passed = fluid.PerformChecks(err_msg);
    if (!checks_passed) {
      return false;
    }

    checks_passed = bed.PerformChecks(adsorbent_directory, err_msg);
    if (!checks_passed) {
      return false;
    }

    checks_passed = ic.PerformChecks(fluid, err_msg);
    if (!checks_passed) {
      return false;
    }

    checks_passed = cycle.PerformChecks(fluid, err_msg);
    if (!checks_passed) {
      return false;
    }

    return true;
  }
  
  void UserData::Initialise() {
    fluid.Initialise(ic);
    bed.Initialise(adsorbent_directory, fluid);
    state.Initialise(fluid, bed, adsorbent_directory);

    // Parameters the residual evaluates per cell are compiled once per run.
    for (const auto& [id, ads] : adsorbent_directory.GetAll())
      adsorbent_directory.Get(id).adsorption_source_.CompileForRun(&variable_registry);

    int N_CELLS = bed.GetNumberOfCells();
    state.Resize(N_CELLS);

    // Set T, Tw, P and yi with user defined initial conditions
    state.SetInitialConditions(ic, bed.equations_.thermal_equilibrium_mode_);
    // Isothermal run: pin fluid, solid and wall temperature at T_iso.
    state.HoldIsothermalTemperature(bed.equations_.GetIsothermalTemperature(ic.temperature));

    // Update the isotherms to get qi_sat 
    state.UpdateIndependentEquations(this);
    state.UpdateDiagnosticEquations(this);

    // Set qi = qi_sat 
    state.SetAdsorbedPhaseInitialConditions();

    // Initialize step start time to 0 for first step
    step_start_time_ = 0.0;
  }

  void UserData::RemoveAdsorbent(AdsorbentID id) {
    adsorbent_directory.Remove(id);
    for (auto& layer : bed.layers_.layers_) layer.RemoveAdsorbent(id);
    cost_estimator.adsorbent_cost_per_kg.erase(id);
  }

  size_t UserData::DuplicateStep(size_t index) {
    const size_t new_index = cycle.DuplicateStep(index);
    if (new_index != static_cast<size_t>(-1))
      kpi_config.OnStepInserted(static_cast<int>(new_index));
    return new_index;
  }

  void UserData::RemoveStep(size_t index) {
    auto& steps = cycle.step_list_;
    if (index >= steps.size()) return;
    // By index rather than Cycle::RemoveStep(name) so a duplicated name can
    // never remove a different step than the one whose refs are dropped.
    DisconnectBoundaryCondition(steps, steps[index].top_boundary_condition_.get());
    DisconnectBoundaryCondition(steps, steps[index].bottom_boundary_condition_.get());
    steps.erase(steps.begin() + static_cast<std::ptrdiff_t>(index));
    kpi_config.OnStepRemoved(static_cast<int>(index));
  }

  void UserData::SwapSteps(size_t a, size_t b) {
    auto& steps = cycle.step_list_;
    if (a >= steps.size() || b >= steps.size() || a == b) return;
    std::swap(steps[a], steps[b]);
    kpi_config.OnStepsSwapped(static_cast<int>(a), static_cast<int>(b));
  }

  void UserData::RenameStep(size_t index, const std::string& new_name) {
    auto& steps = cycle.step_list_;
    if (index >= steps.size()) return;
    const std::string old_name = steps[index].step_name_;
    if (old_name == new_name) return;
    // A buffer-tank source resolves its step by name to the FIRST step with
    // that name (Cycle::FindStep), so it follows the rename only if that is
    // this step (PerformChecks rejects duplicate step names anyway).
    const bool tanks_follow = cycle.FindStep(old_name) == &steps[index];
    steps[index].step_name_ = new_name;
    kpi_config.OnStepRenamed(static_cast<int>(index), old_name, new_name);
    if (tanks_follow) cycle.RenameBufferTankSourceStep(old_name, new_name);
  }

  void UserData::SwapStepBoundaries(size_t index) {
    auto& steps = cycle.step_list_;
    if (index >= steps.size()) return;
    Step& step = steps[index];
    std::swap(step.top_boundary_condition_, step.bottom_boundary_condition_);
    // location_ is otherwise only re-stamped when the GUI draws a side's tab
    // (or at Step::Initialise); until then GetEndpoints / Prune and the
    // "same location" check in Cycle::PerformChecks would see the old sides.
    if (step.top_boundary_condition_)
      step.top_boundary_condition_->location_ = BoundaryLocation::TOP;
    if (step.bottom_boundary_condition_)
      step.bottom_boundary_condition_->location_ = BoundaryLocation::BOTTOM;
    kpi_config.OnStepBoundariesSwapped(static_cast<int>(index), step.step_name_);
    // Buffer-tank sources name an endpoint by stream name or the side-less
    // kPressureBoundaryEndpointName, so the swap leaves them valid.
  }

  void UserData::RenameStream(size_t step_index, BoundaryLocation side, size_t stream_index,
                              const std::string& new_name) {
    auto& steps = cycle.step_list_;
    if (step_index >= steps.size()) return;
    Step& step = steps[step_index];
    auto& slot = (side == BoundaryLocation::TOP) ? step.top_boundary_condition_
                                                 : step.bottom_boundary_condition_;
    auto* fbc = dynamic_cast<FlowBoundaryCondition*>(slot.get());
    if (!fbc || stream_index >= fbc->streams_.size() || !fbc->streams_[stream_index]) return;
    IBoundaryStream& stream = *fbc->streams_[stream_index];
    const std::string old_name = stream.name_;
    if (old_name == new_name) return;

    // A buffer-tank source resolves its endpoint to the FIRST same-named one
    // on the source step (Step::FindSummaryForEndpoint), so it follows the
    // rename only if that is this stream.
    const bool tanks_follow = step.FindSummaryForEndpoint(old_name) == &stream.cycle_summary_;
    stream.name_ = new_name;

    // KPI refs are keyed by (step, name) and cache the endpoint's bed end, and
    // stream names need not be unique within a step. If no other endpoint on
    // the step still carries old_name, every ref to it follows the stream
    // (Prune would drop it otherwise). If one does at the OTHER end, only the
    // refs at this stream's end follow. If one does at the SAME end the refs
    // cannot be told apart, so they stay with the endpoint that kept the name.
    bool kept_here = false, kept_other_end = false;
    for (const auto& ep : step.GetEndpoints()) {
      if (ep.name != old_name) continue;
      (ep.location == side ? kept_here : kept_other_end) = true;
    }
    if (!kept_here) {
      std::optional<BoundaryLocation> only_at;
      if (kept_other_end) only_at = side;
      kpi_config.OnStreamRenamed(static_cast<int>(step_index), old_name, new_name, only_at);
    }
    if (tanks_follow) cycle.RenameBufferTankSourceEndpoint(step.step_name_, old_name, new_name);
  }

  bool UserData::EvaluateAllExpressions(std::string& err_msg) {
    EvalContext ctx;
    ctx.registry = &variable_registry;

    // 1. Resolve user-variable expressions in dependency order
    if (!ResolveVariableExpressions(variable_registry, ctx, err_msg)) {
      return false;
    }

    // 2. Cycle (timing + cascade to BCs in every step).
    cycle.EvaluateExpressions(ctx);
    for (auto& step : cycle.step_list_) {
      if (step.top_boundary_condition_)
        step.top_boundary_condition_->EvaluateExpressions(ctx);
      if (step.bottom_boundary_condition_)
        step.bottom_boundary_condition_->EvaluateExpressions(ctx);
    }

    // 3. Bed (shape, wall, layers, ambient temperature, layer equations).
    bed.EvaluateExpressions(ctx);

    // 4. Adsorbent directory (shared material library).
    for (auto& [id, ads] : const_cast<std::unordered_map<AdsorbentID, Adsorbent>&>(
             adsorbent_directory.GetAll())) {
      ads.EvaluateExpressions(ctx);
    }

    // 5. Initial conditions and solver settings.
    ic.EvaluateExpressions(ctx);
    settings.EvaluateExpressions(ctx);

    // 6. Fluid-owned expressions (per-component adsorbed heat capacities).
    fluid.EvaluateExpressions(ctx);

    return true;
  }

} // end adsim namespace
