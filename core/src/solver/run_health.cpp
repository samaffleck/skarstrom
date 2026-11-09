/*
 * Skarstrom - pressure and temperature swing adsorption simulator
 * Copyright (c) 2024-2026 Sam Affleck and contributors
 *
 * Licensed under the terms in the LICENSE file at the repository root.
 */
#include "../pch.h"
#include "adsim/solver/run_health.h"
#include "adsim/equations/reaction_kinds.h"
#include "adsim/objects/reaction.h"
#include "adsim/objects/user_data.h"
#include "adsim/solver/coupled_solver.h"
#include "adsim/solver/cyclic_steady_state.h"
#include "adsim/solver/kpi_history.h"
#include "adsim/solver/kpis.h"
#include "adsim/solver/log.h"

#include <cmath>
#include <iomanip>
#include <sstream>

namespace adsim {

  namespace {

    using T = RunHealthThresholds;

    std::string Num(double v, int precision = 3) {
      std::ostringstream ss;
      ss << std::setprecision(precision) << v;
      return ss.str();
    }

    std::string SchemeName(AdvectionScheme s) {
      switch (s) {
        case AdvectionScheme::Upwind: return "upwind";
        case AdvectionScheme::WENO:   return "weno";
        case AdvectionScheme::VanLeer: return "van_leer";
      }
      return "?";
    }

    std::string StepName(const UserData& data, int step) {
      const auto& steps = data.cycle.step_list_;
      return (step >= 0 && step < static_cast<int>(steps.size())) ? steps[step].step_name_ : std::string{};
    }

    int TotalCells(const UserData& data) {
      int n = 0;
      for (const auto& layer : data.bed.layers_.layers_) n += layer.domain_.n_cells;
      return n;
    }

    bool Enabled(const WarningSettings& on, RunHealthCheck check) {
      switch (check) {
        case RunHealthCheck::MoleFractionSum:  return on.mole_fraction_sum;
        case RunHealthCheck::MassBalanceStuck: return on.mass_balance_stuck;
        case RunHealthCheck::StepRetries:      return on.step_retries;
        case RunHealthCheck::HistoryTruncated: return on.history_truncated;
        case RunHealthCheck::ReactionStiff:    return on.reaction_stiff;
      }
      return true;
    }

    // Largest |sum y - 1| of one logged cycle over the column's cells (packed
    // and void, no ghosts) and the step it is in.
    void MoleFractionSum(const EquationLog& log, int cycle, RunHealthCycle& out) {
      const auto& sim_log = log.GetSimulationLog();
      const auto it = sim_log.find(cycle);
      if (it == sim_log.end()) return;
      const auto& dom = log.GetDomain();
      std::vector<int> cells = dom.bed_cell_indices;
      if (dom.bottom_void_cell_index > 0) cells.push_back(dom.bottom_void_cell_index);
      if (dom.top_void_cell_index > 0) cells.push_back(dom.top_void_cell_index);
      double worst = std::numeric_limits<double>::quiet_NaN();
      for (const auto& [step, entry] : it->second) {
        const auto v = entry.variables.find(VariableType::TOTAL_MOLE_FRACTION);
        if (v == entry.variables.end()) continue;
        for (const auto& row : v->second.values)
          for (int i : cells) {
            if (i < 0 || i >= static_cast<int>(row.size())) continue;
            const double e = std::abs(row[i] - 1.0);
            if (!std::isfinite(e)) continue;
            if (std::isnan(worst) || e > worst) { worst = e; out.mole_fraction_sum_step = step; }
          }
      }
      out.mole_fraction_sum_error = worst;
    }

    // REACTION_STIFF: which reaction and component of `cell` to name (see the
    // header comment). False when no reaction is active in the cell.
    struct ReactionSuspect {
      std::string reaction, solid, component;
      double y = 0.0;          // the component's smallest mole fraction over the reaction's cells
      int cell = -1;           // where
      double exponent = 0.0;   // its smallest power in the reaction's rate
    };
    bool SuspectReaction(const UserData& data, int cell, const std::string& error_component,
                         ReactionSuspect& out) {
      const auto& state = data.state;
      if (cell < 0 || cell >= static_cast<int>(state.domain_vectors.cell_adsorbents.size())) return false;
      const auto& here = state.domain_vectors.cell_adsorbents[cell];
      struct Candidate { const Reaction* r; AdsorbentID solid_id; std::string solid; double extent; std::string c; double y, m; };
      std::vector<Candidate> all;
      for (const auto& field : state.reaction_fields) {
        bool present = false;
        for (const auto& ca : here)
          if (ca.id == field->adsorbent) { present = true; break; }
        if (!present) continue;
        const Adsorbent& ads = data.adsorbent_directory.Get(field->adsorbent);
        const Reaction* r = ads.reaction_source_.Find(field->reaction);
        if (!r || !r->rate_) continue;
        const auto& ext = field->extent_rate.var;
        const double extent = cell < static_cast<int>(ext.size()) ? std::abs(ext[cell]) : 0.0;
        for (const auto& c : data.fluid.components_) {
          const double m = r->rate_->MinExponent(c, r->nu_);
          if (!std::isfinite(m) && r->nu_.count(c) == 0) continue;
          auto it = state.mole_fraction.find(c);
          if (it == state.mole_fraction.end() || cell >= static_cast<int>(it->second.var.size())) continue;
          all.push_back({r, field->adsorbent, ads.name_, extent, c, it->second.var[cell], m});
        }
      }
      if (all.empty()) return false;
      // The component: the error's own if a power below 1 involves it, else
      // the smallest mole fraction among the powers below 1, else among all.
      auto pick_component = [&](auto keep) -> std::string {
        const Candidate* best = nullptr;
        for (const auto& k : all)
          if (keep(k) && (!best || k.y < best->y)) best = &k;
        return best ? best->c : std::string{};
      };
      std::string c;
      for (const auto& k : all)
        if (k.c == error_component && k.m < 1.0) c = k.c;
      if (c.empty()) c = pick_component([](const Candidate& k) { return k.m < 1.0; });
      if (c.empty()) c = pick_component([](const Candidate&) { return true; });
      // The reaction: of those involving c, the largest rate in the cell.
      const Candidate* best = nullptr;
      for (const auto& k : all)
        if (k.c == c && (!best || k.extent > best->extent)) best = &k;
      out.reaction = best->r->name_;
      out.solid = best->solid;
      out.component = c;
      out.exponent = best->m;
      // Where it is nearest zero: the smallest mole fraction over the cells
      // that carry the reaction (the front ahead of the largest error).
      out.y = best->y;
      out.cell = cell;
      const AdsorbentID solid_id = best->solid_id;
      const auto& y = state.mole_fraction.at(c).var;
      const auto& cells = state.domain_vectors.cell_adsorbents;
      for (int n = 0; n < static_cast<int>(std::min(y.size(), cells.size())); ++n) {
        bool present = false;
        for (const auto& ca : cells[n])
          if (ca.id == solid_id) { present = true; break; }
        if (present && y[n] < out.y) { out.y = y[n]; out.cell = n; }
      }
      out.y = std::max(out.y, 0.0);  // an iterate can sit a round-off below 0
      return true;
    }

  } // namespace

  const char* RunHealthCode(RunHealthCheck check) {
    return kRunHealthCodeList[static_cast<int>(check)];
  }

  void RunHealthMonitor::Begin(const UserData& data) {
    report_ = {};
    report_.components = data.fluid.components_;
    finished_ = false;
  }

  void RunHealthMonitor::Offer(RunHealthWarning w) {
    w.code = RunHealthCode(w.check);
    for (auto& existing : report_.warnings) {
      if (existing.check != w.check) continue;
      if (w.value > existing.value) existing = std::move(w);
      return;
    }
    report_.warnings.push_back(std::move(w));
  }

  void RunHealthMonitor::OnCycle(const UserData& data, const EquationLog& eq_log, const BoundaryLog& b_log,
                                 int cycle, const IntegrationStats& at_start, const IntegrationStats& at_end) {
    RunHealthCycle m;
    m.cycle = cycle;
    MoleFractionSum(eq_log, cycle, m);
    const double change = CyclicSteadyState::EndOfCycleStateChange(eq_log, cycle);
    m.state_change = change >= 0.0 ? change : std::numeric_limits<double>::quiet_NaN();
    for (const auto& c : report_.components) {
      // The reported KPIs' balance: the quadratures, each endpoint falling
      // back to the trapezoid when it has none (KPIIntegration::QUADRATURE).
      m.mass_balance_error.push_back(CycleMassBalanceError(eq_log, b_log, cycle, c));
    }
    m.steps = at_end.num_steps - at_start.num_steps;
    m.rejected = (at_end.num_err_test_fails - at_start.num_err_test_fails)
               + (at_end.num_conv_fails - at_start.num_conv_fails)
               + (at_end.recovered_steps - at_start.recovered_steps);
    report_.cycles.push_back(m);

    const WarningSettings& on = data.settings.warnings_;

    // (a) Sum of the mole fractions.
    if (Enabled(on, RunHealthCheck::MoleFractionSum) && m.mole_fraction_sum_error > T::kMoleFractionSum) {
      RunHealthWarning w;
      w.check = RunHealthCheck::MoleFractionSum;
      w.value = m.mole_fraction_sum_error;
      w.threshold = T::kMoleFractionSum;
      w.cycle = cycle;
      w.step = m.mole_fraction_sum_step;
      w.step_name = StepName(data, w.step);
      const std::string scheme = SchemeName(data.settings.advection_scheme_);
      w.suspects = { "advection scheme (settings.advection_scheme = " + scheme + ")",
                     "grid (" + std::to_string(TotalCells(data)) + " cells)",
                     "step \"" + w.step_name + "\" and its outlet" };
      w.hint = "Re-run with settings.advection_scheme = " + std::string(scheme == "upwind" ? "weno" : "upwind")
               + " and with more cells to compare; a sum that moves with them is numerical, one that does not is "
                 "in the model.";
      w.message = "Warning: the gas mole fractions sum to 1 +/- " + Num(w.value) + " in step \"" + w.step_name
                  + "\" of cycle " + std::to_string(cycle + 1) + " (more than " + Num(w.threshold)
                  + "): the composition is not conserved. " + w.hint;
      Offer(std::move(w));
    }

    // (b) A cycle mole balance that stays open while the state has settled.
    const int W = T::kStuckWindow;
    const int n = static_cast<int>(report_.cycles.size());
    if (Enabled(on, RunHealthCheck::MassBalanceStuck) && n >= W && m.state_change <= T::kStuckStateChange) {
      for (size_t c = 0; c < report_.components.size(); ++c) {
        bool stuck = true;
        for (int k = n - W; k < n && stuck; ++k) {
          const auto& row = report_.cycles[k].mass_balance_error;
          stuck = c < row.size() && std::abs(row[c]) > T::kStuckBalance;
        }
        if (!stuck) continue;
        const double first = std::abs(report_.cycles[n - W].mass_balance_error[c]);
        const double last = std::abs(m.mass_balance_error[c]);
        const double rate = std::pow(last / first, 1.0 / (W - 1));
        if (!(rate >= T::kStuckRate)) continue;
        RunHealthWarning w;
        w.check = RunHealthCheck::MassBalanceStuck;
        w.value = last;
        w.threshold = T::kStuckBalance;
        w.cycle = cycle;
        w.component = report_.components[c];
        w.suspects = { "advection scheme (settings.advection_scheme = "
                           + SchemeName(data.settings.advection_scheme_) + ")",
                       "boundary streams and connections that move " + w.component,
                       "reactions or sources of " + w.component };
        w.hint = "A balance that does not close as the run settles is a model or numerics error, not slow CSS: "
                 "check every endpoint that carries " + w.component + ", then re-run with "
                 "settings.advection_scheme = upwind and more cells to compare.";
        w.message = "Warning: the cycle mole balance of " + w.component + " stays at " + Num(m.mass_balance_error[c])
                    + " in cycle " + std::to_string(cycle + 1) + " (above " + Num(w.threshold) + " for "
                    + std::to_string(W) + " cycles, falling by less than "
                    + Num(100.0 * (1.0 - T::kStuckRate), 2) + " % per cycle) while the end-of-cycle state changes by "
                    + Num(m.state_change) + ": moles are gained or lost. " + w.hint;
        Offer(std::move(w));
      }
    }

    // (c) Rejected integrator attempts.
    if (Enabled(on, RunHealthCheck::StepRetries) && m.steps >= T::kStepRetryMinSteps) {
      const double ratio = static_cast<double>(m.rejected) / static_cast<double>(m.steps);
      if (ratio > T::kStepRetryRatio) {
        RunHealthWarning w;
        w.check = RunHealthCheck::StepRetries;
        w.value = ratio;
        w.threshold = T::kStepRetryRatio;
        w.cycle = cycle;
        w.suspects = { "solver tolerances (settings.reltol, settings.abstol)",
                       "abrupt boundary conditions (constant pressure steps, valve or flow switches)",
                       "grid (" + std::to_string(TotalCells(data)) + " cells)" };
        w.hint = "Check the tolerances first (settings.reltol 1e-6 is the default), then smooth abrupt boundary "
                 "values with a ramp; run at log level high to see where the integrator fails.";
        w.message = "Warning: the integrator rejected " + std::to_string(m.rejected) + " attempts for "
                    + std::to_string(m.steps) + " accepted steps in cycle " + std::to_string(cycle + 1) + " ("
                    + Num(ratio) + " per step, more than " + Num(w.threshold) + "). " + w.hint;
        Offer(std::move(w));
      }
    }
  }

  void RunHealthMonitor::OnStepStart() {
    short_steps_ = 0;
    reaction_named_ = false;
  }

  void RunHealthMonitor::OnShortStep(const UserData& data, int cycle, int step, double t, double h, int cell,
                                     const std::string& component, bool failed) {
    if (!Enabled(data.settings.warnings_, RunHealthCheck::ReactionStiff)) return;
    if (!(std::abs(h) < T::kReactionStiffStep)) return;
    ReactionSuspect suspect;
    if (!SuspectReaction(data, cell, component, suspect)) return;
    if (!failed) ++short_steps_;
    if (reaction_named_) {
      // Named at the threshold (the suspect as the collapse sets in); later
      // short steps of the same cycle step only raise its count.
      for (auto& w : report_.warnings)
        if (w.check == RunHealthCheck::ReactionStiff && w.cycle == cycle && w.step == step)
          w.value = std::max(w.value, static_cast<double>(short_steps_));
      return;
    }
    if (!failed && short_steps_ < T::kReactionStiffMinSteps) return;
    reaction_named_ = true;

    RunHealthWarning w;
    w.check = RunHealthCheck::ReactionStiff;
    w.value = static_cast<double>(short_steps_);
    w.threshold = static_cast<double>(T::kReactionStiffMinSteps);
    w.cycle = cycle;
    w.step = step;
    w.step_name = StepName(data, step);
    w.component = suspect.component;
    w.reaction = suspect.reaction;
    const std::string power = suspect.exponent < 0.0 ? "a negative power (" + Num(suspect.exponent) + ")"
                            : suspect.exponent < 1.0 ? "a power below 1 (" + Num(suspect.exponent) + ")"
                                                     : "power " + Num(suspect.exponent);
    w.suspects = { "reaction \"" + suspect.reaction + "\" on \"" + suspect.solid + "\" (its rate raises "
                       + suspect.component + " to " + power + ")",
                   "the mole fraction of " + suspect.component + " in the reacting cells (down to "
                       + Num(suspect.y) + ", cell " + std::to_string(suspect.cell) + ")",
                   "the initial gas composition (initial.y) and step \"" + w.step_name + "\"" };
    w.hint = "A reaction that relaxes much faster than the flow holds the step at its own time scale, typically "
             "where the gas holds almost no " + suspect.component + " yet. Start the bed with some "
             + suspect.component + " in initial.y (or at the gas-phase equilibrium of the feed), add a trace to "
               "the feed, or check the rate constants and orders of \"" + suspect.reaction
             + "\"; see docs/reactions.md, \"Start-up of fast reactions\".";
    std::ostringstream msg;
    msg << "Warning: the integrator's step collapsed at a reacting cell in step \"" << w.step_name << "\" of cycle "
        << cycle + 1 << " (" << short_steps_ << " steps shorter than " << Num(T::kReactionStiffStep) << " s by t = "
        << Num(t, 4) << " s" << (failed ? ", then the step failed" : "") << "): reaction \"" << suspect.reaction
        << "\" with " << suspect.component << " down to y = " << Num(suspect.y) << " in the reacting cells. "
        << w.hint;
    w.message = msg.str();
    Offer(std::move(w));
  }

  void RunHealthMonitor::Finish(const UserData& data, const EquationLog& eq_log) {
    if (finished_) return;
    finished_ = true;
    // (d) The log kept fewer cycles than the diagnostics read back.
    const int run = static_cast<int>(report_.cycles.size());
    const int kept = static_cast<int>(eq_log.GetSimulationLog().size());
    const int needed = std::min(run, T::kHistoryCycles);
    if (Enabled(data.settings.warnings_, RunHealthCheck::HistoryTruncated) && kept < needed) {
      RunHealthWarning w;
      w.check = RunHealthCheck::HistoryTruncated;
      w.value = kept;
      w.threshold = needed;
      w.suspects = { "settings.max_cycles_to_keep = " + std::to_string(data.settings.max_cycles_to_keep_) };
      w.hint = "Set settings.max_cycles_to_keep to at least " + std::to_string(T::kHistoryCycles)
               + " (default 10) to plot or inspect more of the run. KPIs and the CSS history are kept for every "
                 "cycle; only the logged series and profiles are limited.";
      w.message = "Warning: the log kept " + std::to_string(kept) + " of the " + std::to_string(run)
                  + " cycles run, fewer than " + std::to_string(needed)
                  + ": series and profiles read after the run cover only the logged cycles. " + w.hint;
      Offer(std::move(w));
    }
  }

} // namespace adsim
