/*
 * Skarstrom - pressure and temperature swing adsorption simulator
 * Copyright (c) 2024-2026 Sam Affleck and contributors
 *
 * Licensed under the terms in the LICENSE file at the repository root.
 */
#include "../pch.h"
#include "adsim/solver/run_health.h"
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
      double e = CycleMassBalanceError(eq_log, b_log, cycle, c, KPIIntegration::QUADRATURE);
      if (!std::isfinite(e)) e = CycleMassBalanceError(eq_log, b_log, cycle, c, KPIIntegration::TRAPEZOID);
      m.mass_balance_error.push_back(e);
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
