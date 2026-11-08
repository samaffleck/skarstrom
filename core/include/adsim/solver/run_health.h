/*
 * Skarstrom - pressure and temperature swing adsorption simulator
 * Copyright (c) 2024-2026 Sam Affleck and contributors
 *
 * Licensed under the terms in the LICENSE file at the repository root.
 */
#pragma once

// Run-health checks: cheap tests after every cycle of a run for symptoms of a
// run that completes but should not be trusted (a model or numerics error, a
// grid or tolerance the deck cannot carry, history the diagnostics need that
// was not kept). They never stop a run or change its status or results; they
// add structured warnings (RunResult::health_warnings, RunResult::healthy) and
// the solver logs each once at the end of the run (LogLevel::LOW).
//
//  * MOLE_FRACTION_SUM: the largest |sum_i y_i - 1| over the cycle's logged
//    times and the column's cells (packed and void cells; the logged total
//    mole fraction). Every component balance is solved for its own mole
//    fraction, so nothing forces the sum to 1: it drifts when a scheme or a
//    source does not conserve the gas composition.
//  * MASS_BALANCE_STUCK: the cycle mole balance of a component
//    (CycleMassBalanceError, kpi_history.h; the quadrature balance, else the
//    trapezoid one) is large and has stopped falling while the end-of-cycle
//    state change (CyclicSteadyState::EndOfCycleStateChange) is already
//    small. The balance compares the moles through the endpoints with the
//    change of the column inventory, so a conservative run closes it on every
//    cycle up to the integration error; one that stays open as the run
//    settles is a model error, not slow convergence. The rule, at cycle k:
//    the state change at k is at most kStuckStateChange, and over the
//    kStuckWindow cycles ending at k every |balance| exceeds kStuckBalance
//    and the balance fell by less than (1 - kStuckRate) per cycle:
//    (|e_k| / |e_(k-W+1)|)^(1 / (W - 1)) >= kStuckRate. Needs the KPI
//    quadratures (IDASettings::kpi_quadratures_, on by default): without
//    them the column inventory is not recorded and the check is skipped.
//  * STEP_RETRIES: the integrator rejected too many attempts in a cycle:
//    (error-test failures + nonlinear convergence failures + the solver's own
//    step recoveries, IntegrationStats::recovered_steps) / accepted internal
//    steps above kStepRetryRatio, over a cycle of at least kStepRetryMinSteps
//    accepted steps.
//  * HISTORY_TRUNCATED: the run outlived its log (settings.max_cycles_to_keep)
//    and the log kept fewer cycles than the diagnostics read back
//    (kHistoryCycles: the CSS distance estimate's window, and this check's
//    own balance window), so the series and profiles read after the run cover
//    only the logged cycles. KPIs and the CSS history are kept for every cycle.
//    The checks above run as each cycle completes and are not affected.
//
// Each code gives at most one warning per run, at the cycle where its value
// was largest. The thresholds come from the engine's own decks (every example
// deck and paper case is silent with margin); docs/run-health.md lists the
// measured distribution. Each warning has a switch, settings.warn_<code in
// lower case> (warning_settings.h); a switched-off check is not reported and
// does not make the run unhealthy.

#include <limits>
#include <string>
#include <vector>

namespace adsim {

  class UserData;
  class EquationLog;
  struct BoundaryLog;
  struct IntegrationStats;

  enum class RunHealthCheck {
    MoleFractionSum,
    MassBalanceStuck,
    StepRetries,
    HistoryTruncated
  };

  // The codes, in RunHealthCheck order.
  inline constexpr const char* kRunHealthCodeList[] = {
    "MOLE_FRACTION_SUM",
    "MASS_BALANCE_STUCK",
    "STEP_RETRIES",
    "HISTORY_TRUNCATED"
  };
  const char* RunHealthCode(RunHealthCheck check);

  // Thresholds (see the header comment and docs/run-health.md).
  struct RunHealthThresholds {
    static constexpr double kMoleFractionSum = 1e-3;      // max |sum y - 1| [-]
    static constexpr int kStuckWindow = 5;                // cycles
    static constexpr double kStuckBalance = 2e-3;         // |cycle mole balance| [mol/mol]
    static constexpr double kStuckStateChange = 1e-3;     // end-of-cycle state change [-]
    static constexpr double kStuckRate = 0.9;             // per-cycle ratio of |balance|
    static constexpr double kStepRetryRatio = 0.25;        // rejected attempts per accepted step
    static constexpr long kStepRetryMinSteps = 50;        // accepted steps in the cycle
    // Cycles the diagnostics read back: the CSS distance window (6 cycles)
    // and the balance window.
    static constexpr int kHistoryCycles = 6;
  };

  struct RunHealthWarning {
    RunHealthCheck check = RunHealthCheck::MoleFractionSum;
    std::string code;                 // kRunHealthCodeList entry
    double value = 0.0;               // the measure (see the header comment)
    double threshold = 0.0;           // the threshold it passed
    int cycle = -1;                   // 0-based cycle; -1 for the run as a whole
    int step = -1;                    // 0-based step index; -1 when not one step
    std::string step_name;
    std::string component;            // MASS_BALANCE_STUCK: the component
    std::vector<std::string> suspects;  // what to look at first
    std::string hint;                 // what to try
    std::string message;              // one line for the log ("Warning: ...")
  };

  // The measures of one completed cycle (kept for every cycle of the run).
  struct RunHealthCycle {
    int cycle = 0;
    double mole_fraction_sum_error = std::numeric_limits<double>::quiet_NaN();  // max |sum y - 1|
    int mole_fraction_sum_step = -1;
    double state_change = std::numeric_limits<double>::quiet_NaN();  // NaN: not available
    // Per component (RunHealthReport::components order): the cycle mole
    // balance (quadrature, else trapezoid); NaN when not available.
    std::vector<double> mass_balance_error;
    long steps = 0;                   // accepted internal steps
    long rejected = 0;                // error-test + convergence failures + recoveries
  };

  struct RunHealthReport {
    std::vector<std::string> components;
    std::vector<RunHealthCycle> cycles;
    std::vector<RunHealthWarning> warnings;
    bool Healthy() const { return warnings.empty(); }
  };

  // Collects the measures cycle by cycle and turns them into warnings. The
  // solvers own one per run.
  class RunHealthMonitor {
  public:
    // Starts a run of `data` (clears everything).
    void Begin(const UserData& data);
    // After a cycle completes (before the CSS check stops the run), with the
    // integrator counters at the cycle's start and end (their recovered_steps
    // is the run statistics' count, so the warning and the stats agree).
    void OnCycle(const UserData& data, const EquationLog& eq_log, const BoundaryLog& b_log, int cycle,
                 const IntegrationStats& at_start, const IntegrationStats& at_end);
    // At the end of a run that completed (or reached CSS): the run-level
    // checks. Idempotent.
    void Finish(const UserData& data, const EquationLog& eq_log);
    const RunHealthReport& Report() const { return report_; }

  private:
    void Offer(RunHealthWarning w);
    RunHealthReport report_;
    bool finished_ = false;
  };

} // namespace adsim
