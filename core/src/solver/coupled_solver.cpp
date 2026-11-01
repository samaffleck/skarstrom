/*
 * Skarstrom - pressure and temperature swing adsorption simulator
 * Copyright (c) 2024-2026 Sam Affleck and contributors
 *
 * Licensed under the terms in the LICENSE file at the repository root.
 */
#include "../pch.h"
#include "adsim/solver/coupled_solver.h"
#include "adsim/solver/cyclic_steady_state.h"
#include "adsim/solver/log.h"
#include "adsim/solver/ida_settings.h"
#include "adsim/solver/deck_warnings.h"
#include "adsim/solver/operating_limits.h"
#include <chrono>
#include <optional>
#include <sstream>
#include <iomanip>

namespace adsim {

  namespace {
    // Without a tight max-step cap IDA can take steps much longer than the
    // logging resolution. KPIs and the conservation checks integrate the logged
    // boundary flows, so when an accepted step (t_prev, t_curr] is longer than
    // `interval`, log evenly spaced points inside it from IDA's interpolating
    // polynomial (accurate to the integrator's order). log_at(t) records the
    // state already loaded into user_data.state.
    template <class LogAt>
    void LogFillInPoints(SolverData& solver_data, UserData& user_data, double t_prev,
                         double t_curr, double interval, LogAt log_at) {
      if (!(interval > 0.0) || t_curr - t_prev <= interval) return;
      const int n_extra = static_cast<int>(std::ceil((t_curr - t_prev) / interval)) - 1;
      for (int k = 1; k <= n_extra; ++k) {
        const double t = t_prev + (t_curr - t_prev) * k / (n_extra + 1);
        if (solver_data.UnpackInterpolated(user_data, t)) log_at(t);
      }
    }

    // Logging resolution of a cycle step: at least max_log_points_per_step_
    // points (the CSV output resolution), whatever IDA's step size.
    double LogInterval(const Step& step, const IDASettings& settings) {
      return settings.max_log_points_per_step_ > 0
                 ? step.step_time_ / settings.max_log_points_per_step_ : 0.0;
    }

    // " Residual evaluation threw: <what>." when UpdateResiduals caught an
    // exception since the last call (it then returns an IDA error code), else
    // empty. Consumes the stored message.
    std::string ResidualErrorSuffix(UserData& user_data) {
      const std::string what = user_data.TakeResidualError();
      return what.empty() ? std::string{} : " Residual evaluation threw: " + what + ".";
    }
  }

  CoupledSolver::CoupledSolver(UserData& data) : user_data_(data) {}

  void CoupledSolver::LogMessage(LogLevel level, const std::string& msg) {
    if (message_callback_ && level <= log_level_) {
      message_callback_(level, msg);
    }
  }

  /*
    Initialising should only be done once at the start of the simulation
  */
  bool CoupledSolver::Initialise(Step& first_step, std::string& error_message) {
    auto& state = user_data_.state;      
    
    // Initialise member variables 
    user_data_.Initialise();
    equation_log_.Initialise(user_data_);
    boundary_log_.Initialise(user_data_);
    
    // Set current step for initial conditions
    state.SetCurrentStep(&first_step);
    user_data_.step_start_time_ = 0.0;
    first_step.Initialise(state);
    
    int flag = solver_data_.Initialise(user_data_, error_message);
    if (flag != IDA_SUCCESS) {
      return false;
    }

    return true;
  }

  int CoupledSolver::IntegrateStep(Step& step, bool verbose, int cycle_index, int step_index) {
    step.Initialise(user_data_.state);
    user_data_.state.domain_vectors.PopulateEnergySources(step.energy_sources_, user_data_.bed);

    // Get current time from IDA (to continue from where previous step ended)
    double t_start =  solver_data_.GetTime();
    
    // Set step start time for boundary conditions (they expect time relative to step start)
    user_data_.step_start_time_ = t_start;
    
    // Calculate final time for this step (current time + step duration)
    double t_final = t_start + step.step_time_;
    int flag = solver_data_.SetStopTime(t_final);
    solver_data_.OnStepStart(cycle_index);

    double t_curr = t_start;
    
    auto& profile = user_data_.profile_;
    ++profile.cycle_steps;

    // Log initial state
    {
      ProfileTimer log_timer(profile.enabled, profile.log_ns);
      user_data_.state.UpdateIndependentEquations(&user_data_);
      user_data_.state.UpdateDiagnosticEquations(&user_data_);
      equation_log_.LogState(user_data_.state, t_curr, step_index, cycle_index);
      double t_relative = t_curr - user_data_.step_start_time_;
      step.LogOutletConditions(user_data_.state, t_relative);
    }
  
    if (verbose) {
      std::cout << "Solving Step: " << step.step_name_ << " (from t=" << t_start << " to t=" << t_final << ", step time: 0 to " << step.step_time_ << ")" << std::endl;
    }

    const long int max_iters = user_data_.settings.max_iters_per_step;
    long int iter_count = 0;

    // Retry state. When IDA exhausts its own internal step-shrink budget and
    // returns a recoverable failure, rather than aborting the whole run we
    // re-attempt the step a few times (see SolverData::PrepareStepRetry: before
    // the first accepted step this restarts at half the failed step; afterwards
    // IDA continues from its own reduced step with a fresh failure budget).
    // Consecutive failures are counted; a successful step resets it.
    constexpr int kMaxStepRetries = 4;
    const double min_step = user_data_.settings.min_step.Get();
    int retry_count = 0;

    // Profiling: error-test failures in the first kBoundaryWindow accepted steps
    // after the BC switch, and steps throttled by max_step.
    const double max_step = user_data_.settings.max_step.Get();
    long accepted_in_step = 0;
    long etf_at_step_start = profile.enabled ? solver_data_.GetStatistics().num_err_test_fails : 0;
    bool boundary_window_open = profile.enabled;
    auto close_boundary_window = [&]() {
      if (!boundary_window_open) return;
      profile.boundary_err_test_fails +=
          solver_data_.GetStatistics().num_err_test_fails - etf_at_step_start;
      boundary_window_open = false;
    };

    const double log_interval = LogInterval(step, user_data_.settings);

    while (true) {
      const double t_prev = t_curr;
      // A band audit found a coupling outside the probed matrix band: widen it
      // to the full reach and restart from the current point (rare; the probe
      // at the initial state covers every example deck).
      if (solver_data_.BandNeedsWidening()) {
        if (!solver_data_.WidenBand(user_data_)) {
          flag = IDA_LSETUP_FAIL;
          break;
        }
        LogMessage(LogLevel::MEDIUM,
                   "  Step \"" + step.step_name_ + "\": Jacobian coupling outside the probed band at t="
                   + std::to_string(t_curr) + "; widened the band to the full reach.");
      }
      // flag = IDASolve(ida_mem_, t_final, &t_curr, y_, yp_, IDA_ONE_STEP);
      {
        ProfileTimer step_timer(profile.enabled, profile.one_step_ns);
        flag = solver_data_.OneStep(t_final, t_curr);
      }

      if (flag < 0) {
        // Recoverable failures where a retry may be worth attempting (as opposed
        // to fatal input/memory errors, which will not improve with a smaller
        // step and should abort immediately).
        const bool recoverable =
            (flag == IDA_ERR_FAIL      // repeated error-test failures
             || flag == IDA_CONV_FAIL  // repeated nonlinear convergence failures
             || flag == IDA_REP_RES_ERR // repeated recoverable residual errors
             || flag == IDA_LSETUP_FAIL // linear solver setup failure
             || flag == IDA_LSOLVE_FAIL // linear solver solve failure
             || flag == IDA_TOO_MUCH_WORK);

        SolverData::StepRetry retry;
        if (recoverable && retry_count < kMaxStepRetries)
          retry = solver_data_.PrepareStepRetry(flag, min_step);

        if (retry.retry) {
          ++retry_count;
          std::ostringstream ss;
          ss << "  Step \"" << step.step_name_ << "\" IDA failure (code " << flag
             << ") at t=" << t_curr << "; retry " << retry_count << "/" << kMaxStepRetries
             << (retry.consistent_ic ? " restarting from a consistent initial condition"
                 : retry.restarted ? " restarting at step " : " continuing from IDA's step ");
          if (!retry.consistent_ic) ss << retry.h << "s";
          ss << ".";
          LogMessage(LogLevel::MEDIUM, ss.str());
          continue;
        }

        // Either a fatal (non-recoverable) code, the retry budget is exhausted,
        // or a retry could not change anything: fail with a specific message
        // instead of aborting silently.
        if (recoverable) {
          std::ostringstream ss;
          ss << "  Step \"" << step.step_name_ << "\" failed to converge after "
             << retry_count << " retr" << (retry_count == 1 ? "y" : "ies")
             << " (last IDA code " << flag << ", min step " << min_step << "s) at t="
             << t_curr << ".";
          LogMessage(LogLevel::LOW, ss.str());
        }
        break;
      }

      // Accepted step: reset the consecutive-failure counter.
      retry_count = 0;

      if (profile.enabled) {
        if (max_step > 0.0 && std::abs(solver_data_.GetLastStep()) >= max_step * (1.0 - 1e-9))
          ++profile.steps_at_max_step;
        if (++accepted_in_step >= SolverProfile::kBoundaryWindow) close_boundary_window();
      }

      if (verbose) {
        LogTime(t_curr);
      }

      // Log state if logging is enabled
      {
        ProfileTimer log_timer(profile.enabled, profile.log_ns);
        LogFillInPoints(solver_data_, user_data_, t_prev, t_curr, log_interval, [&](double t) {
          equation_log_.LogState(user_data_.state, t, step_index, cycle_index);
          step.LogOutletConditions(user_data_.state, t - user_data_.step_start_time_);
        });
      }
      Log(step, t_curr, step_index, cycle_index);

      // Stream per-step convergence metrics (time step + local error norm) so a
      // GUI monitor can track solver progress and see where it slows down.
      if (metrics_callback_) {
        metrics_callback_(t_curr, solver_data_.GetLastStep(), solver_data_.GetErrorNorm(),
                          cycle_index, step_index);
      }

      // Stop conditions
      if (flag == IDA_TSTOP_RETURN || t_curr >= t_final) {
        break;
      }

      if (stop_requested_.load()) {
        flag = 99;
        break;
      }

      if (max_iters > 0 && ++iter_count >= max_iters) {
        LogMessage(LogLevel::LOW,
                   "  Step \"" + step.step_name_ + "\" exceeded max OneStep iterations ("
                   + std::to_string(max_iters) + "); aborting.");
        flag = IDA_TOO_MUCH_WORK;
        break;
      }

      // Handle rootfinding events
      if (flag == IDA_ROOT_RETURN) {
        continue;
      }
    }

    close_boundary_window();

    {
      ProfileTimer log_timer(profile.enabled, profile.log_ns);
      boundary_log_.LogEndOfStep(step, cycle_index, step_index);
    }

    // Snapshot the just-finished step's outlet aggregates so any downstream
    // BufferTankSource consumer (later this cycle or any future cycle) reads
    // the most recent run's averages.
    step.UpdateCycleSummary(cycle_index);

    return flag;
  }

  void CoupledSolver::Log(Step& step, double t_curr, int step_index, int cycle_index) {
    ProfileTimer log_timer(user_data_.profile_.enabled, user_data_.profile_.log_ns);
    solver_data_.Unpack(user_data_);
    equation_log_.LogState(user_data_.state, t_curr, step_index, cycle_index);
    double t_relative = t_curr - user_data_.step_start_time_;
    step.LogOutletConditions(user_data_.state, t_relative);
  }

  int CoupledSolver::IntegrateCycle(Cycle& cycle, bool verbose) {
    auto wall_start = std::chrono::steady_clock::now();
    user_data_.profile_.Reset();
    css_cycle_index_ = -1;

    user_data_.TakeResidualError(); // drop any message from a previous run
    std::string error_message;
    bool checks_successful = user_data_.PerformChecks(error_message);
    if (!checks_successful) {
      LogMessage(LogLevel::LOW, "Checks failed: " + error_message);
      return -1;
    }

    // Invalidate the buffer-tank cycle summaries (which are also serialised
    // with the deck) so cycle 1 always uses the first-cycle assumption. Must
    // precede Initialise: the IDA consistent-IC solve already evaluates the
    // first step's inlet flows.
    for (auto& step : cycle.step_list_)
      step.ClearCycleSummary();

    Step& first_step = cycle.step_list_[0];
    bool init_successful = Initialise(first_step, error_message);
    if (!init_successful) {
      LogMessage(LogLevel::LOW, "Initialisation failed: " + error_message);
      return -1;
    }

    // Clear stream_data_ for every step so connected boundary streams don't
    // read stale data from a previous simulation's last cycle.
    for (auto& step : cycle.step_list_)
      step.ClearOutletConditions();

    // Bind every consumer's BufferTankSource.summary_ pointer to the matching
    // source endpoint summary. Must happen after PerformChecks (which already
    // calls this) but is repeated here because the solver can be re-entered
    // without a fresh PerformChecks call.
    cycle.ResolveBufferTankSources();

    int sys_size = user_data_.state.GetSystemSize();
    int total_steps = static_cast<int>(cycle.step_list_.size());
    int total_cycles = cycle.number_of_cycles_;

    {
      std::ostringstream ss;
      ss << "Solver initialised (system size: " << sys_size
         << ", steps per cycle: " << total_steps
         << ", max cycles: " << total_cycles << ")";
      LogMessage(LogLevel::LOW, ss.str());
    }
    // Advisory, once per run: a constant pressure boundary far from the
    // pressure its step starts at (see deck_warnings.h).
    const auto& warn = user_data_.settings.warnings_;
    if (warn.pressure_step_change)
      for (const auto& w : BoundaryPressureStepChangeWarnings(cycle, user_data_.ic.pressure.Get()))
        LogMessage(LogLevel::LOW, w);
    if (warn.wall_thickness)
      for (const auto& w : WallDesignWarnings(user_data_.bed))
        LogMessage(LogLevel::LOW, w);
    for (const auto& w : MixingRuleWarnings(user_data_.adsorbent_directory, user_data_.fluid))
      LogMessage(LogLevel::LOW, w);

    stop_requested_ = false;

    int flag = IDA_SUCCESS;
    int cycle_index = 0;
    for (int i = 1; i <= total_cycles; i++) {
      if (stop_requested_) return 99; // Interrupted by user

      cycle_index = i - 1; // Convert to 0-based index

      LogMessage(LogLevel::LOW, "Cycle " + std::to_string(i) + " of " + std::to_string(total_cycles));

      int step_index = 0;
      for (auto& step : cycle.step_list_) {
        if (stop_requested_) return 99; // Interrupted by user
        // Clear logs for each step (streams will initialize on first use if needed)
        step.ClearOutletConditions();

        LogMessage(LogLevel::MEDIUM, "  Step " + std::to_string(step_index + 1)
                   + "/" + std::to_string(total_steps) + ": \"" + step.step_name_ + "\"");

        auto step_wall_start = std::chrono::steady_clock::now();
        auto ida_steps_before = solver_data_.GetStatistics().num_steps;

        flag = IntegrateStep(step, verbose, cycle_index, step_index);

        if (flag == 99) return 99; // Interrupted by user

        if (flag < 0) {
          // Build human-readable failure message
          ChunkStatus temp_status;
          temp_status.Failed("", flag);
          std::string ida_msg = temp_status.GetFlagMessage();
          LogMessage(LogLevel::LOW, "  Step \"" + step.step_name_ + "\" failed: " + ida_msg
                     + " Largest local error: " + solver_data_.DescribeLargestLocalError(user_data_.state) + "."
                     + ResidualErrorSuffix(user_data_));
          break;
        }

        // HIGH: per-step timing and IDA step count
        auto step_wall_end = std::chrono::steady_clock::now();
        double step_wall_secs = std::chrono::duration<double>(step_wall_end - step_wall_start).count();
        auto ida_steps_after = solver_data_.GetStatistics().num_steps;

        {
          std::ostringstream ss;
          ss << "  Step \"" << step.step_name_ << "\" completed ("
             << std::fixed << std::setprecision(3) << step_wall_secs << " s, "
             << (ida_steps_after - ida_steps_before) << " IDA steps)";
          LogMessage(LogLevel::HIGH, ss.str());
        }

        if (progress_callback_) {
          progress_callback_(cycle_index, step_index,
                             total_cycles, total_steps);
        }

        step_index++;
      }

      if (flag < 0) {
        break; // Stop if a step failed
      }

      std::optional<ProfileTimer> css_timer;  // CSS check + acceleration
      css_timer.emplace(user_data_.profile_.enabled, user_data_.profile_.css_ns);

      auto& css = user_data_.settings.css;
      CSSResult css_result = css.CheckCSS(equation_log_, cycle_index);

      if (css_result.converged) {
        std::ostringstream ss;
        ss << "CSS reached at cycle " << i
           << " (" << css_result.method_name
           << ", error: " << std::scientific << std::setprecision(3) << css_result.error
           << ", tolerance: " << css_result.tolerance << ")";
        LogMessage(LogLevel::LOW, ss.str());
        css_cycle_index_ = cycle_index;
        break;
      } else if (css.method_ != CSSMethod::NONE && css_result.tolerance > 0.0) {
        std::ostringstream ss;
        ss << "  CSS check: not converged (" << css_result.method_name
           << ", error: " << std::scientific << std::setprecision(3) << css_result.error
           << ", tolerance: " << css_result.tolerance << ")";
        LogMessage(LogLevel::MEDIUM, ss.str());
      }

      auto& state = user_data_.state;
      auto& acc_css = user_data_.settings.acc_css;
      bool accelerated = acc_css.Accelerate(state, cycle_index);
      if (accelerated) {
        // Set up the first step's boundary conditions so the re-seed is validated
        // and committed under the context the next cycle will actually run in.
        // That includes the step start time: IDA sits at the cycle boundary, which
        // is where the next cycle's first step starts, so time-varying boundary
        // values must be evaluated at t_rel = 0 during the probe and the live
        // IDACalcIC (not relative to the start of the previous cycle's last step).
        Step& first_step = cycle.step_list_[0];
        user_data_.step_start_time_ = solver_data_.GetTime();
        first_step.Initialise(state);
        state.domain_vectors.PopulateEnergySources(first_step.energy_sources_, user_data_.bed);

        // Validate the extrapolated state in a throwaway integrator before
        // touching the live one. Only if it re-initialises and integrates cleanly
        // is the main integrator re-seeded; otherwise the extrapolation is
        // discarded and the run continues, uninterrupted, from the unaccelerated
        // state. This makes acceleration strictly an optimisation: a re-seed that
        // IDA cannot integrate (which does happen at some cycle-step boundaries)
        // costs a little work but never fails the run.
        const double probe_dt = std::min(2.0, 0.1 * first_step.step_time_);
        std::string why;
        if (solver_data_.ValidateAndCommitReseed(user_data_, probe_dt, &why)) {
          LogMessage(LogLevel::MEDIUM, "  CSS acceleration applied");
        } else {
          LogMessage(LogLevel::MEDIUM, "  CSS acceleration skipped (" + why + ")");
        }
      }
    }

    // Wall time
    auto wall_end = std::chrono::steady_clock::now();
    double wall_secs = std::chrono::duration<double>(wall_end - wall_start).count();

    if (flag >= 0) {
      std::ostringstream ss;
      ss << "Solver completed (" << cycle_index + 1 << " cycle"
         << (cycle_index > 0 ? "s" : "") << ", "
         << std::fixed << std::setprecision(2) << wall_secs << " s)";
      LogMessage(LogLevel::LOW, ss.str());
      // Advisory, once per run: the last cycle against the wall design and
      // fluidisation (see operating_limits.h).
      for (const auto& w : CheckOperatingLimits(user_data_, equation_log_).warnings)
        LogMessage(LogLevel::LOW, w.message);
    } else {
      std::ostringstream ss;
      ss << "Solver failed after " << std::fixed << std::setprecision(2) << wall_secs << " s";
      LogMessage(LogLevel::LOW, ss.str());
    }

    return flag;
  }

  CoupledSolverChuncked::CoupledSolverChuncked(UserData& data) : user_data_(data) {}

  void CoupledSolverChuncked::LogMessage(LogLevel level, const std::string& msg) {
    if (message_callback_ && level <= log_level_) {
      message_callback_(level, msg);
    }
  }

  void CoupledSolverChuncked::LogRunSummary(bool success) {
    double wall_secs = std::chrono::duration<double>(
        std::chrono::steady_clock::now() - run_wall_start_).count();

    if (success) {
      std::ostringstream ss;
      ss << "Solver completed (" << std::fixed << std::setprecision(2)
         << wall_secs << " s)";
      LogMessage(LogLevel::LOW, ss.str());
      for (const auto& w : CheckOperatingLimits(user_data_, equation_log_).warnings)
        LogMessage(LogLevel::LOW, w.message);
    } else {
      std::ostringstream ss;
      ss << "Solver failed after " << std::fixed << std::setprecision(2)
         << wall_secs << " s";
      LogMessage(LogLevel::LOW, ss.str());
    }
  }

  void CoupledSolverChuncked::Stop() {
    status_.Interupt();
  }

  void CoupledSolverChuncked::Start() {
    auto& cycle = user_data_.cycle;
    auto& state = user_data_.state;

    if (!status_.IsActive()) { // Solver has not started
      user_data_.TakeResidualError(); // drop any message from a previous run
      std::string error_message;
      bool checks_successful = user_data_.PerformChecks(error_message);
      if (!checks_successful) {
        LogMessage(LogLevel::LOW, "Checks failed: " + error_message);
        status_.Failed("User data checks failed! " + error_message, IDA_ILL_INPUT);
        return;
      }

      // Initialise member variables
      user_data_.Initialise();
      equation_log_.Initialise(user_data_);
      boundary_log_.Initialise(user_data_);

      // Clear stream_data_ for every step so connected boundary streams don't
      // read stale data from a previous simulation's last cycle, and invalidate
      // the buffer-tank cycle summaries (which are also serialised with the
      // deck) so cycle 1 always uses the first-cycle assumption.
      for (auto& step : cycle.step_list_) {
        step.ClearOutletConditions();
        step.ClearCycleSummary();
      }

      // Set current step for initial conditions
      Step& first_step = cycle.step_list_[0];
      state.SetCurrentStep(&first_step);
      user_data_.step_start_time_ = 0.0;
      first_step.Initialise(state);

      int flag = solver_data_.Initialise(user_data_, error_message);
      if (flag != IDA_SUCCESS) {
        LogMessage(LogLevel::LOW, "Initialisation failed: " + error_message);
        status_.Failed("Initialisation failed! " + error_message, flag);
        return;
      }

      // Set activate to true
      status_.Activate();
      run_wall_start_ = std::chrono::steady_clock::now();

      {
        std::ostringstream ss;
        ss << "Solver initialised (system size: " << state.GetSystemSize()
           << ", steps per cycle: " << cycle.step_list_.size()
           << ", max cycles: " << cycle.number_of_cycles_ << ")";
        LogMessage(LogLevel::LOW, ss.str());
      }
      const auto& warn = user_data_.settings.warnings_;
      if (warn.pressure_step_change)
        for (const auto& w : BoundaryPressureStepChangeWarnings(cycle, user_data_.ic.pressure.Get()))
          LogMessage(LogLevel::LOW, w);
      if (warn.wall_thickness)
        for (const auto& w : WallDesignWarnings(user_data_.bed))
          LogMessage(LogLevel::LOW, w);
      for (const auto& w : MixingRuleWarnings(user_data_.adsorbent_directory, user_data_.fluid))
        LogMessage(LogLevel::LOW, w);
    }
  }

  void CoupledSolverChuncked::Integrate(int max_time_ms) {    
    if (!status_.IsActive()) {
      status_.Failed("Solver is not initialised!", IDA_ILL_INPUT);
      return;
    }

    auto& cycle = user_data_.cycle;    
    auto& css = user_data_.settings.css;
    auto& step_list = cycle.step_list_;
    int total_cycles = cycle.number_of_cycles_;    
    int total_steps = int(step_list.size());

    // Process steps until time limit is reached or simulation completes
    while (status_.GetCycleIndex() < total_cycles) {
      
      // Process step
      ProcessStep(max_time_ms);

      // Check if we've completed all steps in current cycle
      if (status_.GetStepIndex() >= total_steps) {

        // Check for cyclic steady state before moving to next cycle
        bool is_at_css = css.IsAtCSS(equation_log_, status_.GetCycleIndex());
        if (is_at_css) {
          LogMessage(LogLevel::LOW, "CSS reached at cycle "
                     + std::to_string(status_.GetCycleIndex() + 1));
          LogRunSummary(true);
          status_.ReachedCSS(status_.GetCycleIndex(), IDA_SUCCESS);
          return;
        }

        auto& acc_css = user_data_.settings.acc_css;
        bool accelerated = acc_css.Accelerate(user_data_.state, status_.GetCycleIndex());
        if (accelerated) {
          // Set up the first step's BC context, then validate the extrapolated
          // state in a throwaway integrator; only commit it to the live solver if
          // it integrates cleanly (see CoupledSolver::IntegrateCycle). This keeps
          // acceleration from ever aborting a run. The step start time is the
          // cycle boundary IDA sits at, so time-varying boundary values are
          // evaluated at t_rel = 0 during the probe and the live IDACalcIC.
          Step& first_step = cycle.step_list_[0];
          user_data_.step_start_time_ = solver_data_.GetTime();
          first_step.Initialise(user_data_.state);
          user_data_.state.domain_vectors.PopulateEnergySources(first_step.energy_sources_, user_data_.bed);
          const double probe_dt = std::min(2.0, 0.1 * first_step.step_time_);
          std::string why;
          if (solver_data_.ValidateAndCommitReseed(user_data_, probe_dt, &why)) {
            LogMessage(LogLevel::MEDIUM, "  CSS acceleration applied");
          } else {
            LogMessage(LogLevel::MEDIUM, "  CSS acceleration skipped (" + why + ")");
          }
        }

        status_.CycleComplete();

        if (status_.GetCycleIndex() >= total_cycles) {
          LogRunSummary(true);
          status_.Completed("Simulation completed successfully.", IDA_SUCCESS);
          return;
        }
      }

      return; // return after each process

    } // end while loop

    // All cycles completed
    LogRunSummary(true);
    status_.Completed("Simulation completed successfully.", IDA_SUCCESS);

  }

  void CoupledSolverChuncked::ProcessStep(int max_time_ms) {
    auto& state = user_data_.state;
    auto& cycle = user_data_.cycle;
    auto& step_list = cycle.step_list_;
    int cycle_index = status_.GetCycleIndex();
    int step_index = status_.GetStepIndex();
    int total_steps = int(step_list.size());

    if (status_.GetStepIndex() >= total_steps) return;
    
    Step& step = step_list[step_index];
    
    // Initialize step if needed
    if (!status_.IsStepInitialised()) {
      if (!InitialiseStep()) {
        return;
      }
    }
    
    // Get target time for this step
    double t_final = status_.GetEndTime();
    double t_start = status_.GetStartTime();
    double t_curr = solver_data_.GetTime();
    
    // Start timer for this chunk
    auto start_time = std::chrono::steady_clock::now();
    auto max_duration = std::chrono::milliseconds(max_time_ms);
    
    // Retry state (mirrors CoupledSolver::IntegrateStep): on a recoverable IDA
    // failure, retry the step a few times (SolverData::PrepareStepRetry) before
    // giving up with a specific error.
    constexpr int kMaxStepRetries = 4;
    const double min_step = user_data_.settings.min_step.Get();
    int retry_count = 0;

    auto complete_step = [&]() {
      // HIGH: per-step wall-clock and IDA step count (step spans many chunks).
      double step_wall_secs = std::chrono::duration<double>(
          std::chrono::steady_clock::now() - step_wall_start_).count();
      long int ida_steps = solver_data_.GetStatistics().num_steps - ida_steps_before_step_;
      std::ostringstream ss;
      ss << "  Step \"" << step.step_name_ << "\" completed ("
         << std::fixed << std::setprecision(3) << step_wall_secs << " s, "
         << ida_steps << " IDA steps)";
      LogMessage(LogLevel::HIGH, ss.str());

      status_.StepComplete();
      boundary_log_.LogEndOfStep(step, cycle_index, step_index);

      // Snapshot the just-finished step's outlet aggregates for any
      // downstream BufferTankSource consumer (mirrors
      // CoupledSolver::IntegrateStep). cycle_index / step_index were read
      // before StepComplete(), so they still name the step that finished.
      step.UpdateCycleSummary(cycle_index);
    };

    // Defensive: a zero-length step (rejected by Cycle::PerformChecks) would never
    // enter the loop below, so it would never complete and every chunk would
    // return IN_PROGRESS forever. Complete it immediately instead.
    if (t_curr >= t_final) {
      complete_step();
      return;
    }

    const double log_interval = LogInterval(step, user_data_.settings);

    // Integrate step with time limit
    while (t_curr < t_final) {
      const double t_prev = t_curr;

      // Exit condition
      auto elapsed = std::chrono::steady_clock::now() - start_time;
      if (elapsed >= max_duration) {
        return;
      }

      // Take one step
      int flag = solver_data_.OneStep(t_final, t_curr);

      if (flag < 0) {
        const bool recoverable =
            (flag == IDA_ERR_FAIL || flag == IDA_CONV_FAIL
             || flag == IDA_REP_RES_ERR || flag == IDA_LSETUP_FAIL
             || flag == IDA_LSOLVE_FAIL || flag == IDA_TOO_MUCH_WORK);

        if (recoverable && retry_count < kMaxStepRetries
            && solver_data_.PrepareStepRetry(flag, min_step).retry) {
          ++retry_count;
          continue;
        }

        ChunkStatus temp_status;
        temp_status.Failed("", flag);
        const std::string residual_error = ResidualErrorSuffix(user_data_);
        LogMessage(LogLevel::LOW, "  Step \"" + step.step_name_ + "\" failed: "
                   + temp_status.GetFlagMessage() + " Largest local error: "
                   + solver_data_.DescribeLargestLocalError(user_data_.state) + "."
                   + residual_error);
        status_.Failed(
            (recoverable
                ? ("IDASolve failed to converge after " + std::to_string(retry_count)
                   + (retry_count == 1 ? " retry" : " retries")
                   + " (last code " + std::to_string(flag) + ")")
                : ("IDASolve failed with code: " + std::to_string(flag)))
            + residual_error,
            flag);
        LogRunSummary(false);
        return;
      }

      // Accepted step: reset the consecutive-failure counter.
      retry_count = 0;

      LogFillInPoints(solver_data_, user_data_, t_prev, t_curr, log_interval, [&](double t) {
        equation_log_.LogState(user_data_.state, t, status_.GetStepIndex(), cycle_index);
        step.LogOutletConditions(user_data_.state, t - user_data_.step_start_time_);
      });
      Log(step, t_curr, status_.GetStepIndex(), cycle_index);

      // Check if step completed
      if (flag == IDA_TSTOP_RETURN || t_curr >= t_final) {
        complete_step();
        break;
      }
      
      // Handle rootfinding events
      if (flag == IDA_ROOT_RETURN) {
        continue;
      }
    }
    
    return;
  }

  bool CoupledSolverChuncked::InitialiseStep() {
    auto& state = user_data_.state;
    auto& cycle = user_data_.cycle;
    auto& step_list = cycle.step_list_;
    int cycle_index = status_.GetCycleIndex();
    int step_index = status_.GetStepIndex();
    int total_steps = int(step_list.size());
    Step& step = step_list[step_index];

    // A new cycle begins when we initialise its first step.
    if (step_index == 0) {
      LogMessage(LogLevel::LOW, "Cycle " + std::to_string(cycle_index + 1)
                 + " of " + std::to_string(user_data_.cycle.number_of_cycles_));
    }
    LogMessage(LogLevel::MEDIUM, "  Step " + std::to_string(step_index + 1)
               + "/" + std::to_string(total_steps) + ": \"" + step.step_name_ + "\"");

    // Capture markers for the HIGH per-step timing line (read at step completion).
    step_wall_start_ = std::chrono::steady_clock::now();
    ida_steps_before_step_ = solver_data_.GetStatistics().num_steps;

    step.Initialise(state);
    state.domain_vectors.PopulateEnergySources(step.energy_sources_, user_data_.bed);

    // Get current time from IDA
    double t_start = solver_data_.GetTime();

    // Set for boundary conditions to calculate offset
    user_data_.step_start_time_ = t_start;

    // Set the final integration time
    double t_final = t_start + step.step_time_;
    int flag = solver_data_.SetStopTime(t_final);

    if (flag != IDA_SUCCESS) {
      LogMessage(LogLevel::LOW, "  Step \"" + step.step_name_
                 + "\" failed: IDASetStopTime code " + std::to_string(flag));
      status_.Failed("IDASetStopTime failed with code: " + std::to_string(flag), flag);
      LogRunSummary(false);
      return false;
    }
    solver_data_.OnStepStart(cycle_index);

    // Clear logs for this step ??
    step.ClearOutletConditions();

    // Log initial state
    Log(step, t_start, step_index, cycle_index);
    
    status_.StepInitialised(t_start, t_final);
    return true;
  }

  double CoupledSolverChuncked::GetProgressPercent() const {
    SolverStatus solver_status = status_.GetStatus();
    if (solver_status == SolverStatus::COMPLETED || solver_status == SolverStatus::AT_CSS) {
      return 100.0;
    } 
    
    if (!status_.IsActive()) {
      return 0.0;
    }
    
    double current_time = solver_data_.GetTime();
    
    auto& cycle = user_data_.cycle;
    double total_time = cycle.GetTotalSimulationTime();

    if (total_time <= 0.0) {
      return 0.0;
    }

    double progress_percent = (current_time / total_time) * 100.0;    
    return progress_percent;
  }

  void CoupledSolverChuncked::Log(Step& step, double t_curr, int step_index, int cycle_index) {
    solver_data_.Unpack(user_data_);
    equation_log_.LogState(user_data_.state, t_curr, step_index, cycle_index);
    double t_relative = t_curr - user_data_.step_start_time_;
    step.LogOutletConditions(user_data_.state, t_relative);
  }

  SolverData::SolverData() = default;

  SolverData::~SolverData() {
    Free();
    // The SUNContext is released afterwards by context_'s destructor, once
    // jacobian_ has also freed its vectors.
  }

  void SolverData::Free() {
    if (ida_mem_) {
      IDAFree(&ida_mem_);
      ida_mem_ = nullptr;
    }
    if (LS_) {
      SUNLinSolFree(LS_);
      LS_ = nullptr;
    }
    if (A_) {
      SUNMatDestroy(A_);
      A_ = nullptr;
    }
    if (y_) {
      N_VDestroy_Serial(y_);
      y_ = nullptr;
    }
    if (yp_) {
      N_VDestroy_Serial(yp_);
      yp_ = nullptr;
    }
    if (id_) {
      N_VDestroy_Serial(id_);
      id_ = nullptr;
    }
    if (abstol_vec_) {
      N_VDestroy_Serial(abstol_vec_);
      abstol_vec_ = nullptr;
    }
    if (ewt_) {
      N_VDestroy_Serial(ewt_);
      ewt_ = nullptr;
    }
    if (ele_) {
      N_VDestroy_Serial(ele_);
      ele_ = nullptr;
    }
    for (N_Vector* v : { &y_interp_, &yp_interp_ }) {
      if (*v) {
        N_VDestroy_Serial(*v);
        *v = nullptr;
      }
    }
  }

  int SolverData::Initialise(UserData& user_data, std::string& error_message) {
    // Release the previous run's integrator: IntegrateCycle re-initialises the
    // same SolverData on every call (e.g. once per optimiser evaluation).
    Free();

    auto& state = user_data.state;
    auto& settings = user_data.settings;

    int sys_size = state.GetSystemSize(); 
    if (sys_size <= 0) {
      error_message += ("Invalid system size: " + std::to_string(sys_size));
      // `false` here would be 0 == IDA_SUCCESS, so callers testing
      // `flag != IDA_SUCCESS` would proceed with nothing allocated.
      return IDA_ILL_INPUT;
    }

    y_ = N_VNew_Serial(sys_size, sunctx_);
    yp_ = N_VNew_Serial(sys_size, sunctx_);
    id_ = N_VNew_Serial(sys_size, sunctx_);
    abstol_vec_ = N_VNew_Serial(sys_size, sunctx_);
    ewt_ = N_VNew_Serial(sys_size, sunctx_); // scratch for convergence metrics
    ele_ = N_VNew_Serial(sys_size, sunctx_);
    y_interp_ = N_VNew_Serial(sys_size, sunctx_);
    yp_interp_ = N_VNew_Serial(sys_size, sunctx_);
    
    // Pack initial conditions
    sunrealtype* y_data = N_VGetArrayPointer(y_);
    sunrealtype* yp_data = N_VGetArrayPointer(yp_);
    sunrealtype* id_data = N_VGetArrayPointer(id_);
    sunrealtype* abstol_data = N_VGetArrayPointer(abstol_vec_);
    
    state.PackBanded(y_data, yp_data);
    state.SetIDBanded(id_data);
    state.SetToleranceVectorBanded(abstol_data);
    
    // Create and configure IDA
    ida_mem_ = IDACreate(sunctx_);
    if (!ida_mem_) {
      error_message += ("Failed to create IDA memory");
      return IDA_MEM_NULL;
    }

    int retval = IDAInit(ida_mem_, UpdateResiduals, 0.0, y_, yp_);
    if (retval != IDA_SUCCESS) {
      error_message += ("IDAInit failed with code: " + std::to_string(retval));
      return retval;
    }

    retval = IDASetUserData(ida_mem_, &user_data);
    if (retval != IDA_SUCCESS) {
      error_message += ("IDASetUserData failed with code: " + std::to_string(retval));
      return retval;
    }

    retval = IDASVtolerances(ida_mem_, settings.reltol, abstol_vec_);
    if (retval != IDA_SUCCESS) {
      error_message += ("IDASVtolerances failed with code: " + std::to_string(retval));
      return retval;
    }

    retval = IDASetId(ida_mem_, id_);
    if (retval != IDA_SUCCESS) {
      error_message += ("IDASetId failed with code: " + std::to_string(retval));
      return retval;
    }

    // The banded layout is cell-major. Cross-variable couplings reach one cell
    // (|i - j| < 2*neq); the WENO / van Leer face stencils reach two cells but
    // only within the same mole fraction, i.e. exactly 2*neq. So 2*neq covers
    // every scheme (pinned by JacobianBand.* in test_jacobian_band.cpp) and is
    // the reach the band audits scan. The matrix itself is sized from the
    // reach a DQ probe measures at the initial state (about 12 / 16 instead of
    // 20 / 20 at neq = 10); a coupling later found outside it widens the band
    // to the full reach (WidenBand, called from the one-step loop).
    size_t num_equations = state.GetEquations().size();
    const sunindextype reach = static_cast<sunindextype>(num_equations * 2);
    sunindextype mu = reach; // Upper bandwidth
    sunindextype ml = reach; // Lower bandwidth
    if (settings.band_from_probe_) {
      sunindextype mu_p = 0, ml_p = 0;
      if (BandedJacobian::ProbeBandwidth(state, &user_data, sunctx_, y_, yp_, abstol_vec_, reach, mu_p, ml_p)
          && mu_p > 0 && ml_p > 0) {
        mu = SUNMIN(reach, mu_p);
        ml = SUNMIN(reach, ml_p);
        if (settings.band_probe_shrink_for_testing_ > 0) {
          mu = SUNMAX(1, mu - settings.band_probe_shrink_for_testing_);
          ml = SUNMAX(1, ml - settings.band_probe_shrink_for_testing_);
        }
      }
    }
    band_mu_ = mu;
    band_ml_ = ml;
    band_reach_ = reach;
    band_widened_ = false;
    
    A_ = SUNBandMatrix(sys_size, mu, ml, sunctx_);
    if (!A_) {
      error_message += ("Failed to create banded matrix");
      return IDA_ILL_INPUT;
    }
    
    LS_ = SUNLinSol_Band(y_, A_, sunctx_);
    if (!LS_) {
      error_message += ("Failed to create banded linear solver");
      return IDA_ILL_INPUT;
    }

    retval = IDASetLinearSolver(ida_mem_, LS_, A_);
    if (retval != IDA_SUCCESS) {
      error_message += ("IDASetLinearSolver failed with code: " + std::to_string(retval));
      return retval;
    }

    custom_jacobian_ = settings.jacobian_mode_ != JacobianMode::IdaInternal;
    user_data.jacobian_ = nullptr;
    if (custom_jacobian_) {
      jacobian_.Initialise(ida_mem_, state, sunctx_, mu, ml, settings.jacobian_mode_,
                           settings.record_jacobian_pattern_, reach, reach);
      user_data.jacobian_ = &jacobian_;
      retval = IDASetJacFn(ida_mem_, BandedJacobianFn);
      if (retval != IDA_SUCCESS) {
        error_message += ("IDASetJacFn failed with code: " + std::to_string(retval));
        return retval;
      }
    }

    // Apply IDA settings
    retval = ApplyIDASettings(ida_mem_, settings);
    if (retval != IDA_SUCCESS) {
      error_message += ("Applying the solver settings failed with code: " + std::to_string(retval)
                        + ". Please check the solver settings (e.g. the maximum BDF order must be 1-5).");
      return retval;
    }

    // Compute consistent initial conditions
    retval = IDACalcIC(ida_mem_, IDA_YA_YDP_INIT, 1e-9);
    if (retval != IDA_SUCCESS) {
      ChunkStatus temp;
      temp.Failed("", retval);
      error_message += "Initial condition calculation failed. " + temp.GetFlagMessage()
                    + "Please check that your initial conditions are consistent with "
                    + "the boundary conditions of the first step."
                    + ResidualErrorSuffix(user_data);
      return retval;
    }

    // IDACalcIC leaves the corrected values in IDA's history array, not in
    // y_ / yp_ (which still hold the initial guess). Fetch them (only possible
    // before the first IDASolve) so y_ / yp_ and the state logged at t = 0 are
    // the consistent initial condition IDA actually integrates from.
    retval = IDAGetConsistentIC(ida_mem_, y_, yp_);
    if (retval != IDA_SUCCESS) {
      error_message += ("IDAGetConsistentIC failed with code: " + std::to_string(retval));
      return retval;
    }
    state.UnpackBanded(y_data, yp_data);

    return IDA_SUCCESS;
  }

  double SolverData::GetTime() const {
    double time = 0.0;
    int flag = IDAGetCurrentTime(ida_mem_, &time);
    return time;
  }

  int SolverData::SetStopTime(double t_final) {
    int flag = IDASetStopTime(ida_mem_, t_final);
    return flag;  
  }

  int SolverData::OneStep(double t_final, double& t_curr) {
    int flag = IDASolve(ida_mem_, t_final, &t_curr, y_, yp_, IDA_ONE_STEP);
    return flag;
  }

  int SolverData::SetInitStep(double h) {
    if (!ida_mem_ || h <= 0.0) return IDA_ILL_INPUT;
    return IDASetInitStep(ida_mem_, h);
  }

  void SolverData::ResetInitStep(const IDASettings& settings) {
    if (!ida_mem_) return;
    const double h0 = settings.init_step.Get();
    IDASetInitStep(ida_mem_, h0 > 0.0 ? h0 : 0.0);
  }

  double SolverData::GetActualInitStep() const {
    double h0 = 0.0;
    if (ida_mem_) IDAGetActualInitStep(ida_mem_, &h0);
    return h0;
  }

  SolverData::StepRetry SolverData::PrepareStepRetry(int flag, double min_step) {
    StepRetry r;
    if (!ida_mem_) return r;

    // The step IDA attempted last (and failed with): IDA has already shrunk it
    // through its own error-test / convergence-failure handling.
    double h_cur = 0.0;
    IDAGetCurrentStep(ida_mem_, &h_cur);
    h_cur = std::abs(h_cur);
    const double floor = (min_step > 0.0 ? min_step : 0.0);
    const bool at_floor = !(h_cur > floor * (1.0 + 1e-6));

    if (GetNumSteps() == 0) {
      // Before the first accepted step IDA restarts from hin on every IDASolve
      // call, and re-entering without a re-init would also re-scale its
      // already-scaled phi[1] (the predicted derivative). Restart cleanly at the
      // current point (y_/yp_ hold the initial values IDA returned) with hin
      // set to half the failed step.
      if (at_floor) return r;
      r.h = std::max(0.5 * h_cur, floor);
      IDASetInitStep(ida_mem_, r.h);
      if (IDAReInit(ida_mem_, GetTime(), y_, yp_) != IDA_SUCCESS) return StepRetry{};
      r.retry = true;
      r.restarted = true;
      return r;
    }

    // After an accepted step IDA ignores hin, but resets its per-step failure
    // counters on every IDASolve call: re-entering continues from its reduced
    // step with a fresh error-test / convergence-failure budget. That cannot
    // help once a step-size failure has driven the step to the min_step floor.
    const bool step_size_failure =
        (flag == IDA_ERR_FAIL || flag == IDA_CONV_FAIL || flag == IDA_REP_RES_ERR);
    if (step_size_failure && at_floor) {
      // A step-size collapse with history is usually a discontinuity that IDA's
      // multistep history cannot absorb, typically a boundary-condition switch
      // at a cycle-step boundary right after a short step (see
      // RestartFromConsistentIC). Restart once from the last accepted point;
      // if it collapses again before an accepted step, the nst == 0 branch
      // above refuses the retry.
      if (!RestartFromConsistentIC()) return r;
      r.retry = true;
      r.restarted = true;
      r.consistent_ic = true;
      r.h = GetInitStepSetting();
      return r;
    }
    r.retry = true;
    r.h = h_cur;
    return r;
  }

  double SolverData::GetInitStepSetting() const {
    void* ud = nullptr;
    if (!ida_mem_ || IDAGetUserData(ida_mem_, &ud) != IDA_SUCCESS || !ud) return 0.0;
    const double h0 = static_cast<UserData*>(ud)->settings.init_step.Get();
    return h0 > 0.0 ? h0 : 0.0;
  }

  bool SolverData::RestartFromConsistentIC() {
    // Steps run back to back without re-initialising IDA, so a cycle-step
    // boundary that switches a boundary condition (e.g. a flow inlet to a
    // pressure outlet) makes the algebraic ghost-cell values jump. IDA can
    // only absorb that jump with a step whose local error estimate, which
    // scales like h / h_prev times the jump, passes the error test. The first
    // step across therefore has to be a tiny fraction of the last step before
    // the boundary (~1e-8 for a 2 bar ghost-pressure jump at reltol 1e-8), and
    // the step after it smaller still, because the history then reads the
    // jump as an algebraic derivative of order jump / h. When
    // the last step before the boundary is a short remainder clipped by tstop,
    // the required step falls below min_step and the run fails at |h| = hmin.
    //
    // Restart from the last accepted point (which already satisfies the new
    // boundary conditions or lies just before the jump) with fresh history:
    // IDACalcIC makes the algebraic variables and differential derivatives
    // consistent. Algebraic derivatives are not part of the residual, so
    // IDACalcIC leaves them as given; zero them, because the values left from
    // the failed history (~1e17 Pa/s for a ghost pressure) would otherwise
    // put that jump straight back into the first predictor.
    void* ud = nullptr;
    if (!ida_mem_ || IDAGetUserData(ida_mem_, &ud) != IDA_SUCCESS || !ud) return false;
    UserData& user_data = *static_cast<UserData*>(ud);

    const double t = GetTime();
    const sunrealtype* id = N_VGetArrayPointer(id_);
    sunrealtype* yp = N_VGetArrayPointer(yp_);
    const sunindextype n = N_VGetLength(yp_);
    for (sunindextype i = 0; i < n; ++i) {
      if (id[i] == 0.0) yp[i] = 0.0;
    }

    ResetInitStep(user_data.settings);
    int flag = IDAReInit(ida_mem_, t, y_, yp_);
    if (flag == IDA_SUCCESS) flag = IDACalcIC(ida_mem_, IDA_YA_YDP_INIT, t + 1e-8);
    if (flag == IDA_SUCCESS) flag = IDAGetConsistentIC(ida_mem_, y_, yp_);
    if (flag != IDA_SUCCESS) {
      // The run fails with the original IDA error; a residual exception thrown
      // during this attempt is not its cause.
      user_data.TakeResidualError();
      return false;
    }
    user_data.state.UnpackBanded(N_VGetArrayPointer(y_), N_VGetArrayPointer(yp_));
    return true;
  }

  bool SolverData::BandNeedsWidening() const {
    return custom_jacobian_ && A_ && jacobian_.NumOutOfBand() > 0
        && (band_mu_ < band_reach_ || band_ml_ < band_reach_);
  }

  bool SolverData::WidenBand(UserData& user_data) {
    if (!ida_mem_ || !A_ || !LS_) return false;
    const sunindextype n = N_VGetLength(y_);
    SUNMatrix A = SUNBandMatrix(n, band_reach_, band_reach_, sunctx_);
    SUNLinearSolver LS = A ? SUNLinSol_Band(y_, A, sunctx_) : nullptr;
    if (!A || !LS) {
      if (LS) SUNLinSolFree(LS);
      if (A) SUNMatDestroy(A);
      return false;
    }
    // IDASetLinearSolver releases IDA's interface to the old solver; the
    // objects themselves are ours to free. The next IDASolve after the
    // IDAReInit below runs IDA's initial setup again, which initialises the
    // new linear solver and forces a fresh Jacobian.
    if (IDASetLinearSolver(ida_mem_, LS, A) != IDA_SUCCESS) {
      SUNLinSolFree(LS);
      SUNMatDestroy(A);
      return false;
    }
    SUNLinSolFree(LS_);
    SUNMatDestroy(A_);
    LS_ = LS;
    A_ = A;
    band_mu_ = band_reach_;
    band_ml_ = band_reach_;
    band_widened_ = true;
    jacobian_.ResizeBand(band_reach_, band_reach_);
    if (IDASetJacFn(ida_mem_, BandedJacobianFn) != IDA_SUCCESS) return false;
    (void)user_data;
    return RestartFromConsistentIC();
  }

  int SolverData::Reinitialise(UserData& user_data) {
    double t_curr = GetTime();

    // Save current solver state in case reinitialisation fails
    N_Vector y_backup = N_VClone(y_);
    N_Vector yp_backup = N_VClone(yp_);
    N_VScale(1.0, y_, y_backup);
    N_VScale(1.0, yp_, yp_backup);

    // Pack accelerated state
    sunrealtype* y_data = N_VGetArrayPointer(y_);
    sunrealtype* yp_data = N_VGetArrayPointer(yp_);
    user_data.state.PackBanded(y_data, yp_data);

    ResetInitStep(user_data.settings);
    int flag = IDAReInit(ida_mem_, t_curr, y_, yp_);
    if (flag == IDA_SUCCESS) {
      flag = IDACalcIC(ida_mem_, IDA_YA_YDP_INIT, t_curr + 1e-8);
    }

    if (flag != IDA_SUCCESS) {
      // The run carries on from the restored state; don't let a residual
      // exception from the abandoned re-seed be blamed for a later failure.
      user_data.TakeResidualError();
      // Restore original solver state
      N_VScale(1.0, y_backup, y_);
      N_VScale(1.0, yp_backup, yp_);
      IDAReInit(ida_mem_, t_curr, y_, yp_);
      // Restore state vectors
      user_data.state.UnpackBanded(N_VGetArrayPointer(y_), N_VGetArrayPointer(yp_));
    }

    N_VDestroy(y_backup);
    N_VDestroy(yp_backup);
    return IDA_SUCCESS;
  }

  bool SolverData::ValidateAndCommitReseed(UserData& user_data, double probe_dt,
                                           std::string* reason) {
    // Precondition: user_data.state holds the extrapolated (accelerated) state and
    // the first step's BC context is set; the live integrator (ida_mem_, y_, yp_)
    // is still at the cycle boundary and has NOT been re-initialised.
    auto& state = user_data.state;
    auto& settings = user_data.settings;
    const double t0 = GetTime();
    const int sys_size = state.GetSystemSize();
    const sunindextype bw = static_cast<sunindextype>(state.GetEquations().size() * 2);

    // Buffer the extrapolated state so it can be re-applied after validation
    // (residual evaluations during the scratch solve overwrite user_data.state).
    N_Vector y_acc = N_VClone(y_);
    N_Vector yp_acc = N_VClone(yp_);
    state.PackBanded(N_VGetArrayPointer(y_acc), N_VGetArrayPointer(yp_acc));

    // Working copies for the throwaway integrator.
    N_Vector ys = N_VClone(y_acc);
    N_Vector yps = N_VClone(yp_acc);
    N_VScale(1.0, y_acc, ys);
    N_VScale(1.0, yp_acc, yps);

    void* scratch = IDACreate(sunctx_);
    SUNMatrix As = nullptr;
    SUNLinearSolver LSs = nullptr;
    bool ok = (scratch != nullptr);
    if (ok) ok = (IDAInit(scratch, UpdateResiduals, t0, ys, yps) == IDA_SUCCESS);
    if (ok) {
      IDASetUserData(scratch, &user_data);
      IDASVtolerances(scratch, settings.reltol, abstol_vec_);
      IDASetId(scratch, id_);
      As = SUNBandMatrix(sys_size, bw, bw, sunctx_);
      LSs = As ? SUNLinSol_Band(ys, As, sunctx_) : nullptr;
      ok = (As && LSs);
      if (ok) {
        IDASetLinearSolver(scratch, LSs, As);
        ok = (ApplyIDASettings(scratch, settings) == IDA_SUCCESS)
          && (IDACalcIC(scratch, IDA_YA_YDP_INIT, t0 + 1e-8) == IDA_SUCCESS);
      }
    }
    if (ok) {
      // March a short interval: a non-integrable re-seed trips the error test
      // repeatedly down to h = hmin within the first few internal steps.
      const double target = t0 + probe_dt;
      IDASetStopTime(scratch, target);
      double tc = t0;
      for (int i = 0; i < 100; ++i) {
        int f = IDASolve(scratch, target, &tc, ys, yps, IDA_ONE_STEP);
        if (f < 0) { ok = false; break; }
        if (f == IDA_TSTOP_RETURN || tc >= target) break;
      }
    }

    if (scratch) IDAFree(&scratch);
    if (LSs) SUNLinSolFree(LSs);
    if (As) SUNMatDestroy(As);
    N_VDestroy(ys);
    N_VDestroy(yps);

    std::string why = ok ? "" : "re-seed not integrable";
    if (ok) {
      // Commit: re-seed the live integrator with the (validated) extrapolated
      // state. The same IDACalcIC succeeded on the scratch copy, so it should
      // succeed here too, but check: IDAReInit has already discarded the live
      // history, so on failure fall back to the pre-reseed state below.
      N_Vector y_pre = N_VClone(y_);
      N_Vector yp_pre = N_VClone(yp_);
      N_VScale(1.0, y_, y_pre);
      N_VScale(1.0, yp_, yp_pre);

      N_VScale(1.0, y_acc, y_);
      N_VScale(1.0, yp_acc, yp_);
      ResetInitStep(settings);
      int flag = IDAReInit(ida_mem_, t0, y_, yp_);
      if (flag != IDA_SUCCESS) {
        why = "live IDAReInit failed with code " + std::to_string(flag);
      } else {
        flag = IDACalcIC(ida_mem_, IDA_YA_YDP_INIT, t0 + 1e-8);
        if (fail_next_reseed_commit_for_testing_) flag = IDA_NO_RECOVERY;
        if (flag != IDA_SUCCESS) {
          why = "live IDACalcIC failed with code " + std::to_string(flag);
        } else {
          // The corrected IC lives in IDA's history array, not y_ / yp_: fetch it
          // so the state (and the first logged row of the accelerated cycle) is
          // the one IDA integrates from, not the uncorrected extrapolation.
          flag = IDAGetConsistentIC(ida_mem_, y_, yp_);
          if (flag != IDA_SUCCESS) {
            why = "IDAGetConsistentIC failed with code " + std::to_string(flag);
          }
        }
      }
      fail_next_reseed_commit_for_testing_ = false;

      if (flag == IDA_SUCCESS) {
        state.UnpackBanded(N_VGetArrayPointer(y_), N_VGetArrayPointer(yp_));
      } else {
        // The live integrator is half re-initialised. Restart it from the
        // pre-reseed state (its multistep history is gone, so it restarts at
        // order 1) and continue the un-accelerated trajectory.
        ok = false;
        N_VScale(1.0, y_pre, y_);
        N_VScale(1.0, yp_pre, yp_);
        ResetInitStep(settings);
        flag = IDAReInit(ida_mem_, t0, y_, yp_);
        if (flag != IDA_SUCCESS) {
          why += "; restoring the pre-reseed state failed with code " + std::to_string(flag);
        }
      }
      N_VDestroy(y_pre);
      N_VDestroy(yp_pre);
    }
    if (!ok) {
      // Reject: the live integrator holds the good pre-acceleration state (with
      // its continuous history, unless a failed commit forced a restart from it).
      // Restore user_data.state from it and carry on as if no acceleration had
      // been attempted.
      state.UnpackBanded(N_VGetArrayPointer(y_), N_VGetArrayPointer(yp_));
      user_data.TakeResidualError(); // a scratch-solve residual exception is not the run's failure
    }
    if (reason) *reason = why;

    N_VDestroy(y_acc);
    N_VDestroy(yp_acc);
    return ok;
  }


  void SolverData::Unpack(UserData& user_data) {
    sunrealtype* y_data = N_VGetArrayPointer(y_);
    sunrealtype* yp_data = N_VGetArrayPointer(yp_);
    user_data.state.UnpackBanded(y_data, yp_data);
    user_data.state.UpdateIndependentEquations(&user_data);
    user_data.state.UpdateDiagnosticEquations(&user_data);
  }

  bool SolverData::UnpackInterpolated(UserData& user_data, double t) {
    if (!ida_mem_ || !y_interp_ || !yp_interp_) return false;
    if (IDAGetDky(ida_mem_, t, 0, y_interp_) != IDA_SUCCESS) return false;
    if (IDAGetDky(ida_mem_, t, 1, yp_interp_) != IDA_SUCCESS) return false;
    user_data.state.UnpackBanded(N_VGetArrayPointer(y_interp_), N_VGetArrayPointer(yp_interp_));
    user_data.state.UpdateIndependentEquations(&user_data);
    user_data.state.UpdateDiagnosticEquations(&user_data);
    return true;
  }

  IntegrationStats SolverData::GetStatistics() const {
    IntegrationStats stats;
    if (ida_mem_) {
      IDAGetNumSteps(ida_mem_, &stats.num_steps);
      IDAGetNumResEvals(ida_mem_, &stats.num_res_evals);
      IDAGetNumNonlinSolvIters(ida_mem_, &stats.num_nonlin_iters);
      IDAGetNumErrTestFails(ida_mem_, &stats.num_err_test_fails);
      IDAGetNumNonlinSolvConvFails(ida_mem_, &stats.num_conv_fails);
      IDAGetNumLinResEvals(ida_mem_, &stats.num_lin_res_evals);
      IDAGetNumJacEvals(ida_mem_, &stats.num_jac_evals);
      IDAGetNumLinSolvSetups(ida_mem_, &stats.num_lin_setups);
      if (custom_jacobian_) stats.num_lin_res_evals += jacobian_.NumResidualEvals();
    }
    return stats;
  }

  double SolverData::GetLastStep() const {
    double h = 0.0;
    if (ida_mem_) IDAGetLastStep(ida_mem_, &h);
    return h;
  }

  // WRMS norm of IDA's estimated local error, weighted by the error weights
  // ewt = 1/(reltol*|y| + abstol). This is the quantity IDA drives below 1 in
  // its error test, so it is the natural per-step "error" for a convergence
  // monitor: when it approaches 1, the step size is throttled back.
  double SolverData::GetErrorNorm() const {
    if (!ida_mem_ || !ewt_ || !ele_) return 0.0;
    if (IDAGetErrWeights(ida_mem_, ewt_) != IDA_SUCCESS) return 0.0;
    if (IDAGetEstLocalErrors(ida_mem_, ele_) != IDA_SUCCESS) return 0.0;
    return N_VWrmsNorm(ele_, ewt_);
  }

  std::string SolverData::DescribeLargestLocalError(const Equations& state) const {
    if (!ida_mem_ || !ewt_ || !ele_) return "?";
    if (IDAGetErrWeights(ida_mem_, ewt_) != IDA_SUCCESS) return "?";
    if (IDAGetEstLocalErrors(ida_mem_, ele_) != IDA_SUCCESS) return "?";
    const sunrealtype* e = N_VGetArrayPointer(ele_);
    const sunrealtype* w = N_VGetArrayPointer(ewt_);
    const sunindextype n = N_VGetLength(ele_);
    sunindextype worst = 0;
    double worst_value = -1.0;
    for (sunindextype i = 0; i < n; ++i) {
      const double v = std::abs(e[i] * w[i]);
      if (v > worst_value) { worst_value = v; worst = i; }
    }
    return state.DescribeBandedIndex(static_cast<size_t>(worst));
  }

  long int SolverData::GetNumSteps() const {
    long int n = 0;
    if (ida_mem_) IDAGetNumSteps(ida_mem_, &n);
    return n;
  }

  void CoupledSolver::PrintStatistics() const {
    auto stats = solver_data_.GetStatistics();
    std::cout << "Statistics:" << std::endl;
    std::cout << "  Number of steps: " << stats.num_steps << std::endl;
    std::cout << "  Number of residual evals: " << stats.num_res_evals << std::endl;
    std::cout << "  Number of nonlinear iterations: " << stats.num_nonlin_iters << std::endl;
    std::cout << "  Number of error test failures: " << stats.num_err_test_fails << std::endl;
    std::cout << "  Number of convergence failures: " << stats.num_conv_fails << std::endl;
    std::cout << "  Number of Jacobian evals: " << stats.num_jac_evals << std::endl;
    std::cout << "  Number of linear solver setups: " << stats.num_lin_setups << std::endl;
    std::cout << "  Number of Jacobian (DQ) residual evals: " << stats.num_lin_res_evals << std::endl;
  }

} // end adsim namespace
