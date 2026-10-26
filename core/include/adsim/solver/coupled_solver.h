/*
 * Skarstrom - pressure and temperature swing adsorption simulator
 * Copyright (c) 2024-2026 Sam Affleck and contributors
 *
 * Licensed under the terms in the LICENSE file at the repository root.
 */
#pragma once 

// STL includes
#include <atomic>
#include <functional>
#include <string>
#include <chrono>

// adsim includes
#include "adsim/objects/user_data.h"
#include "adsim/solver/log.h"

// SUNDIALS includes
#include <ida/ida.h>
#include <nvector/nvector_serial.h>
#include <sunmatrix/sunmatrix_dense.h>
#include <sunlinsol/sunlinsol_dense.h>
#include <sunmatrix/sunmatrix_band.h>
#include <sunlinsol/sunlinsol_band.h>
#include <sunnonlinsol/sunnonlinsol_newton.h>


namespace adsim {

  enum class LogLevel {
    LOW,      // Essential: init, cycle start, CSS, failures, completion, statistics
    MEDIUM,   // Detailed: per-step progress, CSS convergence metrics, acceleration
    HIGH      // Diagnostic: IDA step counts, per-step wall-clock time
  };

  static const char* LogLevelList[] = { "Low", "Medium", "High" };

  // Callback for solver messages: (level, message)
  using SolverMessageCallback = std::function<void(LogLevel, const std::string&)>;

  struct IntegrationStats {
    long int num_steps = 0;
    long int num_res_evals = 0;
    long int num_nonlin_iters = 0;
    long int num_err_test_fails = 0;
    long int num_conv_fails = 0;
    long int num_lin_res_evals = 0; // residual evals spent building the DQ Jacobian
    long int num_jac_evals = 0;
    long int num_lin_setups = 0;
  };

  enum class SolverStatus {
    NOT_STARTED,      // Solver hasn't been started yet
    IN_PROGRESS,      // Solver is running (chunk completed, more work to do)
    COMPLETED,        // All cycles completed successfully
    FAILED,           // Solver encountered an error
    AT_CSS
  };

  class ChunkStatus {
  public:
    bool IsActive() const { return active; }
    bool IsStepInitialised() const { return step_initialised; }
    double GetStartTime() const { return step_start_time; }
    double GetEndTime() const { return step_end_time; }
    SolverStatus GetStatus() const { return status; }
    const std::string& GetMessage() const { return message; }

    std::string GetFlagMessage() const noexcept { 
      switch(ida_flag) {
        case IDA_SUCCESS:
          return "Solver completed successfully. ";
        case IDA_TSTOP_RETURN:
          return "Solver reached stop time. ";
        case IDA_ROOT_RETURN:
          return "Solver found a root. ";
        case IDA_WARNING:
          return "Solver warning. ";
        case IDA_TOO_MUCH_WORK:
          return "Solver failed: too much work (exceeded maximum number of steps). ";
        case IDA_TOO_MUCH_ACC:
          return "Solver failed: too much accuracy requested. ";
        case IDA_ERR_FAIL:
          return "Solver failed: error test failure. ";
        case IDA_CONV_FAIL:
          return "Solver failed: convergence test failure (may indicate initial condition calculation failure). ";
        case IDA_LINIT_FAIL:
          return "Solver failed: linear solver initialization failure (may occur during initial condition calculation). ";
        case IDA_LSETUP_FAIL:
          return "Solver failed: linear solver setup failure (may occur during initial condition calculation). ";
        case IDA_LSOLVE_FAIL:
          return "Solver failed: linear solver solve failure (may occur during initial condition calculation). ";
        case IDA_RES_FAIL:
          return "Solver failed: residual function failure (may occur during initial condition calculation). ";
        case IDA_REP_RES_ERR:
          return "Solver failed: repeated residual function errors. ";
        case IDA_RTFUNC_FAIL:
          return "Solver failed: rootfinding function failure. ";
        case IDA_CONSTR_FAIL:
          return "Solver failed: constraint violation (may indicate initial condition calculation failure). ";
        case IDA_FIRST_RES_FAIL:
          return "Solver failed: residual function failed at first call (may indicate initial condition calculation failure). ";
        case IDA_LINESEARCH_FAIL:
          return "Solver failed: line search failure (may indicate initial condition calculation failure). ";
        case IDA_NO_RECOVERY:
          return "Solver failed: no recovery possible (may occur during initial condition calculation). ";
        case IDA_NLS_INIT_FAIL:
          return "Solver failed: nonlinear solver initialization failure. ";
        case IDA_NLS_SETUP_FAIL:
          return "Solver failed: nonlinear solver setup failure. ";
        case IDA_NLS_FAIL:
          return "Solver failed: nonlinear solver failure. ";
        case IDA_MEM_NULL:
          return "Solver failed: memory is NULL. ";
        case IDA_MEM_FAIL:
          return "Solver failed: memory allocation failure. ";
        case IDA_ILL_INPUT:
          return "Solver failed: illegal input (may occur during initial condition calculation with bad icopt, tout1, or id). ";
        case IDA_NO_MALLOC:
          return "Solver failed: memory not allocated. ";
        case IDA_BAD_EWT:
          return "Solver failed: bad error weight vector (may occur during initial condition calculation). ";
        case IDA_BAD_K:
          return "Solver failed: bad value of k. ";
        case IDA_BAD_T:
          return "Solver failed: bad value of t. ";
        case IDA_BAD_DKY:
          return "Solver failed: bad derivative vector. ";
        case IDA_VECTOROP_ERR:
          return "Solver failed: vector operation error. ";
        case IDA_CONTEXT_ERR:
          return "Solver failed: context error. ";
        case IDA_UNRECOGNIZED_ERROR:
          return "Solver failed: unrecognized error. ";
        default:
          return "Solver failed. ";
      }
    }

    void Activate() {
      status = SolverStatus::IN_PROGRESS;
      cycle_index = 0;
      step_index = 0;
      active = true;
      step_initialised = false;
      step_start_time = 0.0;
      step_end_time = 0.0;
      ida_flag = IDA_SUCCESS;
    }

    void StepInitialised(double t_start, double t_final) { 
      step_initialised = true; 
      step_start_time = t_start;
      step_end_time = t_final;
    }

    void Failed(const std::string& message_, int flag) {
      status = SolverStatus::FAILED;
      message = message_;
      cycle_index = 0;
      step_index = 0;
      active = false;
      step_initialised = false;
      step_start_time = 0.0;
      step_end_time = 0.0;
      ida_flag = flag;
    }

    void Interupt() {
      status = SolverStatus::NOT_STARTED;
      message = "Solver interupted. ";
      cycle_index = 0;
      step_index = 0;
      active = false;
      step_initialised = false;
      step_start_time = 0.0;
      step_end_time = 0.0;
    }

    void Completed(const std::string& message_, int flag) {
      status = SolverStatus::COMPLETED;
      message = message_;
      cycle_index = 0;
      step_index = 0;
      active = false;
      step_initialised = false;
      step_start_time = 0.0;
      step_end_time = 0.0;
      ida_flag = flag;
    }

    void ReachedCSS(int at_cycle, int flag) {
      status = SolverStatus::AT_CSS;
      message = "Cyclic steady state reached at cycle " + std::to_string(at_cycle + 1);
      cycle_index = 0;
      step_index = 0;
      active = false;
      step_initialised = false;
      step_start_time = 0.0;
      step_end_time = 0.0;
      ida_flag = flag;
    }

    void CycleComplete() {
      cycle_index++;
      step_index = 0;
      step_initialised = false;
      step_start_time = 0.0;
      step_end_time = 0.0;
    }

    void StepComplete() {
      step_index++;
      step_initialised = false;
      step_start_time = 0.0;
      step_end_time = 0.0;
    }
    
    int GetCycleIndex() const { return cycle_index; }
    int GetStepIndex() const { return step_index; }

  private:
    std::string message; 
    double step_start_time = 0.0;
    double step_end_time = 0.0;
    int cycle_index = 0;
    int step_index = 0; 
    int ida_flag = 0;
    bool active = false;
    bool step_initialised = false;
    SolverStatus status = SolverStatus::NOT_STARTED;
  };

  class SolverData {
  public:
    SolverData();
    ~SolverData();
    // Owns raw SUNDIALS handles; copying would double-free them.
    SolverData(const SolverData&) = delete;
    SolverData& operator=(const SolverData&) = delete;

    // (Re)allocates the integrator. Safe to call repeatedly: anything left from
    // a previous Initialise is released first.
    int Initialise(UserData& user_data, std::string& error_message);
    int Reinitialise(UserData& user_data);
    // Validate a CSS-acceleration re-seed (the current user_data.state) in a
    // throwaway integrator, then commit it to the live integrator only if it
    // re-initialises and integrates cleanly for probe_dt seconds. Returns true if
    // committed (state re-seeded), false if rejected (the live integrator is left
    // untouched and user_data.state is restored from it). Precondition: the live
    // integrator is at the cycle boundary and has NOT been re-initialised.
    // On commit, user_data.state holds IDA's consistent IC (from IDACalcIC). If
    // the live re-initialisation itself fails, the integrator is restarted from
    // the pre-reseed state and false is returned. If `reason` is given it
    // receives why a re-seed was rejected (empty on commit).
    bool ValidateAndCommitReseed(UserData& user_data, double probe_dt,
                                 std::string* reason = nullptr);
    // Test hook: make the next ValidateAndCommitReseed treat its live IDACalcIC
    // as failed (after the scratch validation passed), to exercise the commit
    // path's fallback. Cleared after one use; never set in production.
    void FailNextReseedCommitForTesting() { fail_next_reseed_commit_for_testing_ = true; }
    double GetTime() const;
    int SetStopTime(double t_final);
    int OneStep(double t_final, double& t_curr);
    // Sets IDA's initial step size hin [s]. IDA only reads hin on the first
    // internal step after IDAInit/IDAReInit (nst == 0); later it is ignored.
    int SetInitStep(double h);

    // How PrepareStepRetry set up a retry after a recoverable IDA failure.
    struct StepRetry {
      bool retry = false;     // false: re-entering IDASolve cannot change the outcome
      bool restarted = false; // true: re-initialised at the current point with hin = h
      bool consistent_ic = false; // true: the restart also recomputed a consistent IC
      double h = 0.0;         // step size the retry attempts first [s] (0 = IDA's estimate)
    };
    // Prepares the integrator to retry after IDASolve returned the recoverable
    // failure `flag`. Before the first accepted step (nst == 0) IDA restarts
    // from hin, so the integrator is re-initialised at the (unchanged) current
    // point with hin = half the failed step, clamped to min_step; no retry is
    // offered once a step-size failure already sits at the min_step floor.
    // After an accepted step IDA ignores hin: re-entering IDASolve continues
    // from IDA's own (already reduced) step with a fresh failure budget, so
    // nothing is changed. A step-size failure that has reached the floor there
    // instead restarts the integrator from a consistent initial condition at
    // the last accepted point (RestartFromConsistentIC), with hin = init_step,
    // which discards a history that cannot cross a discontinuity.
    StepRetry PrepareStepRetry(int flag, double min_step);
    // Step size IDA used for its first step since the last IDAInit/IDAReInit.
    double GetActualInitStep() const;
    // Notifies the custom Jacobian (if any) that a new cycle step begins.
    void OnStepStart(int cycle_index) { if (custom_jacobian_) jacobian_.OnStepStart(cycle_index); }
    void Unpack(UserData& user_data);
    // Loads IDA's interpolated state at time t (which must lie within the last
    // accepted step) into user_data.state, with all derived fields updated.
    bool UnpackInterpolated(UserData& user_data, double t);
    IntegrationStats GetStatistics() const;

    // Convergence monitoring (valid after at least one accepted step).
    double GetLastStep() const;   // size of the last internal step [s]
    double GetErrorNorm() const;  // WRMS norm of IDA's estimated local error
    // Which state entry dominates IDA's weighted local error estimate (the
    // variable limiting the step size), as "<name> @ cell <i>".
    std::string DescribeLargestLocalError(const Equations& state) const;
    long int GetNumSteps() const; // cumulative internal steps taken
    const BandedJacobian& GetJacobian() const { return jacobian_; }
    // True when a band audit found a coupling outside the probed matrix band
    // (IDASettings::band_from_probe_) and the band is not yet the full reach.
    bool BandNeedsWidening() const;
    // Rebuilds the band matrix and linear solver at the full 2 neq reach and
    // restarts the integrator from the current point (RestartFromConsistentIC).
    bool WidenBand(UserData& user_data);
    bool BandWidened() const { return band_widened_; }
    sunindextype BandUpper() const { return band_mu_; }
    sunindextype BandLower() const { return band_ml_; }

  private:
    // Releases the IDA memory, N_Vectors, band matrix and linear solver and
    // nulls every handle. Idempotent. Does not touch the SUNContext.
    void Free();
    // Restores hin to the configured initial step (0 = IDA's own estimate).
    // Called before every IDAReInit so a retry's reduced hin cannot leak into a
    // later re-initialisation.
    void ResetInitStep(const IDASettings& settings);
    // Re-initialises IDA at its current (last accepted) point with fresh
    // history: algebraic derivatives zeroed, then IDACalcIC (IDA_YA_YDP_INIT)
    // and IDAGetConsistentIC into y_ / yp_ and the state. Returns false (the
    // integrator is then unusable) if any of these fail.
    bool RestartFromConsistentIC();
    // The configured initial step (0 = IDA's own estimate).
    double GetInitStepSetting() const;

    // Owns the SUNContext. Declared first so it is destroyed last: every
    // SUNDIALS object below (including those held by jacobian_) is created
    // against it and must be released while it is still alive.
    struct ContextOwner {
      SUNContext ctx = nullptr;
      ContextOwner() { SUNContext_Create(SUN_COMM_NULL, &ctx); }
      ~ContextOwner() { if (ctx) SUNContext_Free(&ctx); }
      ContextOwner(const ContextOwner&) = delete;
      ContextOwner& operator=(const ContextOwner&) = delete;
    };
    ContextOwner context_;
    SUNContext sunctx_ = context_.ctx;  // non-owning alias
    void* ida_mem_ = nullptr;
    N_Vector y_ = nullptr;
    N_Vector yp_ = nullptr;
    N_Vector id_ = nullptr;
    N_Vector abstol_vec_ = nullptr;
    N_Vector ewt_ = nullptr;   // scratch: error weights (for GetErrorNorm)
    N_Vector ele_ = nullptr;   // scratch: estimated local errors
    N_Vector y_interp_ = nullptr;   // scratch: interpolated state for fill-in logging
    N_Vector yp_interp_ = nullptr;
    SUNMatrix A_ = nullptr;
    SUNLinearSolver LS_ = nullptr;
    BandedJacobian jacobian_;
    // Half-bandwidths of A_ and the full 2 neq reach (see Initialise / WidenBand).
    sunindextype band_mu_ = 0;
    sunindextype band_ml_ = 0;
    sunindextype band_reach_ = 0;
    bool band_widened_ = false;
    bool custom_jacobian_ = false;
    bool fail_next_reseed_commit_for_testing_ = false;

  };

  // Callback: (cycle_index, step_index, total_cycles, total_steps)
  using SolverProgressCallback = std::function<void(int, int, int, int)>;

  // Per-step convergence metrics: (sim_time, dt, error_norm, cycle_index, step_index)
  using SolverMetricsCallback = std::function<void(double, double, double, int, int)>;

  class CoupledSolver {
  public:
    explicit CoupledSolver(UserData& data);
    int IntegrateCycle(Cycle& cycle, bool verbose = false);
    void PrintStatistics() const;
    IntegrationStats GetStatistics() const { return solver_data_.GetStatistics(); }
    const BandedJacobian& GetJacobian() const { return solver_data_.GetJacobian(); }
    const SolverData& GetSolverData() const { return solver_data_; }
    std::string DescribeLargestLocalError() const { return solver_data_.DescribeLargestLocalError(user_data_.state); }

    void SetProgressCallback(SolverProgressCallback cb) { progress_callback_ = std::move(cb); }
    void SetMessageCallback(SolverMessageCallback cb) { message_callback_ = std::move(cb); }
    void SetMetricsCallback(SolverMetricsCallback cb) { metrics_callback_ = std::move(cb); }
    void SetLogLevel(LogLevel level) { log_level_ = level; }

    void RequestStop() { stop_requested_ = true; }
    bool IsStopRequested() const { return stop_requested_; }

    // CSS outcome of the most recent IntegrateCycle call: the 0-based cycle index
    // at which cyclic steady state was reached (the run stops there, so it may be
    // before number_of_cycles_ - 1), or -1 if the run did not reach CSS.
    int GetCSSCycleIndex() const { return css_cycle_index_; }
    bool ReachedCSS() const { return css_cycle_index_ >= 0; }

    void SaveLogs(const std::string& directory = "Results") const { equation_log_.SaveToFile(directory); }
    const EquationLog& GetLog() const { return equation_log_; }
    const BoundaryLog& GetBoundaryLog() const { return boundary_log_; }
    EquationLog ExtractLog() { return std::move(equation_log_); }
    BoundaryLog ExtractBoundaryLog() { return std::move(boundary_log_); }

  private:
    bool Initialise(Step& first_step, std::string& error_message);
    int IntegrateStep(Step& step, bool verbose = false, int cycle_index = 0, int step_index = 0);
    void Log(Step& step, double t_curr, int step_index, int cycle_index);
    void LogMessage(LogLevel level, const std::string& msg);

    UserData& user_data_;
    SolverData solver_data_;
    EquationLog equation_log_;
    BoundaryLog boundary_log_;
    SolverProgressCallback progress_callback_;
    SolverMessageCallback message_callback_;
    SolverMetricsCallback metrics_callback_;
    LogLevel log_level_ = LogLevel::LOW;
    std::atomic<bool> stop_requested_{false};
    int css_cycle_index_ = -1;

  };

  class CoupledSolverChuncked {
  public:
    explicit CoupledSolverChuncked(UserData& data);
    void Integrate(int max_time_ms);
    void Stop();
    void Start();

    void SetMessageCallback(SolverMessageCallback cb) { message_callback_ = std::move(cb); }
    void SetLogLevel(LogLevel level) { log_level_ = level; }

    bool IsActive() const { return status_.IsActive(); }
    ChunkStatus GetStatus() const { return status_; }
    double GetProgressPercent() const;

    void SaveLogs(const std::string& directory = "Results") const { equation_log_.SaveToFile(directory); }
    const EquationLog& GetLog() const { return equation_log_; }
    const BoundaryLog& GetBoundaryLog() const { return boundary_log_; }
    EquationLog ExtractLog() { return std::move(equation_log_); }
    BoundaryLog ExtractBoundaryLog() { return std::move(boundary_log_); }

  private:
    // bool Initialise(Step& first_step, std::string& error_message);
    void Log(Step& step, double t_curr, int step_index, int cycle_index);
    void ProcessStep(int max_time_ms);
    bool InitialiseStep();
    void LogMessage(LogLevel level, const std::string& msg);
    // Emit the end-of-run wall-clock line plus IDA statistics (mirrors the
    // non-chunked CoupledSolver). Called once when the run reaches CSS, finishes,
    // or fails.
    void LogRunSummary(bool success);

    UserData& user_data_;
    SolverData solver_data_;
    EquationLog equation_log_;
    BoundaryLog boundary_log_;
    ChunkStatus status_;

    SolverMessageCallback message_callback_;
    LogLevel log_level_ = LogLevel::LOW;

    // Timing state for HIGH-level diagnostics. A step spans multiple chunks, so
    // the start markers are captured once (in InitialiseStep) and read back when
    // the step completes (in ProcessStep). run_wall_start_ times the whole run.
    std::chrono::steady_clock::time_point run_wall_start_{};
    std::chrono::steady_clock::time_point step_wall_start_{};
    long int ida_steps_before_step_ = 0;
  };

} // end adsim namespace
