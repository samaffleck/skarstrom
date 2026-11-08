/*
 * Skarstrom - pressure and temperature swing adsorption simulator
 * Copyright (c) 2024-2026 Sam Affleck and contributors
 *
 * Licensed under the terms in the LICENSE file at the repository root.
 */
#pragma once

// Simulation: a thin facade over one deck and its solver.
//
//   adsim::Simulation sim = adsim::Simulation::FromFile("deck.json");
//   sim.SetParameter("cycle.number_of_cycles", 5);
//   sim.SetParameter("step[Feed].time", 120.0);          // alias or ID path
//   adsim::RunResult r = sim.Run();
//   for (const auto& k : sim.KPIs()) ...                    // per cycle, per component
//   adsim::TimeSeries y = sim.Series("mole_fraction", "CarbonDioxide",
//                                    adsim::SeriesLocation::Top);
//
// Load or build a deck, set parameters through the ParameterRegistry, run, then
// read typed results. Two ways to run:
//
//   * Run(): blocking, with the full CoupledSolver (what the tests, the desktop
//     app and the Python package use). An optional progress callback is told
//     after every step and can cancel; a cheap cancel poll runs after every
//     internal solver step.
//   * Start() / Advance(ms) / Stop(): cooperative chunks with
//     CoupledSolverChuncked, for hosts that must yield (the browser worker).
//
// Results come from whichever solver ran last. They stay valid until the next
// Load*, Run or Start. The deck lives on the heap, so a Simulation can be
// moved while the solver holds a reference to the deck.

#include <functional>
#include <limits>
#include <memory>
#include <optional>
#include <string>
#include <string_view>
#include <vector>

#include "adsim/objects/user_data.h"
#include "adsim/parameters/parameter_registry.h"
#include "adsim/solver/coupled_solver.h"
#include "adsim/solver/cycle_state.h"
#include "adsim/solver/kpi_history.h"  // KPIRecord
#include "adsim/solver/log.h"
#include "adsim/solver/operating_limits.h"
#include "adsim/solver/run_health.h"
#include "adsim/objects/wall_design.h"

namespace adsim {

  enum class RunStatus {
    NotRun,      // nothing has run yet (or the deck was reloaded)
    Running,     // a chunked run is in progress
    Completed,   // every configured cycle ran without reaching CSS
    ReachedCSS,  // stopped early at cyclic steady state
    Failed,      // checks, initialisation or integration failed (see message)
    Cancelled    // a progress callback, a cancel poll or Stop() ended the run
  };

  const char* RunStatusName(RunStatus status);

  struct RunProgress {
    int cycle = 0;         // 0-based cycle that just finished a step
    int step = 0;          // 0-based step that just finished
    int total_cycles = 0;  // configured number of cycles (an upper bound: CSS may stop earlier)
    int total_steps = 0;   // steps per cycle
    double fraction = 0.0; // completed steps / (total_cycles * total_steps)
  };

  // Return false to cancel the run.
  using RunProgressCallback = std::function<bool(const RunProgress&)>;
  // Polled after every internal integrator step; return true to cancel. Must
  // be cheap.
  using CancelPoll = std::function<bool()>;

  struct RunOptions {
    // Overrides the deck's number of cycles for this run (and leaves it set).
    std::optional<int> cycles;
    RunProgressCallback progress;
    CancelPoll cancel;
    // Solver messages up to this level are kept (Simulation::Messages) and
    // passed to on_message.
    LogLevel log_level = LogLevel::LOW;
    SolverMessageCallback on_message;
    // Called after every accepted internal step with (time, step size, error
    // norm, cycle, step) -- for diagnostics (bench_solver --trace).
    SolverMetricsCallback on_step;
    // Start this run from a state (cycle_state.h) instead of the deck's
    // initial conditions; overrides SetInitialState for this run only. The
    // state is adapted to the deck first (PrepareCycleState: remapped onto
    // other cells per layer, histories stretched onto changed step times);
    // a state that does not fit fails the run with the reason. Must outlive
    // the Run call. nullptr: the warm-start state if set, else the initial
    // conditions.
    const CycleState* initial_state = nullptr;
  };

  struct SolverMessage {
    LogLevel level = LogLevel::LOW;
    std::string text;
  };

  struct RunResult {
    RunStatus status = RunStatus::NotRun;
    std::string message;       // failure / completion text
    int return_code = 0;       // CoupledSolver::IntegrateCycle's return value (blocking runs)
    bool css_reached = false;
    int css_cycle = -1;        // 0-based cycle at which CSS was reached, or -1
    int cycles_run = 0;        // cycles that completed
    double wall_time_s = 0.0;
    IntegrationStats stats{};  // IDA counters (blocking runs)
    // Run-health checks (run_health.h): false when any fired. They never change
    // the status or the results; the warnings say what to look at.
    bool healthy = true;
    std::vector<RunHealthWarning> health_warnings;

    bool Ok() const { return status == RunStatus::Completed || status == RunStatus::ReachedCSS; }
  };

  // Cyclic-steady-state check for one cycle (the deck's CSS method).
  struct CSSRecord {
    int cycle = 0;
    bool converged = false;
    double error = 0.0;
    double tolerance = 0.0;
    std::string method;
    // The end-of-cycle state change and its tolerance (mass-balance methods;
    // tolerance 0 = not checked), and the consecutive passing cycles ending here
    // out of the number required (CSSResult).
    double state_change = 0.0;
    double state_tolerance = 0.0;
    int passing_cycles = 0;
    int required_cycles = 1;
    // Diagnostic (never part of the stopping rule): the estimated remaining
    // distance to the cyclic steady state, in the norm of state_change, its
    // ratio to state_tolerance and the convergence rate it came from
    // (CyclicSteadyState::EstimateDistance). NaN when not estimated; then
    // distance_note says why.
    bool distance_estimated = false;
    double convergence_rate = std::numeric_limits<double>::quiet_NaN();
    double distance = std::numeric_limits<double>::quiet_NaN();
    double distance_ratio = std::numeric_limits<double>::quiet_NaN();
    std::string distance_note;
  };

  struct TimeSeries {
    std::string variable;
    std::string component;  // empty for a bulk variable
    std::string adsorbent;  // the adsorbent's label for a per-adsorbent series, else empty
    std::string reaction;   // the reaction's name for a reaction series, else empty
    std::string unit;
    std::vector<double> time;    // [s]
    std::vector<double> values;
  };

  // Axial profiles of one variable over one step: values is row-major
  // [time][point], with rows == time.size() and cols == z.size(). A cell
  // variable has one point per bed cell, at its centre (the ghost cells that
  // carry the boundary conditions are left out); a face variable (velocity,
  // volumetric and molar flow rate, approach to fluidisation) one per face,
  // from the bottom of the column (z = 0) to the top (z = L).
  struct AxialProfile {
    std::string variable;
    std::string component;
    std::string adsorbent;  // the adsorbent's label for a per-adsorbent profile, else empty
    std::string reaction;   // the reaction's name for a reaction profile, else empty
    std::string unit;
    int cycle = 0;
    int step = 0;
    bool faces = false;          // true: z are face positions, false: cell centres
    std::vector<double> z;       // [m] from the bottom of the column
    std::vector<double> time;    // [s], absolute solver time
    std::vector<double> values;
    size_t rows = 0;
    size_t cols = 0;
  };

  // One endpoint logged by the BoundaryLog (a flow stream, or a pressure
  // boundary named "<step> Pressure Boundary Top|Bottom").
  struct StreamInfo {
    int step = 0;              // step index
    std::string step_name;
    std::string name;          // the BoundaryLog name StreamSeries takes
    size_t endpoint_id = 0;    // the endpoint's stable ID ("endpoint[<id>]" paths)
    BoundaryLocation location = BoundaryLocation::BOTTOM;
  };

  enum class SeriesLocation { Top, Bottom, Position };

  // Names accepted by Series / Profile (snake_case of the log variables, e.g.
  // "pressure", "temperature", "velocity", "mole_fraction",
  // "adsorbed_concentration") and by StreamSeries ("temperature",
  // "molar_flow_rate", "compressor_power", "pressure", "mole_fraction",
  // "molar_flow").
  std::vector<std::string> BulkVariableNames();
  std::vector<std::string> ComponentVariableNames();
  std::vector<std::string> StreamVariableNames();
  std::vector<std::string> StreamComponentVariableNames();
  // The component variables Series / Profile also give per adsorbent (their
  // `adsorbent` argument): "adsorbed_concentration",
  // "equilibrium_adsorbed_concentration", "mass_transfer_coefficient",
  // "slow_adsorbed_concentration". The aggregate (no adsorbent) is the
  // mass-weighted sum over the solids of each cell; a per-adsorbent series is
  // one solid's own value [per kg of it], 0 in the cells without it.
  std::vector<std::string> PerAdsorbentVariableNames();
  // The variables Series / Profile give per reaction (their `reaction`
  // argument): "reaction_rate", eta r of one reaction on one solid
  // [mol/(kg of that solid s)], 0 in the cells without the solid. The
  // reactions' total effect is "reaction_heat" (bulk) and "reaction_source"
  // (per component), which take no reaction.
  std::vector<std::string> ReactionVariableNames();

  // One reaction on one solid whose rate the last run logged.
  struct SeriesReaction {
    size_t id = 0;               // the reaction's ID (reaction[ID] paths)
    std::string name;            // the reaction's name: the `reaction` argument of Series / Profile
    size_t adsorbent_id = 0;     // the solid that carries it
    std::string adsorbent;       // that solid's series label (SeriesAdsorbent::label)
    std::string variable_name;   // the log name, "<reaction> Rate (<adsorbent>)"
  };

  // One adsorbent of the last run's bed, as the per-adsorbent series name it.
  struct SeriesAdsorbent {
    size_t id = 0;       // the adsorbent's deck ID
    std::string label;   // its name (" #<id>" added when two share a name): the
                         // `adsorbent` argument of Series / Profile and the
                         // bracket in the log names, "<component> Adsorbed
                         // Concentration (<label>)"
    bool logged = false; // per-adsorbent series are in the log (the bed has a
                         // mixed layer); otherwise they are read off the aggregate
  };

  // Loads a deck from JSON text in either root shape ("simulation", written by
  // SaveToString, or cereal's default "value0", written by SaveToJson). Throws
  // and leaves `data` unchanged on failure.
  void LoadFromJsonText(UserData& data, const std::string& json);

  class Simulation {
  public:
    Simulation();
    ~Simulation();
    Simulation(Simulation&&) noexcept;
    Simulation& operator=(Simulation&&) noexcept;
    Simulation(const Simulation&) = delete;
    Simulation& operator=(const Simulation&) = delete;

    static Simulation FromFile(const std::string& path);
    static Simulation FromJson(const std::string& json);

    // ---- deck -------------------------------------------------------------
    void Load(const std::string& path);
    void LoadJson(const std::string& json);
    // Replaces the deck (e.g. one built in code or from the web document).
    // The results of the previous run stay readable until the next run.
    // Throws std::logic_error during a chunked run.
    void SetData(UserData&& data);
    // Writes the deck as cereal JSON (the "value0" file format).
    void Save(const std::string& path);
    // The deck as cereal JSON text (the "simulation" string format).
    std::string ToJson();

    UserData& Data() { return *data_; }
    const UserData& Data() const { return *data_; }

    // Re-evaluates every expression and runs the deck checks. Returns false
    // with the reason in `error` when the deck cannot run.
    bool Check(std::string& error);

    // ---- parameters ---------------------------------------------------------
    // The registry for the current deck (rebuilt lazily after structural
    // changes made through this class; call RefreshParameters after editing
    // Data() structurally).
    ParameterRegistry& Parameters();
    void RefreshParameters();
    double GetParameter(std::string_view path);
    std::string GetParameterText(std::string_view path);
    // Numeric set; re-evaluates expressions so dependent values follow.
    void SetParameter(std::string_view path, double value);
    // Text set: an expression ("2*feed_flow", "5 [bar]"), a number, or text
    // for a String parameter. Throws std::invalid_argument when an expression
    // does not evaluate.
    void SetParameterText(std::string_view path, const std::string& text);

    // ---- run ----------------------------------------------------------------
    RunResult Run(const RunOptions& options = {});
    // Asks a blocking Run on another thread to stop after the current internal
    // step (thread-safe).
    void RequestStop();

    // Chunked run. Start() returns false (and sets the result to Failed) when
    // the deck fails its checks.
    bool Start(LogLevel log_level = LogLevel::HIGH, SolverMessageCallback on_message = {});
    // Integrates for up to max_time_ms of wall clock. Returns the run status.
    RunStatus Advance(int max_time_ms);
    void Stop();
    bool IsRunning() const;
    double ProgressPercent() const;

    // ---- cycle state (cycle_state.h) ---------------------------------------
    // The state at the end of the last completed blocking run (Run or
    // RunCycle that finished or reached CSS): the DAE vector, the endpoint
    // histories and the buffer-tank summaries, enough to restart from it.
    // Kept for the session only; it is never saved in the deck. Cleared by the
    // next run, Load* and a failed or cancelled run.
    bool HasFinalState() const { return final_state_ != nullptr; }
    // Throws std::logic_error when there is none.
    const CycleState& FinalState() const;

    // Warm start: every later blocking run (Run, RunCycle without a state)
    // starts from `state` instead of the deck's initial conditions, until
    // ClearInitialState or Load*. The deck's CSS minimum number of cycles
    // counts the cycles behind the state (CycleState::cycles), so a run from a
    // state at CSS can stop as soon as the consecutive-cycles rule holds.
    // Throws std::invalid_argument with the reason when the state does not fit
    // the deck (CheckState); the layout is checked again when the run starts.
    // The chunked run (Start / Advance) always starts from the initial
    // conditions.
    void SetInitialState(CycleState state);
    void ClearInitialState() { initial_state_.reset(); }
    bool HasInitialState() const { return initial_state_ != nullptr; }
    const CycleState& InitialState() const;
    // Why `state` cannot start a run of this deck, or an empty string. Checks
    // the components, layers and void regions, and the step and endpoint
    // structure. Other cell counts and step times are fine: a run adapts the
    // state (PrepareState).
    std::string CheckState(const CycleState& state);
    // The state as a run would start from it: remapped onto this deck's cells
    // per layer and with its histories stretched onto the deck's step times
    // (PrepareCycleState in cycle_state.h). `note`, when given, says what was
    // done. Throws std::invalid_argument when the state does not fit.
    CycleState PrepareState(const CycleState& state, std::string* note = nullptr);
    // The state transferred onto this deck's cells per layer (RemapCycleState
    // in cycle_state.h: a conservative per-layer remap that keeps every
    // component's inventory in every layer). Runs started from a state on
    // another grid do this themselves. Throws std::invalid_argument when the
    // state does not fit the deck.
    CycleState RemapState(const CycleState& state);

    // Runs exactly one cycle, from `from`, without changing the deck's number
    // of cycles; FinalState() is then the cycle map applied to `from`. Without
    // a state: from the last run's final state if there is one, else from the
    // warm-start state if set, else from the initial conditions -- so repeated
    // RunCycle() calls step the cycle map; options.initial_state, when set,
    // comes first. options.cycles is ignored.
    RunResult RunCycle(const RunOptions& options = {});
    RunResult RunCycle(const CycleState& from, const RunOptions& options = {});

    const RunResult& LastResult() const { return result_; }
    const std::vector<SolverMessage>& Messages() const { return messages_; }
    // The solver of the last run: the chunked one when it was chunked, the
    // blocking one otherwise (else nullptr). For diagnostics such as the
    // Jacobian statistics bench_solver prints.
    CoupledSolverChuncked* ChunkedSolver() { return chunked_.get(); }
    CoupledSolver* BlockingSolver() { return solver_.get(); }

    // ---- results --------------------------------------------------------------
    bool HasResults() const;
    const EquationLog& Log() const;
    const BoundaryLog& Boundaries() const;
    // Logged cycle indices (0-based, ascending; old cycles are evicted beyond
    // settings.max_cycles_to_keep, or beyond the cycles the CSS check reads if
    // that is more: CyclicSteadyState::CyclesToLog).
    std::vector<int> Cycles() const;
    int LastCycle() const;  // -1 when nothing is logged
    std::vector<std::string> Components() const;
    std::vector<std::string> StepNames() const;
    std::vector<StreamInfo> Streams() const;
    // Centres of the bed cells [m] from the bottom of the column (no ghost
    // cells): the z of a cell variable's Profile.
    std::vector<double> Positions() const;
    // Cell-face positions [m], from 0 (bottom) to the column length (top):
    // the z of a face variable's Profile.
    std::vector<double> FacePositions() const;

    // KPIs of every tracked component at every cycle of the run (cycle-major),
    // not only the logged ones: a cycle still in the log is computed from it,
    // an evicted one comes from the record the solver made when the cycle
    // completed (KPIHistory, with the run's KPI configuration).
    std::vector<KPIRecord> KPIs() const;
    // KPIs of every tracked component at one cycle (-1 = the last logged one).
    // Throws std::out_of_range for a cycle that is neither logged nor recorded.
    std::vector<KPIRecord> KPIs(int cycle) const;
    // The cycles KPIs() covers: every completed cycle plus the logged ones
    // (ascending; a run stopped mid-cycle logs that cycle without completing it).
    std::vector<int> KPICycles() const;
    // The KPIs the solver recorded at the end of each completed cycle.
    const KPIHistory& KPIHistoryRecords() const;
    // The deck's CSS check evaluated at every logged cycle. The consecutive-
    // cycles rule only sees logged cycles, so the earliest logged cycles can
    // show fewer passing cycles than the run counted.
    std::vector<CSSRecord> CSSHistory();
    // Pressure, wall temperature and fluidisation extremes of one cycle (-1 =
    // last) with the warnings where it leaves the wall design or lifts the bed
    // (operating_limits.h). Empty (cycle -1) without results.
    OperatingLimitsReport OperatingLimits(int cycle = -1) const;
    // The run-health measures of every completed cycle of the last run and its
    // warnings (run_health.h). Empty without a run.
    const RunHealthReport& RunHealth() const;
    // The wall's minimum thickness for its pressure design (wall_design.h).
    // Evaluates the deck's expressions first.
    WallDesignResult WallDesign();

    // Time series of a log variable at one point of the bed over one cycle
    // (-1 = last), concatenated across its steps, with time measured from the
    // start of that cycle. `component` is required for component variables and
    // must be empty for bulk ones. z [m] (from the bottom of the column) is
    // used with SeriesLocation::Position. Top and Bottom are the physical
    // boundary values (bed_series.h): for pressure, temperature and mole
    // fraction the endpoint's logged value where the step has an endpoint at
    // that end (at an inlet, the feed composition and temperature or the set
    // pressure; at an outlet, the bed-end face value; the same as
    // StreamSeries), otherwise, and for every other variable, the face value.
    // `adsorbent` (a label of SeriesAdsorbents(), or a deck ID as text) gives
    // that adsorbent's own series of a PerAdsorbentVariableNames() variable;
    // empty gives the aggregate. `reaction` (a reaction's name or its ID as
    // text, SeriesReactions()) is required for a ReactionVariableNames()
    // variable and refused for the others; `adsorbent` then picks the solid
    // when the reaction is on more than one. A reaction rate has no
    // boundary value: Top and Bottom give the last and first packed cell.
    TimeSeries Series(const std::string& variable, const std::string& component = {},
                      SeriesLocation where = SeriesLocation::Top, double z = 0.0, int cycle = -1,
                      const std::string& adsorbent = {}, const std::string& reaction = {}) const;
    // Axial profiles of a log variable over one step of one cycle (-1 = last).
    AxialProfile Profile(const std::string& variable, const std::string& component = {},
                         int cycle = -1, int step = 0, const std::string& adsorbent = {},
                         const std::string& reaction = {}) const;
    // The reactions whose rates the last run logged (empty without results
    // or without reactions).
    std::vector<SeriesReaction> SeriesReactions() const;
    // The adsorbents of the last run's bed with their series labels (empty
    // without results).
    std::vector<SeriesAdsorbent> SeriesAdsorbents() const;
    // A boundary-stream log variable over one step of one cycle (-1 = last).
    // Time is measured from the start of the step.
    TimeSeries StreamSeries(int step, const std::string& stream, const std::string& variable,
                            const std::string& component = {}, int cycle = -1) const;
    // The same, addressing the endpoint by its stable ID (StreamInfo::endpoint_id).
    TimeSeries EndpointSeries(size_t endpoint_id, const std::string& variable,
                              const std::string& component = {}, int cycle = -1) const;

    // Writes the equation log as CSV files into `directory`.
    void SaveLogs(const std::string& directory) const;

  private:
    void ResetResults();
    // Run from `start` (nullptr: initial conditions).
    RunResult RunFrom(const RunOptions& options, const CycleState* start);
    int ResolveCycle(int cycle) const;
    // The adsorbent ID named by `adsorbent` (a label of SeriesAdsorbents() or
    // a deck ID); throws std::invalid_argument when it names none.
    size_t ResolveAdsorbent(const std::string& adsorbent) const;
    // The logged reaction named by `reaction` (name or ID) on `adsorbent`
    // (optional when only one solid carries it); throws std::invalid_argument.
    SeriesReaction ResolveReaction(const std::string& reaction, const std::string& adsorbent) const;

    std::unique_ptr<UserData> data_;
    std::unique_ptr<ParameterRegistry> registry_;
    std::unique_ptr<CoupledSolver> solver_;
    std::unique_ptr<CoupledSolverChuncked> chunked_;
    RunResult result_;
    std::vector<SolverMessage> messages_;
    LogLevel chunked_log_level_ = LogLevel::HIGH;
    SolverMessageCallback chunked_on_message_;
    std::unique_ptr<CycleState> final_state_;
    std::unique_ptr<CycleState> initial_state_;
  };

}  // namespace adsim
