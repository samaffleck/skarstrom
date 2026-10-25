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
#include "adsim/solver/log.h"
#include "adsim/solver/operating_limits.h"
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

    bool Ok() const { return status == RunStatus::Completed || status == RunStatus::ReachedCSS; }
  };

  // KPIs of one tracked component at one logged cycle. A KPI that is not
  // defined for the run (no product drawn, zero denominator) is NaN, and so is
  // every KPI when the component's KPI configuration is incomplete (`valid`
  // false: no feed or no product endpoint).
  struct KPIRecord {
    int cycle = 0;
    std::string component;
    bool valid = false;
    double purity = std::numeric_limits<double>::quiet_NaN();               // [mol/mol]
    double recovery = std::numeric_limits<double>::quiet_NaN();             // [mol/mol]
    double productivity_mass = std::numeric_limits<double>::quiet_NaN();    // [kg/(kg s)]
    double productivity_volume = std::numeric_limits<double>::quiet_NaN();  // [kg/(m3 s)]
    double energy_consumption = std::numeric_limits<double>::quiet_NaN();   // [J/kg]
    double production_rate = std::numeric_limits<double>::quiet_NaN();      // [mol/s]
  };

  // Cyclic-steady-state check for one logged cycle (the deck's CSS method).
  struct CSSRecord {
    int cycle = 0;
    bool converged = false;
    double error = 0.0;
    double tolerance = 0.0;
    std::string method;
  };

  struct TimeSeries {
    std::string variable;
    std::string component;  // empty for a bulk variable
    std::string unit;
    std::vector<double> time;    // [s]
    std::vector<double> values;
  };

  // Axial profiles of one variable over one step: values is row-major
  // [time][cell], with rows == time.size() and cols == z.size().
  struct AxialProfile {
    std::string variable;
    std::string component;
    std::string unit;
    int cycle = 0;
    int step = 0;
    std::vector<double> z;       // cell-centre positions [m]
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
  // "molar_flow_rate", "compressor_power", "mole_fraction", "molar_flow").
  std::vector<std::string> BulkVariableNames();
  std::vector<std::string> ComponentVariableNames();
  std::vector<std::string> StreamVariableNames();
  std::vector<std::string> StreamComponentVariableNames();

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
    // settings.max_cycles_to_keep).
    std::vector<int> Cycles() const;
    int LastCycle() const;  // -1 when nothing is logged
    std::vector<std::string> Components() const;
    std::vector<std::string> StepNames() const;
    std::vector<StreamInfo> Streams() const;
    // Cell-centre positions [m].
    std::vector<double> Positions() const;

    // KPIs of every tracked component at every logged cycle (cycle-major).
    std::vector<KPIRecord> KPIs() const;
    // KPIs of every tracked component at one cycle (-1 = last).
    std::vector<KPIRecord> KPIs(int cycle) const;
    // The deck's CSS check evaluated at every logged cycle.
    std::vector<CSSRecord> CSSHistory();
    // Pressure, wall temperature and fluidisation extremes of one cycle (-1 =
    // last) with the warnings where it leaves the wall design or lifts the bed
    // (operating_limits.h). Empty (cycle -1) without results.
    OperatingLimitsReport OperatingLimits(int cycle = -1) const;
    // The wall's minimum thickness for its pressure design (wall_design.h).
    // Evaluates the deck's expressions first.
    WallDesignResult WallDesign();

    // Time series of a log variable at one point of the bed over one cycle
    // (-1 = last), concatenated across its steps, with time measured from the
    // start of that cycle. `component` is required for component variables and
    // must be empty for bulk ones. z is used with SeriesLocation::Position.
    TimeSeries Series(const std::string& variable, const std::string& component = {},
                      SeriesLocation where = SeriesLocation::Top, double z = 0.0, int cycle = -1) const;
    // Axial profiles of a log variable over one step of one cycle (-1 = last).
    AxialProfile Profile(const std::string& variable, const std::string& component = {},
                         int cycle = -1, int step = 0) const;
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
    int ResolveCycle(int cycle) const;

    std::unique_ptr<UserData> data_;
    std::unique_ptr<ParameterRegistry> registry_;
    std::unique_ptr<CoupledSolver> solver_;
    std::unique_ptr<CoupledSolverChuncked> chunked_;
    RunResult result_;
    std::vector<SolverMessage> messages_;
    LogLevel chunked_log_level_ = LogLevel::HIGH;
    SolverMessageCallback chunked_on_message_;
  };

}  // namespace adsim
