/*
 * Skarstrom - pressure and temperature swing adsorption simulator
 * Copyright (c) 2024-2026 Sam Affleck and contributors
 *
 * Licensed under the terms in the LICENSE file at the repository root.
 */
#pragma once

// Calibration: fit model parameters to measured data.
//
// A calibration problem is a set of CASES (a deck, the datasets of that deck
// to reproduce, how to run it) and a set of PARAMETERS (registry paths with
// bounds). A parameter is SHARED (one value for every case: an isotherm
// capacity, a wall temperature) or PER CASE (one value in each case: a film
// coefficient that depends on the flow). The datasets are the deck's own
// experimental data, each bound to the quantity it measures
// (ObservableBinding, experimental_data.h); the simulated counterpart is read
// with SimulatedSeries (parameterisation.h).
//
//   using namespace adsim::calibration;
//   std::vector<Case> cases = {Case::FromFile("run_50ccm.json"), Case::FromFile("run_100ccm.json")};
//   std::vector<Parameter> params = {
//     Shared("adsorbent[1].isotherm[CarbonDioxide].m2", 1.0, 9.0),
//     PerCase("step[1].medium.heat_transfer_coefficient", 1.0, 500.0, Scale::Log)};
//   Problem problem(cases, params);                  // validates against every deck
//   Result r = Calibrate(problem);                   // least squares by default
//   for (const auto& p : r.parameters) ... p.estimate, p.standard_error ...
//
// ---- the objective -----------------------------------------------------------
//
// Weighted least squares on the residual vector, one entry per measured point:
//
//   r_k = w_d * (y_sim(x_k) - y_k) / s_k,     w_d = sqrt(weight_d / n_d)
//
// for point k of dataset d, where
//   * y_sim(x_k) is the simulated counterpart at the measured x (linear
//     interpolation of a time series; a KPI matched by cycle number), in the
//     dataset's own y unit;
//   * s_k is the point's uncertainty (ObservableBinding::uncertainty) when the
//     dataset gives one, else the range (max - min) of the dataset's y values,
//     so datasets in different units are comparable (1 for a flat dataset);
//   * weight_d is the dataset's weight (DatasetRef::weight, else
//     ObservableBinding::weight), and n_d its number of points: dividing by
//     n_d makes a dataset's share of the objective independent of how densely
//     it was sampled.
// The objective is the sum of squares S = sum_k r_k^2.
//
// The points of a time series that the simulated series does not cover in the
// first successful evaluation (the starting point, for least squares) -
// measured after the end of the simulated cycle - are left out for the whole
// fit and counted (DatasetFit::skipped); the point set never changes during a
// fit. A KPI measured at a cycle the run did not reach (it
// stopped at cyclic steady state earlier) takes the last cycle's value.
//
// ---- the solvers -------------------------------------------------------------
//
// LeastSquares (the default): a bounded Levenberg-Marquardt method on the
// residual vector. Each parameter value is worked in its scale (the value, or
// its logarithm for Scale::Log) normalised by its bounds to [0, 1]; steps are
// projected onto the box and a parameter held at a bound by the gradient is
// left out of the step (a projected, active-set LM), so bounds are respected
// exactly and the standard errors stay finite at a bound. The Jacobian is by
// finite differences (forward or central, the step a fraction of the
// normalised range, one-sided inwards at a bound). Each Jacobian column is one
// simulation per case the parameter affects (a per-case parameter only moves
// its own case), and the columns run in parallel threads. For strongly
// correlated parameters (two that the data mostly fix in combination, such as
// an isotherm capacity and affinity) use central differences: the forward
// difference's error can stall the fit in the narrow valley.
//
// Global: dlib's MaxLIPO + trust-region search (GlobalOptimizer,
// optimisation.h) on S over the same normalised box, the cases of one
// candidate in parallel. With Options::refine the best point found is then
// polished by LeastSquares with the remaining budget.
//
// Every run is made on a fresh copy of its case's deck. With
// Options::warm_start a run starts from the cyclic steady state of the
// nearest point evaluated so far for that case (CycleStateCache,
// state_cache.h; only CSS states are cached, so a breakthrough deck always
// runs cold). A warm start moves where a run stops within the CSS tolerance,
// which is noise in a finite difference: for cyclic decks tighten the
// tolerance in the case's RunPolicy::settings.
//
// A run that fails (or does not reach CSS when its policy requires it) never
// throws mid-fit: it is counted (Result::failed_simulations), a trial step
// that fails is rejected, a Jacobian column that fails is retried on the
// other side and zeroed if that fails too (with a note), and the global
// search scores it kOptimisationPenaltyCost. Only a failed run at the
// starting point stops a least-squares fit (StopReason::StartFailed).
//
// Determinism: the same problem, options and thread count give the same
// result. The runs of a batch are fixed before it starts, every warm start is
// chosen from the cache as it was before the batch, and the states are cached
// after it in batch order; nothing depends on time or on which thread ran
// what. (Thread count only matters through the warm starts: a batch is the
// same whatever the count, so in practice the results agree across counts.)
//
// ---- standard errors ---------------------------------------------------------
//
// From the Jacobian J at the estimate (in each parameter's scale):
//   cov = sigma^2 (J^T J)^-1,   sigma^2 = S / (N - p)
// with N residuals and p estimated values. sigma^2 is estimated from the fit
// even when the data carry uncertainties (it then tests them: sigma^2 near 1
// means the residuals match the stated uncertainty). For a Scale::Log
// parameter the covariance is in ln(value); its standard error in model units
// is value * se(ln value) (the delta method), and its correlations are those
// of ln(value). When N <= p, or J^T J is singular (a parameter that does not
// change the residuals, or two that change them in the same way), the
// standard errors and correlations are NaN and Result::notes says why. A
// parameter at a bound has a standard error from the same formula, which is
// only a local statement: check at_lower / at_upper.

#include <cmath>
#include <cstddef>
#include <functional>
#include <limits>
#include <string>
#include <vector>

#include "adsim/simulation.h"                 // CancelPoll, RunStatus
#include "adsim/solver/parameterisation.h"    // SimulatedCounterpart
#include "adsim/solver/sensitivity_study.h"   // sensitivity::RunPolicy
#include "adsim/solver/state_cache.h"         // CycleStateCache

namespace adsim {
namespace calibration {

inline constexpr double kNaN = std::numeric_limits<double>::quiet_NaN();

// ---- the problem -------------------------------------------------------------

enum class Scale { Linear, Log };
enum class Sharing { Shared, PerCase };

const char* ScaleName(Scale s);      // "linear", "log"
const char* SharingName(Sharing s);  // "shared", "per_case"

// A dataset of a case's deck, by name (the first dataset with that name) or,
// when `name` is empty, by its 0-based index in the deck's experimental data.
// `weight` (> 0) overrides the dataset's own ObservableBinding::weight; NaN
// keeps it.
struct DatasetRef {
  std::string name;
  int index = -1;
  double weight = kNaN;
};

DatasetRef Dataset(std::string name, double weight = kNaN);
DatasetRef Dataset(int index, double weight = kNaN);

// One experiment: a deck and the datasets in it to reproduce.
struct Case {
  // A label for messages and results (default: "case <i>", 1-based).
  std::string name;
  // The deck, as JSON text (the "value0" file form or the "simulation"
  // string form).
  std::string deck_json;
  // The datasets to fit; empty: every bound dataset of the deck.
  std::vector<DatasetRef> datasets;
  // How each run is made: a cycle cap (0: the deck's number of cycles),
  // whether a run must reach cyclic steady state to count, and registry text
  // sets applied before the parameters. Unlike a sensitivity study, CSS is
  // not required by default: calibration data are often breakthrough runs.
  sensitivity::RunPolicy policy = {false, 0, {}};

  // A case from a deck file. Throws std::runtime_error when it cannot be read.
  static Case FromFile(const std::string& path, std::vector<DatasetRef> datasets = {});
  // A case from a deck in memory (serialised now; later edits to `data` do
  // not reach the case). Saving bumps data.info.sim_version, as every save does.
  static Case FromData(UserData& data, std::vector<DatasetRef> datasets = {});
};

// Bounds for one case of a per-case parameter (NaN: the parameter's own).
struct Bounds {
  double lower = kNaN;
  double upper = kNaN;
};

// One parameter to estimate: a tunable registry path or alias
// (parameters/parameter_registry.h), its bounds and initial value in the
// parameter's SI unit.
struct Parameter {
  std::string path;
  double lower = kNaN;  // finite, lower < upper (and lower > 0 for Scale::Log)
  double upper = kNaN;
  Scale scale = Scale::Linear;
  Sharing sharing = Sharing::Shared;
  // The starting value; NaN: the deck's value (the first case's for a shared
  // parameter, each case's own for a per-case one), moved inside the bounds
  // when the deck's value is outside them (with a note).
  double initial = kNaN;
  // PerCase only, optional: one entry per case (empty: the bounds and initial
  // value above for every case). A NaN entry falls back to the value above.
  std::vector<Bounds> case_bounds;
  std::vector<double> case_initial;
};

Parameter Shared(std::string path, double lower, double upper, Scale scale = Scale::Linear,
                 double initial = kNaN);
Parameter PerCase(std::string path, double lower, double upper, Scale scale = Scale::Linear,
                  double initial = kNaN);

// One estimated value: a shared parameter, or a per-case parameter in one case.
struct Unknown {
  size_t parameter = 0;     // index into Problem::Parameters()
  int case_index = -1;      // -1: shared
  std::string path;         // the parameter's path as given
  std::string label;        // the path, with " [<case name>]" for a per-case value
  double lower = kNaN;      // bounds and start in the parameter's SI unit
  double upper = kNaN;
  double initial = kNaN;
  Scale scale = Scale::Linear;
};

// A dataset as the problem uses it: a copy of the deck's dataset with the
// weight that applies.
struct ProblemDataset {
  size_t case_index = 0;
  int deck_index = 0;       // index in the deck's experimental data
  XYData data;              // the points and binding (observable.weight is the weight in use)
  double scale = 1.0;       // s_k when the dataset gives no uncertainty: its y range (1 when flat)
};

// A validated calibration problem. The constructor loads every case's deck,
// applies its policy settings and checks: at least one case and parameter;
// every parameter is a tunable parameter of every case's deck (an unknown path
// names the case), no path twice, finite bounds lower < upper (> 0 for a log
// scale) and a start inside them; per-case bounds and starts of the right
// length; every dataset exists in its deck, is bound to an observable, has
// points and a positive weight. Throws std::invalid_argument with the reason.
class Problem {
 public:
  Problem(std::vector<Case> cases, std::vector<Parameter> parameters);

  const std::vector<Case>& Cases() const { return cases_; }
  const std::vector<Parameter>& Parameters() const { return parameters_; }
  // The values estimated, in order: for each parameter, one entry if shared,
  // else one per case.
  const std::vector<Unknown>& Unknowns() const { return unknowns_; }
  const std::vector<ProblemDataset>& Datasets() const { return datasets_; }
  // Notes from validation (a deck value moved inside the bounds).
  const std::vector<std::string>& Notes() const { return notes_; }

  // The index into Unknowns() of parameter `p` in case `c`.
  size_t UnknownIndex(size_t p, size_t c) const;

 private:
  std::vector<Case> cases_;
  std::vector<Parameter> parameters_;
  std::vector<Unknown> unknowns_;
  std::vector<ProblemDataset> datasets_;
  std::vector<std::string> notes_;
  std::vector<std::vector<size_t>> unknown_of_;  // [parameter][case]
};

// ---- options -----------------------------------------------------------------

enum class Solver { LeastSquares, Global };
enum class Difference { Forward, Central };

const char* SolverName(Solver s);          // "least_squares", "global"
const char* DifferenceName(Difference d);  // "forward", "central"

// Reported after every iteration of the least-squares solver (an accepted or
// rejected step) and after every candidate of the global search.
struct Progress {
  Solver solver = Solver::LeastSquares;
  int iteration = 0;            // 1-based within the solver
  double objective = kNaN;      // S at the current (best) point
  std::vector<double> values;   // that point, per Unknown, in SI units
  int simulations = 0;          // runs so far (all solvers)
};
using ProgressCallback = std::function<void(const Progress&)>;

struct Options {
  Solver solver = Solver::LeastSquares;
  // Global only: polish the best point with LeastSquares afterwards.
  bool refine = false;

  // The run budget, in simulations (one case run is one simulation). A solver
  // stops before a batch that would exceed it. The Jacobian for the standard
  // errors (at a global optimum, or after a least-squares stop on the budget
  // right after an accepted step) is made even then, and counted.
  int max_simulations = 500;
  // LeastSquares: iterations (accepted or rejected steps).
  int max_iterations = 50;
  // LeastSquares: stop when an accepted step lowers S by less than
  // ftol * S, or moves no normalised value by more than xtol.
  double ftol = 1e-6;
  double xtol = 1e-6;
  Difference difference = Difference::Forward;
  // Finite-difference step as a fraction of the normalised range (the bounds,
  // in the parameter's scale).
  double fd_step = 1e-3;

  // Global: the candidates of the search (each is one run per case).
  int global_evaluations = 30;

  // Threads for the runs of a batch (Jacobian columns, the cases of a
  // candidate). 0: the hardware's concurrency - 1, at least 1. Where threads
  // cannot be started (a WebAssembly build without threads) the runs are made
  // one after another.
  int threads = 0;
  // Start runs from the nearest cached cyclic steady state of their case.
  bool warm_start = true;
  size_t cache_capacity = CycleStateCache::kDefaultCapacity;

  // Make a Jacobian at the end when the solver has none at its estimate.
  bool standard_errors = true;

  // Polled by every run (from several threads at once: it must be
  // thread-safe, e.g. read an atomic flag) and between batches. A cancelled
  // fit returns the best point so far.
  CancelPoll cancel;
  ProgressCallback progress;
};

// ---- the result --------------------------------------------------------------

enum class StopReason {
  Converged,      // S changed by less than ftol relative
  SmallStep,      // the step was below xtol
  NoImprovement,  // no step lowers S (the damping reached its limit)
  Iterations,     // max_iterations
  Budget,         // max_simulations
  GlobalDone,     // the global search used its candidates (and no refine)
  Cancelled,
  StartFailed     // a run at the starting point failed
};

const char* StopReasonName(StopReason r);  // "converged", "small_step", ...

struct ParameterEstimate {
  std::string path;
  std::string label;
  int case_index = -1;        // -1: shared
  double estimate = kNaN;     // SI units
  double standard_error = kNaN;
  double initial = kNaN;
  double lower = kNaN;
  double upper = kNaN;
  bool at_lower = false;      // within 1e-6 of the normalised range of a bound
  bool at_upper = false;
  Scale scale = Scale::Linear;
};

struct DatasetFit {
  size_t case_index = 0;
  std::string case_name;
  std::string name;
  double weight = 1.0;
  size_t points = 0;          // points in the objective
  size_t skipped = 0;         // points left out (not covered by the simulation)
  double rmse = kNaN;         // in the dataset's y unit
  double nrmse = kNaN;        // rmse / (max - min of the measured y); NaN for a flat dataset
  double max_abs_error = kNaN;
  double objective = kNaN;    // this dataset's share of S
  // At the points used: the measured x and y, the simulated y and the
  // residual y_sim - y (dataset units) and r_k (the objective's entries).
  std::vector<double> x;
  std::vector<double> measured;
  std::vector<double> simulated;
  std::vector<double> residuals;
  std::vector<double> weighted_residuals;
  // The whole simulated series at the estimate (SimulatedSeries).
  SimulatedCounterpart series;
};

struct TracePoint {
  Solver solver = Solver::LeastSquares;
  int iteration = 0;
  double objective = kNaN;      // S at the current (best) point
  std::vector<double> values;   // per Unknown, SI units
  int simulations = 0;
};

struct Result {
  std::vector<ParameterEstimate> parameters;     // per Unknown
  std::vector<std::vector<double>> correlation;  // per Unknown x Unknown (NaN when unavailable)
  std::vector<DatasetFit> datasets;              // per Problem dataset
  std::vector<TracePoint> trace;
  double objective = kNaN;      // S at the estimate
  double sigma2 = kNaN;         // S / (N - p)
  size_t residuals = 0;         // N
  int simulations = 0;
  int failed_simulations = 0;
  int iterations = 0;           // least-squares iterations + global candidates
  double wall_time_s = 0.0;
  Solver solver = Solver::LeastSquares;
  bool refined = false;         // a global search followed by least squares
  StopReason stop = StopReason::Converged;
  std::string message;          // why it stopped, in words
  std::vector<std::string> notes;  // what to know about the result
  bool cancelled = false;
  // One deck per case with the estimates applied (JSON text in the
  // "simulation" string form; the case's policy settings are not applied).
  std::vector<std::string> calibrated_decks;

  // True when the fit has an evaluated estimate (it did not fail or get
  // cancelled at the start).
  bool Ok() const { return stop != StopReason::StartFailed && std::isfinite(objective); }
};

// Runs the calibration. Never throws for a failed run; throws
// std::invalid_argument only for invalid options.
Result Calibrate(const Problem& problem, const Options& options = {});

// The residual vector and objective at a point (per Unknown, SI units),
// without fitting: one run per case. Useful to check a problem's set-up.
// `ok` is false when a run failed (then the residuals are empty).
struct Evaluation {
  bool ok = false;
  double objective = kNaN;
  std::vector<double> residuals;
  std::string message;
};
Evaluation Evaluate(const Problem& problem, const std::vector<double>& values, const Options& options = {});

}  // namespace calibration
}  // namespace adsim
