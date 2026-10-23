/*
 * Skarstrom - pressure and temperature swing adsorption simulator
 * Copyright (c) 2024-2026 Sam Affleck and contributors
 *
 * Licensed under the terms in the LICENSE file at the repository root.
 */
#pragma once

// Sensitivity analysis and uncertainty propagation for a deck's KPIs. The
// implementation is in core/src/solver/sensitivity.cpp.
//
// Two studies, both over a list of uncertain inputs (registry paths, each with
// a probability distribution):
//
//   1. Screening (Morris elementary effects). Which inputs move the KPIs?
//      r trajectories of k+1 runs each; along a trajectory the inputs change
//      one at a time, in a random order, by half of their probability range.
//      Per input and KPI: mu* (mean absolute effect: the ranking), mu (mean
//      signed effect) and sigma (spread of the effect: interactions or
//      non-linearity), with a bootstrap interval on mu*.
//
//   2. Propagation (Monte Carlo with Latin hypercube sampling). How uncertain
//      are the KPIs? All inputs vary together; per KPI: mean, standard
//      deviation and the 5th / 50th / 95th percentiles (the 90 % band is the
//      error bar) with bootstrap intervals and a convergence check, and each
//      input's share of the spread from a rank regression (with its R^2).
//      Runs come in batches (independent Latin hypercubes), so a caller adds
//      batches until AnalyseMonteCarlo reports convergence.
//
// The engine plans the runs and analyses their results; the caller executes
// them (one EvaluateSample per run, from the deck's JSON), so runs can go to
// several processes or web workers. Plans are reproducible across platforms:
// the generator is std::mt19937_64 (fully specified by the standard) and every
// draw from it is done here, not by the implementation-defined std::
// distributions.
//
// Inputs are sampled in probability space: u in (0, 1) maps to a value with
// the distribution's quantile function, so Uniform, Triangular and Normal
// inputs are handled alike. Elementary effects are therefore per unit of u:
// mu* is, in KPI units, the typical change of the KPI when the input moves
// across its whole range. Morris levels are the centres of `levels`
// equiprobable cells (u = (j + 1/2) / levels), so an unbounded Normal input
// never maps to an infinite value.
//
// Reference: A. Saltelli et al., Global Sensitivity Analysis: The Primer,
// Wiley (2008), doi:10.1002/9780470725184 (chapter 3 for the screening).

#include <cstdint>
#include <functional>
#include <limits>
#include <string>
#include <string_view>
#include <utility>
#include <vector>

#include "adsim/simulation.h"

namespace adsim {
namespace sensitivity {

inline constexpr double kNaN = std::numeric_limits<double>::quiet_NaN();

// ---- inputs ------------------------------------------------------------------

enum class Distribution {
  Uniform,     // lower .. upper
  Triangular,  // lower .. upper, peak at mode
  Normal       // mean, sd; truncated to [lower, upper] where those are set
};

const char* DistributionName(Distribution d);  // "uniform", "triangular", "normal"
// Throws std::invalid_argument for an unknown name.
Distribution DistributionFromName(std::string_view name);

// One uncertain input. `path` is a registry path or alias (see
// parameters/parameter_registry.h); the parameter must be tunable. Values are
// in the parameter's SI unit.
struct Input {
  std::string path;
  Distribution distribution = Distribution::Uniform;
  double lower = kNaN;  // Uniform, Triangular: required. Normal: optional truncation
  double upper = kNaN;
  double mode = kNaN;   // Triangular
  double mean = kNaN;   // Normal
  double sd = kNaN;     // Normal

  // The value at probability u in (0, 1) (the inverse CDF).
  double Quantile(double u) const;
  // Throws std::invalid_argument naming the path when the description is
  // incomplete or inconsistent (lower >= upper, mode outside, sd <= 0, ...).
  void Validate() const;
};

Input Uniform(std::string path, double lower, double upper);
// Uniform over nominal * (1 -/+ fraction), e.g. fraction 0.1 for +-10 %.
Input Relative(std::string path, double nominal, double fraction);
Input Triangular(std::string path, double lower, double mode, double upper);
Input Normal(std::string path, double mean, double sd, double lower = kNaN, double upper = kNaN);

// ---- plans -------------------------------------------------------------------

enum class Method { Morris, MonteCarlo };

const char* MethodName(Method m);  // "morris", "monte_carlo"

// The runs of a study. Row r of `u` / `values` is run r, one column per input.
struct Plan {
  Method method = Method::Morris;
  std::vector<Input> inputs;
  std::uint64_t seed = 1;
  // Morris: trajectories and grid levels. Run t*(k+1) starts trajectory t;
  // changed[r] is the input that changed from run r-1 to run r, or -1 for
  // the first run of a trajectory.
  int trajectories = 0;
  int levels = 0;
  std::vector<int> changed;
  // Monte Carlo: the batch index of this plan (each batch is an independent
  // Latin hypercube drawn from the seed and the batch index).
  int batch = 0;

  std::vector<std::vector<double>> u;       // probability coordinates in (0, 1)
  std::vector<std::vector<double>> values;  // input values, Quantile(u)

  size_t Runs() const { return values.size(); }
};

// Morris trajectories: trajectories * (k + 1) runs. `levels` must be even and
// >= 2 (4 is the usual choice); every step moves one input by half of the
// levels (|du| = 1/2). Throws std::invalid_argument for bad settings or inputs.
Plan MorrisPlan(const std::vector<Input>& inputs, int trajectories = 10, int levels = 4,
                std::uint64_t seed = 1);

// One Latin hypercube batch of `runs` runs: each input's range is cut into
// `runs` equiprobable strata and every stratum is sampled once.
Plan MonteCarloPlan(const std::vector<Input>& inputs, int runs, std::uint64_t seed = 1,
                    int batch = 0);

// ---- running one sample ------------------------------------------------------

// The KPI fields a study can watch (KPIRecord members): purity, recovery,
// productivity_mass, productivity_volume, energy_consumption, production_rate.
const std::vector<std::string>& KPINames();

// One watched value: a tracked component's KPI at the last logged cycle.
struct Output {
  std::string component;  // CoolProp name, as in the deck's KPI configuration
  std::string kpi;        // one of KPINames()

  std::string Name() const { return component + "." + kpi; }
};

// Parses "Component.kpi" (the component may itself contain dots: the KPI is
// the text after the last one). Throws std::invalid_argument.
Output ParseOutput(std::string_view name);

// Purity, recovery, mass productivity and specific energy of every component
// the deck tracks.
std::vector<Output> DefaultOutputs(const UserData& data);

// How each sample is run.
struct RunPolicy {
  // A run must reach cyclic steady state to count; one that stops at the
  // cycle limit is treated as failed. Turn off only for screening at a fixed
  // number of cycles, where every run is equally (un)converged.
  bool require_css = true;
  // Overrides the deck's number of cycles (0: keep the deck's).
  int cycles = 0;
  // Applied to every run before the sample values, as registry text sets
  // (path, text): typically a tighter CSS tolerance for the study.
  std::vector<std::pair<std::string, std::string>> settings;
};

struct SampleOutcome {
  RunStatus status = RunStatus::NotRun;
  std::string message;
  // The run counts: it finished and (with require_css) reached CSS.
  bool accepted = false;
  // One value per Output; NaN when the run was not accepted or the KPI is
  // undefined for it (e.g. no product drawn).
  std::vector<double> outputs;
  double wall_time_s = 0.0;
};

// Checks a study before any run: loads `deck_json`, applies policy.settings
// and throws std::invalid_argument with the reason when an input is invalid
// or not a tunable parameter of the deck (or listed twice), an output's KPI
// is unknown or its component not tracked, or the policy requires CSS but the
// deck has no CSS check.
void CheckStudy(const std::string& deck_json, const std::vector<Input>& inputs,
                const std::vector<Output>& outputs, const RunPolicy& policy = {});

// Runs one sample on a fresh copy of the deck: loads `deck_json`, applies
// policy.settings, sets inputs[i] to values[i] and runs. `values` may be
// empty for the nominal (unperturbed) run. Never throws for a failed run:
// the failure is in the outcome. Throws std::invalid_argument for a deck,
// setting or input that cannot be applied.
SampleOutcome EvaluateSample(const std::string& deck_json, const std::vector<Input>& inputs,
                             const std::vector<double>& values, const std::vector<Output>& outputs,
                             const RunPolicy& policy = {}, const CancelPoll& cancel = {});

// ---- analysis ----------------------------------------------------------------

struct Interval {
  double low = kNaN;
  double high = kNaN;
};

struct MorrisEffect {
  std::string input;          // the input's path
  double mu = kNaN;           // mean elementary effect
  double mu_star = kNaN;      // mean absolute elementary effect
  double sigma = kNaN;        // sample standard deviation of the effects
  Interval mu_star_interval;  // 95 % bootstrap interval of mu*
  int effects = 0;            // elementary effects used (both runs finite)
};

struct MorrisOutput {
  std::string output;
  std::vector<MorrisEffect> effects;  // one per input, in input order
  int missing = 0;                    // runs whose value is NaN
};

struct MorrisResult {
  std::vector<MorrisOutput> outputs;
  int runs = 0;
};

// `y` is [run][output] in plan order (NaN for a failed run or undefined KPI);
// `output_names` labels the columns. Elementary effects that touch a NaN are
// skipped and counted. Throws std::invalid_argument when the shapes disagree
// or the plan is not a Morris plan.
MorrisResult AnalyseMorris(const Plan& plan, const std::vector<std::string>& output_names,
                           const std::vector<std::vector<double>>& y, int bootstrap = 1000);

struct InputShare {
  std::string input;
  double spearman = kNaN;  // rank correlation of the input with the output
  double srrc = kNaN;      // standardised rank regression coefficient
  double share = kNaN;     // srrc^2: fraction of the output's (rank) variance
};

struct Spread {
  std::string output;
  int valid = 0;    // runs with a finite value
  int missing = 0;  // runs whose value is NaN
  double mean = kNaN, sd = kNaN, min = kNaN, max = kNaN;
  double p05 = kNaN, p50 = kNaN, p95 = kNaN;  // percentiles (linear interpolation)
  Interval p05_interval, p50_interval, p95_interval;  // 95 % bootstrap intervals
  // Every percentile interval is narrower than rel_tol * (p95 - p05) on each
  // side (or the output does not vary at all).
  bool converged = false;
  // R^2 of the rank regression of the output on all inputs. Shares are only
  // a full account of the spread when it is close to 1; a low R^2 means the
  // inputs act through interactions or non-monotonic effects.
  double r2 = kNaN;
  std::vector<InputShare> shares;  // one per input, in input order
};

struct MonteCarloResult {
  std::vector<Spread> outputs;
  int runs = 0;
  bool converged = false;  // every output converged
};

// `values` is [run][input] (all batches stacked), `y` is [run][output].
// Throws std::invalid_argument when the shapes disagree.
MonteCarloResult AnalyseMonteCarlo(const std::vector<Input>& inputs,
                                   const std::vector<std::vector<double>>& values,
                                   const std::vector<std::string>& output_names,
                                   const std::vector<std::vector<double>>& y,
                                   std::uint64_t seed = 1, int bootstrap = 1000,
                                   double rel_tol = 0.1);

// Percentile of finite values by linear interpolation between order
// statistics (numpy's default): q in [0, 1]. NaN when there are none.
double Percentile(std::vector<double> x, double q);

}  // namespace sensitivity
}  // namespace adsim
