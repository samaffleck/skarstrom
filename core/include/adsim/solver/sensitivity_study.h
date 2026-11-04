/*
 * Skarstrom - pressure and temperature swing adsorption simulator
 * Copyright (c) 2024-2026 Sam Affleck and contributors
 *
 * Licensed under the terms in the LICENSE file at the repository root.
 */
#pragma once

// The set-up of a sensitivity study, saved with the deck (UserData's
// Sensitivity_Study): its uncertain inputs, the KPIs it watches, how each run
// is made and how large the screening and the propagation are. Only the
// set-up is saved; a study's runs and results are not (they belong to the
// front end or script that ran it). The study itself (plans, running one
// sample, analysis) is in adsim/solver/sensitivity.h, which includes this
// header; this one has no dependency on the simulation, so UserData can hold
// a Study.
//
// What is saved is kept as written: an input whose parameter has since been
// removed from the deck, or an output whose component is no longer tracked,
// stays in the study and is reported by CheckStudy, never dropped silently.
// Inputs are saved by their registry path (the ID form, e.g. "step[2].time",
// survives renaming a step; an alias does not).

#include <cstdint>
#include <limits>
#include <string>
#include <string_view>
#include <utility>
#include <vector>

#include <cereal/cereal.hpp>
#include <cereal/types/string.hpp>
#include <cereal/types/utility.hpp>
#include <cereal/types/vector.hpp>

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

  // The distribution is saved by name; an unset number is saved as NaN.
  template <class Archive>
  void save(Archive& ar, std::uint32_t const /*version*/) const {
    ar(cereal::make_nvp("Path", path), cereal::make_nvp("Distribution", std::string(DistributionName(distribution))),
       cereal::make_nvp("Lower", lower), cereal::make_nvp("Upper", upper), cereal::make_nvp("Mode", mode),
       cereal::make_nvp("Mean", mean), cereal::make_nvp("SD", sd));
  }
  template <class Archive>
  void load(Archive& ar, std::uint32_t const /*version*/) {
    std::string name;
    ar(cereal::make_nvp("Path", path), cereal::make_nvp("Distribution", name), cereal::make_nvp("Lower", lower),
       cereal::make_nvp("Upper", upper), cereal::make_nvp("Mode", mode), cereal::make_nvp("Mean", mean),
       cereal::make_nvp("SD", sd));
    distribution = DistributionFromName(name);
  }
};

Input Uniform(std::string path, double lower, double upper);
// Uniform over nominal * (1 -/+ fraction), e.g. fraction 0.1 for +-10 %.
Input Relative(std::string path, double nominal, double fraction);
Input Triangular(std::string path, double lower, double mode, double upper);
Input Normal(std::string path, double mean, double sd, double lower = kNaN, double upper = kNaN);

// ---- how each sample is run ----------------------------------------------------

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

  template <class Archive>
  void serialize(Archive& ar, std::uint32_t const /*version*/) {
    ar(cereal::make_nvp("Require_CSS", require_css), cereal::make_nvp("Cycles", cycles),
       cereal::make_nvp("Settings", settings));
  }
};

// ---- the study ---------------------------------------------------------------

// A deck's sensitivity study: the screening (Morris) and the propagation
// (Latin hypercube Monte Carlo) over the same inputs. The defaults are those
// of the Python Study and the `skarstrom sensitivity` command.
struct Study {
  std::vector<Input> inputs;
  // "Component.kpi" names (sensitivity::ParseOutput); empty means the deck's
  // DefaultOutputs at the time the study runs.
  std::vector<std::string> outputs;
  RunPolicy policy;
  std::uint64_t seed = 1;

  // Screening: Morris trajectories and grid levels (even, >= 2).
  int trajectories = 10;
  int levels = 4;
  // Screening then propagation: the inputs propagated are those whose mu*
  // reaches this fraction of the largest mu* of some output.
  double threshold = 0.1;

  // Propagation: the first batch of runs, the runs added per batch until
  // every output's percentiles are known to within rel_tol of its 5-95 %
  // band, and the run limit.
  int runs = 100;
  int batch = 50;
  int max_runs = 400;
  double rel_tol = 0.1;

  // Bootstrap resamples for the intervals of both analyses.
  int bootstrap = 1000;

  // A deck without a study has no inputs.
  bool Empty() const { return inputs.empty(); }

  template <class Archive>
  void serialize(Archive& ar, std::uint32_t const /*version*/) {
    ar(cereal::make_nvp("Inputs", inputs), cereal::make_nvp("Outputs", outputs), cereal::make_nvp("Policy", policy),
       cereal::make_nvp("Seed", seed), cereal::make_nvp("Trajectories", trajectories),
       cereal::make_nvp("Levels", levels), cereal::make_nvp("Threshold", threshold), cereal::make_nvp("Runs", runs),
       cereal::make_nvp("Batch", batch), cereal::make_nvp("Max_Runs", max_runs), cereal::make_nvp("Rel_Tol", rel_tol),
       cereal::make_nvp("Bootstrap", bootstrap));
  }
};

}  // namespace sensitivity
}  // namespace adsim

CEREAL_CLASS_VERSION(adsim::sensitivity::Input, 0);
CEREAL_CLASS_VERSION(adsim::sensitivity::RunPolicy, 0);
CEREAL_CLASS_VERSION(adsim::sensitivity::Study, 0);
