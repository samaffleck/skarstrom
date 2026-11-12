/*
 * Skarstrom - pressure and temperature swing adsorption simulator
 * Copyright (c) 2024-2026 Sam Affleck and contributors
 *
 * Licensed under the terms in the LICENSE file at the repository root.
 */
#include "adsim/solver/calibration.h"

#include <algorithm>
#include <atomic>
#include <chrono>
#include <cmath>
#include <fstream>
#include <memory>
#include <sstream>
#include <stdexcept>
#include <system_error>
#include <thread>
#include <utility>

#include <Eigen/Dense>

#include "adsim/solver/optimisation.h"  // GlobalOptimizer, kOptimisationPenaltyCost

namespace adsim {
namespace calibration {

// ---- names -------------------------------------------------------------------

const char* ScaleName(Scale s) { return s == Scale::Log ? "log" : "linear"; }
const char* SharingName(Sharing s) { return s == Sharing::PerCase ? "per_case" : "shared"; }
const char* SolverName(Solver s) { return s == Solver::Global ? "global" : "least_squares"; }
const char* DifferenceName(Difference d) { return d == Difference::Central ? "central" : "forward"; }

const char* StopReasonName(StopReason r) {
  switch (r) {
    case StopReason::Converged: return "converged";
    case StopReason::SmallStep: return "small_step";
    case StopReason::NoImprovement: return "no_improvement";
    case StopReason::Iterations: return "iterations";
    case StopReason::Budget: return "budget";
    case StopReason::GlobalDone: return "global_done";
    case StopReason::Cancelled: return "cancelled";
    case StopReason::StartFailed: return "start_failed";
  }
  return "unknown";
}

// ---- building the problem ----------------------------------------------------

DatasetRef Dataset(std::string name, double weight) {
  DatasetRef r;
  r.name = std::move(name);
  r.weight = weight;
  return r;
}

DatasetRef Dataset(int index, double weight) {
  DatasetRef r;
  r.index = index;
  r.weight = weight;
  return r;
}

Case Case::FromFile(const std::string& path, std::vector<DatasetRef> datasets) {
  std::ifstream in(path, std::ios::binary);
  if (!in) throw std::runtime_error("cannot read the deck " + path);
  std::ostringstream text;
  text << in.rdbuf();
  Case c;
  c.name = path;
  const auto slash = path.find_last_of("/\\");
  if (slash != std::string::npos) c.name = path.substr(slash + 1);
  c.deck_json = text.str();
  c.datasets = std::move(datasets);
  return c;
}

Case Case::FromData(UserData& data, std::vector<DatasetRef> datasets) {
  Case c;
  SaveToString(data, c.deck_json);
  c.datasets = std::move(datasets);
  return c;
}

Parameter Shared(std::string path, double lower, double upper, Scale scale, double initial) {
  Parameter p;
  p.path = std::move(path);
  p.lower = lower;
  p.upper = upper;
  p.scale = scale;
  p.initial = initial;
  p.sharing = Sharing::Shared;
  return p;
}

Parameter PerCase(std::string path, double lower, double upper, Scale scale, double initial) {
  Parameter p = Shared(std::move(path), lower, upper, scale, initial);
  p.sharing = Sharing::PerCase;
  return p;
}

namespace {

constexpr double kBoundTolerance = 1e-6;  // normalised distance that counts as "at a bound"

std::string Quote(const std::string& s) { return "'" + s + "'"; }

void LoadCase(Simulation& sim, const Case& c) {
  try {
    sim.LoadJson(c.deck_json);
  } catch (const std::exception& e) {
    throw std::invalid_argument("case " + Quote(c.name) + ": the deck does not load: " + e.what());
  }
  for (const auto& [path, text] : c.policy.settings) {
    try {
      sim.SetParameterText(path, text);
    } catch (const std::exception& e) {
      throw std::invalid_argument("case " + Quote(c.name) + ": setting " + Quote(path) + ": " + e.what());
    }
  }
}

}  // namespace

Problem::Problem(std::vector<Case> cases, std::vector<Parameter> parameters)
    : cases_(std::move(cases)), parameters_(std::move(parameters)) {
  if (cases_.empty()) throw std::invalid_argument("a calibration needs at least one case");
  if (parameters_.empty()) throw std::invalid_argument("a calibration needs at least one parameter");
  for (size_t c = 0; c < cases_.size(); ++c)
    if (cases_[c].name.empty()) cases_[c].name = "case " + std::to_string(c + 1);
  const size_t n_cases = cases_.size();

  for (const auto& p : parameters_) {
    const std::string who = "parameter " + Quote(p.path);
    auto check_bounds = [&](double lo, double hi, const std::string& where) {
      if (!std::isfinite(lo) || !std::isfinite(hi) || !(lo < hi))
        throw std::invalid_argument(who + where + " needs finite bounds with lower < upper");
      if (p.scale == Scale::Log && !(lo > 0.0))
        throw std::invalid_argument(who + where + " has a log scale: its lower bound must be > 0");
    };
    check_bounds(p.lower, p.upper, "");
    if (p.sharing == Sharing::Shared) {
      if (!p.case_bounds.empty() || !p.case_initial.empty())
        throw std::invalid_argument(who + " is shared: per-case bounds or starts apply only to a per-case parameter");
    } else {
      if (!p.case_bounds.empty() && p.case_bounds.size() != n_cases)
        throw std::invalid_argument(who + ": " + std::to_string(p.case_bounds.size()) + " per-case bounds for " +
                                    std::to_string(n_cases) + " cases");
      if (!p.case_initial.empty() && p.case_initial.size() != n_cases)
        throw std::invalid_argument(who + ": " + std::to_string(p.case_initial.size()) + " per-case starts for " +
                                    std::to_string(n_cases) + " cases");
    }
  }

  // Check every case's deck and read the deck values (the default starts).
  std::vector<std::vector<double>> deck_value(parameters_.size(), std::vector<double>(n_cases, kNaN));
  for (size_t c = 0; c < n_cases; ++c) {
    const Case& cs = cases_[c];
    const std::string who = "case " + Quote(cs.name);
    if (cs.policy.cycles < 0) throw std::invalid_argument(who + ": the cycle cap must be >= 0");
    Simulation sim;
    LoadCase(sim, cs);
    ParameterRegistry& registry = sim.Parameters();
    std::vector<std::string> seen;
    for (size_t p = 0; p < parameters_.size(); ++p) {
      ::adsim::Parameter* rp = registry.Find(parameters_[p].path);
      if (!rp) throw std::invalid_argument(who + ": no parameter " + Quote(parameters_[p].path) + " in the deck");
      if (!rp->Info().tunable)
        throw std::invalid_argument(who + ": " + Quote(parameters_[p].path) + " is not a tunable parameter");
      const std::string& resolved = rp->Info().path;
      if (std::find(seen.begin(), seen.end(), resolved) != seen.end())
        throw std::invalid_argument("parameter " + Quote(parameters_[p].path) + " is listed twice (" + resolved + ")");
      seen.push_back(resolved);
      deck_value[p][c] = sim.GetParameter(parameters_[p].path);
    }

    // The datasets.
    const std::vector<XYData>& all = sim.Data().exp_data_.GetAllData();
    std::vector<std::pair<int, double>> chosen;  // deck index, weight override
    if (cs.datasets.empty()) {
      for (size_t i = 0; i < all.size(); ++i)
        if (all[i].observable.Bound()) chosen.emplace_back(static_cast<int>(i), kNaN);
      if (chosen.empty()) throw std::invalid_argument(who + ": the deck has no dataset bound to an observable");
    } else {
      for (const auto& ref : cs.datasets) {
        int index = -1;
        if (!ref.name.empty()) {
          for (size_t i = 0; i < all.size(); ++i)
            if (all[i].name == ref.name) {
              index = static_cast<int>(i);
              break;
            }
          if (index < 0) throw std::invalid_argument(who + ": no dataset " + Quote(ref.name) + " in the deck");
        } else {
          if (ref.index < 0 || static_cast<size_t>(ref.index) >= all.size())
            throw std::invalid_argument(who + ": no dataset at index " + std::to_string(ref.index) + " (the deck has " +
                                        std::to_string(all.size()) + ")");
          index = ref.index;
        }
        for (const auto& ch : chosen)
          if (ch.first == index)
            throw std::invalid_argument(who + ": dataset " + Quote(all[static_cast<size_t>(index)].name) +
                                        " is listed twice");
        chosen.emplace_back(index, ref.weight);
      }
    }
    for (const auto& [index, weight] : chosen) {
      ProblemDataset d;
      d.case_index = c;
      d.deck_index = index;
      d.data = all[static_cast<size_t>(index)];
      const std::string dwho = who + ", dataset " + Quote(d.data.name);
      if (!d.data.observable.Bound()) throw std::invalid_argument(dwho + " is not bound to an observable");
      if (d.data.x.empty() || d.data.x.size() != d.data.y.size())
        throw std::invalid_argument(dwho + " needs points with as many x as y values");
      if (std::isfinite(weight)) d.data.observable.weight = weight;
      if (!(d.data.observable.weight > 0.0) || !std::isfinite(d.data.observable.weight))
        throw std::invalid_argument(dwho + ": the weight must be > 0");
      if (d.data.observable.kind == ObservableKind::Kpi) {
        bool tracked = false;
        for (const auto& k : sim.Data().kpi_config.per_component)
          if (k.tracked && k.product_component == d.data.observable.component) tracked = true;
        if (!tracked)
          throw std::invalid_argument(dwho + " measures a KPI of " + Quote(d.data.observable.component) +
                                      ", which the deck's KPI configuration does not track");
      }
      const auto& u = d.data.observable.uncertainty;
      if (!u.empty()) {
        if (u.size() != d.data.y.size())
          throw std::invalid_argument(dwho + ": the uncertainty needs one value per point");
        for (double s : u)
          if (!(s > 0.0) || !std::isfinite(s)) throw std::invalid_argument(dwho + ": every uncertainty must be > 0");
      }
      double lo = std::numeric_limits<double>::infinity(), hi = -lo;
      for (double y : d.data.y)
        if (std::isfinite(y)) {
          lo = std::min(lo, y);
          hi = std::max(hi, y);
        }
      if (!(hi >= lo)) throw std::invalid_argument(dwho + " has no finite y value");
      d.scale = hi > lo ? hi - lo : 1.0;
      datasets_.push_back(std::move(d));
    }
  }

  // The unknowns.
  unknown_of_.assign(parameters_.size(), std::vector<size_t>(n_cases, 0));
  for (size_t p = 0; p < parameters_.size(); ++p) {
    const Parameter& par = parameters_[p];
    auto make = [&](int c) {
      Unknown u;
      u.parameter = p;
      u.case_index = c;
      u.path = par.path;
      u.scale = par.scale;
      u.lower = par.lower;
      u.upper = par.upper;
      u.initial = par.initial;
      u.label = par.path;
      if (c >= 0) {
        const auto cu = static_cast<size_t>(c);
        u.label += " [" + cases_[cu].name + "]";
        if (!par.case_bounds.empty()) {
          if (std::isfinite(par.case_bounds[cu].lower)) u.lower = par.case_bounds[cu].lower;
          if (std::isfinite(par.case_bounds[cu].upper)) u.upper = par.case_bounds[cu].upper;
          if (!std::isfinite(u.lower) || !std::isfinite(u.upper) || !(u.lower < u.upper))
            throw std::invalid_argument("parameter " + Quote(par.path) + " in case " + Quote(cases_[cu].name) +
                                        " needs finite bounds with lower < upper");
          if (par.scale == Scale::Log && !(u.lower > 0.0))
            throw std::invalid_argument("parameter " + Quote(par.path) + " in case " + Quote(cases_[cu].name) +
                                        " has a log scale: its lower bound must be > 0");
        }
        if (!par.case_initial.empty() && std::isfinite(par.case_initial[cu])) u.initial = par.case_initial[cu];
      }
      if (std::isfinite(u.initial)) {
        if (u.initial < u.lower || u.initial > u.upper)
          throw std::invalid_argument("the start of " + u.label + " (" + std::to_string(u.initial) +
                                      ") is outside its bounds");
      } else {
        const double v = deck_value[p][c < 0 ? 0 : static_cast<size_t>(c)];
        u.initial = std::clamp(v, u.lower, u.upper);
        if (!(v >= u.lower && v <= u.upper)) {
          if (!std::isfinite(v)) u.initial = 0.5 * (u.lower + u.upper);
          notes_.push_back("the deck's value of " + u.label + " (" + std::to_string(v) +
                           ") is outside the bounds; the fit starts from " + std::to_string(u.initial));
        }
      }
      return u;
    };
    if (par.sharing == Sharing::Shared) {
      unknowns_.push_back(make(-1));
      for (size_t c = 0; c < n_cases; ++c) unknown_of_[p][c] = unknowns_.size() - 1;
    } else {
      for (size_t c = 0; c < n_cases; ++c) {
        unknowns_.push_back(make(static_cast<int>(c)));
        unknown_of_[p][c] = unknowns_.size() - 1;
      }
    }
  }
}

size_t Problem::UnknownIndex(size_t p, size_t c) const { return unknown_of_.at(p).at(c); }

// ---- the engine --------------------------------------------------------------

namespace {

using Vector = Eigen::VectorXd;
using Matrix = Eigen::MatrixXd;

// The scaled value s (x, or ln x) and its normalised form z in [0, 1].
double ToScaled(const Unknown& u, double x) { return u.scale == Scale::Log ? std::log(x) : x; }
double ScaledLower(const Unknown& u) { return ToScaled(u, u.lower); }
double ScaledRange(const Unknown& u) { return ToScaled(u, u.upper) - ToScaled(u, u.lower); }
double ToNormalised(const Unknown& u, double x) {
  return std::clamp((ToScaled(u, x) - ScaledLower(u)) / ScaledRange(u), 0.0, 1.0);
}
double FromNormalised(const Unknown& u, double z) {
  const double s = ScaledLower(u) + z * ScaledRange(u);
  const double x = u.scale == Scale::Log ? std::exp(s) : s;
  return std::clamp(x, u.lower, u.upper);
}

// One run of one case.
struct Job {
  size_t case_index = 0;
  std::vector<double> z;           // the whole point (normalised)
  const CycleState* start = nullptr;
  double distance = kNaN;
};

struct RunOutcome {
  bool ok = false;
  bool cancelled = false;
  std::string message;
  std::vector<SimulatedCounterpart> series;  // per dataset of the case (problem order)
  std::unique_ptr<CycleState> final_state;   // a CSS state to cache (warm starts only)
};

// The residuals of a point: every case run.
struct PointEval {
  bool ok = false;
  bool cancelled = false;
  std::string message;
  double objective = std::numeric_limits<double>::infinity();
  Vector r;
  std::vector<std::shared_ptr<RunOutcome>> cases;
  std::vector<double> z;
};

class Engine {
 public:
  Engine(const Problem& problem, const Options& options, Result& result)
      : problem_(problem), options_(options), result_(result) {
    const size_t n_cases = problem.Cases().size();
    if (options.warm_start)
      for (size_t c = 0; c < n_cases; ++c) caches_.emplace_back(std::make_unique<CycleStateCache>(options.cache_capacity));
    case_datasets_.resize(n_cases);
    for (size_t d = 0; d < problem.Datasets().size(); ++d) case_datasets_[problem.Datasets()[d].case_index].push_back(d);
    threads_ = options.threads > 0 ? options.threads
                                   : std::max(1, static_cast<int>(std::thread::hardware_concurrency()) - 1);
    // Which unknowns move each case.
    affects_.assign(problem.Unknowns().size(), std::vector<size_t>{});
    for (size_t j = 0; j < problem.Unknowns().size(); ++j) {
      const Unknown& u = problem.Unknowns()[j];
      if (u.case_index < 0)
        for (size_t c = 0; c < n_cases; ++c) affects_[j].push_back(c);
      else
        affects_[j].push_back(static_cast<size_t>(u.case_index));
    }
  }

  size_t NUnknowns() const { return problem_.Unknowns().size(); }
  size_t NCases() const { return problem_.Cases().size(); }
  bool Cancelled() const { return cancelled_; }
  int Simulations() const { return result_.simulations; }
  bool HasActiveSets() const { return have_active_; }

  bool PollCancel() {
    if (!cancelled_ && options_.cancel && options_.cancel()) cancelled_ = true;
    return cancelled_;
  }

  std::vector<double> Values(const std::vector<double>& z) const {
    std::vector<double> x(z.size());
    for (size_t j = 0; j < z.size(); ++j) x[j] = FromNormalised(problem_.Unknowns()[j], z[j]);
    return x;
  }

  // The case-local key of a point (its unknowns for that case).
  std::vector<double> Key(size_t c, const std::vector<double>& z) const {
    std::vector<double> key(problem_.Parameters().size());
    for (size_t p = 0; p < key.size(); ++p) key[p] = z[problem_.UnknownIndex(p, c)];
    return key;
  }

  // Runs a batch of jobs: starts are chosen before, states cached after, in
  // job order.
  std::vector<std::shared_ptr<RunOutcome>> RunBatch(std::vector<Job> jobs) {
    for (auto& job : jobs) {
      if (!options_.warm_start) continue;
      CycleStateCache& cache = *caches_[job.case_index];
      if (!cache.Empty()) job.start = cache.Nearest(Key(job.case_index, job.z), &job.distance);
    }
    std::vector<std::shared_ptr<RunOutcome>> out(jobs.size());
    std::atomic<size_t> next{0};
    auto work = [&]() {
      for (size_t i = next++; i < jobs.size(); i = next++) out[i] = RunOne(jobs[i]);
    };
    const size_t n_threads = std::min(static_cast<size_t>(std::max(1, threads_)), jobs.size());
    std::vector<std::thread> pool;
    if (n_threads > 1) {
      try {
        for (size_t t = 1; t < n_threads; ++t) pool.emplace_back(work);
      } catch (const std::system_error&) {
        // No threads here (e.g. WebAssembly without them): the rest run below.
      }
    }
    work();
    for (auto& t : pool) t.join();

    for (size_t i = 0; i < jobs.size(); ++i) {
      ++result_.simulations;
      if (!out[i]->ok) ++result_.failed_simulations;
      if (out[i]->cancelled) cancelled_ = true;
      if (options_.warm_start && out[i]->final_state) {
        try {
          caches_[jobs[i].case_index]->Add(Key(jobs[i].case_index, jobs[i].z), std::move(*out[i]->final_state));
        } catch (const std::exception&) {
          // Not cached: the next run starts from an older state or cold.
        }
        out[i]->final_state.reset();
      }
    }
    PollCancel();
    return out;
  }

  // The residual entries of one case from its run.
  bool CaseResiduals(size_t c, const RunOutcome& run, std::vector<double>& r, std::string& message) const {
    r.clear();
    const auto& ids = case_datasets_[c];
    for (size_t k = 0; k < ids.size(); ++k) {
      const ProblemDataset& d = problem_.Datasets()[ids[k]];
      const SimulatedCounterpart& s = run.series[k];
      const auto& active = active_[ids[k]];
      for (size_t i : active) {
        const double y_sim = SimulatedAt(d, s, d.data.x[i]);
        if (!std::isfinite(y_sim)) {
          message = "case " + Quote(problem_.Cases()[c].name) + ", dataset " + Quote(d.data.name) +
                    ": the simulated value at x = " + std::to_string(d.data.x[i]) + " is not finite";
          return false;
        }
        r.push_back(Weighted(ids[k], i, y_sim));
      }
    }
    return true;
  }

  // Evaluates every case at z (one batch).
  PointEval EvaluatePoint(const std::vector<double>& z) {
    std::vector<Job> jobs;
    for (size_t c = 0; c < NCases(); ++c) jobs.push_back(Job{c, z, nullptr, kNaN});
    PointEval e;
    e.z = z;
    e.cases = RunBatch(std::move(jobs));
    Assemble(e);
    return e;
  }

  void Assemble(PointEval& e) {
    for (size_t c = 0; c < NCases(); ++c) {
      const RunOutcome& run = *e.cases[c];
      if (run.cancelled) e.cancelled = true;
      if (!run.ok) {
        e.ok = false;
        e.message = "case " + Quote(problem_.Cases()[c].name) + ": " + run.message;
        return;
      }
    }
    if (!have_active_) SetActiveSets(e);
    if (!have_active_) {
      e.ok = false;
      if (e.message.empty()) e.message = "no measured point is covered by the simulation";
      return;
    }
    std::vector<double> all, part;
    for (size_t c = 0; c < NCases(); ++c) {
      if (!CaseResiduals(c, *e.cases[c], part, e.message)) {
        e.ok = false;
        return;
      }
      all.insert(all.end(), part.begin(), part.end());
    }
    e.r = Eigen::Map<Vector>(all.data(), static_cast<Eigen::Index>(all.size()));
    e.objective = e.r.squaredNorm();
    e.ok = std::isfinite(e.objective);
    if (!e.ok) e.message = "the objective is not finite";
  }

  // The offsets of each case's residuals in the full vector.
  std::vector<size_t> CaseOffsets() const {
    std::vector<size_t> off(NCases() + 1, 0);
    for (size_t c = 0; c < NCases(); ++c) {
      size_t n = 0;
      for (size_t d : case_datasets_[c]) n += active_[d].size();
      off[c + 1] = off[c] + n;
    }
    return off;
  }

  // The number of runs a Jacobian takes.
  int JacobianCost() const {
    int n = 0;
    for (const auto& a : affects_) n += static_cast<int>(a.size());
    return options_.difference == Difference::Central ? 2 * n : n;
  }

  // The finite-difference Jacobian at `centre` (normalised variables). Returns
  // false when cancelled.
  bool Jacobian(const PointEval& centre, Matrix& J) {
    const size_t p = NUnknowns();
    const std::vector<size_t> off = CaseOffsets();
    J = Matrix::Zero(centre.r.size(), static_cast<Eigen::Index>(p));
    const double h = options_.fd_step;
    const bool central = options_.difference == Difference::Central;

    struct Column {
      double hp = 0.0, hm = 0.0;  // steps used (hm 0: one-sided)
    };
    std::vector<Column> cols(p);
    std::vector<Job> jobs;
    struct Slot {
      size_t j, c;
      bool plus;
    };
    std::vector<Slot> slots;
    auto add = [&](size_t j, double step, bool plus) {
      std::vector<double> z = centre.z;
      z[j] += step;
      for (size_t c : affects_[j]) {
        jobs.push_back(Job{c, z, nullptr, kNaN});
        slots.push_back(Slot{j, c, plus});
      }
    };
    for (size_t j = 0; j < p; ++j) {
      const double zj = centre.z[j];
      if (central && zj - h >= 0.0 && zj + h <= 1.0) {
        cols[j].hp = h;
        cols[j].hm = -h;
        add(j, h, true);
        add(j, -h, false);
      } else {
        cols[j].hp = zj + h <= 1.0 ? h : -h;  // inwards at the upper bound
        add(j, cols[j].hp, true);
      }
    }
    auto outcomes = RunBatch(jobs);
    if (cancelled_) return false;

    // Per column and case: the residual segments on each side.
    std::vector<std::vector<std::vector<double>>> plus(p, std::vector<std::vector<double>>(NCases()));
    std::vector<std::vector<std::vector<double>>> minus = plus;
    std::vector<std::vector<char>> plus_ok(p, std::vector<char>(NCases(), 0)), minus_ok = plus_ok;
    std::string msg;
    for (size_t i = 0; i < slots.size(); ++i) {
      const Slot& s = slots[i];
      std::vector<double> seg;
      const bool ok = outcomes[i]->ok && CaseResiduals(s.c, *outcomes[i], seg, msg);
      (s.plus ? plus : minus)[s.j][s.c] = std::move(seg);
      (s.plus ? plus_ok : minus_ok)[s.j][s.c] = ok ? 1 : 0;
    }

    // Retry failed one-sided columns on the other side; a failed central side
    // falls back to the other one.
    std::vector<Job> retry;
    std::vector<Slot> retry_slots;
    for (size_t j = 0; j < p; ++j) {
      if (cols[j].hm != 0.0) continue;
      for (size_t c : affects_[j]) {
        if (plus_ok[j][c]) continue;
        const double other = -cols[j].hp;
        const double zj = centre.z[j] + other;
        if (zj < 0.0 || zj > 1.0) continue;
        std::vector<double> z = centre.z;
        z[j] = zj;
        retry.push_back(Job{c, z, nullptr, kNaN});
        retry_slots.push_back(Slot{j, c, false});
      }
    }
    if (!retry.empty()) {
      auto again = RunBatch(retry);
      if (cancelled_) return false;
      for (size_t i = 0; i < retry_slots.size(); ++i) {
        const Slot& s = retry_slots[i];
        std::vector<double> seg;
        const bool ok = again[i]->ok && CaseResiduals(s.c, *again[i], seg, msg);
        minus[s.j][s.c] = std::move(seg);
        minus_ok[s.j][s.c] = ok ? 1 : 0;
      }
    }

    for (size_t j = 0; j < p; ++j) {
      for (size_t c : affects_[j]) {
        const auto o = static_cast<Eigen::Index>(off[c]);
        const auto n = static_cast<Eigen::Index>(off[c + 1] - off[c]);
        const Vector base = centre.r.segment(o, n);
        auto seg = [&](const std::vector<double>& v) {
          return Eigen::Map<const Vector>(v.data(), static_cast<Eigen::Index>(v.size()));
        };
        const bool p_ok = plus_ok[j][c] != 0, m_ok = minus_ok[j][c] != 0;
        if (cols[j].hm != 0.0 && p_ok && m_ok) {
          J.col(static_cast<Eigen::Index>(j)).segment(o, n) = (seg(plus[j][c]) - seg(minus[j][c])) / (2.0 * h);
        } else if (p_ok) {
          J.col(static_cast<Eigen::Index>(j)).segment(o, n) = (seg(plus[j][c]) - base) / cols[j].hp;
        } else if (m_ok) {
          const double hm = cols[j].hm != 0.0 ? cols[j].hm : -cols[j].hp;
          J.col(static_cast<Eigen::Index>(j)).segment(o, n) = (seg(minus[j][c]) - base) / hm;
        } else {
          Note("the derivative of the residuals of case " + Quote(problem_.Cases()[c].name) + " with respect to " +
               problem_.Unknowns()[j].label + " could not be made (the perturbed runs failed); taken as zero");
        }
      }
    }
    return true;
  }

  void Note(const std::string& text) {
    if (std::find(result_.notes.begin(), result_.notes.end(), text) == result_.notes.end())
      result_.notes.push_back(text);
  }

  void Record(Solver solver, int iteration, double objective, const std::vector<double>& z) {
    TracePoint t;
    t.solver = solver;
    t.iteration = iteration;
    t.objective = objective;
    t.values = Values(z);
    t.simulations = result_.simulations;
    result_.trace.push_back(t);
    if (options_.progress) {
      Progress pr;
      pr.solver = solver;
      pr.iteration = iteration;
      pr.objective = objective;
      pr.values = t.values;
      pr.simulations = t.simulations;
      options_.progress(pr);
    }
  }

  // Per dataset, the indices of the points in the objective.
  const std::vector<std::vector<size_t>>& ActiveSets() const { return active_; }
  const std::vector<std::vector<size_t>>& CaseDatasets() const { return case_datasets_; }

  double SimulatedAt(const ProblemDataset& d, const SimulatedCounterpart& s, double x) const {
    if (s.x.empty()) return kNaN;
    if (d.data.observable.kind == ObservableKind::Kpi) {
      const double cyc = std::round(x);
      size_t last = 0;
      for (size_t j = 0; j < s.x.size(); ++j) {
        if (s.x[j] == cyc) return s.y[j];
        if (s.x[j] > s.x[last]) last = j;
      }
      return cyc > s.x[last] ? s.y[last] : kNaN;  // the run stopped (at CSS) before that cycle
    }
    return InterpLinear(s.x, s.y, x);
  }

  double Sigma(size_t d, size_t i) const {
    const ProblemDataset& ds = problem_.Datasets()[d];
    return ds.data.observable.uncertainty.empty() ? ds.scale : ds.data.observable.uncertainty[i];
  }

  double Weight(size_t d) const {
    const double n = static_cast<double>(std::max<size_t>(1, active_[d].size()));
    return std::sqrt(problem_.Datasets()[d].data.observable.weight / n);
  }

  double Weighted(size_t d, size_t i, double y_sim) const {
    return Weight(d) * (y_sim - problem_.Datasets()[d].data.y[i]) / Sigma(d, i);
  }

 private:
  std::shared_ptr<RunOutcome> RunOne(const Job& job) const {
    auto out = std::make_shared<RunOutcome>();
    const Case& cs = problem_.Cases()[job.case_index];
    try {
      Simulation sim;
      LoadCase(sim, cs);
      for (size_t p = 0; p < problem_.Parameters().size(); ++p) {
        const size_t u = problem_.UnknownIndex(p, job.case_index);
        sim.SetParameter(problem_.Parameters()[p].path, FromNormalised(problem_.Unknowns()[u], job.z[u]));
      }
      RunOptions ro;
      if (cs.policy.cycles > 0) ro.cycles = cs.policy.cycles;
      ro.cancel = options_.cancel;
      CachedRunInfo info;
      const RunResult r = RunFromState(sim, job.start, ro, &info, job.distance);
      if (r.status == RunStatus::Cancelled) {
        out->cancelled = true;
        out->message = "cancelled";
        return out;
      }
      if (!r.Ok()) {
        out->message = "the run failed: " + r.message;
        return out;
      }
      if (cs.policy.require_css && !r.css_reached) {
        out->message = "the run did not reach cyclic steady state in " + std::to_string(r.cycles_run) + " cycles";
        return out;
      }
      const BedBoundarySources sources{&sim.Boundaries(), &sim.Data().cycle.step_list_};
      const std::vector<KPIRecord> kpis = sim.KPIs();
      for (size_t d : case_datasets_[job.case_index])
        out->series.push_back(SimulatedSeries(problem_.Datasets()[d].data, sim.Log(), kpis, sources));
      if (options_.warm_start && r.css_reached && sim.HasFinalState())
        out->final_state = std::make_unique<CycleState>(sim.FinalState());
      out->ok = true;
    } catch (const std::exception& e) {
      out->ok = false;
      out->message = e.what();
    }
    return out;
  }

  void SetActiveSets(PointEval& e) {
    active_.assign(problem_.Datasets().size(), {});
    for (size_t c = 0; c < NCases(); ++c) {
      const auto& ids = case_datasets_[c];
      for (size_t k = 0; k < ids.size(); ++k) {
        const ProblemDataset& d = problem_.Datasets()[ids[k]];
        const SimulatedCounterpart& s = e.cases[c]->series[k];
        const bool kpi = d.data.observable.kind == ObservableKind::Kpi;
        double lo = 0.0, hi = 0.0;
        if (!s.x.empty()) {
          lo = *std::min_element(s.x.begin(), s.x.end());
          hi = *std::max_element(s.x.begin(), s.x.end());
        }
        for (size_t i = 0; i < d.data.x.size(); ++i) {
          if (!std::isfinite(d.data.y[i]) || !std::isfinite(d.data.x[i]) || s.x.empty()) continue;
          if (!kpi && (d.data.x[i] < lo || d.data.x[i] > hi)) continue;
          if (!std::isfinite(SimulatedAt(d, s, d.data.x[i]))) continue;
          active_[ids[k]].push_back(i);
        }
        if (active_[ids[k]].empty()) {
          e.message = "case " + Quote(problem_.Cases()[c].name) + ", dataset " + Quote(d.data.name) +
                      ": no measured point is covered by the simulation";
          active_.clear();
          return;
        }
        const size_t skipped = d.data.x.size() - active_[ids[k]].size();
        if (skipped > 0)
          Note("case " + Quote(problem_.Cases()[c].name) + ", dataset " + Quote(d.data.name) + ": " +
               std::to_string(skipped) + " of " + std::to_string(d.data.x.size()) +
               " points are not covered by the simulation and are left out");
      }
    }
    have_active_ = true;
  }

  const Problem& problem_;
  const Options& options_;
  Result& result_;
  std::vector<std::unique_ptr<CycleStateCache>> caches_;
  std::vector<std::vector<size_t>> case_datasets_;  // [case] -> problem dataset indices
  std::vector<std::vector<size_t>> affects_;        // [unknown] -> cases
  std::vector<std::vector<size_t>> active_;         // [dataset] -> point indices
  bool have_active_ = false;
  bool cancelled_ = false;
  int threads_ = 1;
};

// ---- the solvers -------------------------------------------------------------

struct SolverState {
  PointEval best;        // the current (best) point
  Matrix J;              // the Jacobian at best.z when has_jacobian
  bool has_jacobian = false;
  StopReason stop = StopReason::Converged;
  std::string message;
};

void LeastSquares(Engine& eng, const Options& opt, Result& result, SolverState& st) {
  const size_t p = eng.NUnknowns();
  const int n_cases = static_cast<int>(eng.NCases());
  auto budget_left = [&](int runs) { return eng.Simulations() + runs <= opt.max_simulations; };
  auto stop = [&](StopReason r, std::string m) {
    st.stop = r;
    st.message = std::move(m);
  };

  if (!st.best.ok) {
    st.best = eng.EvaluatePoint(st.best.z);
    if (st.best.cancelled || eng.Cancelled()) return stop(StopReason::Cancelled, "cancelled at the starting point");
    if (!st.best.ok) return stop(StopReason::StartFailed, "the run at the starting point failed: " + st.best.message);
    eng.Record(Solver::LeastSquares, 0, st.best.objective, st.best.z);
  }
  st.has_jacobian = false;
  if (!budget_left(eng.JacobianCost()))
    return stop(StopReason::Budget, "the budget does not allow a Jacobian (" + std::to_string(eng.JacobianCost()) +
                                        " runs)");
  if (!eng.Jacobian(st.best, st.J)) return stop(StopReason::Cancelled, "cancelled");
  st.has_jacobian = true;

  double lambda = -1.0, nu = 2.0;
  bool rejected = false;  // the last step was rejected
  const auto zvec = [](const std::vector<double>& z) {
    return Eigen::Map<const Vector>(z.data(), static_cast<Eigen::Index>(z.size()));
  };
  for (int iter = 1;; ++iter) {
    if (iter > opt.max_iterations)
      return stop(StopReason::Iterations, "reached " + std::to_string(opt.max_iterations) + " iterations");
    if (eng.PollCancel()) return stop(StopReason::Cancelled, "cancelled");
    const Vector z = zvec(st.best.z);
    const Vector g = st.J.transpose() * st.best.r;
    std::vector<Eigen::Index> free;
    for (size_t j = 0; j < p; ++j) {
      const auto jj = static_cast<Eigen::Index>(j);
      const bool held = (z(jj) <= 0.0 && g(jj) > 0.0) || (z(jj) >= 1.0 && g(jj) < 0.0);
      if (!held) free.push_back(jj);
    }
    if (free.empty()) return stop(StopReason::Converged, "every parameter is held at a bound");
    const auto nf = static_cast<Eigen::Index>(free.size());
    Matrix JF(st.J.rows(), nf);
    Vector gF(nf);
    for (Eigen::Index k = 0; k < nf; ++k) {
      JF.col(k) = st.J.col(free[static_cast<size_t>(k)]);
      gF(k) = g(free[static_cast<size_t>(k)]);
    }
    const Matrix A = JF.transpose() * JF;
    const double dmax = A.diagonal().maxCoeff();
    if (!(dmax > 0.0)) return stop(StopReason::NoImprovement, "the residuals do not depend on the parameters");
    if (lambda < 0.0) lambda = 1e-3;
    Vector D = A.diagonal().cwiseMax(1e-10 * dmax);
    Matrix M = A;
    M.diagonal() += lambda * D;
    const Vector dF = M.ldlt().solve(-gF);
    Vector d = Vector::Zero(static_cast<Eigen::Index>(p));
    for (Eigen::Index k = 0; k < nf; ++k) d(free[static_cast<size_t>(k)]) = dF(k);
    const Vector zt = (z + d).cwiseMax(0.0).cwiseMin(1.0);
    d = zt - z;
    if (d.lpNorm<Eigen::Infinity>() < opt.xtol) {
      // After rejected steps the damping is what made the step small.
      if (rejected) return stop(StopReason::NoImprovement, "no step lowers the objective (the steps fell below xtol)");
      return stop(StopReason::SmallStep, "the step is below xtol");
    }
    if (!budget_left(n_cases)) return stop(StopReason::Budget, "reached the budget of " +
                                                                   std::to_string(opt.max_simulations) + " runs");
    std::vector<double> ztv(zt.data(), zt.data() + zt.size());
    PointEval trial = eng.EvaluatePoint(ztv);
    if (trial.cancelled || eng.Cancelled()) return stop(StopReason::Cancelled, "cancelled");
    const double predicted = st.best.objective - (st.best.r + st.J * d).squaredNorm();
    const double actual = trial.ok ? st.best.objective - trial.objective : -1.0;
    const double rho = predicted > 0.0 ? actual / predicted : -1.0;
    ++result.iterations;
    if (trial.ok && actual > 0.0 && rho > 1e-4) {
      const double rel = actual / std::max(st.best.objective, std::numeric_limits<double>::min());
      st.best = std::move(trial);
      st.has_jacobian = false;
      lambda *= std::max(1.0 / 3.0, 1.0 - std::pow(2.0 * rho - 1.0, 3));
      nu = 2.0;
      rejected = false;
      eng.Record(Solver::LeastSquares, iter, st.best.objective, st.best.z);
      if (rel < opt.ftol) return stop(StopReason::Converged, "the objective changed by less than ftol");
      if (st.best.objective == 0.0) return stop(StopReason::Converged, "the objective is zero");
      if (!budget_left(eng.JacobianCost()))
        return stop(StopReason::Budget, "reached the budget of " + std::to_string(opt.max_simulations) + " runs");
      if (!eng.Jacobian(st.best, st.J)) return stop(StopReason::Cancelled, "cancelled");
      st.has_jacobian = true;
    } else {
      lambda *= nu;
      nu *= 2.0;
      rejected = true;
      eng.Record(Solver::LeastSquares, iter, st.best.objective, st.best.z);
      if (lambda > 1e16) return stop(StopReason::NoImprovement, "no step lowers the objective");
    }
  }
}

void Global(Engine& eng, const Options& opt, Result& result, SolverState& st) {
  const size_t p = eng.NUnknowns();
  const int n_cases = static_cast<int>(eng.NCases());
  const int evaluations = std::min(opt.global_evaluations, opt.max_simulations / std::max(1, n_cases));
  if (evaluations < 1) {
    st.stop = StopReason::Budget;
    st.message = "the budget does not allow one candidate";
    return;
  }
  int index = 0;
  auto objective = [&](const std::vector<double>& z) -> double {
    if (eng.Cancelled() || index >= evaluations) return kOptimisationPenaltyCost;  // dlib asks for its whole budget
    ++index;
    ++result.iterations;
    PointEval e = eng.EvaluatePoint(z);
    const double cost = e.ok ? e.objective : kOptimisationPenaltyCost;
    if (e.ok && (!st.best.ok || e.objective < st.best.objective)) st.best = std::move(e);
    eng.Record(Solver::Global, index, st.best.ok ? st.best.objective : kNaN, st.best.ok ? st.best.z : z);
    return cost;
  };
  OptimizerSettings search;
  search.n_init_samples = 0;
  search.n_iterations = evaluations;
  GlobalOptimizer optimizer(std::vector<double>(p, 0.0), std::vector<double>(p, 1.0), objective, search);
  optimizer.OptimizeParameters();
  st.has_jacobian = false;
  if (eng.Cancelled()) {
    st.stop = StopReason::Cancelled;
    st.message = "cancelled";
  } else if (!st.best.ok) {
    st.stop = StopReason::StartFailed;
    st.message = "every run of the global search failed";
  } else {
    st.stop = StopReason::GlobalDone;
    st.message = "the global search used its " + std::to_string(evaluations) + " candidates";
  }
}

// ---- the result --------------------------------------------------------------

void StandardErrors(const Problem& problem, const SolverState& st, Result& result) {
  const size_t p = problem.Unknowns().size();
  const size_t n = static_cast<size_t>(st.best.r.size());
  result.correlation.assign(p, std::vector<double>(p, kNaN));
  if (!st.has_jacobian) {
    result.notes.push_back("no Jacobian at the estimate: the standard errors are not available");
    return;
  }
  if (n <= p) {
    result.notes.push_back("the standard errors need more residuals (" + std::to_string(n) + ") than parameters (" +
                           std::to_string(p) + ")");
    return;
  }
  result.sigma2 = st.best.objective / static_cast<double>(n - p);
  // The Jacobian in each parameter's scale (x, or ln x).
  Matrix Js = st.J;
  for (size_t j = 0; j < p; ++j) Js.col(static_cast<Eigen::Index>(j)) /= ScaledRange(problem.Unknowns()[j]);
  const Vector norms = Js.colwise().norm();
  for (size_t j = 0; j < p; ++j) {
    if (!(norms(static_cast<Eigen::Index>(j)) > 0.0)) {
      result.notes.push_back("the residuals do not depend on " + problem.Unknowns()[j].label +
                             ": J^T J is singular and the standard errors are not available");
      return;
    }
  }
  // Normalise the columns before the inverse (the scales differ by orders).
  const Matrix C = norms.cwiseInverse().asDiagonal() * (Js.transpose() * Js) * norms.cwiseInverse().asDiagonal();
  Eigen::SelfAdjointEigenSolver<Matrix> eig(C);
  const double emin = eig.eigenvalues().minCoeff(), emax = eig.eigenvalues().maxCoeff();
  if (!(emin > 1e-12 * emax)) {
    result.notes.push_back("J^T J is rank-deficient (condition number " + std::to_string(emax / std::max(emin, 0.0)) +
                           "): some parameters change the residuals in the same way, and the standard errors are "
                           "not available");
    return;
  }
  const Matrix Cinv = eig.eigenvectors() * eig.eigenvalues().cwiseInverse().asDiagonal() *
                      eig.eigenvectors().transpose();
  const Matrix cov = result.sigma2 * (norms.cwiseInverse().asDiagonal() * Cinv * norms.cwiseInverse().asDiagonal());
  for (size_t i = 0; i < p; ++i) {
    const auto ii = static_cast<Eigen::Index>(i);
    const double se_s = std::sqrt(cov(ii, ii));
    const Unknown& u = problem.Unknowns()[i];
    result.parameters[i].standard_error = u.scale == Scale::Log ? result.parameters[i].estimate * se_s : se_s;
    for (size_t j = 0; j < p; ++j) {
      const auto jj = static_cast<Eigen::Index>(j);
      result.correlation[i][j] = cov(ii, jj) / std::sqrt(cov(ii, ii) * cov(jj, jj));
    }
  }
}

void Finish(const Problem& problem, Engine& eng, const SolverState& st, Result& result) {
  const size_t p = problem.Unknowns().size();
  const std::vector<double> z = st.best.z;
  const std::vector<double> x = eng.Values(z);
  result.parameters.assign(p, ParameterEstimate{});
  for (size_t j = 0; j < p; ++j) {
    const Unknown& u = problem.Unknowns()[j];
    ParameterEstimate& e = result.parameters[j];
    e.path = u.path;
    e.label = u.label;
    e.case_index = u.case_index;
    e.estimate = x[j];
    e.initial = u.initial;
    e.lower = u.lower;
    e.upper = u.upper;
    e.scale = u.scale;
    e.at_lower = z[j] <= kBoundTolerance;
    e.at_upper = z[j] >= 1.0 - kBoundTolerance;
  }
  result.stop = st.stop;
  result.message = st.message;
  result.cancelled = st.stop == StopReason::Cancelled;
  if (!st.best.ok) {
    result.correlation.assign(p, std::vector<double>(p, kNaN));
    return;
  }
  result.objective = st.best.objective;
  result.residuals = static_cast<size_t>(st.best.r.size());
  StandardErrors(problem, st, result);

  // Per dataset, at the estimate.
  const auto& active = eng.ActiveSets();
  for (size_t c = 0; c < problem.Cases().size(); ++c) {
    const auto& ids = eng.CaseDatasets()[c];
    for (size_t k = 0; k < ids.size(); ++k) {
      const size_t d = ids[k];
      const ProblemDataset& ds = problem.Datasets()[d];
      DatasetFit f;
      f.case_index = c;
      f.case_name = problem.Cases()[c].name;
      f.name = ds.data.name;
      f.weight = ds.data.observable.weight;
      f.series = st.best.cases[c]->series[k];
      f.points = active[d].size();
      f.skipped = ds.data.x.size() - f.points;
      double sse = 0.0, wsse = 0.0, max_abs = 0.0;
      for (size_t i : active[d]) {
        const double ys = eng.SimulatedAt(ds, f.series, ds.data.x[i]);
        const double res = ys - ds.data.y[i];
        const double w = eng.Weighted(d, i, ys);
        f.x.push_back(ds.data.x[i]);
        f.measured.push_back(ds.data.y[i]);
        f.simulated.push_back(ys);
        f.residuals.push_back(res);
        f.weighted_residuals.push_back(w);
        sse += res * res;
        wsse += w * w;
        max_abs = std::max(max_abs, std::abs(res));
      }
      if (f.points > 0) {
        f.rmse = std::sqrt(sse / static_cast<double>(f.points));
        f.max_abs_error = max_abs;
        f.objective = wsse;
        double lo = std::numeric_limits<double>::infinity(), hi = -lo;
        for (double y : ds.data.y)
          if (std::isfinite(y)) {
            lo = std::min(lo, y);
            hi = std::max(hi, y);
          }
        if (hi > lo) f.nrmse = f.rmse / (hi - lo);
      }
      result.datasets.push_back(std::move(f));
    }
  }

  // The calibrated decks.
  for (size_t c = 0; c < problem.Cases().size(); ++c) {
    try {
      Simulation sim;
      sim.LoadJson(problem.Cases()[c].deck_json);
      for (size_t q = 0; q < problem.Parameters().size(); ++q)
        sim.SetParameter(problem.Parameters()[q].path, x[problem.UnknownIndex(q, c)]);
      result.calibrated_decks.push_back(sim.ToJson());
    } catch (const std::exception& e) {
      result.calibrated_decks.emplace_back();
      result.notes.push_back("case " + Quote(problem.Cases()[c].name) + ": the calibrated deck was not made: " +
                             e.what());
    }
  }
}

void CheckOptions(const Options& o) {
  if (o.max_simulations < 1) throw std::invalid_argument("max_simulations must be >= 1");
  if (o.max_iterations < 0) throw std::invalid_argument("max_iterations must be >= 0");
  if (!(o.ftol >= 0.0) || !(o.xtol >= 0.0)) throw std::invalid_argument("ftol and xtol must be >= 0");
  if (!(o.fd_step > 0.0 && o.fd_step < 0.5)) throw std::invalid_argument("fd_step must be in (0, 0.5)");
  if (o.threads < 0) throw std::invalid_argument("threads must be >= 0 (0: automatic)");
  if (o.solver == Solver::Global && o.global_evaluations < 1)
    throw std::invalid_argument("global_evaluations must be >= 1");
  if (o.cache_capacity < 1) throw std::invalid_argument("cache_capacity must be >= 1");
}

std::vector<double> StartPoint(const Problem& problem) {
  std::vector<double> z;
  for (const auto& u : problem.Unknowns()) z.push_back(ToNormalised(u, u.initial));
  return z;
}

}  // namespace

Result Calibrate(const Problem& problem, const Options& options) {
  CheckOptions(options);
  const auto t0 = std::chrono::steady_clock::now();
  Result result;
  result.solver = options.solver;
  result.notes = problem.Notes();
  Engine eng(problem, options, result);
  SolverState st;
  st.best.z = StartPoint(problem);

  if (options.solver == Solver::Global) {
    Global(eng, options, result, st);
    if (options.refine && st.stop == StopReason::GlobalDone) {
      result.refined = true;
      LeastSquares(eng, options, result, st);
    }
  } else {
    LeastSquares(eng, options, result, st);
  }

  // The Jacobian for the standard errors, when the solver has none at its estimate.
  if (options.standard_errors && st.best.ok && !st.has_jacobian && st.stop != StopReason::Cancelled &&
      !eng.Cancelled()) {
    if (eng.Jacobian(st.best, st.J)) {
      st.has_jacobian = true;
    } else {
      st.stop = StopReason::Cancelled;
      st.message = "cancelled";
    }
  }
  Finish(problem, eng, st, result);
  result.wall_time_s = std::chrono::duration<double>(std::chrono::steady_clock::now() - t0).count();
  return result;
}

Evaluation Evaluate(const Problem& problem, const std::vector<double>& values, const Options& options) {
  CheckOptions(options);
  if (values.size() != problem.Unknowns().size())
    throw std::invalid_argument("Evaluate: " + std::to_string(values.size()) + " values for " +
                                std::to_string(problem.Unknowns().size()) + " unknowns");
  std::vector<double> z;
  for (size_t j = 0; j < values.size(); ++j) {
    const Unknown& u = problem.Unknowns()[j];
    if (!(values[j] >= u.lower && values[j] <= u.upper))
      throw std::invalid_argument("Evaluate: the value of " + u.label + " is outside its bounds");
    z.push_back(ToNormalised(u, values[j]));
  }
  Result scratch;
  Options o = options;
  o.warm_start = false;
  Engine eng(problem, o, scratch);
  PointEval e = eng.EvaluatePoint(z);
  Evaluation out;
  out.ok = e.ok;
  out.message = e.message;
  if (e.ok) {
    out.objective = e.objective;
    out.residuals.assign(e.r.data(), e.r.data() + e.r.size());
  }
  return out;
}

}  // namespace calibration
}  // namespace adsim
