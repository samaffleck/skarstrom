/*
 * Skarstrom - pressure and temperature swing adsorption simulator
 * Copyright (c) 2024-2026 Sam Affleck and contributors
 *
 * Licensed under the terms in the LICENSE file at the repository root.
 */
#include "adsim/solver/sensitivity.h"

#include <algorithm>
#include <cmath>
#include <numeric>
#include <random>
#include <stdexcept>

#include <Eigen/Dense>

namespace adsim {
namespace sensitivity {

namespace {

// ---- random numbers ------------------------------------------------------------
// Every draw is made here from the raw 64-bit output of std::mt19937_64, whose
// sequence the standard fixes, so a seed gives the same plan on every
// platform (native, Python wheels, WebAssembly).

std::uint64_t SplitMix64(std::uint64_t x) {
  x += 0x9e3779b97f4a7c15ULL;
  x = (x ^ (x >> 30)) * 0xbf58476d1ce4e5b9ULL;
  x = (x ^ (x >> 27)) * 0x94d049bb133111ebULL;
  return x ^ (x >> 31);
}

// A seed for one independent stream (a Monte Carlo batch, one bootstrap).
std::uint64_t StreamSeed(std::uint64_t seed, std::uint64_t stream) {
  return SplitMix64(SplitMix64(seed) ^ SplitMix64(stream + 0x5851f42d4c957f2dULL));
}

class Rng {
public:
  explicit Rng(std::uint64_t seed) : gen_(seed) {}

  // Uniform in the open interval (0, 1).
  double Open01() { return (static_cast<double>(gen_() >> 11) + 0.5) * 0x1.0p-53; }

  // Uniform integer in [0, n), without modulo bias.
  std::uint64_t Below(std::uint64_t n) {
    const std::uint64_t limit = std::numeric_limits<std::uint64_t>::max() -
                                std::numeric_limits<std::uint64_t>::max() % n;
    std::uint64_t x;
    do {
      x = gen_();
    } while (x >= limit);
    return x % n;
  }

  // Fisher-Yates shuffle of 0 .. n-1.
  std::vector<int> Permutation(int n) {
    std::vector<int> p(static_cast<size_t>(n));
    std::iota(p.begin(), p.end(), 0);
    for (int i = n - 1; i > 0; --i) {
      std::swap(p[static_cast<size_t>(i)], p[Below(static_cast<std::uint64_t>(i) + 1)]);
    }
    return p;
  }

private:
  std::mt19937_64 gen_;
};

// ---- the normal distribution ------------------------------------------------------

double NormalCdf(double x) { return 0.5 * std::erfc(-x / std::sqrt(2.0)); }

// Inverse of NormalCdf: Acklam's rational approximation (relative error
// 1.15e-9) refined by one Halley step, which brings it to double precision.
double NormalQuantile(double p) {
  static const double a[] = {-3.969683028665376e+01, 2.209460984245205e+02, -2.759285104469687e+02,
                             1.383577518672690e+02, -3.066479806614716e+01, 2.506628277459239e+00};
  static const double b[] = {-5.447609879822406e+01, 1.615858368580409e+02, -1.556989798598866e+02,
                             6.680131188771972e+01, -1.328068155288572e+01};
  static const double c[] = {-7.784894002430293e-03, -3.223964580411365e-01, -2.400758277161838e+00,
                             -2.549732539343734e+00, 4.374664141464968e+00, 2.938163982698783e+00};
  static const double d[] = {7.784695709041462e-03, 3.224671290700398e-01, 2.445134137142996e+00,
                             3.754408661907416e+00};
  const double plow = 0.02425;
  double x;
  if (p < plow) {
    const double q = std::sqrt(-2.0 * std::log(p));
    x = (((((c[0] * q + c[1]) * q + c[2]) * q + c[3]) * q + c[4]) * q + c[5]) /
        ((((d[0] * q + d[1]) * q + d[2]) * q + d[3]) * q + 1.0);
  } else if (p <= 1.0 - plow) {
    const double q = p - 0.5;
    const double r = q * q;
    x = (((((a[0] * r + a[1]) * r + a[2]) * r + a[3]) * r + a[4]) * r + a[5]) * q /
        (((((b[0] * r + b[1]) * r + b[2]) * r + b[3]) * r + b[4]) * r + 1.0);
  } else {
    const double q = std::sqrt(-2.0 * std::log(1.0 - p));
    x = -(((((c[0] * q + c[1]) * q + c[2]) * q + c[3]) * q + c[4]) * q + c[5]) /
        ((((d[0] * q + d[1]) * q + d[2]) * q + d[3]) * q + 1.0);
  }
  const double e = NormalCdf(x) - p;
  const double u = e * 2.5066282746310002 * std::exp(0.5 * x * x);
  return x - u / (1.0 + 0.5 * x * u);
}

// ---- statistics -----------------------------------------------------------------

double Mean(const std::vector<double>& x) {
  return x.empty() ? kNaN : std::accumulate(x.begin(), x.end(), 0.0) / static_cast<double>(x.size());
}

double SampleSd(const std::vector<double>& x) {
  if (x.size() < 2) return kNaN;
  const double m = Mean(x);
  double s = 0.0;
  for (double v : x) s += (v - m) * (v - m);
  return std::sqrt(s / static_cast<double>(x.size() - 1));
}

// Percentile of an already sorted, finite vector.
double SortedPercentile(const std::vector<double>& s, double q) {
  if (s.empty()) return kNaN;
  const double pos = q * static_cast<double>(s.size() - 1);
  const size_t lo = static_cast<size_t>(std::floor(pos));
  const size_t hi = std::min(lo + 1, s.size() - 1);
  return s[lo] + (pos - static_cast<double>(lo)) * (s[hi] - s[lo]);
}

// Ranks 1..n, ties get their average rank.
std::vector<double> Ranks(const std::vector<double>& x) {
  const size_t n = x.size();
  std::vector<size_t> order(n);
  std::iota(order.begin(), order.end(), size_t{0});
  std::stable_sort(order.begin(), order.end(), [&](size_t a, size_t b) { return x[a] < x[b]; });
  std::vector<double> r(n);
  for (size_t i = 0; i < n;) {
    size_t j = i;
    while (j + 1 < n && x[order[j + 1]] == x[order[i]]) ++j;
    const double rank = 0.5 * static_cast<double>(i + j) + 1.0;
    for (size_t t = i; t <= j; ++t) r[order[t]] = rank;
    i = j + 1;
  }
  return r;
}

// Centred and scaled to unit sample standard deviation; false when constant.
bool Standardise(std::vector<double>& x) {
  const double m = Mean(x);
  const double sd = SampleSd(x);
  if (!(sd > 0.0)) return false;
  for (double& v : x) v = (v - m) / sd;
  return true;
}

double Correlation(const std::vector<double>& a, const std::vector<double>& b) {
  const double ma = Mean(a), mb = Mean(b);
  double sab = 0.0, saa = 0.0, sbb = 0.0;
  for (size_t i = 0; i < a.size(); ++i) {
    sab += (a[i] - ma) * (b[i] - mb);
    saa += (a[i] - ma) * (a[i] - ma);
    sbb += (b[i] - mb) * (b[i] - mb);
  }
  return (saa > 0.0 && sbb > 0.0) ? sab / std::sqrt(saa * sbb) : kNaN;
}

// 95 % percentile-bootstrap interval of a statistic of `x`.
template <typename Stat>
Interval Bootstrap(const std::vector<double>& x, int resamples, std::uint64_t seed, Stat stat) {
  Interval out;
  if (x.empty() || resamples < 1) return out;
  Rng rng(seed);
  std::vector<double> stats;
  stats.reserve(static_cast<size_t>(resamples));
  std::vector<double> sample(x.size());
  for (int b = 0; b < resamples; ++b) {
    for (double& v : sample) v = x[rng.Below(x.size())];
    stats.push_back(stat(sample));
  }
  std::sort(stats.begin(), stats.end());
  out.low = SortedPercentile(stats, 0.025);
  out.high = SortedPercentile(stats, 0.975);
  return out;
}

void CheckInputs(const std::vector<Input>& inputs) {
  if (inputs.empty()) throw std::invalid_argument("a sensitivity study needs at least one input");
  for (const auto& in : inputs) in.Validate();
}

void CheckShape(const std::vector<std::vector<double>>& m, size_t rows, size_t cols, const char* what) {
  if (m.size() != rows) {
    throw std::invalid_argument(std::string(what) + ": expected " + std::to_string(rows) + " rows, got " +
                                std::to_string(m.size()));
  }
  for (const auto& row : m) {
    if (row.size() != cols) {
      throw std::invalid_argument(std::string(what) + ": expected " + std::to_string(cols) +
                                  " columns, got " + std::to_string(row.size()));
    }
  }
}

double KPIValue(const KPIRecord& k, const std::string& kpi) {
  if (kpi == "purity") return k.purity;
  if (kpi == "recovery") return k.recovery;
  if (kpi == "productivity_mass") return k.productivity_mass;
  if (kpi == "productivity_volume") return k.productivity_volume;
  if (kpi == "energy_consumption") return k.energy_consumption;
  if (kpi == "production_rate") return k.production_rate;
  throw std::invalid_argument("unknown KPI '" + kpi + "'");
}

bool IsTracked(const UserData& data, const std::string& component) {
  for (const auto& cfg : data.kpi_config.per_component) {
    if (cfg.tracked && cfg.product_component == component) return true;
  }
  return false;
}

std::string Join(const std::vector<std::string>& v) {
  std::string s;
  for (const auto& x : v) s += (s.empty() ? "" : ", ") + x;
  return s;
}

}  // namespace

// ---- inputs -----------------------------------------------------------------------

const char* DistributionName(Distribution d) {
  switch (d) {
    case Distribution::Uniform: return "uniform";
    case Distribution::Triangular: return "triangular";
    case Distribution::Normal: return "normal";
  }
  return "uniform";
}

Distribution DistributionFromName(std::string_view name) {
  if (name == "uniform") return Distribution::Uniform;
  if (name == "triangular") return Distribution::Triangular;
  if (name == "normal") return Distribution::Normal;
  throw std::invalid_argument("unknown distribution '" + std::string(name) +
                              "' (expected uniform, triangular or normal)");
}

double Input::Quantile(double u) const {
  switch (distribution) {
    case Distribution::Uniform:
      return lower + u * (upper - lower);
    case Distribution::Triangular: {
      const double width = upper - lower;
      const double c = (mode - lower) / width;
      return u < c ? lower + std::sqrt(u * width * (mode - lower))
                   : upper - std::sqrt((1.0 - u) * width * (upper - mode));
    }
    case Distribution::Normal: {
      const double pa = std::isnan(lower) ? 0.0 : NormalCdf((lower - mean) / sd);
      const double pb = std::isnan(upper) ? 1.0 : NormalCdf((upper - mean) / sd);
      const double x = mean + sd * NormalQuantile(pa + u * (pb - pa));
      // Round-off at an extreme truncation point must not leave the interval.
      return std::clamp(x, std::isnan(lower) ? x : lower, std::isnan(upper) ? x : upper);
    }
  }
  return kNaN;
}

void Input::Validate() const {
  const std::string who = "input '" + path + "': ";
  if (path.empty()) throw std::invalid_argument("an input has no parameter path");
  auto finite = [](double v) { return std::isfinite(v); };
  switch (distribution) {
    case Distribution::Uniform:
      if (!finite(lower) || !finite(upper) || !(lower < upper)) {
        throw std::invalid_argument(who + "a uniform range needs finite lower < upper");
      }
      break;
    case Distribution::Triangular:
      if (!finite(lower) || !finite(upper) || !finite(mode) || !(lower < upper) || mode < lower ||
          mode > upper) {
        throw std::invalid_argument(who + "a triangular distribution needs finite lower <= mode <= upper, lower < upper");
      }
      break;
    case Distribution::Normal: {
      if (!finite(mean) || !finite(sd) || !(sd > 0.0)) {
        throw std::invalid_argument(who + "a normal distribution needs a finite mean and sd > 0");
      }
      if ((!std::isnan(lower) && !finite(lower)) || (!std::isnan(upper) && !finite(upper))) {
        throw std::invalid_argument(who + "truncation bounds must be finite (or left unset)");
      }
      if (!std::isnan(lower) && !std::isnan(upper) && !(lower < upper)) {
        throw std::invalid_argument(who + "truncation needs lower < upper");
      }
      const double pa = std::isnan(lower) ? 0.0 : NormalCdf((lower - mean) / sd);
      const double pb = std::isnan(upper) ? 1.0 : NormalCdf((upper - mean) / sd);
      if (!(pb - pa > 1e-12)) {
        throw std::invalid_argument(who + "the truncation bounds leave no probability (they lie far in one tail)");
      }
      break;
    }
  }
}

Input Uniform(std::string path, double lower, double upper) {
  Input in;
  in.path = std::move(path);
  in.distribution = Distribution::Uniform;
  in.lower = lower;
  in.upper = upper;
  return in;
}

Input Relative(std::string path, double nominal, double fraction) {
  if (!(fraction > 0.0) || !std::isfinite(fraction)) {
    throw std::invalid_argument("input '" + path + "': the relative range must be > 0");
  }
  if (!(nominal != 0.0) || !std::isfinite(nominal)) {
    throw std::invalid_argument("input '" + path + "': a relative range needs a finite, non-zero nominal value; give bounds instead");
  }
  const double a = nominal * (1.0 - fraction);
  const double b = nominal * (1.0 + fraction);
  return Uniform(std::move(path), std::min(a, b), std::max(a, b));
}

Input Triangular(std::string path, double lower, double mode, double upper) {
  Input in;
  in.path = std::move(path);
  in.distribution = Distribution::Triangular;
  in.lower = lower;
  in.mode = mode;
  in.upper = upper;
  return in;
}

Input Normal(std::string path, double mean, double sd, double lower, double upper) {
  Input in;
  in.path = std::move(path);
  in.distribution = Distribution::Normal;
  in.mean = mean;
  in.sd = sd;
  in.lower = lower;
  in.upper = upper;
  return in;
}

// ---- plans --------------------------------------------------------------------------

const char* MethodName(Method m) { return m == Method::Morris ? "morris" : "monte_carlo"; }

Plan MorrisPlan(const std::vector<Input>& inputs, int trajectories, int levels, std::uint64_t seed) {
  CheckInputs(inputs);
  if (trajectories < 1) throw std::invalid_argument("Morris screening needs at least one trajectory");
  if (levels < 2 || levels % 2 != 0) throw std::invalid_argument("Morris levels must be even and >= 2");

  Plan plan;
  plan.method = Method::Morris;
  plan.inputs = inputs;
  plan.seed = seed;
  plan.trajectories = trajectories;
  plan.levels = levels;

  const int k = static_cast<int>(inputs.size());
  const int jump = levels / 2;  // |du| = jump / levels = 1/2
  Rng rng(StreamSeed(seed, 0x4d6f72726973ULL));  // "Morris"
  auto push = [&](const std::vector<int>& level, int changed) {
    std::vector<double> u(level.size()), v(level.size());
    for (size_t i = 0; i < level.size(); ++i) {
      u[i] = (static_cast<double>(level[i]) + 0.5) / static_cast<double>(levels);
      v[i] = inputs[i].Quantile(u[i]);
    }
    plan.u.push_back(std::move(u));
    plan.values.push_back(std::move(v));
    plan.changed.push_back(changed);
  };
  for (int t = 0; t < trajectories; ++t) {
    std::vector<int> level(static_cast<size_t>(k));
    for (int& l : level) l = static_cast<int>(rng.Below(static_cast<std::uint64_t>(levels)));
    const std::vector<int> order = rng.Permutation(k);
    push(level, -1);
    for (int i : order) {
      int& l = level[static_cast<size_t>(i)];
      l = (l < jump) ? l + jump : l - jump;
      push(level, i);
    }
  }
  return plan;
}

Plan MonteCarloPlan(const std::vector<Input>& inputs, int runs, std::uint64_t seed, int batch) {
  CheckInputs(inputs);
  if (runs < 1) throw std::invalid_argument("a Monte Carlo batch needs at least one run");
  if (batch < 0) throw std::invalid_argument("the batch index must be >= 0");

  Plan plan;
  plan.method = Method::MonteCarlo;
  plan.inputs = inputs;
  plan.seed = seed;
  plan.batch = batch;

  const size_t k = inputs.size();
  const size_t n = static_cast<size_t>(runs);
  plan.u.assign(n, std::vector<double>(k));
  plan.values.assign(n, std::vector<double>(k));
  Rng rng(StreamSeed(seed, static_cast<std::uint64_t>(batch)));
  for (size_t j = 0; j < k; ++j) {
    const std::vector<int> stratum = rng.Permutation(runs);
    for (size_t r = 0; r < n; ++r) {
      const double u = (static_cast<double>(stratum[r]) + rng.Open01()) / static_cast<double>(n);
      plan.u[r][j] = u;
      plan.values[r][j] = inputs[j].Quantile(u);
    }
  }
  return plan;
}

// ---- running one sample ----------------------------------------------------------------

const std::vector<std::string>& KPINames() {
  static const std::vector<std::string> names = {"purity",           "recovery",
                                                 "productivity_mass", "productivity_volume",
                                                 "energy_consumption", "production_rate"};
  return names;
}

Output ParseOutput(std::string_view name) {
  const size_t dot = name.rfind('.');
  if (dot == std::string_view::npos || dot == 0 || dot + 1 == name.size()) {
    throw std::invalid_argument("output '" + std::string(name) + "' must be Component.kpi (e.g. CarbonDioxide.purity)");
  }
  Output out{std::string(name.substr(0, dot)), std::string(name.substr(dot + 1))};
  const auto& kpis = KPINames();
  if (std::find(kpis.begin(), kpis.end(), out.kpi) == kpis.end()) {
    throw std::invalid_argument("output '" + std::string(name) + "': unknown KPI '" + out.kpi + "' (expected one of " +
                                Join(kpis) + ")");
  }
  return out;
}

std::vector<Output> DefaultOutputs(const UserData& data) {
  std::vector<Output> out;
  for (const auto& cfg : data.kpi_config.per_component) {
    if (!cfg.tracked) continue;
    for (const char* kpi : {"purity", "recovery", "productivity_mass", "energy_consumption"}) {
      out.push_back({cfg.product_component, kpi});
    }
  }
  return out;
}

namespace {

// Applies the policy's settings to a loaded deck.
void ApplySettings(Simulation& sim, const RunPolicy& policy) {
  for (const auto& [path, text] : policy.settings) {
    try {
      sim.SetParameterText(path, text);
    } catch (const std::exception& e) {
      throw std::invalid_argument("study setting " + path + " = '" + text + "': " + e.what());
    }
  }
}

}  // namespace

void CheckStudy(const std::string& deck_json, const std::vector<Input>& inputs, const std::vector<Output>& outputs,
                const RunPolicy& policy) {
  CheckInputs(inputs);
  Simulation sim;
  try {
    sim.LoadJson(deck_json);
  } catch (const std::exception& e) {
    throw std::invalid_argument(std::string("the deck does not load: ") + e.what());
  }
  ApplySettings(sim, policy);
  if (policy.require_css && sim.Data().settings.css.method_ == CSSMethod::NONE) {
    throw std::invalid_argument("the study requires cyclic steady state but the deck has no CSS check (css.method is "
                                "none); choose a CSS method or turn the CSS requirement off");
  }
  ParameterRegistry& registry = sim.Parameters();
  std::vector<std::string> seen;
  for (const auto& in : inputs) {
    const Parameter* p = nullptr;
    try {
      p = registry.Find(in.path);
    } catch (const std::exception& e) {
      throw std::invalid_argument("input '" + in.path + "': " + e.what());
    }
    if (!p) throw std::invalid_argument("input '" + in.path + "': no such parameter in the deck");
    if (!p->Info().tunable) {
      throw std::invalid_argument(
          "input '" + in.path +
          "' cannot be varied: only tunable values can (not solver settings, integers, text, or a value "
          "written as a formula of user variables; vary the variable instead)");
    }
    if (std::find(seen.begin(), seen.end(), p->Path()) != seen.end()) {
      throw std::invalid_argument("input '" + in.path + "' is listed twice");
    }
    seen.push_back(p->Path());
  }
  if (outputs.empty()) {
    throw std::invalid_argument("a sensitivity study needs at least one output (the deck tracks no KPIs; set up the "
                                "KPI configuration first)");
  }
  const auto& kpis = KPINames();
  for (const auto& o : outputs) {
    if (std::find(kpis.begin(), kpis.end(), o.kpi) == kpis.end()) {
      throw std::invalid_argument("output '" + o.Name() + "': unknown KPI '" + o.kpi + "'");
    }
    if (!IsTracked(sim.Data(), o.component)) {
      throw std::invalid_argument("output '" + o.Name() + "': " + o.component +
                                  " is not tracked in the deck's KPI configuration");
    }
  }
}

SampleOutcome EvaluateSample(const std::string& deck_json, const std::vector<Input>& inputs,
                             const std::vector<double>& values, const std::vector<Output>& outputs,
                             const RunPolicy& policy, const CancelPoll& cancel) {
  if (!values.empty() && values.size() != inputs.size()) {
    throw std::invalid_argument("a sample needs one value per input (" + std::to_string(inputs.size()) + "), got " +
                                std::to_string(values.size()));
  }
  Simulation sim;
  try {
    sim.LoadJson(deck_json);
  } catch (const std::exception& e) {
    throw std::invalid_argument(std::string("the deck does not load: ") + e.what());
  }
  ApplySettings(sim, policy);
  for (size_t i = 0; i < values.size(); ++i) {
    try {
      sim.SetParameter(inputs[i].path, values[i]);
    } catch (const std::exception& e) {
      throw std::invalid_argument("input '" + inputs[i].path + "': " + e.what());
    }
  }

  RunOptions options;
  if (policy.cycles > 0) options.cycles = policy.cycles;
  options.cancel = cancel;
  const RunResult r = sim.Run(options);

  SampleOutcome out;
  out.status = r.status;
  out.message = r.message;
  out.wall_time_s = r.wall_time_s;
  out.outputs.assign(outputs.size(), kNaN);
  out.accepted = r.Ok() && (!policy.require_css || r.css_reached);
  if (r.Ok() && policy.require_css && !r.css_reached) {
    out.message = "did not reach cyclic steady state in " + std::to_string(r.cycles_run) + " cycles";
  }
  if (!out.accepted) return out;

  const std::vector<KPIRecord> kpis = sim.KPIs(-1);
  for (size_t o = 0; o < outputs.size(); ++o) {
    for (const auto& k : kpis) {
      if (k.component == outputs[o].component) {
        out.outputs[o] = KPIValue(k, outputs[o].kpi);
        break;
      }
    }
  }
  return out;
}

// ---- analysis ----------------------------------------------------------------------------

double Percentile(std::vector<double> x, double q) {
  x.erase(std::remove_if(x.begin(), x.end(), [](double v) { return !std::isfinite(v); }), x.end());
  std::sort(x.begin(), x.end());
  return SortedPercentile(x, q);
}

MorrisResult AnalyseMorris(const Plan& plan, const std::vector<std::string>& output_names,
                           const std::vector<std::vector<double>>& y, int bootstrap) {
  if (plan.method != Method::Morris) throw std::invalid_argument("AnalyseMorris needs a Morris plan");
  const size_t runs = plan.Runs();
  CheckShape(y, runs, output_names.size(), "Morris results");
  if (plan.changed.size() != runs) throw std::invalid_argument("the Morris plan is inconsistent");

  MorrisResult result;
  result.runs = static_cast<int>(runs);
  const size_t k = plan.inputs.size();
  for (size_t o = 0; o < output_names.size(); ++o) {
    MorrisOutput mo;
    mo.output = output_names[o];
    for (size_t r = 0; r < runs; ++r) mo.missing += std::isfinite(y[r][o]) ? 0 : 1;

    std::vector<std::vector<double>> ee(k);
    for (size_t r = 1; r < runs; ++r) {
      const int i = plan.changed[r];
      if (i < 0) continue;
      const double a = y[r - 1][o], b = y[r][o];
      if (!std::isfinite(a) || !std::isfinite(b)) continue;
      const double du = plan.u[r][static_cast<size_t>(i)] - plan.u[r - 1][static_cast<size_t>(i)];
      ee[static_cast<size_t>(i)].push_back((b - a) / du);
    }
    for (size_t i = 0; i < k; ++i) {
      MorrisEffect e;
      e.input = plan.inputs[i].path;
      e.effects = static_cast<int>(ee[i].size());
      if (!ee[i].empty()) {
        std::vector<double> abs_ee(ee[i].size());
        std::transform(ee[i].begin(), ee[i].end(), abs_ee.begin(), [](double v) { return std::fabs(v); });
        e.mu = Mean(ee[i]);
        e.mu_star = Mean(abs_ee);
        e.sigma = SampleSd(ee[i]);
        e.mu_star_interval = Bootstrap(abs_ee, bootstrap, StreamSeed(plan.seed, 0x10000 + o * 4096 + i), Mean);
      }
      mo.effects.push_back(std::move(e));
    }
    result.outputs.push_back(std::move(mo));
  }
  return result;
}

MonteCarloResult AnalyseMonteCarlo(const std::vector<Input>& inputs, const std::vector<std::vector<double>>& values,
                                   const std::vector<std::string>& output_names,
                                   const std::vector<std::vector<double>>& y, std::uint64_t seed, int bootstrap,
                                   double rel_tol) {
  const size_t runs = values.size();
  const size_t k = inputs.size();
  CheckShape(values, runs, k, "Monte Carlo inputs");
  CheckShape(y, runs, output_names.size(), "Monte Carlo results");

  MonteCarloResult result;
  result.runs = static_cast<int>(runs);
  bool any = false, all = true;
  for (size_t o = 0; o < output_names.size(); ++o) {
    Spread s;
    s.output = output_names[o];
    std::vector<size_t> rows;
    for (size_t r = 0; r < runs; ++r) {
      if (std::isfinite(y[r][o])) rows.push_back(r);
    }
    s.valid = static_cast<int>(rows.size());
    s.missing = static_cast<int>(runs - rows.size());
    for (size_t i = 0; i < k; ++i) s.shares.push_back({inputs[i].path});
    if (rows.empty()) {
      result.outputs.push_back(std::move(s));
      continue;
    }

    std::vector<double> v;
    for (size_t r : rows) v.push_back(y[r][o]);
    std::vector<double> sorted = v;
    std::sort(sorted.begin(), sorted.end());
    s.mean = Mean(v);
    s.sd = SampleSd(v);
    s.min = sorted.front();
    s.max = sorted.back();
    s.p05 = SortedPercentile(sorted, 0.05);
    s.p50 = SortedPercentile(sorted, 0.50);
    s.p95 = SortedPercentile(sorted, 0.95);
    auto percentile = [](double q) {
      return [q](std::vector<double> x) {
        std::sort(x.begin(), x.end());
        return SortedPercentile(x, q);
      };
    };
    const std::uint64_t base = 0x20000 + o * 8;
    s.p05_interval = Bootstrap(v, bootstrap, StreamSeed(seed, base + 0), percentile(0.05));
    s.p50_interval = Bootstrap(v, bootstrap, StreamSeed(seed, base + 1), percentile(0.50));
    s.p95_interval = Bootstrap(v, bootstrap, StreamSeed(seed, base + 2), percentile(0.95));

    const double band = s.p95 - s.p05;
    if (!(band > 1e-12 * std::max(1.0, std::fabs(s.p50)))) {
      s.converged = s.valid >= 2;  // the output does not vary
    } else {
      auto narrow = [&](double p, const Interval& ci) {
        return std::max(p - ci.low, ci.high - p) <= rel_tol * band;
      };
      s.converged = narrow(s.p05, s.p05_interval) && narrow(s.p50, s.p50_interval) && narrow(s.p95, s.p95_interval);
    }
    any = true;
    all = all && s.converged;

    // Rank correlation and rank regression (on the runs with a value).
    std::vector<double> ry = Ranks(v);
    const size_t n = rows.size();
    std::vector<std::vector<double>> rx(k);
    for (size_t i = 0; i < k; ++i) {
      std::vector<double> col;
      for (size_t r : rows) col.push_back(values[r][i]);
      rx[i] = Ranks(col);
      s.shares[i].spearman = Correlation(rx[i], ry);
    }
    if (n > k + 1 && Standardise(ry)) {
      Eigen::MatrixXd X(static_cast<Eigen::Index>(n), static_cast<Eigen::Index>(k));
      Eigen::VectorXd Y(static_cast<Eigen::Index>(n));
      bool ok = true;
      for (size_t i = 0; i < k && ok; ++i) {
        ok = Standardise(rx[i]);
        for (size_t r = 0; r < n; ++r) X(static_cast<Eigen::Index>(r), static_cast<Eigen::Index>(i)) = rx[i][r];
      }
      for (size_t r = 0; r < n; ++r) Y(static_cast<Eigen::Index>(r)) = ry[r];
      if (ok) {
        const Eigen::VectorXd beta = X.colPivHouseholderQr().solve(Y);
        const double ss_res = (Y - X * beta).squaredNorm();
        const double ss_tot = Y.squaredNorm();  // Y is centred
        s.r2 = 1.0 - ss_res / ss_tot;
        for (size_t i = 0; i < k; ++i) {
          s.shares[i].srrc = beta(static_cast<Eigen::Index>(i));
          s.shares[i].share = s.shares[i].srrc * s.shares[i].srrc;
        }
      }
    }
    result.outputs.push_back(std::move(s));
  }
  result.converged = any && all;
  return result;
}

}  // namespace sensitivity
}  // namespace adsim
