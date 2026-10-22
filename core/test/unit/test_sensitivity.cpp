// Sensitivity analysis (adsim/solver/sensitivity.h): distributions, the
// Morris and Latin hypercube plans, and the analyses checked against
// functions whose answers are known in closed form.

#include <gtest/gtest.h>

#include <algorithm>
#include <cmath>
#include <fstream>
#include <set>
#include <sstream>
#include <string>
#include <vector>

#include "adsim/solver/sensitivity.h"

namespace adsim {
namespace sensitivity {
namespace {

std::vector<std::vector<double>> Evaluate(const Plan& plan,
                                          double (*f)(const std::vector<double>&)) {
  std::vector<std::vector<double>> y;
  for (const auto& x : plan.values) y.push_back({f(x)});
  return y;
}

// The Sobol g-function on U(0, 1)^6: input i matters less as a[i] grows.
double GFunction(const std::vector<double>& x) {
  static const double a[] = {0.0, 1.0, 4.5, 9.0, 99.0, 99.0};
  double g = 1.0;
  for (size_t i = 0; i < x.size(); ++i) g *= (std::fabs(4.0 * x[i] - 2.0) + a[i]) / (1.0 + a[i]);
  return g;
}

double Linear(const std::vector<double>& x) { return 3.0 * x[0] - 2.0 * x[1] + 0.0 * x[2]; }

double TwoNormals(const std::vector<double>& x) { return 2.0 * x[0] + x[1]; }

}  // namespace

TEST(Sensitivity, Quantiles) {
  EXPECT_DOUBLE_EQ(Uniform("a", 0.0, 10.0).Quantile(0.25), 2.5);

  const Input tri = Triangular("t", 0.0, 2.0, 10.0);  // CDF at the mode is 0.2
  EXPECT_NEAR(tri.Quantile(0.2), 2.0, 1e-12);
  EXPECT_NEAR(tri.Quantile(0.5), 10.0 - std::sqrt(40.0), 1e-12);

  const Input n = Normal("n", 5.0, 2.0);
  EXPECT_NEAR(n.Quantile(0.5), 5.0, 1e-12);
  EXPECT_NEAR(n.Quantile(0.975), 5.0 + 2.0 * 1.959963984540054, 1e-9);
  EXPECT_NEAR(n.Quantile(1e-6), 5.0 - 2.0 * 4.753424308822899, 1e-7);

  // Half-normal: truncated at the mean, the median is the 75th percentile.
  const Input half = Normal("h", 0.0, 1.0, 0.0);
  EXPECT_NEAR(half.Quantile(0.5), 0.6744897501960817, 1e-9);
  EXPECT_GE(half.Quantile(1e-9), 0.0);

  const Input rel = Relative("r", -2.0, 0.1);
  EXPECT_DOUBLE_EQ(rel.lower, -2.2);
  EXPECT_DOUBLE_EQ(rel.upper, -1.8);
}

TEST(Sensitivity, InvalidInputsAreRejected) {
  EXPECT_THROW(Uniform("a", 1.0, 1.0).Validate(), std::invalid_argument);
  EXPECT_THROW(Uniform("", 0.0, 1.0).Validate(), std::invalid_argument);
  EXPECT_THROW(Triangular("t", 0.0, 2.0, 1.0).Validate(), std::invalid_argument);
  EXPECT_THROW(Normal("n", 0.0, 0.0).Validate(), std::invalid_argument);
  EXPECT_THROW(Normal("n", 0.0, 1.0, 50.0, 60.0).Validate(), std::invalid_argument);
  EXPECT_THROW(Relative("r", 0.0, 0.1), std::invalid_argument);
  EXPECT_THROW(Relative("r", 1.0, 0.0), std::invalid_argument);
  EXPECT_THROW(DistributionFromName("lognormal"), std::invalid_argument);
  EXPECT_EQ(DistributionFromName(DistributionName(Distribution::Triangular)), Distribution::Triangular);
  EXPECT_THROW(MorrisPlan({}, 10), std::invalid_argument);
  EXPECT_THROW(MorrisPlan({Uniform("a", 0, 1)}, 10, 3), std::invalid_argument);
  EXPECT_THROW(MonteCarloPlan({Uniform("a", 0, 1)}, 0), std::invalid_argument);
}

TEST(Sensitivity, MorrisPlanIsOneAtATimeTrajectories) {
  const std::vector<Input> inputs = {Uniform("a", 0, 1), Uniform("b", 0, 1), Normal("c", 0, 1)};
  const int r = 5, k = 3;
  const Plan plan = MorrisPlan(inputs, r, 4, 7);
  ASSERT_EQ(plan.Runs(), static_cast<size_t>(r * (k + 1)));
  for (int t = 0; t < r; ++t) {
    const size_t start = static_cast<size_t>(t * (k + 1));
    EXPECT_EQ(plan.changed[start], -1);
    std::set<int> moved;
    for (size_t s = start + 1; s < start + k + 1; ++s) {
      const int i = plan.changed[s];
      moved.insert(i);
      for (int j = 0; j < k; ++j) {
        const double du = plan.u[s][j] - plan.u[s - 1][j];
        EXPECT_NEAR(std::fabs(du), j == i ? 0.5 : 0.0, 1e-15);
      }
    }
    EXPECT_EQ(moved.size(), static_cast<size_t>(k));  // each input moves exactly once
  }
  for (const auto& row : plan.u) {
    for (double u : row) {
      const double cell = u * 4 - 0.5;
      EXPECT_NEAR(cell, std::round(cell), 1e-12);  // a level centre
    }
  }
  EXPECT_EQ(MorrisPlan(inputs, r, 4, 7).values, plan.values);
  EXPECT_NE(MorrisPlan(inputs, r, 4, 8).values, plan.values);
}

TEST(Sensitivity, MonteCarloPlanIsALatinHypercube) {
  const std::vector<Input> inputs = {Uniform("a", 0, 1), Triangular("b", 0, 0.3, 1), Normal("c", 1, 0.1)};
  const int n = 50;
  const Plan plan = MonteCarloPlan(inputs, n, 3, 0);
  ASSERT_EQ(plan.Runs(), static_cast<size_t>(n));
  for (size_t j = 0; j < inputs.size(); ++j) {
    std::set<int> strata;
    for (const auto& row : plan.u) strata.insert(static_cast<int>(std::floor(row[j] * n)));
    EXPECT_EQ(strata.size(), static_cast<size_t>(n)) << "input " << j;
  }
  for (size_t r = 0; r < plan.Runs(); ++r) {
    EXPECT_DOUBLE_EQ(plan.values[r][1], inputs[1].Quantile(plan.u[r][1]));
  }
  EXPECT_EQ(MonteCarloPlan(inputs, n, 3, 0).values, plan.values);
  EXPECT_NE(MonteCarloPlan(inputs, n, 3, 1).values, plan.values);  // another batch
  EXPECT_NE(MonteCarloPlan(inputs, n, 4, 0).values, plan.values);  // another seed
}

TEST(Sensitivity, MorrisEffectsOfALinearFunctionAreExact) {
  // y = 3 x1 - 2 x2 + 0 x3: an effect per unit of probability is the slope
  // times the range width, the same at every point (sigma = 0).
  const std::vector<Input> inputs = {Uniform("x1", 0, 1), Uniform("x2", 0, 2), Uniform("x3", -1, 1)};
  const Plan plan = MorrisPlan(inputs, 8, 4, 1);
  const MorrisResult m = AnalyseMorris(plan, {"y"}, Evaluate(plan, Linear));
  ASSERT_EQ(m.outputs.size(), 1u);
  const auto& e = m.outputs[0].effects;
  EXPECT_NEAR(e[0].mu_star, 3.0, 1e-12);
  EXPECT_NEAR(e[0].mu, 3.0, 1e-12);
  EXPECT_NEAR(e[1].mu_star, 4.0, 1e-12);
  EXPECT_NEAR(e[1].mu, -4.0, 1e-12);
  EXPECT_NEAR(e[2].mu_star, 0.0, 1e-12);
  for (const auto& x : e) {
    EXPECT_EQ(x.effects, 8);
    EXPECT_NEAR(x.sigma, 0.0, 1e-12);
    EXPECT_LE(x.mu_star_interval.low, x.mu_star + 1e-12);
    EXPECT_GE(x.mu_star_interval.high, x.mu_star - 1e-12);
  }
  EXPECT_EQ(m.outputs[0].missing, 0);
}

TEST(Sensitivity, MorrisRanksTheGFunction) {
  std::vector<Input> inputs;
  for (int i = 0; i < 6; ++i) inputs.push_back(Uniform("x" + std::to_string(i + 1), 0, 1));
  const Plan plan = MorrisPlan(inputs, 40, 4, 11);
  const MorrisResult m = AnalyseMorris(plan, {"g"}, Evaluate(plan, GFunction));
  const auto& e = m.outputs[0].effects;
  EXPECT_GT(e[0].mu_star, e[1].mu_star);
  EXPECT_GT(e[1].mu_star, e[2].mu_star);
  EXPECT_GT(e[2].mu_star, e[3].mu_star);
  EXPECT_LT(e[4].mu_star, 0.05 * e[0].mu_star);
  EXPECT_LT(e[5].mu_star, 0.05 * e[0].mu_star);
  // The inputs interact through the product: sigma is not small.
  EXPECT_GT(e[0].sigma, 0.1 * e[0].mu_star);
}

TEST(Sensitivity, MorrisSkipsEffectsThatTouchAFailedRun) {
  const std::vector<Input> inputs = {Uniform("x1", 0, 1), Uniform("x2", 0, 2), Uniform("x3", -1, 1)};
  const Plan plan = MorrisPlan(inputs, 4, 4, 1);
  auto y = Evaluate(plan, Linear);
  y[5][0] = std::nan("");  // run 5 is the second run of trajectory 1
  const MorrisResult m = AnalyseMorris(plan, {"y"}, y);
  EXPECT_EQ(m.outputs[0].missing, 1);
  int effects = 0;
  for (const auto& e : m.outputs[0].effects) effects += e.effects;
  EXPECT_EQ(effects, 4 * 3 - 2);  // the effects into and out of run 5
  EXPECT_NEAR(m.outputs[0].effects[0].mu_star, 3.0, 1e-12);
  EXPECT_THROW(AnalyseMorris(plan, {"y", "z"}, y), std::invalid_argument);
}

TEST(Sensitivity, MonteCarloRecoversAKnownSpread) {
  // y = 2 x1 + x2 with x ~ N(0, 1): y ~ N(0, 5); four stacked batches.
  const std::vector<Input> inputs = {Normal("x1", 0, 1), Normal("x2", 0, 1)};
  std::vector<std::vector<double>> values, y;
  for (int b = 0; b < 4; ++b) {
    const Plan plan = MonteCarloPlan(inputs, 500, 5, b);
    for (const auto& row : plan.values) {
      values.push_back(row);
      y.push_back({TwoNormals(row), 1.0, std::nan("")});
    }
  }
  const MonteCarloResult mc = AnalyseMonteCarlo(inputs, values, {"y", "constant", "undefined"}, y, 5);
  ASSERT_EQ(mc.outputs.size(), 3u);
  const Spread& s = mc.outputs[0];
  const double sd = std::sqrt(5.0);
  EXPECT_EQ(s.valid, 2000);
  EXPECT_NEAR(s.mean, 0.0, 0.1);
  EXPECT_NEAR(s.sd, sd, 0.08);
  EXPECT_NEAR(s.p05, -1.6448536269514722 * sd, 0.2);
  EXPECT_NEAR(s.p95, 1.6448536269514722 * sd, 0.2);
  EXPECT_NEAR(s.p50, 0.0, 0.15);
  EXPECT_LT(s.p05_interval.low, s.p05);
  EXPECT_GT(s.p95_interval.high, s.p95);
  EXPECT_TRUE(s.converged);
  // Rank shares: about 0.8 and 0.2 of the spread (exactly so for the linear
  // regression; ranks bend it slightly), and the regression explains it all.
  EXPECT_NEAR(s.shares[0].share, 0.8, 0.06);
  EXPECT_NEAR(s.shares[1].share, 0.2, 0.06);
  EXPECT_GT(s.r2, 0.97);
  EXPECT_NEAR(s.shares[0].spearman, 6.0 / std::acos(-1.0) * std::asin(1.0 / std::sqrt(5.0)), 0.02);

  const Spread& c = mc.outputs[1];  // a constant output is converged and has no shares
  EXPECT_TRUE(c.converged);
  EXPECT_DOUBLE_EQ(c.p05, 1.0);
  EXPECT_TRUE(std::isnan(c.r2));

  const Spread& u = mc.outputs[2];  // never defined: reported, not counted for convergence
  EXPECT_EQ(u.valid, 0);
  EXPECT_EQ(u.missing, 2000);
  EXPECT_FALSE(u.converged);
  EXPECT_TRUE(mc.converged);
}

TEST(Sensitivity, MonteCarloIsNotConvergedWithFewRuns) {
  const std::vector<Input> inputs = {Normal("x1", 0, 1), Normal("x2", 0, 1)};
  const Plan plan = MonteCarloPlan(inputs, 12, 2, 0);
  const MonteCarloResult mc = AnalyseMonteCarlo(inputs, plan.values, {"y"}, Evaluate(plan, TwoNormals), 2);
  EXPECT_FALSE(mc.outputs[0].converged);
  EXPECT_FALSE(mc.converged);
  EXPECT_THROW(AnalyseMonteCarlo(inputs, plan.values, {"y"}, {}, 2), std::invalid_argument);
}

TEST(Sensitivity, PercentileMatchesNumpy) {
  EXPECT_DOUBLE_EQ(Percentile({4, 1, 3, 2}, 0.5), 2.5);
  EXPECT_DOUBLE_EQ(Percentile({4, 1, 3, 2}, 0.05), 1.15);
  EXPECT_DOUBLE_EQ(Percentile({4, 1, std::nan(""), 2}, 1.0), 4.0);
  EXPECT_TRUE(std::isnan(Percentile({}, 0.5)));
}

TEST(Sensitivity, Outputs) {
  const Output o = ParseOutput("CarbonDioxide.purity");
  EXPECT_EQ(o.component, "CarbonDioxide");
  EXPECT_EQ(o.kpi, "purity");
  EXPECT_EQ(o.Name(), "CarbonDioxide.purity");
  EXPECT_THROW(ParseOutput("CarbonDioxide"), std::invalid_argument);
  EXPECT_THROW(ParseOutput("CarbonDioxide.yield"), std::invalid_argument);
}

// ---- a real deck ------------------------------------------------------------------

namespace {

std::string ReadDeck(const std::string& name) {
  std::ifstream f("../examples/decks/" + name);
  std::stringstream s;
  s << f.rdbuf();
  return s.str();
}

std::string FirstStepTime(Simulation& sim) {
  for (Parameter* p : sim.Parameters().Tunable()) {
    if (p->Path().rfind("step[", 0) == 0 && p->Path().size() > 5 &&
        p->Path().compare(p->Path().size() - 5, 5, ".time") == 0) {
      return p->Path();
    }
  }
  return {};
}

}  // namespace

TEST(Sensitivity, ChecksAStudyAgainstTheDeck) {
  const std::string deck = ReadDeck("Arora2021_LiLSX_PSA.json");
  ASSERT_FALSE(deck.empty());
  Simulation sim = Simulation::FromJson(deck);
  const std::string time = FirstStepTime(sim);
  ASSERT_FALSE(time.empty());
  const std::vector<Output> outputs = DefaultOutputs(sim.Data());
  ASSERT_EQ(outputs.size(), 4u);  // Oxygen is tracked
  EXPECT_EQ(outputs[0].Name(), "Oxygen.purity");

  const double t = sim.GetParameter(time);
  const std::vector<Input> one = {Relative(time, t, 0.1)};
  // The Arora deck has no CSS check (css.method none): a study that requires
  // CSS must switch one on, or drop the requirement.
  EXPECT_THROW(CheckStudy(deck, one, outputs), std::invalid_argument);
  RunPolicy with_css;
  with_css.settings = {{"css.method", "3"}};
  EXPECT_NO_THROW(CheckStudy(deck, one, outputs, with_css));
  RunPolicy no_css;
  no_css.require_css = false;
  EXPECT_NO_THROW(CheckStudy(deck, one, outputs, no_css));

  EXPECT_THROW(CheckStudy(deck, {Uniform("step[999].time", 1, 2)}, outputs, no_css), std::invalid_argument);
  EXPECT_THROW(CheckStudy(deck, {Uniform("settings.reltol", 1e-6, 1e-4)}, outputs, no_css), std::invalid_argument);
  EXPECT_THROW(CheckStudy(deck, {Relative(time, t, 0.1), Relative(time, t, 0.2)}, outputs, no_css),
               std::invalid_argument);
  EXPECT_THROW(CheckStudy(deck, one, {{"Nitrogen", "purity"}}, no_css), std::invalid_argument);
  EXPECT_THROW(CheckStudy(deck, one, {}, no_css), std::invalid_argument);
  RunPolicy bad = no_css;
  bad.settings = {{"no.such.path", "1"}};
  EXPECT_THROW(CheckStudy(deck, one, outputs, bad), std::invalid_argument);
}

TEST(Sensitivity, EvaluatesSamplesOfARealDeck) {
  const std::string deck = ReadDeck("Arora2021_LiLSX_PSA.json");
  Simulation sim = Simulation::FromJson(deck);
  const std::string time = FirstStepTime(sim);
  const std::vector<Input> inputs = {Relative(time, sim.GetParameter(time), 0.2)};
  const std::vector<Output> outputs = {{"Oxygen", "purity"}, {"Oxygen", "recovery"}};

  // Two cycles: far from CSS, so the CSS requirement rejects the run...
  RunPolicy policy;
  policy.cycles = 2;
  policy.settings = {{"css.method", "3"}};
  const SampleOutcome strict = EvaluateSample(deck, inputs, {}, outputs, policy);
  EXPECT_TRUE(strict.status == RunStatus::Completed || strict.status == RunStatus::ReachedCSS);
  if (strict.status == RunStatus::Completed) {
    EXPECT_FALSE(strict.accepted);
    EXPECT_NE(strict.message.find("cyclic steady state"), std::string::npos);
    EXPECT_TRUE(std::isnan(strict.outputs[0]));
  }

  // ...and without it the KPIs of the last cycle come back, and a longer
  // feed step changes them.
  policy.require_css = false;
  const SampleOutcome base = EvaluateSample(deck, inputs, {}, outputs, policy);
  const SampleOutcome longer = EvaluateSample(deck, inputs, {inputs[0].upper}, outputs, policy);
  ASSERT_TRUE(base.accepted) << base.message;
  ASSERT_TRUE(longer.accepted) << longer.message;
  for (double v : base.outputs) EXPECT_TRUE(v > 0.0 && v <= 1.0) << v;
  EXPECT_NE(base.outputs[0], longer.outputs[0]);

  EXPECT_THROW(EvaluateSample(deck, inputs, {1.0, 2.0}, outputs, policy), std::invalid_argument);
  RunPolicy bad = policy;
  bad.settings = {{"no.such.path", "1"}};
  EXPECT_THROW(EvaluateSample(deck, inputs, {}, outputs, bad), std::invalid_argument);
}

}  // namespace sensitivity
}  // namespace adsim
