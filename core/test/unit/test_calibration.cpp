// test_calibration.cpp
//
// The calibration object (adsim/solver/calibration.h) on the diluted-bed
// breakthrough deck (CO2 into a He-filled bed, one 3000 s step) cut to 12
// cells: a run takes well under a second. The data are synthetic: the deck
// is run at known parameter values, its outlet CO2 is read at the measured
// times and a seeded Gaussian noise of known standard deviation is added, so
// the fit must recover the known values within their standard errors. The
// tests check the solver and the bookkeeping, not the deck's results.

#include <atomic>
#include <cmath>
#include <cstdint>
#include <filesystem>
#include <random>
#include <string>
#include <vector>

#include <gtest/gtest.h>

#include "adsim/simulation.h"
#include "adsim/solver/calibration.h"

namespace adsim {
namespace {

namespace cal = calibration;

const std::string kDeck =
    (std::filesystem::current_path() / ".." / "examples" / "decks" / "Illustrative_Diluted_Bed.json").string();
constexpr int kCells = 12;
const std::string kOutlet = "outlet_CO2_mole_fraction[-]_Qin=100[ccm]";
const std::string kCapacity = "adsorbent[1].isotherm[CarbonDioxide].m2";       // 6.04 mol/kg in the deck
const std::string kAffinity = "adsorbent[1].isotherm[CarbonDioxide].b2";      // 3.32e-6 m3/mol
const std::string kFlow = "endpoint[1].volumetric_flow";                       // 1.667e-6 m3/s
constexpr double kNoise = 0.01;                                                // [mol/mol]

// A deterministic standard normal draw (Box-Muller on std::mt19937_64, whose
// sequence the standard fixes; std::normal_distribution is not portable).
double Gaussian(std::mt19937_64& rng) {
  auto u = [&]() { return (static_cast<double>(rng() >> 11) + 0.5) * 0x1.0p-53; };
  const double u1 = u(), u2 = u();
  return std::sqrt(-2.0 * std::log(u1)) * std::cos(6.283185307179586 * u2);
}

// The deck on kCells cells with `values` set, its outlet CO2 dataset replaced
// by the simulated curve at the measured times plus noise (with that noise as
// its uncertainty), as a case fitting only that dataset.
cal::Case SyntheticCase(const std::vector<std::pair<std::string, double>>& values, std::uint64_t seed,
                        const std::string& name) {
  Simulation sim = Simulation::FromFile(kDeck);
  for (auto& layer : sim.Data().bed.layers_.layers_) layer.domain_.n_cells = kCells;
  for (const auto& [path, v] : values) sim.SetParameter(path, v);
  const RunResult r = sim.Run();
  EXPECT_TRUE(r.Ok()) << r.message;
  XYData* data = sim.Data().exp_data_.Find(kOutlet);
  EXPECT_NE(data, nullptr);
  const SimulatedCounterpart s =
      SimulatedSeries(*data, sim.Log(), sim.KPIs(), BedBoundarySources{&sim.Boundaries(), &sim.Data().cycle.step_list_});
  std::mt19937_64 rng(seed);
  data->observable.uncertainty.clear();
  for (size_t i = 0; i < data->x.size(); ++i) {
    data->y[i] = InterpLinear(s.x, s.y, data->x[i]) + kNoise * Gaussian(rng);
    data->observable.uncertainty.push_back(kNoise);
  }
  cal::Case c = cal::Case::FromData(sim.Data(), {cal::Dataset(kOutlet)});
  c.name = name;
  return c;
}

cal::Options Quick() {
  cal::Options o;
  o.threads = 3;
  o.max_simulations = 120;
  return o;
}

// ---- validation (no runs) ------------------------------------------------------

TEST(Calibration, ProblemValidatesAgainstEveryDeck) {
  cal::Case a = cal::Case::FromFile(kDeck);
  cal::Case b = cal::Case::FromFile(kDeck);
  b.name = "second";
  // Every bound dataset by default (the deck binds three).
  cal::Problem p({a, b}, {cal::Shared(kCapacity, 4.0, 8.0), cal::PerCase(kAffinity, 1e-6, 2e-5, cal::Scale::Log)});
  EXPECT_EQ(p.Datasets().size(), 6u);
  ASSERT_EQ(p.Unknowns().size(), 3u);  // one shared, two per-case
  EXPECT_EQ(p.Unknowns()[0].case_index, -1);
  EXPECT_DOUBLE_EQ(p.Unknowns()[0].initial, 6.04);  // the deck's value
  EXPECT_EQ(p.Unknowns()[2].label, kAffinity + " [second]");
  EXPECT_EQ(p.UnknownIndex(1, 1), 2u);

  auto message = [](auto&& make) {
    try {
      make();
    } catch (const std::invalid_argument& e) {
      return std::string(e.what());
    }
    return std::string("no error");
  };
  EXPECT_NE(message([&] { cal::Problem({a, b}, {cal::Shared("layer[1].no_such_thing", 0, 1)}); }).find("case"),
            std::string::npos);
  EXPECT_NE(message([&] { cal::Problem({a}, {cal::Shared(kCapacity, 0.0, 8.0, cal::Scale::Log)}); }).find("log"),
            std::string::npos);
  EXPECT_NE(message([&] { cal::Problem({a}, {cal::Shared(kCapacity, 8.0, 4.0)}); }).find("bounds"),
            std::string::npos);
  EXPECT_NE(message([&] { cal::Problem({a}, {cal::Shared(kCapacity, 4.0, 8.0, cal::Scale::Linear, 9.0)}); })
                .find("outside"),
            std::string::npos);
  EXPECT_NE(message([&] { cal::Problem({a}, {cal::Shared(kCapacity, 4, 8), cal::Shared(kCapacity, 4, 8)}); })
                .find("twice"),
            std::string::npos);
  cal::Case missing = a;
  missing.datasets = {cal::Dataset("no such dataset")};
  EXPECT_NE(message([&] { cal::Problem({missing}, {cal::Shared(kCapacity, 4, 8)}); }).find("no dataset"),
            std::string::npos);
  cal::Case bad_weight = a;
  bad_weight.datasets = {cal::Dataset(0, -1.0)};
  EXPECT_NE(message([&] { cal::Problem({bad_weight}, {cal::Shared(kCapacity, 4, 8)}); }).find("weight"),
            std::string::npos);
  cal::Parameter wrong = cal::PerCase(kAffinity, 1e-6, 2e-5);
  wrong.case_bounds = {cal::Bounds{2e-6, 1e-5}};  // one entry, two cases
  EXPECT_NE(message([&] { cal::Problem({a, b}, {wrong}); }).find("per-case bounds"), std::string::npos);

  // Unbound data cannot be fitted.
  Simulation sim = Simulation::FromFile(kDeck);
  for (auto& d : sim.Data().exp_data_.GetAllData()) d.observable = ObservableBinding{};
  cal::Case unbound = cal::Case::FromData(sim.Data());
  EXPECT_NE(message([&] { cal::Problem({unbound}, {cal::Shared(kCapacity, 4, 8)}); }).find("bound"),
            std::string::npos);
}

TEST(Calibration, EvaluateGivesTheDocumentedObjective) {
  // At the values the data were made from, each residual is noise / sigma
  // scaled by sqrt(weight / n): S is near weight * mean((noise / sigma)^2) ~ 1.
  cal::Case c = SyntheticCase({{kCapacity, 6.04}}, 7, "truth");
  c.datasets = {cal::Dataset(kOutlet, 2.0)};
  cal::Problem p({c}, {cal::Shared(kCapacity, 4.0, 8.0)});
  const cal::Evaluation e = cal::Evaluate(p, {6.04});
  ASSERT_TRUE(e.ok) << e.message;
  const XYData& d = p.Datasets()[0].data;
  EXPECT_GE(e.residuals.size() + 1, d.x.size());  // the simulated step covers the measured times
  EXPECT_NEAR(e.objective, 2.0, 1.2);
  EXPECT_THROW(cal::Evaluate(p, {9.0}), std::invalid_argument);
}

// ---- (a) recovering known values -------------------------------------------------

TEST(Calibration, LeastSquaresRecoversKnownValues) {
  const double m2 = 5.5, b2 = 4.5e-6;
  cal::Case c = SyntheticCase({{kCapacity, m2}, {kAffinity, b2}}, 11, "synthetic");
  cal::Problem p({c}, {cal::Shared(kCapacity, 3.0, 9.0, cal::Scale::Linear, 6.6),
                       cal::Shared(kAffinity, 1e-6, 2e-5, cal::Scale::Log, 3.32e-6)});
  std::vector<cal::Progress> progress;
  cal::Options o = Quick();
  // The two isotherm parameters are strongly correlated (the outlet curve
  // mostly fixes their product): central differences follow the valley.
  o.difference = cal::Difference::Central;
  o.max_simulations = 200;
  o.progress = [&](const cal::Progress& pr) { progress.push_back(pr); };
  const cal::Result r = cal::Calibrate(p, o);
  ASSERT_TRUE(r.Ok()) << r.message;
  EXPECT_NE(r.stop, cal::StopReason::Budget) << r.message;
  EXPECT_EQ(r.failed_simulations, 0);
  ASSERT_EQ(r.parameters.size(), 2u);
  const double truth[2] = {m2, b2};
  for (size_t i = 0; i < 2; ++i) {
    const auto& e = r.parameters[i];
    SCOPED_TRACE(e.label);
    ASSERT_TRUE(std::isfinite(e.standard_error));
    EXPECT_GT(e.standard_error, 0.0);
    EXPECT_LT(e.standard_error, 0.2 * e.estimate);  // the data determine both
    EXPECT_LT(std::abs(e.estimate - truth[i]), 3.0 * e.standard_error);  // ~2 se, with margin for one seed
    EXPECT_FALSE(e.at_lower || e.at_upper);
  }
  EXPECT_EQ(r.parameters[1].scale, cal::Scale::Log);
  // sigma^2 (N - p) = S, and S ~ 1 when the uncertainty matches the noise:
  // the scale of the standard errors is right.
  EXPECT_NEAR(r.sigma2 * static_cast<double>(r.residuals - 2), r.objective, 1e-12);
  EXPECT_GT(r.objective, 0.3);
  EXPECT_LT(r.objective, 2.0);
  // The truth is worse than the estimate by about chi^2(2) / n in S (n the
  // points: each r_k carries 1 / sqrt(n)), not by orders.
  const cal::Evaluation at_truth = cal::Evaluate(p, {m2, b2});
  ASSERT_TRUE(at_truth.ok);
  const double n = static_cast<double>(r.residuals);
  EXPECT_GE(at_truth.objective, r.objective - 1e-9);
  EXPECT_LT(n * (at_truth.objective - r.objective), 15.0);
  // Correlations: unit diagonal, symmetric, inside [-1, 1].
  ASSERT_EQ(r.correlation.size(), 2u);
  EXPECT_NEAR(r.correlation[0][0], 1.0, 1e-12);
  EXPECT_NEAR(r.correlation[0][1], r.correlation[1][0], 1e-12);
  EXPECT_LE(std::abs(r.correlation[0][1]), 1.0);
  // Per dataset: the noise level comes back.
  ASSERT_EQ(r.datasets.size(), 1u);
  const auto& f = r.datasets[0];
  EXPECT_EQ(f.points, p.Datasets()[0].data.x.size());
  EXPECT_EQ(f.skipped, 0u);
  EXPECT_GT(f.rmse, 0.4 * kNoise);
  EXPECT_LT(f.rmse, 1.6 * kNoise);
  EXPECT_NEAR(f.objective, r.objective, 1e-12);
  EXPECT_EQ(f.residuals.size(), f.points);
  EXPECT_FALSE(f.series.x.empty());
  // The trace starts at the start and falls; progress was reported for it.
  ASSERT_GE(r.trace.size(), 2u);
  EXPECT_NEAR(r.trace.front().values[0], 6.6, 1e-12);
  for (size_t i = 1; i < r.trace.size(); ++i) EXPECT_LE(r.trace[i].objective, r.trace[i - 1].objective);
  EXPECT_EQ(progress.size(), r.trace.size());
  EXPECT_GE(r.simulations, r.trace.back().simulations);  // + a final Jacobian after a budget stop
  // The calibrated deck carries the estimates.
  ASSERT_EQ(r.calibrated_decks.size(), 1u);
  Simulation fitted = Simulation::FromJson(r.calibrated_decks[0]);
  EXPECT_DOUBLE_EQ(fitted.GetParameter(kCapacity), r.parameters[0].estimate);
  EXPECT_DOUBLE_EQ(fitted.GetParameter(kAffinity), r.parameters[1].estimate);
}

TEST(Calibration, TheSameInputsGiveTheSameResult) {
  cal::Case c = SyntheticCase({{kCapacity, 5.5}}, 3, "synthetic");
  cal::Problem p({c}, {cal::Shared(kCapacity, 3.0, 9.0, cal::Scale::Linear, 6.6)});
  cal::Options o = Quick();
  o.max_iterations = 3;
  o.difference = cal::Difference::Central;
  o.threads = 1;
  const cal::Result one = cal::Calibrate(p, o);
  o.threads = 4;
  const cal::Result four = cal::Calibrate(p, o);
  ASSERT_TRUE(one.Ok());
  EXPECT_EQ(one.parameters[0].estimate, four.parameters[0].estimate);
  EXPECT_EQ(one.parameters[0].standard_error, four.parameters[0].standard_error);
  EXPECT_EQ(one.objective, four.objective);
  EXPECT_EQ(one.simulations, four.simulations);
}

// ---- (b) shared and per-case parameters --------------------------------------------

TEST(Calibration, SharedAndPerCaseParameters) {
  // Two experiments at different flows: one capacity, an affinity each (synthetic).
  const double m2 = 5.6, b2_a = 2.5e-6, b2_b = 5.0e-6, flow_b = 2.5e-6;
  cal::Case a = SyntheticCase({{kCapacity, m2}, {kAffinity, b2_a}}, 21, "low flow");
  cal::Case b = SyntheticCase({{kFlow, flow_b}, {kCapacity, m2}, {kAffinity, b2_b}}, 22, "high flow");
  cal::Parameter b2 = cal::PerCase(kAffinity, 1e-6, 2e-5, cal::Scale::Log);
  b2.case_initial = {3.32e-6, 3.32e-6};
  cal::Problem p({a, b}, {cal::Shared(kCapacity, 3.0, 9.0, cal::Scale::Linear, 6.4), b2});
  ASSERT_EQ(p.Unknowns().size(), 3u);
  const cal::Result r = cal::Calibrate(p, Quick());
  ASSERT_TRUE(r.Ok()) << r.message;
  const double truth[3] = {m2, b2_a, b2_b};
  for (size_t i = 0; i < 3; ++i) {
    const auto& e = r.parameters[i];
    SCOPED_TRACE(e.label);
    ASSERT_TRUE(std::isfinite(e.standard_error));
    EXPECT_LT(std::abs(e.estimate - truth[i]), std::max(3.0 * e.standard_error, 0.02 * truth[i]));
  }
  EXPECT_EQ(r.parameters[1].case_index, 0);
  EXPECT_EQ(r.parameters[2].case_index, 1);
  EXPECT_EQ(r.datasets.size(), 2u);
  EXPECT_EQ(r.datasets[1].case_name, "high flow");
  // Each calibrated deck has the shared value and its own affinity, and
  // keeps its own flow.
  ASSERT_EQ(r.calibrated_decks.size(), 2u);
  Simulation da = Simulation::FromJson(r.calibrated_decks[0]);
  Simulation db = Simulation::FromJson(r.calibrated_decks[1]);
  EXPECT_DOUBLE_EQ(da.GetParameter(kCapacity), db.GetParameter(kCapacity));
  EXPECT_DOUBLE_EQ(da.GetParameter(kAffinity), r.parameters[1].estimate);
  EXPECT_DOUBLE_EQ(db.GetParameter(kAffinity), r.parameters[2].estimate);
  EXPECT_DOUBLE_EQ(db.GetParameter(kFlow), flow_b);
}

// ---- (c) the global solver ------------------------------------------------------------

TEST(Calibration, GlobalSolverWithATinyBudget) {
  cal::Case c = SyntheticCase({{kCapacity, 5.5}, {kAffinity, 4.5e-6}}, 11, "synthetic");
  cal::Problem p({c}, {cal::Shared(kCapacity, 3.0, 9.0), cal::Shared(kAffinity, 1e-6, 2e-5, cal::Scale::Log)});
  cal::Options o = Quick();
  o.solver = cal::Solver::Global;
  o.global_evaluations = 5;
  const cal::Result r = cal::Calibrate(p, o);
  ASSERT_TRUE(r.Ok()) << r.message;
  EXPECT_EQ(r.solver, cal::Solver::Global);
  EXPECT_FALSE(r.refined);
  EXPECT_EQ(r.stop, cal::StopReason::GlobalDone);
  EXPECT_EQ(r.trace.size(), 5u);
  EXPECT_TRUE(std::isfinite(r.objective));
  EXPECT_EQ(r.objective, r.trace.back().objective);
  // Five candidates plus the Jacobian for the standard errors.
  EXPECT_EQ(r.simulations, 5 + 2);
  for (const auto& e : r.parameters) {
    EXPECT_GE(e.estimate, e.lower);
    EXPECT_LE(e.estimate, e.upper);
  }

  // Refined by least squares: lower still.
  o.refine = true;
  const cal::Result refined = cal::Calibrate(p, o);
  ASSERT_TRUE(refined.Ok()) << refined.message;
  EXPECT_TRUE(refined.refined);
  EXPECT_LT(refined.objective, r.objective);
  EXPECT_EQ(refined.trace.front().solver, cal::Solver::Global);
  EXPECT_EQ(refined.trace.back().solver, cal::Solver::LeastSquares);
}

// ---- (d) cancelling -------------------------------------------------------------------

TEST(Calibration, CancelReturnsAPartialResult) {
  cal::Case c = SyntheticCase({{kCapacity, 5.5}, {kAffinity, 4.5e-6}}, 11, "synthetic");
  cal::Problem p({c}, {cal::Shared(kCapacity, 3.0, 9.0, cal::Scale::Linear, 6.6),
                       cal::Shared(kAffinity, 1e-6, 2e-5, cal::Scale::Log, 3.32e-6)});
  std::atomic<bool> stop{false};
  cal::Options o = Quick();
  o.cancel = [&]() { return stop.load(); };
  o.progress = [&](const cal::Progress& pr) {
    if (pr.iteration >= 1) stop = true;  // after the first step
  };
  const cal::Result r = cal::Calibrate(p, o);
  EXPECT_TRUE(r.cancelled);
  EXPECT_EQ(r.stop, cal::StopReason::Cancelled);
  EXPECT_TRUE(r.Ok());
  EXPECT_LE(r.iterations, 2);
  ASSERT_EQ(r.parameters.size(), 2u);
  EXPECT_TRUE(std::isfinite(r.parameters[0].estimate));
  EXPECT_TRUE(std::isfinite(r.objective));
  EXPECT_EQ(r.datasets.size(), 1u);

  // Cancelled before anything ran: no estimate beyond the start.
  cal::Options now = Quick();
  now.cancel = [] { return true; };
  const cal::Result none = cal::Calibrate(p, now);
  EXPECT_TRUE(none.cancelled);
  EXPECT_NEAR(none.parameters[0].estimate, 6.6, 1e-12);
  EXPECT_FALSE(none.Ok());
}

// ---- failed runs ----------------------------------------------------------------------

TEST(Calibration, AFailedStartIsReportedNotThrown) {
  cal::Case c = SyntheticCase({{kCapacity, 5.5}}, 3, "synthetic");
  c.policy.require_css = true;  // a breakthrough never reaches CSS: every run is rejected
  cal::Problem p({c}, {cal::Shared(kCapacity, 3.0, 9.0)});
  const cal::Result r = cal::Calibrate(p, Quick());
  EXPECT_FALSE(r.Ok());
  EXPECT_EQ(r.stop, cal::StopReason::StartFailed);
  EXPECT_NE(r.message.find("cyclic steady state"), std::string::npos) << r.message;
  EXPECT_EQ(r.failed_simulations, 1);
  EXPECT_TRUE(std::isnan(r.objective));
}

}  // namespace
}  // namespace adsim
