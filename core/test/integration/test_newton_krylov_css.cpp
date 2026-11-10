/*
 * Skarstrom - pressure and temperature swing adsorption simulator
 * Copyright (c) 2024-2026 Sam Affleck and contributors
 *
 * Licensed under the terms in the LICENSE file at the repository root.
 */

// Newton-Krylov cyclic steady state (AccelerationMethod::NEWTON_KRYLOV):
// GMRES, the damping rule and the fixed-point Newton solver on synthetic maps
// (fast), the accelerator on a coarse Perez 2019 A (fast), and the comparison
// with long plain runs on the shipped decks (NewtonKrylovCSSReference, slow
// tier through the "*CSSReference.*" filter of core/test/CMakeLists.txt).

#include <algorithm>
#include <cmath>
#include <cstdio>
#include <cstdlib>
#include <random>
#include <string>
#include <vector>

#include <gtest/gtest.h>

#include "adsim/simulation.h"
#include "adsim/solver/cyclic_steady_state.h"
#include "adsim/solver/newton_krylov.h"

namespace adsim {

  namespace {
    std::string Deck(const std::string& name) { return "../examples/decks/" + name + ".json"; }

    // The Newton-Krylov solves of a run, oldest first.
    std::vector<NewtonKrylovRecord> Solves(Simulation& sim) {
      std::vector<NewtonKrylovRecord> out;
      for (const auto& [cycle, record] : sim.Data().settings.acc_css.NewtonRecords()) out.push_back(record);
      return out;
    }

    void PrintSolves(const char* label, Simulation& sim) {
      for (const auto& r : Solves(sim)) {
        std::printf("  %s: solve after cycle %d: %d Newton steps, %d extra cycles, Krylov", label, r.cycle,
                    r.newton_iterations, r.cycles);
        for (int k : r.krylov_vectors) std::printf(" %d", k);
        std::printf(", residuals");
        for (double x : r.residuals) std::printf(" %.2e", x);
        std::printf(", multiplier %.4f (%s)\n", r.floquet_multiplier, r.outcome.c_str());
      }
    }

    int ExtraCycles(Simulation& sim) {
      int total = 0;
      for (const auto& r : Solves(sim)) total += r.cycles;
      return total;
    }

    const KPIRecord* Find(const std::vector<KPIRecord>& kpis, const std::string& component) {
      for (const auto& k : kpis)
        if (k.component == component) return &k;
      return nullptr;
    }

    // The state test's distance between two cycle states (the largest field
    // RMS difference, scaled as EndOfCycleStateChange scales it).
    double StateDistance(const CycleState& a, const CycleState& b) {
      auto is_loading = [](const CycleStateField& f) {
        return f.name.find("Adsorbed Concentration") != std::string::npos && f.name.find("Slow") == std::string::npos;
      };
      // RMS total loading over the cells (the loading scale), floored at
      // 1e-8 mol/kg.
      std::vector<double> total;
      for (const auto& f : a.fields) {
        if (!is_loading(f)) continue;
        if (total.size() < f.index.size()) total.resize(f.index.size(), 0.0);
        for (size_t j = 0; j < f.index.size(); ++j) total[j] += a.y[f.index[j]];
      }
      double loading_sq = 0.0;
      for (double q : total) loading_sq += q * q;
      const double loading_scale =
          std::max(total.empty() ? 0.0 : std::sqrt(loading_sq / static_cast<double>(total.size())), 1e-8);
      double worst = 0.0;
      for (size_t k = 0; k < a.fields.size() && k < b.fields.size(); ++k) {
        const auto& f = a.fields[k];
        const auto& g = b.fields[k];
        if (!f.differential || f.index.empty() || g.index.size() != f.index.size()) continue;
        double scale = 0.0;
        if (f.name.find("Mole Fraction") != std::string::npos) {
          scale = 1.0;
        } else if (f.name.find("Adsorbed Concentration") != std::string::npos) {
          scale = loading_scale;
        } else if (f.name.find("Temperature") != std::string::npos) {
          double sq = 0.0;
          for (size_t i : f.index) sq += a.y[i] * a.y[i];
          scale = std::sqrt(sq / static_cast<double>(f.index.size()));
        } else {
          continue;  // the total concentration: set by the pressure, as in the state test
        }
        double sq = 0.0;
        for (size_t j = 0; j < f.index.size(); ++j) {
          const double d = a.y[f.index[j]] - b.y[g.index[j]];
          sq += d * d;
        }
        worst = std::max(worst, std::sqrt(sq / static_cast<double>(f.index.size())) / scale);
      }
      return worst;
    }

    // Perez 2019 A on a coarse grid with short logs: a fast-tier deck that
    // still has the four-step cycle and its slow approach to CSS.
    Simulation CoarsePerezA(int cells = 20) {
      Simulation sim = Simulation::FromFile(Deck("Perez2019_A"));
      for (auto& layer : sim.Data().bed.layers_.layers_) layer.domain_.n_cells = cells;
      sim.Data().settings.max_cycles_to_keep_ = 4;
      return sim;
    }
  }

  // ---- GMRES, damping and the Floquet estimate ------------------------------

  TEST(NewtonKrylov, GmresSolvesASmallSystem) {
    std::mt19937 rng(7);
    std::normal_distribution<double> nd;
    const int n = 12;
    Eigen::MatrixXd A = Eigen::MatrixXd::Identity(n, n) * 3.0;
    for (int i = 0; i < n; ++i)
      for (int j = 0; j < n; ++j) A(i, j) += 0.3 * nd(rng);
    Eigen::VectorXd b(n);
    for (int i = 0; i < n; ++i) b[i] = nd(rng);
    int applications = 0;
    const GmresResult r = Gmres([&](const Eigen::VectorXd& v, Eigen::VectorXd& av) {
      ++applications;
      av = A * v;
      return true;
    }, b, 1e-12 * b.norm(), n + 5);
    EXPECT_TRUE(r.converged);
    EXPECT_LE(r.iterations, n);
    EXPECT_EQ(applications, r.iterations);
    EXPECT_LT((A * r.x - b).norm(), 1e-10 * b.norm());
    // The reported residual is the true one.
    EXPECT_NEAR(r.residual, (A * r.x - b).norm(), 1e-9 * b.norm());
  }

  TEST(NewtonKrylov, GmresStopsAtTheToleranceAndOnFailure) {
    const int n = 30;
    Eigen::VectorXd d(n);
    for (int i = 0; i < n; ++i) d[i] = 1.0 + i;  // condition number 30
    const Eigen::VectorXd b = Eigen::VectorXd::Ones(n);
    const GmresResult loose = Gmres([&](const Eigen::VectorXd& v, Eigen::VectorXd& av) {
      av = d.cwiseProduct(v);
      return true;
    }, b, 0.1 * b.norm(), n);
    EXPECT_TRUE(loose.converged);
    EXPECT_LT(loose.iterations, n);
    EXPECT_LE((d.cwiseProduct(loose.x) - b).norm(), 0.1 * b.norm() * (1 + 1e-12));

    int calls = 0;
    const GmresResult failed = Gmres([&](const Eigen::VectorXd& v, Eigen::VectorXd& av) {
      if (++calls > 3) return false;
      av = d.cwiseProduct(v);
      return true;
    }, b, 1e-14, n);
    EXPECT_TRUE(failed.failed);
    EXPECT_EQ(failed.iterations, 3);
    EXPECT_LT((d.cwiseProduct(failed.x) - b).norm(), b.norm());  // still the best of 3 vectors
  }

  TEST(NewtonKrylov, DominantMultiplierFromTheKrylovBasis) {
    // A = I - J with J having eigenvalues 0.97 and 0.9 and a cluster near 0.
    const int n = 40;
    std::mt19937 rng(3);
    std::normal_distribution<double> nd;
    Eigen::MatrixXd Q = Eigen::MatrixXd::NullaryExpr(n, n, [&]() { return nd(rng); });
    Q = Eigen::HouseholderQR<Eigen::MatrixXd>(Q).householderQ();
    Eigen::VectorXd mu = Eigen::VectorXd::Constant(n, 0.05);
    mu[0] = 0.97;
    mu[1] = 0.9;
    const Eigen::MatrixXd J = Q * mu.asDiagonal() * Q.transpose();
    const Eigen::MatrixXd A = Eigen::MatrixXd::Identity(n, n) - J;
    const GmresResult r = Gmres([&](const Eigen::VectorXd& v, Eigen::VectorXd& av) {
      av = A * v;
      return true;
    }, Eigen::VectorXd::Ones(n), 1e-10 * std::sqrt(double(n)), n);
    EXPECT_TRUE(r.converged);
    EXPECT_LE(r.iterations, 6);  // two slow modes and one cluster
    EXPECT_NEAR(DominantMultiplierFromHessenberg(r.hessenberg), 0.97, 1e-6);
    EXPECT_TRUE(std::isnan(DominantMultiplierFromHessenberg(Eigen::MatrixXd())));
  }

  TEST(NewtonKrylov, GloverDampingHalvesTheStepToZero) {
    EXPECT_DOUBLE_EQ(GloverDamping({1.0, 2.0}, {0.5, -1.0}), 1.0);       // nothing goes negative
    EXPECT_DOUBLE_EQ(GloverDamping({1.0, 2.0}, {-2.0, -1.0}), 0.25);     // x0 would reach -1: half of 1/2
    EXPECT_DOUBLE_EQ(GloverDamping({1.0, 2.0}, {-2.0, -8.0}), 0.125);    // the tightest entry wins
    EXPECT_DOUBLE_EQ(GloverDamping({0.0, 1.0}, {-1.0, 0.0}), 1.0);       // already at zero: ignored
  }

  // ---- the fixed-point solver on a synthetic cycle map ---------------------

  // G(x) = x* + J (x - x*) + c (x - x*)^2 (elementwise), with J symmetric,
  // multipliers 0.98 and 0.95 and the rest between 0 and 0.3, as the cycle
  // maps of the PSA decks: a few slow Floquet modes over a fast cluster.
  TEST(NewtonKrylov, SyntheticCycleMapConvergesSuperlinearly) {
    const int n = 60;
    std::mt19937 rng(11);
    std::normal_distribution<double> nd;
    std::uniform_real_distribution<double> ud(0.0, 0.3);
    Eigen::MatrixXd Q = Eigen::MatrixXd::NullaryExpr(n, n, [&]() { return nd(rng); });
    Q = Eigen::HouseholderQR<Eigen::MatrixXd>(Q).householderQ();
    Eigen::VectorXd mu(n);
    for (int i = 0; i < n; ++i) mu[i] = ud(rng);
    mu[0] = 0.98;
    mu[1] = 0.95;
    const Eigen::MatrixXd J = Q * mu.asDiagonal() * Q.transpose();
    std::vector<double> xs(n), x0(n);
    for (int i = 0; i < n; ++i) {
      xs[i] = 1.0 + 0.5 * std::sin(0.3 * i);
      x0[i] = xs[i] + 0.05 * nd(rng);
    }
    int evaluations = 0;
    auto map = [&](const std::vector<double>& x, std::vector<double>& gx) {
      ++evaluations;
      Eigen::VectorXd e(n);
      for (int i = 0; i < n; ++i) e[i] = x[i] - xs[i];
      const Eigen::VectorXd je = J * e;
      gx.resize(n);
      for (int i = 0; i < n; ++i) gx[i] = xs[i] + je[i] + 0.3 * e[i] * e[i];
      return true;
    };
    auto rms = [](const std::vector<double>& f) {
      double s = 0.0;
      for (double v : f) s += v * v;
      return std::sqrt(s / static_cast<double>(f.size()));
    };
    FixedPointNewtonOptions options;
    options.tolerance = 1e-11;
    options.fd_relative_step = 1e-7;  // an exact map: only rounding in the difference
    options.max_newton_iterations = 10;
    const FixedPointNewtonResult r = SolveFixedPointNewtonKrylov(
        map, x0, std::vector<double>(n, 1.0), rms, {}, options);
    ASSERT_TRUE(r.converged) << r.outcome;
    EXPECT_EQ(r.map_evaluations, evaluations);
    double error = 0.0;
    for (int i = 0; i < n; ++i) error = std::max(error, std::abs(r.x[i] - xs[i]));
    EXPECT_LT(error, 1e-9);
    EXPECT_NEAR(r.dominant_multiplier, 0.98, 0.01);
    // Superlinear: each of the last steps reduces the residual by more than
    // the one before (and by far more than plain iteration's 0.98).
    const auto& h = r.residual_history;
    ASSERT_GE(h.size(), 3u);
    for (size_t k = 1; k < h.size(); ++k) EXPECT_LT(h[k], 0.2 * h[k - 1]) << "step " << k;
    EXPECT_LT(h.back() / h[h.size() - 2], h[h.size() - 2] / h[h.size() - 3]);
    // Plain iteration would need log(1e-11 / 0.05) / log(0.98), over 1000 maps.
    EXPECT_LT(r.map_evaluations, 80);
    std::printf("  synthetic map: %d evaluations, Krylov", r.map_evaluations);
    for (int k : r.krylov_vectors) std::printf(" %d", k);
    std::printf(", residuals");
    for (double x : h) std::printf(" %.1e", x);
    std::printf("\n");
  }

  TEST(NewtonKrylov, FailedStepStopsTheSolve) {
    // A map whose iterates fail away from the start: the first Newton step
    // cannot be evaluated, so the solve stops with the start as its best.
    int calls = 0;
    auto map = [&](const std::vector<double>& x, std::vector<double>& gx) {
      ++calls;
      if (calls > 2 && std::abs(x[0] - 1.0) > 1e-3) return false;
      gx = {0.5 * x[0] + 1.0};
      return true;
    };
    auto norm = [](const std::vector<double>& f) { return std::abs(f[0]); };
    FixedPointNewtonOptions options;
    options.fd_relative_step = 1e-6;
    const FixedPointNewtonResult r = SolveFixedPointNewtonKrylov(map, {1.0}, {1.0}, norm, {}, options);
    EXPECT_FALSE(r.converged);
    EXPECT_EQ(r.newton_iterations, 0);
    EXPECT_DOUBLE_EQ(r.x[0], 1.0);
    EXPECT_EQ(r.outcome, "the Newton step did not reduce the residual");
  }

  // ---- the accelerator on a deck (fast tier) -------------------------------

  TEST(NewtonKrylovCSS, SelectableByName) {
    Simulation sim = Simulation::FromFile(Deck("Perez2019_A"));
    sim.SetParameterText("settings.acceleration_method", "newton_krylov");
    EXPECT_EQ(sim.Data().settings.acc_css.method_, AccelerationMethod::NEWTON_KRYLOV);
    // Stored by number: appended after Anderson, never renumbered.
    EXPECT_EQ(static_cast<int>(AccelerationMethod::NEWTON_KRYLOV), 3);
    EXPECT_STREQ(AccelerationMethodList[3], "Newton-Krylov");
  }

  TEST(NewtonKrylovCSS, CoarseDeckReachesCssInAFewNewtonSteps) {
    Simulation plain = CoarsePerezA();
    plain.Data().cycle.number_of_cycles_ = 150;
    const RunResult pr = plain.Run();
    ASSERT_TRUE(pr.css_reached) << pr.message;

    Simulation sim = CoarsePerezA();
    sim.Data().settings.acc_css.method_ = AccelerationMethod::NEWTON_KRYLOV;
    sim.Data().cycle.number_of_cycles_ = 150;
    int newton_lines = 0;
    RunOptions options;
    options.log_level = LogLevel::MEDIUM;
    options.on_message = [&newton_lines](LogLevel, const std::string& m) {
      if (m.find("Newton-Krylov: |F|") != std::string::npos) ++newton_lines;
    };
    const RunResult r = sim.Run(options);
    ASSERT_TRUE(r.Ok()) << r.message;
    ASSERT_TRUE(r.css_reached) << r.message;
    PrintSolves("coarse Perez 2019 A", sim);
    const auto solves = Solves(sim);
    ASSERT_FALSE(solves.empty());
    EXPECT_EQ(newton_lines, static_cast<int>(solves.size()));
    EXPECT_TRUE(std::any_of(solves.begin(), solves.end(), [](const auto& s) { return s.converged; }));
    const int cost = r.cycles_run + ExtraCycles(sim);
    std::printf("  plain stop %d; Newton-Krylov stop %d, %d cycle equivalents\n", pr.css_cycle, r.css_cycle, cost);
    EXPECT_LT(r.css_cycle, pr.css_cycle);

    // The run's CSS history carries the solve.
    bool found = false;
    for (const auto& rec : sim.CSSHistory()) {
      if (rec.newton_iterations < 0) continue;
      found = true;
      EXPECT_FALSE(rec.newton_residuals.empty());
      EXPECT_GT(rec.newton_cycles, 0);
    }
    EXPECT_TRUE(found);

    // Same KPIs as the plain run, within the plain stop's own distance from
    // its CSS (CSSRobustDeck holds a plain stop to 1e-4).
    const auto pk = plain.KPIs(-1), nk = sim.KPIs(-1);
    const KPIRecord* a = Find(pk, "CarbonDioxide");
    const KPIRecord* b = Find(nk, "CarbonDioxide");
    ASSERT_TRUE(a && b);
    EXPECT_NEAR(a->purity, b->purity, 1e-4);
    EXPECT_NEAR(a->recovery, b->recovery, 1e-4);
  }

  TEST(NewtonKrylovCSS, WarmStartNeedsOneOrTwoNewtonSteps) {
    // The optimiser's warm start: the CSS of one design starts a run of a
    // neighbouring design (feed step 5 % longer).
    Simulation base = CoarsePerezA();
    base.Data().settings.acc_css.method_ = AccelerationMethod::NEWTON_KRYLOV;
    base.Data().cycle.number_of_cycles_ = 150;
    ASSERT_TRUE(base.Run().css_reached);
    const CycleState css = base.FinalState();

    Simulation next = CoarsePerezA();
    next.Data().settings.acc_css.method_ = AccelerationMethod::NEWTON_KRYLOV;
    next.Data().cycle.number_of_cycles_ = 150;
    const double t = next.GetParameter("step[1].time");
    next.SetParameter("step[1].time", 1.05 * t);
    next.SetInitialState(css);
    const RunResult r = next.Run();
    ASSERT_TRUE(r.css_reached) << r.message;
    PrintSolves("warm start", next);
    const auto solves = Solves(next);
    ASSERT_FALSE(solves.empty());
    EXPECT_LE(solves.front().newton_iterations, 2);
    EXPECT_TRUE(solves.front().converged) << solves.front().outcome;
  }

  TEST(NewtonKrylovCSS, WithoutACycleMapItIsAnderson) {
    // An engine-level run (CoupledSolver directly, no Simulation) has no
    // cycle map: the accelerator falls back to Anderson.
    Simulation sim = CoarsePerezA(10);
    UserData& data = sim.Data();
    data.settings.acc_css.method_ = AccelerationMethod::NEWTON_KRYLOV;
    data.cycle.number_of_cycles_ = 40;
    std::string err;
    ASSERT_TRUE(data.EvaluateAllExpressions(err)) << err;
    data.settings.acc_css.ResetRun();
    CoupledSolver solver(data);
    int applied = 0;
    solver.SetLogLevel(LogLevel::MEDIUM);
    solver.SetMessageCallback([&applied](LogLevel, const std::string& m) {
      if (m.find("CSS acceleration applied") != std::string::npos) ++applied;
    });
    ASSERT_GE(solver.IntegrateCycle(data.cycle, false), 0);
    EXPECT_GT(applied, 0);
    EXPECT_TRUE(data.settings.acc_css.NewtonRecords().empty());
  }


  // ---- comparison with plain iteration (slow tier) --------------------------
  //
  // Each deck is run to CSS with Newton-Krylov from its own start, with its
  // own CSS rule. Then:
  //  * fixed point: plain cycles continued from the Newton-Krylov end state
  //    (CSS check off) must barely move it (state test norm, in multiples of
  //    the CSS tolerance eps);
  //  * reference (Perez 2019 A and E): as CSSReference builds its reference,
  //    the plain stop continued with the CSS check off for two segments of
  //    max(100, plain stop + 1) cycles. The drift between the two ends must
  //    be at most 1.5 eps (else the reference is not converged enough to judge
  //    by), and the Newton-Krylov end state and CO2 purity and recovery must
  //    be within c eps of the reference with the CSSReference bound c <= 100.
  // Nguyen 2023 X has no feasible plain reference: its slowest Floquet
  // multiplier is about 0.998, so a plain run still drifts after a thousand
  // cycles; it gets the fixed-point check only.
  namespace {
    constexpr double kEps = 2e-5;      // CyclicSteadyState::kDefaultStateTolerance on every deck
    constexpr double kBound = 100.0;   // test_css_reference.cpp's bound on c

    struct EndPoint {
      CycleState state;
      double purity = 0.0, recovery = 0.0;
    };

    EndPoint Capture(Simulation& sim) {
      EndPoint e;
      e.state = sim.FinalState();
      const auto kpis = sim.KPIs(-1);
      if (const KPIRecord* k = Find(kpis, "CarbonDioxide")) {
        e.purity = k->purity;
        e.recovery = k->recovery;
      }
      return e;
    }

    // Plain cycles from `from`, CSS check off.
    EndPoint Continue(const std::string& deck, const CycleState& from, int cycles) {
      Simulation sim = Simulation::FromFile(Deck(deck));
      sim.Data().settings.css.method_ = CSSMethod::NONE;
      sim.Data().settings.acc_css.method_ = AccelerationMethod::NONE;
      sim.Data().cycle.number_of_cycles_ = cycles;
      RunOptions options;
      options.initial_state = &from;
      const RunResult r = sim.Run(options);
      EXPECT_TRUE(r.Ok()) << r.message;
      return Capture(sim);
    }

    void CheckNewtonKrylov(const std::string& deck, bool with_reference, int max_cost, int fixed_point_cycles,
                           double max_fixed_point_move) {
      Simulation sim = Simulation::FromFile(Deck(deck));
      ASSERT_NE(sim.Data().settings.css.method_, CSSMethod::NONE);
      sim.Data().settings.acc_css.method_ = AccelerationMethod::NEWTON_KRYLOV;
      sim.Data().cycle.number_of_cycles_ = 1000;
      const RunResult r = sim.Run();
      ASSERT_TRUE(r.css_reached) << deck << ": " << r.message;
      PrintSolves(deck.c_str(), sim);
      const int cost = r.cycles_run + ExtraCycles(sim);
      const EndPoint nk = Capture(sim);

      const EndPoint after = Continue(deck, nk.state, fixed_point_cycles);
      const double c_fixed = StateDistance(nk.state, after.state) / kEps;
      std::printf("  %s: Newton-Krylov stop %d, %d cycle equivalents; %d plain cycles from its end move it by "
                  "%.3f eps (purity %.1e, recovery %.1e)\n",
                  deck.c_str(), r.css_cycle, cost, fixed_point_cycles, c_fixed, std::abs(after.purity - nk.purity),
                  std::abs(after.recovery - nk.recovery));
      EXPECT_LE(cost, max_cost) << deck;
      EXPECT_LE(c_fixed, max_fixed_point_move) << deck << ": the Newton-Krylov CSS is not a fixed point of plain cycles";

      if (!with_reference) return;
      Simulation plain = Simulation::FromFile(Deck(deck));
      plain.Data().cycle.number_of_cycles_ = 2000;
      const RunResult pr = plain.Run();
      ASSERT_TRUE(pr.css_reached) << deck << ": " << pr.message;
      const EndPoint plain_stop = Capture(plain);
      const int segment = std::max(100, pr.css_cycle + 1);
      const EndPoint mid = Continue(deck, plain_stop.state, segment);
      const EndPoint ref = Continue(deck, mid.state, segment);
      const double drift = StateDistance(mid.state, ref.state) / kEps;
      const double c_state = StateDistance(nk.state, ref.state) / kEps;
      const double c_kpi = std::max(std::abs(nk.purity - ref.purity), std::abs(nk.recovery - ref.recovery)) / kEps;
      const double c_plain_state = StateDistance(plain_stop.state, ref.state) / kEps;
      const double c_plain_kpi =
          std::max(std::abs(plain_stop.purity - ref.purity), std::abs(plain_stop.recovery - ref.recovery)) / kEps;
      std::printf("  %s: reference drift %.3f eps; c state %.3f, c KPIs %.3f (plain stop %d: %.2f, %.2f)\n",
                  deck.c_str(), drift, c_state, c_kpi, pr.css_cycle, c_plain_state, c_plain_kpi);
      EXPECT_LE(drift, 1.5) << deck;  // CSSReference's drift bound
      EXPECT_LE(c_state, kBound) << deck;
      EXPECT_LE(c_kpi, kBound) << deck;
    }
  }  // namespace

  // Measured (cycle equivalents; plain cycles' move from the end, c state and
  // c KPIs against the reference, in eps): Perez 2019 A 35, 0.43, 0.43, 0.38
  // (plain stop 4.4, 1.2); Perez 2019 E 48, 3.5, 3.7, 3.4 (plain stop 36,
  // 3.4); Nguyen 2023 X 218, 0.47. A solve stops at a residual r and leaves
  // about r / (1 - mu) along the slowest mode, hence E's larger move. The bounds on cost and move are regression bounds from
  // these measurements; c <= 100 and the drift <= 1.5 eps are CSSReference's.
  TEST(NewtonKrylovCSSReference, Perez2019A) { CheckNewtonKrylov("Perez2019_A", true, 60, 50, 1.0); }
  TEST(NewtonKrylovCSSReference, Perez2019E) { CheckNewtonKrylov("Perez2019_E", true, 80, 50, 10.0); }
  TEST(NewtonKrylovCSSReference, Nguyen2023ExpX) { CheckNewtonKrylov("Nguyen2023_Exp_X", false, 300, 50, 2.0); }

  // Below the integrator's noise floor at the deck's tolerance: with reltol
  // 1e-12, one solve takes the residual from 5e-3 to the floor the absolute
  // tolerances set (about 1e-8 in the state norm: the loadings' 1e-8 mol/kg
  // and the temperatures' 1e-4 K are fixed), three orders below the stopping
  // rule's tolerance eps = 2e-5, with factors that shrink while clear of the
  // floor. Measured first solves (57 to 62 extra cycles, ~45 s on Linux):
  //   Linux (GCC):          4.97e-3, 5.01e-4, 3.20e-5, 1.71e-6, 1.81e-8,
  //                         1.22e-8, 6.03e-9 (converged)
  //   macOS (Apple clang):  4.97e-3, 5.00e-4, 3.19e-5, 1.10e-6, 9.36e-8,
  //                         4.80e-8, 1.59e-8, 1.20e-8, 7.59e-9 (converged)
  //   Windows (MSVC):       4.97e-3, 5.16e-4, 3.70e-5, 1.35e-6, 1.02e-7,
  //                         1.02e-8, then a step that does not reduce it
  // Where the floor lies differs with the platform's libm and FMA use: the
  // Newton tolerance of 1e-8 sits on it, so whether the solve reports
  // "converged" (Linux, macOS) or stops one step short with "the Newton step
  // did not reduce the residual" at 1.02e-8 (Windows) is rounding, not
  // Newton-Krylov. So the test asks for either ending and for a best
  // residual within kFloorMargin of the floor (3x the 1e-8 estimate, about 3x the
  // largest measured 1.02e-8, still ~700x below eps; plain iteration, contracting by
  // a constant factor near 1, would need hundreds of cycles to get there).
  //
  // The contraction factors: only steps that start AND end above kFloorBand
  // (30x the floor) are measured. A step that ends nearer the floor carries
  // its noise (and its finite-difference J v products are noise-dominated
  // there): the step to the floor gives 0.011 (Linux), 0.085 (macOS); the
  // Windows step to 1.02e-7, 10x the floor, gives 0.076. Clear of the band
  // the three platforms give 0.101, 0.064, 0.053 (Linux), 0.101, 0.064,
  // 0.035 (macOS) and 0.104, 0.072, 0.036 (Windows): each well below plain
  // iteration's and the last at most three quarters of the first (0.53,
  // 0.34, 0.35), so the test still fails if the solve stops contracting
  // superlinearly (a constant factor) or stalls.
  TEST(NewtonKrylovCSSReference, ConvergesBelowTheNoiseFloor) {
    Simulation sim = CoarsePerezA();
    sim.Data().settings.reltol = 1e-12;
    sim.Data().settings.acc_css.method_ = AccelerationMethod::NEWTON_KRYLOV;
    sim.Data().settings.acc_css.newton_tolerance_ = 1e-8;
    sim.Data().cycle.number_of_cycles_ = 60;
    const RunResult r = sim.Run();
    ASSERT_TRUE(r.css_reached) << r.message;
    PrintSolves("coarse Perez 2019 A, reltol 1e-12", sim);
    const auto solves = Solves(sim);
    ASSERT_FALSE(solves.empty());
    const auto& first = solves.front();
    const auto& h = first.residuals;
    ASSERT_FALSE(h.empty());
    // Ended at the floor: converged there, or the next step could not reduce
    // the residual further (a Krylov failure or the iteration limit is not
    // accepted).
    EXPECT_TRUE(first.converged || first.outcome == "the Newton step did not reduce the residual") << first.outcome;
    constexpr double kFloor = 1e-8;  // the absolute tolerances' floor in the state norm
    constexpr double kFloorMargin = 3.0;
    const double best = *std::min_element(h.begin(), h.end());
    std::printf("  best residual %.3e (%.2f x the 1e-8 floor, %.1e eps)\n", best, best / kFloor, best / kEps);
    EXPECT_LE(best, kFloorMargin * kFloor);
    constexpr double kFloorBand = 30.0 * kFloor;
    std::vector<double> factors;
    for (size_t k = 1; k < h.size(); ++k)
      if (h[k - 1] > kFloorBand && h[k] > kFloorBand) factors.push_back(h[k] / h[k - 1]);
    std::string sequence;
    for (double q : factors) sequence += " " + std::to_string(q);
    std::printf("  factors above the floor band:%s\n", sequence.c_str());
    ASSERT_GE(factors.size(), 3u);
    for (double q : factors) EXPECT_LT(q, 0.25);
    EXPECT_LT(factors.back(), 0.75 * factors.front());
  }

  // ---- measurement helper (not run by default) ------------------------------
  // ADSIM_NK_DECK=<deck> [ADSIM_NK_RELTOL=r] [ADSIM_NK_METHOD=0..3] [ADSIM_NK_CELLS=n]
  // [ADSIM_NK_TOL=newton tolerance] [ADSIM_NK_FACTOR=tightened reltol factor] [ADSIM_NK_VERBOSE=1]
  //   ./AdsorptionSimulatorTests --gtest_also_run_disabled_tests
  //   --gtest_filter=NewtonKrylovCSS.DISABLED_Measure
  TEST(NewtonKrylovCSS, DISABLED_Measure) {
    const char* deck = std::getenv("ADSIM_NK_DECK");
    if (!deck) GTEST_SKIP() << "set ADSIM_NK_DECK";
    Simulation sim = Simulation::FromFile(Deck(deck));
    const char* method = std::getenv("ADSIM_NK_METHOD");
    sim.Data().settings.acc_css.method_ = static_cast<AccelerationMethod>(method ? std::atoi(method) : 3);
    if (const char* rt = std::getenv("ADSIM_NK_RELTOL")) sim.Data().settings.reltol = std::atof(rt);
    if (const char* f = std::getenv("ADSIM_NK_FACTOR")) sim.Data().settings.acc_css.newton_reltol_factor_ = std::atof(f);
    if (const char* tol = std::getenv("ADSIM_NK_TOL")) sim.Data().settings.acc_css.newton_tolerance_ = std::atof(tol);
    if (const char* cells = std::getenv("ADSIM_NK_CELLS"))
      for (auto& layer : sim.Data().bed.layers_.layers_) layer.domain_.n_cells = std::atoi(cells);
    sim.Data().cycle.number_of_cycles_ = 1000;
    RunOptions options;
    options.log_level = LogLevel::MEDIUM;
    if (std::getenv("ADSIM_NK_VERBOSE"))
      options.on_message = [](LogLevel, const std::string& m) { std::printf("%s\n", m.c_str()); };
    const RunResult r = sim.Run(options);
    PrintSolves(deck, sim);
    std::printf("  %s: %s, stop %d, cycles run %d, extra %d, wall %.1f s\n", deck, r.message.c_str(), r.css_cycle,
                r.cycles_run, ExtraCycles(sim), r.wall_time_s);
    for (const auto& k : sim.KPIs(-1))
      std::printf("  %s purity %.8f recovery %.8f\n", k.component.c_str(), k.purity, k.recovery);
  }

}  // namespace adsim
