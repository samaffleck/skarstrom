// test_cycle_sensitivity.cpp
//
// The forward-sensitivity prototype (adsim/solver/cycle_sensitivity.h): one
// cycle with IDAS forward sensitivities along directions v in the
// differential initial state gives the cycle map's Jacobian-vector product
// Phi_x v. Checked against a centred finite difference of the plain cycle,
// (G(x + h v) - G(x - h v)) / 2h, whose own error is estimated from the steps
// h and 2h.
//
// Where it works: a cycle whose residual stays smooth along the run. The
// single-step Casas 2012 breakthrough deck, started after one warm-up cycle
// (so the measured cycle starts from a flowing bed, not from a bed at rest),
// with temperature-bump directions. Where it does not: across the step
// switches of a PSA cycle (Perez 2019 A), where the velocities cross zero and
// the residual is not differentiable; see docs/cyclic-steady-state.md,
// "Cycle-map Jacobian by forward sensitivities".
//
//  * CycleSensitivity.MatchesCentredDifferenceOnASmoothCycle (fast tier):
//    Casas 2012 on 10 cells, one temperature bump, both correctors.
//  * CycleSensitivity.RejectsConnectedDecks (fast tier).
//  * CycleSensitivity.MeasureCostAndAccuracy (slow tier): the cost per
//    direction for 1, 4 and 16 directions with both correctors against a
//    plain cycle and the accuracy against a centred difference at reltol
//    1e-10, on Casas 2012 at its own grid; then what happens on Perez 2019 A.
//    Prints the tables the documentation reports. Runs that fail are part
//    of the measurement (printed, not asserted); a run that completes must
//    agree with a reliable difference.

#include <algorithm>
#include <cmath>
#include <cstdio>
#include <random>
#include <string>
#include <vector>

#include <gtest/gtest.h>

#include "adsim/objects/user_data.h"
#include "adsim/solver/cycle_sensitivity.h"

namespace adsim {

  namespace {

    using Vec = std::vector<double>;
    using Options = CycleSensitivityOptions;

    constexpr const char* kSmoothDeck = "Casas2012_Fig4_50_50_15bar_25C.json";
    constexpr const char* kPsaDeck = "Perez2019_A.json";

    UserData Load(const char* deck, int cells) {
      UserData data;
      LoadFromJson(data, std::string("../examples/decks/") + deck);
      if (cells > 0)
        for (auto& layer : data.bed.layers_.layers_) layer.domain_.n_cells = cells;
      return data;
    }

    // Per-component scale of the state, |x_j| + abstol_j: errors are measured
    // in units of it, so concentrations, temperatures and loadings count alike.
    Vec Scale(const CycleSensitivityResult& r) {
      Vec s(r.x.size());
      for (size_t j = 0; j < s.size(); ++j) s[j] = std::abs(r.x[j]) + r.abstol[j];
      return s;
    }

    double ScaledNorm(const Vec& a, const Vec& scale) {
      double s = 0.0;
      for (size_t j = 0; j < a.size(); ++j) s += (a[j] / scale[j]) * (a[j] / scale[j]);
      return std::sqrt(s / static_cast<double>(a.size()));
    }

    double RelErr(const Vec& a, const Vec& ref, const Vec& scale) {
      Vec d(a.size());
      for (size_t j = 0; j < a.size(); ++j) d[j] = a[j] - ref[j];
      return ScaledNorm(d, scale) / ScaledNorm(ref, scale);
    }

    int CellOf(const std::string& description) {
      const auto at = description.rfind("cell ");
      return at == std::string::npos ? -1 : std::stoi(description.substr(at + 5));
    }

    // +1 K on the fluid and solid temperatures, a smooth (sin^2) bump of
    // width w centred at the bed fraction zc.
    Vec TemperatureBump(const CycleSensitivityResult& r, double zc, double w) {
      int max_cell = 0;
      for (const auto& d : r.description) max_cell = std::max(max_cell, CellOf(d));
      Vec v(r.x.size(), 0.0);
      for (size_t j = 0; j < v.size(); ++j) {
        const auto& d = r.description[j];
        if (d.find("Fluid Temperature") == std::string::npos && d.find("Solid Temperature") == std::string::npos)
          continue;
        const int cell = CellOf(d);
        if (cell < 0) continue;
        const double z = (cell + 0.5) / (max_cell + 1.0);
        if (std::abs(z - zc) < 0.5 * w) {
          const double s = std::cos(M_PI * (z - zc) / w);
          v[j] = s * s;
        }
      }
      return v;
    }

    // A random direction of scaled RMS norm `size`.
    Vec RandomDirection(const Vec& scale, unsigned seed, double size) {
      std::mt19937 gen(seed);
      std::normal_distribution<double> dist(0.0, 1.0);
      Vec v(scale.size());
      for (size_t j = 0; j < v.size(); ++j) v[j] = dist(gen) * scale[j];
      const double n = ScaledNorm(v, scale);
      for (double& x : v) x *= size / n;
      return v;
    }

    Options Base(int warmup) {
      Options o;
      o.warmup_cycles = warmup;
      o.sens_reltol = 1e-6;
      return o;
    }

    CycleSensitivityResult Plain(UserData& data, int warmup, const Vec& shift, double reltol) {
      Options o = Base(warmup);
      o.sensitivities = false;
      o.shift = shift;
      o.reltol = reltol;
      return CycleJacobianVector(data, {}, o);
    }

    // Centred difference of G along v with step h (in units of v).
    Vec CentredDifference(UserData& data, int warmup, const Vec& v, double h, double reltol) {
      Vec plus(v.size()), minus(v.size());
      for (size_t j = 0; j < v.size(); ++j) { plus[j] = h * v[j]; minus[j] = -h * v[j]; }
      const auto gp = Plain(data, warmup, plus, reltol);
      const auto gm = Plain(data, warmup, minus, reltol);
      EXPECT_TRUE(gp.ok) << gp.message;
      EXPECT_TRUE(gm.ok) << gm.message;
      Vec d(v.size(), 0.0);
      if (!gp.ok || !gm.ok) return d;
      for (size_t j = 0; j < v.size(); ++j) d[j] = (gp.g[j] - gm.g[j]) / (2.0 * h);
      return d;
    }

    const char* Name(Options::Corrector c) {
      return c == Options::Corrector::Staggered ? "staggered" : "simultaneous";
    }

  }  // namespace

  TEST(CycleSensitivity, MatchesCentredDifferenceOnASmoothCycle) {
    UserData data = Load(kSmoothDeck, 10);
    const auto base = Plain(data, 1, {}, 0.0);
    ASSERT_TRUE(base.ok) << base.message;
    const Vec scale = Scale(base);
    const Vec v = TemperatureBump(base, 0.5, 0.4);
    ASSERT_GT(ScaledNorm(v, scale), 0.0);
    const Vec fd = CentredDifference(data, 1, v, 1e-2, 1e-10);
    const Vec fd2 = CentredDifference(data, 1, v, 2e-2, 1e-10);
    const double fd_err = RelErr(fd2, fd, scale);

    for (auto corrector : { Options::Corrector::Staggered, Options::Corrector::Simultaneous }) {
      Options o = Base(1);
      o.corrector = corrector;
      const auto r = CycleJacobianVector(data, { v }, o);
      ASSERT_TRUE(r.ok) << Name(corrector) << ": " << r.message;
      ASSERT_EQ(r.jv.size(), 1u);
      EXPECT_GE(r.reinits, 1) << "the cycle start restarts the sensitivities";
      const double err = RelErr(r.jv[0], fd, scale);
      std::printf("  %s: |Phi v - FD| / |FD| = %.2e (FD error estimate %.2e)\n", Name(corrector), err, fd_err);
      EXPECT_LT(err, 1e-2) << Name(corrector);
    }
  }

  TEST(CycleSensitivity, RejectsConnectedDecks) {
    UserData data = Load("FiveStepProductPurge.json", 0);
    const auto r = CycleJacobianVector(data, {}, Options{});
    EXPECT_FALSE(r.ok);
    EXPECT_NE(r.message.find("connected"), std::string::npos) << r.message;
  }

  TEST(CycleSensitivity, MeasureCostAndAccuracy) {
    // --- Casas 2012 (single step), the deck's own grid, one warm-up cycle.
    UserData data = Load(kSmoothDeck, 0);
    const int warmup = 1;
    const auto base = Plain(data, warmup, {}, 0.0);
    ASSERT_TRUE(base.ok) << base.message;
    const Vec scale = Scale(base);
    std::printf("%s: %zu differential unknowns, reltol %.0e\n", kSmoothDeck, base.x.size(),
                static_cast<double>(data.settings.reltol));

    double plain_wall = 1e300;
    long plain_res = 0, plain_steps = 0;
    for (int rep = 0; rep < 5; ++rep) {
      const auto p = Plain(data, warmup, {}, 0.0);
      ASSERT_TRUE(p.ok) << p.message;
      plain_wall = std::min(plain_wall, p.wall_s);
      plain_res = p.residual_calls;
      plain_steps = p.steps;
    }
    std::printf("plain cycle: %.4f s, %ld residual calls, %ld steps\n", plain_wall, plain_res, plain_steps);
    std::printf("| corrector | DQ | directions | wall ratio | extra wall per direction [cycles] |"
                " extra residual calls per direction [cycles] | steps | sens. iterations | sens. conv. failures |\n");
    for (auto corrector : { Options::Corrector::Staggered, Options::Corrector::Simultaneous }) {
      for (auto dq : { Options::DifferenceQuotient::Centred, Options::DifferenceQuotient::Forward }) {
        for (int nv : { 1, 4, 16 }) {
          std::vector<Vec> dirs;
          for (int k = 0; k < nv; ++k) dirs.push_back(TemperatureBump(base, (k + 0.5) / nv, std::max(0.125, 1.0 / nv)));
          Options o = Base(warmup);
          o.corrector = corrector;
          o.dq = dq;
          double wall = 1e300;
          CycleSensitivityResult r;
          for (int rep = 0; rep < 3; ++rep) {
            r = CycleJacobianVector(data, dirs, o);
            if (!r.ok) break;
            wall = std::min(wall, r.wall_s);
          }
          const char* dq_name = dq == Options::DifferenceQuotient::Forward ? "forward" : "centred";
          if (!r.ok) {
            std::printf("| %s | %s | %d | failed: %s |\n", Name(corrector), dq_name, nv, r.message.c_str());
            continue;
          }
          std::printf("| %s | %s | %d | %.2f | %.2f | %.2f | %ld | %ld | %ld |\n", Name(corrector), dq_name, nv,
                      wall / plain_wall, (wall / plain_wall - 1.0) / nv,
                      (static_cast<double>(r.residual_calls) / plain_res - 1.0) / nv, r.steps, r.sens_nonlin_iters,
                      r.sens_conv_fails);
        }
      }
    }

    std::printf("| direction | FD error estimate | Phi v (deck reltol) vs FD | Phi v (reltol 1e-10) vs FD |\n");
    std::vector<std::pair<std::string, Vec>> dirs = {
        { "temperature bump, middle", TemperatureBump(base, 0.5, 0.4) },
        { "temperature bump, inlet quarter", TemperatureBump(base, 0.125, 0.25) },
        { "temperature bump, outlet quarter", TemperatureBump(base, 0.875, 0.25) },
        { "random (all unknowns, 1e-3 of the state)", RandomDirection(scale, 7, 1e-3) },
    };
    for (const auto& [label, v] : dirs) {
      const Vec fd1 = CentredDifference(data, warmup, v, 1e-2, 1e-10);
      const Vec fd2 = CentredDifference(data, warmup, v, 2e-2, 1e-10);
      Options o = Base(warmup);
      const auto r1 = CycleJacobianVector(data, { v }, o);
      o.reltol = 1e-10;
      const auto r2 = CycleJacobianVector(data, { v }, o);
      auto cell = [&](const CycleSensitivityResult& r) {
        return r.ok ? [&] { char b[32]; std::snprintf(b, sizeof b, "%.1e", RelErr(r.jv[0], fd1, scale)); return std::string(b); }()
                    : "failed (" + r.message + ")";
      };
      std::printf("| %s | %.1e | %s | %s |\n", label.c_str(), RelErr(fd2, fd1, scale), cell(r1).c_str(),
                  cell(r2).c_str());
      // A completed run must agree with the difference where the
      // difference itself is reliable; failures are the measurement.
      if (r2.ok && RelErr(fd2, fd1, scale) < 1e-2) EXPECT_LT(RelErr(r2.jv[0], fd1, scale), 2e-2) << label;
    }

    // --- Perez 2019 A (four-step PSA): sensitivities across step switches.
    UserData psa = Load(kPsaDeck, 0);
    const auto pbase = Plain(psa, 2, {}, 0.0);
    ASSERT_TRUE(pbase.ok) << pbase.message;
    const Vec pscale = Scale(pbase);
    for (auto corrector : { Options::Corrector::Staggered, Options::Corrector::Simultaneous }) {
      for (const auto& [label, v] : std::vector<std::pair<std::string, Vec>>{
               { "temperature bump", TemperatureBump(pbase, 0.5, 0.4) },
               { "random", RandomDirection(pscale, 11, 1e-3) } }) {
        Options o = Base(2);
        o.corrector = corrector;
        const auto r = CycleJacobianVector(psa, { v }, o);
        std::printf("%s, %s, %s: %s\n", kPsaDeck, Name(corrector), label.c_str(),
                    r.ok ? "completed" : r.message.c_str());
      }
    }
  }

}  // namespace adsim
