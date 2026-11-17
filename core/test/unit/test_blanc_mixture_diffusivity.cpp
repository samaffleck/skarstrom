// Blanc's law with a smooth floor on the partner mole fractions
// (BlancMixtureDiffusivity, binary_diffusion_model.h): agreement with plain
// Blanc away from the pure limit, a finite value at the pure limit, bounds,
// independence of the partner order, smoothness, and a derivative that the
// solver's difference quotient resolves.
#include <algorithm>
#include <cmath>
#include <numeric>
#include <random>
#include <vector>

#include <gtest/gtest.h>

#include "adsim/equations/binary_diffusion_model.h"

namespace adsim {

  namespace {

    double Blanc(const std::vector<double>& y, const std::vector<double>& D) {
      return BlancMixtureDiffusivity(y.data(), D.data(), y.size());
    }

    // Blanc's law without the floor: sum y_j / sum (y_j / D_ij).
    double PlainBlanc(const std::vector<double>& y, const std::vector<double>& D) {
      double num = 0.0, den = 0.0;
      for (size_t j = 0; j < y.size(); ++j) { num += y[j]; den += y[j] / D[j]; }
      return num / den;
    }

    // D_ij [m^2/s] of N2 with H2 and with CO2 at 300 K, 1 bar (Chapman-Enskog).
    const std::vector<double> kN2Partners = { 7.568e-5, 1.585e-5 };

  } // namespace

  TEST(BlancMixtureDiffusivity, SmoothFloorValues) {
    const double d = kBlancMoleFractionFloor;
    EXPECT_DOUBLE_EQ(SmoothMoleFractionFloor(0.0, d), 0.5 * d);
    EXPECT_NEAR(SmoothMoleFractionFloor(0.3, d), 0.3, 1e-12);
    // Both branches are one function: equal at y = 0 and positive for y < 0.
    EXPECT_NEAR(SmoothMoleFractionFloor(-1e-300, d), SmoothMoleFractionFloor(1e-300, d), 1e-20);
    EXPECT_GT(SmoothMoleFractionFloor(-0.5, d), 0.0);
    EXPECT_NEAR(SmoothMoleFractionFloor(-0.5, d), d * d / (4 * 0.5), 1e-6 * d * d);
  }

  // Away from the pure limit the floor does not change Blanc's law (measured:
  // 1.8e-7 with every partner >= 1e-3, 1.6e-5 with every partner >= 1e-4).
  TEST(BlancMixtureDiffusivity, EqualsPlainBlancAwayFromThePureLimit) {
    std::mt19937 rng(7);
    std::uniform_real_distribution<double> logy(-3.0, 0.0), logD(-5.5, -3.5);
    std::uniform_int_distribution<int> count(1, 5);
    double worst = 0.0;
    for (int trial = 0; trial < 5000; ++trial) {
      const int k = count(rng);
      std::vector<double> y(k), D(k);
      for (int j = 0; j < k; ++j) { y[j] = std::pow(10.0, logy(rng)); D[j] = std::pow(10.0, logD(rng)); }
      worst = std::max(worst, std::abs(Blanc(y, D) / PlainBlanc(y, D) - 1.0));
    }
    EXPECT_LT(worst, 1e-6);
  }

  // A trace of i in the others: Blanc's law is the exact Stefan-Maxwell limit.
  TEST(BlancMixtureDiffusivity, TraceLimit) {
    const std::vector<double> D = { 6.337e-5, 1.585e-5 };   // CO2 with H2 and with N2
    const std::vector<double> y = { 0.5, 0.5 };
    const double stefan_maxwell = 1.0 / (0.5 / D[0] + 0.5 / D[1]);
    EXPECT_NEAR(Blanc(y, D), stefan_maxwell, 1e-6 * stefan_maxwell);
  }

  // A cell of pure N2 in an H2/CO2/N2 deck. Without the floor this was 0, and
  // every diffusion-based mass-transfer coefficient of N2 was 0.
  TEST(BlancMixtureDiffusivity, PureLimitIsTheHarmonicMeanOfThePartners) {
    const auto& D = kN2Partners;
    const double harmonic = 2.0 / (1.0 / D[0] + 1.0 / D[1]);
    EXPECT_NEAR(Blanc({ 0.0, 0.0 }, D), harmonic, 1e-12 * harmonic);
    // A solver overshoot below zero stays finite and positive.
    const double overshoot = Blanc({ -1e-9, -3e-9 }, D);
    EXPECT_GT(overshoot, std::min(D[0], D[1]));
    EXPECT_LT(overshoot, std::max(D[0], D[1]));
  }

  // A weighted harmonic mean with positive weights: between min and max D_ij
  // for every composition, negative and tiny fractions included.
  TEST(BlancMixtureDiffusivity, BoundedByThePartnerDiffusivities) {
    std::mt19937 rng(11);
    std::uniform_real_distribution<double> u(-1.0, 1.0), logD(-5.5, -3.5);
    std::uniform_real_distribution<double> logy(-14.0, 0.0);
    for (int trial = 0; trial < 5000; ++trial) {
      const int k = 1 + trial % 5;
      std::vector<double> y(k), D(k);
      for (int j = 0; j < k; ++j) {
        const double mag = std::pow(10.0, logy(rng));
        y[j] = u(rng) < -0.8 ? -mag * 1e-6 : mag;   // some small negative overshoots
        D[j] = std::pow(10.0, logD(rng));
      }
      const double v = Blanc(y, D);
      ASSERT_TRUE(std::isfinite(v));
      EXPECT_GE(v, *std::min_element(D.begin(), D.end()) * (1 - 1e-12));
      EXPECT_LE(v, *std::max_element(D.begin(), D.end()) * (1 + 1e-12));
    }
  }

  TEST(BlancMixtureDiffusivity, IndependentOfThePartnerOrder) {
    std::vector<int> order = { 0, 1, 2, 3 };
    const std::vector<double> y = { 0.6, 1e-9, 0.25, 0.0 };
    const std::vector<double> D = { 7.6e-5, 1.6e-5, 2.1e-5, 6.3e-5 };
    const double ref = Blanc(y, D);
    do {
      std::vector<double> yp, Dp;
      for (int j : order) { yp.push_back(y[j]); Dp.push_back(D[j]); }
      EXPECT_NEAR(Blanc(yp, Dp), ref, 1e-14 * ref);
    } while (std::next_permutation(order.begin(), order.end()));
  }

  // The derivative in a partner fraction is continuous through y = 0 (the
  // switch between the two branches of the floor) and at the scale of the
  // floor: central differences at two steps agree.
  TEST(BlancMixtureDiffusivity, SmoothInThePartnerFractions) {
    const auto& D = kN2Partners;
    auto f = [&](double y_h2) { return Blanc({ y_h2, 0.0 }, D); };
    auto dfdy = [&](double y, double h) { return (f(y + h) - f(y - h)) / (2 * h); };
    const double d = kBlancMoleFractionFloor;
    EXPECT_NEAR(dfdy(-1e-14, 1e-12), dfdy(1e-14, 1e-12), 1e-4 * std::abs(dfdy(0.0, 1e-12)));
    for (double y = -5 * d; y <= 5 * d; y += 0.05 * d) {
      const double a = dfdy(y, 1e-3 * d), b = dfdy(y, 1e-4 * d);
      EXPECT_NEAR(a, b, 1e-4 * std::abs(b) + 1e-12) << "y = " << y;
    }
  }

  // The solver's difference-quotient Jacobian perturbs a mole fraction near
  // zero by about its absolute tolerance, 1e-8 (jacobian.cpp, Increment). Over
  // that step the secant equals the derivative within 1 % (measured 0.57 %),
  // so the Jacobian sees the floor's transition, not a jump.
  TEST(BlancMixtureDiffusivity, SolverDifferenceQuotientResolvesTheFloor) {
    const auto& D = kN2Partners;
    auto f = [&](double y_h2) { return Blanc({ y_h2, 0.0 }, D); };
    const double h_solver = 1e-8;
    double worst = 0.0;
    for (double y = -2e-5; y <= 2e-5; y += 1e-8) {
      const double derivative = (f(y + 1e-11) - f(y - 1e-11)) / 2e-11;
      if (std::abs(derivative) < 1e-3) continue;   // flat far from the pure limit
      const double secant = (f(y + h_solver) - f(y)) / h_solver;
      worst = std::max(worst, std::abs(secant / derivative - 1.0));
    }
    EXPECT_LT(worst, 0.01);
  }

} // namespace adsim
