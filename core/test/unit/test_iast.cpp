/*
 * Skarstrom - pressure and temperature swing adsorption simulator
 * Copyright (c) 2024-2026 Sam Affleck and contributors
 *
 * Licensed under the terms in the LICENSE file at the repository root.
 */
// Ideal Adsorbed Solution Theory (adsim/equations/iast.h) and the adsorbent
// mixing rule that selects it (AdsorptionSource::mixing_rule_).
//
//  * the reduced spreading pressures: closed forms against a brute-force
//    integral, the quadrature kinds against their closed-form special cases,
//    and the inversion psi -> P° round trip for every kind;
//  * the pure loading of the prepared isotherm matches the engine's own
//    EvalIsotherm for every kind at its default parameters;
//  * the solve: the pure limit, the classic identity (equal-capacity Langmuir
//    IAST == extended Langmuir), a binary against a brute-force bisection,
//    the Raoult / Gibbs relations in a ternary with an inert and an absent
//    component, and the warm start;
//  * the deck: the mixing rule in the registry, in a save / load round trip,
//    the old-deck default, the deck warning, and a breakthrough run under IAST
//    that reproduces the extended-Langmuir run with equal capacities.

#include <algorithm>
#include <cmath>
#include <filesystem>
#include <map>
#include <sstream>
#include <string>
#include <vector>

#include <gtest/gtest.h>

#include "../src/pch.h"
#include "adsim/equations/iast.h"
#include "adsim/equations/isotherm.h"
#include "adsim/equations/isotherm_kinds.h"
#include "adsim/math/constants.h"
#include "adsim/objects/deck_builder.h"
#include "adsim/objects/deck_io.h"
#include "adsim/objects/porous_media.h"
#include "adsim/simulation.h"
#include "adsim/solver/deck_warnings.h"

namespace adsim {
namespace {

constexpr double kT = 298.15;
const std::string kComp = "A";

// A pure isotherm of `kind` with its default parameters, some overridden.
std::unique_ptr<IIsotherm> Make(IsothermType kind, const std::map<std::string, double>& values = {}) {
  auto iso = CreateIsotherm(kind, kComp);
  for (auto& [name, p] : IsothermFields(iso.get(), kComp)) {
    auto it = values.find(name);
    if (it != values.end()) p->Get() = it->second;
  }
  return iso;
}

// Brute-force ∫_0^P q(p)/p dp = ∫ q(e^u) du by the trapezoid rule in u, from a
// pressure 1e-12 P below which the isotherm is a power law q ∝ p^n (n = 1 in
// the Henry regime, 1/v for Freundlich) whose tail integrates to q(p_min)/n.
double BruteForcePsi(const PureIsotherm& iso, double P) {
  const int n = 400000;
  const double u1 = std::log(P), u0 = u1 - std::log(1e12);
  const double h = (u1 - u0) / n;
  double sum = 0.5 * (PureLoading(iso, std::exp(u0)) + PureLoading(iso, P));
  for (int i = 1; i < n; ++i) sum += PureLoading(iso, std::exp(u0 + i * h));
  const double q_min = PureLoading(iso, std::exp(u0));
  const double slope = std::log(PureLoading(iso, std::exp(u0 + 0.01)) / q_min) / 0.01;
  return h * sum + q_min / slope;
}

// Every kind with the parameter values the tests use (defaults are either
// tiny b's meant for Pa or placeholder magnitudes; these give loadings of
// order 1 mol/kg between 0.1 and 10 bar).
struct Case {
  IsothermType kind;
  std::map<std::string, double> values;
};

std::vector<Case> AllCases() {
  return {
      {IsothermType::HENRY, {{"k", 1e-6}, {"Q", 5e3}}},
      {IsothermType::SINGLE_SITE_LANGMUIR_P, {{"m1", 3.0}, {"b1", 2e-9}, {"q1", 2.5e4}}},
      {IsothermType::SINGLE_SITE_LANGMUIR_C, {{"m1", 3.0}, {"b1", 2e-5}, {"q1", 2.0e4}}},
      {IsothermType::DUAL_SITE_LANGMUIR_P, {{"m1", 2.0}, {"b1", 3e-9}, {"q1", 2.5e4}, {"m2", 1.5}, {"b2", 1e-10}, {"q2", 3.0e4}}},
      {IsothermType::DUAL_SITE_LANGMUIR_C, {{"m1", 2.0}, {"b1", 3e-5}, {"q1", 2.0e4}, {"m2", 1.5}, {"b2", 1e-6}, {"q2", 2.5e4}}},
      {IsothermType::EXTENDED_SINGLE_SITE_LANGMUIR_P, {{"m1", 3.0}, {"b1", 2e-9}, {"q1", 2.5e4}}},
      {IsothermType::EXTENDED_DUAL_SITE_LANGMUIR_P, {{"m1", 2.0}, {"b1", 3e-9}, {"q1", 2.5e4}, {"m2", 1.5}, {"b2", 1e-10}, {"q2", 3.0e4}}},
      {IsothermType::EXTENDED_SINGLE_SITE_LANGMUIR_C, {{"m1", 3.0}, {"b1", 2e-5}, {"q1", 2.0e4}}},
      {IsothermType::EXTENDED_DUAL_SITE_LANGMUIR_C, {{"m1", 2.0}, {"b1", 3e-5}, {"q1", 2.0e4}, {"m2", 1.5}, {"b2", 1e-6}, {"q2", 2.5e4}}},
      {IsothermType::FREUNDLICH, {{"a", 1e-5}, {"v", 2.0}}},
      {IsothermType::SIPS, {{"q_sat", 4.0}, {"b", 3e-5}, {"v", 1.4}}},
      {IsothermType::SIPS_2, {{"A1", 4.0}, {"A2", -0.6}, {"B1", 3.6e-9}, {"B2", 1.5e4}, {"C1", 0.8}, {"C2", 0.1}, {"T_ref", 298.0}}},
      {IsothermType::EXTENDED_SIPS_2, {{"A1", 4.0}, {"A2", -0.6}, {"B1", 3.6e-9}, {"B2", 1.5e4}, {"C1", 0.8}, {"C2", 0.1}, {"T_ref", 298.0}}},
      {IsothermType::DUAL_SITE_SIPS, {{"q_sat1", 2.0}, {"b1", 3e-5}, {"v1", 1.3}, {"q_sat2", 1.5}, {"b2", 2e-6}, {"v2", 0.9}}},
      {IsothermType::REDLICH_PETERSON, {{"a", 1e-4}, {"b", 2e-5}, {"v", 0.9}}},
      {IsothermType::TOTH, {{"q_sat", 4.0}, {"b", 2e-5}, {"v", 0.6}}},
      {IsothermType::TOTH_2, {{"q_sat", 3.715}, {"x", 0.1}, {"T_0", 324.0}, {"b0", 2e-5}, {"dH0", 2e4}, {"w0", 0.5}, {"alpha", 0.3}}},
      {IsothermType::ARRHENIUS_SIPS, {{"n_inf_a", 1.0}, {"n_inf_b", 3e3}, {"k_A", 1e-8}, {"k_B", 2e4}, {"c_ref", 0.9}, {"c_amp", 0.1}, {"c_beta", 0.05}, {"c_Tref", 298.15}}},
      {IsothermType::ANTI_LANGMUIR, {{"a", 1e-6}, {"b", 5e-7}}},
  };
}

}  // namespace

// ---- spreading pressures ------------------------------------------------------------

TEST(Iast, PureLoadingMatchesTheEngineIsotherm) {
  for (const auto& c : AllCases()) {
    SCOPED_TRACE(IsothermKindString(c.kind));
    auto iso = Make(c.kind, c.values);
    for (double T : {273.15, 298.15, 350.0}) {
      const PureIsotherm pure = PreparePureIsotherm(iso.get(), kComp, T);
      for (double P : {1e3, 2e4, 1e5, 5e5, 1e6}) {
        const double expected = EvalIsotherm(iso.get(), kComp, T, P, 1.0);
        EXPECT_NEAR(PureLoading(pure, P), expected, 1e-12 * std::abs(expected) + 1e-300) << "T " << T << " P " << P;
      }
    }
  }
  // Inert and a component without parameters are inert.
  auto inert = Make(IsothermType::INERT);
  EXPECT_EQ(PreparePureIsotherm(inert.get(), kComp, kT).kind, IsothermType::INERT);
  auto langmuir = Make(IsothermType::SINGLE_SITE_LANGMUIR_P);
  EXPECT_EQ(PreparePureIsotherm(langmuir.get(), "not_a_component", kT).kind, IsothermType::INERT);
  EXPECT_EQ(PreparePureIsotherm(nullptr, kComp, kT).kind, IsothermType::INERT);
}

TEST(Iast, SpreadingPressuresMatchABruteForceIntegral) {
  for (const auto& c : AllCases()) {
    SCOPED_TRACE(IsothermKindString(c.kind));
    auto iso = Make(c.kind, c.values);
    const PureIsotherm pure = PreparePureIsotherm(iso.get(), kComp, kT);
    for (double P : {2e4, 1e5, 8e5}) {
      const double psi = ReducedSpreadingPressure(pure, P);
      const double ref = BruteForcePsi(pure, P);
      EXPECT_GT(psi, 0.0);
      EXPECT_NEAR(psi, ref, 1e-7 * ref) << "P " << P;
    }
  }
  const PureIsotherm any = PreparePureIsotherm(Make(IsothermType::TOTH).get(), kComp, kT);
  EXPECT_EQ(ReducedSpreadingPressure(any, 0.0), 0.0);
  EXPECT_TRUE(std::isnan(ReducedSpreadingPressure(any, std::nan(""))));
}

TEST(Iast, QuadratureKindsReduceToTheirClosedFormSpecialCases) {
  // Toth with v = 1 is Langmuir: psi = q_sat ln(1 + b P).
  const PureIsotherm toth = PreparePureIsotherm(Make(IsothermType::TOTH, {{"q_sat", 4.0}, {"b", 2e-5}, {"v", 1.0}}).get(), kComp, kT);
  // Redlich-Peterson with v = 1 is Langmuir with m = a / b.
  const PureIsotherm rp = PreparePureIsotherm(Make(IsothermType::REDLICH_PETERSON, {{"a", 1e-4}, {"b", 2e-5}, {"v", 1.0}}).get(), kComp, kT);
  for (double P : {1e2, 1e4, 1e5, 1e6, 1e7}) {
    EXPECT_NEAR(ReducedSpreadingPressure(toth, P), 4.0 * std::log1p(2e-5 * P), 1e-11 * 4.0 * std::log1p(2e-5 * P)) << P;
    EXPECT_NEAR(ReducedSpreadingPressure(rp, P), 5.0 * std::log1p(2e-5 * P), 1e-11 * 5.0 * std::log1p(2e-5 * P)) << P;
  }
}

TEST(Iast, InvertingTheSpreadingPressureRoundTrips) {
  for (const auto& c : AllCases()) {
    SCOPED_TRACE(IsothermKindString(c.kind));
    auto iso = Make(c.kind, c.values);
    const PureIsotherm pure = PreparePureIsotherm(iso.get(), kComp, kT);
    for (double P : {5e2, 2e4, 1e5, 8e5}) {
      const double psi = ReducedSpreadingPressure(pure, P);
      // Cold (Henry estimate), warm, and a wrong guess all land on P.
      EXPECT_NEAR(InvertSpreadingPressure(pure, psi), P, 1e-9 * P) << "cold, P " << P;
      EXPECT_NEAR(InvertSpreadingPressure(pure, psi, P), P, 1e-9 * P) << "warm, P " << P;
      EXPECT_NEAR(InvertSpreadingPressure(pure, psi, 50.0 * P), P, 1e-9 * P) << "high guess, P " << P;
    }
    EXPECT_EQ(InvertSpreadingPressure(pure, 0.0), 0.0);
  }
}

// ---- the solve ------------------------------------------------------------------------

TEST(Iast, OneAdsorbingComponentIsItsPureIsothermAtThePartialPressure) {
  auto a = Make(IsothermType::TOTH, {{"q_sat", 4.0}, {"b", 2e-5}, {"v", 0.6}});
  auto inert = Make(IsothermType::INERT);
  std::vector<IIsotherm*> isos{a.get(), inert.get()};
  std::vector<std::string> names{kComp, kComp};
  std::vector<double> q;
  const double P = 2e5, y = 0.3;
  const auto r = EvalIastMixture(isos, names, kT, P, {y, 1.0 - y}, q);
  ASSERT_TRUE(r.converged);
  EXPECT_NEAR(q[0], EvalIsotherm(a.get(), kComp, kT, P, y), 1e-13);
  EXPECT_EQ(q[1], 0.0);

  // Nothing adsorbing, or no pressure: zero loadings; a bad input gives NaN.
  EXPECT_TRUE(EvalIastMixture(isos, names, kT, 0.0, {y, 1.0 - y}, q).converged);
  EXPECT_EQ(q[0], 0.0);
  EXPECT_TRUE(EvalIastMixture(isos, names, kT, P, {0.0, 1.0}, q).converged);
  EXPECT_EQ(q[0], 0.0);
  EXPECT_FALSE(EvalIastMixture(isos, names, kT, std::nan(""), {y, 1.0 - y}, q).converged);
  EXPECT_TRUE(std::isnan(q[0]));
}

// The textbook identity: with equal saturation capacities, IAST over Langmuir
// isotherms is exactly the extended Langmuir model.
TEST(Iast, EqualCapacityLangmuirIsExtendedLangmuir) {
  const double m = 3.2;
  const std::vector<double> b{2e-9, 3e-10, 8e-11}, Q{2.5e4, 1.8e4, 1.5e4};
  std::vector<std::unique_ptr<IIsotherm>> owned;
  std::vector<IIsotherm*> isos;
  std::vector<std::string> names;
  for (size_t i = 0; i < 3; ++i) {
    owned.push_back(Make(IsothermType::SINGLE_SITE_LANGMUIR_P, {{"m1", m}, {"b1", b[i]}, {"q1", Q[i]}}));
    isos.push_back(owned.back().get());
    names.push_back(kComp);
  }
  std::vector<double> q;
  for (double T : {280.0, 298.15, 330.0}) {
    for (double P : {5e4, 1e5, 6e5}) {
      for (const std::vector<double>& y : {std::vector<double>{0.15, 0.80, 0.05}, {0.5, 0.25, 0.25}, {0.02, 0.9, 0.08}}) {
        const auto r = EvalIastMixture(isos, names, T, P, y, q);
        ASSERT_TRUE(r.converged) << "T " << T << " P " << P;
        double denominator = 1.0;
        std::vector<double> B(3);
        for (size_t i = 0; i < 3; ++i) {
          B[i] = b[i] * std::exp(Q[i] / (constants::kIdealGasConstant * T));
          denominator += B[i] * P * y[i];
        }
        for (size_t i = 0; i < 3; ++i) {
          const double expected = m * B[i] * P * y[i] / denominator;
          EXPECT_NEAR(q[i], expected, 1e-10 * expected) << "T " << T << " P " << P << " i " << i;
        }
      }
    }
  }
}

// A binary with unequal capacities against a bisection on psi written out
// with the Langmuir closed forms, and the IAST relations themselves.
TEST(Iast, BinaryLangmuirAgainstABruteForceSolve) {
  const double m1 = 4.0, m2 = 2.5, B1 = 3e-5, B2 = 4e-6;  // q = m B P / (1 + B P) at this T (Q = 0)
  auto a = Make(IsothermType::SINGLE_SITE_LANGMUIR_P, {{"m1", m1}, {"b1", B1}, {"q1", 0.0}});
  auto bb = Make(IsothermType::SINGLE_SITE_LANGMUIR_P, {{"m1", m2}, {"b1", B2}, {"q1", 0.0}});
  std::vector<IIsotherm*> isos{a.get(), bb.get()};
  std::vector<std::string> names{kComp, kComp};
  const double P = 1.5e5, y1 = 0.2, y2 = 0.8;

  // Bisection on psi: F(psi) = y1 P / P1°(psi) + y2 P / P2°(psi) - 1.
  auto P0 = [](double psi, double m, double B) { return std::expm1(psi / m) / B; };
  double lo = 0.0, hi = 100.0;
  for (int it = 0; it < 300; ++it) {
    const double mid = 0.5 * (lo + hi);
    const double F = y1 * P / P0(mid, m1, B1) + y2 * P / P0(mid, m2, B2) - 1.0;
    (F > 0.0 ? lo : hi) = mid;
  }
  const double psi = 0.5 * (lo + hi);
  const double P1 = P0(psi, m1, B1), P2 = P0(psi, m2, B2);
  const double x1 = y1 * P / P1, x2 = y2 * P / P2;
  const double q1p = m1 * B1 * P1 / (1 + B1 * P1), q2p = m2 * B2 * P2 / (1 + B2 * P2);
  const double qT = 1.0 / (x1 / q1p + x2 / q2p);

  std::vector<double> q;
  const auto r = EvalIastMixture(isos, names, kT, P, {y1, y2}, q);
  ASSERT_TRUE(r.converged);
  EXPECT_NEAR(x1 + x2, 1.0, 1e-12);
  EXPECT_NEAR(r.psi, psi, 1e-10 * psi);
  EXPECT_NEAR(q[0], x1 * qT, 1e-10 * q[0]);
  EXPECT_NEAR(q[1], x2 * qT, 1e-10 * q[1]);
  // The strongly adsorbed component is enriched in the adsorbed phase; each
  // loads less than it would alone (the strong one at the total pressure, the
  // weak one at its partial pressure).
  EXPECT_GT(x1, y1);
  EXPECT_LT(q[0], EvalIsotherm(a.get(), kComp, kT, P, 1.0));
  EXPECT_LT(q[1], EvalIsotherm(bb.get(), kComp, kT, P, y2));
}

// Mixed kinds, an inert carrier and a component absent from the gas: the
// loadings satisfy the IAST relations (equal spreading pressures at the
// fictitious pure pressures, the Raoult analogue, the total-loading mixing
// rule) and the absent and inert components load nothing.
TEST(Iast, TernaryWithMixedKindsInertAndAbsentComponent) {
  auto a = Make(IsothermType::TOTH, {{"q_sat", 4.0}, {"b", 2e-5}, {"v", 0.6}});
  auto b = Make(IsothermType::DUAL_SITE_LANGMUIR_P, {{"m1", 2.0}, {"b1", 3e-9}, {"q1", 2.5e4}, {"m2", 1.5}, {"b2", 1e-10}, {"q2", 3.0e4}});
  auto c = Make(IsothermType::SIPS, {{"q_sat", 3.0}, {"b", 1e-5}, {"v", 1.3}});
  auto d = Make(IsothermType::SINGLE_SITE_LANGMUIR_P, {{"m1", 3.0}, {"b1", 2e-9}, {"q1", 2.5e4}});
  auto inert = Make(IsothermType::INERT);
  std::vector<IIsotherm*> isos{a.get(), b.get(), c.get(), d.get(), inert.get()};
  std::vector<std::string> names(5, kComp);
  const double P = 3e5;
  const std::vector<double> y{0.2, 0.3, 0.1, 0.0, 0.4};
  std::vector<double> q;
  const auto r = EvalIastMixture(isos, names, kT, P, y, q);
  ASSERT_TRUE(r.converged);
  EXPECT_EQ(q[3], 0.0);
  EXPECT_EQ(q[4], 0.0);

  double qT = 0.0;
  for (int i = 0; i < 3; ++i) qT += q[i];
  double inv_qT = 0.0;
  for (int i = 0; i < 3; ++i) {
    const double x = q[i] / qT;
    const double Pi = y[i] * P / x;  // Raoult's analogue
    const PureIsotherm pure = PreparePureIsotherm(isos[i], kComp, kT);
    EXPECT_NEAR(ReducedSpreadingPressure(pure, Pi), r.psi, 1e-9 * r.psi) << i;
    inv_qT += x / PureLoading(pure, Pi);
    EXPECT_LT(q[i], PureLoading(pure, P));
  }
  EXPECT_NEAR(1.0 / inv_qT, qT, 1e-9 * qT);
}

// Ahead of a front the integrator produces mole fractions of order ±1e-12:
// the loading must be smooth through y = 0 (continuous, with a continuous
// slope), or every residual call there sees a kink and the run crawls.
TEST(Iast, LoadingIsSmoothThroughZeroMoleFraction) {
  auto a = Make(IsothermType::TOTH, {{"q_sat", 4.0}, {"b", 2e-5}, {"v", 0.6}});
  auto b = Make(IsothermType::DUAL_SITE_LANGMUIR_P, {{"m1", 2.0}, {"b1", 3e-9}, {"q1", 2.5e4}, {"m2", 1.5}, {"b2", 1e-10}, {"q2", 3.0e4}});
  auto inert = Make(IsothermType::INERT);
  std::vector<IIsotherm*> isos{a.get(), b.get(), inert.get()};
  std::vector<std::string> names(3, kComp);
  const double P = 2e5, ya = 0.3;
  auto qb = [&](double yb) {
    std::vector<double> q;
    EXPECT_TRUE(EvalIastMixture(isos, names, kT, P, {ya, yb, 1.0 - ya - yb}, q).converged) << yb;
    return q;
  };
  const double h = 1e-9;
  const auto minus = qb(-h), zero = qb(0.0), plus = qb(h), plus2 = qb(2 * h);
  EXPECT_EQ(zero[1], 0.0);
  EXPECT_LT(minus[1], 0.0);
  EXPECT_GT(plus[1], 0.0);
  // Equal slopes either side of zero (the dilute limit is the first-order
  // limit of the full solve), and the other component's loading is continuous.
  const double slope_left = (zero[1] - minus[1]) / h, slope_right = (plus2[1] - plus[1]) / h;
  EXPECT_NEAR(slope_left, slope_right, 1e-5 * slope_right);
  EXPECT_NEAR(minus[0], plus[0], 1e-7 * plus[0]);
  EXPECT_EQ(minus[2], 0.0);

  // Nothing in the gas: every adsorbing component is dilute, q = y P K.
  std::vector<double> q;
  ASSERT_TRUE(EvalIastMixture(isos, names, kT, P, {-1e-12, -2e-12, 1.0}, q).converged);
  const double Ka = DiluteSlope(PreparePureIsotherm(a.get(), kComp, kT));
  const double Kb = DiluteSlope(PreparePureIsotherm(b.get(), kComp, kT));
  EXPECT_NEAR(q[0], -1e-12 * P * Ka, 1e-12 * std::abs(q[0]));
  EXPECT_NEAR(q[1], -2e-12 * P * Kb, 1e-12 * std::abs(q[1]));
  EXPECT_EQ(q[2], 0.0);
}

TEST(Iast, AWarmStartConvergesInAStepOrTwo) {
  auto a = Make(IsothermType::DUAL_SITE_LANGMUIR_P, {{"m1", 2.0}, {"b1", 3e-9}, {"q1", 2.5e4}, {"m2", 1.5}, {"b2", 1e-10}, {"q2", 3.0e4}});
  auto b = Make(IsothermType::TOTH, {{"q_sat", 4.0}, {"b", 2e-7}, {"v", 0.7}});
  PureIsotherm pure[2] = {PreparePureIsotherm(a.get(), kComp, kT), PreparePureIsotherm(b.get(), kComp, kT)};
  const double y[2] = {0.15, 0.85};
  double q[2];
  IastScratch scratch;
  const auto cold = SolveIast(pure, 2, 1e5, y, q, scratch);
  ASSERT_TRUE(cold.converged);
  const double q_cold[2] = {q[0], q[1]};
  const auto warm = SolveIast(pure, 2, 1e5, y, q, scratch, cold.psi);
  ASSERT_TRUE(warm.converged);
  EXPECT_LE(warm.iterations, 2);
  EXPECT_LE(warm.iterations, cold.iterations);
  EXPECT_NEAR(q[0], q_cold[0], 1e-12 * q_cold[0]);
  EXPECT_NEAR(q[1], q_cold[1], 1e-12 * q_cold[1]);
  // A nearby state (the next cell) from the last psi: still a short iteration.
  const double y2[2] = {0.16, 0.84};
  const auto next = SolveIast(pure, 2, 1.02e5, y2, q, scratch, cold.psi);
  ASSERT_TRUE(next.converged);
  EXPECT_LE(next.iterations, 4);
}

// ---- the deck -----------------------------------------------------------------------

namespace {

// A one-step CO2/N2 breakthrough on a single-site Langmuir adsorbent: the
// caller chooses the isotherm kind and the mixing rule.
struct Breakthrough {
  Simulation sim;
  AdsorbentID ads = 0;
};

void BuildBreakthrough(Breakthrough& b, const std::string& isotherm_kind, const std::string& mixing_rule) {
  UserData& d = b.sim.Data();
  build::AddComponent(d, "CarbonDioxide");
  build::AddComponent(d, "Nitrogen");
  b.ads = build::AddAdsorbent(d, "13X");
  build::SetIsotherm(d, b.ads, isotherm_kind, {"CarbonDioxide", "Nitrogen"});
  build::SetMixingRule(d, b.ads, mixing_rule);
  build::SetRateModel(d, b.ads, "CarbonDioxide", "ldf", "constant");
  build::SetRateModel(d, b.ads, "Nitrogen", "ldf", "constant");
  build::SetBedShape(d, "cylinder");
  build::SetWallMaterial(d, "ss316l");
  build::SetWallBoundary(d, "zero_flux");
  build::SetThermalEquilibrium(d, "lte");
  const size_t layer = build::AddLayer(d, "13X");
  build::SetLayerComposition(d, layer, {{b.ads, 1.0}});
  build::SetMomentum(d, layer, "ergun");
  build::SetDispersion(d, layer, "constant");
  build::SetHeatDispersion(d, layer, "constant");
  build::SetInternalHeatTransfer(d, layer, "constant");
  build::SetExternalHeatTransfer(d, layer, "constant");
  const size_t feed = build::AddStep(d, "Feed");
  build::SetFlowBoundary(d, feed, "bottom", "inlet");
  build::AddStream(d, feed, "bottom", "Feed inlet", "velocity", "constant");
  build::SetPressureBoundary(d, feed, "top", "outlet", "constant");
  b.sim.RefreshParameters();

  // Equal capacities (IAST over Langmuir == extended Langmuir), small enough
  // that CO2 breaks through the 0.1 m bed after about 90 s of the 150 s feed.
  const double m = 1.0;
  const std::vector<std::pair<std::string, std::string>> values = {
      {"bed.diameter", "0.032"},
      {"wall.thickness", "0.0031"},
      {"adsorbent[13X].solid_density", "1138"},
      {"adsorbent[13X].heat_capacity", "900"},
      {"adsorbent[13X].thermal_conductivity", "0.63"},
      {"adsorbent[13X].particle_void_fraction", "0.35"},
      {"adsorbent[13X].isotherm[CarbonDioxide].m1", std::to_string(m)},
      {"adsorbent[13X].isotherm[CarbonDioxide].b1", "2e-9"},
      {"adsorbent[13X].isotherm[CarbonDioxide].q1", "25000"},
      {"adsorbent[13X].isotherm[Nitrogen].m1", std::to_string(m)},
      {"adsorbent[13X].isotherm[Nitrogen].b1", "3e-10"},
      {"adsorbent[13X].isotherm[Nitrogen].q1", "15000"},
      {"adsorbent[13X].kinetics[CarbonDioxide].k", "0.3"},
      {"adsorbent[13X].kinetics[Nitrogen].k", "0.5"},
      {"adsorbent[13X].heat_of_adsorption[CarbonDioxide]", "30000"},
      {"adsorbent[13X].heat_of_adsorption[Nitrogen]", "15000"},
      {"layer[13X].length", "0.1"},
      {"layer[13X].cells", "15"},
      {"layer[13X].bed_void_fraction", "0.41"},
      {"layer[13X].particle_diameter", "0.0015"},
      {"layer[13X].dispersion_coefficient", "1e-4"},
      {"initial.pressure", "1e5"},
      {"initial.temperature", "298"},
      {"initial.y[Nitrogen]", "1"},
      {"initial.y[CarbonDioxide]", "0"},
      {"step[Feed].time", "150"},
      {"step[Feed].bottom[Feed inlet].velocity", "0.05"},
      {"step[Feed].bottom[Feed inlet].temperature", "298"},
      {"step[Feed].bottom[Feed inlet].y[CarbonDioxide]", "0.15"},
      {"step[Feed].bottom[Feed inlet].y[Nitrogen]", "0.85"},
      {"step[Feed].top.pressure", "1e5"},
  };
  for (const auto& [path, text] : values) b.sim.SetParameterText(path, text);
}

}  // namespace

TEST(IastDeck, MixingRuleIsARegistryParameterAndSurvivesSaveAndLoad) {
  Breakthrough b;
  BuildBreakthrough(b, "single_site_langmuir_p", "none");
  EXPECT_EQ(b.sim.GetParameter("adsorbent[13X].mixing_rule"), 0.0);
  EXPECT_EQ(b.sim.Data().adsorbent_directory.Get(b.ads).adsorption_source_.mixing_rule_, MixingRule::None);
  b.sim.SetParameterText("adsorbent[13X].mixing_rule", "iast");
  EXPECT_EQ(b.sim.GetParameter("adsorbent[13X].mixing_rule"), 1.0);
  EXPECT_EQ(b.sim.Data().adsorbent_directory.Get(b.ads).adsorption_source_.mixing_rule_, MixingRule::IAST);
  EXPECT_THROW(b.sim.SetParameterText("adsorbent[13X].mixing_rule", "rast"), std::exception);
  EXPECT_THROW(build::SetMixingRule(b.sim.Data(), b.ads, "rast"), std::invalid_argument);
  EXPECT_EQ(build::MixingRuleKinds(), (std::vector<std::string>{"none", "iast"}));

  // Save / load round trip through cereal's JSON archive, and an archive from
  // before the mixing rule existed (AdsorptionSource version 0) loads as None.
  std::ostringstream oss;
  {
    cereal::JSONOutputArchive ar(oss);
    ar(b.sim.Data());
  }
  const std::string json = oss.str();
  EXPECT_NE(json.find("\"Mixing_Rule\": 1"), std::string::npos);
  UserData loaded;
  ASSERT_NO_THROW(LoadDeckFlexible(loaded, json));
  EXPECT_EQ(loaded.adsorbent_directory.Get(b.ads).adsorption_source_.mixing_rule_, MixingRule::IAST);

  std::string old = json;
  const std::string version_line = "\"cereal_class_version\": 1,";
  const auto source = old.find("\"Adsorption_Source\"");
  ASSERT_NE(source, std::string::npos);
  const auto version = old.find(version_line, source);
  ASSERT_NE(version, std::string::npos);
  old.replace(version, version_line.size(), "\"cereal_class_version\": 0,");
  const auto rule = old.find("\"Mixing_Rule\": 1");
  ASSERT_NE(rule, std::string::npos);
  const auto rule_start = old.rfind(',', rule);  // the comma that ends the previous field
  old.erase(rule_start, rule + std::string("\"Mixing_Rule\": 1").size() - rule_start);
  UserData from_old;
  ASSERT_NO_THROW(LoadDeckFlexible(from_old, old));
  EXPECT_EQ(from_old.adsorbent_directory.Get(b.ads).adsorption_source_.mixing_rule_, MixingRule::None);
}

TEST(IastDeck, OldDecksLoadWithTheNoneMixingRule) {
  const auto decks = std::filesystem::current_path() / ".." / "examples" / "decks";
  UserData d;
  ASSERT_NO_THROW(LoadFromJson(d, (decks / "Perez2019_A.json").string()));
  for (const auto& [id, ads] : d.adsorbent_directory.GetAll())
    EXPECT_EQ(ads.adsorption_source_.mixing_rule_, MixingRule::None) << ads.name_;
}

TEST(IastDeck, WarnsWhenIastIsCombinedWithACompetitiveKind) {
  Breakthrough b;
  BuildBreakthrough(b, "extended_single_site_langmuir_p", "iast");
  const auto& d = b.sim.Data();
  const auto warnings = MixingRuleWarnings(d.adsorbent_directory, d.fluid);
  ASSERT_EQ(warnings.size(), 1u);
  EXPECT_NE(warnings[0].find("13X"), std::string::npos);
  EXPECT_NE(warnings[0].find("extended_single_site_langmuir_p"), std::string::npos);
  EXPECT_NE(warnings[0].find("CarbonDioxide, Nitrogen"), std::string::npos);

  Breakthrough plain;
  BuildBreakthrough(plain, "single_site_langmuir_p", "iast");
  EXPECT_TRUE(MixingRuleWarnings(plain.sim.Data().adsorbent_directory, plain.sim.Data().fluid).empty());
  Breakthrough none;
  BuildBreakthrough(none, "extended_single_site_langmuir_p", "none");
  EXPECT_TRUE(MixingRuleWarnings(none.sim.Data().adsorbent_directory, none.sim.Data().fluid).empty());
}

// The end-to-end check: with equal saturation capacities IAST over Langmuir is
// the extended Langmuir model, so the two breakthrough runs give the same
// outlet curves (the IAST path solves what the other writes in closed form).
TEST(IastDeck, BreakthroughUnderIastReproducesExtendedLangmuirWithEqualCapacities) {
  Breakthrough extended, iast;
  BuildBreakthrough(extended, "extended_single_site_langmuir_p", "none");
  BuildBreakthrough(iast, "single_site_langmuir_p", "iast");
  std::string error;
  ASSERT_TRUE(extended.sim.Check(error)) << error;
  ASSERT_TRUE(iast.sim.Check(error)) << error;

  RunOptions options;
  options.cycles = 1;
  const RunResult re = extended.sim.Run(options);
  ASSERT_TRUE(re.Ok()) << re.message;
  const RunResult ri = iast.sim.Run(options);
  ASSERT_TRUE(ri.Ok()) << ri.message;

  // The IAST loadings are smooth (including through zero mole fraction ahead
  // of the front), so the integrator takes the same path: a residual with a
  // kink or solver noise would show up here as many more steps.
  EXPECT_LE(ri.stats.num_steps, re.stats.num_steps * 11 / 10 + 5);
  // The runs log at their own integrator steps: compare on the first run's times.
  auto at = [](const TimeSeries& s, double t) {
    if (t <= s.time.front()) return s.values.front();
    if (t >= s.time.back()) return s.values.back();
    const auto it = std::upper_bound(s.time.begin(), s.time.end(), t);
    const size_t i = static_cast<size_t>(it - s.time.begin());
    const double w = (t - s.time[i - 1]) / (s.time[i] - s.time[i - 1]);
    return s.values[i - 1] + w * (s.values[i] - s.values[i - 1]);
  };
  for (const char* variable : {"mole_fraction", "adsorbed_concentration"}) {
    for (const char* component : {"CarbonDioxide", "Nitrogen"}) {
      const TimeSeries a = extended.sim.Series(variable, component, SeriesLocation::Top);
      const TimeSeries b = iast.sim.Series(variable, component, SeriesLocation::Top);
      ASSERT_GT(a.values.size(), 10u);
      ASSERT_GT(b.values.size(), 10u);
      double max_diff = 0.0, max_val = 0.0;
      for (size_t i = 0; i < a.values.size(); ++i) {
        max_diff = std::max(max_diff, std::abs(a.values[i] - at(b, a.time[i])));
        max_val = std::max(max_val, std::abs(a.values[i]));
      }
      EXPECT_LT(max_diff, 1e-6 * max_val) << variable << " " << component;
    }
  }
  EXPECT_GT(extended.sim.Series("adsorbed_concentration", "CarbonDioxide", SeriesLocation::Top).values.back(), 0.1);

  // A rerun of the same simulation starts its IAST iterations cold and
  // reproduces the first run exactly (the warm-start state is per run).
  const std::vector<double> first = iast.sim.Series("mole_fraction", "CarbonDioxide", SeriesLocation::Top).values;
  const RunResult again = iast.sim.Run(options);
  ASSERT_TRUE(again.Ok()) << again.message;
  EXPECT_EQ(iast.sim.Series("mole_fraction", "CarbonDioxide", SeriesLocation::Top).values, first);
}

}  // namespace adsim
