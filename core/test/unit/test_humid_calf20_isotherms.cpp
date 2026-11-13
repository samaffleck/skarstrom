/*
 * Skarstrom - pressure and temperature swing adsorption simulator
 * Copyright (c) 2024-2026 Sam Affleck and contributors
 *
 * Licensed under the terms in the LICENSE file at the repository root.
 */
// Water on CALF-20 (Nguyen et al., Ind. Eng. Chem. Res. 63 (2024) 3265-3281):
//   * cubic_langmuir: the type V water isotherm (the paper's Cubic-Langmuir,
//     the Quadratic-Langmuir of T. T. T. Nguyen's thesis, Univ. of Alberta 2021),
//   * the optional water factors (lambda1, lambda2, lambda3, gamma2) of the
//     dual-site Langmuir kinds, which lower the CO2 loading with the adsorbed
//     water loading.

#include <algorithm>
#include <cmath>
#include <filesystem>
#include <fstream>
#include <iostream>
#include <sstream>
#include <stdexcept>
#include <string>
#include <tuple>
#include <vector>

#include <gtest/gtest.h>

#include "adsim/equations/equation.h"
#include "adsim/equations/iast.h"
#include "adsim/equations/isotherm.h"
#include "adsim/equations/isotherm_kinds.h"
#include "adsim/math/constants.h"
#include "adsim/objects/deck_builder.h"
#include "adsim/objects/deck_io.h"
#include "adsim/objects/model_kinds.h"
#include "adsim/objects/porous_media.h"
#include "adsim/simulation.h"
#include "adsim/solver/cyclic_steady_state.h"
#include "adsim/solver/isotherm_fit.h"

namespace adsim {
namespace {

const std::string kCO2 = "CarbonDioxide";
const std::string kN2 = "Nitrogen";
const std::string kWater = kCoadsorbedWaterComponent;

double Rel(double a, double b) { return std::fabs(a / b - 1.0); }

double Cubic(const CubicLangmuirIsothermParameters& ips, double T, double C) {
  double q = 0.0;
  IsothermFunction::CubicLangmuir(q, T, C * R * T, 1.0, ips);
  return q;
}

// The CO2 dual-site Langmuir of the paper's SI (Table S2, the same as Nguyen
// et al., Chem. Eng. J. 442 (2022) 136263): q_sb 2.387 mol/kg, b0 5.519e-7
// m3/mol, dU_b -35.06 kJ/mol; q_sd 3.271 mol/kg, d0 5.187e-8 m3/mol, dU_d
// -28.95 kJ/mol (the engine's q1, q2 are -dU). The per-site water factors of
// the paper (lambda1 0.0066, lambda2 2.483, lambda3 3.55, q_w in mol/kg).
DualSiteLangmuirIsothermParametersConcentration Calf20CO2(bool water) {
  DualSiteLangmuirIsothermParametersConcentration p;
  p.m1 = 2.387; p.b1 = 5.519e-7; p.q1 = 35060.0;
  p.m2 = 3.271; p.b2 = 5.187e-8; p.q2 = 28950.0;
  if (water) { p.lambda1 = 0.0066; p.lambda2 = 2.483; p.lambda3 = 3.55; }
  return p;
}

template <class Map, class... Args>
void AddCells(Map& map, const std::string& c, size_t n, double v, Args&&... args) {
  auto [it, ok] = map.emplace(std::piecewise_construct, std::forward_as_tuple(c),
                              std::forward_as_tuple(c, std::forward<Args>(args)...));
  it->second.var.assign(n, v);
}

}  // namespace

// ---- cubic_langmuir ------------------------------------------------------------

// Hand values with the defaults (the paper's adsorption branch), worked
// independently in Python with R = 8.31446261815324.
TEST(HumidCalf20, CubicLangmuirHandValues) {
  const CubicLangmuirIsothermParameters ips;
  const std::vector<std::tuple<double, double, double>> cases = {
      {295.15, 0.1, 1.0440147055608748}, {295.15, 0.5, 7.136870366311029},  {295.15, 1.0, 8.175142871766676},
      {323.15, 2.0, 5.525058304516076},  {348.15, 5.0, 4.728444878809606},  {373.15, 10.0, 4.088209276911139},
      {295.15, -0.01, -0.048120159590595096}};
  for (const auto& [T, C, q] : cases) EXPECT_LT(Rel(Cubic(ips, T, C), q), 1e-12) << T << " K, " << C << " mol/m3";
}

// The loading is the Henry line below C = 0 (C1 at 0) and is floored at 0
// where b0_1 < 0 makes the raw form negative.
TEST(HumidCalf20, CubicLangmuirContinuationAndFloor) {
  const CubicLangmuirIsothermParameters ips;
  const double T = 295.15;
  const auto k = CubicLangmuirConstantsAt(T, ips);
  const double slope = ips.q_sb * k.b1 + ips.q_sd * k.d;
  EXPECT_DOUBLE_EQ(Cubic(ips, T, 0.0), 0.0);
  EXPECT_NEAR(Cubic(ips, T, 1e-9) / 1e-9, slope, 1e-6 * std::fabs(slope));
  EXPECT_DOUBLE_EQ(Cubic(ips, T, -1e-3), -1e-3 * slope);

  CubicLangmuirIsothermParameters neg = ips;
  neg.q_sd = 0.0;  // the anti-Langmuir start alone: negative raw loading at small C
  const auto kn = CubicLangmuirConstantsAt(T, neg);
  ASSERT_LT(CubicLangmuirRawLoading(0.05, neg.q_sb, neg.q_sd, kn), 0.0);
  EXPECT_EQ(Cubic(neg, T, 0.05), 0.0);
  EXPECT_GT(Cubic(neg, T, 0.8), 0.0);
}

// The denominator minimum is exact (stationary points) against a fine scan,
// and positive for the paper's parameters up to saturation at 22-100 C.
TEST(HumidCalf20, CubicLangmuirDenominatorMinimum) {
  CubicLangmuirIsothermParameters ips;
  for (double T : {295.15, 323.15, 348.15, 373.15}) {
    const double C_sat = SaturationPressure(kWater, T) / (R * T);
    const double m = CubicLangmuirDenominatorMin(ips, T, C_sat);
    EXPECT_GT(m, 0.0) << T;
    const auto k = CubicLangmuirConstantsAt(T, ips);
    double scan = 1.0;
    for (int i = 0; i <= 200000; ++i) {
      const double C = C_sat * i / 200000.0;
      scan = std::min(scan, 1.0 + C * (k.b1 + C * (k.b2 + C * k.b3)));
    }
    EXPECT_NEAR(m, scan, 1e-9) << T;
  }
  ips.b0_2 = 0.0;
  ips.b0_3 = 1e-19;  // b1 C dominates: the denominator crosses zero
  double C_at = 0.0;
  EXPECT_LT(CubicLangmuirDenominatorMin(ips, 295.15, 1.1, &C_at), 0.0);
  EXPECT_GT(C_at, 0.0);
}

// IAST: the closed spreading pressure equals a numeric quadrature of q/p, and
// a binary with an inert gas gives the pure loading back.
TEST(HumidCalf20, CubicLangmuirSpreadingPressureAndIast) {
  CubicLangmuirIsotherm iso(kWater);
  for (double T : {295.15, 323.15, 373.15}) {
    const double p_sat = SaturationPressure(kWater, T);
    for (double x : {0.05, 0.2, 0.5, 0.9}) {
      const double P = x * p_sat;
      // Simpson in u = ln p over [ln P - 30, ln P]; below, q = K p adds q(p_min).
      const int n = 20000;
      const double u1 = std::log(P), u0 = u1 - 30.0, h = (u1 - u0) / n;
      double s = 0.0;
      for (int i = 0; i <= n; ++i) {
        const double w = (i == 0 || i == n) ? 1.0 : (i % 2 ? 4.0 : 2.0);
        s += w * EvalIsotherm(&iso, kWater, T, std::exp(u0 + i * h), 1.0);
      }
      const double quad = s * h / 3.0 + EvalIsotherm(&iso, kWater, T, std::exp(u0), 1.0);
      const double psi = ReducedSpreadingPressure(&iso, kWater, T, P);
      EXPECT_LT(Rel(psi, quad), 1e-9) << T << " K, x " << x;
    }
  }
  InertIsotherm inert(kN2);
  const double T = 295.15, P = 1e5, y = 0.5 * SaturationPressure(kWater, T) / P;
  std::vector<double> q;
  const auto r = EvalIastMixture({&iso, &inert}, {kWater, kN2}, T, P, {y, 1.0 - y}, q);
  ASSERT_TRUE(r.converged);
  EXPECT_LT(Rel(q[0], EvalIsotherm(&iso, kWater, T, P, y)), 1e-9);
  EXPECT_EQ(q[1], 0.0);
}

// The paper's pure-water data on the structured sorbent (SI of Nguyen et al.,
// Ind. Eng. Chem. Res. 63 (2024) 3265, Table S4: T [C], RH [%], q [mmol/g of
// sorbent], gravimetric; p = RH/100 p_sat(T)), against the defaults (the
// paper's adsorption branch). The paper's fit reaches RMSE 0.21 mmol/g over
// all temperatures (with the engine's p_sat 0.204), worst at 22 C (0.34).
TEST(HumidCalf20, CubicLangmuirReproducesTableS4) {
  struct Row { double t_c, rh, q; };
  const std::vector<Row> rows = {
      {22, 0.00, 0.00},  {22, 8.00, 0.68},  {22, 15.30, 2.27}, {22, 32.00, 5.74}, {22, 50.00, 6.91},
      {22, 75.00, 7.65}, {22, 79.00, 7.70}, {22, 100.00, 8.30},
      {30, 5.10, 0.38},  {30, 9.76, 0.97},  {30, 20.42, 3.26}, {30, 31.90, 5.61}, {30, 47.85, 6.64},
      {30, 50.40, 6.68}, {30, 63.80, 7.35},
      {40, 2.93, 0.22},  {40, 5.61, 0.51},  {40, 11.73, 1.31}, {40, 18.33, 2.74}, {40, 27.49, 4.49},
      {40, 28.96, 4.67}, {40, 36.65, 5.73},
      {50, 1.75, 0.13},  {50, 3.35, 0.30},  {50, 7.00, 0.70},  {50, 10.94, 1.23}, {50, 16.40, 2.11},
      {50, 17.28, 2.34}, {50, 21.87, 3.16},
      {60, 1.08, 0.06},  {60, 2.07, 0.21},  {60, 4.32, 0.41},  {60, 6.75, 0.68},  {60, 10.13, 1.05},
      {60, 10.67, 1.10}, {60, 13.51, 1.58},
      {80, 0.45, 0.03},  {80, 0.86, 0.09},  {80, 1.80, 0.20},  {80, 2.82, 0.29},  {80, 4.23, 0.39},
      {80, 4.45, 0.37},  {80, 5.64, 0.49},
      {100, 0.21, 0.01}, {100, 0.40, 0.06}, {100, 0.83, 0.14}, {100, 1.30, 0.19}, {100, 1.95, 0.21},
      {100, 2.06, 0.19}, {100, 2.61, 0.21}};
  ASSERT_EQ(rows.size(), 50u);
  CubicLangmuirIsotherm iso(kWater);
  double sse = 0.0, sse22 = 0.0;
  int n22 = 0;
  for (const auto& r : rows) {
    const double T = r.t_c + 273.15;
    const double p = r.rh / 100.0 * SaturationPressure(kWater, T);
    const double e = EvalIsotherm(&iso, kWater, T, p, 1.0) - r.q;  // mol/kg = mmol/g
    sse += e * e;
    if (r.t_c == 22) { sse22 += e * e; ++n22; }
  }
  const double rmse = std::sqrt(sse / rows.size()), rmse22 = std::sqrt(sse22 / n22);
  EXPECT_LT(rmse, 0.25) << "all temperatures";
  EXPECT_LT(rmse22, 0.40) << "22 C";
  std::cout << "Table S4: RMSE " << rmse << " mmol/g (all), " << rmse22 << " (22 C)\n";
}

// The isosteric heat from the isotherm's temperature dependence is finite at
// every humidity. It is negative at low humidity, where the anti-Langmuir start
// (b0_1 < 0, which grows in magnitude as T falls) dominates, and positive,
// up to about 87 kJ/mol, above about 10 % RH (a property of the published set).
TEST(HumidCalf20, CubicLangmuirIsostericHeat) {
  CubicLangmuirIsotherm iso(kWater);
  const double T = 295.15, p_sat = SaturationPressure(kWater, T);
  for (double x : {0.02, 0.05, 0.1, 0.3, 0.5, 0.9}) {
    const double qst = IsostericHeat(&iso, kWater, T, x * p_sat, 1.0);
    EXPECT_TRUE(std::isfinite(qst)) << x;
    if (x >= 0.15) {
      EXPECT_GT(qst, 4.0e4) << x;
      EXPECT_LT(qst, 1.0e5) << x;
    }
  }
  EXPECT_LT(IsostericHeat(&iso, kWater, T, 0.05 * p_sat, 1.0), 0.0);
}

// ---- the dual-site Langmuir water factors --------------------------------------

// The paper's modified dual-site Langmuir at 22 C in CO2 (0.95 bar) at the
// seven measured water loadings: hand values (Python, same R).
TEST(HumidCalf20, DslWaterFactorsHandValues) {
  const auto ips = Calf20CO2(true);
  const std::vector<std::pair<double, double>> cases = {
      {0.137, 3.0065978722834403}, {0.22, 3.0001814472547452}, {0.45, 2.9195082289340863},
      {1.07, 2.4656430152326543},  {2.08, 2.1499961176555225}, {7.85, 0.21264957613732155},
      {9.01, 0.13515711552461565}};
  for (const auto& [qw, q] : cases) {
    double v = 0.0;
    IsothermFunction::DualSiteLangmuir_Concentration(v, 295.15, 0.95e5, 1.0, qw, ips);
    EXPECT_LT(Rel(v, q), 1e-12) << qw;
  }
  double dry = 0.0;
  IsothermFunction::DualSiteLangmuir_Concentration(dry, 295.15, 0.95e5, 1.0, ips);
  EXPECT_LT(Rel(dry, 3.0080817645502913), 1e-12);
}

// The competitive CO2 loadings measured in pure CO2 at 0.97 bar with water
// (Nguyen et al., ibid., Table 4; T [C], RH [%], q_H2O, q_CO2 [mmol/g]): the
// dry dual-site Langmuir of the SI (Table S2) at each point's temperature and
// CO2 partial pressure, with the paper's per-site factors at the measured
// water loading. The paper reports RMSE 0.08 mol/kg; the bound is 0.1.
TEST(HumidCalf20, DslWaterFactorsReproduceTable4) {
  struct Row { double t_c, rh, q_w, q_co2; };
  const std::vector<Row> rows = {{22.2, 13, 0.137, 2.95}, {22.5, 18, 0.22, 2.93}, {24.2, 23, 0.45, 2.90},
                                 {22.7, 32, 1.07, 2.51},  {23.1, 47, 2.08, 2.00}, {23.1, 72, 7.85, 0.26},
                                 {22.5, 87, 9.01, 0.03}};
  const auto ips = Calf20CO2(true);
  const double P = 0.97e5;
  double sse = 0.0;
  for (const auto& r : rows) {
    const double T = r.t_c + 273.15;
    const double y_w = r.rh / 100.0 * SaturationPressure(kWater, T) / P;
    double q = 0.0;
    IsothermFunction::DualSiteLangmuir_Concentration(q, T, P, 1.0 - y_w, r.q_w, ips);
    sse += (q - r.q_co2) * (q - r.q_co2);
  }
  const double rmse = std::sqrt(sse / rows.size());
  EXPECT_LT(rmse, 0.1);
  std::cout << "Table 4: RMSE " << rmse << " mmol/g\n";
  // The dry loading at 22 C and 0.97 bar: the paper's text gives 3.00 mmol/g.
  double dry = 0.0;
  IsothermFunction::DualSiteLangmuir_Concentration(dry, 295.15, P, 1.0, ips);
  EXPECT_NEAR(dry, 3.00, 0.05);
}

// All factors 0: every path is the dry isotherm, bit for bit, at any q_w.
TEST(HumidCalf20, DslZeroFactorsAreTheDryKindExactly) {
  const auto c = Calf20CO2(false);
  DualSiteLangmuirIsothermParameters p;
  p.m1 = 2.0; p.b1 = 3e-10; p.q1 = 3e4; p.m2 = 1.0; p.b2 = 2e-11; p.q2 = 2e4;
  for (double T : {280.0, 300.0, 350.0})
    for (double P : {1e3, 1e5})
      for (double qw : {0.0, 0.5, 5.0}) {
        double a = 0.0, b = 0.0;
        IsothermFunction::DualSiteLangmuir_Concentration(a, T, P, 0.3, c);
        IsothermFunction::DualSiteLangmuir_Concentration(b, T, P, 0.3, qw, c);
        EXPECT_EQ(a, b);
        IsothermFunction::DualSiteLangmuir(a, T, P, 0.3, p);
        IsothermFunction::DualSiteLangmuir(b, T, P, 0.3, qw, p);
        EXPECT_EQ(a, b);
      }
  // lambda3 alone changes nothing.
  auto c3 = c;
  c3.lambda3 = 3.55;
  EXPECT_FALSE(DslWaterFactors::Active(c3));
}

// The extended (competitive) kind: CO2 takes its factors, N2 (factors 0)
// keeps the shared denominator, which carries CO2's affinity factor. Update
// with the adsorbent's loadings matches the hand formula cell by cell; with
// zero factors Update is the dry pass bit for bit.
TEST(HumidCalf20, ExtendedDslWithWaterFactorsUpdate) {
  const size_t n = 4;
  ComponentEquation yi, q_map;
  ComponentIndependentEquation qi_sat, dry_sat;
  std::vector<double> T(n), P(n);
  for (const std::string& c : {kCO2, kN2, kWater}) {
    AddCells(yi, c, n, 0.0, CellType::CellCentre, EquationType::Differential, [](UserData*) {}, 1e-8);
    AddCells(q_map, c, n, 0.0, CellType::CellCentre, EquationType::Differential, [](UserData*) {}, 1e-8);
    AddCells(qi_sat, c, n, -1.0, CellType::CellCentre, [](UserData*) {});
    AddCells(dry_sat, c, n, -1.0, CellType::CellCentre, [](UserData*) {});
  }
  for (size_t i = 0; i < n; ++i) {
    T[i] = 295.15 + 10.0 * i;
    P[i] = 1.0e5;
    yi.at(kCO2).var[i] = 0.15;
    yi.at(kN2).var[i] = 0.8;
    yi.at(kWater).var[i] = 0.05;
    q_map.at(kWater).var[i] = 1.5 * i - 0.5;  // cell 0 negative: dry
  }
  ExtendedDualSiteLangmuirIsotherm_Concentration iso(kCO2);
  iso.params_[kCO2] = Calf20CO2(true);
  iso.params_[kCO2].gamma2 = 0.2;
  auto& n2 = iso.params_[kN2];
  n2.m1 = 2.387; n2.b1 = 8.139e-7; n2.q1 = 17960.0; n2.m2 = 3.271; n2.b2 = 8.139e-7; n2.q2 = 17960.0;

  IsothermUpdateParams with(yi, qi_sat, T, P, kCO2, 0, static_cast<int>(n) - 1, nullptr, &q_map);
  iso.UpdateAll(with, {kCO2, kN2});
  const auto& c = iso.params_[kCO2];
  for (size_t i = 0; i < n; ++i) {
    const double qw = std::max(q_map.at(kWater).var[i], 0.0);
    const double RT = R * T[i];
    const double cc = P[i] * 0.15 / RT, cn = P[i] * 0.8 / RT;
    const double aff = std::exp(-0.2 * qw), x = std::pow(qw, 3.55);
    const double B1c = 5.519e-7 * std::exp(35060.0 / RT) * aff, B2c = 5.187e-8 * std::exp(28950.0 / RT) * aff;
    const double Bn = 8.139e-7 * std::exp(17960.0 / RT);
    const double d1 = 1 + B1c * cc + Bn * cn, d2 = 1 + B2c * cc + Bn * cn;
    const double q_co2 = 2.387 * B1c * cc / d1 / (1 + 0.0066 * x) + 3.271 * B2c * cc / d2 / (1 + 2.483 * x);
    const double q_n2 = 2.387 * Bn * cn / d1 + 3.271 * Bn * cn / d2;
    EXPECT_LT(Rel(qi_sat.at(kCO2).var[i], q_co2), 1e-12) << i;
    EXPECT_LT(Rel(qi_sat.at(kN2).var[i], q_n2), 1e-12) << i;
  }
  (void)c;

  // Zero factors: with or without the loadings, the dry pass bit for bit.
  iso.params_[kCO2] = Calf20CO2(false);
  IsothermUpdateParams without(yi, dry_sat, T, P, kCO2, 0, static_cast<int>(n) - 1);
  iso.UpdateAll(without, {kCO2, kN2});
  iso.UpdateAll(with, {kCO2, kN2});
  for (size_t i = 0; i < n; ++i)
    for (const std::string& comp : {kCO2, kN2}) EXPECT_EQ(qi_sat.at(comp).var[i], dry_sat.at(comp).var[i]);
}

// Kind layout, registry, deck round trip; a deck whose factors are 0 writes no
// factor fields, so it saves exactly as before.
TEST(HumidCalf20, KindsRegistryAndDeckRoundTrip) {
  IsothermType t{};
  ASSERT_TRUE(KindFromString("cubic_langmuir", t));
  EXPECT_EQ(t, IsothermType::CUBIC_LANGMUIR);
  EXPECT_EQ(static_cast<int>(IsothermType::CUBIC_LANGMUIR), 25);  // appended
  EXPECT_TRUE(IsSupportedIsothermKind(t));
  EXPECT_FALSE(IsCompetitiveIsothermKind(t));
  EXPECT_FALSE(IsWaterCoupledIsothermKind(t));
  EXPECT_NE(IsothermKind(t).iast, IastFamily::None);
  for (auto k : {IsothermType::DUAL_SITE_LANGMUIR_C, IsothermType::EXTENDED_DUAL_SITE_LANGMUIR_P}) {
    const auto& ps = IsothermKind(k).parameters;
    ASSERT_GE(ps.size(), 10u);
    EXPECT_EQ(ps[6].name, "lambda1");
    EXPECT_EQ(ps[9].name, "gamma2");
    EXPECT_EQ(ps[9].default_value, 0.0);
  }

  Simulation sim;
  UserData& d = sim.Data();
  build::AddComponent(d, kCO2);
  build::AddComponent(d, kN2);
  build::AddComponent(d, kWater);
  const AdsorbentID a = build::AddAdsorbent(d, "CALF-20");
  build::SetIsotherm(d, a, "extended_dual_site_langmuir_c", {kCO2, kN2});
  build::SetIsotherm(d, a, "cubic_langmuir", {kWater});
  sim.RefreshParameters();

  std::string dry_text;
  SaveToString(d, dry_text);
  EXPECT_EQ(dry_text.find("IP_lambda1"), std::string::npos);
  EXPECT_NE(dry_text.find("IP_q_sb"), std::string::npos);

  const std::string base = "adsorbent[CALF-20].isotherm[CarbonDioxide].";
  sim.SetParameter(base + "lambda1", 0.0066);
  sim.SetParameter(base + "lambda2", 2.483);
  sim.SetParameter(base + "lambda3", 3.55);
  sim.SetParameterText("adsorbent[CALF-20].isotherm[Water].dU1", "-97.99 [kJ/mol]");
  std::string text;
  SaveToString(d, text);
  EXPECT_NE(text.find("IP_lambda1"), std::string::npos);
  UserData back;
  LoadDeckFlexible(back, text);
  auto* iso = dynamic_cast<ExtendedDualSiteLangmuirIsotherm_Concentration*>(
      back.adsorbent_directory.Get(a).adsorption_source_.GetIsotherm(kCO2).get());
  ASSERT_NE(iso, nullptr);
  EXPECT_DOUBLE_EQ(iso->params_[kCO2].lambda2.Get(), 2.483);
  EXPECT_DOUBLE_EQ(iso->params_[kCO2].lambda3.Get(), 3.55);
  EXPECT_DOUBLE_EQ(iso->params_[kN2].lambda1.Get(), 0.0);  // N2 wrote none: read back as 0
  auto* w = dynamic_cast<CubicLangmuirIsotherm*>(
      back.adsorbent_directory.Get(a).adsorption_source_.GetIsotherm(kWater).get());
  ASSERT_NE(w, nullptr);
  EXPECT_DOUBLE_EQ(w->params_[kWater].dU1.Get(), -97990.0);

  const Parameter* p = sim.Parameters().Find("adsorbent[" + std::to_string(a) + "].isotherm[Water].b0_2");
  ASSERT_NE(p, nullptr);
  EXPECT_EQ(p->Info().unit, "(m3/mol)^2");
  EXPECT_TRUE(IsothermUsesWaterLoading(iso, kCO2));
  EXPECT_FALSE(IsothermUsesWaterLoading(iso, kN2));
}

// Every bundled deck that uses a dual-site Langmuir kind saves without any
// water factor field: the decks are unchanged by the new parameters.
TEST(HumidCalf20, BundledDecksWriteNoWaterFactors) {
  const std::filesystem::path dir = std::filesystem::current_path() / "../examples/decks";  // ctest runs in core/
  if (!std::filesystem::exists(dir)) GTEST_SKIP() << "no example decks at " << dir;
  int checked = 0;
  for (const auto& e : std::filesystem::directory_iterator(dir)) {
    if (e.path().extension() != ".json") continue;
    std::ifstream in(e.path());
    std::stringstream ss;
    ss << in.rdbuf();
    const std::string original = ss.str();
    if (original.find("DualSiteLangmuirIsotherm") == std::string::npos) continue;
    UserData data;
    LoadDeckFlexible(data, original);
    std::string saved;
    SaveToString(data, saved);
    EXPECT_EQ(saved.find("IP_lambda1"), std::string::npos) << e.path().filename();
    ++checked;
  }
  EXPECT_GT(checked, 0);
}

// The deck checks: IAST refuses a dual-site Langmuir with water factors (and
// accepts it with the factors at 0), a missing water component is a warning,
// and a Cubic-Langmuir whose denominator crosses zero below saturation at the
// deck's temperatures is an error.
TEST(HumidCalf20, DeckChecks) {
  Simulation sim;
  UserData& d = sim.Data();
  build::AddComponent(d, kCO2);
  build::AddComponent(d, kN2);
  const AdsorbentID a = build::AddAdsorbent(d, "CALF-20");
  build::SetIsotherm(d, a, "extended_dual_site_langmuir_c", {kCO2, kN2});
  sim.RefreshParameters();
  MixedAdsorbentLayer layer;
  layer.name_ = "L";
  layer.AddAdsorbent(a, 1.0);
  using S = LayerCompositionIssue::Severity;
  auto has = [&](S s, const char* text) {
    for (const auto& i : CheckLayerComposition(layer, d.adsorbent_directory))
      if (i.severity == s && i.message.find(text) != std::string::npos) return true;
    return false;
  };
  build::SetMixingRule(d, a, "iast");
  EXPECT_FALSE(has(S::Error, "cannot mix"));
  sim.SetParameter("adsorbent[CALF-20].isotherm[CarbonDioxide].lambda2", 2.483);
  EXPECT_TRUE(has(S::Error, "cannot mix"));
  EXPECT_TRUE(has(S::Warning, "has no such component"));
  IIsotherm* iso = d.adsorbent_directory.Get(a).adsorption_source_.GetIsotherm(kCO2).get();
  EXPECT_THROW(ResolvePureIsotherm(iso, kCO2), std::invalid_argument);
  build::SetMixingRule(d, a, "none");

  build::AddComponent(d, kWater);
  build::SetIsotherm(d, a, "cubic_langmuir", {kWater});
  sim.RefreshParameters();
  EXPECT_FALSE(has(S::Warning, "has no such component"));
  std::string err;
  EXPECT_TRUE(d.CubicLangmuirChecks(err)) << err;
  sim.SetParameter("adsorbent[CALF-20].isotherm[Water].b0_2", 0.0);
  sim.SetParameter("adsorbent[CALF-20].isotherm[Water].b0_3", 1e-19);
  err.clear();
  EXPECT_FALSE(d.CubicLangmuirChecks(err));
  EXPECT_NE(err.find("denominator"), std::string::npos) << err;
}

// Fitter: Cubic-Langmuir data from the defaults over 22-100 C are recovered;
// dry dual-site data hold the water factors at 0, humid data fit them.
TEST(HumidCalf20, Fitter) {
  {
    const auto& kinds = isofit::FittableKinds();
    EXPECT_NE(std::find(kinds.begin(), kinds.end(), IsothermType::CUBIC_LANGMUIR), kinds.end());
    isofit::Model truth(IsothermType::CUBIC_LANGMUIR, kWater);
    std::vector<isofit::Point> pts;
    for (double T : {295.15, 323.15, 348.15, 373.15}) {
      const double ps = SaturationPressure(kWater, T);
      for (double x : {0.03, 0.08, 0.15, 0.25, 0.35, 0.5, 0.65, 0.8, 0.9})
        pts.push_back({x * ps, truth.Eval(T, x * ps, 1.0), T, 1.0, 0.0});
    }
    std::vector<isofit::ParamSpec> specs;
    double t_ref = 0.0;
    std::string err;
    ASSERT_TRUE(isofit::DefaultSpecs(IsothermType::CUBIC_LANGMUIR, pts, specs, t_ref, err)) << err;
    EXPECT_EQ(specs.size(), truth.Names().size());
    // Start at the truth: checks the set-up and the sign of the pairing.
    for (auto& s : specs) s.value = truth.Get(s.name);
    isofit::FitOptions opt;
    opt.max_global_evals = 1;
    const auto r = isofit::Fit(IsothermType::CUBIC_LANGMUIR, pts, specs, t_ref, opt);
    ASSERT_TRUE(r.ok) << r.error;
    EXPECT_LT(r.rmse, 1e-6);
  }
  {
    isofit::Model truth(IsothermType::DUAL_SITE_LANGMUIR_C);
    for (const auto& [n, v] : std::vector<std::pair<std::string, double>>{
             {"m1", 2.387}, {"b1", 5.519e-7}, {"q1", 35060.0}, {"m2", 3.271}, {"b2", 5.187e-8}, {"q2", 28950.0},
             {"lambda1", 0.0066}, {"lambda2", 2.483}, {"lambda3", 3.55}})
      truth.Set(n, v);
    std::vector<isofit::Point> dry, wet;
    for (double p : {1e3, 5e3, 2e4, 5e4, 1e5}) {
      dry.push_back({p, truth.Eval(295.15, p, 1.0, 0.0), 295.15, 1.0, 0.0});
      for (double qw : {0.45, 1.07, 2.08, 7.85}) wet.push_back({p, truth.Eval(295.15, p, 1.0, qw), 295.15, 1.0, qw});
    }
    std::vector<isofit::ParamSpec> specs;
    double t_ref = 0.0;
    std::string err;
    ASSERT_TRUE(isofit::DefaultSpecs(IsothermType::DUAL_SITE_LANGMUIR_C, dry, specs, t_ref, err)) << err;
    EXPECT_EQ(specs.size(), truth.Names().size());
    for (const auto& s : specs)
      if (s.name.rfind("lambda", 0) == 0 || s.name == "gamma2") {
        EXPECT_TRUE(s.fixed) << s.name;
        EXPECT_EQ(s.value, 0.0) << s.name;
      }
    ASSERT_TRUE(isofit::DefaultSpecs(IsothermType::DUAL_SITE_LANGMUIR_C, wet, specs, t_ref, err)) << err;
    for (auto& s : specs) {
      if (s.name == "lambda1" || s.name == "lambda2" || s.name == "lambda3") {
        EXPECT_FALSE(s.fixed) << s.name;
      } else {
        s.fixed = true;
        s.value = truth.Get(s.name);
        s.pair_energy.clear();
      }
    }
    isofit::FitOptions opt;
    opt.max_global_evals = 200;
    const auto r = isofit::Fit(IsothermType::DUAL_SITE_LANGMUIR_C, wet, specs, t_ref, opt);
    ASSERT_TRUE(r.ok) << r.error;
    EXPECT_LT(r.rmse, 1e-4);
    for (const auto& pr : r.params)
      if (pr.name == "lambda2") EXPECT_NEAR(pr.value, 2.483, 0.01);
  }
}

// A small humid VSA cycle on CALF-20: the CO2 dual-site Langmuir with the
// paper's water factors, N2 competing in the same extended isotherm, water on
// cubic_langmuir with the paper's water LDF (0.005 1/s, SI Table S3). Feed of
// 15 % CO2 at 30 % relative humidity and 22 C, evacuation to 0.1 bar, feed
// pressurisation, from a bed in equilibrium with the feed humidity.
static void BuildHumidCycle(Simulation& sim) {
  UserData& d = sim.Data();
  for (const auto& c : {kCO2, kN2, kWater}) build::AddComponent(d, c);
  const AdsorbentID ads = build::AddAdsorbent(d, "CALF-20");
  build::SetIsotherm(d, ads, "extended_dual_site_langmuir_c", {kCO2, kN2});
  build::SetIsotherm(d, ads, "cubic_langmuir", {kWater});
  for (const auto& c : {kCO2, kN2, kWater}) build::SetRateModel(d, ads, c, "ldf", "constant");
  build::SetBedShape(d, "cylinder");
  build::SetWallMaterial(d, "ss316l");
  build::SetWallBoundary(d, "zero_flux");
  build::SetThermalEquilibrium(d, "lte");
  const size_t layer = build::AddLayer(d, "Bed");
  build::SetLayerComposition(d, layer, {{ads, 1.0}});
  build::SetMomentum(d, layer, "kozeny_carman");
  build::SetDispersion(d, layer, "constant");
  build::SetHeatDispersion(d, layer, "constant");
  build::SetInternalHeatTransfer(d, layer, "constant");
  build::SetExternalHeatTransfer(d, layer, "constant");
  const size_t feed = build::AddStep(d, "Feed");
  build::SetFlowBoundary(d, feed, "bottom", "inlet");
  const size_t feed_in = build::AddStream(d, feed, "bottom", "Feed inlet", "velocity", "constant");
  build::SetPressureBoundary(d, feed, "top", "outlet", "constant");
  const size_t evac = build::AddStep(d, "Evacuation");
  const size_t evac_end = build::SetPressureBoundary(d, evac, "bottom", "outlet", "exponential");
  const size_t press = build::AddStep(d, "Pressurisation");
  const size_t press_end = build::SetPressureBoundary(d, press, "bottom", "inlet", "linear");
  build::SetKPI(d, kCO2, {feed_in, press_end}, {evac_end});
  sim.RefreshParameters();

  const double T0 = 295.15;
  const double y_w = 0.3 * SaturationPressure(kWater, T0) / 1e5;
  auto num = [](double v) { std::ostringstream ss; ss.precision(17); ss << v; return ss.str(); };
  const std::string co2 = "adsorbent[CALF-20].isotherm[CarbonDioxide].";
  const std::string n2 = "adsorbent[CALF-20].isotherm[Nitrogen].";
  const std::vector<std::pair<std::string, std::string>> values = {
      {"bed.diameter", "0.021"},
      {"wall.thickness", "0.002"},
      {"adsorbent[CALF-20].solid_density", "570"},
      {"adsorbent[CALF-20].heat_capacity", "1350"},
      {"adsorbent[CALF-20].thermal_conductivity", "0.63"},
      {"adsorbent[CALF-20].particle_void_fraction", "0.35"},
      {co2 + "m1", "2.387"}, {co2 + "b1", "5.519e-7"}, {co2 + "q1", "35060"},
      {co2 + "m2", "3.271"}, {co2 + "b2", "5.187e-8"}, {co2 + "q2", "28950"},
      {co2 + "lambda1", "0.0066"}, {co2 + "lambda2", "2.483"}, {co2 + "lambda3", "3.55"},
      {n2 + "m1", "2.387"}, {n2 + "b1", "8.139e-7"}, {n2 + "q1", "17960"},
      {n2 + "m2", "3.271"}, {n2 + "b2", "8.139e-7"}, {n2 + "q2", "17960"},
      {"adsorbent[CALF-20].kinetics[CarbonDioxide].k", "0.1"},
      {"adsorbent[CALF-20].kinetics[Nitrogen].k", "0.5"},
      {"adsorbent[CALF-20].kinetics[Water].k", "0.005"},  // SI Table S3
      {"adsorbent[CALF-20].heat_of_adsorption[CarbonDioxide]", "38000"},
      {"adsorbent[CALF-20].heat_of_adsorption[Nitrogen]", "15000"},
      {"adsorbent[CALF-20].heat_of_adsorption[Water]", "50000"},
      {"layer[Bed].length", "0.1"},
      {"layer[Bed].cells", "10"},
      {"layer[Bed].bed_void_fraction", "0.4"},
      {"layer[Bed].particle_diameter", "0.001"},
      {"layer[Bed].dispersion_coefficient", "1e-4"},
      {"initial.pressure", "1e5"},
      {"initial.temperature", num(T0)},
      {"initial.wall_temperature", num(T0)},
      {"initial.y[Water]", num(y_w)},  // the bed starts wet, in equilibrium with the feed humidity
      {"initial.y[Nitrogen]", num(1.0 - y_w)},
      {"step[Feed].time", "60"},
      {"step[Feed].bottom[Feed inlet].velocity", "0.1"},
      {"step[Feed].bottom[Feed inlet].temperature", num(T0)},
      {"step[Feed].bottom[Feed inlet].y[CarbonDioxide]", "0.15"},
      {"step[Feed].bottom[Feed inlet].y[Water]", num(y_w)},
      {"step[Feed].bottom[Feed inlet].y[Nitrogen]", num(0.85 - y_w)},
      {"step[Feed].top.pressure", "1e5"},
      {"step[Evacuation].time", "60"},
      {"step[Evacuation].bottom.pressure_exponential_end", "1e4"},
      {"step[Pressurisation].time", "20"},
      {"step[Pressurisation].bottom.pressure_linear_end", "1e5"},
      {"step[Pressurisation].bottom.temperature", num(T0)},
      {"step[Pressurisation].bottom.y[CarbonDioxide]", "0.15"},
      {"step[Pressurisation].bottom.y[Water]", num(y_w)},
      {"step[Pressurisation].bottom.y[Nitrogen]", num(0.85 - y_w)},
  };
  for (const auto& [path, v] : values) sim.SetParameterText(path, v);
  // The molar-profile test (the end-of-feed gas profile): the per-component
  // balance of water converges only over hundreds of cycles here, because
  // the water held on the sorbent is about 10^3 times the water in the gas.
}

// Fast tier: three cycles run, water is on the sorbent and the CO2
// equilibrium follows the adsorbed water.
TEST(HumidCalf20, SmallHumidCycleRuns) {
  Simulation sim;
  BuildHumidCycle(sim);
  UserData& d = sim.Data();
  d.settings.css.method_ = CSSMethod::NONE;
  d.cycle.number_of_cycles_ = 3;
  std::string error;
  ASSERT_TRUE(sim.Check(error)) << error;
  const RunResult r = sim.Run(RunOptions{});
  ASSERT_TRUE(r.Ok()) << r.message;
  EXPECT_EQ(r.cycles_run, 3);
  const auto qw = sim.Series("adsorbed_concentration", kWater, SeriesLocation::Bottom).values;
  ASSERT_FALSE(qw.empty());
  EXPECT_GT(qw.back(), 0.5);  // about 6 mol/kg at 30 % RH, less after the evacuation
  for (double v : sim.Series("adsorbed_concentration", kCO2, SeriesLocation::Bottom).values)
    EXPECT_TRUE(std::isfinite(v));
  // The CO2 equilibrium loading at the end of the run reads the adsorbed water.
  const auto ads = d.adsorbent_directory.GetAll().begin()->first;
  IIsotherm* iso = d.adsorbent_directory.Get(ads).adsorption_source_.GetIsotherm(kCO2).get();
  const auto& s = d.state;
  const int i = d.bed.GetFirstBedCellIndex();
  const double q_w = s.ads_concentration.at(ads).at(kWater).var[i];
  const double T = s.solid_temperature.var[i], P = s.pressure.var[i], y = s.mole_fraction.at(kCO2).var[i];
  ASSERT_GT(q_w, 0.5);
  EXPECT_LT(EvalIsotherm(iso, kCO2, T, P, y, q_w), EvalIsotherm(iso, kCO2, T, P, y));
}

// Slow tier (~11 s): the same cycle reaches cyclic steady state on the
// per-component mass balance (161 cycles: the water held on the sorbent is
// about 10^3 times the water in the gas, so its profile settles slowly).
TEST(HumidCalf20, SmallHumidCycleReachesCss) {
  Simulation sim;
  BuildHumidCycle(sim);
  UserData& d = sim.Data();
  d.settings.css.method_ = CSSMethod::MASS_BALANCE_COMPONENT;
  d.settings.css.min_cycles_ = 3;
  d.cycle.number_of_cycles_ = 300;
  std::string error;
  ASSERT_TRUE(sim.Check(error)) << error;
  const RunResult r = sim.Run(RunOptions{});
  ASSERT_TRUE(r.Ok()) << r.message;
  EXPECT_EQ(r.status, RunStatus::ReachedCSS) << r.message;
  std::cout << "humid cycle: CSS after " << r.cycles_run << " cycles\n";
  const auto qw = sim.Series("adsorbed_concentration", kWater, SeriesLocation::Bottom).values;
  ASSERT_FALSE(qw.empty());
  EXPECT_GT(qw.back(), 0.5);
}

}  // namespace adsim
