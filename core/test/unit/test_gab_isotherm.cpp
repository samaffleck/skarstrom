/*
 * Skarstrom - pressure and temperature swing adsorption simulator
 * Copyright (c) 2024-2026 Sam Affleck and contributors
 *
 * Licensed under the terms in the LICENSE file at the repository root.
 */
// GAB (Guggenheim-Anderson-de Boer) water isotherm kinds, gab and gab_young:
// closed-form values, the BET limit, the Henry slope at x -> 0, the
// continuation near the pole K x -> 1, the saturation pressure, an independent
// cross-check against pyDAC's GAB, the solver's Update path, the isosteric heat
// it implies, and the kind's plumbing (registry paths, deck round trip).

#include <cmath>
#include <string>
#include <vector>

#include <gtest/gtest.h>

#include "CoolProp.h"
#include "adsim/equations/equation.h"
#include "adsim/equations/isotherm.h"
#include "adsim/equations/isotherm_kinds.h"
#include "adsim/math/constants.h"
#include "adsim/objects/deck_builder.h"
#include "adsim/objects/deck_io.h"
#include "adsim/simulation.h"
#include "adsim/solver/isotherm_fit.h"

namespace adsim {
namespace {

const std::string kWater = "Water";

double Rel(double a, double b) { return std::fabs(a / b - 1.0); }

GabIsothermParameters HandParams() {
  GabIsothermParameters p;
  p.q_m0 = 2.5;
  p.dH_m = 500.0;
  p.C0 = 0.8;
  p.dH_C = 6000.0;
  p.K0 = 0.75;
  p.dH_K = 150.0;
  return p;
}

// q at relative pressure x (P = x p_sat, y = 1) through the public function.
double GabAt(const GabIsothermParameters& p, double x, double T) {
  const double ps = SaturationPressure(kWater, T);
  double q = 0.0;
  IsothermFunction::Gab(q, T, x * ps, 1.0, ps, p);
  return q;
}

}  // namespace

// Values worked independently (Python, R = 8.31446261815324):
//   q_m = 2.5 exp(500/RT), C = 0.8 exp(6000/RT), K = 0.75 exp(150/RT),
//   q = q_m C K x / ((1 - K x)(1 - K x + C K x)).
TEST(GabIsotherm, ClosedFormValues) {
  const auto p = HandParams();
  const GabConstants g = GabConstantsAt(298.15, p);
  EXPECT_NEAR(g.C, 9.0000414366, 1e-9);
  EXPECT_NEAR(g.K, 0.796783122821, 1e-11);
  EXPECT_NEAR(GabAt(p, 0.1, 298.15), 1.45552110761, 1e-10);
  EXPECT_NEAR(GabAt(p, 0.5, 298.15), 4.35370023938, 1e-10);
  EXPECT_NEAR(GabAt(p, 0.9, 323.15), 9.98368691404, 1e-10);
  EXPECT_NEAR(GabAt(p, 0.3, 273.15), 3.20088487523, 1e-10);
  // The relative pressure is P y / p_sat: the same partial pressure at another
  // total pressure gives the same loading.
  const double ps = SaturationPressure(kWater, 298.15);
  double q1 = 0.0, q2 = 0.0;
  IsothermFunction::Gab(q1, 298.15, 1.0e5, 0.5 * ps / 1.0e5, ps, p);
  IsothermFunction::Gab(q2, 298.15, 2.0e5, 0.25 * ps / 1.0e5, ps, p);
  EXPECT_NEAR(q1, 4.35370023938, 1e-10);
  EXPECT_NEAR(q2, q1, 1e-12);
}

// BET is GAB with K = 1: q = q_m C x / ((1 - x)(1 - x + C x)).
TEST(GabIsotherm, BetLimit) {
  EXPECT_NEAR(GabLoading(2.0, 10.0, 1.0, 0.5), 3.63636363636, 1e-10);
  for (double x : {0.01, 0.1, 0.3, 0.6, 0.9, 0.94}) {
    const double bet = 2.0 * 10.0 * x / ((1.0 - x) * (1.0 - x + 10.0 * x));
    EXPECT_NEAR(GabLoading(2.0, 10.0, 1.0, x), bet, 1e-12 * bet) << x;
  }
  GabIsothermParameters p;
  p.q_m0 = 2.0; p.dH_m = 0.0; p.C0 = 10.0; p.dH_C = 0.0; p.K0 = 1.0; p.dH_K = 0.0;
  EXPECT_NEAR(GabAt(p, 0.5, 300.0), 3.63636363636, 1e-10);
}

// At x -> 0 the GAB is linear with slope q_m C K; for x < 0 (a transiently
// negative mole fraction) it stays on that line, so q is C1 through zero.
TEST(GabIsotherm, HenrySlopeAtZero) {
  const double q_m = 2.5, C = 9.0, K = 0.8, A = q_m * C * K;
  EXPECT_EQ(GabLoading(q_m, C, K, 0.0), 0.0);
  EXPECT_NEAR(GabLoading(q_m, C, K, 1e-9) / 1e-9, A, 1e-7 * A);
  EXPECT_DOUBLE_EQ(GabLoading(q_m, C, K, -1e-3), -A * 1e-3);
  EXPECT_DOUBLE_EQ(GabLoading(q_m, C, K, -0.5), -A * 0.5);
  // Slope continuity through zero.
  const double h = 1e-7;
  const double right = (GabLoading(q_m, C, K, h) - GabLoading(q_m, C, K, 0.0)) / h;
  const double left = (GabLoading(q_m, C, K, 0.0) - GabLoading(q_m, C, K, -h)) / h;
  EXPECT_NEAR(right, left, 1e-5 * A);
}

// The pole K x = 1. With K <= 0.95 the GAB is exact up to x = 1 and continued
// linearly above; with K near or above 1 (BET) it is exact up to K x = 0.95.
// Either way q is finite, continuous, C1 and increasing.
TEST(GabIsotherm, ContinuationNearThePole) {
  struct Case { double q_m, C, K, x_c; };
  for (const Case c : {Case{2.0, 10.0, 1.0, 0.95}, Case{2.0, 10.0, 1.2, 0.95 / 1.2},
                       Case{3.2, 1.74, 0.84, 1.0}, Case{36.48, 0.1489, 0.5751, 1.0}}) {
    SCOPED_TRACE("K = " + std::to_string(c.K) + ", C = " + std::to_string(c.C));
    auto q = [&](double x) { return GabLoading(c.q_m, c.C, c.K, x); };
    const double exact_below = c.q_m * c.C * c.K * c.x_c /
                               ((1 - c.K * c.x_c) * (1 - c.K * c.x_c + c.C * c.K * c.x_c));
    EXPECT_NEAR(q(c.x_c), exact_below, 1e-12 * exact_below);
    const double h = 1e-6;
    EXPECT_NEAR(q(c.x_c + h), q(c.x_c - h), 1e-4 * exact_below);         // continuous
    const double s_left = (q(c.x_c) - q(c.x_c - h)) / h;
    const double s_right = (q(c.x_c + h) - q(c.x_c)) / h;
    EXPECT_NEAR(s_right, s_left, 1e-3 * s_left);                           // C1
    double prev = q(0.0);
    for (double x = 0.01; x <= 3.0; x += 0.01) {                            // finite, increasing
      const double v = q(x);
      ASSERT_TRUE(std::isfinite(v)) << x;
      ASSERT_GT(v, prev) << x;
      prev = v;
    }
    // Above x_c the continuation is the tangent line.
    EXPECT_NEAR(q(c.x_c + 0.5) - q(c.x_c), 0.5 * s_right, 1e-3 * 0.5 * s_right);
  }
  // K = 0.8: exact at saturation itself.
  const double at1 = 2.0 * 10.0 * 0.8 / ((1 - 0.8) * (1 - 0.8 + 8.0));
  EXPECT_NEAR(GabLoading(2.0, 10.0, 0.8, 1.0), at1, 1e-12 * at1);
}

// p_sat from CoolProp's saturation ancillary: within 2e-4 of the full
// Helmholtz-energy saturation pressure and of the IAPWS correlation pyDAC uses
// (Wagner and Pruss, 1993); p_c above T_c; NaN for an unknown fluid.
TEST(GabIsotherm, SaturationPressure) {
  // pyDAC water_saturation_pressure(T) at these temperatures (see the
  // cross-check below for its source).
  const std::vector<std::pair<double, double>> wagner = {
      {288.15, 1705.676858}, {298.15, 3169.824241}, {313.15, 7385.109965},
      {343.15, 31202.19682}, {373.15, 101417.9895}};
  for (const auto& [T, ps] : wagner) {
    EXPECT_LT(Rel(SaturationPressure(kWater, T), ps), 2e-4) << T;
    EXPECT_LT(Rel(SaturationPressure(kWater, T), CoolProp::PropsSI("P", "T", T, "Q", 1, kWater)), 2e-4) << T;
  }
  EXPECT_NEAR(SaturationPressure(kWater, 700.0), CoolProp::PropsSI("pcrit", "", 0, "", 0, kWater), 1.0);
  EXPECT_GT(SaturationPressure(kWater, 253.15), 0.0);  // below the triple point: extrapolated, finite
  EXPECT_TRUE(std::isnan(SaturationPressure("NotAFluid", 300.0)));
}

// Independent cross-check against pyDAC 3.1.0 (MIT, PyPI pydac3), whose H2O
// model is this layer-energy GAB (src/pydac/core/isotherm_models.py,
// GABIsotherm.isotherm, with its fitted q_m, C, D, F, G and E10 = 57220 -
// 44.38 T). The reference loadings were computed from pyDAC's equations
// transcribed into a separate script (pyDAC itself was not run), at the
// relative pressures shown; the engine is evaluated at the same relative
// pressure, so the comparison does not depend on the p_sat correlation.
TEST(GabIsotherm, CrossCheckAgainstPyDac) {
  GabYoungIsothermParameters p;
  p.q_m = 3.1970616886895318;
  p.C_E = 46779.95249302809;
  p.D = 0.02356068511609903;
  p.F = 56233.72466565154;
  p.G = -42.738276947397196;
  p.E10_A = 57220.0;
  p.E10_B = -44.38;
  const std::vector<double> xs = {0.05, 0.2, 0.4, 0.6, 0.8};
  struct Row { double T, c, k; std::vector<double> q; };
  const std::vector<Row> rows = {
      {288.15, 1.839332219, 0.8071768847, {0.239206642572, 0.996981339204, 2.20609497768, 3.92646191797, 6.95132518505}},
      {298.15, 1.959749883, 0.8183991607, {0.257218603605, 1.0597574693, 2.3202669092, 4.1085646286, 7.29595790045}},
      {313.15, 2.040649758, 0.8341454443, {0.272132957665, 1.1130972894, 2.42485206583, 4.2981812843, 7.72312287378}},
      {343.15, 1.718017565, 0.8622322066, {0.240033529764, 1.01845102548, 2.31772303432, 4.2926892387, 8.16787517627}},
      {373.15, 0.8625194259, 0.8865313486, {0.128685376066, 0.609148818468, 1.59279492548, 3.38083933155, 7.45250157088}}};
  for (const auto& r : rows) {
    const GabConstants g = GabConstantsAt(r.T, p);
    EXPECT_LT(Rel(g.C, r.c), 1e-9) << r.T;
    EXPECT_LT(Rel(g.K, r.k), 1e-9) << r.T;
    const double ps = SaturationPressure(kWater, r.T);
    for (size_t i = 0; i < xs.size(); ++i) {
      double q = 0.0;
      IsothermFunction::GabYoung(q, r.T, xs[i] * ps, 1.0, ps, p);
      EXPECT_LT(Rel(q, r.q[i]), 1e-9) << "T = " << r.T << ", x = " << xs[i];
    }
  }
}

// The layer-energy form with the printed Low et al. (2023, J. Chem. Eng. Data
// 68, Table 9) Lewatit VP OC 1065 values, the energies taken in J/mol: at
// 298.15 K C = 1.7395 and K = 0.8436 (worked independently).
TEST(GabIsotherm, LayerEnergyFormWithLowLewatitValues) {
  GabYoungIsothermParameters p;
  p.q_m = 3.23; p.C_E = 46830.0; p.D = 0.02446; p.F = 53310.0; p.G = -32.68;
  const GabConstants g = GabConstantsAt(298.15, p);
  EXPECT_NEAR(g.C, 1.73946278, 1e-8);
  EXPECT_NEAR(g.K, 0.8435891628, 1e-9);
  EXPECT_NEAR(GabLoading(g.q_m, g.C, g.K, 0.5), 3.12417468, 1e-7);
}

// The solver's per-cell Update gives exactly what EvalIsotherm gives, cell by
// cell, for both kinds; the water loading ignores the other components.
TEST(GabIsotherm, UpdateMatchesEvalIsotherm) {
  const size_t n = 6;
  ComponentEquation yi;
  ComponentIndependentEquation qi_sat;
  std::vector<double> T(n), P(n);
  for (const std::string c : {kWater, std::string("Nitrogen")}) {
    auto [it, ok] = yi.emplace(std::piecewise_construct, std::forward_as_tuple(c),
                               std::forward_as_tuple(c, CellType::CellCentre, EquationType::Differential,
                                                     [](UserData*) {}, 1e-8));
    it->second.var.resize(n);
    auto [jt, ok2] = qi_sat.emplace(std::piecewise_construct, std::forward_as_tuple(c),
                                    std::forward_as_tuple(c, CellType::CellCentre, [](UserData*) {}));
    jt->second.var.assign(n, -1.0);
  }
  for (size_t i = 0; i < n; ++i) {
    T[i] = 283.15 + 15.0 * i;
    P[i] = 1.0e5 + 2.0e4 * i;
    yi.at(kWater).var[i] = -1e-6 + 0.004 * i;
    yi.at("Nitrogen").var[i] = 1.0 - yi.at(kWater).var[i];
  }
  // Cell 0 has a negative mole fraction; cell 1 (298 K, 1.2 bar, y = 0.05) is
  // supersaturated, x = 1.9: the continuation above x_c.
  yi.at(kWater).var[1] = 0.05;
  GabIsotherm gab(kWater);
  gab.params_[kWater] = HandParams();
  GabYoungIsotherm young(kWater);
  for (IIsotherm* iso : std::vector<IIsotherm*>{&gab, &young}) {
    SCOPED_TRACE(IsothermKindString(iso->GetIsothermType()));
    IsothermUpdateParams params(yi, qi_sat, T, P, kWater, 0, static_cast<int>(n) - 1);
    iso->Update(params);
    for (size_t i = 0; i < n; ++i) {
      const double expected = EvalIsotherm(iso, kWater, T[i], P[i], yi.at(kWater).var[i]);
      EXPECT_DOUBLE_EQ(qi_sat.at(kWater).var[i], expected) << "cell " << i;
      EXPECT_TRUE(std::isfinite(expected)) << "cell " << i;
    }
    EXPECT_LT(qi_sat.at(kWater).var[0], 0.0);  // negative mole fraction: on the Henry line
  }
  EXPECT_DOUBLE_EQ(qi_sat.at("Nitrogen").var[0], -1.0);  // untouched
}

// The isosteric heat implied by the temperature dependence. Langmuir with
// b = b0 exp(Q / RT) gives Q at every loading (a check of IsostericHeat
// itself). For GAB at low loading, q = q_m C K p / p_sat, so
//   q_st = dH_vap + dH_m + dH_C + dH_K,  dH_vap = R T^2 d ln p_sat / dT
// (the Clausius-Clapeyron heat of p_sat).
TEST(GabIsotherm, IsostericHeat) {
  LangmuirIsotherm lang("CarbonDioxide");
  lang.params_["CarbonDioxide"].m1 = 3.0;
  lang.params_["CarbonDioxide"].b1 = 1e-10;
  lang.params_["CarbonDioxide"].q1 = 30000.0;
  for (double P : {1e3, 1e5, 1e7})
    EXPECT_NEAR(IsostericHeat(&lang, "CarbonDioxide", 300.0, P, 1.0), 30000.0, 0.5) << P;

  GabIsotherm gab(kWater);
  gab.params_[kWater] = HandParams();
  const double T = 298.15, R = constants::kIdealGasConstant;
  const double h = 0.01;
  const double dH_vap = R * T * T *
      (std::log(SaturationPressure(kWater, T + h)) - std::log(SaturationPressure(kWater, T - h))) / (2 * h);
  EXPECT_NEAR(dH_vap, 43990.0, 100.0);  // water at 25 degC: 43.99 kJ/mol
  const double ps = SaturationPressure(kWater, T);
  EXPECT_NEAR(IsostericHeat(&gab, kWater, T, 1e-6 * ps, 1.0), dH_vap + 500.0 + 6000.0 + 150.0, 20.0);
  // At higher loading the multilayer approaches condensation: q_st falls
  // towards dH_vap + dH_m + dH_K (the C term no longer acts).
  const double q_st_high = IsostericHeat(&gab, kWater, T, 0.9 * ps, 1.0);
  EXPECT_LT(q_st_high, dH_vap + 6650.0);
  EXPECT_GT(q_st_high, dH_vap);
  EXPECT_TRUE(std::isnan(IsostericHeat(&gab, kWater, T, 0.0, 1.0)));
}

// Plumbing: kind tags, field layouts, registry paths, builder and a deck round
// trip through the serialiser for both kinds.
TEST(GabIsotherm, KindsRegistryAndDeckRoundTrip) {
  IsothermType t{};
  ASSERT_TRUE(KindFromString("gab", t));
  EXPECT_EQ(t, IsothermType::GAB);
  ASSERT_TRUE(KindFromString("gab_young", t));
  EXPECT_EQ(t, IsothermType::GAB_YOUNG);
  EXPECT_TRUE(IsSupportedIsothermKind(IsothermType::GAB));
  EXPECT_FALSE(IsCompetitiveIsothermKind(IsothermType::GAB));
  EXPECT_FALSE(IsCompetitiveIsothermKind(IsothermType::GAB_YOUNG));
  IsothermType paired{};
  EXPECT_FALSE(PairedIsothermKind(IsothermType::GAB, paired));
  EXPECT_STREQ(isothermModelNames[static_cast<size_t>(IsothermType::GAB)], "GAB");

  Simulation sim;
  UserData& d = sim.Data();
  build::AddComponent(d, kWater);
  build::AddComponent(d, "Nitrogen");
  const AdsorbentID a = build::AddAdsorbent(d, "Amine A");
  build::SetIsotherm(d, a, "gab", {kWater});
  const AdsorbentID b = build::AddAdsorbent(d, "Amine B");
  build::SetIsotherm(d, b, "gab_young", {kWater});
  sim.RefreshParameters();

  const std::vector<std::pair<std::string, double>> gab_values = {
      {"q_m0", 3.51e0}, {"dH_m", 0.0}, {"C0", 2.73}, {"dH_C", 0.0}, {"K0", 0.867}, {"dH_K", 0.0}};
  const std::vector<std::pair<std::string, double>> young_values = {
      {"q_m", 3.23}, {"C_E", 46830.0}, {"D", 0.02446}, {"F", 53310.0}, {"G", -32.68},
      {"E10_A", 57220.0}, {"E10_B", -44.38}};
  for (const auto& [name, v] : gab_values) sim.SetParameter("adsorbent[Amine A].isotherm[Water]." + name, v);
  for (const auto& [name, v] : young_values) sim.SetParameter("adsorbent[Amine B].isotherm[Water]." + name, v);

  std::string text;
  SaveToString(d, text);
  UserData back;
  LoadDeckFlexible(back, text);

  auto* ga = dynamic_cast<GabIsotherm*>(back.adsorbent_directory.Get(a).adsorption_source_.GetIsotherm(kWater).get());
  ASSERT_NE(ga, nullptr);
  auto fields_a = IsothermFields(ga, kWater);
  ASSERT_EQ(fields_a.size(), gab_values.size());
  for (size_t i = 0; i < fields_a.size(); ++i) {
    EXPECT_EQ(fields_a[i].first, gab_values[i].first);
    EXPECT_DOUBLE_EQ(fields_a[i].second->Get(), gab_values[i].second) << gab_values[i].first;
  }
  auto* gb = dynamic_cast<GabYoungIsotherm*>(back.adsorbent_directory.Get(b).adsorption_source_.GetIsotherm(kWater).get());
  ASSERT_NE(gb, nullptr);
  auto fields_b = IsothermFields(gb, kWater);
  ASSERT_EQ(fields_b.size(), young_values.size());
  for (size_t i = 0; i < fields_b.size(); ++i) {
    EXPECT_EQ(fields_b[i].first, young_values[i].first);
    EXPECT_DOUBLE_EQ(fields_b[i].second->Get(), young_values[i].second) << young_values[i].first;
  }

  // The registry gives each path its unit.
  EXPECT_DOUBLE_EQ(sim.GetParameter("adsorbent[Amine A].isotherm[Water].K0"), 0.867);
  EXPECT_DOUBLE_EQ(sim.GetParameter("adsorbent[Amine B].isotherm[Water].G"), -32.68);

  // Grossmann et al. (2025, Table 5) values in the constant form: 1.76 mol/kg at RH 0.5.
  const double q = EvalIsotherm(ga, kWater, 298.15, 0.5 * SaturationPressure(kWater, 298.15), 1.0);
  const double expected = 3.51 * 2.73 * 0.867 * 0.5 / ((1 - 0.867 * 0.5) * (1 - 0.867 * 0.5 + 2.73 * 0.867 * 0.5));
  EXPECT_NEAR(q, expected, 1e-9);
}

// The fitter's evaluation path (what skarstrom.evaluate_isotherm and
// isotherm_parameters use): a GAB model takes water's saturation pressure.
TEST(GabIsotherm, FitterModelUsesWater) {
  isofit::Model m(IsothermType::GAB);
  ASSERT_TRUE(m.valid());
  EXPECT_EQ(m.Names(), (std::vector<std::string>{"q_m0", "dH_m", "C0", "dH_C", "K0", "dH_K"}));
  const double ps = SaturationPressure(kWater, 298.15);
  std::vector<isofit::Point> pts = {{0.2 * ps, 0.0, 298.15, 1.0}, {0.7 * ps, 0.0, 298.15, 1.0}};
  std::vector<double> q;
  std::string err;
  ASSERT_TRUE(isofit::Evaluate(IsothermType::GAB, {{"q_m0", 3.51}, {"dH_C", 0.0}, {"C0", 2.73}, {"K0", 0.867}},
                               pts, q, err)) << err;
  for (size_t i = 0; i < pts.size(); ++i) {
    const double x = pts[i].P / ps;
    const double expected = 3.51 * 2.73 * 0.867 * x / ((1 - 0.867 * x) * (1 - 0.867 * x + 2.73 * 0.867 * x));
    EXPECT_NEAR(q[i], expected, 1e-9) << i;
  }
  isofit::Model ethanol(IsothermType::GAB_YOUNG, "Ethanol");  // another vapour, by name
  EXPECT_TRUE(std::isfinite(ethanol.Eval(298.15, 1000.0, 1.0)));
}

}  // namespace adsim
