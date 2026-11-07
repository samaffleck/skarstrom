/*
 * Skarstrom - pressure and temperature swing adsorption simulator
 * Copyright (c) 2024-2026 Sam Affleck and contributors
 *
 * Licensed under the terms in the LICENSE file at the repository root.
 */
// The CO2 isotherm kinds whose loading depends on the adsorbed water loading:
// toth_wadst (Young et al. 2021, Energy Environ. Sci. 14, 5377, Eq. 14) and
// toth_humid (Stampi-Bombelli, Storione and Mazzotti 2020, Adsorption 26,
// 1183, Eqs. 11-14, corrected 2022).

#include <algorithm>
#include <cmath>
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
#include "adsim/objects/porous_media.h"
#include "adsim/simulation.h"
#include "adsim/solver/isotherm_fit.h"

namespace adsim {
namespace {

const std::string kCO2 = "CarbonDioxide";
const std::string kWater = kCoadsorbedWaterComponent;

double Rel(double a, double b) { return std::fabs(a / b - 1.0); }

struct Case { double p, T, q_w, q; };

double Wadst(const TothWadstIsothermParameters& ips, double p, double T, double q_w) {
  double q = 0.0;
  IsothermFunction::TothWadst(q, T, p, 1.0, q_w, ips);
  return q;
}

double Humid(const TothHumidIsothermParameters& ips, double p, double T, double q_w) {
  double q = 0.0;
  IsothermFunction::TothHumid(q, T, p, 1.0, q_w, ips);
  return q;
}

// Young et al. (2021)'s refit of the humid Toth to Lewatit VP OC 1065: the dry
// Toth of the WADST defaults written in the toth_2 form (b0 at T_0) with
// gamma = -0.137 and beta = 5.612 kg/mol.
TothHumidIsothermParameters YoungLewatitHumid() {
  TothHumidIsothermParameters p;
  p.q_sat = 4.86;
  p.x = 0.0;
  p.T_0 = 298.15;
  p.b0 = 2.85e-21 * std::exp(117798.0 / (constants::kIdealGasConstant * 298.15));
  p.dH0 = 117798.0;
  p.w0 = 0.209;
  p.alpha = 0.523;
  p.gamma = -0.137;
  p.beta = 5.612;
  return p;
}

// A component equation map with `n` cells for each name.
template <class Map, class... Args>
void AddCells(Map& map, const std::string& c, size_t n, double v, Args&&... args) {
  auto [it, ok] = map.emplace(std::piecewise_construct, std::forward_as_tuple(c),
                              std::forward_as_tuple(c, std::forward<Args>(args)...));
  it->second.var.assign(n, v);
}

}  // namespace

// Hand values (worked independently in Python, R = 8.31446261815324) with the
// defaults, which are the paper's Lewatit VP OC 1065 parameters:
// T_0 = 298.15 K; dry q_inf = 4.86 mol/kg, b0 = 2.85e-21 1/Pa,
// dH = 117798 J/mol, tau0 = 0.209, alpha = 0.523; wet q_inf = 9.035,
// b0 = 1.230e-18, dH = 203687, tau0 = alpha = 0.053; chi = 0; A = 1.532.
TEST(CoadsorptionIsotherms, WadstHandValues) {
  const TothWadstIsothermParameters ips;  // the defaults
  const std::vector<Case> cases = {
      {40, 298.15, 0.0, 0.8420245338722876},   {40, 298.15, 1.0, 1.0185996655819423},
      {40, 298.15, 3.0, 1.3323551030122467},   {1000, 313.15, 2.0, 1.719685426802942},
      {10, 288.15, 0.5, 0.7291403985743092},   {40, 373.15, 1.0, 0.1802687833911343},
      {1e4, 298.15, 5.0, 2.53738869009161}};
  for (const auto& c : cases)
    EXPECT_LT(Rel(Wadst(ips, c.p, c.T, c.q_w), c.q), 1e-12) << c.p << " Pa, " << c.T << " K, q_w " << c.q_w;
}

// At q_w -> 0 the weight exp(-A / q_w) vanishes: the kind is exactly the dry
// site, which is toth_2 with b0' = b0 exp(dH / (R T_0)). Negative water
// loadings (Newton iterates) are taken as 0.
TEST(CoadsorptionIsotherms, WadstDryLimitIsTheDryToth) {
  const TothWadstIsothermParameters ips;
  Toth2IsothermParameters dry;
  dry.q_sat = 4.86;
  dry.x = 0.0;
  dry.T_0 = 298.15;
  dry.b0 = 2.85e-21 * std::exp(117798.0 / (constants::kIdealGasConstant * 298.15));
  dry.dH0 = 117798.0;
  dry.w0 = 0.209;
  dry.alpha = 0.523;
  for (double T : {288.15, 298.15, 333.15, 373.15}) {
    for (double p : {1.0, 40.0, 1e3, 1e5}) {
      double q2 = 0.0;
      IsothermFunction::Toth2(q2, T, p, 1.0, dry);
      const double q0 = Wadst(ips, p, T, 0.0);
      EXPECT_LT(Rel(q0, q2), 1e-12) << T << " K, " << p << " Pa";
      EXPECT_EQ(Wadst(ips, p, T, -0.3), q0);
      EXPECT_EQ(Wadst(ips, p, T, 1e-4), q0);  // exp(-15320) == 0: no wet share yet
      EXPECT_LT(Rel(Wadst(ips, p, T, 0.02), q0), 1e-12);  // exp(-76.6): a negligible wet share
    }
  }
}

// Stampi-Bombelli et al. (2020), corrected 2022, APDES-NFC values (the
// defaults): T_0 = 296 K, n_s0 = 2.38 mol/kg, b0 = 70.7 1/kPa, t0 = 0.4148,
// dH0 = 57047 J/mol, alpha = -1.606, chi = 0, gamma = 0.0016 kg/mol,
// beta = 59.1 kg/mol. Then Young et al. (2021)'s Lewatit refit.
TEST(CoadsorptionIsotherms, HumidTothHandValues) {
  const TothHumidIsothermParameters ips;  // the defaults
  const std::vector<Case> apdes = {
      {40, 296, 0.0, 0.7119412653033451},     {40, 296, 0.5, 1.6748044208814166},
      {40, 296, 2.0, 1.9421735317365478},     {100, 313.15, 1.0, 1.3321524141099492},
      {10, 288.15, 0.2, 1.4077526705613985},  {40, 353.15, 1.0, 0.05229459454346384}};
  for (const auto& c : apdes)
    EXPECT_LT(Rel(Humid(ips, c.p, c.T, c.q_w), c.q), 1e-12) << c.p << " Pa, " << c.T << " K, q_w " << c.q_w;

  const auto young = YoungLewatitHumid();
  const std::vector<Case> lewatit = {
      {40, 298.15, 0.0, 0.8420245338722876}, {40, 298.15, 1.0, 1.2263823507502378},
      {40, 298.15, 3.0, 1.2196463858084956}, {1000, 323.15, 2.0, 1.681003440993059}};
  for (const auto& c : lewatit)
    EXPECT_LT(Rel(Humid(young, c.p, c.T, c.q_w), c.q), 1e-12) << c.p << " Pa, " << c.T << " K, q_w " << c.q_w;

  // b falls with temperature (dH0 > 0) at fixed water loading.
  EXPECT_GT(Humid(ips, 40, 288.15, 0.5), Humid(ips, 40, 313.15, 0.5));
}

// At q_w = 0 the humid Toth is toth_2 with the same parameters, bit for bit.
TEST(CoadsorptionIsotherms, HumidDryLimitIsToth2) {
  const TothHumidIsothermParameters ips;
  Toth2IsothermParameters t2;
  t2.q_sat = ips.q_sat.Get(); t2.x = ips.x.Get(); t2.T_0 = ips.T_0.Get(); t2.b0 = ips.b0.Get();
  t2.dH0 = ips.dH0.Get(); t2.w0 = ips.w0.Get(); t2.alpha = ips.alpha.Get();
  for (double T : {283.15, 296.0, 330.0})
    for (double p : {0.5, 40.0, 2e3}) {
      double q2 = 0.0;
      IsothermFunction::Toth2(q2, T, p, 1.0, t2);
      EXPECT_EQ(Humid(ips, p, T, 0.0), q2) << T << " " << p;
      EXPECT_EQ(Humid(ips, p, T, -1.0), q2);
    }
}

// Both kinds increase with the CO2 pressure at fixed (T, q_w).
TEST(CoadsorptionIsotherms, MonotonicInPressure) {
  const TothWadstIsothermParameters wadst;
  const TothHumidIsothermParameters humid;
  const auto young = YoungLewatitHumid();
  for (double T : {288.15, 298.15, 323.15, 363.15}) {
    for (double q_w : {0.0, 0.3, 1.0, 3.0, 6.0}) {
      double prev_w = 0.0, prev_h = 0.0, prev_y = 0.0;
      for (double p = 0.1; p <= 1e5; p *= 1.5) {
        const double qw = Wadst(wadst, p, T, q_w), qh = Humid(humid, p, T, q_w), qy = Humid(young, p, T, q_w);
        EXPECT_GT(qw, prev_w) << T << " K, q_w " << q_w << ", " << p << " Pa";
        EXPECT_GT(qh, prev_h) << T << " K, q_w " << q_w << ", " << p << " Pa";
        EXPECT_GT(qy, prev_y) << T << " K, q_w " << q_w << ", " << p << " Pa";
        prev_w = qw;
        prev_h = qh;
        prev_y = qy;
      }
    }
  }
}

// pyDAC (MIT licence, PyPI pydac3 3.1.0, src/pydac/core/isotherm_models.py,
// TothIsotherm.co_isotherm_wadst with its "Cleaned-Toth" parameters) writes the
// same WADST with each site's affinity at a reference temperature,
// b = b0 exp(DH0 / (R T0) (1 - T0 / T)) with DH0 < 0, and constant capacity;
// in this kind's parameters b0 -> b0 exp(DH0 / (R T0)), dH = -DH0, chi = 0.
// Reference values from a transcription of that function (scipy's R).
TEST(CoadsorptionIsotherms, WadstCrossCheckAgainstPyDac) {
  const double R = constants::kIdealGasConstant, T0 = 288.0;
  TothWadstIsothermParameters p;
  p.T_0 = T0;
  p.q_inf_dry = 4.077685;
  p.chi_dry = 0.0;
  p.b0_dry = 3.07587735 * std::exp(-104935.08 / (R * T0));
  p.dH_dry = 104935.08;
  p.tau0_dry = 0.222508;
  p.alpha_dry = 0.548347;
  p.q_inf_wet = 4.586366;
  p.chi_wet = 0.0;
  p.b0_wet = 6500.42 * std::exp(-130435.87 / (R * T0));
  p.dH_wet = 130435.87;
  p.tau0_wet = 0.148360;
  p.alpha_wet = 0.202008;
  p.A = 0.221114;
  const std::vector<Case> cases = {
      {40, 288.15, 0.01, 1.0816674098478305}, {40, 298.15, 0.5, 1.2613825059653363},
      {40, 298.15, 2.0, 1.412848760290231},   {1000, 313.15, 1.0, 1.8306437747608122},
      {10, 353.15, 3.0, 0.16796428493845794}, {400, 333.15, 0.2, 0.8362360206964081}};
  // The gas constants differ by 2e-11 relative (scipy's 8.314462618).
  for (const auto& c : cases) EXPECT_LT(Rel(Wadst(p, c.p, c.T, c.q_w), c.q), 1e-7) << c.p << " " << c.T;
}

// Update reads the water loading from the adsorbent's loadings (q_map) cell by
// cell and matches EvalIsotherm; without a loading map, or without a water
// entry in it, the kind is its dry isotherm.
TEST(CoadsorptionIsotherms, UpdateReadsTheAdsorbedWaterLoading) {
  const size_t n = 5;
  ComponentEquation yi, q_map;
  ComponentIndependentEquation qi_sat;
  std::vector<double> T(n), P(n);
  for (const std::string& c : {kCO2, kWater}) {
    AddCells(yi, c, n, 0.0, CellType::CellCentre, EquationType::Differential, [](UserData*) {}, 1e-8);
    AddCells(q_map, c, n, 0.0, CellType::CellCentre, EquationType::Differential, [](UserData*) {}, 1e-8);
    AddCells(qi_sat, c, n, -1.0, CellType::CellCentre, [](UserData*) {});
  }
  for (size_t i = 0; i < n; ++i) {
    T[i] = 288.15 + 10.0 * i;
    P[i] = 1.0e5;
    yi.at(kCO2).var[i] = 4e-4 * (i + 1);
    yi.at(kWater).var[i] = 0.01;
    q_map.at(kWater).var[i] = 0.8 * i - 0.4;  // cell 0 negative: dry
    q_map.at(kCO2).var[i] = 99.0;             // not read
  }
  TothWadstIsotherm wadst(kCO2);
  TothHumidIsotherm humid(kCO2);
  for (IIsotherm* iso : std::vector<IIsotherm*>{&wadst, &humid}) {
    SCOPED_TRACE(IsothermKindString(iso->GetIsothermType()));
    IsothermUpdateParams with(yi, qi_sat, T, P, kCO2, 0, static_cast<int>(n) - 1, nullptr, &q_map);
    iso->Update(with);
    for (size_t i = 0; i < n; ++i) {
      const double expected = EvalIsotherm(iso, kCO2, T[i], P[i], yi.at(kCO2).var[i], q_map.at(kWater).var[i]);
      EXPECT_DOUBLE_EQ(qi_sat.at(kCO2).var[i], expected) << "cell " << i;
      if (i > 1) EXPECT_NE(expected, EvalIsotherm(iso, kCO2, T[i], P[i], yi.at(kCO2).var[i])) << "cell " << i;
    }
    EXPECT_DOUBLE_EQ(qi_sat.at(kCO2).var[0], EvalIsotherm(iso, kCO2, T[0], P[0], yi.at(kCO2).var[0]));
    EXPECT_DOUBLE_EQ(qi_sat.at(kWater).var[0], -1.0);  // untouched

    IsothermUpdateParams without(yi, qi_sat, T, P, kCO2, 0, static_cast<int>(n) - 1);
    iso->Update(without);
    for (size_t i = 0; i < n; ++i)
      EXPECT_DOUBLE_EQ(qi_sat.at(kCO2).var[i], EvalIsotherm(iso, kCO2, T[i], P[i], yi.at(kCO2).var[i]));
  }
}

// Plumbing: kind tags, field layouts, registry paths with units, builder and a
// deck round trip through the serialiser for both kinds.
TEST(CoadsorptionIsotherms, KindsRegistryAndDeckRoundTrip) {
  for (const auto& [tag, type] : {std::pair<const char*, IsothermType>{"toth_wadst", IsothermType::TOTH_WADST},
                                  {"toth_humid", IsothermType::TOTH_HUMID}}) {
    IsothermType t{};
    ASSERT_TRUE(KindFromString(tag, t));
    EXPECT_EQ(t, type);
    EXPECT_STREQ(IsothermKindString(type), tag);
    EXPECT_TRUE(IsSupportedIsothermKind(type));
    EXPECT_FALSE(IsCompetitiveIsothermKind(type));
    EXPECT_TRUE(IsWaterCoupledIsothermKind(type));
    IsothermType paired{};
    EXPECT_FALSE(PairedIsothermKind(type, paired));
  }
  EXPECT_FALSE(IsWaterCoupledIsothermKind(IsothermType::TOTH_2));
  EXPECT_FALSE(IsWaterCoupledIsothermKind(IsothermType::GAB));
  // Appended after the GAB kinds: existing enum values are unchanged.
  EXPECT_EQ(static_cast<int>(IsothermType::GAB_YOUNG), 22);
  EXPECT_EQ(static_cast<int>(IsothermType::TOTH_WADST), 23);
  EXPECT_EQ(static_cast<int>(IsothermType::TOTH_HUMID), 24);

  Simulation sim;
  UserData& d = sim.Data();
  build::AddComponent(d, kCO2);
  build::AddComponent(d, kWater);
  build::AddComponent(d, "Nitrogen");
  const AdsorbentID a = build::AddAdsorbent(d, "Amine A");
  build::SetIsotherm(d, a, "toth_wadst", {kCO2});
  build::SetIsotherm(d, a, "gab_young", {kWater});
  const AdsorbentID b = build::AddAdsorbent(d, "Amine B");
  build::SetIsotherm(d, b, "toth_humid", {kCO2});
  sim.RefreshParameters();

  const std::vector<std::pair<std::string, double>> wadst_values = {
      {"T_0", 300.0},      {"q_inf_dry", 4.5},  {"chi_dry", 0.1},   {"b0_dry", 3e-21}, {"dH_dry", 1.1e5},
      {"tau0_dry", 0.25},  {"alpha_dry", 0.4},  {"q_inf_wet", 8.0}, {"chi_wet", 0.2},  {"b0_wet", 2e-18},
      {"dH_wet", 2e5},     {"tau0_wet", 0.06},  {"alpha_wet", 0.07}, {"A", 1.4}};
  const std::vector<std::pair<std::string, double>> humid_values = {
      {"q_sat", 2.4}, {"x", 0.3},  {"T_0", 295.0},  {"b0", 0.07},  {"dH0", 57000.0},
      {"w0", 0.41},   {"alpha", -1.5}, {"gamma", 0.002}, {"beta", 60.0}};
  for (const auto& [name, v] : wadst_values)
    sim.SetParameter("adsorbent[Amine A].isotherm[CarbonDioxide]." + name, v);
  for (const auto& [name, v] : humid_values)
    sim.SetParameter("adsorbent[Amine B].isotherm[CarbonDioxide]." + name, v);

  std::string text;
  SaveToString(d, text);
  UserData back;
  LoadDeckFlexible(back, text);

  auto check = [&](AdsorbentID id, IsothermType type, const std::vector<std::pair<std::string, double>>& values) {
    IIsotherm* iso = back.adsorbent_directory.Get(id).adsorption_source_.GetIsotherm(kCO2).get();
    ASSERT_NE(iso, nullptr);
    EXPECT_EQ(iso->GetIsothermType(), type);
    const auto fields = IsothermFields(iso, kCO2);
    ASSERT_EQ(fields.size(), values.size());
    for (size_t i = 0; i < fields.size(); ++i) {
      EXPECT_EQ(fields[i].first, values[i].first);
      EXPECT_DOUBLE_EQ(fields[i].second->Get(), values[i].second) << values[i].first;
    }
  };
  check(a, IsothermType::TOTH_WADST, wadst_values);
  check(b, IsothermType::TOTH_HUMID, humid_values);

  auto unit_of = [&](const std::string& path) {
    const Parameter* p = sim.Parameters().Find(path);
    return p ? p->Info().unit : std::string("<missing>");
  };
  const std::string pa = "adsorbent[" + std::to_string(a) + "].isotherm[CarbonDioxide].";
  const std::string pb = "adsorbent[" + std::to_string(b) + "].isotherm[CarbonDioxide].";
  EXPECT_EQ(unit_of(pa + "b0_dry"), "1/Pa");
  EXPECT_EQ(unit_of(pa + "A"), "mol/kg");
  EXPECT_EQ(unit_of(pb + "beta"), "kg/mol");
  EXPECT_DOUBLE_EQ(sim.GetParameter("adsorbent[Amine B].isotherm[CarbonDioxide].gamma"), 0.002);
}

// The deck checks: IAST cannot mix these kinds (an error, not a silently inert
// component), and without a water component they are only their dry form (a
// warning). The IAST code itself also refuses them.
TEST(CoadsorptionIsotherms, DeckChecks) {
  Simulation sim;
  UserData& d = sim.Data();
  build::AddComponent(d, kCO2);
  build::AddComponent(d, "Nitrogen");
  const AdsorbentID a = build::AddAdsorbent(d, "Amine");
  build::SetIsotherm(d, a, "toth_wadst", {kCO2});
  MixedAdsorbentLayer layer;
  layer.name_ = "L";
  layer.AddAdsorbent(a, 1.0);

  auto issues = [&] { return CheckLayerComposition(layer, d.adsorbent_directory); };
  auto has = [](const std::vector<LayerCompositionIssue>& v, LayerCompositionIssue::Severity s, const char* text) {
    for (const auto& i : v)
      if (i.severity == s && i.message.find(text) != std::string::npos) return true;
    return false;
  };
  using S = LayerCompositionIssue::Severity;
  auto v = issues();
  EXPECT_TRUE(has(v, S::Warning, "has no such component")) << (v.empty() ? "" : v[0].message);
  EXPECT_FALSE(has(v, S::Error, "IAST"));

  build::AddComponent(d, kWater);
  build::SetIsotherm(d, a, "gab", {kWater});
  v = issues();
  EXPECT_FALSE(has(v, S::Warning, "has no such component"));
  EXPECT_TRUE(v.empty()) << v[0].message;

  build::SetMixingRule(d, a, "iast");
  v = issues();
  EXPECT_TRUE(has(v, S::Error, "cannot mix")) << (v.empty() ? "" : v[0].message);

  TothWadstIsotherm iso(kCO2);
  EXPECT_THROW(ResolvePureIsotherm(&iso, kCO2), std::invalid_argument);
}

// In a run: humid air through a Lewatit-like bed (toth_wadst for CO2, the
// Young GAB for water, nitrogen inert). At the start every cell's CO2 loading
// is in equilibrium with the water loading already on the sorbent (the second
// equilibrium pass), and a short feed step runs; the solver passes the
// adsorbent's loadings to the isotherm.
TEST(CoadsorptionIsotherms, HumidBreakthroughStartsInEquilibriumWithTheWater) {
  Simulation sim;
  UserData& d = sim.Data();
  build::AddComponent(d, kCO2);
  build::AddComponent(d, kWater);
  build::AddComponent(d, "Nitrogen");
  const AdsorbentID ads = build::AddAdsorbent(d, "Amine");
  build::SetIsotherm(d, ads, "toth_wadst", {kCO2});
  build::SetIsotherm(d, ads, "gab_young", {kWater});
  for (const auto& c : {kCO2, kWater, std::string("Nitrogen")}) build::SetRateModel(d, ads, c, "ldf", "constant");
  build::SetBedShape(d, "cylinder");
  build::SetWallMaterial(d, "ss316l");
  build::SetWallBoundary(d, "zero_flux");
  build::SetThermalEquilibrium(d, "lte");
  const size_t layer = build::AddLayer(d, "Amine");
  build::SetLayerComposition(d, layer, {{ads, 1.0}});
  build::SetMomentum(d, layer, "ergun");
  build::SetDispersion(d, layer, "constant");
  build::SetHeatDispersion(d, layer, "constant");
  build::SetInternalHeatTransfer(d, layer, "constant");
  build::SetExternalHeatTransfer(d, layer, "constant");
  const size_t feed = build::AddStep(d, "Feed");
  build::SetFlowBoundary(d, feed, "bottom", "inlet");
  build::AddStream(d, feed, "bottom", "Feed inlet", "velocity", "constant");
  build::SetPressureBoundary(d, feed, "top", "outlet", "constant");
  sim.RefreshParameters();
  const std::vector<std::pair<std::string, std::string>> values = {
      {"bed.diameter", "0.032"},
      {"wall.thickness", "0.0031"},
      {"adsorbent[Amine].solid_density", "1100"},
      {"adsorbent[Amine].heat_capacity", "1500"},
      {"adsorbent[Amine].thermal_conductivity", "0.2"},
      {"adsorbent[Amine].particle_void_fraction", "0.3"},
      {"adsorbent[Amine].kinetics[CarbonDioxide].k", "0.01"},
      {"adsorbent[Amine].kinetics[Water].k", "0.002"},
      {"adsorbent[Amine].kinetics[Nitrogen].k", "0.1"},
      {"adsorbent[Amine].heat_of_adsorption[CarbonDioxide]", "70000"},
      {"adsorbent[Amine].heat_of_adsorption[Water]", "46000"},
      {"layer[Amine].length", "0.05"},
      {"layer[Amine].cells", "10"},
      {"layer[Amine].bed_void_fraction", "0.4"},
      {"layer[Amine].particle_diameter", "0.0005"},
      {"layer[Amine].dispersion_coefficient", "1e-4"},
      {"initial.pressure", "1e5"},
      {"initial.temperature", "298.15"},
      {"initial.y[CarbonDioxide]", "4e-4"},
      {"initial.y[Water]", "0.01"},
      {"initial.y[Nitrogen]", "0.9896"},
      {"step[Feed].time", "30"},
      {"step[Feed].bottom[Feed inlet].velocity", "0.05"},
      {"step[Feed].bottom[Feed inlet].temperature", "298.15"},
      {"step[Feed].bottom[Feed inlet].y[CarbonDioxide]", "4e-4"},
      {"step[Feed].bottom[Feed inlet].y[Water]", "0.02"},
      {"step[Feed].bottom[Feed inlet].y[Nitrogen]", "0.9796"},
      {"step[Feed].top.pressure", "1e5"},
  };
  for (const auto& [path, v] : values) sim.SetParameterText(path, v);
  std::string error;
  ASSERT_TRUE(sim.Check(error)) << error;

  // Initial equilibrium, cell by cell (Run initialises again).
  {
    UserData& copy = d;
    copy.Initialise();
    IIsotherm* iso = copy.adsorbent_directory.Get(ads).adsorption_source_.GetIsotherm(kCO2).get();
    const auto& q = copy.state.ads_concentration.at(ads);
    const auto& s = copy.state;
    const int first = copy.bed.GetFirstBedCellIndex(), last = copy.bed.GetLastBedCellIndex();
    ASSERT_LT(first, last);
    for (int i = first; i <= last; ++i) {
      const double q_w = q.at(kWater).var[i];
      ASSERT_GT(q_w, 0.5) << i;  // GAB at 31 % relative humidity
      const double T = s.solid_temperature.var[i], P = s.pressure.var[i];
      const double y = s.mole_fraction.at(kCO2).var[i];
      const double expected = EvalIsotherm(iso, kCO2, T, P, y, q_w);
      EXPECT_NEAR(q.at(kCO2).var[i], expected, 1e-12 * expected) << i;
      EXPECT_GT(expected, EvalIsotherm(iso, kCO2, T, P, y) * 1.01) << i;  // above the dry loading
    }
  }

  RunOptions options;
  options.cycles = 1;
  const RunResult r = sim.Run(options);
  ASSERT_TRUE(r.Ok()) << r.message;
  const auto qw_bottom = sim.Series("adsorbed_concentration", kWater, SeriesLocation::Bottom).values;
  ASSERT_FALSE(qw_bottom.empty());
  EXPECT_GT(qw_bottom.back(), qw_bottom.front());  // the wetter feed loads water at the inlet
  for (double v : sim.Series("adsorbed_concentration", kCO2, SeriesLocation::Bottom).values)
    EXPECT_TRUE(std::isfinite(v));
}

// Fitter: dry data fit the dry block and hold the water terms; humid data at
// fixed q_w fit the water terms with the dry block held. Synthetic data from
// the defaults are recovered.
TEST(CoadsorptionIsotherms, FitterSetUpAndRecovery) {
  for (auto type : {IsothermType::TOTH_WADST, IsothermType::TOTH_HUMID}) {
    SCOPED_TRACE(IsothermKindString(type));
    const auto& kinds = isofit::FittableKinds();
    EXPECT_NE(std::find(kinds.begin(), kinds.end(), type), kinds.end());
    EXPECT_TRUE(isofit::IsTemperatureDependent(type));

    isofit::Model truth(type);
    std::vector<isofit::Point> dry, wet;
    for (double T : {288.15, 298.15, 313.15})
      for (double p : {5.0, 20.0, 40.0, 100.0, 400.0, 1000.0}) {
        dry.push_back({p, truth.Eval(T, p, 1.0, 0.0), T, 1.0, 0.0});
        for (double q_w : {0.5, 1.5, 3.0}) wet.push_back({p, truth.Eval(T, p, 1.0, q_w), T, 1.0, q_w});
      }

    std::vector<isofit::ParamSpec> specs;
    double t_ref = 0.0;
    std::string err;
    ASSERT_TRUE(isofit::DefaultSpecs(type, dry, specs, t_ref, err)) << err;
    const std::vector<std::string> water_terms = type == IsothermType::TOTH_WADST
        ? std::vector<std::string>{"q_inf_wet", "chi_wet", "b0_wet", "dH_wet", "tau0_wet", "alpha_wet", "A"}
        : std::vector<std::string>{"gamma", "beta"};
    auto spec = [&](const std::string& name) -> const isofit::ParamSpec& {
      for (const auto& s : specs) if (s.name == name) return s;
      ADD_FAILURE() << "no spec " << name;
      return specs.front();
    };
    for (const auto& name : water_terms) EXPECT_TRUE(spec(name).fixed) << name;
    EXPECT_EQ(specs.size(), truth.Names().size());

    // Humid data only: the dry block is held at the given values (here the
    // truth, as the defaults), the water terms are fitted.
    ASSERT_TRUE(isofit::DefaultSpecs(type, wet, specs, t_ref, err)) << err;
    for (const auto& name : water_terms) EXPECT_FALSE(spec(name).fixed) << name;
    const std::string dry_capacity = type == IsothermType::TOTH_WADST ? "q_inf_dry" : "q_sat";
    EXPECT_TRUE(spec(dry_capacity).fixed);
    EXPECT_TRUE(spec("T_0").fixed);
    EXPECT_DOUBLE_EQ(spec("T_0").value, truth.Get("T_0"));
    // Keep the test short: fit only the main water terms, the rest held at truth.
    for (auto& s : specs) {
      const bool keep = type == IsothermType::TOTH_WADST ? (s.name == "A" || s.name == "q_inf_wet")
                                                         : (s.name == "gamma" || s.name == "beta");
      if (!keep) { s.fixed = true; s.value = truth.Get(s.name); s.pair_energy.clear(); }
    }
    isofit::FitOptions opt;
    opt.max_global_evals = 60;
    const auto r = isofit::Fit(type, wet, specs, t_ref, opt);
    ASSERT_TRUE(r.ok) << r.error;
    EXPECT_LT(r.rmse, 1e-6);
    for (const auto& pr : r.params)
      if (!pr.fixed) EXPECT_LT(Rel(pr.value, truth.Get(pr.name)), 1e-4) << pr.name;
  }
}

}  // namespace adsim
