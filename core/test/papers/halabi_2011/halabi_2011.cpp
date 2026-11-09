/*
 * Skarstrom - pressure and temperature swing adsorption simulator
 * Copyright (c) 2024-2026 Sam Affleck and contributors
 *
 * Licensed under the terms in the LICENSE file at the repository root.
 */
// Halabi (2011), sorption-enhanced steam methane reforming over Rh/Ce0.6Zr0.4O2
// mixed with K2CO3-promoted hydrotalcite (PhD thesis, TU Eindhoven,
// doi:10.6100/IR709035). See README.md in this folder for the provenance of
// every value, the model-form audit and the attribution runs.
//
//   GasPhaseEquilibrium      the chapter-8 rate model on an inert bed: outlet
//                            CH4 conversion against the equilibria the thesis
//                            prints (ch. 9, Fig. 7 legend and Table 6).
//   SorptionEnhanced_1p5bar  Table 6 at 1.5 bar, S/C 3-6: pre-breakthrough CH4
//                            conversion and dry, Ar-free H2 purity.
//   SorptionEnhanced_HighPressure_RecordOnly
//                            the 4.5 and 6 bar rows of Table 6, recorded only
//                            (README: the Freundlich fit underestimates the
//                            sorbent capacity at these conditions).
//
// The decks are built from kind names (deck_builder.h) and the registry.

#include <chrono>
#include <cmath>
#include <cstdio>
#include <iostream>
#include <map>
#include <string>
#include <utility>
#include <vector>

#include <gtest/gtest.h>

#include "adsim/math/constants.h"
#include "adsim/objects/deck_builder.h"
#include "adsim/simulation.h"

#include "../../fixtures/validation_metrics.h"

namespace adsim {
namespace {

const std::string CH4 = "Methane", H2O = "Water", H2 = "Hydrogen", CO = "CarbonMonoxide", CO2 = "CarbonDioxide",
                  AR = "Argon";
constexpr const char* kPaper = "halabi_2011";

constexpr double kR = constants::kIdealGasConstant;  // [J/(mol K)]
constexpr double kT = 723.15;             // 450 C, every case [K]

// ---- thesis constants (ch. 8, Tables 2-4), converted from kPa to bar ----------
//
// The engine's partial-pressure basis is bar. A constant multiplying a product
// of partial pressures of total degree s converts as C_bar = C_kPa * 100^s.

// Reaction equilibrium constants K = A exp(-theta / T) (Table 2).
constexpr double kKI_A = 1.198e17 / 1e4, kKI_theta = 26830.0;    // (kPa)^2 -> bar^2
constexpr double kKII_A = 1.767e-2, kKII_theta = -4400.0;        // dimensionless
constexpr double kKIII_A = 2.117e15 / 1e4, kKIII_theta = 22430.0;  // (kPa)^2 -> bar^2

double KI(double T) { return kKI_A * std::exp(-kKI_theta / T); }
double KII(double T) { return kKII_A * std::exp(-kKII_theta / T); }

// One factor of a site balance: K0 [bar^-s], van't Hoff enthalpy dH [J/mol],
// exponents of the partial pressures.
struct SiteTerm {
  double K0;
  double dH;
  std::map<std::string, double> exponent;
};

// Omega_l = 1 / (1 + K_CH4 p_CH4 / p_H2^0.5 + K_CO p_CO + K_CO2 p_CO2 + K_H2 p_H2),
// Omega_S = 1 / (1 + K_H2O p_H2O / p_H2 + K_H2 p_H2)   (Table 1, Eq. 4; Table 4).
const std::vector<SiteTerm> kSiteL = {
    {1.49e-8 * 10.0, -98.8e3, {{CH4, 1.0}, {H2, -0.5}}},   // degree 0.5
    {2.34e-8 * 100.0, -111.2e3, {{CO, 1.0}}},              // degree 1
    {8.33e-10 * 100.0, -115.6e3, {{CO2, 1.0}}},            // degree 1
    {6.11e-6 * 100.0, -60.8e3, {{H2, 1.0}}},               // degree 1
};
const std::vector<SiteTerm> kSiteS = {
    {3.14e8, 126.9e3, {{H2O, 1.0}, {H2, -1.0}}},  // degree 0
    {6.11e-6 * 100.0, -60.8e3, {{H2, 1.0}}},      // degree 1
};

// Omega_l Omega_S = 1 / [(1 + sum_l)(1 + sum_S)] = 1 / (1 + sum_l + sum_S + sum_l x sum_S):
// one LHHW denominator of 4 + 2 + 8 = 14 terms (exact).
std::vector<SiteTerm> ExpandedDenominator() {
  std::vector<SiteTerm> out(kSiteL);
  out.insert(out.end(), kSiteS.begin(), kSiteS.end());
  for (const auto& l : kSiteL)
    for (const auto& s : kSiteS) {
      SiteTerm t{l.K0 * s.K0, l.dH + s.dH, l.exponent};
      for (const auto& [c, e] : s.exponent) t.exponent[c] += e;
      out.push_back(t);
    }
  return out;
}

// Rate laws (Table 1, Eqs. 1-3; Table 3; Table 2 for K).
struct RateLaw {
  std::string name;
  std::map<std::string, double> nu;
  double k0;                             // [mol/(kg_cat s) bar^-degree]
  double E;                              // [J/mol]
  std::map<std::string, double> orders;  // overrides of the default orders
  double KA, Ktheta;                     // equilibrium constant, engine units
};

const std::vector<RateLaw> kRates = {
    // R1 = k1 / p_H2^2.5 (p_CH4 p_H2O - p_H2^3 p_CO / K_I) Omega_l Omega_S; degree -0.5
    {"R1", {{CH4, -1}, {H2O, -1}, {CO, 1}, {H2, 3}}, 1.62e8 / 10.0, 83.8e3, {{H2, -2.5}}, kKI_A, kKI_theta},
    // R2 = k2 / p_H2 (p_CO p_H2O - p_H2 p_CO2 / K_II) Omega_l Omega_S; degree 1
    {"R2", {{CO, -1}, {H2O, -1}, {CO2, 1}, {H2, 1}}, 2.34e5 * 100.0, 15.1e3, {{H2, -1.0}}, kKII_A, kKII_theta},
    // R3 = k3 / p_H2^3.5 (p_CH4 p_H2O^2 - p_H2^4 p_CO2 / K_III) Omega_l Omega_S; degree -0.5
    {"R3", {{CH4, -1}, {H2O, -2}, {CO2, 1}, {H2, 4}}, 4.55e8 / 10.0, 89.2e3, {{H2, -3.5}}, kKIII_A, kKIII_theta},
};

// ---- feeds -----------------------------------------------------------------------

// Mole fractions of a CH4 / steam / H2 / Ar feed.
std::map<std::string, double> Feed(double y_ch4, double steam_to_carbon, double y_h2) {
  const double y_h2o = steam_to_carbon * y_ch4;
  return {{CH4, y_ch4}, {H2O, y_h2o}, {H2, y_h2}, {CO, 0.0}, {CO2, 0.0}, {AR, 1.0 - y_ch4 - y_h2o - y_h2}};
}

// Ch. 9, Table 6 feed: 3.94 vol% CH4, 0.1 vol% H2, Ar balance (p. 203).
std::map<std::string, double> Table6Feed(double sc) { return Feed(0.0394, sc, 0.001); }
// Ch. 9, Fig. 7 feed: 3 vol% CH4, 0.1 vol% H2, S/C 4 (p. 205 and the caption).
std::map<std::string, double> Fig7Feed() { return Feed(0.03, 4.0, 0.001); }

// 75 nml/min (0 C, 1 atm; both feeds, pp. 203 and 205) as a molar flow [mol/s].
double FeedMolarFlow() { return 101325.0 * 75e-6 / 60.0 / (kR * 273.15); }

// ---- gas-phase equilibrium (R1 and R2; R3 is their sum) ------------------------
//
// Moles per mole of feed after extents a (R1) and b (R2), as mole fractions.
std::map<std::string, double> Advance(const std::map<std::string, double>& y0, double a, double b) {
  std::map<std::string, double> n = y0;
  n[CH4] -= a;
  n[H2O] -= a + b;
  n[CO] += a - b;
  n[H2] += 3.0 * a + b;
  n[CO2] += b;
  double total = 0.0;
  for (const auto& [c, v] : n) total += v;
  for (auto& [c, v] : n) v /= total;
  return n;
}

// Bisection on a monotone function f with f(lo) < 0 < f(hi).
template <class F>
double Bisect(F f, double lo, double hi) {
  for (int i = 0; i < 200; ++i) {
    const double m = 0.5 * (lo + hi);
    (f(m) > 0.0 ? hi : lo) = m;
  }
  return 0.5 * (lo + hi);
}

// The equilibrium of `y0` at (T, P [bar]) under the thesis's K_I and K_II.
std::map<std::string, double> Equilibrium(const std::map<std::string, double>& y0, double T, double P) {
  const double ch4 = y0.at(CH4);
  auto wgs_extent = [&](double a) {
    auto g = [&](double b) {
      auto y = Advance(y0, a, b);
      return std::log(y[CO2] * y[H2] / (y[CO] * y[H2O])) - std::log(KII(T));
    };
    return Bisect(g, 1e-300, a * (1.0 - 1e-15));
  };
  auto h = [&](double a) {
    auto y = Advance(y0, a, wgs_extent(a));
    return std::log(y[CO] * std::pow(y[H2] * P, 3) * P / (y[CH4] * P * y[H2O] * P)) - std::log(KI(T));
  };
  const double a = Bisect(h, ch4 * 1e-12, ch4 * (1.0 - 1e-12));
  return Advance(y0, a, wgs_extent(a));
}

// ---- deck ------------------------------------------------------------------------

std::string Num(double v) {
  char buf[40];
  std::snprintf(buf, sizeof(buf), "%.17g", v);
  return buf;
}

struct Case {
  std::map<std::string, double> feed;
  double pressure_bar = 1.5;
  bool sorbent_active = true;  // false: the sorbent is inert (gas-phase equilibrium case)
  double step_time = 60.0;     // [s]; the read time (README, "Reading")
  int cells = 40;
};

struct Outcome {
  bool ok = false;
  std::string message;
  bool healthy = true;
  std::map<std::string, double> y_out;  // outlet mole fractions at the end of the step
  double wall_s = 0.0;
};

// Bed (ch. 9, pp. 192-195, 201): 9.0 mm ID quartz tube, 4.2 cm bed of 2.0 g
// K-HTC and 0.120 g Rh/Ce0.6Zr0.4O2, 0.3 mm pellets, 450 C.
constexpr double kDiameter = 9.0e-3;          // [m]
constexpr double kLength = 0.042;             // [m]
constexpr double kMassSorbent = 2.0e-3;       // [kg]
constexpr double kMassCatalyst = 0.120e-3;    // [kg]
constexpr double kPelletDensity = 2000.0;     // [kg/m3] ch. 8, Table 12 (assumed for both solids)
constexpr double kParticleDiameter = 0.3e-3;  // [m]

double BedVoidFraction() {
  const double volume = M_PI / 4.0 * kDiameter * kDiameter * kLength;
  return 1.0 - (kMassSorbent + kMassCatalyst) / volume / kPelletDensity;  // 0.603
}

Outcome RunCase(const Case& c) {
  Simulation sim;
  UserData& d = sim.Data();
  for (const auto& comp : {CH4, H2O, H2, CO, CO2, AR}) build::AddComponent(d, comp);
  const AdsorbentID htc = build::AddAdsorbent(d, "K-HTC");
  const AdsorbentID cat = build::AddAdsorbent(d, "Rh catalyst");
  if (c.sorbent_active) {
    build::SetIsotherm(d, htc, "freundlich", {CO2});
    build::SetRateModel(d, htc, CO2, "ldf", "series_resistance");
  }
  build::SetBedShape(d, "cylinder");
  build::SetWallMaterial(d, "ss316l");
  build::SetWallBoundary(d, "zero_flux");
  const size_t layer = build::AddLayer(d, "Bed");
  const double w_cat = kMassCatalyst / (kMassSorbent + kMassCatalyst);
  build::SetLayerComposition(d, layer, {{htc, 1.0 - w_cat}, {cat, w_cat}});
  build::SetMomentum(d, layer, "ergun");
  build::SetDispersion(d, layer, "ruthven");

  const size_t feed = build::AddStep(d, "Feed");
  build::SetFlowBoundary(d, feed, "bottom", "inlet");
  build::AddStream(d, feed, "bottom", "Feed inlet", "molar_flow", "constant");
  build::SetPressureBoundary(d, feed, "top", "outlet", "constant");

  const auto terms = ExpandedDenominator();
  for (const auto& r : kRates) {
    const ReactionID id = build::AddReaction(d, cat, r.name, r.nu, "lhhw");
    for (size_t j = 0; j < terms.size(); ++j) build::AddAdsorptionTerm(d, cat, id);
  }
  sim.RefreshParameters();

  const double P = c.pressure_bar * 1e5;
  std::vector<std::pair<std::string, std::string>> v = {
      {"bed.diameter", Num(kDiameter)},
      {"adsorbent[K-HTC].solid_density", Num(kPelletDensity)},
      // Table 2: pore volume 0.037 cm3/g x 2000 kg/m3.
      {"adsorbent[K-HTC].particle_void_fraction", Num(0.037e-3 * kPelletDensity)},
      {"adsorbent[Rh catalyst].solid_density", Num(kPelletDensity)},
      // Ch. 9 Table 1: pore volume 0.20 cm3/g x 2000 kg/m3.
      {"adsorbent[Rh catalyst].particle_void_fraction", Num(0.20e-3 * kPelletDensity)},
      {"layer[Bed].length", Num(kLength)},
      {"layer[Bed].cells", Num(c.cells)},
      {"layer[Bed].bed_void_fraction", Num(BedVoidFraction())},
      {"layer[Bed].particle_diameter", Num(kParticleDiameter)},
      {"initial.pressure", Num(P)},
      {"initial.temperature", Num(kT)},
      {"step[Feed].time", Num(c.step_time)},
      {"step[Feed].bottom[Feed inlet].molar_flow", Num(FeedMolarFlow())},
      {"step[Feed].bottom[Feed inlet].temperature", Num(kT)},
      {"step[Feed].top.pressure", Num(P)},
  };
  for (const auto& [comp, y] : c.feed) v.push_back({"step[Feed].bottom[Feed inlet].y[" + comp + "]", Num(y)});

  // Start-up: the bed starts in the gas-phase equilibrium of the feed at the
  // run pressure. With an active sorbent its CO2 is replaced by Ar, so the
  // sorbent starts clean (README, "Start-up").
  auto y0 = Equilibrium(c.feed, kT, c.pressure_bar);
  if (c.sorbent_active) {
    y0[AR] += y0[CO2];
    y0[CO2] = 0.0;
  }
  for (const auto& [comp, y] : y0) v.push_back({"initial.y[" + comp + "]", Num(y)});

  if (c.sorbent_active) {
    // Wet Freundlich q = k p^(1/n) [mol/kg, p in bar] (ch. 9, Table 4) as the
    // engine's q = (a p)^(1/v) [p in Pa]: a = k^n 1e-5, v = n.
    const double k = 1.11, n = 1.76;
    v.push_back({"adsorbent[K-HTC].isotherm[CarbonDioxide].a", Num(std::pow(k, n) * 1e-5)});
    v.push_back({"adsorbent[K-HTC].isotherm[CarbonDioxide].v", Num(n)});
  }

  const std::string rx = "adsorbent[Rh catalyst].reaction[";
  for (const auto& r : kRates) {
    const std::string h = rx + r.name + "].";
    v.push_back({h + "basis", "partial_pressure"});
    v.push_back({h + "rate.pre_exp_factor", Num(r.k0)});
    v.push_back({h + "rate.activation_energy", Num(r.E)});
    for (const auto& [comp, o] : r.orders) v.push_back({h + "rate.order[" + comp + "]", Num(o)});
    v.push_back({h + "rate.equilibrium_pre_exp_factor", Num(r.KA)});
    v.push_back({h + "rate.equilibrium_energy", Num(r.Ktheta * kR)});
    v.push_back({h + "rate.denominator_exponent", "1"});
    for (size_t j = 0; j < terms.size(); ++j) {
      const std::string t = h + "rate.adsorption[" + std::to_string(j + 1) + "].";
      v.push_back({t + "pre_exp_factor", Num(terms[j].K0)});
      v.push_back({t + "energy", Num(terms[j].dH)});
      for (const auto& [comp, e] : terms[j].exponent) v.push_back({t + "exponent[" + comp + "]", Num(e)});
    }
    // Isothermal run: the heat of reaction does not enter.
    v.push_back({h + "heat_of_reaction", "0"});
  }
  for (const auto& [path, text] : v) sim.SetParameterText(path, text);
  for (const auto& r : kRates) {
    sim.SetParameterText(rx + r.name + "].rate.reversible", "1");
    sim.SetParameterText(rx + r.name + "].heat_mode", "constant");
  }
  if (c.sorbent_active) {
    // Ch. 8, Eq. 6: k = 15 eps_p D_p / (r_p^2 K), K = rho_p R T dq*/dp, with
    // D_p = 1.1e-6 m2/s, eps_p = 0.4 and tortuosity 1 (ch. 8, Table 12): the
    // macropore term alone with the isotherm-slope capacity factor.
    const std::string kp = "adsorbent[K-HTC].kinetics[CarbonDioxide].";
    sim.SetParameterText(kp + "resistances.macropore", "1");
    sim.SetParameterText(kp + "resistances.lumped", "0");
    sim.RefreshParameters();
    sim.SetParameterText(kp + "capacity_factor_mode", "slope");
    sim.RefreshParameters();
    sim.SetParameterText(kp + "tortuosity", "1");
    sim.SetParameterText(kp + "macropore_void_fraction", "0.4");
    sim.SetParameterText(kp + "molecular_diffusion_coefficient", "1.1e-6");
  }

  d.settings.css.method_ = CSSMethod::NONE;
  d.cycle.number_of_cycles_ = 1;
  d.bed.equations_.energy_balance_mode_ = EnergyBalanceMode::ISOTHERMAL;

  RunOptions options;
  options.cycles = 1;
  options.on_message = [](LogLevel, const std::string&) {};
  const auto t0 = std::chrono::steady_clock::now();
  const RunResult res = sim.Run(options);
  Outcome out;
  out.wall_s = std::chrono::duration<double>(std::chrono::steady_clock::now() - t0).count();
  out.ok = res.Ok();
  out.message = res.message;
  out.healthy = res.healthy;
  if (!out.ok) return out;
  for (const auto& comp : {CH4, H2O, H2, CO, CO2, AR}) {
    const auto e = sim.Log().GetComponentFaceVariable(0, 0, ComponentVariableType::MOLE_FRACTION, comp);
    if (!e.values.empty() && !e.values.back().empty()) out.y_out[comp] = e.values.back().back();
  }
  return out;
}

// CH4 conversion with Ar as the internal standard (ch. 9, p. 193) [%].
double Conversion(const std::map<std::string, double>& feed, const std::map<std::string, double>& y) {
  return 100.0 * (1.0 - (y.at(CH4) / y.at(AR)) / (feed.at(CH4) / feed.at(AR)));
}

// H2 on the dry, Ar-free basis [%].
double DryH2(const std::map<std::string, double>& y) {
  return 100.0 * y.at(H2) / (y.at(H2) + y.at(CH4) + y.at(CO) + y.at(CO2));
}

std::string Tag(double sc, double p) {
  char buf[32];
  std::snprintf(buf, sizeof(buf), "SC%g_P%g", sc, p);
  return buf;
}

// Tolerances (regression bounds, README "Measured errors"): each covers one
// quantity over the rows listed and is set from the largest measured error.
constexpr double kTolFig7Equilibrium = 1.0;    // pp; measured 0.46 (6 bar)
constexpr double kTolTable6Equilibrium = 2.5;  // pp; measured 1.95 (S/C 3, 1.5 bar)
constexpr double kTolSeConversion = 1.0;       // pp; measured <= 0.7 (S/C 3)
constexpr double kTolSeH2 = 0.3;               // pp; measured <= 0.2

}  // namespace

class Halabi2011Test : public ::testing::Test {};

// Case 1: inert sorbent. The outlet of the 4.2 cm catalyst bed reaches the
// gas-phase equilibrium of the chapter-8 rate model; the thesis prints the
// equilibrium conversions it computed for the same feeds.
TEST_F(Halabi2011Test, GasPhaseEquilibrium) {
  validation::CaseReport report(kPaper, "GasPhaseEquilibrium");
  for (double p : {1.5, 4.5, 6.0}) {
    Case c;
    c.feed = Fig7Feed();
    c.pressure_bar = p;
    c.sorbent_active = false;
    c.step_time = 20.0;
    const Outcome o = RunCase(c);
    ASSERT_TRUE(o.ok) << o.message;
    EXPECT_TRUE(o.healthy);
    const double x = Conversion(c.feed, o.y_out);
    char q[64];
    std::snprintf(q, sizeof(q), "Fig7_P%g", p);
    std::printf("[halabi_2011] Fig. 7 feed, %.1f bar: X = %.2f %% (%.2f s)\n", p, x, o.wall_s);
    report.CheckAbs(std::string("CH4 equilibrium conversion, Fig. 7 feed, ") + Num(p) + " bar", x,
                    validation::Reference(kPaper, "GasPhaseEquilibrium", q), kTolFig7Equilibrium);
  }
  for (double sc : {3.0, 4.0, 5.0, 6.0})
    for (double p : {1.5, 4.5, 6.0}) {
      Case c;
      c.feed = Table6Feed(sc);
      c.pressure_bar = p;
      c.sorbent_active = false;
      c.step_time = 20.0;
      const Outcome o = RunCase(c);
      ASSERT_TRUE(o.ok) << o.message;
      EXPECT_TRUE(o.healthy);
      const double x = Conversion(c.feed, o.y_out);
      std::printf("[halabi_2011] Table 6 feed S/C %g, %.1f bar: X = %.2f %% (%.2f s)\n", sc, p, x, o.wall_s);
      report.CheckAbs("CH4 equilibrium conversion, Table 6, S/C " + Num(sc) + ", " + Num(p) + " bar", x,
                      validation::Reference(kPaper, "GasPhaseEquilibrium", "Table6_" + Tag(sc, p)),
                      kTolTable6Equilibrium);
    }
}

// Case 2: the sorption-enhanced runs of Table 6 at 1.5 bar, read before CO2
// breakthrough (t = 60 s).
TEST_F(Halabi2011Test, SorptionEnhanced_1p5bar) {
  validation::CaseReport report(kPaper, "SorptionEnhanced_1p5bar");
  double wall = 0.0;
  for (double sc : {3.0, 4.0, 5.0, 6.0}) {
    Case c;
    c.feed = Table6Feed(sc);
    c.pressure_bar = 1.5;
    const Outcome o = RunCase(c);
    ASSERT_TRUE(o.ok) << o.message;
    EXPECT_TRUE(o.healthy);
    wall += o.wall_s;
    const double x = Conversion(c.feed, o.y_out), h2 = DryH2(o.y_out);
    std::printf("[halabi_2011] SE S/C %g, 1.5 bar: X = %.2f %%, H2 = %.2f %%, CO+CO2 = %.0f ppm dry (%.2f s)\n", sc,
                x, h2, 1e6 * (o.y_out.at(CO) + o.y_out.at(CO2)) / (o.y_out.at(H2) + o.y_out.at(CH4) +
                                                                   o.y_out.at(CO) + o.y_out.at(CO2)),
                o.wall_s);
    const std::string tag = Tag(sc, 1.5);
    report.CheckAbs("CH4 conversion, S/C " + Num(sc) + ", 1.5 bar", x,
                    validation::Reference(kPaper, "SorptionEnhanced", "X_" + tag), kTolSeConversion);
    report.CheckAbs("H2 dry Ar-free, S/C " + Num(sc) + ", 1.5 bar", h2,
                    validation::Reference(kPaper, "SorptionEnhanced", "H2_" + tag), kTolSeH2);
  }
  std::printf("[halabi_2011] four 1.5 bar runs: %.1f s\n", wall);
}

// The 4.5 and 6 bar rows: recorded, not asserted (README, "Attribution").
TEST_F(Halabi2011Test, SorptionEnhanced_HighPressure_RecordOnly) {
  validation::CaseReport report(kPaper, "SorptionEnhanced_HighPressure");
  for (double p : {4.5, 6.0})
    for (double sc : {3.0, 4.0, 5.0, 6.0}) {
      Case c;
      c.feed = Table6Feed(sc);
      c.pressure_bar = p;
      const Outcome o = RunCase(c);
      ASSERT_TRUE(o.ok) << o.message;
      const double x = Conversion(c.feed, o.y_out), h2 = DryH2(o.y_out);
      std::printf("[halabi_2011] SE S/C %g, %.1f bar: X = %.2f %%, H2 = %.2f %% (%.2f s)\n", sc, p, x, h2, o.wall_s);
      const std::string tag = Tag(sc, p);
      const std::string note = "recorded only: the Freundlich fit lies below the Table 3 capacities (README)";
      report.Info("CH4 conversion, S/C " + Num(sc) + ", " + Num(p) + " bar", x,
                  validation::Reference(kPaper, "SorptionEnhanced", "X_" + tag), "abs_error", note);
      report.Info("H2 dry Ar-free, S/C " + Num(sc) + ", " + Num(p) + " bar", h2,
                  validation::Reference(kPaper, "SorptionEnhanced", "H2_" + tag), "abs_error", note);
    }
}

}  // namespace adsim
