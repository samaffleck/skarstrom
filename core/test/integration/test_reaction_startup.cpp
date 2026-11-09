/*
 * Skarstrom - pressure and temperature swing adsorption simulator
 * Copyright (c) 2024-2026 Sam Affleck and contributors
 *
 * Licensed under the terms in the LICENSE file at the repository root.
 */
// Start-up of a very fast reversible reaction from a product-free bed.
//
// An isothermal packed bed of a steam-reforming catalyst with the LHHW rates
// of Xu and Froment (1989, AIChE J. 35:88; three reversible reactions, H2
// orders -2.5, -1 and -3.5, a four-term denominator squared), fed with
// 3.94 % CH4, steam to carbon 4, 0.1 % H2 in argon at 450 C and 1.5 bar
// (75 nml/min, 9 mm ID, 4.2 cm bed). At these conditions the catalyst can
// convert some 10^4 times the methane fed: the gas reaches equilibrium within
// microseconds of meeting the catalyst. The bed starts with steam and argon
// only (no product), so the H2-poor front that moves into the bed carries
// rates with p_H2^-2.5 near zero H2.
//
//  (a) the run from the product-free bed reaches 1200 s, and its outlet CH4
//      conversion before breakthrough of the feed front equals that of the
//      run from a bed that starts at the gas-phase equilibrium of the feed,
//      and the equilibrium conversion itself, to < 1 point;
//  (b) the reaction start-up hint (REACTION_STIFF) names the reaction and the
//      near-zero component when the step size collapses at a reacting cell
//      (an artificially stiff case: the Jacobian's reaction treatment off,
//      BandedJacobian::SetReactionIncrementForTesting), and stays silent on
//      the runs of (a).

#include <chrono>
#include <cmath>
#include <cstdio>
#include <iostream>
#include <map>
#include <string>
#include <vector>

#include <gtest/gtest.h>

#include "adsim/equations/reaction_kinds.h"
#include "adsim/objects/deck_builder.h"
#include "adsim/objects/reaction.h"
#include "adsim/simulation.h"
#include "adsim/solver/jacobian.h"
#include "adsim/solver/run_health.h"

namespace adsim {
namespace {

const std::string CH4 = "Methane", H2O = "Water", CO = "CarbonMonoxide", H2 = "Hydrogen", CO2 = "CarbonDioxide",
                  AR = "Argon";

constexpr double kR = 8.314462618;
constexpr double kT = 723.15;   // [K]
constexpr double kP = 1.5e5;    // [Pa]

std::string Num(double v) {
  char buf[40];
  std::snprintf(buf, sizeof(buf), "%.17g", v);
  return buf;
}

// Feed: 3.94 % CH4, S/C = 4, 0.1 % H2, argon balance.
const std::map<std::string, double> kFeed = {
    {CH4, 0.0394}, {H2O, 4 * 0.0394}, {H2, 0.001}, {AR, 1.0 - 5 * 0.0394 - 0.001}};

// Xu-Froment equilibrium constants, K = exp(a + b / T) = A exp(-E / (R T)):
//   K1 (CH4 + H2O = CO + 3 H2)  exp(30.114 - 26830 / T)  [bar^2]
//   K2 (CO + H2O = CO2 + H2)    exp(-4.036 + 4400 / T)   [-]
//   K3 = K1 K2                                           [bar^2]
double K1(double T) { return std::exp(30.114 - 26830.0 / T); }
double K2(double T) { return std::exp(-4.036 + 4400.0 / T); }

// The gas-phase equilibrium of a feed (mole fractions) at (T, P): two
// independent extents per mole of feed, by a damped Newton on the logs of the
// two equilibrium conditions.
std::map<std::string, double> Equilibrium(const std::map<std::string, double>& feed, double T, double P) {
  auto y0 = [&](const std::string& c) { auto it = feed.find(c); return it == feed.end() ? 0.0 : it->second; };
  const double pbar = P * 1e-5;
  auto composition = [&](double x1, double x2) {
    std::map<std::string, double> n = {{CH4, y0(CH4) - x1}, {H2O, y0(H2O) - x1 - x2}, {CO, y0(CO) + x1 - x2},
                                       {H2, y0(H2) + 3 * x1 + x2}, {CO2, y0(CO2) + x2}, {AR, y0(AR)}};
    const double total = 1.0 + 2 * x1;
    for (auto& [c, v] : n) v /= total;
    return n;
  };
  auto residual = [&](double x1, double x2, double& f1, double& f2) {
    auto y = composition(x1, x2);
    f1 = std::log(y[CO]) + 3 * std::log(y[H2]) + 2 * std::log(pbar) - std::log(y[CH4]) - std::log(y[H2O]) -
         std::log(K1(T));
    f2 = std::log(y[CO2]) + std::log(y[H2]) - std::log(y[CO]) - std::log(y[H2O]) - std::log(K2(T));
  };
  double x1 = 0.5 * y0(CH4), x2 = 0.25 * y0(CH4);
  for (int it = 0; it < 200; ++it) {
    double f1, f2;
    residual(x1, x2, f1, f2);
    if (std::abs(f1) + std::abs(f2) < 1e-13) break;
    const double h = 1e-9;
    double a1, a2, b1, b2;
    residual(x1 + h, x2, a1, a2);
    residual(x1, x2 + h, b1, b2);
    const double j11 = (a1 - f1) / h, j21 = (a2 - f2) / h, j12 = (b1 - f1) / h, j22 = (b2 - f2) / h;
    const double det = j11 * j22 - j12 * j21;
    double d1 = -(f1 * j22 - f2 * j12) / det, d2 = -(j11 * f2 - j21 * f1) / det;
    // Keep every amount positive.
    double s = 1.0;
    for (int k = 0; k < 60; ++k) {
      auto y = composition(x1 + s * d1, x2 + s * d2);
      bool ok = true;
      for (const auto& [c, v] : y) ok = ok && v > 0.0;
      if (ok) break;
      s *= 0.5;
    }
    x1 += s * d1;
    x2 += s * d2;
  }
  return composition(x1, x2);
}

struct Options {
  bool equilibrium_start = false;   // bed gas at the feed's equilibrium, else steam/argon only
  double rate_scale = 1.0;          // multiplies the three rate constants
  int cells = 30;
  double step_time = 1200.0;        // [s]
  double velocity = 0.0;            // superficial [m/s]; 0: the 75 nml/min of the feed
};

// 75 nml/min (273.15 K, 1.01325 bar) as a superficial velocity at (T, P).
double FeedVelocity() {
  const double area = 0.25 * M_PI * 0.009 * 0.009;
  const double q_normal = 75e-6 / 60.0;  // [m3/s]
  return q_normal * (kT / 273.15) * (1.01325e5 / kP) / area;
}

struct Reformer {
  Simulation sim;
  AdsorbentID cat = 0;
  ReactionID smr = 0, wgs = 0, glob = 0;
};

void Build(Reformer& r, const Options& o) {
  UserData& d = r.sim.Data();
  for (const auto& c : {CH4, H2O, CO, H2, CO2, AR}) build::AddComponent(d, c);
  r.cat = build::AddAdsorbent(d, "Catalyst");  // inert isotherms
  build::SetBedShape(d, "cylinder");
  build::SetWallMaterial(d, "ss316l");
  build::SetWallBoundary(d, "zero_flux");
  const size_t layer = build::AddLayer(d, "Bed");
  build::SetLayerComposition(d, layer, {{r.cat, 1.0}});
  build::SetMomentum(d, layer, "ergun");
  build::SetDispersion(d, layer, "constant");
  build::SetHeatDispersion(d, layer, "constant");
  build::SetInternalHeatTransfer(d, layer, "constant");
  build::SetExternalHeatTransfer(d, layer, "constant");

  const size_t feed = build::AddStep(d, "Feed");
  build::SetFlowBoundary(d, feed, "bottom", "inlet");
  build::AddStream(d, feed, "bottom", "Feed inlet", "velocity", "constant");
  build::SetPressureBoundary(d, feed, "top", "outlet", "constant");

  r.smr = build::AddReaction(d, r.cat, "SMR", {{CH4, -1}, {H2O, -1}, {CO, 1}, {H2, 3}}, "lhhw");
  r.wgs = build::AddReaction(d, r.cat, "WGS", {{CO, -1}, {H2O, -1}, {CO2, 1}, {H2, 1}}, "lhhw");
  r.glob = build::AddReaction(d, r.cat, "Global", {{CH4, -1}, {H2O, -2}, {CO2, 1}, {H2, 4}}, "lhhw");
  for (ReactionID id : {r.smr, r.wgs, r.glob}) {
    build::SetReactionBasis(d, r.cat, id, "partial_pressure");
    for (int j = 0; j < 4; ++j) build::AddAdsorptionTerm(d, r.cat, id);
  }
  r.sim.RefreshParameters();

  const double u = o.velocity > 0.0 ? o.velocity : FeedVelocity();
  std::vector<std::pair<std::string, std::string>> values = {
      {"bed.diameter", "0.009"},
      {"adsorbent[Catalyst].solid_density", "1500"},
      {"adsorbent[Catalyst].heat_capacity", "900"},
      {"adsorbent[Catalyst].particle_void_fraction", "0.4"},
      {"layer[Bed].length", "0.042"},
      {"layer[Bed].cells", Num(o.cells)},
      {"layer[Bed].bed_void_fraction", "0.4"},
      {"layer[Bed].particle_diameter", "3e-4"},
      {"initial.pressure", Num(kP)},
      {"initial.temperature", Num(kT)},
      {"step[Feed].time", Num(o.step_time)},
      {"step[Feed].bottom[Feed inlet].velocity", Num(u)},
      {"step[Feed].bottom[Feed inlet].temperature", Num(kT)},
      {"step[Feed].top.pressure", Num(kP)},
  };
  for (const auto& [c, y] : kFeed) values.push_back({"step[Feed].bottom[Feed inlet].y[" + c + "]", Num(y)});
  if (o.equilibrium_start) {
    for (const auto& [c, y] : Equilibrium(kFeed, kT, kP)) values.push_back({"initial.y[" + c + "]", Num(y)});
  } else {
    values.push_back({"initial.y[" + H2O + "]", Num(kFeed.at(H2O))});
    values.push_back({"initial.y[" + AR + "]", Num(1.0 - kFeed.at(H2O))});
  }

  // Xu and Froment (1989), Table 3 (rate constants, kmol/(kg h) -> mol/(kg s))
  // and Table 4 (adsorption constants); E in J/mol.
  const double h = 1000.0 / 3600.0;
  struct Rate { const char* name; double A, E, A_K, E_K; std::map<std::string, double> order; };
  const double A_K1 = std::exp(30.114), E_K1 = 26830.0 * kR;
  const double A_K2 = std::exp(-4.036), E_K2 = -4400.0 * kR;
  const std::vector<Rate> rates = {
      {"SMR", 4.225e15 * h, 240.1e3, A_K1, E_K1, {{CH4, 1}, {H2O, 1}, {H2, -2.5}, {CO, 0}}},
      {"WGS", 1.955e6 * h, 67.13e3, A_K2, E_K2, {{CO, 1}, {H2O, 1}, {H2, -1}, {CO2, 0}}},
      {"Global", 1.020e15 * h, 243.9e3, A_K1 * A_K2, E_K1 + E_K2, {{CH4, 1}, {H2O, 2}, {H2, -3.5}, {CO2, 0}}},
  };
  // DEN = 1 + K_CO p_CO + K_H2 p_H2 + K_CH4 p_CH4 + K_H2O p_H2O / p_H2.
  struct Term { double A, E; std::map<std::string, double> exponent; };
  const std::vector<Term> terms = {
      {8.23e-5, -70.65e3, {{CO, 1}}},
      {6.12e-9, -82.90e3, {{H2, 1}}},
      {6.65e-4, -38.28e3, {{CH4, 1}}},
      {1.77e5, 88.68e3, {{H2O, 1}, {H2, -1}}},
  };
  for (const auto& rate : rates) {
    const std::string head = std::string("adsorbent[Catalyst].reaction[") + rate.name + "].";
    values.push_back({head + "rate.pre_exp_factor", Num(rate.A * o.rate_scale)});
    values.push_back({head + "rate.activation_energy", Num(rate.E)});
    values.push_back({head + "rate.equilibrium_pre_exp_factor", Num(rate.A_K)});
    values.push_back({head + "rate.equilibrium_energy", Num(rate.E_K)});
    values.push_back({head + "rate.denominator_exponent", "2"});
    for (const auto& [c, n] : rate.order) values.push_back({head + "rate.order[" + c + "]", Num(n)});
    for (size_t j = 0; j < terms.size(); ++j) {
      const std::string t = head + "rate.adsorption[" + std::to_string(j + 1) + "].";
      values.push_back({t + "pre_exp_factor", Num(terms[j].A)});
      values.push_back({t + "energy", Num(terms[j].E)});
      for (const auto& [c, e] : terms[j].exponent) values.push_back({t + "exponent[" + c + "]", Num(e)});
    }
  }
  for (const auto& [path, text] : values) r.sim.SetParameterText(path, text);

  d.settings.advection_scheme_ = AdvectionScheme::Upwind;
  d.settings.css.method_ = CSSMethod::NONE;
  d.cycle.number_of_cycles_ = 1;
  d.bed.equations_.energy_balance_mode_ = EnergyBalanceMode::ISOTHERMAL;
}

struct Outcome {
  RunResult result;
  double seconds = 0.0;
  std::vector<double> time;          // logged times of the step [s]
  std::vector<double> conversion;    // outlet CH4 conversion at those times [-]
};

Outcome RunReformer(Reformer& r) {
  RunOptions options;
  options.cycles = 1;
  options.on_message = [](LogLevel, const std::string&) {};
  Outcome out;
  const auto t0 = std::chrono::steady_clock::now();
  out.result = r.sim.Run(options);
  out.seconds = std::chrono::duration<double>(std::chrono::steady_clock::now() - t0).count();
  if (!out.result.Ok()) return out;
  // Outlet conversion from the carbon in the gas: (CO + CO2) / (CH4 + CO + CO2).
  const auto ych4 = r.sim.Log().GetComponentVariable(0, 0, ComponentVariableType::MOLE_FRACTION, CH4);
  const auto yco = r.sim.Log().GetComponentVariable(0, 0, ComponentVariableType::MOLE_FRACTION, CO);
  const auto yco2 = r.sim.Log().GetComponentVariable(0, 0, ComponentVariableType::MOLE_FRACTION, CO2);
  out.time = ych4.times;
  for (size_t k = 0; k < ych4.values.size(); ++k) {
    const size_t n = ych4.values[k].size() - 2;  // the last bed cell
    const double c = ych4.values[k][n], a = yco.values[k][n], b = yco2.values[k][n];
    out.conversion.push_back(c + a + b > 1e-6 ? (a + b) / (c + a + b) : std::nan(""));
  }
  return out;
}

// Mean outlet conversion over [t0, t1].
double MeanConversion(const Outcome& o, double t0, double t1) {
  double sum = 0.0;
  int n = 0;
  for (size_t k = 0; k < o.time.size(); ++k)
    if (o.time[k] >= t0 && o.time[k] <= t1 && std::isfinite(o.conversion[k])) { sum += o.conversion[k]; ++n; }
  return n ? sum / n : std::nan("");
}

void Report(const char* label, const Outcome& o) {
  const auto& s = o.result.stats;
  std::cout << label << ": status " << static_cast<int>(o.result.status) << " (" << o.result.message << "), "
            << o.seconds << " s wall, steps " << s.num_steps << ", res " << s.num_res_evals << ", jac "
            << s.num_jac_evals << ", conv fails " << s.num_conv_fails << ", err fails " << s.num_err_test_fails
            << ", recovered " << s.recovered_steps << "\n";
  for (const auto& w : o.result.health_warnings) std::cout << "  " << w.message << "\n";
}

const RunHealthWarning* FindWarning(const RunResult& r, RunHealthCheck check) {
  for (const auto& w : r.health_warnings)
    if (w.check == check) return &w;
  return nullptr;
}

// Restores the reaction increment at scope exit.
struct ReactionIncrementOff {
  ReactionIncrementOff() { BandedJacobian::SetReactionIncrementForTesting(false); }
  ~ReactionIncrementOff() { BandedJacobian::SetReactionIncrementForTesting(true); }
};

// ---- the powers the Jacobian reads -------------------------------------------------

// MinExponent over orders, reverse exponents and LHHW denominator exponents,
// on the deck's Xu-Froment rates; the activity floor's mole fraction.
TEST(ReactionStartup, MinExponentAndActivityFloor) {
  Reformer r;
  Build(r, Options{});
  const auto& src = r.sim.Data().adsorbent_directory.Get(r.cat).reaction_source_;
  const Reaction& smr = src.Get(r.smr);
  const Reaction& wgs = src.Get(r.wgs);
  // SMR: H2 order -2.5, reverse 0.5, denominator -1 (K_H2O p_H2O / p_H2).
  EXPECT_DOUBLE_EQ(smr.rate_->MinExponent(H2, smr.nu_), -2.5);
  // CO: reverse 1, denominator 1; H2O: order 1, denominator 1; Ar: absent.
  EXPECT_DOUBLE_EQ(smr.rate_->MinExponent(CO, smr.nu_), 1.0);
  EXPECT_DOUBLE_EQ(smr.rate_->MinExponent(H2O, smr.nu_), 1.0);
  EXPECT_TRUE(std::isinf(smr.rate_->MinExponent(AR, smr.nu_)));
  // WGS: H2 order -1, reverse 0 (skipped), denominator -1.
  EXPECT_DOUBLE_EQ(wgs.rate_->MinExponent(H2, wgs.nu_), -1.0);
  // Power law: orders only; the expression kind: 1 when named.
  PowerLawRate power;
  power.order_.emplace(CH4, 0.5);
  EXPECT_DOUBLE_EQ(power.MinExponent(CH4, {}), 0.5);
  EXPECT_TRUE(std::isinf(power.MinExponent(H2, {})));
  ExpressionRate expr;
  expr.expression_.SetExpr("1e-3 * p_Methane / 1e5");
  EXPECT_DOUBLE_EQ(expr.MinExponent(CH4, {}), 1.0);
  EXPECT_TRUE(std::isinf(expr.MinExponent(H2, {})));
  // 1e-10 bar at 1.5 bar; 1e-6 mol/m3 at 1.5 bar and 723.15 K.
  EXPECT_NEAR(ActivityFloorMoleFraction(RateBasis::PartialPressure, kT, kP), 1e-10 / 1.5, 1e-25);
  EXPECT_NEAR(ActivityFloorMoleFraction(RateBasis::Concentration, kT, kP), 1e-6 * kR * kT / kP, 1e-6 * 1e-6 * kR * kT / kP);
}

// ---- (a) -----------------------------------------------------------------------

TEST(ReactionStartup, ProductFreeBedMatchesTheEquilibriumStart) {
  Reformer cold, warm;
  Options o;
  Build(cold, o);
  o.equilibrium_start = true;
  Build(warm, o);
  const Outcome a = RunReformer(cold);
  Report("product-free bed", a);
  const Outcome b = RunReformer(warm);
  Report("equilibrium start", b);
  ASSERT_TRUE(a.result.Ok()) << a.result.message;
  ASSERT_TRUE(b.result.Ok()) << b.result.message;

  const auto eq = Equilibrium(kFeed, kT, kP);
  const double x_eq = (eq.at(CO) + eq.at(CO2)) / (eq.at(CH4) + eq.at(CO) + eq.at(CO2));
  const double xa = MeanConversion(a, 600.0, 1200.0), xb = MeanConversion(b, 600.0, 1200.0);
  std::cout << "outlet CH4 conversion 600-1200 s: product-free " << xa << ", equilibrium start " << xb
            << ", equilibrium " << x_eq << "\n";
  EXPECT_NEAR(xa, xb, 0.01);
  EXPECT_NEAR(xa, x_eq, 0.01);
  // Both runs integrate at the flow's time scale: no step collapse.
  EXPECT_EQ(FindWarning(a.result, RunHealthCheck::ReactionStiff), nullptr);
  EXPECT_EQ(FindWarning(b.result, RunHealthCheck::ReactionStiff), nullptr);
  EXPECT_LT(a.result.stats.num_steps, 2000);
}

// ---- (b) -----------------------------------------------------------------------

// The product-free start with IDA's own increment on the reaction columns
// (the treatment off) and ten times faster rates: the step collapses to the
// reaction's time scale, and the hint names a reaction and H2.
TEST(ReactionStartup, StepCollapseNamesTheReactionAndTheNearZeroComponent) {
  ReactionIncrementOff off;
  Reformer r;
  Options o;
  o.rate_scale = 10.0;
  o.step_time = 0.2;
  Build(r, o);
  const Outcome a = RunReformer(r);
  Report("stiff, reaction increment off", a);
  const RunHealthWarning* w = FindWarning(a.result, RunHealthCheck::ReactionStiff);
  ASSERT_NE(w, nullptr);
  EXPECT_STREQ(w->code.c_str(), "REACTION_STIFF");
  EXPECT_EQ(w->component, H2);
  EXPECT_TRUE(w->reaction == "SMR" || w->reaction == "WGS" || w->reaction == "Global") << w->reaction;
  EXPECT_GE(w->value, w->threshold);
  EXPECT_EQ(w->step_name, "Feed");
  EXPECT_NE(w->message.find(w->reaction), std::string::npos);
  EXPECT_NE(w->message.find(H2), std::string::npos);
  EXPECT_FALSE(a.result.healthy);
}

// The same collapse at a hundred times the rates, with a low step limit:
// the run fails with "too much work" (the failure of the original report),
// and the failed run carries the hint.
TEST(ReactionStartup, FailedRunCarriesTheHint) {
  ReactionIncrementOff off;
  Reformer r;
  Options o;
  o.rate_scale = 100.0;
  o.step_time = 10.0;
  Build(r, o);
  r.sim.Data().settings.max_iters_per_step = 500;
  const Outcome a = RunReformer(r);
  Report("stiff x100, reaction increment off, 500 steps", a);
  EXPECT_FALSE(a.result.Ok());
  const RunHealthWarning* w = FindWarning(a.result, RunHealthCheck::ReactionStiff);
  ASSERT_NE(w, nullptr);
  EXPECT_EQ(w->component, H2);
  EXPECT_FALSE(w->reaction.empty());
}

// Timings for the record (the PR and docs/reactions.md): the product-free
// start with the reaction treatment off and on, at the published rates.
TEST(ReactionStartup, DISABLED_TimingsWithAndWithoutTheReactionIncrement) {
  for (double scale : {1.0, 10.0, 100.0}) {
    for (bool on : {false, true}) {
      BandedJacobian::SetReactionIncrementForTesting(on);
      Reformer r;
      Options o;
      o.rate_scale = scale;
      Build(r, o);
      const Outcome a = RunReformer(r);
      const std::string label = "rates x" + std::to_string(static_cast<int>(scale)) +
                                (on ? ", reaction increment on" : ", reaction increment off");
      Report(label.c_str(), a);
      if (a.result.Ok()) std::cout << "  conversion 600-1200 s " << MeanConversion(a, 600.0, 1200.0) << "\n";
    }
  }
  BandedJacobian::SetReactionIncrementForTesting(true);
}

}  // namespace
}  // namespace adsim
