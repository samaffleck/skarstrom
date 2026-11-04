/*
 * Skarstrom - pressure and temperature swing adsorption simulator
 * Copyright (c) 2024-2026 Sam Affleck and contributors
 *
 * Licensed under the terms in the LICENSE file at the repository root.
 */

// Cross-code regression: single-column breakthrough and purge cases run in two
// other open-source simulators, frozen as data, and replayed here.
//
// Each case is a deck in decks/ (written by the WP18 harness in
// skarstrom-gui benchmarks/open-source-psa/matrix/, saved before the run) and
// the reference tool's outlet curve in data/ (README.md: tools, versions,
// licences, the case specification and how every difference is explained).
// The test runs the deck, computes the same front metrics as the harness
// (metrics.py) and asserts that Skarstrom - reference stays within a bound set
// from the MEASURED, attributed difference:
//
//   f(t) = (y(t) - y_init) / (y_feed - y_init)            0 -> 1 for adsorption and purge
//   t_X  = first time f reaches X (X = 0.05, 0.5, 0.95)
//   rms  = RMS of f_sk - f_ref on a 1 s grid where either is in (0.01, 0.99)
//   T_ext = max (adsorption) / min (purge) outlet temperature
//
// Bounds (regression bounds on a measured, attributed difference, not
// acceptance criteria; the papers/ convention), from the harness run of the
// same deck on engine main:
//  * where the codes are taken to agree, |Skarstrom - reference| <=
//    1.25 |measured| + floor;
//  * where the difference is attributed to a named property choice (the
//    adiabatic cases, the viscosity in BT08), Skarstrom - reference =
//    measured +- (0.25 |measured| + floor);
// floors 0.25 s, 0.05 K; RMS <= 1.25 measured + 0.002; rounded up to 0.1 s /
// 0.01 K / 0.001. Every t50 band is checked to be below 0.5 % of t50: a 1 %
// error in the feed flux, the bed capacity or the voidage moves t50 by about
// 1 % and fails the test (the guard that the comparison is sensitive).
//
// What the measured differences are (README.md has the full attribution):
//  * isothermal cases, same grid and a WENO-type scheme in both codes: within
//    0.6 s on t50 and 3 s on t05/t95; the codes differ in their grids
//    (cell-centred finite volumes against the reference's nodes) and in the
//    viscosity (Skarstrom's polynomial; it matters only through the pressure
//    drop: BT08, +3.4 s);
//  * adiabatic cases: Skarstrom's CoolProp gas heat capacity and its constant
//    heat of adsorption -dU + R T0 (the reference: 30.7 J/(mol K) and the local
//    -(dU - R T)) delay the fronts by up to 14.5 s (BT06); with the reference's
//    property choices (skarstrom#45 switches) every case agrees within 4 s
//    at the same grid.

#include <gtest/gtest.h>

#include <algorithm>
#include <cmath>
#include <fstream>
#include <limits>
#include <sstream>
#include <stdexcept>
#include <string>
#include <vector>

#include "adsim/objects/deck_io.h"
#include "adsim/objects/user_data.h"
#include "adsim/solver/coupled_solver.h"
#include "adsim/solver/log.h"

namespace adsim::crosscode {
namespace {

struct Curve {
  std::vector<double> t, y, T;
};

struct Front {
  double t05, t50, t95, T_ext;
};

// Skarstrom - reference is expected to be `offset` within +-`band`
// (NaN: not asserted, e.g. t95 of a front still open at the end).
struct Bound {
  double offset, band;
};

constexpr double kNaN = std::numeric_limits<double>::quiet_NaN();

struct CrossCase {
  const char* id;
  const char* reference;   // data/<file>
  double y_init, y_feed;   // component 1 (CarbonDioxide)
  bool adsorption;
  Bound t05, t50, t95, T_ext;
  double max_rms;
};

const std::string kDir = "test/crosscode/";   // ctest runs in core/

std::string ReadFile(const std::string& path) {
  std::ifstream in(path);
  EXPECT_TRUE(in.good()) << "cannot open " << path;
  std::stringstream ss;
  ss << in.rdbuf();
  return ss.str();
}

Curve ReadCurve(const std::string& path) {
  Curve c;
  std::ifstream in(path);
  EXPECT_TRUE(in.good()) << "cannot open " << path;
  std::string line;
  std::getline(in, line);  // header t_s,y1_out,T_out_K
  while (std::getline(in, line)) {
    double a, b, d;
    char c1, c2;
    std::istringstream ls(line);
    if (ls >> a >> c1 >> b >> c2 >> d) {
      c.t.push_back(a);
      c.y.push_back(b);
      c.T.push_back(d);
    }
  }
  return c;
}

double Interp(const std::vector<double>& x, const std::vector<double>& y, double xi) {
  if (xi <= x.front()) return y.front();
  if (xi >= x.back()) return y.back();
  const auto it = std::upper_bound(x.begin(), x.end(), xi);
  const size_t i = static_cast<size_t>(it - x.begin());
  const double w = (xi - x[i - 1]) / (x[i] - x[i - 1]);
  return y[i - 1] + w * (y[i] - y[i - 1]);
}

std::vector<double> Fraction(const CrossCase& k, const std::vector<double>& y) {
  std::vector<double> f(y.size());
  for (size_t i = 0; i < y.size(); ++i) f[i] = (y[i] - k.y_init) / (k.y_feed - k.y_init);
  return f;
}

double Cross(const std::vector<double>& t, const std::vector<double>& f, double X) {
  for (size_t i = 0; i < f.size(); ++i) {
    if (f[i] >= X) {
      if (i == 0) return t[0];
      return t[i - 1] + (X - f[i - 1]) * (t[i] - t[i - 1]) / (f[i] - f[i - 1]);
    }
  }
  return std::nan("");
}

Front Metrics(const CrossCase& k, const Curve& c) {
  const auto f = Fraction(k, c.y);
  const double T_ext = k.adsorption ? *std::max_element(c.T.begin(), c.T.end())
                                    : *std::min_element(c.T.begin(), c.T.end());
  return {Cross(c.t, f, 0.05), Cross(c.t, f, 0.5), Cross(c.t, f, 0.95), T_ext};
}

double Rms(const CrossCase& k, const Curve& a, const Curve& b) {
  const auto fa = Fraction(k, a.y);
  const auto fb = Fraction(k, b.y);
  const double t_end = std::min(a.t.back(), b.t.back());
  double s = 0.0;
  int n = 0;
  for (double t = 0.0; t <= t_end + 1e-9; t += 1.0) {
    const double x = Interp(a.t, fa, t), y = Interp(b.t, fb, t);
    if ((x > 0.01 && x < 0.99) || (y > 0.01 && y < 0.99)) {
      s += (x - y) * (x - y);
      ++n;
    }
  }
  return n ? std::sqrt(s / n) : std::nan("");
}

Curve RunDeck(const std::string& id) {
  UserData data;
  LoadDeckFlexible(data, ReadFile(kDir + "decks/" + id + ".json"));
  std::string error;
  EXPECT_TRUE(data.PerformChecks(error)) << error;
  CoupledSolver solver{data};
  EXPECT_GE(solver.IntegrateCycle(data.cycle, false), 0) << id << ": solver did not complete";
  const auto& log = solver.GetLog();
  const auto y = log.GetComponentFaceVariable(0, 0, ComponentVariableType::MOLE_FRACTION, "CarbonDioxide");
  const auto T = log.GetFaceVariable(0, 0, VariableType::TEMPERATURE);
  Curve c;
  for (size_t i = 0; i < y.times.size(); ++i) {
    if (y.values[i].empty() || T.values[i].empty()) continue;
    if (!c.t.empty() && !(y.times[i] > c.t.back())) continue;
    c.t.push_back(y.times[i]);
    c.y.push_back(y.values[i].back());
    c.T.push_back(T.values[i].back());
  }
  return c;
}

void Check(const CrossCase& k) {
  SCOPED_TRACE(k.id);
  const Curve ref = ReadCurve(kDir + "data/" + k.reference);
  ASSERT_GT(ref.t.size(), 10u);
  const Curve sk = RunDeck(k.id);
  ASSERT_GT(sk.t.size(), 10u);
  ASSERT_NEAR(sk.t.back(), ref.t.back(), 1e-6) << "decks and references cover the same time";

  const Front r = Metrics(k, ref), s = Metrics(k, sk);
  const double rms = Rms(k, sk, ref);
  std::printf("%s  Skarstrom t05 %.2f t50 %.2f t95 %.2f T_ext %.3f | reference t05 %.2f t50 %.2f t95 %.2f "
              "T_ext %.3f | rms %.4f\n", k.id, s.t05, s.t50, s.t95, s.T_ext, r.t05, r.t50, r.t95, r.T_ext, rms);

  // Guards: the front is in the window and the comparison is sensitive to a
  // 1 % change in what sets the front time.
  ASSERT_TRUE(std::isfinite(r.t50) && std::isfinite(s.t50)) << "front not reached";
  EXPECT_LT(k.t50.band, 0.005 * r.t50) << "t50 band too wide to catch a 1 % error";
  EXPECT_GT(r.t50 - r.t05, 0.0);

  auto near = [](double a, double b, const Bound& bd, const char* what) {
    if (std::isnan(bd.band)) return;
    EXPECT_NEAR(a - b, bd.offset, bd.band) << what << ": Skarstrom " << a << ", reference " << b;
  };
  near(s.t05, r.t05, k.t05, "t05");
  near(s.t50, r.t50, k.t50, "t50");
  near(s.t95, r.t95, k.t95, "t95");
  near(s.T_ext, r.T_ext, k.T_ext, "outlet T extreme");
  EXPECT_LE(rms, k.max_rms) << "curve RMS";
}

}  // namespace

// Cases. BT* against the PEESE column model (Yancy-Caballero et al. 2020; MIT):
// same grid in both codes (100 cells; BT08 85, the finest the reference
// completes), Skarstrom with WENO, the reference with its own WENO faces.
// TP* against toPSAil (GPL-3.0; numerical outputs only).
const CrossCase kCases[] = {
  // BT01: isothermal 15 % CO2 (the original benchmark's case); same grid, WENO in both.
  //   measured: t05 -0.45 s, t50 -0.48 s, t95 -0.51 s, T_ext -0.00 K, rms 0.0057
  {"BT01", "BT01_peese_N100.csv", 0.0, 0.15, true,
   {0.0, 0.9}, {0.0, 0.9}, {0.0, 0.9}, {0.0, 0.06}, 0.010},
  // BT02: adiabatic 15 % CO2: t50 +10.8 s = CoolProp gas cp (-6.8 s when set to 30.7) + constant heat -dU + R T0 (-4.3 s when local); both matched (skarstrom#45): -0.5 s.
  //   measured: t05 +10.49 s, t50 +10.80 s, T_ext -0.15 K, rms 0.0227
  {"BT02", "BT02_peese_N100.csv", 0.0, 0.15, true,
   {+10.49, 2.9}, {+10.80, 3.0}, {kNaN, kNaN}, {-0.15, 0.09}, 0.031},
  // BT03: isothermal purge of a 15 % bed with N2 (desorption simple wave).
  //   measured: t05 -0.14 s, t50 +0.15 s, T_ext +0.00 K, rms 0.0000
  {"BT03", "BT03_peese_N100.csv", 0.15, 0.0, false,
   {0.0, 0.5}, {0.0, 0.5}, {kNaN, kNaN}, {0.0, 0.06}, 0.003},
  // BT04: adiabatic purge: outlet T_min -0.43 K = CoolProp gas cp (matched: -0.03 K).
  //   measured: t05 -0.33 s, t50 -0.31 s, T_ext -0.43 K, rms 0.0026
  {"BT04", "BT04_peese_N100.csv", 0.15, 0.0, false,
   {-0.33, 0.4}, {-0.31, 0.4}, {kNaN, kNaN}, {-0.43, 0.16}, 0.006},
  // BT05: isothermal 50 % CO2, pressure-form Langmuir (bulk feed: the velocity falls across the front).
  //   measured: t05 -0.04 s, t50 -0.09 s, t95 -0.14 s, T_ext -0.00 K, rms 0.0028
  {"BT05", "BT05_peese_N100.csv", 0.0, 0.5, true,
   {0.0, 0.4}, {0.0, 0.4}, {0.0, 0.5}, {0.0, 0.06}, 0.006},
  // BT06: adiabatic 50 % CO2: t50 +14.5 s = gas cp (-6.5 s) + heat form (-4.9 s); both matched: +2.9 s at 100 cells (+1.0 s at 400).
  //   measured: t05 +14.45 s, t50 +14.46 s, t95 +15.93 s, T_ext -0.87 K, rms 0.0524
  {"BT06", "BT06_peese_N100.csv", 0.0, 0.5, true,
   {+14.45, 3.9}, {+14.46, 3.9}, {+15.93, 4.3}, {-0.87, 0.27}, 0.068},
  // BT07: Henry isotherm, 1 % CO2 at 393 K: t05 / t95 differ by 3 s, the size of either code's change from 100 to 400 cells.
  //   measured: t05 +3.02 s, t50 +0.41 s, t95 -3.35 s, T_ext +0.00 K, rms 0.0021
  {"BT07", "BT07_peese_N100.csv", 0.0, 0.01, true,
   {0.0, 4.1}, {0.0, 0.8}, {0.0, 4.5}, {0.0, 0.05}, 0.005},
  // BT08: 0.3 mm particles (Ergun drop 25 % of P0), 85 cells (the reference fails at 100): t50 +2.8 s = Skarstrom's viscosity polynomial, 1.85e-5 against 1.72e-5 Pa s (matched: -0.6 s).
  //   measured: t05 +2.89 s, t50 +2.84 s, t95 +2.67 s, T_ext -0.00 K, rms 0.0298
  {"BT08", "BT08_peese_N85.csv", 0.0, 0.15, true,
   {+2.89, 1.0}, {+2.84, 1.0}, {+2.67, 1.0}, {-0.00, 0.06}, 0.040},
  // BT09: LiX, 40 % CO2 in CH4 (biogas), fast LDF (k = 46.7 1/s).
  //   measured: t05 -0.24 s, t50 -0.37 s, t95 -0.59 s, T_ext -0.00 K, rms 0.0301
  {"BT09", "BT09_peese_N100.csv", 0.0, 0.4, true,
   {0.0, 0.6}, {0.0, 0.8}, {0.0, 1.0}, {0.0, 0.06}, 0.040},
  // BT10: LiX biogas, adiabatic: t50 +1.8 s = gas cp (-0.6 s) + heat form (-1.1 s); both matched: +0.2 s.
  //   measured: t05 +2.04 s, t50 +1.83 s, t95 +1.75 s, T_ext +0.00 K, rms 0.0147
  {"BT10", "BT10_peese_N100.csv", 0.0, 0.4, true,
   {+2.04, 0.8}, {+1.83, 0.8}, {+1.75, 0.7}, {+0.00, 0.06}, 0.021},
  // TP01: toPSAil's model choices (no dispersion; CSTRs in series = first-order upwind), 100 cells in both.
  //   measured: t05 -0.51 s, t50 -0.51 s, t95 -0.55 s, T_ext +0.00 K, rms 0.0058
  {"TP01", "TP01_topsail_N100.csv", 0.0, 0.15, true,
   {0.0, 0.9}, {0.0, 0.9}, {0.0, 1.0}, {0.0, 0.05}, 0.010},
  // TP02: toPSAil, steel wall + ambient, 40 cells, h scaled to toPSAil's per-CSTR area: t50 +6.4 s = CoolProp gas cp (matched cp: +0.1 s).
  //   measured: t05 +7.19 s, t50 +6.40 s, T_ext -0.25 K, rms 0.0062
  {"TP02", "TP02_topsail_N40.csv", 0.0, 0.15, true,
   {+7.19, 2.1}, {+6.40, 1.9}, {kNaN, kNaN}, {-0.25, 0.12}, 0.010},

};

const CrossCase& Case(const std::string& id) {
  for (const auto& k : kCases)
    if (id == k.id) return k;
  throw std::out_of_range(id);
}

TEST(CrossCodePeese, BT01_Isothermal15Percent) { Check(Case("BT01")); }
TEST(CrossCodePeese, BT02_Adiabatic15Percent) { Check(Case("BT02")); }
TEST(CrossCodePeese, BT03_IsothermalPurge) { Check(Case("BT03")); }
TEST(CrossCodePeese, BT04_AdiabaticPurge) { Check(Case("BT04")); }
TEST(CrossCodePeese, BT05_IsothermalBulkLangmuirP) { Check(Case("BT05")); }
TEST(CrossCodePeese, BT06_AdiabaticBulkLangmuirP) { Check(Case("BT06")); }
TEST(CrossCodePeese, BT07_HenryTrace) { Check(Case("BT07")); }
TEST(CrossCodePeese, BT08_ErgunPressureDrop) { Check(Case("BT08")); }
TEST(CrossCodePeese, BT09_BiogasLiXIsothermal) { Check(Case("BT09")); }
TEST(CrossCodePeese, BT10_BiogasLiXAdiabatic) { Check(Case("BT10")); }

// ---------------------------------------------------------------------------
// Against toPSAil (GPL-3.0; numerical outputs only). The decks take toPSAil's
// model choices (README.md): no dispersion or conduction, first-order upwind
// on the same cells, one heat per component, negligible pressure drop.
// ---------------------------------------------------------------------------

TEST(CrossCodeToPSAil, TP01_IsothermalSameDiscretisation) { Check(Case("TP01")); }
TEST(CrossCodeToPSAil, TP02_WallAndAmbientHeatTransfer) { Check(Case("TP02")); }

}  // namespace adsim::crosscode
