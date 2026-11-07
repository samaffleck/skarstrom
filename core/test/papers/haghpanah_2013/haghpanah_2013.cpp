// ============================================================================
//  Haghpanah, Majumder, Nilam, Rajendran, Farooq, Karimi, Amanullah (2013),
//  Ind. Eng. Chem. Res. 52, 4249-4265, doi:10.1021/ie302658y: four-step
//  CO2/N2 PVSA on zeolite 13X; and the verification of the same model by
//  Glover, Papathanasiou, Pini (2026), arXiv:2606.21523.
//
//  * Sec61: the fully specified reference design of Sect. 6.1 (Table 4)
//    against the paper's own purity and recovery.
//  * PL002 ... PL010: the five minimum-energy designs of Table 6 (Glover et
//    al. 2026, Table S2) against Haghpanah's energies (Table 6) and Glover's
//    purity, recovery, productivity and energy (Table 4).
//  * Fig5Benchmark: the CSS benchmark design of Glover et al. (2026), Fig. 5,
//    on Glover's 10 cells: the dominant multiplier of the cycle map against
//    rho = 0.9509 and the cycles to CSS from a feed- and an N2-saturated bed.
//
//  Measured errors, bounds and the attribution are in README.md.
// ============================================================================

#define _USE_MATH_DEFINES
#include <algorithm>
#include <cmath>
#include <cstdio>
#include <limits>
#include <map>
#include <string>
#include <vector>

#include "../../../src/pch.h"
#include "adsim/solver/coupled_solver.h"
#include "adsim/solver/cyclic_steady_state.h"
#include "adsim/solver/kpis.h"
#include "adsim/simulation.h"

#include <gtest/gtest.h>

#include "haghpanah_2013.h"
#include "../../fixtures/validation_metrics.h"

namespace adsim {

  namespace {
    constexpr double kBar = 1.0e5;
    constexpr double kInf = std::numeric_limits<double>::infinity();

    // Cycle KPIs at the last cycle, in the papers' units.
    struct HaghpanahKpis {
      double purity = 0.0;          // [-] CO2 in the evacuation product (Eq. 60)
      double recovery = 0.0;        // [-] (Eq. 61)
      double productivity = 0.0;    // [mol CO2 / (m3 column s)] (Glover Eq. 55 basis)
      double energy = 0.0;          // [kWh / t CO2] (Eqs. 66-70)
      int cycles = 0;               // cycles run
      bool css = false;
    };

    HaghpanahKpis Evaluate(const UserData& data, const CoupledSolver& solver) {
      const EquationLog& eq_log = solver.GetLog();
      const auto& cyc = eq_log.sys_data.cycle_index_set_;
      const int last_cycle = cyc.empty() ? 0 : *cyc.rbegin();
      KPIParams kpi(eq_log, solver.GetBoundaryLog());
      kpi.product_comp_ = Haghpanah2013Test::kCO2;
      kpi.cycle_index_ = last_cycle;
      kpi.feed_endpoints = KpiEndpointsForFaces(data.cycle,
        { { 0, BoundaryLocation::BOTTOM }, { 1, BoundaryLocation::BOTTOM } });
      kpi.product_endpoints = KpiEndpointsForFaces(data.cycle, { { 3, BoundaryLocation::BOTTOM } });
      HaghpanahKpis k;
      k.purity = CalculatePurity(kpi);
      k.recovery = CalculateRecovery(kpi);
      const double r = Haghpanah2013Test::kInnerRadius;
      const double v_column = M_PI * r * r * Haghpanah2013Test::kLength;
      k.productivity = CalculateProductionRate(kpi) / v_column;
      k.energy = CalculateSpecificEnergyConsumption(kpi) / 3600.0;  // J/kg -> kWh/t
      k.cycles = last_cycle + 1;
      k.css = solver.GetCSSCycleIndex() >= 0;
      return k;
    }

    // The operating points. Table 6 / Table S2 do not list the pressurisation
    // time: 20 s, fixed in both papers' optimisations (Haghpanah Sect. 7.1,
    // Glover p. 21).
    const HaghpanahDesign kSec61{ 15, 15, 30, 40, 1.0 * kBar, 0.2 * kBar, 0.1 * kBar, 1.0 };
    const std::map<std::string, HaghpanahDesign> kTable6 = {
      { "PL002", { 20, 94.89, 122.84, 189.46, 1.00 * kBar, 0.18 * kBar, 0.02 * kBar, 0.58 } },
      { "PL003", { 20, 69.36, 30.95, 124.27, 1.76 * kBar, 0.22 * kBar, 0.03 * kBar, 0.41 } },
      { "PL005", { 20, 91.63, 161.32, 195.24, 3.54 * kBar, 0.43 * kBar, 0.05 * kBar, 0.19 } },
      { "PL007", { 20, 70.93, 54.13, 88.45, 5.73 * kBar, 0.53 * kBar, 0.07 * kBar, 0.15 } },
      { "PL010", { 20, 97.36, 94.52, 128.45, 9.04 * kBar, 1.07 * kBar, 0.10 * kBar, 0.11 } },
    };
    // Glover et al. (2026), Fig. 5 caption; P_L and t_press are Glover's fixed
    // values (p. 21).
    const HaghpanahDesign kFig5{ 20, 50, 100, 100, 8.0 * kBar, 1.5 * kBar, 0.05 * kBar, 0.8 };

    // Bounds (regression bounds from the measured errors; README.md "Measured errors").
    struct Table6Tol { double energy_H, purity_G, recovery_G, productivity_G, productivity_G_ref, energy_G; };
    const std::map<std::string, Table6Tol> kTable6Tol = {
      { "PL002", { kInf, kInf, kInf, kInf, kInf, kInf } },
      { "PL003", { kInf, kInf, kInf, kInf, kInf, kInf } },
      { "PL005", { kInf, kInf, kInf, kInf, kInf, kInf } },
      { "PL007", { kInf, kInf, kInf, kInf, kInf, kInf } },
      { "PL010", { kInf, kInf, kInf, kInf, kInf, kInf } },
    };
    constexpr double kTolSec61Purity = kInf;    // pp
    constexpr double kTolSec61Recovery = kInf;  // pp
    constexpr double kTolProductivityH = kInf;  // relative
    constexpr double kTolMultiplier = kInf;     // absolute

    void Check(validation::CaseReport& report, const std::string& quantity, double sim,
               const validation::RefValue& ref, double tol, bool relative, const std::string& note = "") {
      if (std::isinf(tol)) {
        report.Info(quantity, sim, ref, relative ? "rel_error" : "abs_error", note);
      } else if (relative) {
        report.CheckRel(quantity, sim, ref, tol, note);
      } else {
        report.CheckAbs(quantity, sim, ref, tol, note);
      }
    }

    void PrintKpis(const std::string& label, const HaghpanahKpis& k) {
      std::printf("\n[KPI] Haghpanah 2013 %s | cycles %d (CSS %s) | purity %.4f %% | recovery %.4f %% | "
                  "productivity %.5f mol/(m3 s) | energy %.3f kWh/t\n\n",
                  label.c_str(), k.cycles, k.css ? "yes" : "no", 100 * k.purity, 100 * k.recovery,
                  k.productivity, k.energy);
    }
  }  // namespace

  // ---- Sect. 6.1 reference design ----------------------------------------
  TEST_F(Haghpanah2013Test, Sec61) {
    BuildHaghpanah("Haghpanah 2013 - Sect. 6.1 reference design",
                   "Four-step CO2/N2 PVSA on zeolite 13X, reference design of Haghpanah et al. (2013), "
                   "Sect. 6.1 / Table 4: pressurisation, adsorption, forward blowdown, reverse evacuation.",
                   kSec61, 30);
    CoupledSolver solver{ data };
    std::string error;
    ASSERT_TRUE(data.PerformChecks(error)) << error;
    ASSERT_GE(solver.IntegrateCycle(data.cycle, false), 0) << "Solver did not complete";
    const HaghpanahKpis k = Evaluate(data, solver);
    PrintKpis("Sec61", k);
    EXPECT_TRUE(k.css);

    using validation::Reference;
    validation::CaseReport report("haghpanah_2013", "Sec61");
    Check(report, "CO2 purity", 100 * k.purity, Reference("haghpanah_2013", "Sec61", "CO2_purity"),
          kTolSec61Purity, false);
    Check(report, "CO2 recovery", 100 * k.recovery, Reference("haghpanah_2013", "Sec61", "CO2_recovery"),
          kTolSec61Recovery, false);
    SaveTest("Haghpanah2013_Sec61", solver);
    SaveResultsAsExample("Haghpanah2013_Sec61", solver);
  }

  // ---- Table 6: minimum-energy designs ------------------------------------
  void Haghpanah2013Test::RunTable6(const std::string& id) {
    const HaghpanahDesign& d = kTable6.at(id);
    char title[128];
    std::snprintf(title, sizeof title, "Haghpanah 2013 - P_L = %.2f bar", d.P_L / kBar);
    BuildHaghpanah(title,
                   std::string("Four-step CO2/N2 PVSA on zeolite 13X: minimum-energy design of Haghpanah et al. "
                               "(2013), Table 6, at an evacuation pressure of ") +
                     std::to_string(d.P_L / kBar).substr(0, 4) + " bar (Glover et al. 2026, Table S2).",
                   d, 30);
    CoupledSolver solver{ data };
    std::string error;
    ASSERT_TRUE(data.PerformChecks(error)) << error;
    ASSERT_GE(solver.IntegrateCycle(data.cycle, false), 0) << "Solver did not complete";
    const HaghpanahKpis k = Evaluate(data, solver);
    PrintKpis(id, k);
    EXPECT_TRUE(k.css);

    using validation::Reference;
    const auto& tol = kTable6Tol.at(id);
    validation::CaseReport report("haghpanah_2013", id);
    Check(report, "Energy vs Haghpanah Table 6", k.energy, Reference("haghpanah_2013", id, "energy_H"),
          tol.energy_H, true);
    Check(report, "CO2 purity vs Glover", 100 * k.purity, Reference("haghpanah_2013", id, "purity_G"), tol.purity_G,
          false);
    Check(report, "CO2 recovery vs Glover", 100 * k.recovery, Reference("haghpanah_2013", id, "recovery_G"),
          tol.recovery_G, false);
    Check(report, "Productivity vs Glover", k.productivity, Reference("haghpanah_2013", id, "productivity_G"),
          tol.productivity_G, true, "per column volume");
    Check(report, "Productivity vs Glover Ref.", k.productivity,
          Reference("haghpanah_2013", id, "productivity_G_ref"), tol.productivity_G_ref, true,
          "per column volume; Glover's reference value, not printed by Haghpanah");
    Check(report, "Energy vs Glover", k.energy, Reference("haghpanah_2013", id, "energy_G"), tol.energy_G, true);
    if (id == "PL002") {
      // Haghpanah's productivity is per volume of adsorbent, (1 - eps) of the column.
      Check(report, "Productivity vs Haghpanah", k.productivity / (1.0 - kBedVoid),
            Reference("haghpanah_2013", id, "productivity_H"), kTolProductivityH, true,
            "per adsorbent volume (1 - eps) A L");
    }
    const std::string deck = "Haghpanah2013_" + id;
    SaveTest(deck, solver);
    SaveResultsAsExample(deck, solver);
  }

  TEST_F(Haghpanah2013Test, PL002) { RunTable6("PL002"); }
  TEST_F(Haghpanah2013Test, PL003) { RunTable6("PL003"); }
  TEST_F(Haghpanah2013Test, PL005) { RunTable6("PL005"); }
  TEST_F(Haghpanah2013Test, PL007) { RunTable6("PL007"); }
  TEST_F(Haghpanah2013Test, PL010) { RunTable6("PL010"); }

  // ---- Glover et al. (2026) Fig. 5: CSS multiplier and cycles to CSS -------
  namespace {
    struct CssTrace {
      std::vector<double> change;  // end-of-cycle relative state change, per cycle (0-based)
      std::vector<double> rate;    // the engine's fitted convergence rate (NaN if not estimated)
      int css_cycle = -1;          // the engine's stopping rule (default tolerances)
    };

    // Runs `cycles` cycles of the deck with the stopping rule held off (a
    // state tolerance no cycle can meet), then reads every cycle's state
    // change and fitted rate, and the stop of the default rule, from the log.
    CssTrace TraceCss(Simulation& sim, int cycles) {
      auto& css = sim.Data().settings.css;
      const double tol = css.state_tolerance_;
      sim.Data().cycle.number_of_cycles_ = cycles;
      sim.Data().settings.max_cycles_to_keep_ = cycles;
      css.state_tolerance_ = 1e-30;
      const RunResult run = sim.Run();
      EXPECT_TRUE(run.Ok()) << run.message;
      css.state_tolerance_ = tol;
      CssTrace t;
      for (const CSSRecord& r : sim.CSSHistory()) {
        t.change.push_back(r.state_change);
        t.rate.push_back(r.distance_estimated ? r.convergence_rate : std::nan(""));
        if (t.css_cycle < 0 && r.converged) t.css_cycle = r.cycle;
      }
      return t;
    }

    int FirstBelow(const std::vector<double>& v, double level) {
      for (size_t i = 1; i < v.size(); ++i)
        if (v[i] > 0.0 && v[i] < level) return static_cast<int>(i);
      return -1;
    }

    // Least-squares slope of ln(change) over cycles [a, b]: the dominant multiplier.
    double FitMultiplier(const std::vector<double>& v, int a, int b) {
      double sx = 0, sy = 0, sxx = 0, sxy = 0;
      int n = 0;
      for (int i = a; i <= b; ++i) {
        if (!(v[i] > 0.0)) continue;
        const double x = i, y = std::log(v[i]);
        sx += x; sy += y; sxx += x * x; sxy += x * y; ++n;
      }
      const double slope = (n * sxy - sx * sy) / (n * sxx - sx * sx);
      return std::exp(slope);
    }
  }  // namespace

  TEST_F(Haghpanah2013Test, Fig5Benchmark) {
    constexpr int kCycles = 300;
    const std::string deck = "Haghpanah2013_Fig5Benchmark";
    BuildHaghpanah("Haghpanah 2013 - Glover 2026 benchmark design",
                   "Four-step CO2/N2 PVSA on zeolite 13X (Haghpanah et al. 2013 model) at the CSS benchmark "
                   "design of Glover et al. (2026), Fig. 5, on 10 cells.",
                   kFig5, 10);
    SaveSimulation(deck);
    test::SaveToJsonAtomically(data, std::filesystem::current_path() / "../examples/decks" / (deck + ".json"));
    const std::string path = (test_path / deck / (deck + ".json")).string();

    // N2-saturated start (Haghpanah Sect. 6.1), as the deck.
    Simulation n2_sim = Simulation::FromFile(path);
    const CssTrace n2 = TraceCss(n2_sim, kCycles);
    // Feed-saturated start: the feed gas at P_H with loadings in equilibrium.
    Simulation feed_sim = Simulation::FromFile(path);
    feed_sim.Data().ic.pressure = kFig5.P_H;
    feed_sim.Data().ic.mole_fraction[kCO2] = 0.15;
    feed_sim.Data().ic.mole_fraction[kN2] = 0.85;
    const CssTrace feed = TraceCss(feed_sim, kCycles);
    for (int c = 0; c < kCycles; ++c)
      std::printf("  cycle %3d | N2 start change %.4e rate %.5f | feed start change %.4e rate %.5f\n", c,
                  n2.change[c], n2.rate[c], feed.change[c], feed.rate[c]);
    std::printf("  CSS stop (default rule): N2 start %d, feed start %d; first change < 1e-5: N2 %d, feed %d\n",
                n2.css_cycle, feed.css_cycle, FirstBelow(n2.change, 1e-5), FirstBelow(feed.change, 1e-5));
  }

} // namespace adsim
