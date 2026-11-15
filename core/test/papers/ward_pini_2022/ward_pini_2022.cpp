// ============================================================================
//  Ward, Pini (2022), Efficient Bayesian optimization of industrial-scale
//  pressure-vacuum swing adsorption processes for CO2 capture, Ind. Eng. Chem.
//  Res. 61, 13650-13668, doi:10.1021/acs.iecr.2c02313.
//
//  Four-step CO2/N2 PVSA with feed pressurisation (adsorption, forward
//  blowdown, reverse evacuation, feed pressurisation). Each case is one of the
//  ten cost-optimal designs of the Supporting Information, Tables S2 (zeolite
//  13X) and S3 (ZIF-36-FRL), at p_L = 0.02, 0.03, 0.05, 0.07 and 0.1 bar. The
//  designs are fully specified (pressures, step times, feed velocity, column
//  length), so the deck is run once at each and compared with the paper's
//  purity, recovery, productivity and energy. 13X_PL003 and ZIF_PL007 are the
//  minimum-cost designs of the main text, Table 7; they are saved as the
//  example decks WardPini2022_13X_PL003 and WardPini2022_ZIF36FRL_PL007.
//
//  Measured errors, bounds and the attribution are in README.md.
// ============================================================================

#define _USE_MATH_DEFINES
#include <cmath>
#include <cstdio>
#include <limits>
#include <map>
#include <string>

#include "../../../src/pch.h"
#include "adsim/solver/coupled_solver.h"
#include "adsim/solver/kpis.h"

#include <gtest/gtest.h>

#include "ward_pini_2022.h"
#include "../../fixtures/validation_metrics.h"

namespace adsim {

  namespace {
    constexpr double kBar = 1.0e5;
    constexpr double kInf = std::numeric_limits<double>::infinity();
    constexpr const char* kPaper = "ward_pini_2022";

    // The paper's model runs WENO on 10 volume elements (Sect. 2.2).
    constexpr int kCells = 10;

    struct WardPiniKpis {
      double purity = 0.0;        // [-] CO2 in the evacuation product (Eq. 7)
      double recovery = 0.0;      // [-] (Eq. 8)
      double productivity = 0.0;  // [mol CO2 / (m3 column s)] (Eq. 9)
      double energy = 0.0;        // [kWh / t CO2] (Eq. 10)
      double step_energy[4] = {}; // [kWh / t CO2] per step: pres, ads, bd, evac
      int cycles = 0;
      bool css = false;
    };

    WardPiniKpis Evaluate(const UserData& data, const CoupledSolver& solver, const WardPiniDesign& d) {
      const EquationLog& eq_log = solver.GetLog();
      const auto& cyc = eq_log.sys_data.cycle_index_set_;
      const int last_cycle = cyc.empty() ? 0 : *cyc.rbegin();
      KPIParams kpi(eq_log, solver.GetBoundaryLog());
      kpi.product_comp_ = WardPini2022Test::kCO2;
      kpi.cycle_index_ = last_cycle;
      kpi.feed_endpoints = KpiEndpointsForFaces(data.cycle,
        { { 0, BoundaryLocation::BOTTOM }, { 1, BoundaryLocation::BOTTOM } });
      kpi.product_endpoints = KpiEndpointsForFaces(data.cycle, { { 3, BoundaryLocation::BOTTOM } });
      WardPiniKpis k;
      k.purity = CalculatePurity(kpi);
      k.recovery = CalculateRecovery(kpi);
      k.productivity = CalculateProductionRate(kpi) / WardPini2022Test::ColumnVolume(d);
      k.energy = CalculateSpecificEnergyConsumption(kpi) / 3600.0;  // J/kg -> kWh/t
      const double product_kg =
        CalculateNetComponentMolesInProduct(kpi) * eq_log.sys_data.molar_mass.at(WardPini2022Test::kCO2);
      for (int s = 0; s < 4; ++s) k.step_energy[s] = CalculateStepEnergy(kpi, s) / product_kg / 3600.0;
      k.cycles = last_cycle + 1;
      k.css = solver.GetCSSCycleIndex() >= 0;
      return k;
    }

    // SI Tables S2 and S3, in the order p_H, p_I, p_L, t_ads, t_bd, t_evac,
    // v_F, L. t_pres = 20 s (Sect. 3.3).
    using M = WardPiniMaterial;
    const std::map<std::string, WardPiniDesign> kDesigns = {
      { "13X_PL002", { M::ZEOLITE_13X, 3.81 * kBar, 0.879 * kBar, 0.02 * kBar, 50.3, 93.6, 111, 1.69, 1.57 } },
      { "13X_PL003", { M::ZEOLITE_13X, 4.72 * kBar, 0.781 * kBar, 0.03 * kBar, 67.5, 40.5, 174, 1.70, 3.56 } },
      { "13X_PL005", { M::ZEOLITE_13X, 8.24 * kBar, 1.22 * kBar, 0.05 * kBar, 51.0, 80.1, 186, 1.87, 4.18 } },
      { "13X_PL007", { M::ZEOLITE_13X, 9.01 * kBar, 1.68 * kBar, 0.07 * kBar, 21.1, 78.9, 113, 2.0, 1.74 } },
      { "13X_PL010", { M::ZEOLITE_13X, 8.32 * kBar, 1.31 * kBar, 0.10 * kBar, 37.7, 166, 135, 1.12, 2.22 } },
      { "ZIF_PL002", { M::ZIF_36_FRL, 3.66 * kBar, 0.765 * kBar, 0.02 * kBar, 99.1, 159, 154, 1.02, 2.35 } },
      { "ZIF_PL003", { M::ZIF_36_FRL, 3.15 * kBar, 0.478 * kBar, 0.03 * kBar, 48.2, 152, 151, 1.88, 2.97 } },
      { "ZIF_PL005", { M::ZIF_36_FRL, 4.65 * kBar, 0.969 * kBar, 0.05 * kBar, 33.1, 189, 159, 2.0, 2.12 } },
      { "ZIF_PL007", { M::ZIF_36_FRL, 5.52 * kBar, 1.00 * kBar, 0.07 * kBar, 44.1, 134, 160, 1.17, 2.39 } },
      { "ZIF_PL010", { M::ZIF_36_FRL, 9.22 * kBar, 1.59 * kBar, 0.10 * kBar, 43.1, 116, 192, 0.848, 2.01 } },
    };

    // The two minimum-cost designs of Table 7, shipped as example decks.
    const std::map<std::string, std::string> kExampleDecks = {
      { "13X_PL003", "WardPini2022_13X_PL003" },
      { "ZIF_PL007", "WardPini2022_ZIF36FRL_PL007" },
    };

    // Bounds: regression bounds from the measured errors (README.md "Measured
    // errors"), asserted where the error is within the 5 % target (5 pp for
    // purity and recovery); kInf records the comparison without asserting (the
    // energy gap and the 13X recovery at 0.02 bar are open, README
    // "Attribution"). Purity and recovery in percentage points, productivity
    // and energy relative.
    struct Tol { double purity, recovery, productivity, energy; };
    const std::map<std::string, Tol> kTol = {
      // measured: purity 0.374 pp, recovery 5.539 pp, productivity 0.0717, energy 0.0069
      { "13X_PL002", { 0.45, kInf, kInf, 0.012 } },
      // measured: 0.047, 2.258, 0.0336, 0.0726
      { "13X_PL003", { 0.1, 2.6, 0.04, kInf } },
      // measured: 0.263, 0.774, 0.0108, 0.1731
      { "13X_PL005", { 0.35, 0.9, 0.015, kInf } },
      // measured: 0.420, 1.776, 0.0203, 0.3669
      { "13X_PL007", { 0.5, 2.0, 0.025, kInf } },
      // measured: 0.401, 0.378, 0.0035, 0.2493
      { "13X_PL010", { 0.5, 0.45, 0.008, kInf } },
      // measured: 0.218, 0.344, 0.0041, 0.0694
      { "ZIF_PL002", { 0.3, 0.45, 0.008, kInf } },
      // measured: 0.164, 0.758, 0.0052, 0.1293
      { "ZIF_PL003", { 0.25, 0.9, 0.009, kInf } },
      // measured: 0.528, 0.108, 0.0025, 0.2488
      { "ZIF_PL005", { 0.65, 0.2, 0.006, kInf } },
      // measured: 0.380, 0.172, 0.0139, 0.2145
      { "ZIF_PL007", { 0.45, 0.25, 0.018, kInf } },
      // measured: 0.608, 1.237, 0.0196, 0.2565
      { "ZIF_PL010", { 0.75, 1.45, 0.024, kInf } },
    };

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
  }  // namespace

  void WardPini2022Test::RunDesign(const std::string& id) {
    const WardPiniDesign& d = kDesigns.at(id);
    const bool zif = d.material == WardPiniMaterial::ZIF_36_FRL;
    const char* material = zif ? "ZIF-36-FRL" : "zeolite 13X";
    char title[160], description[512];
    std::snprintf(title, sizeof title, "Ward & Pini 2022 - %s, p_L = %.2f bar", material, d.P_L / kBar);
    std::snprintf(description, sizeof description,
                  "Four-step CO2/N2 PVSA with feed pressurisation on %s: the cost-optimal design of Ward & "
                  "Pini (2022), SI Table %s, at an evacuation pressure of %.2f bar (p_H %.2f bar, p_I %.3g bar, "
                  "v_F %.3g m/s, L %.3g m).",
                  material, zif ? "S3" : "S2", d.P_L / kBar, d.P_H / kBar, d.P_I / kBar, d.v_feed, d.length);
    BuildWardPini(title, description, d, kCells);
    CoupledSolver solver{ data };
    std::string error;
    ASSERT_TRUE(data.PerformChecks(error)) << error;
    ASSERT_GE(solver.IntegrateCycle(data.cycle, false), 0) << "Solver did not complete";
    const WardPiniKpis k = Evaluate(data, solver, d);
    std::printf("\n[KPI] Ward & Pini 2022 %s | cells %d | cycles %d (CSS %s) | purity %.4f %% | recovery %.4f %% | "
                "productivity %.5f mol/(m3 s) | energy %.3f kWh/t "
                "(pres %.1f, ads %.1f, bd %.1f, evac %.1f)\n\n",
                id.c_str(), kCells, k.cycles, k.css ? "yes" : "no", 100 * k.purity, 100 * k.recovery,
                k.productivity, k.energy, k.step_energy[0], k.step_energy[1], k.step_energy[2], k.step_energy[3]);
    EXPECT_TRUE(k.css);

    using validation::Reference;
    const Tol& tol = kTol.at(id);
    validation::CaseReport report(kPaper, id);
    Check(report, "CO2 purity", 100 * k.purity, Reference(kPaper, id, "CO2_purity"), tol.purity, false);
    Check(report, "CO2 recovery", 100 * k.recovery, Reference(kPaper, id, "CO2_recovery"), tol.recovery, false);
    Check(report, "Productivity", k.productivity, Reference(kPaper, id, "productivity"), tol.productivity, true,
          "per column volume");
    Check(report, "Energy", k.energy, Reference(kPaper, id, "energy"), tol.energy, true);
    if (auto it = kExampleDecks.find(id); it != kExampleDecks.end()) SaveResultsAsExample(it->second, solver);
  }

  TEST_F(WardPini2022Test, Zeolite13X_PL002) { RunDesign("13X_PL002"); }
  TEST_F(WardPini2022Test, Zeolite13X_PL003) { RunDesign("13X_PL003"); }
  TEST_F(WardPini2022Test, Zeolite13X_PL005) { RunDesign("13X_PL005"); }
  TEST_F(WardPini2022Test, Zeolite13X_PL007) { RunDesign("13X_PL007"); }
  TEST_F(WardPini2022Test, Zeolite13X_PL010) { RunDesign("13X_PL010"); }
  TEST_F(WardPini2022Test, ZIF36FRL_PL002) { RunDesign("ZIF_PL002"); }
  TEST_F(WardPini2022Test, ZIF36FRL_PL003) { RunDesign("ZIF_PL003"); }
  TEST_F(WardPini2022Test, ZIF36FRL_PL005) { RunDesign("ZIF_PL005"); }
  TEST_F(WardPini2022Test, ZIF36FRL_PL007) { RunDesign("ZIF_PL007"); }
  TEST_F(WardPini2022Test, ZIF36FRL_PL010) { RunDesign("ZIF_PL010"); }

} // namespace adsim
