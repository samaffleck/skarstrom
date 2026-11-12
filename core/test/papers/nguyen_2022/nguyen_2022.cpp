// STL includes
#define _USE_MATH_DEFINES
#include <cmath>
#include <filesystem>
#include <limits>
#include <stdexcept>
#include <string>
#include <vector>

// AdSim includes
#include "../../../src/pch.h"
#include "adsim/solver/coupled_solver.h"
#include "adsim/equations/isotherm.h"
#include "adsim/boundary_conditions/boundary_conditions.h"
#include "adsim/boundary_conditions/boundary_stream.h"
#include "adsim/boundary_conditions/step.h"
#include "adsim/equations/dispersion_model.h"
#include "adsim/equations/heat_transfer_model.h"

// GTest includes
#include <gtest/gtest.h>

#include "nguyen_2022.h"
#include "../../fixtures/validation_metrics.h"

namespace adsim {

  // ==========================================================================
  //  Nguyen, Lin, Shimizu & Rajendran (2022), "Separation of CO2 and N2 on a
  //  hydrophobic metal organic framework CALF-20", Chem. Eng. J. 442, 136263,
  //  doi:10.1016/j.cej.2022.136263.
  //  The six dynamic column breakthrough experiments of Fig. 4 (CO2/He,
  //  50 cm3/min) and Fig. 5 (CO2/N2, 100 cm3/min), each saved as an example
  //  deck. Every case is a prediction with the paper's parameters (Table 2,
  //  Table 3) and the conditions of the thesis's Table A2.1; nothing is fitted.
  //  Compared with
  //    * the measured markers, re-digitised from the 300 dpi figures
  //      (core/data/nguyen_2022/Fig{4,5}/): outlet y_CO2 (NRMSE / feed y_CO2 and
  //      the front time at half the feed fraction), the thermocouple 1.82 cm
  //      from the outlet (RMSE) and, in Fig. 5, the outlet y_N2 (RMSE);
  //    * the paper's own model curves (core/data/nguyen_2022/paper_model/),
  //      code to code: the y_CO2 NRMSE and front time ("vs paper model").
  //  README "Validation" has the results, the attribution and the open items.
  // ==========================================================================
  namespace {

    constexpr double kNotUsed = std::numeric_limits<double>::infinity();

    struct NguyenTolerances {
      double nrmse_co2, t50, t_rmse, rmse_n2, model_nrmse, model_t50;
    };

    struct NguyenEntry {
      NguyenCase c;
      NguyenTolerances tol;
      std::string description;
    };

    const std::vector<NguyenEntry>& NguyenCases() {
      static const std::vector<NguyenEntry> cases = {
        // Tolerances: measured error x 1.1, rounded up (front times at least 0.002;
        // README "Validation"). The measured errors are in the comment above each case.
        // measured: NRMSE y_CO2, front, T RMSE [K], y_N2 RMSE, NRMSE vs model, front vs model = 0.2944, 0.0652, 7.965, -, 0.1533, 0.01787
        { { "Fig4_100pct", "Fig4", 100, 50, 22.4, 0.97, 7100.0, "Nguyen2022_Fig4_100pct" },
          { 0.33, 0.075, 9, kNotUsed, 0.17, 0.02 },
          "pure CO2 (balance He)" },
        // measured: NRMSE y_CO2, front, T RMSE [K], y_N2 RMSE, NRMSE vs model, front vs model = 0.1179, 0.03247, 5.232, -, 0.08878, 0.01651
        { { "Fig4_70pct", "Fig4", 70, 50, 22.3, 0.97, 7100.0, "Nguyen2022_Fig4" },
          { 0.13, 0.04, 6, kNotUsed, 0.1, 0.02 },
          "70 % CO2 in He" },
        // measured: NRMSE y_CO2, front, T RMSE [K], y_N2 RMSE, NRMSE vs model, front vs model = 0.1875, 0.1086, 6.28, -, 0.1069, 0.03891
        { { "Fig4_50pct", "Fig4", 50, 50, 22.4, 0.97, 7100.0, "Nguyen2022_Fig4_50pct" },
          { 0.21, 0.12, 7, kNotUsed, 0.12, 0.045 },
          "50 % CO2 in He" },
        // measured: NRMSE y_CO2, front, T RMSE [K], y_N2 RMSE, NRMSE vs model, front vs model = 0.1518, 0.07249, 2.659, 0.2822, 0.09945, 0.001455
        { { "Fig5_30pct", "Fig5", 30, 100, 21.8, 0.98, 11800.0, "Nguyen2022_Fig5_30pct" },
          { 0.17, 0.08, 3, 0.32, 0.11, 0.005 },
          "30 % CO2 in N2" },
        // measured: NRMSE y_CO2, front, T RMSE [K], y_N2 RMSE, NRMSE vs model, front vs model = 0.1336, 0.05556, 1.028, 0.1387, 0.08495, 0.04006
        { { "Fig5_15pct", "Fig5", 15, 100, 20.3, 0.97, 11800.0, "Nguyen2022_Fig5" },
          { 0.15, 0.065, 1.25, 0.16, 0.1, 0.045 },
          "15 % CO2 in N2" },
        // measured: NRMSE y_CO2, front, T RMSE [K], y_N2 RMSE, NRMSE vs model, front vs model = 0.09487, 0.09655, 1.721, 0.03553, 0.1299, 0.04139
        { { "Fig5_5pct", "Fig5", 5, 100, 23.4, 0.97, 11800.0, "Nguyen2022_Fig5_5pct" },
          { 0.11, 0.11, 2, 0.04, 0.15, 0.05 },
          "5 % CO2 in N2" },
      };
      return cases;
    }

    const NguyenEntry& NguyenEntryFor(const std::string& name) {
      for (const auto& e : NguyenCases())
        if (e.c.name == name) return e;
      throw std::invalid_argument("no Nguyen 2022 case " + name);
    }

    std::string Pct(int p) { return std::to_string(p) + "%"; }
  }

  class Nguyen2022Test : public Nguyen2022Fixture {
  public:
    void Run(const std::string& name) {
      const NguyenEntry& e = NguyenEntryFor(name);
      const NguyenCase& c = e.c;
      const std::string fig_label = c.fig == "Fig4" ? "Fig. 4" : "Fig. 5";
      data.info.title = "Nguyen 2022 - " + fig_label + ", " + e.description;
      data.info.description =
        "CO2 breakthrough on CALF-20 (Nguyen et al. 2022, " + fig_label + "): " + e.description + " fed at " +
        std::to_string(static_cast<int>(c.flow_ccm)) + " cm3/min into a column full of helium. The measured outlet "
        "composition and the thermocouple 1.82 cm from the outlet, re-digitised from the paper's figure, are bound to "
        "the deck (Experimental Data).";
      data.info.author = "Sam Affleck";
      SetUpNguyen2022(c);

      // Every measured curve of the figure is embedded; the case's own curves are
      // bound to what they measure.
      SaveExperimentalData("data/nguyen_2022/" + c.fig);
      for (auto& xy : data.exp_data_.GetAllData()) {
        ObservableBinding b;
        if (xy.name == "Y_CO2_" + Pct(c.percent)) {
          b.kind = ObservableKind::Component; b.variable = "mole_fraction"; b.component = "CarbonDioxide";
        } else if (xy.name == "Y_N2_" + Pct(c.percent)) {
          b.kind = ObservableKind::Component; b.variable = "mole_fraction"; b.component = "Nitrogen";
        } else if (xy.name == "T_CO2_" + Pct(c.percent)) {
          b.kind = ObservableKind::Bulk; b.variable = "temperature";
          b.location = ObservableLocation::Position; b.z = ThermocoupleZ();
        } else {
          continue;  // another feed of the same figure
        }
        xy.observable = b;
      }

      CoupledSolver solver{data};
      std::string error{};
      EXPECT_TRUE(data.PerformChecks(error)) << "System checks failed: " << error;
      ASSERT_GE(solver.IntegrateCycle(data.cycle, false), 0) << "Solver did not complete";
      EXPECT_TRUE(IsMassConserved(solver.GetLog())) << "Mass is not conserved.";
      EXPECT_TRUE(IsMolesConserved(solver.GetLog())) << "Moles are not conserved.";
      EXPECT_TRUE(IsComponentMolesConserved(solver.GetLog())) << "Component moles not conserved.";
      SaveResultsAsExample(c.deck, solver);

      using namespace validation;
      CaseReport report("nguyen_2022", c.name);
      const double y = c.percent / 100.0;
      const std::string dir = "data/nguyen_2022/" + c.fig + "/";
      const std::string model_dir = "data/nguyen_2022/paper_model/" + c.fig + "/";
      const std::string src = fig_label + " (measured markers, re-digitised)";
      const std::string model_src = fig_label + " (paper's model curve, digitised)";
      const Series co2 = FaceMoleFraction(solver, "CarbonDioxide", 0, 0);

      report.CheckCurve("outlet y_CO2", co2, LoadSeriesCsv(dir + "Y_CO2_" + Pct(c.percent) + ".csv"), e.tol.nrmse_co2, "-", src, y);
      report.CheckCrossingTime("CO2 front time (y_CO2 = y_feed/2)", co2, LoadSeriesCsv(dir + "Y_CO2_" + Pct(c.percent) + ".csv"),
                               0.5 * y, e.tol.t50, src);
      const Series T = CellVariableAt(solver, VariableType::TEMPERATURE, ThermocoupleZ(), BedLength(), 30);
      report.CheckCurve("thermocouple T (1.82 cm from the outlet)", T, LoadSeriesCsv(dir + "T_CO2_" + Pct(c.percent) + ".csv"),
                        e.tol.t_rmse, "K", src);
      if (c.fig == "Fig5")
        report.CheckCurve("outlet y_N2", FaceMoleFraction(solver, "Nitrogen", 0, 0),
                          LoadSeriesCsv(dir + "Y_N2_" + Pct(c.percent) + ".csv"), e.tol.rmse_n2, "-", src);

      report.CheckCurve("outlet y_CO2 vs paper model", co2, LoadSeriesCsv(model_dir + "Y_CO2_" + Pct(c.percent) + ".csv"),
                        e.tol.model_nrmse, "-", model_src, y);
      report.CheckCrossingTime("CO2 front time (y_CO2 = y_feed/2) vs paper model", co2,
                               LoadSeriesCsv(model_dir + "Y_CO2_" + Pct(c.percent) + ".csv"), 0.5 * y, e.tol.model_t50, model_src);
    }
  };

  TEST_F(Nguyen2022Test, Fig4_100pct) { Run("Fig4_100pct"); }
  TEST_F(Nguyen2022Test, Fig4_70pct) { Run("Fig4_70pct"); }
  TEST_F(Nguyen2022Test, Fig4_50pct) { Run("Fig4_50pct"); }
  TEST_F(Nguyen2022Test, Fig5_30pct) { Run("Fig5_30pct"); }
  TEST_F(Nguyen2022Test, Fig5_15pct) { Run("Fig5_15pct"); }
  TEST_F(Nguyen2022Test, Fig5_5pct) { Run("Fig5_5pct"); }

} // namespace adsim
