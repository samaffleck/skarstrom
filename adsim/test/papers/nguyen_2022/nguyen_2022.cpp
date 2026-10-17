// STL includes
#define _USE_MATH_DEFINES
#include <cmath>
#include <filesystem>

// AdSim includes
#include "../../../src/pch.h"
#include "adsim/solver/optimisation.h"
#include "adsim/solver/coupled_solver.h"
#include "adsim/solver/kpis.h"
#include "adsim/math/constants.h"
#include "adsim/equations/isotherm.h"
#include "adsim/boundary_conditions/boundary_conditions.h"
#include "adsim/boundary_conditions/boundary_stream.h"
#include "adsim/boundary_conditions/connection.h"
#include "adsim/boundary_conditions/step.h"
#include "adsim/equations/dispersion_model.h"
#include "adsim/objects/compressor.h"
#include "adsim/equations/heat_transfer_model.h"

// GTest includes
#include <gtest/gtest.h>

#include "nguyen_2022.h"
#include "../../fixtures/validation_metrics.h"
#include <limits>

namespace adsim {

  // ==========================================================================
  //  Nguyen, Lin, Shimizu & Rajendran (2022), "Separation of CO2 and N2 on a
  //  hydrophobic metal organic framework CALF-20", Chem. Eng. J. 442, 136263,
  //  doi:10.1016/j.cej.2022.136263.
  //  Quantitative check: outlet CO2 mole fraction vs the digitised breakthrough
  //  curves in adsim/data/nguyen_2022/Fig{4,5}/Y_CO2_<feed%>.csv, as a
  //  normalised RMSE (RMSE / feed y_CO2). The thermocouple position of the
  //  digitised T_CO2_*.csv curves is not stated in the material available to us
  //  (the paper text is not in the repo), so temperatures are not scored.
  // ==========================================================================
  namespace {

    void ValidateNguyen2022(const CoupledSolver& solver, const std::string& case_name,
                            const std::string& fig, int feed_percent, double tol_nrmse,
                            double tol_t50) {
      using namespace validation;
      CaseReport report("nguyen_2022", case_name);
      const std::string file = "data/nguyen_2022/" + fig + "/Y_CO2_" + std::to_string(feed_percent) + "%.csv";
      report.CheckCurve("outlet y_CO2", FaceMoleFraction(solver, "CarbonDioxide", 0, 0), LoadSeriesCsv(file),
                        tol_nrmse, "-", (fig == "Fig4" ? "Fig. 4" : "Fig. 5") + std::string(" (digitised)"),
                        feed_percent / 100.0);
      report.CheckCrossingTime("CO2 front time (y_CO2 = y_feed/2)", FaceMoleFraction(solver, "CarbonDioxide", 0, 0),
                               LoadSeriesCsv(file), 0.5 * feed_percent / 100.0, tol_t50,
                               (fig == "Fig4" ? "Fig. 4" : "Fig. 5") + std::string(" (digitised)"));
    }
  }

  // Additional feed compositions of Figs. 4 and 5 (same column, same transport
  // parameters, only the inlet CO2 fraction and flow change). Not saved as GUI
  // examples; the two bundled examples are the tests further below.
  void RunNguyen2022Extra(Nguyen2022Test& t, const std::string& case_name, const std::string& fig,
                          int feed_percent, double flow_ccm, double duration, double tol_nrmse,
                          double tol_t50) {
    const std::string CO2 = "CarbonDioxide";
    const std::string N2 = "Nitrogen";
    const std::string He = "Helium";
    t.data.info.title = "Nguyen 2022 - " + case_name;
    t.SetUpNguyen2022();
    t.AddInitialConditions(0.97e5, 295, 295, { { CO2, 0.0 }, { N2, 0.0 }, { He, 1.0 } });
    t.AddStep("Feed", duration);
    const double y = feed_percent / 100.0;
    std::vector<std::unique_ptr<IBoundaryStream>> streams;
    streams.push_back(t.CreateNoValveInletStream(t.CreateVolumeFlowBC(flow_ccm * 1.666666666666667e-8),
                                                 { { CO2, y }, { N2, 0.0 }, { He, 1.0 - y } }, 295.0));
    t.AddFlowBC(std::move(streams), BoundaryLocation::BOTTOM, BoundaryDirection::INLET, "Feed");
    t.AddPressureOutletBC(t.CreateConstantBoundaryValue(0.97e5), BoundaryLocation::TOP, "Feed",
                          { { CO2, 0.0 }, { N2, 0.0 }, { He, 1.0 } });
    t.ConfigureCycle(1);

    CoupledSolver solver{t.data};
    std::string error{};
    EXPECT_TRUE(t.data.PerformChecks(error)) << "System checks failed: " << error;
    ASSERT_GE(solver.IntegrateCycle(t.data.cycle, false), 0) << "Solver did not complete";
    EXPECT_TRUE(t.IsMassConserved(solver.GetLog()));
    EXPECT_TRUE(t.IsMolesConserved(solver.GetLog()));
    EXPECT_TRUE(t.IsComponentMolesConserved(solver.GetLog()));
    t.SaveTest("Nguyen_2022_" + case_name, solver);
    ValidateNguyen2022(solver, case_name, fig, feed_percent, tol_nrmse, tol_t50);
  }

  TEST_F(Nguyen2022Test, Fig5_TernaryBreakthrough) {
    data.info.title = "Nguyen 2022 - Fig. 5";
    data.info.description = "Ternary CO2/N2/He breakthrough experiment with CALF-20. Results correspond Fig 5 from the publication. The example is initialised with an inlet CO2 mole fraction of 0.15. Varying this value will correspond to the other cases in the paper. Experimental data from the publication is included in the example under the Data Plotter window.";
    data.info.author = "Sam Affleck";
    
    const std::string CO2 = "CarbonDioxide";
    const std::string N2 = "Nitrogen";
    const std::string He = "Helium";
    const std::string calf20 = "CALF-20";

    SetUpNguyen2022();

    AddInitialConditions(
      0.97e5, 295, 295,
      {
        { CO2, 0.0 },
        { N2, 0.0 },
        { He, 1.0 }
      }
    );

    AddStep("Feed", 12000.0);

    const double volumetric_flowrate_ccm = 100;
    const double volumetric_flowrate_m3_s = volumetric_flowrate_ccm * 1.666666666666667e-8;

    Step& step_feed = data.cycle.GetStep("Feed");
    auto feed_inlet_stream = CreateNoValveInletStream(
      CreateVolumeFlowBC(volumetric_flowrate_m3_s),
      {
        { CO2, 0.15 },
        { N2, 0.0 },
        { He, 0.85 }
      },
      295.0
    );
    feed_inlet_stream->name_ = "CO2 Inlet";
    
    auto feed_flow_bc = std::make_unique<FlowBoundaryCondition>();
    feed_flow_bc->direction_ = BoundaryDirection::INLET;
    feed_flow_bc->location_ = BoundaryLocation::BOTTOM;
    feed_flow_bc->streams_.push_back(std::move(feed_inlet_stream));
    step_feed.bottom_boundary_condition_ = std::move(feed_flow_bc);

    auto feed_pressure_out = std::make_unique<PressureBoundaryCondition>();
    feed_pressure_out->direction_ = BoundaryDirection::OUTLET;
    feed_pressure_out->location_ = BoundaryLocation::TOP;
    feed_pressure_out->pressure_ = CreateConstantBoundaryValue(0.97e5);
    feed_pressure_out->mole_fraction_[CO2] = 0.0;
    feed_pressure_out->mole_fraction_[N2] = 0.0;
    feed_pressure_out->mole_fraction_[He] = 1.0;
    step_feed.top_boundary_condition_ = std::move(feed_pressure_out);

    ConfigureCycle(1);

    SaveExperimentalData("data/nguyen_2022/Fig5");

    CoupledSolver solver{data};

    std::string error{};
    EXPECT_TRUE(data.PerformChecks(error)) << "System checks failed: " << error;
    
    bool verbose = false;
    int retval = solver.IntegrateCycle(data.cycle, verbose);
    EXPECT_GE(retval, 0) << "Solver did not complete";
    
    EXPECT_TRUE(IsMassConserved(solver.GetLog())) << "Mass is not conserved.";
    EXPECT_TRUE(IsMolesConserved(solver.GetLog())) << "Mass is not conserved.";
    EXPECT_TRUE(IsComponentMolesConserved(solver.GetLog())) << "Component moles not conserved.";

    SaveTest("Nguyen_2022_Fig5", solver);
    SaveResultsAsExample("Nguyen_2022_Fig5", solver);

    ValidateNguyen2022(solver, "Fig5_TernaryBreakthrough", "Fig5", 15, 0.097, 0.059);  // measured NRMSE 0.0842158, front-time error 0.0488449
  }

  TEST_F(Nguyen2022Test, Fig4_BinaryBreakthrough) {
    data.info.title = "Nguyen 2022 - Fig. 4";
    data.info.description = "Binary CO2/He breakthrough experiment with CALF-20. Results correspond Fig 4 from the publication. The example is initialised with an inlet CO2 mole fraction of 0.7. Varying this value will correspond to the other cases in the paper. Experimental data from the publication is included in the example under the Data Plotter window.";
    data.info.author = "Sam Affleck";
    
    const std::string CO2 = "CarbonDioxide";
    const std::string N2 = "Nitrogen";
    const std::string He = "Helium";
    const std::string calf20 = "CALF-20";

    SetUpNguyen2022();

    AddInitialConditions(
      0.97e5, 295, 295,
      {
        { CO2, 0.0 },
        { N2, 0.0 },
        { He, 1.0 }
      }
    );

    AddStep("Feed", 6000.0);

    const double volumetric_flowrate_ccm = 50;
    const double volumetric_flowrate_m3_s = volumetric_flowrate_ccm * 1.666666666666667e-8;

    Step& step_feed = data.cycle.GetStep("Feed");
    auto feed_inlet_stream = CreateNoValveInletStream(
      CreateVolumeFlowBC(volumetric_flowrate_m3_s),
      {
        { CO2, 0.7 },
        { N2, 0.0 },
        { He, 0.3 }
      },
      295.0
    );
    feed_inlet_stream->name_ = "CO2 Inlet";
    
    auto feed_flow_bc = std::make_unique<FlowBoundaryCondition>();
    feed_flow_bc->direction_ = BoundaryDirection::INLET;
    feed_flow_bc->location_ = BoundaryLocation::BOTTOM;
    feed_flow_bc->streams_.push_back(std::move(feed_inlet_stream));
    step_feed.bottom_boundary_condition_ = std::move(feed_flow_bc);

    auto feed_pressure_out = std::make_unique<PressureBoundaryCondition>();
    feed_pressure_out->direction_ = BoundaryDirection::OUTLET;
    feed_pressure_out->location_ = BoundaryLocation::TOP;
    feed_pressure_out->pressure_ = CreateConstantBoundaryValue(0.97e5);
    feed_pressure_out->mole_fraction_[CO2] = 0.0;
    feed_pressure_out->mole_fraction_[N2] = 0.0;
    feed_pressure_out->mole_fraction_[He] = 1.0;
    step_feed.top_boundary_condition_ = std::move(feed_pressure_out);

    ConfigureCycle(1);

    SaveExperimentalData("data/nguyen_2022/Fig4");

    CoupledSolver solver{data};

    std::string error{};
    EXPECT_TRUE(data.PerformChecks(error)) << "System checks failed: " << error;
    
    bool verbose = false;
    int retval = solver.IntegrateCycle(data.cycle, verbose);
    EXPECT_GE(retval, 0) << "Solver did not complete";
    
    EXPECT_TRUE(IsMassConserved(solver.GetLog())) << "Mass is not conserved.";
    EXPECT_TRUE(IsMolesConserved(solver.GetLog())) << "Mass is not conserved.";
    EXPECT_TRUE(IsComponentMolesConserved(solver.GetLog())) << "Component moles not conserved.";

    SaveTest("Nguyen_2022_Fig4", solver);
    SaveResultsAsExample("Nguyen_2022_Fig4", solver);

    ValidateNguyen2022(solver, "Fig4_BinaryBreakthrough", "Fig4", 70, 0.1, 0.027);  // measured NRMSE 0.0861904, front-time error 0.0169452
  }

  // Fig. 4 feeds use 50 ccm (as the 70 % case above); Fig. 5 feeds use 100 ccm
  // (as the 15 % case). Durations cover the digitised record.
  // measured NRMSE 0.121728, front-time error 0.0497634
  TEST_F(Nguyen2022Test, Fig4_Binary_50pct) { RunNguyen2022Extra(*this, "Fig4_Binary_50pct", "Fig4", 50, 50, 7100.0, 0.14, 0.06); }
  // measured NRMSE 0.0993589, front-time error 0.00687957
  TEST_F(Nguyen2022Test, Fig4_Binary_100pct) { RunNguyen2022Extra(*this, "Fig4_Binary_100pct", "Fig4", 100, 50, 7100.0, 0.12, 0.017); }
  // measured NRMSE 0.124475, front-time error 0.0155786
  TEST_F(Nguyen2022Test, Fig5_Ternary_5pct) { RunNguyen2022Extra(*this, "Fig5_Ternary_5pct", "Fig5", 5, 100, 11800.0, 0.15, 0.026); }
  // measured NRMSE 0.104462, front-time error 0.0282569
  TEST_F(Nguyen2022Test, Fig5_Ternary_30pct) { RunNguyen2022Extra(*this, "Fig5_Ternary_30pct", "Fig5", 30, 100, 11800.0, 0.13, 0.039); }

} // namespace adsim
