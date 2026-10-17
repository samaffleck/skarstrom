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

// Paper fixture (Ward et al. 2022) — see README.md in this folder.
#include "../../fixtures/user_data_test.h"
#include "../../fixtures/validation_metrics.h"
#include <limits>

namespace adsim {

  namespace {
    // Validation tolerances = measured error rounded up.
    // Q50ccm: measured NRMSE 0.381727, T RMSE 2.58832 K, front-time error 0.0427429
    constexpr double kTolY_50 = 0.44, kTolT_50 = 3.1, kTolT50_50 = 0.053;
    // Q100ccm: measured NRMSE 0.379455, T RMSE 3.54418 K, front-time error 0.0611066
    constexpr double kTolY_100 = 0.44, kTolT_100 = 4.1, kTolT50_100 = 0.072;
    // Q200ccm: measured NRMSE 0.335826, T RMSE 4.20647 K, front-time error 0.0458177
    constexpr double kTolY_200 = 0.39, kTolT_200 = 4.9, kTolT50_200 = 0.056;
  }

  class Ward2022Test : public UserDataTest {
  public:
    // Quantitative comparison with the digitised experiments of Ward & Pini
    // (2022), Adsorption 28:161-183, doi:10.1007/s10450-022-00361-z (data of
    // their ref. [5]): outlet CO2 mole fraction (feed is pure CO2, so the
    // RMSE is already normalised) and the thermocouple at z = 0.0715 m. The
    // paper's 95 % credibility intervals for the model-predicted breakthrough
    // time and temperature rise (Table 10) are recorded for context only - the
    // indicator definitions (Table 9) did not survive the text extraction.
    void ValidateWard(const CoupledSolver& solver, const std::string& case_name, int flow_ccm,
                      double tol_y_rmse, double tol_T_rmse,
                      double tol_t50) {
      using namespace validation;
      CaseReport report("ward_2022", case_name);
      const std::string q = std::to_string(flow_ccm);
      const std::string dir = "data/ward_2022/Qin=" + q + "/";
      const Series y_sim = FaceMoleFraction(solver, "CarbonDioxide", 0, 0);
      report.CheckCurve("outlet y_CO2", y_sim,
                        LoadSeriesCsv(dir + "outlet_CO2_mole_fraction[-]_Qin=" + q + "[ccm].csv"),
                        tol_y_rmse, "-", "Fig. 3 / data of ref. [5] (digitised)", 1.0);
      report.CheckCrossingTime("CO2 front time (y_CO2 = 0.5)", y_sim,
                               LoadSeriesCsv(dir + "outlet_CO2_mole_fraction[-]_Qin=" + q + "[ccm].csv"),
                               0.5, tol_t50, "Fig. 3 / data of ref. [5] (digitised)");
      const Series T_sim = CellVariableAt(solver, VariableType::TEMPERATURE, 0.0715, 0.286, 60);
      report.CheckCurve("thermocouple T at z = 0.0715 m", T_sim,
                        LoadSeriesCsv(dir + "temperature[K]_z=0.0715[m]_Qin=" + q + "[ccm].csv"),
                        tol_T_rmse, "K", "Fig. 3 / data of ref. [5] (digitised)");

      auto ci_note = [](double v, double lo, double hi) {
        return std::string(v >= lo && v <= hi ? "inside" : "outside") + " the paper's 95% CI [" +
               JsonNumber(lo) + ", " + JsonNumber(hi) + "]";
      };
      const RefValue t_lo = Reference("ward_2022", case_name, "breakthrough_time_CI_low");
      const RefValue t_hi = Reference("ward_2022", case_name, "breakthrough_time_CI_high");
      RefValue t_mid = t_lo; t_mid.value = 0.5 * (t_lo.value + t_hi.value);
      const double t_br = CrossingTime(y_sim, 0.5);
      report.Info("breakthrough time (y_CO2,out = 0.5)", t_br, t_mid, "abs_error",
                  "reference = CI midpoint; simulated value " + ci_note(t_br, t_lo.value, t_hi.value) +
                  "; t_br read as the time y_CO2,out reaches 0.5 (Table 9 only partly legible)");
      const RefValue d_lo = Reference("ward_2022", case_name, "temperature_rise_CI_low");
      const RefValue d_hi = Reference("ward_2022", case_name, "temperature_rise_CI_high");
      RefValue d_mid = d_lo; d_mid.value = 0.5 * (d_lo.value + d_hi.value);
      const double dT = T_sim.empty() ? 0.0 : MaxValue(T_sim) - T_sim.y.front();
      report.Info("temperature rise at z = 0.0715 m", dT, d_mid, "abs_error",
                  "reference = CI midpoint; simulated value " + ci_note(dT, d_lo.value, d_hi.value) +
                  "; simulated = max T at the thermocouple minus initial T");
    }
  };

  TEST_F(Ward2022Test, Q100ccm) {
    data.info.title = "Ward 2022 - 100 ccm Flow Rate";
    data.info.description = "CO2 breakthrough curve test based on Ward et al. 2022 for activated carbon in a lab-scale test rig. "
                            "This is a binary CO2/Air breakthrough simulation where the bed is initially saturated with air and then fed with pure CO2. "
                            "The adsorbent uses a dual-site Langmuir isotherm model and Haynes linear driving force model for mass transfer. "
                            "The simulation includes coupled heat and mass transfer effects with wall heat transfer. "
                            "Experimental breakthrough data is included and can be plotted to validate the solver performance. "
                            "This configuration tests the medium flow rate condition (100 ccm/min) of three flow rates studied.";
    data.info.author = "Sam Affleck";
    
    AddComponents({ "CarbonDioxide", "Air" });

    const double Di = 0.0189;

    auto [layer_id, ads_id] = AddLayer("Activated Carbon");
    auto& AC = GetLayer(layer_id);
    AC.domain_.length = 0.286;
    std::get<CylinderBed>(data.bed.shape_).diameter = Di;
    AC.particle_diameter_ = 0.003;
    AC.eb_ = 0.456;
    AC.domain_.n_cells = 60;

    auto& AC_ads = GetAdsorbent(ads_id);
    AC_ads.solid_density_ = 799;
    AC_ads.heat_capacity_ = 1040;
    AC_ads.thermal_conductivity_ = 0.63;
    AC_ads.ep_ = 0.0; // 0.59

    data.bed.wall_.heat_capacity_ = 902;
    data.bed.wall_.density_ = 2700;
    data.bed.wall_.thermal_conductivity_ = 205;
    data.bed.wall_.thickness_ = 3.25e-3;
    data.bed.ambient_temperature_ = 291.4; // Table 1: T_a = 291.4 K (was 293 = T_F, a transcription slip)
    data.fluid.adsorbed_heat_capacity_["CarbonDioxide"] = 36.96;

    auto dsl_co2 = std::make_unique<DualSiteLangmuirIsotherm_Concentration>("CarbonDioxide");
    dsl_co2->params_["CarbonDioxide"].m1 = 0.418;
    dsl_co2->params_["CarbonDioxide"].b1 = 2.05e-6;
    dsl_co2->params_["CarbonDioxide"].q1 = 29980;
    dsl_co2->params_["CarbonDioxide"].m2 = 6.04;
    dsl_co2->params_["CarbonDioxide"].b2 = 3.32e-6;
    dsl_co2->params_["CarbonDioxide"].q2 = 20367;
    AddIsotherm({ "CarbonDioxide" }, ads_id, std::move(dsl_co2));

    AddLinearDrivingForceModel_Haynes(ads_id, "CarbonDioxide", 3, 0.59);
    AddLinearDrivingForceModel(ads_id, "Air", 0.0);

    SetHOA(ads_id, "CarbonDioxide", 29980.0);
    SetHOA(ads_id, "Air", 0.0);

    auto& layer = GetLayer(layer_id);
    layer.equations_.dispersion_model_ = CreateAxialDispersionModel(AxialDispersionType::RUTHVEN);
    layer.equations_.heat_dispersion_model_ = CreateAxialHeatDispersionModel(AxialHeatDispersionType::YAGI_AND_KUNII);

    AddConstantInternalHeatTransferCoefficientModel(57.9, layer_id);
    auto ihtc = std::make_unique<ExternalHeatTransferModel_Holman>();
    ihtc->c = 0.027;
    ihtc->m = 0.8;
    ihtc->n = 1.0/3.0;
    ihtc->external_velocity = 0.5;
    layer.equations_.external_heat_transfer_model_ = std::move(ihtc);

    AddInitialConditions(
      1.0e5, 291.4, 291.4,
      {
        {"CarbonDioxide", 0.0 },
        {"Air", 1.0 }
      }
    );

    AddStep("Feed", 3000.0);

    const double volumetric_flowrate_ccm = 100;
    const double volumetric_flowrate_m3_s = volumetric_flowrate_ccm * 1.666666666666667e-8;

    Step& step_feed = data.cycle.GetStep("Feed");
    auto feed_inlet_stream = CreateNoValveInletStream(
      CreateVolumeFlowBC(volumetric_flowrate_m3_s),
      {
        {"CarbonDioxide", 1.0 },
        {"Air", 0.0 }
      },
      293.0
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
    feed_pressure_out->pressure_ = CreateConstantBoundaryValue(1e5);
    feed_pressure_out->mole_fraction_["CarbonDioxide"] = 0.0;
    feed_pressure_out->mole_fraction_["Air"] = 1.0;
    step_feed.top_boundary_condition_ = std::move(feed_pressure_out);

    ConfigureCycle(1);

    SaveExperimentalData("data/ward_2022/Qin=100");

    CoupledSolver solver{data};

    std::string error{};
    EXPECT_TRUE(data.PerformChecks(error)) << "System checks failed: " << error;
    
    bool verbose = false;
    int retval = solver.IntegrateCycle(data.cycle, verbose);
    EXPECT_GE(retval, 0) << "Solver did not complete";
    
    EXPECT_TRUE(IsMassConserved(solver.GetLog())) << "Mass is not conserved.";
    EXPECT_TRUE(IsMolesConserved(solver.GetLog())) << "Mass is not conserved.";
    EXPECT_TRUE(IsComponentMolesConserved(solver.GetLog())) << "Component moles not conserved.";


    SaveTest("Ward_2022_Q=100[ccm]", solver);
    SaveResultsAsExample("Ward_2022_Q=100[ccm]", solver);

    ValidateWard(solver, "Q100ccm", 100, kTolY_100, kTolT_100, kTolT50_100);
  }

  TEST_F(Ward2022Test, Q50ccm) {
    data.info.title = "Ward 2022 - 50 ccm Flow Rate";
    data.info.description = "CO2 breakthrough curve test based on Ward et al. 2022 for activated carbon in a lab-scale test rig. "
                            "This is a binary CO2/Air breakthrough simulation where the bed is initially saturated with air and then fed with pure CO2. "
                            "The adsorbent uses a dual-site Langmuir isotherm model and Haynes linear driving force model for mass transfer. "
                            "The simulation includes coupled heat and mass transfer effects with wall heat transfer. "
                            "Experimental breakthrough data is included and can be plotted to validate the solver performance. "
                            "This configuration tests the low flow rate condition (50 ccm/min) of three flow rates studied.";
    data.info.author = "Sam Affleck";
    
    AddComponents({ "CarbonDioxide", "Air" });

    const double Di = 0.0189;

    auto [layer_id, ads_id] = AddLayer("Activated Carbon");
    auto& AC = GetLayer(layer_id);
    AC.domain_.length = 0.286;
    std::get<CylinderBed>(data.bed.shape_).diameter = Di;
    AC.particle_diameter_ = 0.003;
    AC.eb_ = 0.456;
    AC.domain_.n_cells = 60;

    auto& AC_ads = GetAdsorbent(ads_id);
    AC_ads.solid_density_ = 799;
    AC_ads.heat_capacity_ = 1040;
    AC_ads.thermal_conductivity_ = 0.63;
    AC_ads.ep_ = 0.0;

    data.bed.wall_.heat_capacity_ = 902;
    data.bed.wall_.density_ = 2700;
    data.bed.wall_.thermal_conductivity_ = 205;
    data.bed.wall_.thickness_ = 3.25e-3;
    data.bed.ambient_temperature_ = 291.4; // Table 1: T_a = 291.4 K (was 293 = T_F, a transcription slip)
    data.fluid.adsorbed_heat_capacity_["CarbonDioxide"] = 36.96;

    auto dsl_co2 = std::make_unique<DualSiteLangmuirIsotherm_Concentration>("CarbonDioxide");
    dsl_co2->params_["CarbonDioxide"].m1 = 0.418;
    dsl_co2->params_["CarbonDioxide"].b1 = 2.05e-6;
    dsl_co2->params_["CarbonDioxide"].q1 = 29980;
    dsl_co2->params_["CarbonDioxide"].m2 = 6.04;
    dsl_co2->params_["CarbonDioxide"].b2 = 3.32e-6;
    dsl_co2->params_["CarbonDioxide"].q2 = 20367;
    AddIsotherm({ "CarbonDioxide" }, ads_id, std::move(dsl_co2));

    AddLinearDrivingForceModel_Haynes(ads_id, "CarbonDioxide", 3, 0.59);
    AddLinearDrivingForceModel(ads_id, "Air", 0.0);

    SetHOA(ads_id, "CarbonDioxide", 29980.0);
    SetHOA(ads_id, "Air", 0.0);

    auto& layer = GetLayer(layer_id);
    layer.equations_.dispersion_model_ = CreateAxialDispersionModel(AxialDispersionType::RUTHVEN);
    layer.equations_.heat_dispersion_model_ = CreateAxialHeatDispersionModel(AxialHeatDispersionType::YAGI_AND_KUNII);

    AddConstantInternalHeatTransferCoefficientModel(77.5, layer_id);
    auto ihtc = std::make_unique<ExternalHeatTransferModel_Holman>();
    ihtc->c = 0.027;
    ihtc->m = 0.8;
    ihtc->n = 1.0/3.0;
    ihtc->external_velocity = 0.5;
    layer.equations_.external_heat_transfer_model_ = std::move(ihtc);

    AddInitialConditions(
      1.0e5, 291.4, 291.4,
      {
        {"CarbonDioxide", 0.0 },
        {"Air", 1.0 }
      }
    );

    AddStep("Feed", 3000.0);

    const double volumetric_flowrate_ccm = 50;
    const double volumetric_flowrate_m3_s = volumetric_flowrate_ccm * 1.666666666666667e-8;
    
    Step& step_feed = data.cycle.GetStep("Feed");
    auto feed_inlet_stream = CreateNoValveInletStream(
      CreateVolumeFlowBC(volumetric_flowrate_m3_s),
      {
        {"CarbonDioxide", 1.0 },
        {"Air", 0.0 }
      },
      293.0
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
    feed_pressure_out->pressure_ = CreateConstantBoundaryValue(1e5);
    feed_pressure_out->mole_fraction_["CarbonDioxide"] = 0.0;
    feed_pressure_out->mole_fraction_["Air"] = 1.0;
    step_feed.top_boundary_condition_ = std::move(feed_pressure_out);

    ConfigureCycle(1);

    SaveExperimentalData("data/ward_2022/Qin=50");

    CoupledSolver solver{data};

    std::string error{};
    EXPECT_TRUE(data.PerformChecks(error)) << "System checks failed: " << error;
    
    bool verbose = false;
    int retval = solver.IntegrateCycle(data.cycle, verbose);
    EXPECT_GE(retval, 0) << "Solver did not complete";

    EXPECT_TRUE(IsMassConserved(solver.GetLog())) << "Mass is not conserved.";
    EXPECT_TRUE(IsMolesConserved(solver.GetLog())) << "Mass is not conserved.";
    EXPECT_TRUE(IsComponentMolesConserved(solver.GetLog())) << "Component moles not conserved.";

    // solver.PrintStatistics();
    SaveTest("Ward_2022_Q=50[ccm]", solver);
    SaveResultsAsExample("Ward_2022_Q=50[ccm]", solver);

    ValidateWard(solver, "Q50ccm", 50, kTolY_50, kTolT_50, kTolT50_50);
  }

  TEST_F(Ward2022Test, Q200ccm) {
    data.info.title = "Ward 2022 - 200 ccm Flow Rate";
    data.info.description = "CO2 breakthrough curve test based on Ward et al. 2022 for activated carbon in a lab-scale test rig. "
                            "This is a binary CO2/Air breakthrough simulation where the bed is initially saturated with air and then fed with pure CO2. "
                            "The adsorbent uses a dual-site Langmuir isotherm model and Haynes linear driving force model for mass transfer. "
                            "The simulation includes coupled heat and mass transfer effects with wall heat transfer. "
                            "Experimental breakthrough data is included and can be plotted to validate the solver performance. "
                            "This configuration tests the high flow rate condition (200 ccm/min) of three flow rates studied.";
    data.info.author = "Sam Affleck";
    
    // Start with Q100 setup
    AddComponents({ "CarbonDioxide", "Air" });

    const double Di = 0.0189;

    auto [layer_id, ads_id] = AddLayer("Activated Carbon");
    auto& AC = GetLayer(layer_id);
    AC.domain_.length = 0.286;
    std::get<CylinderBed>(data.bed.shape_).diameter = Di;
    AC.particle_diameter_ = 0.003;
    AC.eb_ = 0.456;
    AC.domain_.n_cells = 60;

    auto& AC_ads = GetAdsorbent(ads_id);
    AC_ads.solid_density_ = 799;
    AC_ads.heat_capacity_ = 1040;
    AC_ads.thermal_conductivity_ = 0.63;
    AC_ads.ep_ = 0.0;

    data.bed.wall_.heat_capacity_ = 902;
    data.bed.wall_.density_ = 2700;
    data.bed.wall_.thermal_conductivity_ = 205;
    data.bed.wall_.thickness_ = 3.25e-3;
    data.bed.ambient_temperature_ = 291.4; // Table 1: T_a = 291.4 K (was 293 = T_F, a transcription slip)
    data.fluid.adsorbed_heat_capacity_["CarbonDioxide"] = 36.96;

    auto dsl_co2 = std::make_unique<DualSiteLangmuirIsotherm_Concentration>("CarbonDioxide");
    dsl_co2->params_["CarbonDioxide"].m1 = 0.418;
    dsl_co2->params_["CarbonDioxide"].b1 = 2.05e-6;
    dsl_co2->params_["CarbonDioxide"].q1 = 29980;
    dsl_co2->params_["CarbonDioxide"].m2 = 6.04;
    dsl_co2->params_["CarbonDioxide"].b2 = 3.32e-6;
    dsl_co2->params_["CarbonDioxide"].q2 = 20367;
    AddIsotherm({ "CarbonDioxide" }, ads_id, std::move(dsl_co2));

    AddLinearDrivingForceModel_Haynes(ads_id, "CarbonDioxide", 3, 0.59);
    AddLinearDrivingForceModel(ads_id, "Air", 0.0);

    SetHOA(ads_id, "CarbonDioxide", 29980.0);
    SetHOA(ads_id, "Air", 0.0);

    auto& layer = GetLayer(layer_id);
    layer.equations_.dispersion_model_ = CreateAxialDispersionModel(AxialDispersionType::RUTHVEN);
    layer.equations_.heat_dispersion_model_ = CreateAxialHeatDispersionModel(AxialHeatDispersionType::YAGI_AND_KUNII);

    AddConstantInternalHeatTransferCoefficientModel(45.2, layer_id);
    auto ihtc = std::make_unique<ExternalHeatTransferModel_Holman>();
    ihtc->c = 0.027;
    ihtc->m = 0.8;
    ihtc->n = 1.0/3.0;
    ihtc->external_velocity = 0.5;
    layer.equations_.external_heat_transfer_model_ = std::move(ihtc);

    AddInitialConditions(
      1.0e5, 291.4, 291.4,
      {
        {"CarbonDioxide", 0.0 },
        {"Air", 1.0 }
      }
    );

    AddStep("Feed", 3000.0);

    const double volumetric_flowrate_ccm = 200;
    const double volumetric_flowrate_m3_s = volumetric_flowrate_ccm * 1.666666666666667e-8;
    
    Step& step_feed = data.cycle.GetStep("Feed");
    auto feed_inlet_stream = CreateNoValveInletStream(
      CreateVolumeFlowBC(volumetric_flowrate_m3_s),
      {
        {"CarbonDioxide", 1.0 },
        {"Air", 0.0 }
      },
      293.0
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
    feed_pressure_out->pressure_ = CreateConstantBoundaryValue(1e5);
    feed_pressure_out->mole_fraction_["CarbonDioxide"] = 0.0;
    feed_pressure_out->mole_fraction_["Air"] = 1.0;
    step_feed.top_boundary_condition_ = std::move(feed_pressure_out);

    ConfigureCycle(1);

    SaveExperimentalData("data/ward_2022/Qin=200");

    CoupledSolver solver{data};

    std::string error{};
    EXPECT_TRUE(data.PerformChecks(error)) << "System checks failed: " << error;
    
    bool verbose = false;
    int retval = solver.IntegrateCycle(data.cycle, verbose);
    EXPECT_GE(retval, 0) << "Solver did not complete";

    EXPECT_TRUE(IsMassConserved(solver.GetLog())) << "Mass is not conserved.";
    EXPECT_TRUE(IsMolesConserved(solver.GetLog())) << "Mass is not conserved.";
    EXPECT_TRUE(IsComponentMolesConserved(solver.GetLog())) << "Component moles not conserved.";

    // solver.PrintStatistics();
    SaveTest("Ward_2022_Q=200[ccm]", solver);
    SaveResultsAsExample("Ward_2022_Q=200[ccm]", solver);

    ValidateWard(solver, "Q200ccm", 200, kTolY_200, kTolT_200, kTolT50_200);
  }

  // Reproduces a suspected bug: the total CO2 accumulated in the bed should be
  // independent of the mass transfer coefficient (given sufficient simulation time
  // for both cases to reach equilibrium).  A left-shift in the breakthrough curve
  // when MTC is reduced suggests the equilibrium loading is being affected by the
  // MTC, which is physically incorrect.
  TEST_F(Ward2022Test, MTC_Invariance) {
    // Helper that configures the Ward 2022 Q=80 ccm setup on 'data' with a
    // given CO2 mass transfer coefficient, runs the solver, and returns the
    // total moles of CO2 accumulated (gas + solid) during the feed step.
    auto RunAndGetCO2Accumulated = [&](double co2_mtc, double step_time = 3000.0) -> double {
      // Reset fixture data to a clean state
      data = UserData{};
      SetUp();

      AddComponents({ "CarbonDioxide", "Air" });

      const double Di = 0.0189;
      auto [layer_id, ads_id] = AddLayer("Activated Carbon");
      auto& AC = GetLayer(layer_id);
      AC.domain_.length = 0.286;
      std::get<CylinderBed>(data.bed.shape_).diameter = Di;
      AC.particle_diameter_ = 0.003;
      AC.eb_ = 0.456;
      AC.domain_.n_cells = 60;

      auto& AC_ads = GetAdsorbent(ads_id);
      AC_ads.solid_density_ = 799;
      AC_ads.heat_capacity_ = 1040;
      AC_ads.thermal_conductivity_ = 0.63;
      AC_ads.ep_ = 0.0;

      data.bed.wall_.heat_capacity_ = 902;
      data.bed.wall_.density_ = 2700;
      data.bed.wall_.thermal_conductivity_ = 205;
      data.bed.wall_.thickness_ = 3.25e-3;
      data.bed.ambient_temperature_ = 293;

      auto dsl_co2 = std::make_unique<DualSiteLangmuirIsotherm_Concentration>("CarbonDioxide");
      dsl_co2->params_["CarbonDioxide"].m1 = 0.418;
      dsl_co2->params_["CarbonDioxide"].b1 = 2.05e-6;
      dsl_co2->params_["CarbonDioxide"].q1 = 29980;
      dsl_co2->params_["CarbonDioxide"].m2 = 6.04;
      dsl_co2->params_["CarbonDioxide"].b2 = 3.32e-6;
      dsl_co2->params_["CarbonDioxide"].q2 = 20367;
      AddIsotherm({ "CarbonDioxide" }, ads_id, std::move(dsl_co2));

      AddLinearDrivingForceModel(ads_id, "CarbonDioxide", co2_mtc);
      AddLinearDrivingForceModel(ads_id, "Air", 0.0);

      AddInitialConditions(
        1.0e5, 293.0, 293.0,
        {
          {"CarbonDioxide", 0.0 },
          {"Air", 1.0 }
        }
      );

      AddStep("Feed", step_time);

      const double volumetric_flowrate_m3_s = 80.0 * 1.666666666666667e-8;

      Step& step_feed = data.cycle.GetStep("Feed");
      auto feed_inlet_stream = CreateNoValveInletStream(
        CreateVolumeFlowBC(volumetric_flowrate_m3_s),
        { {"CarbonDioxide", 1.0}, {"Air", 0.0} },
        293.0
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
      feed_pressure_out->pressure_ = CreateConstantBoundaryValue(1e5);
      feed_pressure_out->mole_fraction_["CarbonDioxide"] = 0.0;
      feed_pressure_out->mole_fraction_["Air"] = 1.0;
      step_feed.top_boundary_condition_ = std::move(feed_pressure_out);

      ConfigureCycle(1);

      CoupledSolver solver{data};
      std::string error{};
      EXPECT_TRUE(data.PerformChecks(error)) << "System checks failed: " << error;

      int retval = solver.IntegrateCycle(data.cycle, false);
      EXPECT_GE(retval, 0) << "Solver did not complete (MTC=" << co2_mtc << ")";

      const auto& log = solver.GetLog();
      double gas   = CalculateComponentAccumulation_Volume(log, 0, 0, ComponentVariableType::COMPONENT_MOLAR_CONCENTRATION, "CarbonDioxide");
      double solid = CalculateComponentAccumulation_Mass(log, 0, 0, ComponentVariableType::ADSORBED_CONCENTRATION, "CarbonDioxide");
      double total = gas + solid;

      std::cout << "[MTC=" << co2_mtc << "] CO2 accumulated:"
                << "  gas="   << gas
                << "  solid=" << solid
                << "  total=" << total << "\n";

      std::string test_name = "Ward_2022_Q=80[ccm]_k=" + std::to_string(co2_mtc);
      SaveTest(test_name, solver);

      return total;
    };

    double co2_high_mtc = RunAndGetCO2Accumulated(0.1, 3000.0);
    // Use much longer time for low MTC to allow solid to reach equilibrium.
    // For k=0.001, time constant = 1000 s, so 20000 s = 20 time constants.
    double co2_low_mtc  = RunAndGetCO2Accumulated(0.001, 20000.0);

    // At t=3000 s both cases have been running for hundreds of kinetic half-lives
    // (t_half = ln2/k = 69 s for k=0.01), so both beds should be at the same
    // equilibrium loading.  Any significant difference indicates a solver bug.
    double tolerance = 0.01 * std::max(std::abs(co2_high_mtc), std::abs(co2_low_mtc)); // 1 % relative
    EXPECT_NEAR(co2_high_mtc, co2_low_mtc, tolerance)
      << "Total CO2 accumulated differs between MTC=0.1 (" << co2_high_mtc
      << " mol) and MTC=0.01 (" << co2_low_mtc << " mol)."
      << " Difference = " << std::abs(co2_high_mtc - co2_low_mtc) << " mol.";
  }

} // namespace adsim
