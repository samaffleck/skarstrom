// STL includes
#define _USE_MATH_DEFINES
#include <cmath>
#include <filesystem>
#include <limits>

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

// Paper fixture (Streb et al. 2021) — see README.md in this folder.
#include "../../fixtures/user_data_test.h"
#include "../../fixtures/validation_metrics.h"

// ============================================================================
//  Streb & Mazzotti (2021), Adsorption 27:541-558, doi:10.1007/s10450-021-00306-y
//  Breakthrough experiments on zeolite 13X (Table 2):
//    A1-A7: binary CO2:CH4 = 50:50;  B1-B7: ternary H2:CO2:CH4 = 75:20:5.
//  Digitised outlet compositions and five centre-line thermocouples
//  (z = 10, 35, 60, 85, 110 cm) live in core/data/streb_2021/<case>/.
//
//  Every case is a PREDICTION with the parameters estimated on A3 (Table 4) and
//  the extended-Sips competitive isotherm of Table 3. Note the paper itself
//  finds IAST/RAST predict the binary cases better than extended isotherms
//  (Sect. 4.2.2) and fits the feed velocity per experiment within the MFC
//  uncertainty (Fig. 7, 12); adsim uses the MFC set-point and extended Sips,
//  so some disagreement in the CH4 front at 10-25 bar is expected. The paper
//  also models ~10.5 m of downstream piping to the MS (Table 1); adsim does not,
//  which shifts the simulated outlet curves slightly earlier.
//
//  Metrics (recorded in the validation summary):
//    * outlet mole-fraction NRMSE per component = RMSE / feed mole fraction,
//    * thermocouple RMSE [K] pooled over the five thermocouples.
//  Tolerances: measured value rounded up modestly (see comments per case).
// ============================================================================

namespace adsim {

  class Streb2021Test : public UserDataTest {
  public:
    struct StrebCase {
      std::string name;          // "A1", ...
      double pressure;           // [Pa]
      double temperature;        // [K] feed = initial = thermostat
      double vol_flow;           // [m3/s] at column P and T (Table 2)
      bool ternary;              // false: CO2:CH4 = 50:50, true: H2:CO2:CH4 = 75:20:5
      double duration;           // [s]
      // Regression tolerances (measured error, rounded up; see TEST_F comments)
      double tol_nrmse_co2;
      double tol_nrmse_ch4;
      double tol_nrmse_h2;
      double tol_temperature_rmse;
      double tol_t50_co2;
      double tol_t50_ch4;
    };

    // Shared column / isotherm / transport set-up (Table 1, 3, 4). Kept in the
    // exact order of the original A1/A3 tests so their example decks (which the
    // WASM-parity and coloured-Jacobian tests also load) are unchanged.
    void SetUpStreb(const StrebCase& c) {
      if (c.ternary) AddComponents({ "CarbonDioxide", "Methane", "Hydrogen", "Helium" });
      else           AddComponents({ "CarbonDioxide", "Methane", "Helium" });

      // The paper's model uses "constant, average values" for the heat
      // capacity of the gas phase C_g (model assumptions) and prints no value.
      // Our choice of the average: each component's ideal-gas C_p0 at the
      // case's feed (= initial) temperature, held for the run (the engine's
      // constant gas heat capacity mode; see README).
      data.fluid.gas_heat_capacity_mode_ = GasHeatCapacityMode::CONSTANT;
      for (const auto& component : data.fluid.components_)
        data.fluid.gas_heat_capacity_[component] =
          Fluid::EvalIdealGasHeatCapacity(*Fluid::IdealGasCoeffs(component), c.temperature);

      auto [layer_id, ads_id] = AddLayer("13X");
      double skeletal_density = 2359;
      double solid_density = 1085;
      double bed_density = 595;
      double ep = 1 - solid_density / skeletal_density;
      double eb = 1 - bed_density / solid_density;

      auto& layer = GetLayer(layer_id);
      layer.domain_.length = kLength;
      std::get<CylinderBed>(data.bed.shape_).diameter = 0.025;
      layer.particle_diameter_ = 0.0018;
      layer.eb_ = eb;
      layer.domain_.n_cells = kCells;

      auto& ads_13X = GetAdsorbent(ads_id);
      ads_13X.solid_density_ = solid_density;
      ads_13X.heat_capacity_ = 1100;
      ads_13X.thermal_conductivity_ = 0.15;
      ads_13X.ep_ = ep;

      data.bed.wall_.heat_capacity_ = 500;
      data.bed.wall_.density_ = 8000;
      data.bed.wall_.thermal_conductivity_ = 0; // They don't include wall conductions, so set to 0
      data.bed.wall_.thickness_ = 2.5e-3;
      data.bed.equations_.wall_temperature_boundary_condition_ = CreateWallTemperatureBoundaryCondition(WallTemperatureBoundaryConditionType::ZERO_FLUX);
      data.bed.ambient_temperature_ = c.temperature;

      // Extended Sips2 isotherm, Table 3 (B1 converted from 1e-4 1/bar to 1/Pa,
      // B2 from kJ/mol to J/mol).
      auto sips = std::make_unique<ExtendedSips2Isotherm>();
      sips->params_["CarbonDioxide"].A1 = 7.268;
      sips->params_["CarbonDioxide"].A2 = -0.6168;
      sips->params_["CarbonDioxide"].B1 = 1.129e-9;
      sips->params_["CarbonDioxide"].B2 = 28389;
      sips->params_["CarbonDioxide"].C1 = 0.425;
      sips->params_["CarbonDioxide"].C2 = 0.7238;
      sips->params_["CarbonDioxide"].T_ref = 298;

      sips->params_["Methane"].A1 = 4.473;
      sips->params_["Methane"].A2 = -0.6569;
      sips->params_["Methane"].B1 = 3.605e-9;
      sips->params_["Methane"].B2 = 15561;
      sips->params_["Methane"].C1 = 1.073;
      sips->params_["Methane"].C2 = 0.0;
      sips->params_["Methane"].T_ref = 298;

      if (c.ternary) {
        sips->params_["Hydrogen"].A1 = 5.013;
        sips->params_["Hydrogen"].A2 = 0.0;
        sips->params_["Hydrogen"].B1 = 1.034e-9;
        sips->params_["Hydrogen"].B2 = 9453;
        sips->params_["Hydrogen"].C1 = 1.006;
        sips->params_["Hydrogen"].C2 = 0.0;
        sips->params_["Hydrogen"].T_ref = 298;
        AddIsotherm({ "CarbonDioxide", "Methane", "Hydrogen" }, ads_id, std::move(sips));
      } else {
        AddIsotherm({ "CarbonDioxide", "Methane" }, ads_id, std::move(sips));
      }

      // LDF coefficients estimated on A3 (Table 4); H2 assumed 1 1/s (Table 4, b).
      AddLinearDrivingForceModel(ads_id, "CarbonDioxide", 0.06);
      AddLinearDrivingForceModel(ads_id, "Methane", 0.2);
      if (c.ternary) AddLinearDrivingForceModel(ads_id, "Hydrogen", 1.0);

      // Loading-averaged isosteric heats (Sect. 4.1).
      SetHOA(ads_id, "CarbonDioxide", 36.3e3);
      SetHOA(ads_id, "Methane", 18.8e3);
      if (c.ternary) SetHOA(ads_id, "Hydrogen", 10.7e3);

      layer.equations_.dispersion_model_ = CreateAxialDispersionModel(AxialDispersionType::CONSTANT);
      layer.equations_.heat_dispersion_model_ = CreateAxialHeatDispersionModel(AxialHeatDispersionType::CONSTANT);

      layer.equations_.internal_heat_transfer_model_ = CreateInternalHeatTransferModel(InternalHeatTransferType::DIXON_1996);
      AddConstantExternalHeatTransferCoefficientModel(220.0, layer_id); // [W/(m-K)]

      std::unordered_map<std::string, double> y0 = {
        {"CarbonDioxide", 0.0 }, {"Methane", 0.0 }, {"Helium", 1.0 } };
      if (c.ternary) y0["Hydrogen"] = 0.0;
      AddInitialConditions(c.pressure, c.temperature, c.temperature, y0);

      AddStep("Feed", c.duration);

      Step& step_feed = data.cycle.GetStep("Feed");
      std::unordered_map<std::string, double> feed = c.ternary
        ? std::unordered_map<std::string, double>{ {"CarbonDioxide", 0.20}, {"Methane", 0.05}, {"Hydrogen", 0.75}, {"Helium", 0.0} }
        : std::unordered_map<std::string, double>{ {"CarbonDioxide", 0.5 }, {"Methane", 0.5 }, {"Helium", 0.0 } };
      auto feed_inlet_stream = CreateNoValveInletStream(CreateVolumeFlowBC(c.vol_flow), feed, c.temperature);
      feed_inlet_stream->name_ = c.ternary ? "Ternary Feed" : "Binary Feed";

      auto feed_flow_bc = std::make_unique<FlowBoundaryCondition>();
      feed_flow_bc->direction_ = BoundaryDirection::INLET;
      feed_flow_bc->location_ = BoundaryLocation::BOTTOM;
      feed_flow_bc->streams_.push_back(std::move(feed_inlet_stream));
      step_feed.bottom_boundary_condition_ = std::move(feed_flow_bc);

      auto feed_pressure_out = std::make_unique<PressureBoundaryCondition>();
      feed_pressure_out->direction_ = BoundaryDirection::OUTLET;
      feed_pressure_out->location_ = BoundaryLocation::TOP;
      feed_pressure_out->pressure_ = CreateConstantBoundaryValue(c.pressure);
      feed_pressure_out->mole_fraction_["CarbonDioxide"] = 0.0;
      feed_pressure_out->mole_fraction_["Methane"] = 0.0;
      if (c.ternary) feed_pressure_out->mole_fraction_["Hydrogen"] = 0.0;
      feed_pressure_out->mole_fraction_["Helium"] = 1.0;
      step_feed.top_boundary_condition_ = std::move(feed_pressure_out);

      ConfigureCycle(1);

      SaveExperimentalData("data/streb_2021/" + c.name);
    }

    // Run the case, check conservation and record/assert the validation metrics.
    void RunAndValidate(const StrebCase& c, bool save_example) {
      CoupledSolver solver{data};

      std::string error{};
      EXPECT_TRUE(data.PerformChecks(error)) << "System checks failed: " << error;

      bool verbose = false;
      int retval = solver.IntegrateCycle(data.cycle, verbose);
      ASSERT_GE(retval, 0) << "Solver did not complete";

      EXPECT_TRUE(IsMassConserved(solver.GetLog())) << "Mass is not conserved.";
      EXPECT_TRUE(IsMolesConserved(solver.GetLog())) << "Moles are not conserved.";
      EXPECT_TRUE(IsComponentMolesConserved(solver.GetLog())) << "Component moles not conserved.";

      if (save_example) SaveResultsAsExample("Streb2021_" + c.name, solver);

      Validate(solver, c);
    }

    void Validate(const CoupledSolver& solver, const StrebCase& c) {
      using namespace validation;
      CaseReport report("streb_2021", "Experiment" + c.name);
      const std::string dir = "data/streb_2021/" + c.name + "/" + c.name;
      const std::string fig = c.ternary ? "Fig. 11 (digitised)" : (c.name == "A3" ? "Fig. 6 (digitised)" : "Fig. 7 (digitised)");

      const double y_co2 = c.ternary ? 0.20 : 0.5;
      const double y_ch4 = c.ternary ? 0.05 : 0.5;
      report.CheckCurve("outlet y_CO2", FaceMoleFraction(solver, "CarbonDioxide", 0, 0),
                        LoadSeriesCsv(dir + "_mole_fraction_CO2.csv"), c.tol_nrmse_co2, "-", fig, y_co2);
      report.CheckCurve("outlet y_CH4", FaceMoleFraction(solver, "Methane", 0, 0),
                        LoadSeriesCsv(dir + "_mole_fraction_CH4.csv"), c.tol_nrmse_ch4, "-", fig, y_ch4);
      report.CheckCrossingTime("CO2 front time (y_CO2 = y_feed/2)", FaceMoleFraction(solver, "CarbonDioxide", 0, 0),
                               LoadSeriesCsv(dir + "_mole_fraction_CO2.csv"), 0.5 * y_co2, c.tol_t50_co2, fig);
      report.CheckCrossingTime("CH4 front time (y_CH4 = y_feed/2)", FaceMoleFraction(solver, "Methane", 0, 0),
                               LoadSeriesCsv(dir + "_mole_fraction_CH4.csv"), 0.5 * y_ch4, c.tol_t50_ch4, fig);
      if (c.ternary)
        report.CheckCurve("outlet y_H2", FaceMoleFraction(solver, "Hydrogen", 0, 0),
                          LoadSeriesCsv(dir + "_mole_fraction_H2.csv"), c.tol_nrmse_h2, "-", fig, 0.75);

      // Pooled RMSE over the five centre-line thermocouples.
      double sum_sq = 0.0;
      size_t n_tot = 0;
      for (int z_cm : { 10, 35, 60, 85, 110 }) {
        Series sim = CellVariableAt(solver, VariableType::TEMPERATURE, z_cm / 100.0, kLength, kCells);
        Series exp = LoadSeriesCsv(dir + "_temperature_z=" + std::to_string(z_cm) + "[cm].csv");
        size_t n = 0;
        const double rmse = CurveRMSE(sim, exp, &n);
        if (n > 0 && std::isfinite(rmse)) { sum_sq += rmse * rmse * n; n_tot += n; }
      }
      const double t_rmse = n_tot ? std::sqrt(sum_sq / n_tot) : std::numeric_limits<double>::quiet_NaN();
      report.CheckError("thermocouple T (5 positions, pooled)", t_rmse, "rmse", c.tol_temperature_rmse, "K", fig,
                        "RMSE pooled over " + std::to_string(n_tot) +
                        " digitised points of the centre-line thermocouples at z = 10, 35, 60, 85, 110 cm");

      if (c.name == "A3") {
        // Text-stated (approximate) breakthrough times / temperature rises —
        // recorded for context only (no assertion: 'about', 'after' in the text).
        const Series co2 = FaceMoleFraction(solver, "CarbonDioxide", 0, 0);
        const Series ch4 = FaceMoleFraction(solver, "Methane", 0, 0);
        report.Info("CH4 breakthrough time (y > 0.01)", CrossingTime(ch4, 0.01),
                    Reference("streb_2021", "ExperimentA3", "CH4_breakthrough_time"), "abs_error",
                    "simulated: first time outlet y_CH4 > 0.01");
        report.Info("CO2 breakthrough time (y > 0.01)", CrossingTime(co2, 0.01),
                    Reference("streb_2021", "ExperimentA3", "CO2_breakthrough_time"), "abs_error",
                    "simulated: first time outlet y_CO2 > 0.01");
      }
    }

    static constexpr double kLength = 1.2;
    static constexpr int kCells = 160;  // grid-converged with the default WENO scheme (README.md "Grid")
  };

  namespace {
    // Ternary-only tolerance slot for the binary cases.
    constexpr double kNotUsed = std::numeric_limits<double>::infinity();
  }

  // Reference experiment (Table 2: 25 C, 10 bar, 20 cm3/s) — the one the paper
  // estimated the LDF coefficients on (Fig. 6).
  TEST_F(Streb2021Test, ExperimentA3) {
    data.info.title = "Streb 2021 - Experiment A3";
    data.info.description = "Binary CO2/CH4 separation breakthrough test based on Streb et al. 2021, experiment A3. "
                            "This simulation models a breakthrough experiment for 13X zeolite in a lab-scale test rig at elevated pressure (10 bar). "
                            "The bed is initially saturated with helium and then fed with an equimolar binary mixture of CO2 and CH4. "
                            "The adsorbent uses an Extended Sips2 isotherm model for both components and linear driving force models for mass transfer. "
                            "The simulation includes coupled heat and mass transfer with internal and external heat transfer mechanisms. "
                            "Wall heat transfer is neglected in this configuration. Experimental breakthrough data is included and can be plotted "
                            "to validate the solver performance. This test demonstrates the ability to model competitive adsorption in binary systems "
                            "at high pressure conditions.";
    data.info.author = "Sam Affleck";
    // Tolerances = measured error rounded up, with WENO on 160 cells (README.md "Grid and advection scheme") (measured: NRMSE y_CO2 0.184339, y_CH4 0.603105, T RMSE 12.0664 K, front time CO2 0.046863, CH4 0.030998; upwind on 60 cells: y_CH4 0.5591).
    StrebCase c{ "A3", 10e5, 298, 20e-6, false, 4100.0,
                 0.27, 0.73, kNotUsed, 19, 0.08, 0.041 };
    SetUpStreb(c);
    RunAndValidate(c, true);
  }

  // Table 2: 25 C, 1 bar, 20 cm3/s.
  TEST_F(Streb2021Test, ExperimentA1) {
    data.info.title = "Streb 2021 - Experiment A1";
    data.info.description = "Binary CO2/CH4 separation breakthrough test based on Streb et al. 2021, experiment A1. "
                            "This simulation models a breakthrough experiment for 13X zeolite in a lab-scale test rig at atmospheric pressure (1 bar). "
                            "This is identical to experiment A3 except operated at lower pressure, allowing comparison of pressure effects on separation performance. "
                            "The bed is initially saturated with helium and then fed with an equimolar binary mixture of CO2 and CH4. "
                            "The adsorbent uses an Extended Sips2 isotherm model for both components and linear driving force models for mass transfer. "
                            "The simulation includes coupled heat and mass transfer with internal and external heat transfer mechanisms. "
                            "Wall heat transfer is neglected in this configuration. Experimental breakthrough data is included and can be plotted "
                            "to validate the solver performance. This test demonstrates the ability to model competitive adsorption in binary systems "
                            "at atmospheric pressure conditions.";
    data.info.author = "Sam Affleck";
    // Tolerances = measured error rounded up (measured: NRMSE y_CO2 0.456647, y_CH4 0.721357, T RMSE 21.1045 K, front time CO2 0.0843783, CH4 0.0970627).
    StrebCase c{ "A1", 1e5, 298, 20e-6, false, 4100.0,
                 0.53, 0.83, kNotUsed, 25, 0.098, 0.12 };
    SetUpStreb(c);
    RunAndValidate(c, true);
  }

  // The remaining binary (A2, A4-A7) and ternary (B1, B2) experiments of Table 2
  // for which digitised data exist. Durations cover the digitised record.
  TEST_F(Streb2021Test, ExperimentA2) {
    data.info.title = "Streb 2021 - Experiment A2";
    // Tolerances = measured error rounded up (re-measured after the molar-basis gas heat capacity fix, 2026-11; see README.md) (measured: NRMSE y_CO2 0.350419, y_CH4 0.838922, T RMSE 18.5785 K, front time CO2 0.0546649, CH4 0.110449; before: y_CO2 0.33839).
    StrebCase c{ "A2", 10e5, 298, 10e-6, false, 1500.0,
                 0.37, 0.96, kNotUsed, 22, 0.058, 0.13 };
    SetUpStreb(c);
    RunAndValidate(c, false);
  }

  TEST_F(Streb2021Test, ExperimentA4) {
    data.info.title = "Streb 2021 - Experiment A4";
    // Tolerances = measured error rounded up, with WENO on 160 cells (README.md "Grid and advection scheme") (measured: NRMSE y_CO2 0.290623, y_CH4 0.792805, T RMSE 16.6768 K, front time CO2 0.100904, CH4 0.177891; upwind on 60 cells: y_CO2 0.2865).
    StrebCase c{ "A4", 10e5, 298, 40e-6, false, 1000.0,
                 0.35, 0.91, kNotUsed, 20, 0.105, 0.195 };
    SetUpStreb(c);
    RunAndValidate(c, false);
  }

  TEST_F(Streb2021Test, ExperimentA5) {
    data.info.title = "Streb 2021 - Experiment A5";
    // Tolerances = measured error rounded up. y_CO2 NRMSE and CO2 front time
    // re-measured with the paper's constant gas heat capacity (Cp0 at the feed
    // temperature; 2026-11, owner-approved; see README.md): y_CO2 0.184734,
    // CO2 front 0.0647726 (were 0.1728 and 0.0594 with the held CoolProp value).
    // Earlier (molar-basis gas heat capacity fix): y_CH4 0.659749, T RMSE
    // 16.0716 K, CH4 front 0.148211.
    StrebCase c{ "A5", 25e5, 298, 20e-6, false, 750.0,
                 0.19, 0.76, kNotUsed, 21, 0.066, 0.17 };
    SetUpStreb(c);
    RunAndValidate(c, false);
  }

  TEST_F(Streb2021Test, ExperimentA6) {
    data.info.title = "Streb 2021 - Experiment A6";
    // Tolerances = measured error rounded up. y_CO2 NRMSE and CO2 front time
    // re-measured with the paper's constant gas heat capacity (Cp0 at the feed
    // temperature; 2026-11, owner-approved; see README.md): y_CO2 0.170486,
    // CO2 front 0.0292289 (were 0.1608 and 0.0273 with the held CoolProp value).
    // Earlier (molar-basis gas heat capacity fix): y_CH4 0.749957, T RMSE
    // 11.7501 K, CH4 front 0.122827.
    StrebCase c{ "A6", 10e5, 318.15, 20e-6, false, 1050.0,
                 0.175, 0.85, kNotUsed, 13, 0.030, 0.14 };
    SetUpStreb(c);
    RunAndValidate(c, false);
  }

  TEST_F(Streb2021Test, ExperimentA7) {
    data.info.title = "Streb 2021 - Experiment A7";
    // Tolerances = measured error rounded up (re-measured after the Danckwerts inlet and Darcy flow-face fixes, 2026-11; see README.md) (measured: NRMSE y_CO2 0.0794838, y_CH4 0.692014, T RMSE 11.1461 K, front time CO2 0.0195941, CH4 0.102098; before: CH4 front 0.0726965).
    StrebCase c{ "A7", 10e5, 338.15, 20e-6, false, 1100.0,
                 0.14, 0.72, kNotUsed, 14, 0.04, 0.12 };
    SetUpStreb(c);
    RunAndValidate(c, false);
  }

  TEST_F(Streb2021Test, ExperimentB1) {
    data.info.title = "Streb 2021 - Experiment B1";
    // Tolerances = measured error rounded up, with WENO on 160 cells (README.md "Grid and advection scheme") (measured: NRMSE y_CO2 0.376751, y_CH4 0.126405, y_H2 0.372369, T RMSE 17.6917 K, front time CO2 0.06425, CH4 0.001939; upwind on 60 cells: CH4 front 0.02382).
    StrebCase c{ "B1", 1e5, 298, 20e-6, true, 10400.0,
                 0.45, 0.15, 0.42, 20, 0.081, 0.004 };
    SetUpStreb(c);
    RunAndValidate(c, false);
  }

  TEST_F(Streb2021Test, ExperimentB2) {
    data.info.title = "Streb 2021 - Experiment B2";
    // Tolerances = measured error rounded up (measured: NRMSE y_CO2 0.488768, y_CH4 0.297433, y_H2 0.399014, T RMSE 32.6309 K, front time CO2 0.0912807, CH4 0.121991). After the WP19-1 competitive-Sips trace fix (2026-11): y_CO2 0.44785, y_CH4 0.32043, y_H2 0.39777, T RMSE 27.005 K, front time CO2 0.1072, CH4 0.13975 (before it, on the same main: y_CH4 0.31859, CH4 front 0.13722).
    StrebCase c{ "B2", 10e5, 298, 10e-6, true, 3350.0,
                 0.57, 0.35, 0.46, 38, 0.11, 0.15 };
    SetUpStreb(c);
    RunAndValidate(c, false);
  }

} // namespace adsim
