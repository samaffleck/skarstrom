// STL includes
#define _USE_MATH_DEFINES
#include <cmath>
#include <map>
#include <filesystem>
#include <limits>
#include <stdexcept>
#include <string>
#include <vector>

// AdSim includes
#include "../../../src/pch.h"
#include "adsim/solver/optimisation.h"
#include "adsim/solver/coupled_solver.h"
#include "adsim/solver/kpis.h"
#include "adsim/math/constants.h"
#include "adsim/equations/isotherm.h"
#include "adsim/equations/iast.h"
#include "adsim/objects/porous_media.h"
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
//  CALIBRATED on A3 (README, "Calibration on A3"): the LDF coefficients of CO2
//  and CH4, the wall-to-thermofluid coefficient and the volume of the sampling
//  line to the MS (an inert, nearly empty layer above the bed, run as plug
//  flow) are fitted to the outlet curves and thermocouples of A3 alone, with
//  A3's feed flow held within 2 % of the set point. Every other case is a
//  PREDICTION with those values and the extended-Sips competitive isotherm of
//  Table 3, except its feed flow: like the paper, which fits the feed velocity
//  per experiment (Sect. 3.2.1), each case's flow factor is fitted, here to the
//  speed of its own CO2 temperature front between the thermocouples at 35 and
//  110 cm. That uses none of the outlet data. Six cases are more than 2 % from
//  the set point (A1 +5.1, A2 -6.1, A4 -9.4, A6 -3.5, B1 +5.7, B2 -11.3 %); the
//  paper's bound is 2 % of the unprinted MFC full scale (README). The
//  *_IAST_RecordOnly tests run every case under IAST over the pure
//  Sips terms of the same parameters and record it without assertion (README,
//  "IAST versus extended Sips").
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
      double flow_factor;        // actual / Table 2 flow, from the CO2 thermal-front speed (README "Calibration on A3")
      // Regression tolerances (measured error, rounded up; see TEST_F comments)
      double tol_nrmse_co2;
      double tol_nrmse_ch4;
      double tol_nrmse_h2;
      double tol_temperature_rmse;
      double tol_t50_co2;
      double tol_t50_ch4;
    };

    // Competition model of the adsorbent: false = the extended-Sips isotherm
    // the decks use (Table 3, Eqs. 1-4, 6); true = IAST over the pure Sips
    // terms of the same parameters (the paper's Sect. 4.2.1 alternative;
    // recorded only, README "IAST versus extended Sips").
    bool iast_ = false;

    static std::vector<std::string> StrebAdsorbing(bool ternary) {
      if (ternary) return { "CarbonDioxide", "Methane", "Hydrogen" };
      return { "CarbonDioxide", "Methane" };
    }

    // Table 3 (CO2 from the paper's ref. [14]): A1 [mol/kg], A2 [-], B1 [1/Pa],
    // B2 [J/mol], C1, C2 [-], T_ref [K].
    static void SetTable3Sips(ParamMap<Sips2IsothermParameters>& params, bool ternary) {
      params["CarbonDioxide"].A1 = 7.268;
      params["CarbonDioxide"].A2 = -0.6168;
      params["CarbonDioxide"].B1 = 1.129e-9;
      params["CarbonDioxide"].B2 = 28389;
      params["CarbonDioxide"].C1 = 0.425;
      params["CarbonDioxide"].C2 = 0.7238;
      params["CarbonDioxide"].T_ref = 298;

      params["Methane"].A1 = 4.473;
      params["Methane"].A2 = -0.6569;
      params["Methane"].B1 = 3.605e-9;
      params["Methane"].B2 = 15561;
      params["Methane"].C1 = 1.073;
      params["Methane"].C2 = 0.0;
      params["Methane"].T_ref = 298;

      if (ternary) {
        params["Hydrogen"].A1 = 5.013;
        params["Hydrogen"].A2 = 0.0;
        params["Hydrogen"].B1 = 1.034e-9;
        params["Hydrogen"].B2 = 9453;
        params["Hydrogen"].C1 = 1.006;
        params["Hydrogen"].C2 = 0.0;
        params["Hydrogen"].T_ref = 298;
      }
    }

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

      // Table 3 Sips parameters (B1 converted from 1e-4 1/bar to 1/Pa, B2 from
      // kJ/mol to J/mol), as the extended-Sips competitive isotherm (Eqs. 1-4,
      // 6) or as pure Sips 2 isotherms under the IAST mixing rule.
      if (iast_) {
        auto sips = std::make_unique<Sips2Isotherm>();
        SetTable3Sips(sips->params_, c.ternary);
        AddIsotherm(StrebAdsorbing(c.ternary), ads_id, std::move(sips));
        GetAdsorbent(ads_id).adsorption_source_.mixing_rule_ = MixingRule::IAST;
      } else {
        auto sips = std::make_unique<ExtendedSips2Isotherm>();
        SetTable3Sips(sips->params_, c.ternary);
        AddIsotherm(StrebAdsorbing(c.ternary), ads_id, std::move(sips));
      }

      // LDF coefficients calibrated on A3 (README "Calibration on A3"; Table 4:
      // 0.06 and 0.2 1/s); H2 assumed 1 1/s (Table 4, b).
      AddLinearDrivingForceModel(ads_id, "CarbonDioxide", kLdfCO2);
      AddLinearDrivingForceModel(ads_id, "Methane", kLdfCH4);
      if (c.ternary) AddLinearDrivingForceModel(ads_id, "Hydrogen", 1.0);

      // Loading-averaged isosteric heats (Sect. 4.1).
      SetHOA(ads_id, "CarbonDioxide", 36.3e3);
      SetHOA(ads_id, "Methane", 18.8e3);
      if (c.ternary) SetHOA(ads_id, "Hydrogen", 10.7e3);

      layer.equations_.dispersion_model_ = CreateAxialDispersionModel(AxialDispersionType::CONSTANT);
      layer.equations_.heat_dispersion_model_ = CreateAxialHeatDispersionModel(AxialHeatDispersionType::CONSTANT);

      layer.equations_.internal_heat_transfer_model_ = CreateInternalHeatTransferModel(InternalHeatTransferType::DIXON_1996);
      // Wall-to-thermofluid coefficient calibrated on A3 (Sect. 3.2.1 takes
      // 220 W/(m2 K) from ref. [23]).
      AddConstantExternalHeatTransferCoefficientModel(kExternalHTC, layer_id); // [W/(m2 K)]

      // The sampling line to the MS (Table 1: 10.51 m of 6 mm tube, modelled by
      // the paper as plug flow at column pressure): an inert layer above the bed,
      // nearly empty (bed voidage 0.99), adiabatic and without dispersion, so it
      // carries the outlet gas as plug flow over its cells. Its gas volume
      // (length x area x voidage) is calibrated on A3. The engine's void region
      // is one well-mixed cell, which smears the CH4 front (README).
      {
        auto [line_id, line_ads] = AddLayer("Sampling line");
        auto& line = GetLayer(line_id);
        line.eb_ = kLineVoidFraction;
        line.domain_.length = kLineVolume / (kLineVoidFraction * kArea);
        line.domain_.n_cells = kLineCells;
        line.particle_diameter_ = 0.01;
        std::vector<std::string> line_components = StrebAdsorbing(c.ternary);
        line_components.push_back("Helium");
        auto line_inert = std::make_unique<InertIsotherm>();
        for (const auto& component : line_components) line_inert->params_[component];
        AddIsotherm(line_components, line_ads, std::move(line_inert));
        auto& line_solid = GetAdsorbent(line_ads);
        line_solid.solid_density_ = 1000;
        line_solid.heat_capacity_ = 500;
        line_solid.thermal_conductivity_ = 0.15;
        line_solid.ep_ = 0.0;
        line.equations_.dispersion_model_ = CreateAxialDispersionModel(AxialDispersionType::CONSTANT);
        line.equations_.heat_dispersion_model_ = CreateAxialHeatDispersionModel(AxialHeatDispersionType::CONSTANT);
        auto line_ihtc = std::make_unique<ConstantInternalHeatTransferModel>();
        line_ihtc->heat_transfer_coefficient_ = 0.0;
        line.equations_.internal_heat_transfer_model_ = std::move(line_ihtc);
        AddConstantExternalHeatTransferCoefficientModel(0.0, line_id);
      }

      std::unordered_map<std::string, double> y0 = {
        {"CarbonDioxide", 0.0 }, {"Methane", 0.0 }, {"Helium", 1.0 } };
      if (c.ternary) y0["Hydrogen"] = 0.0;
      AddInitialConditions(c.pressure, c.temperature, c.temperature, y0);

      AddStep("Feed", c.duration);

      Step& step_feed = data.cycle.GetStep("Feed");
      std::unordered_map<std::string, double> feed = c.ternary
        ? std::unordered_map<std::string, double>{ {"CarbonDioxide", 0.20}, {"Methane", 0.05}, {"Hydrogen", 0.75}, {"Helium", 0.0} }
        : std::unordered_map<std::string, double>{ {"CarbonDioxide", 0.5 }, {"Methane", 0.5 }, {"Helium", 0.0 } };
      auto feed_inlet_stream = CreateNoValveInletStream(CreateVolumeFlowBC(c.vol_flow * c.flow_factor), feed, c.temperature);
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
      BindStrebData(c.name);
    }

    // What each digitised dataset measures (experimental_data.h): the outlet
    // mole fractions and the centre-line thermocouples at z [cm] from the feed end.
    void BindStrebData(const std::string& case_name) {
      const std::map<std::string, std::string> species = {
          {"CH4", "Methane"}, {"CO2", "CarbonDioxide"}, {"H2", "Hydrogen"}};
      const std::string y_prefix = case_name + "_mole_fraction_";
      const std::string t_prefix = case_name + "_temperature_z=";
      for (auto& xy : data.exp_data_.GetAllData()) {
        ObservableBinding b;
        if (xy.name.rfind(y_prefix, 0) == 0 && species.count(xy.name.substr(y_prefix.size()))) {
          b.kind = ObservableKind::Component;
          b.variable = "mole_fraction";
          b.component = species.at(xy.name.substr(y_prefix.size()));
        } else if (xy.name.rfind(t_prefix, 0) == 0) {
          b.kind = ObservableKind::Bulk;
          b.variable = "temperature";
          b.location = ObservableLocation::Position;
          b.z = std::stod(xy.name.substr(t_prefix.size())) / 100.0;  // "10[cm]" -> 0.10 m
        } else {
          ADD_FAILURE() << "unbound Streb dataset " << xy.name;
          continue;
        }
        xy.observable = b;
      }
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

    // Under IAST (iast_) the case is recorded as Experiment<name>_IAST without
    // assertion: every tolerance is infinite and the compared quantities are
    // named "... vs experiment (IAST)".
    void Validate(const CoupledSolver& solver, const StrebCase& c) {
      using namespace validation;
      CaseReport report("streb_2021", "Experiment" + c.name + (iast_ ? "_IAST" : ""));
      const std::string tag = iast_ ? " vs experiment (IAST)" : "";
      const std::string dir = "data/streb_2021/" + c.name + "/" + c.name;
      const std::string fig = c.ternary ? "Fig. 11 (digitised)" : (c.name == "A3" ? "Fig. 6 (digitised)" : "Fig. 7 (digitised)");

      const double y_co2 = c.ternary ? 0.20 : 0.5;
      const double y_ch4 = c.ternary ? 0.05 : 0.5;
      report.CheckCurve(std::string("outlet y_CO2") + tag, FaceMoleFraction(solver, "CarbonDioxide", 0, 0),
                        LoadSeriesCsv(dir + "_mole_fraction_CO2.csv"), c.tol_nrmse_co2, "-", fig, y_co2);
      report.CheckCurve(std::string("outlet y_CH4") + tag, FaceMoleFraction(solver, "Methane", 0, 0),
                        LoadSeriesCsv(dir + "_mole_fraction_CH4.csv"), c.tol_nrmse_ch4, "-", fig, y_ch4);
      report.CheckCrossingTime(std::string("CO2 front time (y_CO2 = y_feed/2)") + tag, FaceMoleFraction(solver, "CarbonDioxide", 0, 0),
                               LoadSeriesCsv(dir + "_mole_fraction_CO2.csv"), 0.5 * y_co2, c.tol_t50_co2, fig);
      report.CheckCrossingTime(std::string("CH4 front time (y_CH4 = y_feed/2)") + tag, FaceMoleFraction(solver, "Methane", 0, 0),
                               LoadSeriesCsv(dir + "_mole_fraction_CH4.csv"), 0.5 * y_ch4, c.tol_t50_ch4, fig);
      if (c.ternary)
        report.CheckCurve(std::string("outlet y_H2") + tag, FaceMoleFraction(solver, "Hydrogen", 0, 0),
                          LoadSeriesCsv(dir + "_mole_fraction_H2.csv"), c.tol_nrmse_h2, "-", fig, 0.75);

      // Pooled RMSE over the five centre-line thermocouples.
      double sum_sq = 0.0;
      size_t n_tot = 0;
      for (int z_cm : { 10, 35, 60, 85, 110 }) {
        Series sim = CellVariableAt(solver, VariableType::TEMPERATURE, z_cm / 100.0, kLength, kCells, 0, 0, /*first_cell=*/1);
        Series exp = LoadSeriesCsv(dir + "_temperature_z=" + std::to_string(z_cm) + "[cm].csv");
        size_t n = 0;
        const double rmse = CurveRMSE(sim, exp, &n);
        if (n > 0 && std::isfinite(rmse)) { sum_sq += rmse * rmse * n; n_tot += n; }
      }
      const double t_rmse = n_tot ? std::sqrt(sum_sq / n_tot) : std::numeric_limits<double>::quiet_NaN();
      report.CheckError(std::string("thermocouple T (5 positions, pooled)") + tag, t_rmse, "rmse", c.tol_temperature_rmse, "K", fig,
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
    static constexpr double kArea = 0.25 * 3.14159265358979323846 * 0.025 * 0.025;  // column cross-section [m2]

    // Calibrated on A3 (README "Calibration on A3"; weighted least squares on
    // the outlet y_CO2, y_CH4 and the five thermocouples, A3's flow within 2 %).
    static constexpr double kLdfCO2 = 0.0592;        // [1/s] (+-4.7 %; Table 4: 0.06)
    static constexpr double kLdfCH4 = 0.137;         // [1/s] (+-13 %; Table 4: 0.2)
    static constexpr double kExternalHTC = 49.7;     // [W/(m2 K)] (+-6.7 %; Sect. 3.2.1: 220)
    static constexpr double kLineVolume = 1.42e-4;   // [m3] gas volume of the sampling line (+-0.03e-4)
    static constexpr double kLineVoidFraction = 0.99;
    static constexpr int kLineCells = 40;
    static constexpr int kCells = 160;  // grid-converged with the default WENO scheme (README.md "Grid")
  };

  namespace {
    // Ternary-only tolerance slot for the binary cases.
    constexpr double kNotUsed = std::numeric_limits<double>::infinity();

    // The experiments of Table 2 with digitised data: conditions and the
    // regression tolerances of the extended-Sips decks (measured error rounded
    // up). Durations cover the digitised record.
    Streb2021Test::StrebCase StrebCaseFor(const std::string& name) {
      using StrebCase = Streb2021Test::StrebCase;
      static const std::vector<StrebCase> cases = {
        // Tolerances = measured error x 1.1, rounded up; front times at least 0.002 (calibrated decks, README "Calibration on A3") (measured: NRMSE y_CO2 0.08258, y_CH4 0.1772, T RMSE 5.451 K, front time CO2 0.01426, CH4 0.02393).
        StrebCase{ "A3", 10e5, 298, 20e-6, false, 4100.0, 1.020,
                   0.1, 0.2, kNotUsed, 6, 0.02, 0.03 },
        // Tolerances = measured error x 1.1, rounded up; front times at least 0.002 (calibrated decks, README "Calibration on A3") (measured: NRMSE y_CO2 0.1646, y_CH4 0.2034, T RMSE 6.912 K, front time CO2 0.00639, CH4 0.01299).
        StrebCase{ "A1", 1e5, 298, 20e-6, false, 4100.0, 1.051,
                   0.19, 0.23, kNotUsed, 8, 0.01, 0.015 },
        // Tolerances = measured error x 1.1, rounded up; front times at least 0.002 (calibrated decks, README "Calibration on A3") (measured: NRMSE y_CO2 0.1046, y_CH4 0.264, T RMSE 8.513 K, front time CO2 0.007827, CH4 0.01764).
        StrebCase{ "A2", 10e5, 298, 10e-6, false, 1500.0, 0.939,
                   0.12, 0.3, kNotUsed, 9.5, 0.01, 0.02 },
        // Tolerances = measured error x 1.1, rounded up; front times at least 0.002 (calibrated decks, README "Calibration on A3") (measured: NRMSE y_CO2 0.02341, y_CH4 0.2891, T RMSE 7.072 K, front time CO2 0.003391, CH4 0.006905).
        StrebCase{ "A4", 10e5, 298, 40e-6, false, 1000.0, 0.906,
                   0.03, 0.32, kNotUsed, 8, 0.005, 0.01 },
        // Tolerances = measured error x 1.1, rounded up; front times at least 0.002 (calibrated decks, README "Calibration on A3") (measured: NRMSE y_CO2 0.04989, y_CH4 0.3121, T RMSE 11.4 K, front time CO2 0.01402, CH4 0.0101).
        StrebCase{ "A5", 25e5, 298, 20e-6, false, 750.0, 0.978,
                   0.06, 0.35, kNotUsed, 13, 0.02, 0.015 },
        // Tolerances = measured error x 1.1, rounded up; front times at least 0.002 (calibrated decks, README "Calibration on A3") (measured: NRMSE y_CO2 0.01319, y_CH4 0.1175, T RMSE 7.459 K, front time CO2 0.001863, CH4 0.001715).
        StrebCase{ "A6", 10e5, 318.15, 20e-6, false, 1050.0, 0.965,
                   0.02, 0.13, kNotUsed, 8.5, 0.005, 0.005 },
        // Tolerances = measured error x 1.1, rounded up; front times at least 0.002 (calibrated decks, README "Calibration on A3") (measured: NRMSE y_CO2 0.01381, y_CH4 0.2592, T RMSE 7.16 K, front time CO2 0.004063, CH4 0.002858).
        StrebCase{ "A7", 10e5, 338.15, 20e-6, false, 1100.0, 0.995,
                   0.02, 0.29, kNotUsed, 8, 0.005, 0.005 },
        // Tolerances = measured error x 1.1, rounded up; front times at least 0.002 (calibrated decks, README "Calibration on A3") (measured: NRMSE y_CO2 0.04515, y_CH4 0.06507, y_H2 0.2493, T RMSE 7.367 K, front time CO2 0.0009689, CH4 0.04696).
        StrebCase{ "B1", 1e5, 298, 20e-6, true, 10400.0, 1.057,
                   0.05, 0.08, 0.28, 8.5, 0.005, 0.055 },
        // Tolerances = measured error x 1.1, rounded up; front times at least 0.002 (calibrated decks, README "Calibration on A3") (measured: NRMSE y_CO2 0.04955, y_CH4 0.04676, y_H2 0.2165, T RMSE 9.126 K, front time CO2 0.0006391, CH4 0.01669).
        StrebCase{ "B2", 10e5, 298, 10e-6, true, 3350.0, 0.887,
                   0.06, 0.06, 0.24, 10.5, 0.005, 0.02 },
      };
      for (const auto& c : cases)
        if (c.name == name) return c;
      throw std::invalid_argument("no Streb 2021 case " + name);
    }

    // The same case with every tolerance infinite (recorded only).
    Streb2021Test::StrebCase RecordOnly(Streb2021Test::StrebCase c) {
      c.tol_nrmse_co2 = c.tol_nrmse_ch4 = c.tol_nrmse_h2 = kNotUsed;
      c.tol_temperature_rmse = c.tol_t50_co2 = c.tol_t50_ch4 = kNotUsed;
      return c;
    }
  }  // namespace

  // The IAST variant uses the pure Sips 2 form of the same Table 3 parameters,
  // so at y_i -> 1 it must give the extended-Sips loading exactly: the switch
  // changes only the competition. Checked over the three test temperatures and
  // 1 mbar - 25 bar for every component. The 50:50 CO2/CH4 loadings at
  // 10 bar, 25 C of both models are printed for the README.
  TEST_F(Streb2021Test, IastPureLimitsMatchExtendedSips) {
    Sips2Isotherm pure;
    SetTable3Sips(pure.params_, true);
    const std::vector<std::string> comps = StrebAdsorbing(true);
    std::vector<IIsotherm*> isos(comps.size(), &pure);
    for (double T : { 298.0, 318.15, 338.15 }) {
      for (double P = 1e2; P <= 25e5 * 1.0001; P *= std::sqrt(10.0)) {
        for (size_t i = 0; i < comps.size(); ++i) {
          std::vector<double> y(comps.size(), 0.0), q;
          y[i] = 1.0;
          const IastResult r = EvalIastMixture(isos, comps, T, P, y, q);
          ASSERT_TRUE(r.converged) << comps[i] << " T=" << T << " P=" << P;
          double q_ext = 0.0;
          IsothermFunction::ExtendedSips2(q_ext, T, P, 1.0, pure.params_.at(comps[i]));
          EXPECT_NEAR(q[i], q_ext, 1e-10 * std::max(1.0, q_ext)) << comps[i] << " T=" << T << " P=" << P;
          for (size_t j = 0; j < comps.size(); ++j)
            if (j != i) EXPECT_EQ(q[j], 0.0) << comps[j] << " absent, T=" << T << " P=" << P;
        }
      }
    }

    // 50:50 CO2/CH4 at 10 bar, 25 C: extended Sips (shared denominator) and IAST.
    const double T = 298.0, P = 10e5;
    const std::vector<std::string> bin = StrebAdsorbing(false);
    std::vector<IIsotherm*> bin_isos(bin.size(), &pure);
    std::vector<double> y = { 0.5, 0.5 }, q_iast;
    ASSERT_TRUE(EvalIastMixture(bin_isos, bin, T, P, y, q_iast).converged);
    const double R = constants::kIdealGasConstant;
    double z[2], qmax[2], zsum = 0.0;
    for (int i = 0; i < 2; ++i) {
      const auto& ip = pure.params_.at(bin[i]);
      const double b = ip.B1 * std::exp(ip.B2 / (R * T));
      const double s = ip.C1 + ip.C2 * (T / ip.T_ref - 1);
      z[i] = std::pow(b * P * y[i], s);
      qmax[i] = ip.A1 * std::exp(ip.A2 * (T / ip.T_ref - 1));
      zsum += z[i];
    }
    for (int i = 0; i < 2; ++i)
      std::cout << "[streb iast] 50:50, 10 bar, 298 K: " << bin[i] << " extended Sips "
                << qmax[i] * z[i] / (1 + zsum) << " mol/kg, IAST " << q_iast[i] << " mol/kg" << std::endl;
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
    const StrebCase c = StrebCaseFor("A3");
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
    const StrebCase c = StrebCaseFor("A1");
    SetUpStreb(c);
    RunAndValidate(c, true);
  }

  // The remaining binary (A2, A4-A7) and ternary (B1, B2) experiments of Table 2
  // for which digitised data exist.
  TEST_F(Streb2021Test, ExperimentA2) {
    data.info.title = "Streb 2021 - Experiment A2";
    const StrebCase c = StrebCaseFor("A2");
    SetUpStreb(c);
    RunAndValidate(c, false);
  }

  TEST_F(Streb2021Test, ExperimentA4) {
    data.info.title = "Streb 2021 - Experiment A4";
    const StrebCase c = StrebCaseFor("A4");
    SetUpStreb(c);
    RunAndValidate(c, false);
  }

  TEST_F(Streb2021Test, ExperimentA5) {
    data.info.title = "Streb 2021 - Experiment A5";
    const StrebCase c = StrebCaseFor("A5");
    SetUpStreb(c);
    RunAndValidate(c, false);
  }

  TEST_F(Streb2021Test, ExperimentA6) {
    data.info.title = "Streb 2021 - Experiment A6";
    const StrebCase c = StrebCaseFor("A6");
    SetUpStreb(c);
    RunAndValidate(c, false);
  }

  TEST_F(Streb2021Test, ExperimentA7) {
    data.info.title = "Streb 2021 - Experiment A7";
    const StrebCase c = StrebCaseFor("A7");
    SetUpStreb(c);
    RunAndValidate(c, false);
  }

  TEST_F(Streb2021Test, ExperimentB1) {
    data.info.title = "Streb 2021 - Experiment B1";
    const StrebCase c = StrebCaseFor("B1");
    SetUpStreb(c);
    RunAndValidate(c, false);
  }

  TEST_F(Streb2021Test, ExperimentB2) {
    data.info.title = "Streb 2021 - Experiment B2";
    const StrebCase c = StrebCaseFor("B2");
    SetUpStreb(c);
    RunAndValidate(c, false);
  }

  // IAST over the pure Sips terms of Table 3 (the competition model the paper
  // finds better for the binary breakthroughs, Sect. 4.2.1), every case
  // recorded without assertion. Extended Sips stays the decks' model: README
  // "IAST versus extended Sips" (2026-11-10) has the comparison.
  TEST_F(Streb2021Test, ExperimentA3_IAST_RecordOnly) {
    iast_ = true;
    const StrebCase c = RecordOnly(StrebCaseFor("A3"));
    SetUpStreb(c);
    RunAndValidate(c, false);
  }

  TEST_F(Streb2021Test, ExperimentA1_IAST_RecordOnly) {
    iast_ = true;
    const StrebCase c = RecordOnly(StrebCaseFor("A1"));
    SetUpStreb(c);
    RunAndValidate(c, false);
  }

  TEST_F(Streb2021Test, ExperimentA2_IAST_RecordOnly) {
    iast_ = true;
    const StrebCase c = RecordOnly(StrebCaseFor("A2"));
    SetUpStreb(c);
    RunAndValidate(c, false);
  }

  TEST_F(Streb2021Test, ExperimentA4_IAST_RecordOnly) {
    iast_ = true;
    const StrebCase c = RecordOnly(StrebCaseFor("A4"));
    SetUpStreb(c);
    RunAndValidate(c, false);
  }

  TEST_F(Streb2021Test, ExperimentA5_IAST_RecordOnly) {
    iast_ = true;
    const StrebCase c = RecordOnly(StrebCaseFor("A5"));
    SetUpStreb(c);
    RunAndValidate(c, false);
  }

  TEST_F(Streb2021Test, ExperimentA6_IAST_RecordOnly) {
    iast_ = true;
    const StrebCase c = RecordOnly(StrebCaseFor("A6"));
    SetUpStreb(c);
    RunAndValidate(c, false);
  }

  TEST_F(Streb2021Test, ExperimentA7_IAST_RecordOnly) {
    iast_ = true;
    const StrebCase c = RecordOnly(StrebCaseFor("A7"));
    SetUpStreb(c);
    RunAndValidate(c, false);
  }

  TEST_F(Streb2021Test, ExperimentB1_IAST_RecordOnly) {
    iast_ = true;
    const StrebCase c = RecordOnly(StrebCaseFor("B1"));
    SetUpStreb(c);
    RunAndValidate(c, false);
  }

  TEST_F(Streb2021Test, ExperimentB2_IAST_RecordOnly) {
    iast_ = true;
    const StrebCase c = RecordOnly(StrebCaseFor("B2"));
    SetUpStreb(c);
    RunAndValidate(c, false);
  }

} // namespace adsim
