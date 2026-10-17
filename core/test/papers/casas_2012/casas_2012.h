#pragma once

// ============================================================================
//  Casas et al. (2012) — CO2/H2 breakthrough on AP3-60 activated carbon
//  See README.md in this folder for citations / DOIs. Papers + markdown under
//  ./assets/. Isotherm equilibria come from the companion equilibrium study
//  (Schell et al. 2012); the fixed-bed breakthrough setup comes from Casas 2012.
// ============================================================================
//
//  Isotherm: competitive Sips (Langmuir-Freundlich) with the paper's EXACT
//  temperature-dependence, implemented by adsim's ArrheniusSipsIsotherm:
//
//      q_i   = n∞_i · z_i / (1 + Σ_j z_j),   z_i = (k_i · P · y_i)^{c_i}
//      n∞_i(T) = a_i · exp(b_i / (R·T))                        [Arrhenius]
//      k_i(T)  = A_i · exp(B_i / (R·T))                        [Arrhenius]
//      c_i(T)  = c_ref,i + amp_i·(2/π)·arctan(β_i·(T − Tref,i)) [S-shaped]
//
//  These are Eqs. (5)–(7) of Schell et al. (2012) plus the inverse-tangent
//  exponent for CO2; H2 uses constant n∞ (b=0) and constant c (amp=0). The
//  per-component constants below were obtained by fitting these exact forms to
//  the paper's per-temperature Table 3 values across all five measured
//  temperatures (25–140 °C); the fit recovers the paper's own Table 4 Arrhenius
//  magnitudes to 3–4 significant figures and reproduces Table 3 loadings to
//  within the paper's own Table-3→Table-4 smoothing (≤7 % for CO2, <1 % for H2).
//  See `casas_2012_isotherm.cpp` for the acceptance test and README.md for the
//  full derivation.
// ============================================================================

#define _USE_MATH_DEFINES
#include <cmath>
#include <limits>
#include <string>

#include "../../fixtures/user_data_test.h"
#include "../../fixtures/validation_metrics.h"
#include "adsim/equations/isotherm.h"
#include "adsim/equations/dispersion_model.h"
#include "adsim/equations/heat_transfer_model.h"
#include "adsim/boundary_conditions/boundary_conditions.h"
#include "adsim/boundary_conditions/boundary_stream.h"
#include "adsim/boundary_conditions/step.h"
#include "adsim/solver/coupled_solver.h"
#include "adsim/solver/kpis.h"

#include <gtest/gtest.h>

namespace adsim {

  // Single source of truth for the Schell/Casas (2012) AP3-60 isotherm constants.
  // Obtained by fitting the paper's exact functional forms (Arrhenius n∞ & k,
  // arctan c for CO2) to the per-temperature Sips parameters of Table 3 across
  // 25–140 °C. Used by both the breakthrough fixture and the isotherm acceptance
  // test (`casas_2012_isotherm.cpp`).
  inline void SetSchellSipsParams(ArrheniusSipsIsothermParameters& co2,
                                  ArrheniusSipsIsothermParameters& h2) {
    // CO2: temperature-dependent n∞, k and (S-shaped) exponent c.
    co2.n_inf_a = 1.38286;      // a      [mol/kg]
    co2.n_inf_b = 5625.93;      // b      [J/mol]
    co2.k_A     = 1.6814e-08;   // A      [1/Pa]
    co2.k_B     = 9160.22;      // B      [J/mol]
    co2.c_ref   = 0.832365;     // c_ref  [-]
    co2.c_amp   = 0.116571;     // c_amp  [-]
    co2.c_beta  = 0.0928714;    // c_beta [1/K]
    co2.c_Tref  = 330.129;      // c_Tref [K]

    // H2: temperature-independent n∞ (b=0) and exponent c (amp=0).
    h2.n_inf_a = 6.66;          // [mol/kg]
    h2.n_inf_b = 0.0;           // [J/mol]
    h2.k_A     = 6.6886e-10;    // [1/Pa]
    h2.k_B     = 9954.73;       // [J/mol]
    h2.c_ref   = 0.96;          // [-]
    h2.c_amp   = 0.0;           // [-]
    h2.c_beta  = 0.0;           // [1/K]
    h2.c_Tref  = 298.15;        // [K]
  }

  class Casas2012Test : public UserDataTest {
  public:
    // Physical setup shared by every Casas 2012 breakthrough case.
    // feed_temperature_K sets the (isothermal) column-wall / ambient reference;
    // the isotherm itself is temperature dependent, so the same parameter set
    // serves both the 25 °C and 45 °C experiments.
    void SetUpCasas2012(double feed_temperature_K) {
      const std::string CO2 = "CarbonDioxide";
      const std::string H2 = "Hydrogen";
      const std::string He = "Helium";
      const std::string carbon = "AP3-60";

      AddComponents({ CO2, H2, He });
      auto [layer_id, ads_id] = AddLayer(carbon);

      // --- Column geometry & packing (Casas 2012, Table 2) ---
      const double length = 1.20;              // [m]
      const double inner_diameter = 0.025;     // [m] (2.5 cm)
      const double eb = 0.403;                 // measured bed void fraction [-]
      const double particle_density = 850.0;   // [kg/m3] (manufacturer)
      const double skeletal_density = 1970.0;  // [kg/m3] (He pycnometer, 1.97 g/cm3)
      const double ep = 1.0 - particle_density / skeletal_density; // ≈ 0.5685
      const double particle_diameter = 0.003;  // [m]

      auto& carbon_layer = GetLayer(layer_id);
      carbon_layer.domain_.length = length;
      carbon_layer.particle_diameter_ = particle_diameter;
      carbon_layer.eb_ = eb;
      carbon_layer.domain_.n_cells = 60;
      std::get<CylinderBed>(data.bed.shape_).diameter = inner_diameter;

      auto& carbon_ads = GetAdsorbent(ads_id);
      carbon_ads.solid_density_ = particle_density;
      carbon_ads.heat_capacity_ = 1000.0;      // C_s [J/(kg-K)] (Table 2)
      carbon_ads.thermal_conductivity_ = 0.15; // not tabulated; nominal (see README)
      carbon_ads.ep_ = ep;

      // --- Column wall (lumped heating shells, Casas 2012 Table 2) ---
      // Reported as a lumped volumetric heat capacity C_w = 4e6 J/(m3-K) and a
      // lumped outer radius R_o = 0.02 m (wall + copper/brass shells + heater).
      data.bed.wall_.heat_capacity_ = 500.0;   // J/(kg-K)  }
      data.bed.wall_.density_ = 8000.0;        // kg/m3     } product = 4e6 J/(m3-K)
      data.bed.wall_.thermal_conductivity_ = 0.0; // axial wall conduction neglected
      data.bed.wall_.thickness_ = 0.02 - 0.5 * inner_diameter; // R_o - R_i = 0.0075 m
      data.bed.equations_.wall_temperature_boundary_condition_ =
          CreateWallTemperatureBoundaryCondition(WallTemperatureBoundaryConditionType::ZERO_FLUX);
      data.bed.ambient_temperature_ = feed_temperature_K;

      // --- Competitive Arrhenius Sips isotherm (Schell 2012 exact form) ---
      auto sips = std::make_unique<ArrheniusSipsIsotherm>();
      SetSchellSipsParams(sips->params_[CO2], sips->params_[H2]);
      AddIsotherm({ CO2, H2 }, ads_id, std::move(sips));

      // --- Mass transfer: lumped LDF, fitted coefficients (Casas 2012 Table 3) ---
      AddLinearDrivingForceModel(ads_id, CO2, 0.15); // k_CO2 [1/s]
      AddLinearDrivingForceModel(ads_id, H2, 1.0);   // k_H2  [1/s]

      // --- Isosteric heats of adsorption (Casas 2012 Table 2, Clausius-Clapeyron) ---
      SetHOA(ads_id, CO2, 26.0e3); // [J/mol]
      SetHOA(ads_id, H2, 9.8e3);   // [J/mol]

      // --- Dispersion & heat transfer ---
      carbon_layer.equations_.dispersion_model_ =
          CreateAxialDispersionModel(AxialDispersionType::RUTHVEN);
      carbon_layer.equations_.heat_dispersion_model_ =
          CreateAxialHeatDispersionModel(AxialHeatDispersionType::CONSTANT);
      // Internal (gas→wall) HTC h_L: Leva (1947) correlation with the
      // coefficient/exponent Casas fitted to the reference experiment.
      auto leva = std::make_unique<InternalHeatTransferModel_Leva>();
      leva->eta1 = 41.13;  // fitted (Casas 2012, Table 3; Leva canonical 0.813)
      leva->eta2 = 0.32;   // fitted (Casas 2012, Table 3; Leva canonical 0.90)
      carbon_layer.equations_.internal_heat_transfer_model_ = std::move(leva);
      // Wall→ambient HTC h_W: fitted constant (Casas 2012, Table 3).
      AddConstantExternalHeatTransferCoefficientModel(5.0, layer_id); // h_W [W/(m2-K)]
    }

    // Runs a single CO2/H2 breakthrough: the bed is pre-pressurised with inert
    // Helium at `total_pressure_Pa`, then fed a CO2/H2 mixture at constant
    // superficial volumetric flow (10 cm3/s, Casas 2012) held at the outlet
    // pressure. `y_CO2` is the CO2 mole fraction of the dry CO2/H2 feed.
    //
    // Validation: the H2 and CO2 breakthrough times (first time the outlet mole
    // fraction exceeds 0.002, the paper's Sect. 4.4 definition) are compared with
    // Table 5 of Casas et al. (2012) as a relative error. Table 5 lists the
    // paper's SIMULATED times; the experimental curves (Figs. 4, 6) have not been
    // digitised. The paper's model also includes the column-to-MS piping and a
    // per-experiment fitted feed velocity (Fig. 10), neither of which adsim has.
    void RunBreakthroughCase(double total_pressure_Pa, double feed_temperature_K,
                             double y_CO2, double duration_s,
                             const std::string& save_label,
                             const std::string& case_name = "",
                             double tol_h2_rel = std::numeric_limits<double>::infinity(),
                             double tol_co2_rel = std::numeric_limits<double>::infinity()) {
      const std::string CO2 = "CarbonDioxide";
      const std::string H2 = "Hydrogen";
      const std::string He = "Helium";
      const double y_H2 = 1.0 - y_CO2;

      SetUpCasas2012(feed_temperature_K);

      AddInitialConditions(
          total_pressure_Pa, feed_temperature_K, feed_temperature_K,
          { { CO2, 0.0 }, { H2, 0.0 }, { He, 1.0 } });

      AddStep("Feed", duration_s);

      // Feed velocity is held constant across all experiments at 10 cm3/s
      // (superficial). 1 cm3/s = 1e-6 m3/s.
      const double volumetric_flowrate_m3_s = 10.0e-6;

      Step& step_feed = data.cycle.GetStep("Feed");
      auto feed_inlet_stream = CreateNoValveInletStream(
          CreateVolumeFlowBC(volumetric_flowrate_m3_s),
          { { CO2, y_CO2 }, { H2, y_H2 }, { He, 0.0 } },
          feed_temperature_K);
      feed_inlet_stream->name_ = "CO2/H2 Feed";

      auto feed_flow_bc = std::make_unique<FlowBoundaryCondition>();
      feed_flow_bc->direction_ = BoundaryDirection::INLET;
      feed_flow_bc->location_ = BoundaryLocation::BOTTOM;
      feed_flow_bc->streams_.push_back(std::move(feed_inlet_stream));
      step_feed.bottom_boundary_condition_ = std::move(feed_flow_bc);

      auto feed_pressure_out = std::make_unique<PressureBoundaryCondition>();
      feed_pressure_out->direction_ = BoundaryDirection::OUTLET;
      feed_pressure_out->location_ = BoundaryLocation::TOP;
      feed_pressure_out->pressure_ = CreateConstantBoundaryValue(total_pressure_Pa);
      feed_pressure_out->mole_fraction_[CO2] = 0.0;
      feed_pressure_out->mole_fraction_[H2] = 0.0;
      feed_pressure_out->mole_fraction_[He] = 1.0;
      step_feed.top_boundary_condition_ = std::move(feed_pressure_out);

      ConfigureCycle(1);

      CoupledSolver solver{ data };

      std::string error{};
      EXPECT_TRUE(data.PerformChecks(error)) << "System checks failed: " << error;

      bool verbose = false;
      int retval = solver.IntegrateCycle(data.cycle, verbose);
      EXPECT_GE(retval, 0) << "Solver did not complete";

      EXPECT_TRUE(IsMassConserved(solver.GetLog())) << "Mass is not conserved.";
      EXPECT_TRUE(IsMolesConserved(solver.GetLog())) << "Moles are not conserved.";
      EXPECT_TRUE(IsComponentMolesConserved(solver.GetLog())) << "Component moles not conserved.";

      SaveTest(save_label, solver);
      SaveResultsAsExample(save_label, solver);

      if (!case_name.empty()) {
        using namespace validation;
        CaseReport report("casas_2012", case_name);
        const double t_h2 = CrossingTime(FaceMoleFraction(solver, H2, 0, 0), 0.002);
        const double t_co2 = CrossingTime(FaceMoleFraction(solver, CO2, 0, 0), 0.002);
        report.CheckRel("H2 breakthrough time (y > 0.002)", t_h2,
                        Reference("casas_2012", case_name, "H2_breakthrough_time"), tol_h2_rel);
        report.CheckRel("CO2 breakthrough time (y > 0.002)", t_co2,
                        Reference("casas_2012", case_name, "CO2_breakthrough_time"), tol_co2_rel);
        if (case_name == "Fig4_Reference_50_50_15bar_25C")
          report.Info("H2 breakthrough time vs experiment (approx.)", t_h2,
                      Reference("casas_2012", case_name, "H2_breakthrough_time_experiment"), "rel_error",
                      "text value is approximate ('about 110 s')");
      }
    }
  };

} // namespace adsim
