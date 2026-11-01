#pragma once

// ============================================================================
//  Perez et al. (2019) — CO2/N2 vacuum swing adsorption on Zeolite 13X
//  See README.md in this folder for full citation, DOI and reproduced
//  configurations. Paper / SI live under ./assets/.
//
//  Test fixture + shared bed/isotherm setup for the Perez 2019 validation
//  cases. Kept in its own header so the setup is co-located with the tests
//  that use it (moved out of the shared user_data_test.h god-fixture).
// ============================================================================

#define _USE_MATH_DEFINES
#include <cmath>

#include "../../fixtures/user_data_test.h"

namespace adsim {

  class Perez2019Test : public UserDataTest {
  public:
    // Parameterised cycle runners (defined in perez_2019.cpp) so each Table 3
    // experiment is a one-line TEST_F. Run to CSS and compare CO2 purity/recovery.
    void RunPerez4Step(const std::string& save_name, const std::string& title,
                       double feed_time, double interstitial_feed_vel,
                       double bd_pressure, double evac_pressure,
                       double blowdown_time, double evac_time, double press_time,
                       double target_purity, double target_recovery);
    void RunPerezLPP(const std::string& save_name, const std::string& title,
                     double feed_time, double interstitial_feed_vel,
                     double bd_pressure, double evac_pressure,
                     double blowdown_time, double evac_time, double press_time,
                     double target_purity, double target_recovery);

    void SetUpPerez2019() {
      const std::string CO2 = "CarbonDioxide";
      const std::string N2 = "Nitrogen";
      const std::string z13X = "13X";
   
      AddComponents({ CO2, N2 });
      auto [layer_id, ads_id] = AddLayer(z13X);

      double adsorbent_mass = 0.162; // [kg]
      double length = 0.3;
      double Di = 2 * 0.016;
      double volume = (M_PI / 4.0) * Di * Di * length;
      double bulk_density = adsorbent_mass / volume;
      double eb = 0.41;
      double solid_density = bulk_density / (1 - eb);

      auto& _13X = GetLayer(layer_id);
      _13X.domain_.length = length;
      _13X.particle_diameter_ = 2 * 7.5e-4;
      _13X.eb_ = eb;
      _13X.domain_.n_cells = 30;
      // SI Eq. S5 is the Blake-Kozeny law (constant 150); the deck keeps Darcy /
      // Kozeny-Carman (180): the paper calls the pressure drop negligible and
      // the difference is 0.3 pp of recovery (README). ConfigD1_ConnectedBC
      // runs on Ergun (see there).
      std::get<CylinderBed>(data.bed.shape_).diameter = Di;

      auto& _13X_ads = GetAdsorbent(ads_id);
      _13X_ads.solid_density_ = solid_density;
      _13X_ads.heat_capacity_ = 900;
      _13X_ads.thermal_conductivity_ = 0.63;
      _13X_ads.ep_ = 0.35;
      data.bed.wall_.heat_capacity_ = 502;
      data.bed.wall_.density_ = 7800;
      data.bed.wall_.thermal_conductivity_ = 16;
      data.bed.wall_.thickness_ = 3.1e-3;
      data.bed.ambient_temperature_ = 293;
      // SI Table S2: T_w = T_a at both ends of the wall in every step.
      auto wall_bc = std::make_unique<WallTemperatureBoundaryConditionFixed>();
      wall_bc->temperature_ = 293.0;
      data.bed.equations_.wall_temperature_boundary_condition_ = std::move(wall_bc);

      data.fluid.adsorbed_heat_capacity_[CO2] = 30.7;
      data.fluid.adsorbed_heat_capacity_[N2] = 30.7;

      // Competitive (extended) DSL. Perez et al. (2019), Sect. 4.1: "the description
      // of the competitive behavior of a given gas mixture is critical" - the
      // Wilkins & Rajendran approach (N2 with equal site energies, Table 2) used to
      // predict the competitive loadings. This deck previously used the
      // non-competitive DualSiteLangmuirIsotherm_Concentration.
      auto dsl = std::make_unique<ExtendedDualSiteLangmuirIsotherm_Concentration>(CO2);
      dsl->params_[CO2].m1 = 2.768;
      dsl->params_[CO2].b1 = 3.934e-7;
      dsl->params_[CO2].q1 = 40.98e3;
      dsl->params_[CO2].m2 = 2.758;
      dsl->params_[CO2].b2 = 1.256e-7;
      dsl->params_[CO2].q2 = 34.19e3;
      dsl->params_[N2].m1 = 2.768;
      dsl->params_[N2].b1 = 3.134e-6;
      dsl->params_[N2].q1 = 15.71e3;
      dsl->params_[N2].m2 = 2.758;
      dsl->params_[N2].b2 = 3.134e-6;
      dsl->params_[N2].q2 = 15.72e3;
      AddIsotherm({ CO2, N2 }, ads_id, std::move(dsl));

      // SI Eq. S4: the macropore coefficient times the capacity factor c / (rho_p q*).
      AddLinearDrivingForceModel_MacroporeCapacity(ads_id, CO2, 3.0, 0.35);
      AddLinearDrivingForceModel_MacroporeCapacity(ads_id, N2, 3.0, 0.35);

      SetHOA(ads_id, CO2, 30.85e3);
      SetHOA(ads_id, N2, 15.93e3);

      AddConstantInternalHeatTransferCoefficientModel(8.6, layer_id);
      AddConstantExternalHeatTransferCoefficientModel(10, layer_id);

      auto& layer = GetLayer(layer_id);
      layer.equations_.dispersion_model_ = CreateAxialDispersionModel(AxialDispersionType::RUTHVEN);

      // Effective gas thermal conductivity K_z = 0.0903 J/(m K s) (Table 1). This
      // used to be left at the constant model's default of 0.
      auto eff_thermal_conductivity = std::make_unique<HeatDispersionModelConstant>();
      eff_thermal_conductivity->dispersion_coefficient_ = 0.0903; // [W/(m-K)]
      layer.equations_.heat_dispersion_model_ = std::move(eff_thermal_conductivity);

      AddInitialConditions(
        1.0e5, 298, 298,
        {
          {"CarbonDioxide", 0.0 },
          {"Nitrogen", 1.0 }
        }
      );
    }
  };

} // namespace adsim
