#pragma once

// ============================================================================
//  Nguyen et al. (2022) — CO2 breakthrough on CALF-20
//  See README.md in this folder for citation / DOI. Paper + SI under ./assets/.
// ============================================================================

#define _USE_MATH_DEFINES
#include <cmath>

#include "../../fixtures/user_data_test.h"

namespace adsim {

  class Nguyen2022Test : public UserDataTest {
  public:
    void SetUpNguyen2022() {
      const std::string CO2 = "CarbonDioxide";
      const std::string N2 = "Nitrogen";
      const std::string He = "Helium";
      const std::string calf20 = "CALF-20";
   
      AddComponents({ CO2, N2, He });
      auto [layer_id, ads_id] = AddLayer(calf20);

      double adsorbent_mass = 16.7 / 1000.0; // [kg]
      double Di = 0.0282;
      double Do = 0.0318;
      double tw = (Do - Di) / 2;
      double area = (M_PI / 4.0) * Di * Di;
      double solid_density = 570.0;
      // double bulk_density = adsorbent_mass / volume; // ~389 kg/m3?
      double eb = 0.4;
      double ep = 0.35;
      double rho_b = solid_density * (1 - eb);
      double volume = adsorbent_mass / rho_b;
      double length = volume / area;
      // double volume = area * length;
      // double solid_density = bulk_density / (1 - eb);

      // Column length is derived from the adsorbent mass and rho_b = rho_p (1 - eb)
      // so that the deck holds exactly 16.7 g of CALF-20.

      auto& _calf20 = GetLayer(layer_id);
      _calf20.domain_.length = length;
      // Particle diameter: not stated in the material available for this paper
      // (its text was never added to assets/). 1.5 mm is the deck author's
      // assumption; it enters only through the macropore LDF coefficient and the
      // Darcy/Carman-Kozeny pressure drop. (The follow-up paper, nguyen_2023,
      // reports ~1 mm pellets for structured CALF-20.)
      _calf20.particle_diameter_ = 2.0 * 7.5e-4;
      _calf20.eb_ = eb;
      _calf20.domain_.n_cells = 30;
      std::get<CylinderBed>(data.bed.shape_).diameter = Di;

      auto& _calf20_ads = GetAdsorbent(ads_id);
      _calf20_ads.solid_density_ = solid_density;
      _calf20_ads.heat_capacity_ = 1371;
      _calf20_ads.thermal_conductivity_ = 0.63;
      _calf20_ads.ep_ = ep;
      data.bed.wall_.heat_capacity_ = 502.0;
      data.bed.wall_.density_ = 7800.0;
      data.bed.wall_.thermal_conductivity_ = 16.0;
      data.bed.wall_.thickness_ = tw;
      // Ambient / initial temperature: not stated in the available material; 295 K
      // matches the initial thermocouple readings of the digitised experiments
      // (core/data/nguyen_2022/Fig*/T_CO2_*.csv start at 294.8-296.7 K).
      data.bed.ambient_temperature_ = 295.0;

      data.bed.equations_.wall_temperature_boundary_condition_ = std::make_unique<WallTemperatureBoundaryConditionZeroFlux>();

      data.fluid.adsorbed_heat_capacity_[CO2] = 37.435; // J/(mol-K)
      data.fluid.adsorbed_heat_capacity_[N2] = 29.171; // J/(mol-K)

      // Parameters from ESS model in Table 3
      // Separation of CO2 and N2 on a hydrophobic metal organic framework CALF-20
      auto dsl = std::make_unique<ExtendedDualSiteLangmuirIsotherm_Concentration>(CO2);
      dsl->params_[CO2].m1 = 2.387; // mol/kg
      dsl->params_[CO2].b1 = 5.519e-7; // m3/mol
      dsl->params_[CO2].q1 = 35.06e3; // J/mol
      dsl->params_[CO2].m2 = 3.271; // mol/kg
      dsl->params_[CO2].b2 = 5.187e-8; // m3/mol
      dsl->params_[CO2].q2 = 28.95e3; // J/mol

      dsl->params_[N2].m1 = 2.387; // mol/kg
      dsl->params_[N2].b1 = 8.139e-7; // m3/mol
      dsl->params_[N2].q1 = 17.96e3; // J/mol
      dsl->params_[N2].m2 = 3.271; // mol/kg
      dsl->params_[N2].b2 = 8.139e-7; // m3/mol
      dsl->params_[N2].q2 = 17.96e3; // J/mol
      AddIsotherm({ CO2, N2 }, ads_id, std::move(dsl));

      AddLinearDrivingForceModel_Macropore(ads_id, CO2, 3, 0.35);
      AddLinearDrivingForceModel_Macropore(ads_id, N2, 3, 0.35);

      SetHOA(ads_id, CO2, 38.0e3);
      SetHOA(ads_id, N2, 15.0e3);

      AddConstantInternalHeatTransferCoefficientModel(6.963, layer_id);
      AddConstantExternalHeatTransferCoefficientModel(14.0, layer_id);

      _calf20.equations_.dispersion_model_ = CreateAxialDispersionModel(AxialDispersionType::RUTHVEN);

      auto eff_thermal_conductivity = std::make_unique<HeatDispersionModelConstant>();
      eff_thermal_conductivity->dispersion_coefficient_ = 0.4; // [W/(m-K)]
      _calf20.equations_.heat_dispersion_model_ = std::move(eff_thermal_conductivity);
    }
  };

} // namespace adsim
