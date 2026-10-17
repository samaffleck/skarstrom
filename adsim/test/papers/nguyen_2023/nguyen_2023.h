#pragma once

// ============================================================================
//  Nguyen et al. (2023) — CO2/N2 VPSA cycles on CALF-20
//  DOI: 10.1016/j.cej.2022.139550
//  See README.md in this folder. Paper + SI under ./assets/.
// ============================================================================

#define _USE_MATH_DEFINES
#include <cmath>

#include "../../fixtures/user_data_test.h"

namespace adsim {

  class Nguyen2023Test : public UserDataTest {
  public:
    void SetUpNguyen2023() {
      const std::string CO2 = "CarbonDioxide";
      const std::string N2 = "Nitrogen";
      const std::string calf20 = "CALF-20";
   
      AddComponents({ CO2, N2 });
      auto [layer_id, ads_id] = AddLayer(calf20);

      double adsorbent_mass = 36.0 / 1000.0; // [kg]
      // double length = 0.267; // [m]
      double Di = 0.021;
      double Do = 0.025;
      double tw = (Do - Di) / 2;
      double area = (M_PI / 4.0) * Di * Di;
      double solid_density = 570.0;
      double eb = 0.4;
      double ep = 0.35;
      double rho_b = solid_density * (1 - eb);
      double volume = adsorbent_mass / rho_b;
      double length = volume / area;
      // double volume = area * length;
      // double solid_density = bulk_density / (1 - eb);

      // Table 1 is internally inconsistent: m_ads = 36.0 g, rho_p = 570 kg/m3 and
      // eps = 0.4 in a 0.021 m x 0.267 m column would need rho_p ~ 649 kg/m3.
      // The deck keeps the measured mass and particle density and derives the
      // length from them (L ~ 0.304 m instead of the measured 0.267 m), so the
      // adsorbent inventory - which sets the CO2 capacity - matches the paper.

      auto& _calf20 = GetLayer(layer_id);
      _calf20.domain_.length = length;
      _calf20.particle_diameter_ = 1.0e-3; // ~1 mm pellet diameter (paper Section 2); r_p = 0.5 mm
      _calf20.eb_ = eb;
      _calf20.domain_.n_cells = 100;
      std::get<CylinderBed>(data.bed.shape_).diameter = Di;

      auto& _calf20_ads = GetAdsorbent(ads_id);
      _calf20_ads.solid_density_ = solid_density;
      _calf20_ads.heat_capacity_ = 1350;
      _calf20_ads.thermal_conductivity_ = 0.63;
      _calf20_ads.ep_ = ep;
      data.bed.wall_.heat_capacity_ = 502.0;
      data.bed.wall_.density_ = 7800.0;
      data.bed.wall_.thermal_conductivity_ = 10.0;
      data.bed.wall_.thickness_ = tw;
      // Ambient temperature: not stated in the paper or SI text (Table 1 lists
      // none). 293 K is an assumption; the isotherms were measured from 22 C up.
      data.bed.ambient_temperature_ = 293.0;

      data.fluid.adsorbed_heat_capacity_[CO2] = 37.435; // J/(mol-K)
      data.fluid.adsorbed_heat_capacity_[N2] = 29.171; // J/(mol-K)

      // Competitive (extended) DSL. Nguyen et al. (2023), Sect. 2: "The
      // single-component isotherms and the competitive forms were shown to be
      // well-described by the dual-site Langmuir isotherm [10]" (N2 with the equal
      // energy scheme, SI Table S1) - the same extended form the nguyen_2022 deck
      // uses. This deck previously used the non-competitive
      // DualSiteLangmuirIsotherm_Concentration.
      auto dsl = std::make_unique<ExtendedDualSiteLangmuirIsotherm_Concentration>(CO2);
      dsl->params_[CO2].m1 = 2.687; // mol/kg
      dsl->params_[CO2].b1 = 5.742e-7; // m3/mol
      dsl->params_[CO2].q1 = 35.18e3; // J/mol
      dsl->params_[CO2].m2 = 4.597; // mol/kg
      dsl->params_[CO2].b2 = 3.684e-9; // m3/mol
      dsl->params_[CO2].q2 = 35.03e3; // J/mol

      dsl->params_[N2].m1 = 2.687; // mol/kg
      dsl->params_[N2].b1 = 7.466e-7; // m3/mol
      dsl->params_[N2].q1 = 17.91e3; // J/mol
      dsl->params_[N2].m2 = 4.597; // mol/kg
      dsl->params_[N2].b2 = 7.466e-7; // m3/mol
      dsl->params_[N2].q2 = 17.91e3; // J/mol
      AddIsotherm({ CO2, N2 }, ads_id, std::move(dsl));

      AddLinearDrivingForceModel_Macropore(ads_id, CO2, 3, 0.35);
      AddLinearDrivingForceModel_Macropore(ads_id, N2, 3, 0.35);

      SetHOA(ads_id, CO2, 38.0e3);
      SetHOA(ads_id, N2, 15.0e3);

      AddConstantInternalHeatTransferCoefficientModel(8.0, layer_id);
      AddConstantExternalHeatTransferCoefficientModel(10.0, layer_id);

      _calf20.equations_.dispersion_model_ = CreateAxialDispersionModel(AxialDispersionType::RUTHVEN);

      auto eff_thermal_conductivity = std::make_unique<HeatDispersionModelConstant>();
      eff_thermal_conductivity->dispersion_coefficient_ = 0.4; // [W/(m-K)]
      _calf20.equations_.heat_dispersion_model_ = std::move(eff_thermal_conductivity);

      AddInitialConditions(
        1.0e5, 293, 293,
        {
          {"CarbonDioxide", 0.0 },
          {"Nitrogen", 1.0 }
        }
      );
    }
  };

} // namespace adsim
