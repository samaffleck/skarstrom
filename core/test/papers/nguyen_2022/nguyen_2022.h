#pragma once

// ============================================================================
//  Nguyen et al. (2022) — CO2 breakthrough on CALF-20
//  See README.md in this folder for citation / DOI and the provenance table.
// ============================================================================

#define _USE_MATH_DEFINES
#include <cmath>
#include <string>

#include "../../fixtures/user_data_test.h"

namespace adsim {

  // One dynamic column breakthrough experiment of the paper's Figs. 4 and 5.
  // Flow, temperature and pressure: Table A2.1 of the author's thesis (T.T.T.
  // Nguyen, PhD thesis, University of Alberta), which lists every experiment.
  struct NguyenCase {
    std::string name;      // test case name, e.g. "Fig4_70pct"
    std::string fig;       // "Fig4" (CO2/He) or "Fig5" (CO2/N2)
    int percent;           // feed CO2 [mol %]
    double flow_ccm;       // feed flow [cm3/min], taken at feed conditions
    double temperature_c;  // feed, initial and ambient temperature [C]
    double pressure_bar;   // column pressure [bar]
    double duration;       // feed time [s], covers the figure
    std::string deck;      // example deck name
  };

  class Nguyen2022Fixture : public UserDataTest {
  public:
    // Column, CALF-20 and transport (Table 2 of the paper unless noted); the
    // feed and the initial state of case `c`.
    void SetUpNguyen2022(const NguyenCase& c) {
      const std::string CO2 = "CarbonDioxide";
      const std::string N2 = "Nitrogen";
      const std::string He = "Helium";

      AddComponents({ CO2, N2, He });
      auto [layer_id, ads_id] = AddLayer("CALF-20");

      const double T = c.temperature_c + 273.15;
      const double P = c.pressure_bar * 1e5;

      // Constant gas heat capacities of Table 2 (NIST; per kg: CO2 841, N2 1040,
      // He 5200 J/(kg K)) per mole, and the adsorbed phase equal to the gas
      // ("assumed", Table 2). Helium does not adsorb.
      data.fluid.gas_heat_capacity_mode_ = GasHeatCapacityMode::CONSTANT;
      data.fluid.gas_heat_capacity_[CO2] = kCpCO2;
      data.fluid.gas_heat_capacity_[N2] = kCpN2;
      data.fluid.gas_heat_capacity_[He] = kCpHe;
      data.fluid.adsorbed_heat_capacity_[CO2] = kCpCO2;
      data.fluid.adsorbed_heat_capacity_[N2] = kCpN2;
      data.fluid.adsorbed_heat_capacity_[He] = 0.0;

      const double Di = 0.0282, Do = 0.0318;
      const double area = 0.25 * M_PI * Di * Di;
      auto& layer = GetLayer(layer_id);
      // L = m / (rho_p (1 - eb) A): the bed holds exactly the 16.7 g of Table 2
      // (0.0782 m; Table 2 gives 0.0786 m, 0.5 % longer, README).
      layer.domain_.length = kMass / (kParticleDensity * (1 - kBedVoidage) * area);
      layer.particle_diameter_ = kParticleDiameter;
      layer.eb_ = kBedVoidage;
      layer.domain_.n_cells = 30;                              // Sect. 4.1: 30 cells
      layer.momentum_model_ = MomentumModel::BlakeKozeny;      // Eq. 7 (Darcy, 150)
      std::get<CylinderBed>(data.bed.shape_).diameter = Di;

      auto& calf20 = GetAdsorbent(ads_id);
      calf20.solid_density_ = kParticleDensity;
      calf20.heat_capacity_ = 1371;                            // Table 2, fitted
      calf20.thermal_conductivity_ = 0.63;                     // not in the paper; unused (no solid conduction in a packed bed)
      // The paper's balances (Eqs. 1-3) hold up gas in the bed voidage only;
      // the particle voidage enters the LDF coefficient alone (Eq. 5).
      calf20.ep_ = 0.0;

      data.bed.wall_.heat_capacity_ = 502.0;                   // Table 2
      data.bed.wall_.density_ = 7800.0;
      data.bed.wall_.thermal_conductivity_ = 16.0;
      data.bed.wall_.thickness_ = 0.5 * (Do - Di);
      // The paper does not state the wall-end condition (README).
      data.bed.equations_.wall_temperature_boundary_condition_ = std::make_unique<WallTemperatureBoundaryConditionZeroFlux>();
      data.bed.ambient_temperature_ = T;

      // Extended dual-site Langmuir, concentration form, equal-energy-site N2
      // (Table 3; Eq. 20).
      auto dsl = std::make_unique<ExtendedDualSiteLangmuirIsotherm_Concentration>(CO2);
      dsl->params_[CO2].m1 = 2.387;   dsl->params_[CO2].b1 = 5.519e-7; dsl->params_[CO2].q1 = 35.06e3;
      dsl->params_[CO2].m2 = 3.271;   dsl->params_[CO2].b2 = 5.187e-8; dsl->params_[CO2].q2 = 28.95e3;
      dsl->params_[N2].m1 = 2.387;    dsl->params_[N2].b1 = 8.139e-7;  dsl->params_[N2].q1 = 17.96e3;
      dsl->params_[N2].m2 = 3.271;    dsl->params_[N2].b2 = 8.139e-7;  dsl->params_[N2].q2 = 17.96e3;
      AddIsotherm({ CO2, N2 }, ads_id, std::move(dsl));

      // LDF with the macropore coefficient and the capacity factor (Eq. 5:
      // k = 15 eps_p D_p c / (r_p^2 q*), D_p = D_m / tau), D_m of Table 2.
      for (const std::string& comp : { CO2, N2 }) {
        auto ldf = std::make_unique<LinearDrivingForce>();
        auto mtc = std::make_unique<MassTransferCoefficientMacroporeCapacity>();
        mtc->tortuosity_ = 3.0;
        mtc->macropore_void_fraction_ = 0.35;
        mtc->molecular_diffusion_coefficient_ = kDm;
        ldf->mtc_model_ = ToResistances(std::move(mtc));
        calf20.adsorption_source_.mass_transfer_models_[comp] = std::move(ldf);
      }

      // Isosteric heats of Fig. 2d (CO2 about 37.5, N2 about 15 kJ/mol; read
      // from the figure, README).
      SetHOA(ads_id, CO2, 38.0e3);
      SetHOA(ads_id, N2, 15.0e3);

      AddConstantInternalHeatTransferCoefficientModel(6.963, layer_id);   // Table 2, fitted
      AddConstantExternalHeatTransferCoefficientModel(14.0, layer_id);    // Table 2, fitted

      // Eq. 6: D_L = 0.7 D_m + v r_p (the engine's Ruthven form), D_m of Table 2.
      auto disp = std::make_unique<DispersionModelRuthven>();
      disp->Dm = kDm;
      layer.equations_.dispersion_model_ = std::move(disp);
      auto kz = std::make_unique<HeatDispersionModelConstant>();
      kz->dispersion_coefficient_ = 0.4;                       // Table 2, fitted
      layer.equations_.heat_dispersion_model_ = std::move(kz);

      // Column saturated with helium (Sect. 3.2); feed at the case's conditions.
      AddInitialConditions(P, T, T, { { CO2, 0.0 }, { N2, 0.0 }, { He, 1.0 } });
      AddStep("Feed", c.duration);
      const double y = c.percent / 100.0;
      std::unordered_map<std::string, double> feed = { { CO2, y }, { N2, 0.0 }, { He, 0.0 } };
      feed[c.fig == "Fig4" ? He : N2] = 1.0 - y;               // Fig. 4: balance He; Fig. 5: balance N2
      auto stream = CreateNoValveInletStream(CreateVolumeFlowBC(c.flow_ccm * 1e-6 / 60.0), feed, T);
      stream->name_ = "Feed";
      auto inlet = std::make_unique<FlowBoundaryCondition>();
      inlet->direction_ = BoundaryDirection::INLET;
      inlet->location_ = BoundaryLocation::BOTTOM;
      inlet->streams_.push_back(std::move(stream));
      Step& step = data.cycle.GetStep("Feed");
      step.bottom_boundary_condition_ = std::move(inlet);
      auto outlet = std::make_unique<PressureBoundaryCondition>();
      outlet->direction_ = BoundaryDirection::OUTLET;
      outlet->location_ = BoundaryLocation::TOP;
      outlet->pressure_ = CreateConstantBoundaryValue(P);
      outlet->mole_fraction_[CO2] = 0.0;
      outlet->mole_fraction_[N2] = 0.0;
      outlet->mole_fraction_[He] = 1.0;
      step.top_boundary_condition_ = std::move(outlet);
      ConfigureCycle(1);
    }

    double BedLength() const {
      return kMass / (kParticleDensity * (1 - kBedVoidage) * 0.25 * M_PI * 0.0282 * 0.0282);
    }
    // The thermocouple sits 1.82 cm from the column outlet (Sect. 3.2).
    double ThermocoupleZ() const { return BedLength() - 0.0182; }

    static constexpr double kMass = 16.7e-3;                 // [kg] Table 2
    static constexpr double kParticleDensity = 570.0;        // [kg/m3] Table 2
    static constexpr double kBedVoidage = 0.4;               // Table 2
    static constexpr double kParticleDiameter = 1.5e-3;      // [m] assumption: pellets 1.0-3.0 mm (Sect. 2)
    static constexpr double kDm = 1.496e-5;                  // [m2/s] Table 2
    static constexpr double kCpCO2 = 841.0 * 0.0440095;      // [J/(mol K)] 841 J/(kg K)
    static constexpr double kCpN2 = 1040.0 * 0.0280134;      // [J/(mol K)] 1040 J/(kg K)
    static constexpr double kCpHe = 5200.0 * 0.004002602;    // [J/(mol K)] 5200 J/(kg K)
  };

} // namespace adsim
