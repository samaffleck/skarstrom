#pragma once

// ============================================================================
//  Haghpanah et al. (2013) - four-step CO2/N2 PVSA on zeolite 13X
//  See README.md in this folder for the citation, the provenance table and the
//  measured errors.
//
//  The deck: Haghpanah's column, bed, wall, gas and transport inputs (Table 3),
//  the Krishnamurthy et al. (2012) dual-site Langmuir fit (Table 2), the four
//  steps of Sect. 6.1 with exponential pressure profiles (time constant
//  lambda = 0.5 1/s, Eqs. 57-59) and the loading-dependent LDF coefficient of
//  Eq. 40. One function builds every case: the Sect. 6.1 reference design, the
//  five minimum-energy designs of Table 6 (also Glover et al. 2026, Table S2)
//  and the CSS benchmark design of Glover et al. (2026), Figure 5.
// ============================================================================

#define _USE_MATH_DEFINES
#include <cmath>
#include <memory>
#include <string>
#include <unordered_map>

#include "adsim/boundary_conditions/boundary_conditions.h"
#include "adsim/boundary_conditions/boundary_stream.h"
#include "adsim/equations/dispersion_model.h"
#include "adsim/equations/heat_transfer_model.h"
#include "adsim/equations/isotherm.h"
#include "adsim/equations/mass_transfer_model.h"
#include "adsim/equations/wall_temperature_boundary_condition.h"
#include "adsim/objects/compressor.h"

#include "../../fixtures/user_data_test.h"

namespace adsim {

  // One operating point. Pressures in Pa, times in s, v_feed the interstitial
  // feed velocity in m/s (Table 3: "interstitial feed velocity").
  struct HaghpanahDesign {
    double t_press, t_ads, t_bd, t_evac;
    double P_H, P_I, P_L;
    double v_feed;
  };

  // Rate model of the deck. CAPACITY is Haghpanah's Eq. 40 (the validation
  // deck); CONSTANT_K is Glover et al.'s constant 26.25 1/s (a variant, used
  // only to compare the CSS multiplier).
  enum class HaghpanahKinetics { CAPACITY, CONSTANT_K };

  class Haghpanah2013Test : public UserDataTest {
  public:
    static constexpr const char* kCO2 = "CarbonDioxide";
    static constexpr const char* kN2 = "Nitrogen";

    // Table 3 (p. 4255) and the text around it.
    static constexpr double kLength = 1.0;              // [m]
    static constexpr double kInnerRadius = 0.1445;      // [m]
    static constexpr double kOuterRadius = 0.1620;      // [m]
    static constexpr double kBedVoid = 0.37;            // [-]
    static constexpr double kParticleVoid = 0.35;       // [-]
    static constexpr double kParticleRadius = 1.0e-3;   // [m]
    static constexpr double kTortuosity = 3.0;          // [-]
    static constexpr double kParticleDensity = 1130.0;  // [kg/m3] adsorbent density rho_s
    static constexpr double kMolecularDiffusivity = 1.60e-5;  // [m2/s]
    static constexpr double kGasConductivity = 0.09;    // [W/(m K)] K_z
    static constexpr double kWallConductivity = 16.0;   // [W/(m K)]
    static constexpr double kGasCp = 30.7;              // [J/(mol K)] C_p,g = C_p,a
    static constexpr double kSolidCp = 1070.0;          // [J/(kg K)]
    static constexpr double kWallCp = 502.0;            // [J/(kg K)]
    static constexpr double kWallDensity = 7800.0;      // [kg/m3]
    static constexpr double kViscosity = 1.72e-5;       // [Pa s]
    static constexpr double kInsideHtc = 8.6;           // [W/(m2 K)]
    static constexpr double kOutsideHtc = 2.5;          // [W/(m2 K)]
    static constexpr double kFeedTemperature = 298.15;  // [K] = T_a
    static constexpr double kAtmosphere = 1.0e5;        // [Pa] P_f = P_atm = 1 bar (p. 4260)
    static constexpr double kEfficiency = 0.72;         // compression / evacuation efficiency
    static constexpr double kGamma = 1.4;               // adiabatic constant
    static constexpr double kLambda = 0.5;              // [1/s] pressure-profile time constant (Table 4)
    // Constant heats of adsorption (the paper prints none, Sect. 2.5 of the
    // README): CO2 the dual-site isosteric weighting of the two site energies
    // at the feed partial pressure (0.15 bar, 298.15 K), N2 its single site.
    static constexpr double kHoaCO2 = 35983.8;          // [J/mol]
    static constexpr double kHoaN2 = 15800.0;           // [J/mol]
    // Glover et al. (2026) Table 2: constant LDF coefficient (variant only).
    static constexpr double kGloverK = 26.25;           // [1/s]

    // One Table 6 design to CSS, compared with both papers (haghpanah_2013.cpp).
    void RunTable6(const std::string& id);

    // Builds the whole deck for one design. `cells` is the grid; the bed
    // starts full of N2 at P_L and 298.15 K (Sect. 6.1).
    void BuildHaghpanah(const std::string& title, const std::string& description,
                        const HaghpanahDesign& d, int cells,
                        HaghpanahKinetics kinetics = HaghpanahKinetics::CAPACITY) {
      data.info.title = title;
      data.info.description = description;
      data.info.author = "Skarstrom validation";

      const std::string CO2 = kCO2, N2 = kN2;
      AddComponents({ CO2, N2 });
      auto [layer_id, ads_id] = AddLayer("13X");

      // ---- gas: constant properties (Table 3; Appendix A text: mu constant) ----
      data.fluid.gas_heat_capacity_mode_ = GasHeatCapacityMode::CONSTANT;
      data.fluid.gas_heat_capacity_[CO2] = kGasCp;
      data.fluid.gas_heat_capacity_[N2] = kGasCp;
      data.fluid.adsorbed_heat_capacity_[CO2] = kGasCp;
      data.fluid.adsorbed_heat_capacity_[N2] = kGasCp;
      data.fluid.viscosity_mode_ = ViscosityMode::CONSTANT;
      data.fluid.viscosity_ = kViscosity;
      // Eq. A.8 carries -C_pa T dq/dt next to the heat of adsorption: the
      // internal-energy convention.
      data.fluid.heat_of_adsorption_convention_ = HeatOfAdsorptionConvention::INTERNAL_ENERGY;

      // ---- column and wall ----
      auto& layer = GetLayer(layer_id);
      layer.domain_.length = kLength;
      layer.domain_.n_cells = cells;
      layer.eb_ = kBedVoid;
      layer.particle_diameter_ = 2.0 * kParticleRadius;
      layer.momentum_model_ = MomentumModel::BlakeKozeny;  // Eq. 6: Darcy with the Blake-Kozeny constant 150
      std::get<CylinderBed>(data.bed.shape_).diameter = 2.0 * kInnerRadius;
      data.bed.wall_.thickness_ = kOuterRadius - kInnerRadius;
      data.bed.wall_.density_ = kWallDensity;
      data.bed.wall_.heat_capacity_ = kWallCp;
      data.bed.wall_.thermal_conductivity_ = kWallConductivity;
      data.bed.ambient_temperature_ = kFeedTemperature;
      // Eq. 12 / SI Table 2: T_w = T_a at both ends of the wall.
      auto wall_bc = std::make_unique<WallTemperatureBoundaryConditionFixed>();
      wall_bc->temperature_ = kFeedTemperature;
      data.bed.equations_.wall_temperature_boundary_condition_ = std::move(wall_bc);

      // ---- adsorbent ----
      auto& ads = GetAdsorbent(ads_id);
      ads.solid_density_ = kParticleDensity;
      ads.heat_capacity_ = kSolidCp;
      // The balances (Eqs. A.2-A.3) hold gas in the bed voidage only: no
      // macropore holdup (eps_t = eps_b). eps_p = 0.35 stays in the rate (Eq. 40).
      ads.ep_ = 0.0;
      ads.thermal_conductivity_ = 0.63;  // not used: gas and solid share one temperature

      // Table 2: competitive dual-site Langmuir in concentration form, b(T) =
      // b0 exp(-dU / RT); CO2 two sites, N2 one site (d0 = 0, q_sd = 0).
      auto dsl = std::make_unique<ExtendedDualSiteLangmuirIsotherm_Concentration>(CO2);
      dsl->params_[CO2].m1 = 3.09;
      dsl->params_[CO2].b1 = 8.65e-7;
      dsl->params_[CO2].q1 = 36641.21;
      dsl->params_[CO2].m2 = 2.54;
      dsl->params_[CO2].b2 = 2.63e-8;
      dsl->params_[CO2].q2 = 35690.66;
      dsl->params_[N2].m1 = 5.84;
      dsl->params_[N2].b1 = 2.50e-6;
      dsl->params_[N2].q1 = 15800.0;
      dsl->params_[N2].m2 = 0.0;
      dsl->params_[N2].b2 = 0.0;
      dsl->params_[N2].q2 = 0.0;
      AddIsotherm({ CO2, N2 }, ads_id, std::move(dsl));

      for (const std::string& c : { CO2, N2 }) {
        if (kinetics == HaghpanahKinetics::CAPACITY) {
          // Eq. 40-41: k_i = (c_i / q_i*) 15 eps_p D_p / r_p^2, D_p = D_m / tau.
          auto ldf = std::make_unique<LinearDrivingForce>();
          auto mtc = std::make_unique<MassTransferCoefficientMacroporeCapacity>();
          mtc->tortuosity_ = kTortuosity;
          mtc->macropore_void_fraction_ = kParticleVoid;
          mtc->molecular_diffusion_coefficient_ = kMolecularDiffusivity;
          ldf->mtc_model_ = ToResistances(std::move(mtc));  // generic model, as the example decks store it
          ads.adsorption_source_.mass_transfer_models_[c] = std::move(ldf);
        } else {
          AddLinearDrivingForceModel(ads_id, c, kGloverK);
        }
      }
      SetHOA(ads_id, CO2, kHoaCO2);
      SetHOA(ads_id, N2, kHoaN2);

      // ---- transport: constant per design (Eq. 42, Table 3) ----
      auto dispersion = std::make_unique<DispersionModelConstant>();
      dispersion->dispersion_coefficient_ = 0.7 * kMolecularDiffusivity + d.v_feed * kParticleRadius;
      layer.equations_.dispersion_model_ = std::move(dispersion);
      auto conduction = std::make_unique<HeatDispersionModelConstant>();
      conduction->dispersion_coefficient_ = kGasConductivity;
      layer.equations_.heat_dispersion_model_ = std::move(conduction);
      AddConstantInternalHeatTransferCoefficientModel(kInsideHtc, layer_id);
      AddConstantExternalHeatTransferCoefficientModel(kOutsideHtc, layer_id);

      // ---- initial state: N2 at P_L (Sect. 6.1) ----
      AddInitialConditions(d.P_L, kFeedTemperature, kFeedTemperature, { { CO2, 0.0 }, { N2, 1.0 } });

      // ---- cycle ----
      const std::unordered_map<std::string, double> feed_y = { { CO2, 0.15 }, { N2, 0.85 } };
      AddStep("Pressurisation", d.t_press);
      AddStep("Adsorption", d.t_ads);
      AddStep("Blowdown", d.t_bd);
      AddStep("Evacuation", d.t_evac);

      // Pressurisation with feed gas at the feed end, P_L -> P_H (Eq. 57).
      {
        auto bc = CreatePressureInletBoundary(Ramp(d.P_H), feed_y, kFeedTemperature);
        bc->location_ = BoundaryLocation::BOTTOM;
        bc->direction_ = BoundaryDirection::INLET;
        bc->pump_ = MakeCompressor();
        Step& s = data.cycle.GetStep("Pressurisation");
        s.bottom_boundary_condition_ = std::move(bc);
        AddClosedBC(BoundaryLocation::TOP, "Pressurisation");
      }
      // Adsorption: feed at v_feed (superficial eps v_feed), raffinate at P_H.
      {
        auto stream = CreateNoValveInletStream(CreateVelocityBC(kBedVoid * d.v_feed), feed_y, kFeedTemperature);
        stream->name_ = "Feed inlet";
        stream->pump_ = MakeCompressor();
        auto inlet = std::make_unique<FlowBoundaryCondition>();
        inlet->direction_ = BoundaryDirection::INLET;
        inlet->location_ = BoundaryLocation::BOTTOM;
        inlet->streams_.push_back(std::move(stream));
        data.cycle.GetStep("Adsorption").bottom_boundary_condition_ = std::move(inlet);
        AddPressureOutletBC(CreateConstantBoundaryValue(d.P_H), BoundaryLocation::TOP, "Adsorption",
                            { { CO2, 0.0 }, { N2, 1.0 } });
      }
      // Forward blowdown at the product end, P_H -> P_I (Eq. 58).
      {
        AddClosedBC(BoundaryLocation::BOTTOM, "Blowdown");
        AddPressureOutletBC(Ramp(d.P_I), BoundaryLocation::TOP, "Blowdown", { { CO2, 0.0 }, { N2, 1.0 } });
        static_cast<PressureBoundaryCondition&>(*data.cycle.GetStep("Blowdown").top_boundary_condition_).pump_ = MakeVacuumPump();
      }
      // Reverse evacuation at the feed end, P_I -> P_L (Eq. 59): the CO2 product.
      {
        AddClosedBC(BoundaryLocation::TOP, "Evacuation");
        AddPressureOutletBC(Ramp(d.P_L), BoundaryLocation::BOTTOM, "Evacuation", { { CO2, 0.0 }, { N2, 1.0 } });
        static_cast<PressureBoundaryCondition&>(*data.cycle.GetStep("Evacuation").bottom_boundary_condition_).pump_ = MakeVacuumPump();
      }

      // KPIs (Eqs. 60-61): CO2 fed in pressurisation and adsorption, recovered by evacuation.
      ConfigureKPI(0, { 0, 1 }, { BoundaryLocation::BOTTOM, BoundaryLocation::BOTTOM },
                   { 3 }, { BoundaryLocation::BOTTOM });

      // ---- solver ----
      ConfigureCycle(500);  // Sect. 6.2: up to 500 cycles
      data.settings.css.method_ = CSSMethod::MASS_BALANCE_COMPONENT;
      data.settings.css.min_cycles_ = 5;
      data.settings.max_step = 2.0;
    }

    // exp(-lambda t) approach from the bed pressure at the start of the step to
    // `end` (Eqs. 57-59): the two-term profile with both rates lambda.
    static std::unique_ptr<IBoundaryValue> Ramp(double end) {
      auto bv = std::make_unique<DoubleExponentialBoundaryValue>();
      bv->end_value_ = end;
      bv->rate_1_ = kLambda;
      bv->rate_2_ = kLambda;
      bv->weight_ = 1.0;
      return bv;
    }

    // Eqs. 66-69: isentropic work gamma/(gamma-1) P Q [(P2/P1)^((gamma-1)/gamma) - 1]
    // over the efficiency 0.72; the feed is compressed from 1 bar, the vacuum
    // pumps discharge to 1 bar, and the gas enters the bed at the feed
    // temperature (isothermal outlet model; the power still uses gamma).
    static std::unique_ptr<IPump> MakeCompressor() {
      auto p = CreatePump(PumpType::COMPRESSOR);
      ConfigurePump(*p);
      static_cast<adsim::Compressor&>(*p).suction_pressure_ = kAtmosphere;
      return p;
    }
    static std::unique_ptr<IPump> MakeVacuumPump() {
      auto p = CreatePump(PumpType::VACUUM_PUMP);
      ConfigurePump(*p);
      static_cast<adsim::VacuumPump&>(*p).discharge_pressure_ = kAtmosphere;
      return p;
    }
    static void ConfigurePump(IPump& p) {
      p.use_pump_ = true;
      p.poly_index_ = kGamma;
      p.energy_eq_ = CreateCompressorEnergyModel(CompressorEnergyType::ISOTHERMAL);
      auto eff = std::make_unique<ConstantPumpEfficiencyModel>();
      eff->efficiency_ = kEfficiency;
      p.efficiency_model_ = std::move(eff);
    }
  };

} // namespace adsim
