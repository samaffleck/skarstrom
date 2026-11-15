#pragma once

// ============================================================================
//  Ward & Pini (2022), Ind. Eng. Chem. Res. 61, 13650-13668: four-step CO2/N2
//  PVSA with feed pressurisation on zeolite 13X and ZIF-36-FRL.
//  See README.md in this folder for the citation, the provenance table and the
//  measured errors.
//
//  The deck: the column, bed, wall, gas and transport inputs of Table 3, the
//  extended dual-site Langmuir parameters of Table 1, the boundary conditions
//  of SI Sect. S.2 (exponential pressure profiles, lambda = 0.5 1/s), the LDF
//  coefficient of Eq. 6 and the energy model of SI Sect. S.4 (isentropic
//  compression and evacuation with the efficiency eta(p) of Eq. 11). One
//  function builds every case: the ten cost-optimal designs of SI Tables S2
//  (zeolite 13X) and S3 (ZIF-36-FRL), whose column length is a decision
//  variable at the fixed aspect ratio L / r_in = 6.43 (Sect. 6.1).
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

  enum class WardPiniMaterial { ZEOLITE_13X, ZIF_36_FRL };

  // One operating point of SI Table S2 / S3. Pressures in Pa, times in s,
  // v_feed the interstitial feed velocity in m/s, L the column length in m.
  struct WardPiniDesign {
    WardPiniMaterial material;
    double P_H, P_I, P_L;
    double t_ads, t_bd, t_evac;
    double v_feed;
    double length;
  };

  class WardPini2022Test : public UserDataTest {
  public:
    static constexpr const char* kCO2 = "CarbonDioxide";
    static constexpr const char* kN2 = "Nitrogen";

    // Table 3 (p. 13653).
    static constexpr double kAspectRatio = 6.43;           // [-] L / r_in (Sect. 6.1)
    static constexpr double kRadiusRatio = 0.162 / 0.145;  // [-] r_out / r_in of the Table 3 column
    static constexpr double kBedDensity = 712.0;           // [kg/m3] rho_b
    static constexpr double kBedVoid = 0.37;               // [-]
    static constexpr double kParticleVoid = 0.35;          // [-]
    static constexpr double kParticleRadius = 1.0e-3;      // [m]
    static constexpr double kTortuosity = 3.0;             // [-]
    static constexpr double kMolecularDiffusivity = 1.5e-5;  // [m2/s]
    static constexpr double kGasConductivity = 0.09;       // [W/(m K)] K_z
    static constexpr double kWallConductivity = 16.0;      // [W/(m K)]
    static constexpr double kGasCp = 30.7;                 // [J/(mol K)] C_p,g = C_p,a
    static constexpr double kSolidCp = 1070.0;             // [J/(kg K)]
    static constexpr double kWallCp = 502.0;               // [J/(kg K)]
    static constexpr double kWallDensity = 7800.0;         // [kg/m3]
    static constexpr double kViscosity = 1.72e-5;          // [Pa s]
    static constexpr double kInsideHtc = 8.6;              // [W/(m2 K)]
    static constexpr double kOutsideHtc = 2.5;             // [W/(m2 K)]
    static constexpr double kGamma = 1.4;                  // ratio of heat capacities
    static constexpr double kLambda = 0.5;                 // [1/s] pressure-profile time constant
    // Feed and surroundings (Sect. 2.1, Table 3): 1 bar, 298.15 K; the
    // compressors draw from p_F and the vacuum pumps discharge to p_atm (SI
    // Eqs. 45-48), both 1 bar.
    static constexpr double kFeedPressure = 1.0e5;         // [Pa]
    static constexpr double kFeedTemperature = 298.15;     // [K] = T_a
    static constexpr double kFeedCO2 = 0.15;               // [-]
    static constexpr double kPressurisationTime = 20.0;    // [s] fixed (Sect. 3.3)
    static constexpr double kR = 8.314462618;              // [J/(mol K)]

    // Pellet density: rho_b / (1 - eps), so that the bed holds 712 kg/m3.
    static constexpr double ParticleDensity() { return kBedDensity / (1.0 - kBedVoid); }
    static double InnerRadius(const WardPiniDesign& d) { return d.length / kAspectRatio; }
    static double ColumnVolume(const WardPiniDesign& d) {
      const double r = InnerRadius(d);
      return M_PI * r * r * d.length;
    }

    // One design of SI Table S2 / S3 to CSS, compared with the paper (ward_pini_2022.cpp).
    void RunDesign(const std::string& id);

    // Builds the whole deck for one design. `cells` is the grid; the bed
    // starts full of N2 at P_H and 298.15 K (p. 13653).
    void BuildWardPini(const std::string& title, const std::string& description, const WardPiniDesign& d,
                       int cells) {
      data.info.title = title;
      data.info.description = description;
      data.info.author = "Skarstrom validation";

      const std::string CO2 = kCO2, N2 = kN2;
      const bool zif = d.material == WardPiniMaterial::ZIF_36_FRL;
      AddComponents({ CO2, N2 });
      auto [layer_id, ads_id] = AddLayer(zif ? "ZIF-36-FRL" : "13X");

      // ---- gas: constant properties (Table 3) ----
      data.fluid.gas_heat_capacity_mode_ = GasHeatCapacityMode::CONSTANT;
      data.fluid.gas_heat_capacity_[CO2] = kGasCp;
      data.fluid.gas_heat_capacity_[N2] = kGasCp;
      data.fluid.adsorbed_heat_capacity_[CO2] = kGasCp;
      data.fluid.adsorbed_heat_capacity_[N2] = kGasCp;
      data.fluid.viscosity_mode_ = ViscosityMode::CONSTANT;
      data.fluid.viscosity_ = kViscosity;
      // The column energy balance (Table 2, SI Eq. 8) has the source
      // sum_i (sigma_i - Omega_3 T) dx_i/dtau: heat -dU_b,i + R T per mole
      // adsorbed, less the adsorbed-phase sensible term C_p,a T. That is the
      // internal-energy convention with -dU_b,i + R T_F as the stored heat.
      data.fluid.heat_of_adsorption_convention_ = HeatOfAdsorptionConvention::INTERNAL_ENERGY;

      // ---- column and wall ----
      const double r_in = InnerRadius(d);
      auto& layer = GetLayer(layer_id);
      layer.domain_.length = d.length;
      layer.domain_.n_cells = cells;
      layer.eb_ = kBedVoid;
      layer.particle_diameter_ = 2.0 * kParticleRadius;
      layer.momentum_model_ = MomentumModel::BlakeKozeny;  // Table 2: Darcy, 150 (1 - eps)^2 / eps^2
      std::get<CylinderBed>(data.bed.shape_).diameter = 2.0 * r_in;
      data.bed.wall_.thickness_ = (kRadiusRatio - 1.0) * r_in;
      data.bed.wall_.density_ = kWallDensity;
      data.bed.wall_.heat_capacity_ = kWallCp;
      data.bed.wall_.thermal_conductivity_ = kWallConductivity;
      data.bed.ambient_temperature_ = kFeedTemperature;
      // SI Eq. 18: T_w = T_a at both ends of the wall, in every step.
      auto wall_bc = std::make_unique<WallTemperatureBoundaryConditionFixed>();
      wall_bc->temperature_ = kFeedTemperature;
      data.bed.equations_.wall_temperature_boundary_condition_ = std::move(wall_bc);

      // ---- adsorbent ----
      auto& ads = GetAdsorbent(ads_id);
      ads.solid_density_ = ParticleDensity();
      ads.heat_capacity_ = kSolidCp;
      // The balances (Table 2, SI Eq. 1: psi = R T0 rho_b q_s0 / (p0 eps))
      // hold gas in the bed voidage only: no macropore holdup. eps_p = 0.35
      // stays in the rate (Eq. 6).
      ads.ep_ = 0.0;
      ads.thermal_conductivity_ = 0.63;  // not used: gas and solid share one temperature

      // Table 1: extended dual-site Langmuir in concentration form (Eqs. 1-4),
      // b(T) = b0 exp(-dU / RT); CO2 two sites, N2 one site.
      auto dsl = std::make_unique<ExtendedDualSiteLangmuirIsotherm_Concentration>(CO2);
      double dU_b_CO2 = 0.0, dU_b_N2 = 0.0;
      if (zif) {
        dsl->params_[CO2].m1 = 2.01;
        dsl->params_[CO2].b1 = 5.00e-8;
        dsl->params_[CO2].q1 = 38100.0;
        dsl->params_[CO2].m2 = 1.01;
        dsl->params_[CO2].b2 = 9.10e-7;
        dsl->params_[CO2].q2 = 21900.0;
        dsl->params_[N2].m1 = 3.02;
        dsl->params_[N2].b1 = 6.29e-7;
        dsl->params_[N2].q1 = 16230.0;
        dU_b_CO2 = 38100.0;
        dU_b_N2 = 16230.0;
      } else {
        dsl->params_[CO2].m1 = 3.09;
        dsl->params_[CO2].b1 = 8.65e-7;
        dsl->params_[CO2].q1 = 36600.0;
        dsl->params_[CO2].m2 = 2.54;
        dsl->params_[CO2].b2 = 2.63e-8;
        dsl->params_[CO2].q2 = 35700.0;
        dsl->params_[N2].m1 = 5.84;
        dsl->params_[N2].b1 = 2.50e-6;
        dsl->params_[N2].q1 = 15800.0;
        dU_b_CO2 = 36600.0;
        dU_b_N2 = 15800.0;
      }
      dsl->params_[N2].m2 = 0.0;
      dsl->params_[N2].b2 = 0.0;
      dsl->params_[N2].q2 = 0.0;
      AddIsotherm({ CO2, N2 }, ads_id, std::move(dsl));

      // Eq. 6: k_i = (c_i / q_i*) 15 eps_p D_m / (tau_p r_p^2), macropore
      // diffusion with the Hassan et al. extension (Haynes and Sarma).
      for (const std::string& c : { CO2, N2 }) {
        auto ldf = std::make_unique<LinearDrivingForce>();
        auto mtc = std::make_unique<MassTransferCoefficientMacroporeCapacity>();
        mtc->tortuosity_ = kTortuosity;
        mtc->macropore_void_fraction_ = kParticleVoid;
        mtc->molecular_diffusion_coefficient_ = kMolecularDiffusivity;
        ldf->mtc_model_ = ToResistances(std::move(mtc));
        ads.adsorption_source_.mass_transfer_models_[c] = std::move(ldf);
      }
      // SI Eq. 8: -dU_b,i + R T, with T at the feed temperature.
      SetHOA(ads_id, CO2, dU_b_CO2 + kR * kFeedTemperature);
      SetHOA(ads_id, N2, dU_b_N2 + kR * kFeedTemperature);

      // ---- transport (Eq. 5, Table 3) ----
      auto dispersion = std::make_unique<DispersionModelConstant>();
      dispersion->dispersion_coefficient_ = 0.7 * kMolecularDiffusivity + d.v_feed * kParticleRadius;
      layer.equations_.dispersion_model_ = std::move(dispersion);
      auto conduction = std::make_unique<HeatDispersionModelConstant>();
      conduction->dispersion_coefficient_ = kGasConductivity;
      layer.equations_.heat_dispersion_model_ = std::move(conduction);
      AddConstantInternalHeatTransferCoefficientModel(kInsideHtc, layer_id);
      AddConstantExternalHeatTransferCoefficientModel(kOutsideHtc, layer_id);

      // ---- initial state: N2 at P_H (p. 13653) ----
      AddInitialConditions(d.P_H, kFeedTemperature, kFeedTemperature, { { CO2, 0.0 }, { N2, 1.0 } });

      // ---- cycle (Sect. 2.3, Fig. 2, SI Sect. S.2) ----
      const std::unordered_map<std::string, double> feed_y = { { CO2, kFeedCO2 }, { N2, 1.0 - kFeedCO2 } };
      AddStep("Pressurisation", kPressurisationTime);
      AddStep("Adsorption", d.t_ads);
      AddStep("Blowdown", d.t_bd);
      AddStep("Evacuation", d.t_evac);

      // Feed pressurisation at the feed end, P_L -> P_H (SI Eq. 36).
      {
        auto bc = CreatePressureInletBoundary(Ramp(d.P_H), feed_y, kFeedTemperature);
        bc->location_ = BoundaryLocation::BOTTOM;
        bc->direction_ = BoundaryDirection::INLET;
        bc->pump_ = MakeCompressor();
        data.cycle.GetStep("Pressurisation").bottom_boundary_condition_ = std::move(bc);
        AddClosedBC(BoundaryLocation::TOP, "Pressurisation");
      }
      // Adsorption: feed at v_feed (superficial eps v_feed, SI Eq. 45), light
      // product at P_H (SI Eqs. 12-13).
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
      // Forward blowdown at the product end, P_H -> P_I (SI Eq. 24).
      {
        AddClosedBC(BoundaryLocation::BOTTOM, "Blowdown");
        AddPressureOutletBC(Ramp(d.P_I), BoundaryLocation::TOP, "Blowdown", { { CO2, 0.0 }, { N2, 1.0 } });
        static_cast<PressureBoundaryCondition&>(*data.cycle.GetStep("Blowdown").top_boundary_condition_).pump_ =
          MakeVacuumPump();
      }
      // Reverse evacuation at the feed end, P_I -> P_L (SI Eq. 30): the CO2 product.
      {
        AddClosedBC(BoundaryLocation::TOP, "Evacuation");
        AddPressureOutletBC(Ramp(d.P_L), BoundaryLocation::BOTTOM, "Evacuation", { { CO2, 0.0 }, { N2, 1.0 } });
        static_cast<PressureBoundaryCondition&>(*data.cycle.GetStep("Evacuation").bottom_boundary_condition_).pump_ =
          MakeVacuumPump();
      }

      // KPIs (Eqs. 7-8): CO2 fed in pressurisation and adsorption, recovered by evacuation.
      ConfigureKPI(0, { 0, 1 }, { BoundaryLocation::BOTTOM, BoundaryLocation::BOTTOM }, { 3 },
                   { BoundaryLocation::BOTTOM });

      // ---- solver ----
      ConfigureCycle(500);
      data.settings.css.method_ = CSSMethod::MASS_BALANCE_COMPONENT;
      data.settings.css.min_cycles_ = 5;
      data.settings.max_step = 2.0;
    }

    // exp(-lambda t) approach from the bed pressure at the start of the step
    // to `end` (SI Eqs. 24, 30, 36): the two-term profile with both rates lambda.
    static std::unique_ptr<IBoundaryValue> Ramp(double end) {
      auto bv = std::make_unique<DoubleExponentialBoundaryValue>();
      bv->end_value_ = end;
      bv->rate_1_ = kLambda;
      bv->rate_2_ = kLambda;
      bv->weight_ = 1.0;
      return bv;
    }

    // SI Eqs. 45-48: isentropic work gamma/(gamma-1) P Q [(P2/P1)^((gamma-1)/gamma) - 1]
    // over eta(p) = 0.8 * 19.55 p / (1 + 19.55 p), p in bar (Eq. 11, SI Eq. 49;
    // Table 3 gives eta = eta(p) for every pump). The engine evaluates eta at
    // the suction pressure: the bed for a vacuum pump, p_F = 1 bar for a
    // compressor (eta = 0.761). Work is zero when the pressure ratio is below
    // one (SI Eqs. 46, 48). The gas enters the bed at the feed temperature.
    static std::unique_ptr<IPump> MakeCompressor() {
      auto p = CreatePump(PumpType::COMPRESSOR);
      ConfigurePump(*p);
      static_cast<adsim::Compressor&>(*p).suction_pressure_ = kFeedPressure;
      return p;
    }
    static std::unique_ptr<IPump> MakeVacuumPump() {
      auto p = CreatePump(PumpType::VACUUM_PUMP);
      ConfigurePump(*p);
      static_cast<adsim::VacuumPump&>(*p).discharge_pressure_ = kFeedPressure;
      return p;
    }
    static void ConfigurePump(IPump& p) {
      p.use_pump_ = true;
      p.poly_index_ = kGamma;
      p.energy_eq_ = CreateCompressorEnergyModel(CompressorEnergyType::ISOTHERMAL);
      auto eff = std::make_unique<VacuumPumpEfficiencyModel>();
      eff->eta_max_ = 0.8;
      eff->alpha_ = 19.55;
      p.efficiency_model_ = std::move(eff);
    }
  };

} // namespace adsim
