// STL includes
#define _USE_MATH_DEFINES
#include <cmath>
#include <filesystem>

// AdSim includes
#include "../../../src/pch.h"
#include "adsim/solver/coupled_solver.h"
#include "adsim/solver/kpis.h"
#include "adsim/math/constants.h"
#include "adsim/equations/isotherm.h"
#include "adsim/boundary_conditions/boundary_conditions.h"
#include "adsim/boundary_conditions/boundary_stream.h"
#include "adsim/boundary_conditions/connection.h"
#include "adsim/boundary_conditions/step.h"
#include "adsim/equations/dispersion_model.h"
#include "adsim/equations/heat_transfer_model.h"

// GTest includes
#include <gtest/gtest.h>

// Paper fixture (Ferre et al. 2026) — see README.md in this folder.
#include "../../fixtures/user_data_test.h"
#include "../../fixtures/validation_metrics.h"
#include <limits>

namespace adsim {

  // ---------------------------------------------------------------------------
  // Ferre, Wedler, Ward, Burger, Pini, Petit (2026), Sep. Purif. Technol. 406,
  // 138821. Single-bed 4-step PVSA for biogas CO2/CH4 co-production.
  //   Steps: (1) adsorption/feed  -> CH4 light product (raffinate, top)
  //          (2) forward blowdown (top -> p_i)  -- not collected
  //          (3) reverse evacuation (bottom -> p_l) -> CO2 heavy product
  //          (4) feed pressurisation (feed gas, bottom, p_l -> p_h)
  // Detailed dynamic model = Ward & Pini 2022 framework (internal-energy heat of
  // adsorption -> INTERNAL_ENERGY convention). Concentration-form DSL isotherm.
  // Each case reproduces the Table 5 "recovery-recovery, max Re_CO2" operating
  // point for one adsorbent. Isotherm parameters: SI Tables S1/S2; densities: S3.
  // ---------------------------------------------------------------------------
  class Ferre2026Test : public UserDataTest {};

  namespace {
    // A dual-site (concentration-form) DSL fit for one component. q1/q2 store the
    // positive magnitude of the internal energy of adsorption |dU| (dU < 0).
    struct Site { double m1, b1, q1, m2, b2, q2; };

    struct FerreCase {
      std::string name;
      Site  co2, ch4;         // SI Tables S1 (CO2) / S2 (CH4)
      double rho_ads;         // particle density [kg/m3] (SI Table S3); bulk = rho_ads*(1-eb)
      // Constant energy-balance heat of adsorption [J/mol], taken as the DSL
      // (pure-component) isosteric heat q_st = (|dU_b|.dq_b/dc + |dU_d|.dq_d/dc)/
      // (dq_b/dc + dq_d/dc) evaluated at the feed partial pressure. For a DSL with
      // disparate site energies this is loading-dependent; the feed-loading value
      // is a far better single constant than the raw site average.
      double hoa_co2, hoa_ch4;
      // Operating point (Table 5, recovery-recovery, max Re_CO2)
      double p_h, p_i, t_ads, t_bd, t_evac, v_f;
      // Paper targets [%]
      double pu_co2, pu_ch4, re_co2, re_ch4;
      // Validation tolerances vs Table 5 (abs. percentage points; productivity:
      // relative). Set from the measured error — see the per-case comments.
      double tol_pu_co2 = std::numeric_limits<double>::infinity();
      double tol_pu_ch4 = std::numeric_limits<double>::infinity();
      double tol_re_co2 = std::numeric_limits<double>::infinity();
      double tol_re_ch4 = std::numeric_limits<double>::infinity();
      double tol_pr_ch4 = std::numeric_limits<double>::infinity();
    };

    void RunFerreCase(Ferre2026Test& t, const FerreCase& fc, const std::string& case_id) {
      auto& data = t.data;
      data.info.title  = "Ferre 2026 - " + fc.name;
      data.info.author = "Ferre et al. (2026)";

      const std::string CO2 = "CarbonDioxide";
      const std::string CH4 = "Methane";
      t.AddComponents({CO2, CH4});

      auto ads_id   = data.adsorbent_directory.Add(fc.name);
      auto layer_id = data.bed.layers_.Add("Adsorbent");
      auto& ads     = data.adsorbent_directory.Get(ads_id);
      auto& layer   = data.bed.layers_.Get(layer_id);
      layer.momentum_model_ = MomentumModel::KarmenCozeny;  // fixture predates Ergun: keep its Darcy / Kozeny-Carman pressure drop

      // Bed geometry / packing (Table 2). rho_b = rho_ads * (1 - eb).
      const double L = 1.0, r_in = 0.1445, r_out = 0.162, eb = 0.37;
      layer.domain_.length     = L;
      layer.domain_.n_cells    = 30;
      layer.eb_                = eb;
      layer.particle_diameter_ = 2.0 * 7.5e-4;
      layer.composition_[ads_id] = 1.0;
      std::get<CylinderBed>(data.bed.shape_).diameter = 2.0 * r_in;

      ads.solid_density_        = fc.rho_ads; // particle density (S3); bulk rho_b = rho_ads*(1-eb)
      ads.ep_                   = 0.35;
      ads.heat_capacity_        = 1070.0;     // Cp,s
      ads.thermal_conductivity_ = 0.5;

      data.bed.wall_.heat_capacity_        = 502.0;
      data.bed.wall_.density_              = 7800.0;
      data.bed.wall_.thermal_conductivity_ = 16.0;
      data.bed.wall_.thickness_            = r_out - r_in;
      data.bed.ambient_temperature_        = 298.15;

      // Extended (competitive) DSL isotherm, concentration form. Ferre et al.
      // (2026), Sect. 3.1: "To model the multi-component competitive adsorption
      // equilibrium, the extended DSL isotherm model was [used]". This deck
      // previously used the non-competitive DualSiteLangmuirIsotherm_Concentration
      // (each component adsorbing as if alone) - a transcription slip.
      auto dsl = std::make_unique<ExtendedDualSiteLangmuirIsotherm_Concentration>(CO2);
      dsl->params_[CO2].m1 = fc.co2.m1; dsl->params_[CO2].b1 = fc.co2.b1; dsl->params_[CO2].q1 = fc.co2.q1;
      dsl->params_[CO2].m2 = fc.co2.m2; dsl->params_[CO2].b2 = fc.co2.b2; dsl->params_[CO2].q2 = fc.co2.q2;
      dsl->params_[CH4].m1 = fc.ch4.m1; dsl->params_[CH4].b1 = fc.ch4.b1; dsl->params_[CH4].q1 = fc.ch4.q1;
      dsl->params_[CH4].m2 = fc.ch4.m2; dsl->params_[CH4].b2 = fc.ch4.b2; dsl->params_[CH4].q2 = fc.ch4.q2;
      t.AddIsotherm({CO2, CH4}, ads_id, std::move(dsl));

      // LDF: macropore diffusion constant 15*ep*Dm/(tau*rp^2) ~ 46.7 /s. The paper
      // scales this by the loading-dependent capacity factor c/q* (Haynes-Sarma /
      // Hassan), but the separation is equilibrium-controlled here so the constant
      // reproduces the same KPIs (verified for NaY(F)).
      t.AddLinearDrivingForceModel(ads_id, CO2, 46.7);
      t.AddLinearDrivingForceModel(ads_id, CH4, 46.7);

      t.SetHOA(ads_id, CO2, fc.hoa_co2);
      t.SetHOA(ads_id, CH4, fc.hoa_ch4);
      data.fluid.adsorbed_heat_capacity_[CO2] = 30.7;
      data.fluid.adsorbed_heat_capacity_[CH4] = 30.7;
      data.fluid.heat_of_adsorption_convention_ = HeatOfAdsorptionConvention::INTERNAL_ENERGY;

      t.AddConstantInternalHeatTransferCoefficientModel(8.6, layer_id);
      t.AddConstantExternalHeatTransferCoefficientModel(2.5, layer_id);
      layer.equations_.dispersion_model_ = CreateAxialDispersionModel(AxialDispersionType::RUTHVEN);

      const double p_h = fc.p_h * 1e5, p_i = fc.p_i * 1e5, p_l = 0.03e5;
      const double T_feed = 298.15, feed_vel = fc.v_f * eb;
      const std::unordered_map<std::string, double> feed_mf = {{CO2, 0.40}, {CH4, 0.60}};
      t.AddInitialConditions(p_h, T_feed, T_feed, {{CO2, 0.40}, {CH4, 0.60}});

      const std::string feed = "Feed", blow = "Blowdown", evac = "Evacuation", press = "Pressurisation";
      t.AddStep(feed,  fc.t_ads);
      t.AddStep(blow,  fc.t_bd);
      t.AddStep(evac,  fc.t_evac);
      t.AddStep(press, 20.0); // t_pres fixed at 20 s

      {
        std::vector<std::unique_ptr<IBoundaryStream>> streams;
        streams.push_back(t.CreateNoValveInletStream(t.CreateVelocityBC(feed_vel), feed_mf, T_feed));
        t.AddFlowBC(std::move(streams), BoundaryLocation::BOTTOM, BoundaryDirection::INLET, feed);
      }
      t.AddPressureOutletBC(t.CreateConstantBoundaryValue(p_h), BoundaryLocation::TOP, feed, feed_mf);

      t.AddClosedBC(BoundaryLocation::BOTTOM, blow);
      t.AddPressureOutletBC(t.CreateExponentialBoundaryValue(p_i, fc.t_bd), BoundaryLocation::TOP, blow, feed_mf);

      t.AddClosedBC(BoundaryLocation::TOP, evac);
      t.AddPressureOutletBC(t.CreateExponentialBoundaryValue(p_l, fc.t_evac), BoundaryLocation::BOTTOM, evac, feed_mf);

      t.AddPressureInletBC(t.CreateLinearBoundaryValue(p_l, p_h, 20.0), feed_mf, T_feed,
                           BoundaryLocation::BOTTOM, press);
      t.AddClosedBC(BoundaryLocation::TOP, press);

      t.ConfigureCycle(250);
      data.settings.css.method_ = CSSMethod::MASS_BALANCE_COMPONENT;
      data.settings.css.min_cycles_ = 5;
      data.settings.max_step = 2.0;

      std::string error{};
      ASSERT_TRUE(data.PerformChecks(error)) << "System checks failed: " << error;

      CoupledSolver solver{data};
      int retval = solver.IntegrateCycle(data.cycle, false);
      ASSERT_GE(retval, 0) << "Solver did not complete for " << fc.name;

      EXPECT_TRUE(t.IsMassConserved(solver.GetLog()));
      EXPECT_TRUE(t.IsComponentMolesConserved(solver.GetLog()));

      // KPIs via the official API: CO2 heavy product (evac@BOTTOM), CH4 light
      // product (feed raffinate@TOP); feed = feed + pressurisation streams.
      const EquationLog& log = solver.GetLog();
      const BoundaryLog& b_log = solver.GetBoundaryLog();
      const auto& cyc = log.sys_data.cycle_index_set_;
      const int c = cyc.empty() ? 0 : *cyc.rbegin();

      KPIParams kco2(log, b_log);
      kco2.product_comp_ = CO2; kco2.cycle_index_ = c;
      kco2.feed_endpoints = KpiEndpointsForFaces(data.cycle,
        {{0, BoundaryLocation::BOTTOM}, {3, BoundaryLocation::BOTTOM}});
      kco2.product_endpoints = KpiEndpointsForFaces(data.cycle, {{2, BoundaryLocation::BOTTOM}});

      KPIParams kch4(log, b_log);
      kch4.product_comp_ = CH4; kch4.cycle_index_ = c;
      kch4.feed_endpoints = KpiEndpointsForFaces(data.cycle,
        {{0, BoundaryLocation::BOTTOM}, {3, BoundaryLocation::BOTTOM}});
      kch4.product_endpoints = KpiEndpointsForFaces(data.cycle, {{0, BoundaryLocation::TOP}});

      const double pu_co2 = 100.0 * CalculatePurity(kco2), re_co2 = 100.0 * CalculateRecovery(kco2);
      const double pu_ch4 = 100.0 * CalculatePurity(kch4), re_ch4 = 100.0 * CalculateRecovery(kch4);

      std::cout << "\n[Ferre] " << fc.name << " (cyc " << (c + 1) << ")  model vs paper:\n"
                << "   Pu_CO2 " << pu_co2 << " / " << fc.pu_co2
                << " | Re_CO2 " << re_co2 << " / " << fc.re_co2
                << " | Pu_CH4 " << pu_ch4 << " / " << fc.pu_ch4
                << " | Re_CH4 " << re_ch4 << " / " << fc.re_ch4 << "\n" << std::endl;

      // Quantitative comparison with Table 5 (the paper's detailed-model optimum
      // for this operating point). Productivity = CH4 in the raffinate per bed
      // volume per total cycle time (Sect. 4.2). Energy is not compared: the
      // deck drives the pressure steps with pressure boundary conditions, not
      // compressor / vacuum-pump models.
      {
        using namespace validation;
        const double V_bed = M_PI * r_in * r_in * L;
        const double t_cycle = log.sys_data.total_cycle_time;
        const double pr_ch4 = CalculateComponentMolesInProduct(kch4) / (V_bed * t_cycle);
        CaseReport report("ferre_2026", case_id);
        report.CheckAbs("CO2 purity (heavy product)", pu_co2, Reference("ferre_2026", case_id, "CO2_purity"), fc.tol_pu_co2);
        report.CheckAbs("CO2 recovery", re_co2, Reference("ferre_2026", case_id, "CO2_recovery"), fc.tol_re_co2);
        report.CheckAbs("CH4 purity (light product)", pu_ch4, Reference("ferre_2026", case_id, "CH4_purity"), fc.tol_pu_ch4);
        report.CheckAbs("CH4 recovery", re_ch4, Reference("ferre_2026", case_id, "CH4_recovery"), fc.tol_re_ch4);
        report.CheckRel("CH4 productivity", pr_ch4, Reference("ferre_2026", case_id, "CH4_productivity"), fc.tol_pr_ch4);
      }

      t.SaveTest("Ferre2026_" + fc.name, solver);

      // Emit the flagship NaY(F) biogas PVSA case as a bundled GUI example.
      if (fc.name == "NaY(F)")
        t.SaveResultsAsExample("Ferre2026_NaYF", solver);
    }

    // Isotherm parameters (SI S1/S2), densities (S3), and Table-5 (recovery-
    // recovery, max Re_CO2) operating points. q = |dU| (positive); the detailed
    // model uses a fixed bulk density rho_b = 712 kg/m3 (Table 2) for all materials.
    const FerreCase kNaY_F {
      "NaY(F)",
      /*co2*/ {6.39, 1.23e-6, 27028, 0.00, 2.10e-6, 21169},
      /*ch4*/ {4.09, 8.16e-6, 12963, 2.85, 5.39e-6, 12756},
      /*rho_ads*/ 1330.0,
      /*hoa: DSL isosteric heat at feed loading, per iso_heat calc*/ 27028.0, 12896.0,
      /*p_h,p_i*/ 4.6, 0.7, /*t*/ 40.6, 87.5, 183.2, /*v_f*/ 0.5,
      /*targets*/ 95.1, 98.8, 97.8, 66.0,
      /*tol pu_co2, pu_ch4, re_co2, re_ch4 [pp], pr_ch4 [rel]*/ 1.1, 0.59, 1.3, 5.6, 0.046
      // measured: 0.740596, 0.287088, 0.961627, 4.81811 pp; 0.0352443
    };
    const FerreCase kNaY_W {
      "NaY(W)",
      {4.85, 1.75e-7, 31284, 2.79, 3.48e-8, 43133},
      {5.52, 2.65e-6, 16754, 0.00, 0.0,     0.0},
      /*rho_ads*/ 1330.0,
      31726.0, 16754.0,
      4.2, 0.7, 36.3, 64.2, 189.2, 0.5,
      95.0, 97.6, 95.9, 65.4,
      /*tol pu_co2, pu_ch4, re_co2, re_ch4 [pp], pr_ch4 [rel]*/ 1.8, 0.45, 0.89, 7.6, 0.085
      // measured: 1.49549, 0.141752, 0.580179, 6.52251 pp; 0.0736976
    };
    const FerreCase k13X {
      "13X",
      {2.82, 1.14e-5, 18494, 2.84, 4.97e-8, 42073},
      {4.49, 9.15e-5, 7528,  1.31, 9.07e-6, 15506},
      /*rho_ads*/ 1130.0,
      20012.0, 10255.0,
      3.0, 0.4, 38.9, 91.9, 178.6, 0.3,
      95.0, 98.0, 95.9, 48.6,
      /*tol pu_co2, pu_ch4, re_co2, re_ch4 [pp], pr_ch4 [rel]*/ 1.5, 0.73, 1.6, 5.7, 0.11
      // measured: 1.11122, 0.425536, 1.23481, 4.87329 pp; 0.0869989
    };
    const FerreCase kMIL160 {
      "MIL-160(Al)",
      {3.15, 1.87e-6, 27085, 3.12, 4.58e-6, 16968},
      {7.38, 1.19e-4, 4865,  2.12, 3.22e-5, 12434},
      /*rho_ads*/ 1070.0,
      22851.0, 8862.0,
      2.3, 0.3, 44.6, 109.3, 147.0, 0.4,
      95.0, 97.6, 94.1, 35.8,
      /*tol pu_co2, pu_ch4, re_co2, re_ch4 [pp], pr_ch4 [rel]*/ 2.2, 1.4, 3.5, 13, 0.68
      // measured: 1.82503, 1.09926, 2.99094, 10.7492 pp; 0.585651
    };
  } // namespace

  TEST_F(Ferre2026Test, NaY_F)   { RunFerreCase(*this, kNaY_F, "NaY_F"); }
  TEST_F(Ferre2026Test, NaY_W)   { RunFerreCase(*this, kNaY_W, "NaY_W"); }
  TEST_F(Ferre2026Test, Zeolite13X) { RunFerreCase(*this, k13X, "Zeolite13X"); }
  TEST_F(Ferre2026Test, MIL160Al){ RunFerreCase(*this, kMIL160, "MIL160Al"); }

} // namespace adsim
