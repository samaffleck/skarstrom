// E-09: a single-component fluid must get a non-zero molecular diffusivity
// (the Chapman-Enskog self-diffusion coefficient), otherwise every
// diffusion-based MTC model (Macropore, Macropore+Knudsen, Haynes) evaluates
// to zero and the adsorbed phase is frozen at its initial loading.

// STL includes
#include <cmath>

// Gtest includes
#include <gtest/gtest.h>

// Test helpers
#include "../fixtures/user_data_test.h"

namespace adsim {

  namespace {

    // Independent evaluation of the Chapman-Enskog self-diffusion coefficient
    // [m^2/s] for a pure gas: D_ii = 2.66e-7 T^1.5 / (P[bar] sqrt(M_i) sigma_ii^2 Omega_D),
    // with Omega_D from the Neufeld et al. (1972) correlation and M_AB = M_i.
    double ChapmanEnskogSelfDiffusion(const std::string& component, double T, double P) {
      const ComponentData& cd = Components::GetComponentData(component);
      const double M = CoolProp::Props1SI(component, "molar_mass") * 1e3; // [g/mol]
      const double Ts = T / cd.eps;
      const double omega = 1.06036 / std::pow(Ts, 0.15610)
                         + 0.19300 / std::exp(0.47635 * Ts)
                         + 1.03587 / std::exp(1.52996 * Ts)
                         + 1.76474 / std::exp(3.89411 * Ts);
      return 2.66e-7 * std::pow(T, 1.5) / ((P / 1e5) * std::sqrt(M) * cd.sigma * cd.sigma * omega);
    }

    void AddCO2Langmuir(UserDataTest& t, AdsorbentID ads_id) {
      auto iso = std::make_unique<LangmuirIsotherm>();
      iso->params_["CarbonDioxide"].m1 = 5.0;
      iso->params_["CarbonDioxide"].b1 = 3e-11;
      iso->params_["CarbonDioxide"].q1 = 35000;
      t.AddIsotherm({ "CarbonDioxide" }, ads_id, std::move(iso));
    }

  } // namespace

  enum class DiffusionMTC { Macropore, MacroporeKnudsen, Haynes };

  // Pure-CO2 bed: after initialisation Dmi must equal the Chapman-Enskog
  // self-diffusion value and the diffusion-based MTC must be finite and > 0.
  static void CheckSingleComponentDiffusion(UserDataTest& t, DiffusionMTC model) {
    const double tortuosity = 3.0;
    const double ep = 0.35;

    t.AddComponents({ "CarbonDioxide" });
    auto [layer_id, ads_id] = t.AddLayer("13X");
    AddCO2Langmuir(t, ads_id);
    t.SetHOA(ads_id, "CarbonDioxide", 0.0);
    switch (model) {
      case DiffusionMTC::Macropore:
        t.AddLinearDrivingForceModel_Macropore(ads_id, "CarbonDioxide", tortuosity, ep); break;
      case DiffusionMTC::MacroporeKnudsen:
        t.AddLinearDrivingForceModel_MacroporeKnudsen(ads_id, "CarbonDioxide", 1e-7, tortuosity, ep); break;
      case DiffusionMTC::Haynes:
        t.AddLinearDrivingForceModel_Haynes(ads_id, "CarbonDioxide", tortuosity, ep); break;
    }
    t.AddConstantInternalHeatTransferCoefficientModel(0, layer_id);
    t.AddConstantExternalHeatTransferCoefficientModel(0, layer_id);
    t.AddInitialConditions(2e5, 310, 310, { {"CarbonDioxide", 1.0} });

    t.data.Initialise();

    const auto& state = t.data.state;
    const auto& P = state.pressure.var;
    const auto& T = state.fluid_temperature.var;
    const auto& Dm = state.molecular_diffusion.at("CarbonDioxide").var;
    const auto& mtc = state.mass_transfer_coefficient.at(ads_id).at("CarbonDioxide").var;
    const auto& layer = t.GetLayer(layer_id);
    const double dp = layer.particle_diameter_;

    for (int n = layer.domain_.start_index; n <= layer.domain_.end_index; ++n) {
      const double D_ce = ChapmanEnskogSelfDiffusion("CarbonDioxide", T[n], P[n]);
      ASSERT_GT(D_ce, 1e-6);  // CO2 self-diffusion at 2 bar, 310 K is ~5e-6 m^2/s
      EXPECT_NEAR(Dm[n], D_ce, 1e-12 * D_ce) << "cell " << n;

      EXPECT_TRUE(std::isfinite(mtc[n])) << "cell " << n;
      EXPECT_GT(mtc[n], 0.0) << "cell " << n;
      if (model == DiffusionMTC::Macropore) {
        EXPECT_NEAR(mtc[n], 60.0 * Dm[n] * ep / tortuosity / (dp * dp), 1e-12 * mtc[n]) << "cell " << n;
      }
    }

    // Schmidt number must come from the real diffusivity, not the D_min floor.
    const auto& Sc = state.schmidt_number.var;
    const auto& vis = state.viscosity.var;
    const auto& rho = state.fluid_density.var;
    const int n_mid = (layer.domain_.start_index + layer.domain_.end_index) / 2;
    EXPECT_NEAR(Sc[n_mid], vis[n_mid] / (rho[n_mid] * Dm[n_mid]), 1e-9 * Sc[n_mid]);
  }

  TEST_F(UserDataTest, SingleComponentDiffusion_Macropore) {
    CheckSingleComponentDiffusion(*this, DiffusionMTC::Macropore);
  }

  TEST_F(UserDataTest, SingleComponentDiffusion_MacroporeKnudsen) {
    CheckSingleComponentDiffusion(*this, DiffusionMTC::MacroporeKnudsen);
  }

  TEST_F(UserDataTest, SingleComponentDiffusion_Haynes) {
    CheckSingleComponentDiffusion(*this, DiffusionMTC::Haynes);
  }

  // Pressurise a pure-CO2 bed with a Macropore MTC model: the loading must rise
  // from its initial (equilibrium at the initial pressure) value towards the
  // new equilibrium. Before the fix mtc == 0 and the loading never moved.
  TEST_F(UserDataTest, SingleComponentDiffusion_MacroporeBedAdsorbs) {
    AddComponents({ "CarbonDioxide" });
    auto [layer_id, ads_id] = AddLayer("13X");
    AddCO2Langmuir(*this, ads_id);
    SetHOA(ads_id, "CarbonDioxide", 0.0);
    AddLinearDrivingForceModel_Macropore(ads_id, "CarbonDioxide", 3.0, 0.35);
    AddConstantInternalHeatTransferCoefficientModel(0, layer_id);
    AddConstantExternalHeatTransferCoefficientModel(0, layer_id);
    AddInitialConditions(1e5, 298, 298, { {"CarbonDioxide", 1.0} });

    AddStep("Pressurise", 100);
    Step& step = data.cycle.GetStep("Pressurise");
    auto inlet = CreatePressureInletBoundary(
      CreateLinearBoundaryValue(1e5, 3e5, 10.0), { {"CarbonDioxide", 1.0} }, 298);
    inlet->location_ = BoundaryLocation::BOTTOM;
    inlet->direction_ = BoundaryDirection::INLET;
    step.bottom_boundary_condition_ = std::move(inlet);
    AddClosedBC(BoundaryLocation::TOP, "Pressurise");
    ConfigureCycle(1);

    std::string error{};
    ASSERT_TRUE(data.PerformChecks(error)) << error;

    CoupledSolver solver{data};
    int retval = solver.IntegrateCycle(data.cycle, false);
    ASSERT_GE(retval, 0) << "Solver did not complete";

    const auto& layer = GetLayer(layer_id);
    const auto& q = data.state.ads_concentration.at(ads_id).at("CarbonDioxide").var;
    const auto& q_eq = data.state.eq_ads_concentration.at(ads_id).at("CarbonDioxide").var;
    const auto& mtc = data.state.mass_transfer_coefficient.at(ads_id).at("CarbonDioxide").var;

    // Equilibrium loading at the initial state (1 bar, 298 K, pure CO2).
    LangmuirIsothermParameters p0;
    p0.m1 = 5.0; p0.b1 = 3e-11; p0.q1 = 35000;
    double q0 = 0.0;
    IsothermFunction::Langmuir(q0, 298.0, 1e5, 1.0, p0);

    for (int n = layer.domain_.start_index; n <= layer.domain_.end_index; ++n) {
      EXPECT_GT(mtc[n], 0.0) << "cell " << n;
      EXPECT_GT(q[n], q0 + 0.1) << "loading did not rise in cell " << n;
      // Fast macropore diffusion: after 100 s the loading is close to equilibrium.
      EXPECT_NEAR(q[n], q_eq[n], 0.05 * q_eq[n]) << "cell " << n;
    }
  }

} // namespace adsim
