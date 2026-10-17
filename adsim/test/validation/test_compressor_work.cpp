/*
 * Skarstrom - pressure and temperature swing adsorption simulator
 * Copyright (c) 2024-2026 Sam Affleck and contributors
 *
 * Licensed under the terms in the LICENSE file at the repository root.
 */

// Analytical validation of the compression work models.
//
// Compression duty is the single largest operating cost in a PSA plant, so the
// power model feeds straight into every economic KPI. The polytropic work of a
// reversible compressor is textbook,
//
//     W = (1/eta) * (n/(n-1)) * n_dot R T1 [ (P2/P1)^((n-1)/n) - 1 ]
//
// with two limits worth pinning: n -> 1 must recover the isothermal work
// n_dot R T1 ln(P2/P1) (the bracket and the prefactor both diverge, so this is a
// genuine L'Hopital limit and an easy thing to get wrong), and n = gamma must
// give the reversible adiabatic work, which must equal the enthalpy rise
// n_dot Cp (T2 - T1) computed from the model's own outlet temperature.

#include "validation_helpers.h"

#include <limits>

namespace adsim {

  using namespace adsim::validation;

  namespace {

    // Compressor owns unique_ptr models, so configure in place rather than
    // returning by value.
    void Configure(Compressor& compressor, double suction_pressure,
                   double polytropic_index, double efficiency) {
      compressor.use_pump_ = true;
      compressor.suction_pressure_ = suction_pressure;
      compressor.poly_index_ = polytropic_index;
      auto eff = std::make_unique<ConstantPumpEfficiencyModel>();
      eff->efficiency_ = efficiency;
      compressor.efficiency_model_ = std::move(eff);
    }

  } // namespace

  TEST(CompressorWork, PolytropicPowerMatchesClosedForm) {
    const double P1 = 1.0e5;
    const double P2 = 5.0e5;
    const double T1 = 300.0;
    const double molar_flow = 0.1;  // [mol/s]
    const double eta = 0.72;
    const double n = 1.4;

    Compressor compressor;
    Configure(compressor, P1, n, eta);

    const double expected =
      (1.0 / eta) * (n / (n - 1.0)) * molar_flow * kR * T1 *
      (std::pow(P2 / P1, (n - 1.0) / n) - 1.0);

    EXPECT_NEAR(compressor.CalculatePower(P2, molar_flow, T1), expected, 1e-6 * expected);
  }

  TEST(CompressorWork, PowerTendsToIsothermalWorkAsPolytropicIndexApproachesOne) {
    const double P1 = 1.0e5;
    const double P2 = 4.0e5;
    const double T1 = 310.0;
    const double molar_flow = 0.05;
    const double eta = 1.0;  // reversible, so the comparison is exact

    // Isothermal reversible work: n_dot R T ln(P2/P1).
    const double isothermal = molar_flow * kR * T1 * std::log(P2 / P1);

    // Approach n -> 1 from above; the error is O(n-1), so successive halvings of
    // (n-1) must halve the discrepancy. Checking the *rate* of convergence is a
    // far stronger statement than checking one loose tolerance.
    double previous_error = 0.0;
    for (double delta : {0.08, 0.04, 0.02, 0.01}) {
      Compressor compressor;
      Configure(compressor, P1, 1.0 + delta, eta);
      const double error =
        std::abs(compressor.CalculatePower(P2, molar_flow, T1) - isothermal);
      if (previous_error > 0.0) {
        EXPECT_NEAR(error / previous_error, 0.5, 0.05)
          << "Polytropic power does not converge linearly to the isothermal work "
             "as n -> 1 (delta = " << delta << ")";
      }
      previous_error = error;
    }

    // And the residual at the last step must be small in absolute terms.
    EXPECT_LT(previous_error, 0.01 * isothermal);
  }

  TEST(CompressorWork, ReversibleAdiabaticPowerEqualsEnthalpyRise) {
    const double P1 = 1.0e5;
    const double P2 = 6.0e5;
    const double T1 = 295.0;
    const double molar_flow = 0.2;
    const double gamma = 1.4;  // diatomic ideal gas

    Compressor compressor;
    Configure(compressor, P1, gamma, /*efficiency=*/1.0);
    compressor.energy_eq_ = CreateCompressorEnergyModel(CompressorEnergyType::ADIABATIC);

    // For a reversible adiabatic compression the shaft work is exactly the
    // enthalpy rise, n_dot Cp (T2 - T1), with Cp = gamma R / (gamma - 1).
    const double T2 = compressor.CalculateOutletTemperature(P2, T1);
    const double Cp = gamma * kR / (gamma - 1.0);
    const double enthalpy_rise = molar_flow * Cp * (T2 - T1);

    EXPECT_NEAR(compressor.CalculatePower(P2, molar_flow, T1), enthalpy_rise,
                1e-8 * enthalpy_rise)
      << "Reversible adiabatic compression power does not equal the enthalpy rise";

    // Sanity: the isentropic outlet temperature must satisfy T2/T1 = (P2/P1)^((g-1)/g).
    EXPECT_NEAR(T2, T1 * std::pow(P2 / P1, (gamma - 1.0) / gamma), 1e-9 * T1);
  }

  TEST(CompressorWork, PowerIsMonotonicAndVanishesAtUnityPressureRatio) {
    Compressor compressor;
    Configure(compressor, 1.0e5, 1.4, 0.7);
    const double T1 = 300.0;
    const double molar_flow = 0.1;

    EXPECT_NEAR(compressor.CalculatePower(1.0e5, molar_flow, T1), 0.0, 1e-9)
      << "Compressing to the suction pressure should cost nothing";

    double previous = -1.0;
    for (double P2 : {1.5e5, 2.0e5, 4.0e5, 8.0e5}) {
      const double power = compressor.CalculatePower(P2, molar_flow, T1);
      EXPECT_GT(power, previous) << "Power is not monotonic in discharge pressure";
      previous = power;
    }

    // Expansion is not a turbine: negative work is clipped to zero.
    EXPECT_DOUBLE_EQ(compressor.CalculatePower(0.5e5, molar_flow, T1), 0.0);
  }

  TEST(CompressorWork, MolarFlowIsTheInverseOfPower) {
    Compressor compressor;
    Configure(compressor, 1.0e5, 1.35, 0.68);
    const double P2 = 3.0e5;
    const double T1 = 305.0;
    const double molar_flow = 0.15;

    const double power = compressor.CalculatePower(P2, molar_flow, T1);
    EXPECT_NEAR(compressor.CalculateMolarFlow(P2, power, T1), molar_flow,
                1e-8 * molar_flow)
      << "CalculateMolarFlow does not invert CalculatePower";
  }

  // E-19: the inversion must hold for every efficiency model, not only the
  // constant one. The pressure-dependent models used to be evaluated at the
  // discharge pressure in CalculatePower but at the suction pressure in
  // CalculateMolarFlow, so the flow was not the inverse of the power.
  namespace {

    std::unique_ptr<IPumpEfficiencyModel> MakeEfficiencyModel(PumpEfficiencyType type) {
      auto model = CreatePumpEfficiencyModel(type);
      if (auto* poly = dynamic_cast<PolyPumpEfficiencyModel*>(model.get())) {
        // eta = 0.4 + 1e-6 P (P in Pa): 0.5 at 1 bar, 0.7 at 3 bar.
        poly->a_ = 0.4;
        poly->b_ = 1e-6;
      }
      return model;
    }

  } // namespace

  TEST(CompressorWork, MolarFlowIsTheInverseOfPowerForEveryEfficiencyModel) {
    const double T1 = 305.0;
    const double molar_flow = 0.15;

    for (auto type : {PumpEfficiencyType::CONSTANT, PumpEfficiencyType::VACUUM_PUMP,
                      PumpEfficiencyType::POLYNOMIAL}) {
      SCOPED_TRACE("efficiency model " + std::to_string(static_cast<int>(type)));

      // Compressor: suction 1 bar -> bed 3 bar (explicit-suction overload too).
      Compressor compressor;
      Configure(compressor, 1.0e5, 1.35, 0.68);
      compressor.efficiency_model_ = MakeEfficiencyModel(type);
      const double P2 = 3.0e5;
      const double power = compressor.CalculatePower(P2, molar_flow, T1);
      ASSERT_GT(power, 0.0);
      EXPECT_NEAR(compressor.CalculateMolarFlow(P2, power, T1), molar_flow, 1e-10 * molar_flow);
      EXPECT_NEAR(compressor.CalculateMolarFlow(P2, power, T1, 1.0e5), molar_flow, 1e-10 * molar_flow);

      // Vacuum pump: bed 0.2 bar -> discharge 1 bar.
      VacuumPump pump;
      pump.use_pump_ = true;
      pump.poly_index_ = 1.35;
      pump.discharge_pressure_ = 1.0e5;
      pump.efficiency_model_ = MakeEfficiencyModel(type);
      const double P1 = 0.2e5;
      const double vac_power = pump.CalculatePower(P1, molar_flow, T1);
      ASSERT_GT(vac_power, 0.0);
      EXPECT_NEAR(pump.CalculateMolarFlow(P1, vac_power, T1), molar_flow, 1e-10 * molar_flow);
      EXPECT_NEAR(pump.CalculateMolarFlow(P1, vac_power, T1, 1.0e5), molar_flow, 1e-10 * molar_flow);
    }
  }

  // E-19: with the discharge at the suction pressure the head factor
  // (P2/P1)^((n-1)/n) - 1 vanishes and a fixed power used to give an infinite
  // flow (and an infinite IDA residual at t = 0 of a pressurisation step that
  // starts at the suction pressure). The flow must be finite and continuous
  // across P2 = P1, non-increasing in P2, and must not collapse to zero when
  // the discharge drops below the suction pressure.
  TEST(CompressorWork, MolarFlowIsFiniteAndContinuousAtUnityPressureRatio) {
    const double P_suction = 1.0e5;
    const double T1 = 300.0;
    const double power = 500.0;

    Compressor compressor;
    Configure(compressor, P_suction, 1.4, 0.7);

    const double n_equal = compressor.CalculateMolarFlow(P_suction, power, T1);
    ASSERT_TRUE(std::isfinite(n_equal));
    EXPECT_GT(n_equal, 0.0);

    for (double eps : {1e-12, 1e-9, 1e-6}) {
      const double n_above = compressor.CalculateMolarFlow(P_suction * (1.0 + eps), power, T1);
      const double n_below = compressor.CalculateMolarFlow(P_suction * (1.0 - eps), power, T1);
      ASSERT_TRUE(std::isfinite(n_above));
      ASSERT_TRUE(std::isfinite(n_below));
      EXPECT_NEAR(n_above, n_equal, 1e-3 * n_equal) << "eps = " << eps;
      EXPECT_NEAR(n_below, n_equal, 1e-3 * n_equal) << "eps = " << eps;
    }

    // Monotone (non-increasing) and finite over a sweep straddling P2 = P1,
    // including through the regularisation band just above it.
    double previous = std::numeric_limits<double>::infinity();
    for (double ratio = 0.5; ratio <= 1.5; ratio += 1e-4) {
      const double n = compressor.CalculateMolarFlow(ratio * P_suction, power, T1);
      ASSERT_TRUE(std::isfinite(n)) << "ratio = " << ratio;
      ASSERT_GT(n, 0.0) << "ratio = " << ratio;
      EXPECT_LE(n, previous * (1.0 + 1e-12)) << "ratio = " << ratio;
      previous = n;
    }

    // Outside the small regularisation band the inversion is exact.
    const double P2 = 1.002 * P_suction;
    const double n_ref = 0.3;
    EXPECT_NEAR(compressor.CalculateMolarFlow(P2, compressor.CalculatePower(P2, n_ref, T1), T1),
                n_ref, 1e-10 * n_ref);

    // Same singularity for a vacuum pump with the bed at the discharge pressure.
    VacuumPump pump;
    pump.use_pump_ = true;
    pump.poly_index_ = 1.4;
    pump.discharge_pressure_ = P_suction;
    const double v_equal = pump.CalculateMolarFlow(P_suction, power, T1);
    ASSERT_TRUE(std::isfinite(v_equal));
    EXPECT_GT(v_equal, 0.0);
    for (double eps : {1e-12, 1e-9, 1e-6}) {
      EXPECT_NEAR(pump.CalculateMolarFlow(P_suction * (1.0 + eps), power, T1), v_equal, 1e-3 * v_equal);
      EXPECT_NEAR(pump.CalculateMolarFlow(P_suction * (1.0 - eps), power, T1), v_equal, 1e-3 * v_equal);
    }
    // Bed above the discharge pressure: finite and not zero.
    const double v_above = pump.CalculateMolarFlow(1.2 * P_suction, power, T1);
    EXPECT_TRUE(std::isfinite(v_above));
    EXPECT_GT(v_above, 0.0);
  }

} // namespace adsim
