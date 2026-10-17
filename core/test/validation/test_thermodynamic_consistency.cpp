/*
 * Skarstrom - pressure and temperature swing adsorption simulator
 * Copyright (c) 2024-2026 Sam Affleck and contributors
 *
 * Licensed under the terms in the LICENSE file at the repository root.
 */

// Thermodynamic consistency of the isotherm library.
//
// The energy balance takes the heat of adsorption as a separate user input from
// the isotherm parameters. Nothing in the solver forces the two to agree, yet
// thermodynamics does: the heat released on adsorption is fixed by the
// temperature dependence of the isotherm through the Clausius-Clapeyron
// relation,
//
//     q_st = R T^2 ( d ln P / d T )_q
//
// These tests evaluate that derivative numerically from the isotherm functions
// themselves and check it against the parameter a user would naturally reach for
// when filling in the heat of adsorption. They need no solver run, so they are
// cheap enough to cover every isotherm form.

#include "validation_helpers.h"

#include <functional>

namespace adsim {

  using namespace adsim::validation;

  namespace {

    // Pressure at which the isotherm returns `loading` at temperature T, found
    // by bisection. `isotherm` must be monotonically increasing in P.
    double PressureAtLoading(const std::function<double(double, double)>& isotherm,
                             double loading, double T) {
      double lo = 1.0e-3, hi = 1.0e9;
      for (int i = 0; i < 200; ++i) {
        const double mid = std::sqrt(lo * hi);  // geometric bisection: wide range
        if (isotherm(mid, T) < loading) lo = mid; else hi = mid;
      }
      return std::sqrt(lo * hi);
    }

    // Isosteric heat of adsorption [J/mol] from the Clausius-Clapeyron relation,
    // evaluated by central differences in temperature at fixed loading.
    double IsostericHeat(const std::function<double(double, double)>& isotherm,
                         double loading, double T, double dT = 0.05) {
      const double P_hot = PressureAtLoading(isotherm, loading, T + dT);
      const double P_cold = PressureAtLoading(isotherm, loading, T - dT);
      return kR * T * T * (std::log(P_hot) - std::log(P_cold)) / (2.0 * dT);
    }

  } // namespace

  // ===========================================================================
  // Pressure-basis isotherms: q_st equals the fitted energy parameter
  // ===========================================================================

  TEST(ThermodynamicConsistency, HenryIsostericHeatEqualsQParameter) {
    HenryIsothermParameters params;
    params.k = 1.0e-9;
    params.Q = 28000.0;  // [J/mol]

    auto henry = [&](double P, double T) {
      double q = 0.0;
      IsothermFunction::Henry(q, T, P, 1.0, params);
      return q;
    };

    // A Henry isotherm is linear, so q_st is independent of loading.
    for (double T : {273.0, 300.0, 350.0}) {
      for (double loading : {0.01, 0.1, 1.0}) {
        EXPECT_NEAR(IsostericHeat(henry, loading, T), 28000.0, 1.0)
          << "Henry isosteric heat != Q at T = " << T << " K, q = " << loading;
      }
    }
  }

  TEST(ThermodynamicConsistency, LangmuirIsostericHeatEqualsQ1AtAllCoverages) {
    LangmuirIsothermParameters params;
    params.m1 = 5.0;         // [mol/kg]
    params.b1 = 1.0e-10;     // [1/Pa]
    params.q1 = 32000.0;     // [J/mol]

    auto langmuir = [&](double P, double T) {
      double q = 0.0;
      IsothermFunction::Langmuir(q, T, P, 1.0, params);
      return q;
    };

    // For q = m b(T) p / (1 + b(T) p) with b = b0 exp(q1/RT), the group b*p is
    // fixed at fixed loading, so ln p = const - q1/(R T) exactly: q_st = q1 at
    // every coverage, not just in the Henry limit.
    for (double coverage : {0.05, 0.5, 0.9}) {
      const double loading = coverage * 5.0;
      EXPECT_NEAR(IsostericHeat(langmuir, loading, 300.0), 32000.0, 5.0)
        << "Langmuir isosteric heat != q1 at coverage " << coverage;
    }
  }

  TEST(ThermodynamicConsistency, DualSiteLangmuirIsostericHeatIsBoundedBySiteEnergies) {
    DualSiteLangmuirIsothermParameters params;
    params.m1 = 3.0;
    params.b1 = 1.0e-10;
    params.q1 = 40000.0;  // strong site
    params.m2 = 4.0;
    params.b2 = 1.0e-12;
    params.q2 = 25000.0;  // weak site

    auto dsl = [&](double P, double T) {
      double q = 0.0;
      IsothermFunction::DualSiteLangmuir(q, T, P, 1.0, params);
      return q;
    };

    // A dual-site isotherm has a loading-dependent isosteric heat: it starts at
    // the strong site's energy and decays towards the weak site's as the strong
    // site fills. It must always stay between the two.
    const double q_low = IsostericHeat(dsl, 0.01, 300.0);
    const double q_high = IsostericHeat(dsl, 5.0, 300.0);

    EXPECT_GT(q_low, 25000.0);
    EXPECT_LT(q_low, 40000.0 + 1.0);
    EXPECT_GT(q_high, 25000.0 - 1.0);
    EXPECT_LT(q_high, 40000.0);
    EXPECT_GT(q_low, q_high)
      << "Isosteric heat should fall as the strong site saturates";

    // In the zero-coverage limit the strong site dominates the Henry constant,
    // so q_st tends to the site energy weighted by the site Henry constants.
    const double KH1 = params.m1 * params.b1 * std::exp(params.q1 / (kR * 300.0));
    const double KH2 = params.m2 * params.b2 * std::exp(params.q2 / (kR * 300.0));
    const double q_st_limit = (KH1 * params.q1 + KH2 * params.q2) / (KH1 + KH2);
    EXPECT_NEAR(IsostericHeat(dsl, 1e-4, 300.0), q_st_limit, 10.0)
      << "Zero-coverage isosteric heat does not match the Henry-weighted mean";
  }

  // ===========================================================================
  // Concentration-basis isotherms carry an extra R T
  // ===========================================================================
  //
  // This one is a trap worth codifying. For a concentration-form Langmuir,
  //
  //     q = m b(T) c / (1 + b(T) c),   c = p / (R T),  b = b0 exp(q1/RT)
  //
  // fixed loading fixes b*c, so c ~ exp(-q1/RT) and p = c R T ~ T exp(-q1/RT).
  // Differentiating, q_st = q1 + R T -- NOT q1. At 300 K that is an extra
  // 2.5 kJ/mol, roughly 8 % of a typical CO2 heat of adsorption, which lands
  // straight in the energy balance if the heat of adsorption is filled in with
  // the fitted q1.

  TEST(ThermodynamicConsistency, ConcentrationLangmuirIsostericHeatIsQ1PlusRT) {
    LangmuirIsothermParametersConcentration params;
    params.m1 = 5.0;
    params.b1 = 1.0e-6;
    params.q1 = 30000.0;

    auto langmuir_c = [&](double P, double T) {
      double q = 0.0;
      IsothermFunction::Langmuir_Concentration(q, T, P, 1.0, params);
      return q;
    };

    for (double T : {280.0, 300.0, 340.0}) {
      const double expected = 30000.0 + kR * T;
      EXPECT_NEAR(IsostericHeat(langmuir_c, 1.0, T), expected, 10.0)
        << "Concentration-basis Langmuir isosteric heat should be q1 + RT at T = " << T;
    }
  }

  // ===========================================================================
  // Low-pressure limits
  // ===========================================================================
  //
  // Every physically admissible isotherm must become linear in pressure as
  // p -> 0, with a Henry constant that is analytic in the fitted parameters.
  // A form that does not (Freundlich, for instance) cannot be used near zero
  // partial pressure without care -- worth having pinned explicitly.

  TEST(ThermodynamicConsistency, LangmuirApproachesHenryLawAtLowPressure) {
    LangmuirIsothermParameters params;
    params.m1 = 5.0;
    params.b1 = 1.0e-10;
    params.q1 = 32000.0;

    const double T = 300.0;
    const double KH = params.m1 * params.b1 * std::exp(params.q1 / (kR * T));  // [mol/(kg-Pa)]

    // The deviation from linearity is exactly b*p to leading order, so tolerance
    // is set from the parameters rather than picked arbitrarily. That makes this
    // a check on the *rate* of approach, not just on approach.
    const double b = params.b1 * std::exp(params.q1 / (kR * T));
    for (double P : {1.0, 10.0, 100.0}) {
      double q = 0.0;
      IsothermFunction::Langmuir(q, T, P, 1.0, params);
      EXPECT_NEAR(q, KH * P, 2.0 * b * P * KH * P)
        << "Langmuir does not reduce to Henry's law at P = " << P << " Pa";
    }
  }

  TEST(ThermodynamicConsistency, DualSiteLangmuirHenryConstantIsSumOfSites) {
    DualSiteLangmuirIsothermParameters params;
    params.m1 = 3.0;
    params.b1 = 1.0e-10;
    params.q1 = 40000.0;
    params.m2 = 4.0;
    params.b2 = 1.0e-12;
    params.q2 = 25000.0;

    const double T = 300.0;
    const double KH = params.m1 * params.b1 * std::exp(params.q1 / (kR * T)) +
                      params.m2 * params.b2 * std::exp(params.q2 / (kR * T));

    // Evaluate deep in the linear regime; the residual curvature is O(b1*p).
    const double P = 1.0e-3;
    const double b1 = params.b1 * std::exp(params.q1 / (kR * T));
    double q = 0.0;
    IsothermFunction::DualSiteLangmuir(q, T, P, 1.0, params);
    EXPECT_NEAR(q, KH * P, 2.0 * b1 * P * KH * P)
      << "Dual-site Henry constant is not the sum of the site Henry constants";
  }

  // ===========================================================================
  // Monotonicity and boundedness
  // ===========================================================================

  TEST(ThermodynamicConsistency, IsothermsAreMonotonicAndSaturate) {
    LangmuirIsothermParameters langmuir;
    langmuir.m1 = 5.0;
    langmuir.b1 = 1.0e-8;
    langmuir.q1 = 20000.0;

    double previous = -1.0;
    for (double P = 1.0; P <= 1e8; P *= 2.0) {
      double q = 0.0;
      IsothermFunction::Langmuir(q, 300.0, P, 1.0, langmuir);
      EXPECT_GE(q, previous) << "Langmuir loading decreased with pressure at P = " << P;
      EXPECT_LE(q, static_cast<double>(langmuir.m1) + 1e-9)
        << "Langmuir loading exceeded the saturation capacity at P = " << P;
      previous = q;
    }
    EXPECT_NEAR(previous, static_cast<double>(langmuir.m1), 1e-3 * langmuir.m1)
      << "Langmuir did not saturate at high pressure";

    // Loading must fall with temperature at fixed pressure (adsorption is
    // exothermic), for every form.
    double q_cold = 0.0, q_hot = 0.0;
    IsothermFunction::Langmuir(q_cold, 280.0, 1e5, 1.0, langmuir);
    IsothermFunction::Langmuir(q_hot, 360.0, 1e5, 1.0, langmuir);
    EXPECT_LT(q_hot, q_cold) << "Loading increased with temperature";
  }

} // namespace adsim
