// STL includes
#include <filesystem>

// ads includes
#include "adsim/equations/isotherm.h"

// Gtest includes
#include <gtest/gtest.h>

#include "gtest/gtest.h"
#include "adsim/equations/isotherm.h"
#include "adsim/math/constants.h"


namespace adsim {

  // Henry isotherm: qi = k * P * yi * exp(Q / (R * T))
  TEST(IsothermFunctionTests, Henry) {
    HenryIsothermParameters params{ 1.0, 2000.0 }; // k, Q
    double T = 300.0;  // K
    double P = 2.0;    // bar (or whatever units)
    double yi = 0.5;
    double qi = 0.0;

    IsothermFunction::Henry(qi, T, P, yi, params);

    double expected = params.k * P * yi * std::exp(params.Q / (constants::kIdealGasConstant * T));
    EXPECT_NEAR(qi, expected, 1e-12);
  }

  // Langmuir isotherm: qi = m1*b1*exp(q1/RT)*P*yi / (1 + b1*exp(q1/RT)*P*yi)
  TEST(IsothermFunctionTests, Langmuir) {
    LangmuirIsothermParameters params{ 2.0, 0.2, 1000.0 }; // m1, b1, q1
    double T = 320.0;
    double P = 1.2;
    double yi = 0.8;
    double qi = 0.0;

    IsothermFunction::Langmuir(qi, T, P, yi, params);

    double expq = std::exp(params.q1 / (constants::kIdealGasConstant * T));
    double numerator = params.m1 * params.b1 * expq * P * yi;
    double denominator = 1 + params.b1 * expq * P * yi;
    double expected = numerator / denominator;
    EXPECT_NEAR(qi, expected, 1e-12);
  }

  // AntiLangmuir: qi = a*P*yi / (1 - b*P*yi)
  TEST(IsothermFunctionTests, AntiLangmuir) {
    AntiLangmuirIsothermParameters params{ 1.5, 0.1 }; // a, b
    double T = 300.0;
    double P = 0.8;
    double yi = 0.7;
    double qi = 0.0;

    IsothermFunction::AntiLangmuir(qi, T, P, yi, params);

    double expected = params.a * P * yi / (1 - params.b * P * yi);
    EXPECT_NEAR(qi, expected, 1e-12);
  }

  // Freundlich: qi = pow(a*P*yi, 1/v)
  TEST(IsothermFunctionTests, Freundlich) {
    FreundlichIsothermParameters params{ 3.0, 2.0 }; // a, v
    double T = 350.0;
    double P = 0.9;
    double yi = 0.4;
    double qi = 0.0;

    IsothermFunction::Freundlich(qi, T, P, yi, params);

    double expected = std::pow(params.a * P * yi, 1.0 / params.v);
    EXPECT_NEAR(qi, expected, 1e-12);
  }

  // Sips: qi = q_sat * z / (1 + z), z = pow(b*P*yi, 1/v)
  TEST(IsothermFunctionTests, Sips) {
    SipsIsothermParameters params{ 5.0, 0.8, 2.0 }; // q_sat, b, v
    double T = 298.0;
    double P = 0.5;
    double yi = 0.6;
    double qi = 0.0;

    IsothermFunction::Sips(qi, T, P, yi, params);

    double z = std::pow(params.b * P * yi, 1.0 / params.v);
    double expected = params.q_sat * z / (1 + z);
    EXPECT_NEAR(qi, expected, 1e-12);
  }

  TEST(IsothermFunctionTests, Sips2) {
    Sips2IsothermParameters params;
    double T = 298.0;
    double P = 0.5;
    double yi = 0.6;
    double qi = 0.0;

    IsothermFunction::Sips2(qi, T, P, yi, params);
  }

  // DualSiteSips
  TEST(IsothermFunctionTests, DualSiteSips) {
    DualSiteSipsIsothermParameters params{ 4.0, 0.7, 1.5, 2.5, 0.6, 2.0 }; // q_sat1, b1, v1, q_sat2, b2, v2
    double T = 340.0;
    double P = 1.1;
    double yi = 0.4;
    double qi = 0.0;

    IsothermFunction::DualSiteSips(qi, T, P, yi, params);

    double z1 = std::pow(params.b1 * P * yi, 1.0 / params.v1);
    double z2 = std::pow(params.b2 * P * yi, 1.0 / params.v2);
    double expected = params.q_sat1 * z1 / (1 + z1) + params.q_sat2 * z2 / (1 + z2);
    EXPECT_NEAR(qi, expected, 1e-12);
  }

  // RedlichPeterson: qi = a*P*yi / (1 + pow(b*P*yi, v))
  TEST(IsothermFunctionTests, RedlichPeterson) {
    RedlichPetersonIsothermParameters params{ 2.5, 0.4, 1.2 }; // a, b, v
    double T = 315.0;
    double P = 1.5;
    double yi = 0.5;
    double qi = 0.0;

    IsothermFunction::RedlichPeterson(qi, T, P, yi, params);

    double z = std::pow(params.b * P * yi, params.v);
    double expected = params.a * P * yi / (1 + z);
    EXPECT_NEAR(qi, expected, 1e-12);
  }

  // Toth: qi = q_sat * bp / pow(1 + pow(bp, v), 1/v)
  TEST(IsothermFunctionTests, Toth) {
    TothIsothermParameters params{ 3.5, 0.5, 1.3 }; // q_sat, b, v
    double T = 330.0;
    double P = 0.7;
    double yi = 0.6;
    double qi = 0.0;

    IsothermFunction::Toth(qi, T, P, yi, params);

    double bp = params.b * P * yi;
    double bpv = std::pow(bp, params.v);
    double denom = std::pow(1 + bpv, 1.0 / params.v);
    double expected = params.q_sat * bp / denom;
    EXPECT_NEAR(qi, expected, 1e-12);
  }

  // DualSiteLangmuir
  TEST(IsothermFunctionTests, DualSiteLangmuir) {
    DualSiteLangmuirIsothermParameters params{ 1.0, 0.4, 900.0, 2.0, 0.2, 800.0 }; // m1, b1, q1, m2, b2, q2
    double T = 310.0;
    double P = 1.0;
    double yi = 0.5;
    double qi = 0.0;

    IsothermFunction::DualSiteLangmuir(qi, T, P, yi, params);

    double expq1 = std::exp(params.q1 / (constants::kIdealGasConstant * T));
    double expq2 = std::exp(params.q2 / (constants::kIdealGasConstant * T));
    double part1 = params.m1 * params.b1 * expq1 * P * yi / (1 + params.b1 * expq1 * P * yi);
    double part2 = params.m2 * params.b2 * expq2 * P * yi / (1 + params.b2 * expq2 * P * yi);
    double expected = part1 + part2;
    EXPECT_NEAR(qi, expected, 1e-12);
  }

  // Langmuir_Concentration
  TEST(IsothermFunctionTests, Langmuir_Concentration) {
    LangmuirIsothermParametersConcentration params{ 1.1, 0.33, 600.0 };
    double T = 290.0;
    double P = 1.4;
    double yi = 0.7;
    double qi = 0.0;

    IsothermFunction::Langmuir_Concentration(qi, T, P, yi, params);

    double ci = P * yi / (constants::kIdealGasConstant * T);
    double expq = std::exp(params.q1 / (constants::kIdealGasConstant * T));
    double numerator = params.m1 * params.b1 * expq * ci;
    double denominator = 1 + params.b1 * expq * ci;
    double expected = numerator / denominator;
    EXPECT_NEAR(qi, expected, 1e-12);
  }

  // Toth2: temperature-dependent Toth
  TEST(IsothermFunctionTests, Toth2) {
    Toth2IsothermParameters params;
    params.q_sat = 3.715;
    params.x = 0.0;
    params.T_0 = 324.0;
    params.b0 = 0.027;
    params.dH0 = 108836.0;
    params.w0 = 0.28;
    params.alpha = 0.446;

    double qi = 0.0;
    IsothermFunction::Toth2(qi, 393.15, 99784.005, 1.0, params);
    EXPECT_NEAR(qi, 0.7763031, 1e-3);

    IsothermFunction::Toth2(qi, 298.15, 672.366, 1.0, params);
    EXPECT_NEAR(qi, 1.667238, 1e-3);
  }

  // DualSiteLangmuir_Concentration
  TEST(IsothermFunctionTests, DualSiteLangmuir_Concentration) {
    DualSiteLangmuirIsothermParametersConcentration params{ 1.1, 0.33, 600.0, 0.7, 0.17, 300.0 };
    double T = 310.0;
    double P = 1.8;
    double yi = 0.6;
    double qi = 0.0;

    IsothermFunction::DualSiteLangmuir_Concentration(qi, T, P, yi, params);

    double ci = P * yi / (constants::kIdealGasConstant * T);
    double expq1 = std::exp(params.q1 / (constants::kIdealGasConstant * T));
    double expq2 = std::exp(params.q2 / (constants::kIdealGasConstant * T));
    double part1 = params.m1 * params.b1 * expq1 * ci / (1 + params.b1 * expq1 * ci);
    double part2 = params.m2 * params.b2 * expq2 * ci / (1 + params.b2 * expq2 * ci);
    double expected = part1 + part2;
    EXPECT_NEAR(qi, expected, 1e-12);
  }

}
