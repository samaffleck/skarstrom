// adsim includes
#include "../src/pch.h"
#include "adsim/objects/wall.h"
#include "adsim/objects/wall_design.h"

// Gtest includes
#include <gtest/gtest.h>

#include <cmath>


namespace adsim {

  // ---------------------------------------------------------------------------
  // GetMaxAllowableStress
  // ---------------------------------------------------------------------------

  TEST(WallAllowableStressTest, SS316L_AtRoomTemperature) {
    // ASME IID Table 1A: SS316L at 38°C (100°F) = 115.1 MPa
    const double S = GetMaxAllowableStress(WallMaterial::SS316L, 311.15); // 38°C in K
    EXPECT_NEAR(S, 115.1e6, 1e3);
  }

  TEST(WallAllowableStressTest, SS316L_At300C) {
    // ASME IID Table 1A: SS316L at 149°C ≈ 104.8 MPa, at 204°C ≈ 101.4 MPa
    // At 300°C, interpolate between 260°C (97.2 MPa) and 316°C (93.1 MPa)
    // frac = (300 - 260) / (316 - 260) = 40/56 ≈ 0.7143
    // expected ≈ 97.2e6 + 0.7143 * (93.1e6 - 97.2e6) ≈ 97.2e6 - 2.929e6 ≈ 94.27e6
    const double S = GetMaxAllowableStress(WallMaterial::SS316L, 573.15); // 300°C in K
    EXPECT_GT(S, 93.1e6);
    EXPECT_LT(S, 97.2e6);
  }

  TEST(WallAllowableStressTest, SS316L_At500C) {
    // ASME IID Table 1A: SS316L at 482°C ≈ 62.1 MPa, at 538°C ≈ 48.3 MPa
    // At 500°C, interpolate between 482°C and 538°C
    const double S = GetMaxAllowableStress(WallMaterial::SS316L, 773.15); // 500°C in K
    EXPECT_GT(S, 48.3e6);
    EXPECT_LT(S, 62.1e6);
  }

  TEST(WallAllowableStressTest, SS316L_AboveMaxTemperatureClamps) {
    // Above 566°C, should return the minimum (41.4 MPa)
    const double S = GetMaxAllowableStress(WallMaterial::SS316L, 1273.15); // 1000°C
    EXPECT_NEAR(S, 41.4e6, 1e3);
  }

  TEST(WallAllowableStressTest, SS316L_BelowMinTemperatureClamps) {
    // Below -29°C, should return the maximum (115.1 MPa)
    const double S = GetMaxAllowableStress(WallMaterial::SS316L, 173.15); // -100°C
    EXPECT_NEAR(S, 115.1e6, 1e3);
  }

  TEST(WallAllowableStressTest, SS304_AtRoomTemperature) {
    // ASME IID: SS304 at RT = 129.6 MPa
    const double S = GetMaxAllowableStress(WallMaterial::SS304, 293.15);
    EXPECT_NEAR(S, 129.6e6, 1e3);
  }

  TEST(WallAllowableStressTest, CarbonSteelA106_AtRoomTemperature) {
    // SA-106 Gr.B: constant 120.7 MPa up to ~316°C
    const double S = GetMaxAllowableStress(WallMaterial::CarbonSteel_A106, 293.15);
    EXPECT_NEAR(S, 120.7e6, 1e3);
  }

  TEST(WallAllowableStressTest, CarbonSteelA106_At350C) {
    // SA-106 Gr.B: constant 120.7 MPa up to 316°C, then drops
    // At 350°C, interpolating between 316°C (120.7 MPa) and 371°C (111.0 MPa)
    const double S = GetMaxAllowableStress(WallMaterial::CarbonSteel_A106, 623.15);
    EXPECT_GT(S, 111.0e6);
    EXPECT_LT(S, 120.7e6);
  }

  TEST(WallAllowableStressTest, Hastelloy_AtRoomTemperature) {
    // Hastelloy C-276 at RT = 194.4 MPa
    const double S = GetMaxAllowableStress(WallMaterial::Hastelloy_C276, 293.15);
    EXPECT_NEAR(S, 194.4e6, 1e3);
  }

  TEST(WallAllowableStressTest, Aluminium6061_AtRoomTemperature) {
    const double S = GetMaxAllowableStress(WallMaterial::Aluminium_6061, 293.15);
    EXPECT_NEAR(S, 103.4e6, 1e3);
  }

  TEST(WallAllowableStressTest, Aluminium6061_At200C_IsLow) {
    // Al 6061-T6 loses strength rapidly above ~150°C
    const double S_RT  = GetMaxAllowableStress(WallMaterial::Aluminium_6061, 293.15);
    const double S_200 = GetMaxAllowableStress(WallMaterial::Aluminium_6061, 473.15);
    EXPECT_LT(S_200, S_RT / 5.0); // should be much lower
  }

  TEST(WallAllowableStressTest, TitaniumGr2_AtRoomTemperature) {
    const double S = GetMaxAllowableStress(WallMaterial::Titanium_Gr2, 293.15);
    EXPECT_NEAR(S, 98.6e6, 1e3);
  }

  TEST(WallAllowableStressTest, CustomMaterialReturnsZero) {
    const double S = GetMaxAllowableStress(WallMaterial::Custom, 293.15);
    EXPECT_DOUBLE_EQ(S, 0.0);
  }

  TEST(WallAllowableStressTest, InterpolationIsMonotonicallyDecreasing_SS316L) {
    // Allowable stress must not increase with temperature for SS316L
    double prev = GetMaxAllowableStress(WallMaterial::SS316L, 293.15);
    for (int T_C = 50; T_C <= 550; T_C += 50) {
      const double S = GetMaxAllowableStress(WallMaterial::SS316L, T_C + 273.15);
      EXPECT_LE(S, prev + 1.0); // allow 1 Pa rounding tolerance
      prev = S;
    }
  }

  // ---------------------------------------------------------------------------
  // DesignWall — internal pressure
  // ---------------------------------------------------------------------------

  // 100 mm ID SS316L column at 20 °C, E = 1: S = 115.1 MPa.
  WallDesignInput Column(double P, double c = 0.0) {
    WallDesignInput in;
    in.material = WallMaterial::SS316L;
    in.inner_diameter = 0.1;
    in.length = 1.0;
    in.design_pressure = P;
    in.design_temperature = 293.15;
    in.joint_efficiency = 1.0;
    in.corrosion_allowance = c;
    return in;
  }

  // UG-27(c)(1): t = P R / (S E - 0.6 P) = 1e6 * 0.05 / (115.1e6 - 0.6e6) = 0.4367 mm.
  TEST(WallDesignTest, Ug27HandCalculation) {
    const auto r = DesignWall(Column(10e5));
    ASSERT_TRUE(r.valid) << r.message;
    EXPECT_NEAR(r.allowable_stress, 115.1e6, 1e3);
    EXPECT_NEAR(r.required_thickness, 1e6 * 0.05 / (115.1e6 - 0.6e6), 1e-12);
    EXPECT_FALSE(r.thick_shell);
    EXPECT_FALSE(r.external_governs);
    EXPECT_DOUBLE_EQ(r.external_thickness, 0.0);
  }

  // The corrosion allowance enlarges the radius and is added to the thickness:
  // t = 1e6 * (0.05 + 0.0015) / 114.5e6 + 0.0015.
  TEST(WallDesignTest, CorrosionAllowance) {
    const auto r = DesignWall(Column(10e5, 0.0015));
    ASSERT_TRUE(r.valid) << r.message;
    EXPECT_NEAR(r.required_thickness, 1e6 * 0.0515 / 114.5e6 + 0.0015, 1e-12);
  }

  // Above P = 0.385 S E the thick-shell rule of Appendix 1-2 applies:
  // t = R (sqrt(Z) - 1), Z = (S E + P) / (S E - P).
  TEST(WallDesignTest, Appendix1_2ThickShell) {
    const double P = 50e6;
    const auto r = DesignWall(Column(P));
    ASSERT_TRUE(r.valid) << r.message;
    EXPECT_TRUE(r.thick_shell);
    const double Z = (115.1e6 + P) / (115.1e6 - P);
    EXPECT_NEAR(r.required_thickness, 0.05 * (std::sqrt(Z) - 1.0), 1e-12);
  }

  TEST(WallDesignTest, TrendsWithPressureTemperatureAndJointEfficiency) {
    const double t0 = DesignWall(Column(10e5)).required_thickness;
    EXPECT_GT(DesignWall(Column(50e5)).required_thickness, t0);
    auto hot = Column(10e5);
    hot.design_temperature = 773.15;
    EXPECT_GT(DesignWall(hot).required_thickness, t0);
    auto spot = Column(10e5);
    spot.joint_efficiency = 0.85;
    EXPECT_GT(DesignWall(spot).required_thickness, t0);
  }

  // B31.3 eq. (3b) with Y = 0.4 equals UG-27 for the inside diameter...
  TEST(WallDesignTest, B31_3EqualsUg27InTheThinRange) {
    auto in = Column(10e5, 0.001);
    in.code = WallDesignCode::ASME_B31_3;
    const auto r = DesignWall(in);
    ASSERT_TRUE(r.valid) << r.message;
    EXPECT_NEAR(r.required_thickness, DesignWall(Column(10e5, 0.001)).required_thickness, 1e-12);
  }

  // ...but has no thick-wall rule: past P / (S E) = 0.385 it needs special consideration.
  TEST(WallDesignTest, B31_3ThickWallIsNotValid) {
    auto in = Column(50e6);
    in.code = WallDesignCode::ASME_B31_3;
    const auto r = DesignWall(in);
    EXPECT_FALSE(r.valid);
    EXPECT_NE(r.message.find("304.1.2(b)"), std::string::npos);
  }

  // The MAWP of the required wall is the design pressure, for both codes.
  TEST(WallDesignTest, MaxAllowableInternalPressureInvertsTheDesign) {
    for (auto code : { WallDesignCode::ASME_VIII_1, WallDesignCode::ASME_B31_3 }) {
      for (double P : { 5e5, 40e5, 200e5 }) {
        auto in = Column(P, 0.001);
        in.code = code;
        const auto r = DesignWall(in);
        ASSERT_TRUE(r.valid) << r.message;
        EXPECT_NEAR(MaxAllowableInternalPressure(in, r.required_thickness), P, 1e-6 * P);
      }
    }
    auto thick = Column(50e6);
    EXPECT_NEAR(MaxAllowableInternalPressure(thick, DesignWall(thick).required_thickness), 50e6, 1.0);
  }

  TEST(WallDesignTest, NotValidWithoutAKnownStressOrADesign) {
    auto custom = Column(10e5);
    custom.material = WallMaterial::Custom;
    EXPECT_FALSE(DesignWall(custom).valid);
    EXPECT_FALSE(DesignWall(Column(0.0)).valid);
    EXPECT_FALSE(DesignWall(Column(200e6)).valid);  // P > S E
  }

  // ---------------------------------------------------------------------------
  // External pressure (Code Case 2286 / VIII-2 4.4.5.1)
  // ---------------------------------------------------------------------------

  // A long tube in the elastic range: Ch = 0.55 t/Do, so Fhe = 0.88 E (t/Do)^2,
  // FS = 2 and Pa = 2 (Fhe / 2)(t/Do) = 0.88 E (t/Do)^3 (the classical long-tube
  // collapse pressure 2E/(1-v^2) (t/D)^3 with a 0.8 knock-down, halved).
  TEST(WallExternalPressureTest, LongTubeElasticClosedForm) {
    WallDesignInput in = Column(0.0);
    in.inner_diameter = 0.5;
    in.length = 10.0;
    const double E = GetElasticModulus(WallMaterial::SS316L, 293.15);
    EXPECT_NEAR(E, 195.1e9, 1e6);
    const double t = 0.004;
    const double Do = 0.5 + 2 * t;
    EXPECT_NEAR(AllowableExternalPressure(in, t), 0.88 * E * std::pow(t / Do, 3), 1e-9 * E);

    // Designing for full vacuum solves 0.88 E (t/Do)^3 = 101 325 Pa.
    in.vacuum_service = true;
    const auto r = DesignWall(in);
    ASSERT_TRUE(r.valid) << r.message;
    EXPECT_TRUE(r.external_governs);
    const double k = std::cbrt(kFullVacuumExternalPressure / (0.88 * E));  // t / Do
    EXPECT_NEAR(r.required_thickness, k * 0.5 / (1 - 2 * k), 1e-8);
  }

  // A thick, short tube buckles at yield: Fic = Sy, FS = 1.667, Pa = 2 (Sy / 1.667)(t/Do),
  // with Sy the lower bound 1.5 S for carbon steel.
  TEST(WallExternalPressureTest, YieldGovernedThickTube) {
    WallDesignInput in = Column(0.0);
    in.material = WallMaterial::CarbonSteel_A106;
    in.inner_diameter = 0.05;
    const double t = 0.005, Do = 0.06;
    const double Sy = YieldStrengthLowerBound(WallMaterial::CarbonSteel_A106, 293.15);
    EXPECT_NEAR(Sy, 1.5 * 120.7e6, 1.0);
    EXPECT_NEAR(AllowableExternalPressure(in, t), 2.0 * (Sy / 1.667) * t / Do, 1.0);
  }

  TEST(WallExternalPressureTest, ShorterAndThickerHoldsMore) {
    WallDesignInput in = Column(0.0);
    in.inner_diameter = 0.3;
    double prev = 0.0;
    for (double t : { 0.001, 0.002, 0.004, 0.008 }) {
      const double pa = AllowableExternalPressure(in, t);
      EXPECT_GT(pa, prev);
      prev = pa;
    }
    in.length = 3.0;
    const double long_pa = AllowableExternalPressure(in, 0.002);
    in.length = 0.3;
    EXPECT_GT(AllowableExternalPressure(in, 0.002), long_pa);
  }

  // The previous simplification (1 atm as an internal pressure in UG-27) gave
  // 0.13 mm for a 300 mm column; buckling needs several times that.
  TEST(WallExternalPressureTest, VacuumNeedsMoreThanUg27AtOneAtmosphere) {
    WallDesignInput in = Column(0.0);
    in.inner_diameter = 0.3;
    in.vacuum_service = true;
    const auto r = DesignWall(in);
    ASSERT_TRUE(r.valid) << r.message;
    const double ug27_1atm = 101325.0 * 0.15 / (115.1e6 - 0.6 * 101325.0);
    EXPECT_GT(r.required_thickness, 5.0 * ug27_1atm);
    EXPECT_NEAR(AllowableExternalPressure(in, r.required_thickness), kFullVacuumExternalPressure, 1.0);
  }

  // Vacuum service and a design pressure: the larger requirement governs, and
  // the corrosion allowance is kept in both.
  TEST(WallExternalPressureTest, TheLargerRequirementGoverns) {
    WallDesignInput in = Column(2e5, 0.001);
    in.inner_diameter = 0.3;
    in.vacuum_service = true;
    const auto r = DesignWall(in);
    ASSERT_TRUE(r.valid) << r.message;
    EXPECT_GT(r.internal_thickness, 0.001);
    EXPECT_GT(r.external_thickness, r.internal_thickness);
    EXPECT_DOUBLE_EQ(r.required_thickness, r.external_thickness);
    EXPECT_GE(AllowableExternalPressure(in, r.required_thickness), kFullVacuumExternalPressure);

    in.design_pressure = 50e5;
    const auto high = DesignWall(in);
    EXPECT_FALSE(high.external_governs);
    EXPECT_DOUBLE_EQ(high.required_thickness, high.internal_thickness);
  }

  // ---------------------------------------------------------------------------
  // BedWall
  // ---------------------------------------------------------------------------

  TEST(BedWallTest, DesignInUse) {
    BedWall wall;
    EXPECT_FALSE(wall.DesignInUse());
    wall.design_pressure_ = 5e5;
    EXPECT_TRUE(wall.DesignInUse());
    wall.design_pressure_ = 0.0;
    wall.vacuum_service_ = true;
    EXPECT_TRUE(wall.DesignInUse());
  }

} // end adsim namespace
