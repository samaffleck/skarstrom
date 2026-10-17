// adsim includes
#include "../src/pch.h"
#include "adsim/objects/wall.h"

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
  // CalculateAsmeMinWallThickness — internal pressure
  // ---------------------------------------------------------------------------
  // Hand calculation: 100 mm ID, 10 bar gauge, SS316L at 20°C, E = 1.0
  //   S = 115.1 MPa, R = 0.05 m, P = 1e6 Pa
  //   t = P*R / (S*E - 0.6*P) = 1e6 * 0.05 / (115.1e6 * 1.0 - 0.6 * 1e6)
  //     = 50000 / (115100000 - 600000)
  //     = 50000 / 114500000
  //     ≈ 4.367e-4 m  ≈ 0.437 mm

  TEST(AsmeWallThicknessTest, InternalPressure_BasicHandCalc) {
    const double P  = 10.0 * 1e5;   // 10 bar gauge → Pa
    const double D  = 0.1;          // 100 mm ID
    const double T_K = 293.15;      // 20°C
    const double E  = 1.0;

    const auto result = CalculateAsmeMinWallThickness(D, P, T_K, WallMaterial::SS316L, E, false);

    EXPECT_TRUE(result.is_valid);
    EXPECT_NEAR(result.allowable_stress, 115.1e6, 1e3);

    const double S   = result.allowable_stress;
    const double R   = D / 2.0;
    const double expected = P * R / (S * E - 0.6 * P);
    EXPECT_NEAR(result.min_thickness, expected, 1e-9);
    // Sanity: should be a fraction of a millimetre for 10 bar on a 100mm vessel
    EXPECT_GT(result.min_thickness, 0.0);
    EXPECT_LT(result.min_thickness, 0.01); // < 10 mm
  }

  TEST(AsmeWallThicknessTest, InternalPressure_HigherPressureGivesThickerWall) {
    const auto r1 = CalculateAsmeMinWallThickness(0.1, 10e5, 293.15, WallMaterial::SS316L, 1.0, false);
    const auto r2 = CalculateAsmeMinWallThickness(0.1, 50e5, 293.15, WallMaterial::SS316L, 1.0, false);
    EXPECT_TRUE(r1.is_valid);
    EXPECT_TRUE(r2.is_valid);
    EXPECT_GT(r2.min_thickness, r1.min_thickness);
  }

  TEST(AsmeWallThicknessTest, InternalPressure_HigherTemperatureGivesThickerWall) {
    // Higher temperature → lower S → thicker wall
    const auto r_cold = CalculateAsmeMinWallThickness(0.1, 10e5, 293.15, WallMaterial::SS316L, 1.0, false);
    const auto r_hot  = CalculateAsmeMinWallThickness(0.1, 10e5, 773.15, WallMaterial::SS316L, 1.0, false);
    EXPECT_TRUE(r_cold.is_valid);
    EXPECT_TRUE(r_hot.is_valid);
    EXPECT_GT(r_hot.min_thickness, r_cold.min_thickness);
  }

  TEST(AsmeWallThicknessTest, InternalPressure_LowerWeldEfficiencyGivesThickerWall) {
    const auto r_full = CalculateAsmeMinWallThickness(0.1, 10e5, 293.15, WallMaterial::SS316L, 1.0, false);
    const auto r_spot = CalculateAsmeMinWallThickness(0.1, 10e5, 293.15, WallMaterial::SS316L, 0.85, false);
    EXPECT_TRUE(r_full.is_valid);
    EXPECT_TRUE(r_spot.is_valid);
    EXPECT_GT(r_spot.min_thickness, r_full.min_thickness);
  }

  TEST(AsmeWallThicknessTest, InvalidWhenPressureExceedsSE) {
    // If P >= S*E the denominator is zero or negative — should return is_valid = false
    // Use enormous pressure with Custom material (S = 0)
    const auto result = CalculateAsmeMinWallThickness(0.1, 1e9, 293.15, WallMaterial::Custom, 1.0, false);
    EXPECT_FALSE(result.is_valid);
  }

  // ---------------------------------------------------------------------------
  // CalculateAsmeMinWallThickness — vacuum service
  // ---------------------------------------------------------------------------
  // Vacuum uses P_eff = 101 325 Pa regardless of design_pressure
  // Hand calculation: 100 mm ID, SS316L at 20°C, E = 1.0
  //   P = 101325 Pa, R = 0.05 m, S = 115.1 MPa
  //   t = 101325 * 0.05 / (115.1e6 - 0.6 * 101325)
  //     = 5066.25 / (115100000 - 60795)
  //     = 5066.25 / 115039205
  //     ≈ 4.404e-5 m  ≈ 0.044 mm

  TEST(AsmeWallThicknessTest, VacuumService_UsesAtmosphericPressure) {
    const double D = 0.1;
    const auto result = CalculateAsmeMinWallThickness(D, 0.0, 293.15, WallMaterial::SS316L, 1.0, true);

    EXPECT_TRUE(result.is_valid);

    // Verify: result should match UG-27 formula at P = 101325 Pa
    const double P_vac = 101325.0;
    const double S     = result.allowable_stress;
    const double R     = D / 2.0;
    const double expected = P_vac * R / (S - 0.6 * P_vac);
    EXPECT_NEAR(result.min_thickness, expected, 1e-9);
  }

  TEST(AsmeWallThicknessTest, VacuumService_ThicknessIsPositive) {
    const auto result = CalculateAsmeMinWallThickness(0.3, 0.0, 293.15, WallMaterial::SS316L, 1.0, true);
    EXPECT_TRUE(result.is_valid);
    EXPECT_GT(result.min_thickness, 0.0);
  }

  // ---------------------------------------------------------------------------
  // BedWall convenience method
  // ---------------------------------------------------------------------------

  TEST(BedWallTest, CalculateMinThickness_DelegatesToFreeFunction) {
    BedWall wall;
    wall.ApplyMaterial(WallMaterial::SS316L);
    wall.design_pressure_      = 10e5;   // 10 bar gauge
    wall.design_temperature_   = 293.15; // 20°C
    wall.weld_joint_efficiency_ = 1.0;
    wall.vacuum_service_       = false;

    const double D = 0.1;
    const auto r1 = wall.CalculateMinThickness(D);
    const auto r2 = CalculateAsmeMinWallThickness(D, 10e5, 293.15, WallMaterial::SS316L, 1.0, false);

    EXPECT_EQ(r1.is_valid, r2.is_valid);
    EXPECT_DOUBLE_EQ(r1.min_thickness, r2.min_thickness);
    EXPECT_DOUBLE_EQ(r1.allowable_stress, r2.allowable_stress);
  }

  TEST(BedWallTest, VacuumServiceFlag_AffectsResult) {
    BedWall wall;
    wall.ApplyMaterial(WallMaterial::SS316L);
    wall.design_pressure_      = 5e5;
    wall.design_temperature_   = 293.15;
    wall.weld_joint_efficiency_ = 1.0;
    wall.vacuum_service_       = false;
    const auto r_pos = wall.CalculateMinThickness(0.1);

    wall.vacuum_service_ = true;
    const auto r_vac = wall.CalculateMinThickness(0.1);

    EXPECT_TRUE(r_pos.is_valid);
    EXPECT_TRUE(r_vac.is_valid);
    // Vacuum uses 101325 Pa; 5 bar is larger, so positive design gives thicker wall
    EXPECT_GT(r_pos.min_thickness, r_vac.min_thickness);
  }

} // end adsim namespace
