#include "adsim/math/unit_parser.h"

#include <gtest/gtest.h>
#include <cmath>
#include <string>

namespace adsim {

// ─── Helpers ────────────────────────────────────────────────────────────────

// Convenience: parse and return the converted value (ui * mul + add).
static double ConvertToBase(const std::string& expr, double ui_value) {
    std::string math;
    UnitXform xf;
    if (!ParseUnitSuffix(expr, math, xf)) return std::numeric_limits<double>::quiet_NaN();
    return xf.ToBase(ui_value);
}

// ─── ParseUnitSuffix — no-op cases ──────────────────────────────────────────

TEST(ParseUnitSuffix, NoSuffixReturnsFalse) {
    std::string math = "original";
    UnitXform xf{2.0, 5.0};  // should not be modified
    EXPECT_FALSE(ParseUnitSuffix("10 + 3", math, xf));
    EXPECT_EQ(math, "original");
    EXPECT_DOUBLE_EQ(xf.mul, 2.0);
    EXPECT_DOUBLE_EQ(xf.add, 5.0);
}

TEST(ParseUnitSuffix, UnknownSymbolReturnsFalse) {
    std::string math;
    UnitXform xf;
    EXPECT_FALSE(ParseUnitSuffix("10 [xyz]", math, xf));
    EXPECT_FALSE(ParseUnitSuffix("10 [megafoo]", math, xf));
}

TEST(ParseUnitSuffix, EmptyStringReturnsFalse) {
    std::string math;
    UnitXform xf;
    EXPECT_FALSE(ParseUnitSuffix("", math, xf));
}

TEST(ParseUnitSuffix, PlainNumberReturnsFalse) {
    std::string math;
    UnitXform xf;
    EXPECT_FALSE(ParseUnitSuffix("101325", math, xf));
    EXPECT_FALSE(ParseUnitSuffix("1e6", math, xf));
}

// ─── ParseUnitSuffix — math_out splitting ───────────────────────────────────

TEST(ParseUnitSuffix, SplitsMathFromUnit_Simple) {
    std::string math;
    UnitXform xf;
    EXPECT_TRUE(ParseUnitSuffix("10 [bar]", math, xf));
    // Leading whitespace before '[' is consumed by \s* in the regex;
    // the math part has no trailing space.
    EXPECT_EQ(math, "10");
}

TEST(ParseUnitSuffix, SplitsMathFromUnit_Expression) {
    std::string math;
    UnitXform xf;
    EXPECT_TRUE(ParseUnitSuffix("2 + 8 [kPa]", math, xf));
    EXPECT_EQ(math, "2 + 8");
}

TEST(ParseUnitSuffix, SplitsMathFromUnit_NoSpaceBeforeBracket) {
    std::string math;
    UnitXform xf;
    EXPECT_TRUE(ParseUnitSuffix("10[bar]", math, xf));
    EXPECT_EQ(math, "10");
}

TEST(ParseUnitSuffix, TrailingWhitespaceAfterBracket) {
    std::string math;
    UnitXform xf;
    EXPECT_TRUE(ParseUnitSuffix("10 [bar]  ", math, xf));
    // Whitespace before '[' and after ']' is absorbed by the regex.
    EXPECT_EQ(math, "10");
}

TEST(ParseUnitSuffix, UnitOnlyNoMathPart) {
    std::string math = "unchanged";
    UnitXform xf;
    EXPECT_TRUE(ParseUnitSuffix("[Pa]", math, xf));
    EXPECT_EQ(math, "");   // no math before the bracket
}

// ─── Pressure units ─────────────────────────────────────────────────────────

TEST(ParseUnitSuffix, Pressure_Pa) {
    std::string math; UnitXform xf;
    EXPECT_TRUE(ParseUnitSuffix("101325 [Pa]", math, xf));
    EXPECT_DOUBLE_EQ(xf.ToBase(101325.0), 101325.0);
}

TEST(ParseUnitSuffix, Pressure_kPa) {
    std::string math; UnitXform xf;
    EXPECT_TRUE(ParseUnitSuffix("101.325 [kPa]", math, xf));
    EXPECT_NEAR(xf.ToBase(101.325), 101325.0, 1e-6);
}

TEST(ParseUnitSuffix, Pressure_MPa) {
    std::string math; UnitXform xf;
    EXPECT_TRUE(ParseUnitSuffix("1 [MPa]", math, xf));
    EXPECT_DOUBLE_EQ(xf.ToBase(1.0), 1e6);
}

TEST(ParseUnitSuffix, Pressure_bar) {
    std::string math; UnitXform xf;
    EXPECT_TRUE(ParseUnitSuffix("10 [bar]", math, xf));
    EXPECT_DOUBLE_EQ(xf.ToBase(10.0), 1e6);
}

TEST(ParseUnitSuffix, Pressure_atm) {
    std::string math; UnitXform xf;
    EXPECT_TRUE(ParseUnitSuffix("1 [atm]", math, xf));
    EXPECT_NEAR(xf.ToBase(1.0), 101325.0, 1e-6);
}

TEST(ParseUnitSuffix, Pressure_psi) {
    std::string math; UnitXform xf;
    EXPECT_TRUE(ParseUnitSuffix("1 [psi]", math, xf));
    EXPECT_NEAR(xf.ToBase(1.0), 6894.757293168, 1e-6);
}

TEST(ParseUnitSuffix, Pressure_Torr) {
    std::string math; UnitXform xf;
    EXPECT_TRUE(ParseUnitSuffix("760 [Torr]", math, xf));
    EXPECT_NEAR(xf.ToBase(760.0), 101325.0, 1.0);  // 760 Torr ≈ 1 atm
}

// ─── Temperature units (affine) ─────────────────────────────────────────────

TEST(ParseUnitSuffix, Temperature_K) {
    std::string math; UnitXform xf;
    EXPECT_TRUE(ParseUnitSuffix("300 [K]", math, xf));
    EXPECT_DOUBLE_EQ(xf.ToBase(300.0), 300.0);
}

TEST(ParseUnitSuffix, Temperature_degC_FreezingPoint) {
    std::string math; UnitXform xf;
    EXPECT_TRUE(ParseUnitSuffix("0 [degC]", math, xf));
    EXPECT_NEAR(xf.ToBase(0.0), 273.15, 1e-10);
}

TEST(ParseUnitSuffix, Temperature_degC_BoilingPoint) {
    std::string math; UnitXform xf;
    EXPECT_TRUE(ParseUnitSuffix("100 [degC]", math, xf));
    EXPECT_NEAR(xf.ToBase(100.0), 373.15, 1e-10);
}

TEST(ParseUnitSuffix, Temperature_Unicode_degC) {
    std::string math; UnitXform xf;
    EXPECT_TRUE(ParseUnitSuffix("100 [°C]", math, xf));
    EXPECT_NEAR(xf.ToBase(100.0), 373.15, 1e-10);
}

TEST(ParseUnitSuffix, Temperature_degF_FreezingPoint) {
    std::string math; UnitXform xf;
    EXPECT_TRUE(ParseUnitSuffix("32 [degF]", math, xf));
    EXPECT_NEAR(xf.ToBase(32.0), 273.15, 1e-6);   // 32°F = 0°C = 273.15 K
}

TEST(ParseUnitSuffix, Temperature_degF_BoilingPoint) {
    std::string math; UnitXform xf;
    EXPECT_TRUE(ParseUnitSuffix("212 [degF]", math, xf));
    EXPECT_NEAR(xf.ToBase(212.0), 373.15, 1e-6);  // 212°F = 100°C = 373.15 K
}

TEST(ParseUnitSuffix, Temperature_AliasC) {
    std::string math; UnitXform xf;
    EXPECT_TRUE(ParseUnitSuffix("25 [C]", math, xf));
    EXPECT_NEAR(xf.ToBase(25.0), 298.15, 1e-10);
}

// ─── Length units ───────────────────────────────────────────────────────────

TEST(ParseUnitSuffix, Length_m) {
    std::string math; UnitXform xf;
    EXPECT_TRUE(ParseUnitSuffix("1 [m]", math, xf));
    EXPECT_DOUBLE_EQ(xf.ToBase(1.0), 1.0);
}

TEST(ParseUnitSuffix, Length_cm) {
    std::string math; UnitXform xf;
    EXPECT_TRUE(ParseUnitSuffix("100 [cm]", math, xf));
    EXPECT_DOUBLE_EQ(xf.ToBase(100.0), 1.0);
}

TEST(ParseUnitSuffix, Length_mm) {
    std::string math; UnitXform xf;
    EXPECT_TRUE(ParseUnitSuffix("25 [mm]", math, xf));
    EXPECT_NEAR(xf.ToBase(25.0), 0.025, 1e-15);
}

TEST(ParseUnitSuffix, Length_um_ASCII) {
    std::string math; UnitXform xf;
    EXPECT_TRUE(ParseUnitSuffix("500 [um]", math, xf));
    EXPECT_NEAR(xf.ToBase(500.0), 5e-4, 1e-15);
}

TEST(ParseUnitSuffix, Length_nm) {
    std::string math; UnitXform xf;
    EXPECT_TRUE(ParseUnitSuffix("1 [nm]", math, xf));
    EXPECT_DOUBLE_EQ(xf.ToBase(1.0), 1e-9);
}

TEST(ParseUnitSuffix, Length_Angstrom) {
    std::string math; UnitXform xf;
    EXPECT_TRUE(ParseUnitSuffix("1 [A]", math, xf));
    EXPECT_DOUBLE_EQ(xf.ToBase(1.0), 1e-10);
}

// ─── Time units ─────────────────────────────────────────────────────────────

TEST(ParseUnitSuffix, Time_s) {
    std::string math; UnitXform xf;
    EXPECT_TRUE(ParseUnitSuffix("60 [s]", math, xf));
    EXPECT_DOUBLE_EQ(xf.ToBase(60.0), 60.0);
}

TEST(ParseUnitSuffix, Time_min) {
    std::string math; UnitXform xf;
    EXPECT_TRUE(ParseUnitSuffix("1 [min]", math, xf));
    EXPECT_DOUBLE_EQ(xf.ToBase(1.0), 60.0);
}

TEST(ParseUnitSuffix, Time_h) {
    std::string math; UnitXform xf;
    EXPECT_TRUE(ParseUnitSuffix("1 [h]", math, xf));
    EXPECT_DOUBLE_EQ(xf.ToBase(1.0), 3600.0);
}

TEST(ParseUnitSuffix, Time_ms) {
    std::string math; UnitXform xf;
    EXPECT_TRUE(ParseUnitSuffix("500 [ms]", math, xf));
    EXPECT_NEAR(xf.ToBase(500.0), 0.5, 1e-12);
}

// ─── Energy / specific molar energy ─────────────────────────────────────────

TEST(ParseUnitSuffix, Energy_J) {
    std::string math; UnitXform xf;
    EXPECT_TRUE(ParseUnitSuffix("1 [J]", math, xf));
    EXPECT_DOUBLE_EQ(xf.ToBase(1.0), 1.0);
}

TEST(ParseUnitSuffix, Energy_kJ) {
    std::string math; UnitXform xf;
    EXPECT_TRUE(ParseUnitSuffix("1 [kJ]", math, xf));
    EXPECT_DOUBLE_EQ(xf.ToBase(1.0), 1000.0);
}

TEST(ParseUnitSuffix, SpecificMolarEnergy_kJ_mol) {
    std::string math; UnitXform xf;
    EXPECT_TRUE(ParseUnitSuffix("30 [kJ/mol]", math, xf));
    EXPECT_DOUBLE_EQ(xf.ToBase(30.0), 30000.0);
}

TEST(ParseUnitSuffix, SpecificMolarEnergy_kcal_mol) {
    std::string math; UnitXform xf;
    EXPECT_TRUE(ParseUnitSuffix("1 [kcal/mol]", math, xf));
    EXPECT_NEAR(xf.ToBase(1.0), 4184.0, 1e-9);
}

// ─── Mass and flow units ─────────────────────────────────────────────────────

TEST(ParseUnitSuffix, Mass_g) {
    std::string math; UnitXform xf;
    EXPECT_TRUE(ParseUnitSuffix("1000 [g]", math, xf));
    EXPECT_NEAR(xf.ToBase(1000.0), 1.0, 1e-12);
}

TEST(ParseUnitSuffix, MolarFlow_mol_s) {
    std::string math; UnitXform xf;
    EXPECT_TRUE(ParseUnitSuffix("1 [mol/s]", math, xf));
    EXPECT_DOUBLE_EQ(xf.ToBase(1.0), 1.0);
}

TEST(ParseUnitSuffix, MolarFlow_mol_min) {
    std::string math; UnitXform xf;
    EXPECT_TRUE(ParseUnitSuffix("60 [mol/min]", math, xf));
    EXPECT_DOUBLE_EQ(xf.ToBase(60.0), 1.0);
}

TEST(ParseUnitSuffix, Volume_L) {
    std::string math; UnitXform xf;
    EXPECT_TRUE(ParseUnitSuffix("1 [L]", math, xf));
    EXPECT_DOUBLE_EQ(xf.ToBase(1.0), 1e-3);
}

// ─── Inverse-time (mass transfer coefficient) ────────────────────────────────

TEST(ParseUnitSuffix, InverseTime_1_s) {
    std::string math; UnitXform xf;
    EXPECT_TRUE(ParseUnitSuffix("0.1 [1/s]", math, xf));
    EXPECT_DOUBLE_EQ(xf.ToBase(0.1), 0.1);
}

TEST(ParseUnitSuffix, InverseTime_1_min) {
    std::string math; UnitXform xf;
    EXPECT_TRUE(ParseUnitSuffix("6 [1/min]", math, xf));
    EXPECT_NEAR(xf.ToBase(6.0), 0.1, 1e-14);
}

// ─── Adsorbed concentration ──────────────────────────────────────────────────

TEST(ParseUnitSuffix, AdsorbedConcentration_mmol_g) {
    std::string math; UnitXform xf;
    EXPECT_TRUE(ParseUnitSuffix("5 [mmol/g]", math, xf));
    EXPECT_DOUBLE_EQ(xf.ToBase(5.0), 5.0);   // 1 mmol/g = 1 mol/kg
}

TEST(ParseUnitSuffix, AdsorbedConcentration_mol_kg) {
    std::string math; UnitXform xf;
    EXPECT_TRUE(ParseUnitSuffix("5 [mol/kg]", math, xf));
    EXPECT_DOUBLE_EQ(xf.ToBase(5.0), 5.0);
}

// ─── Mole fraction units ─────────────────────────────────────────────────────

TEST(ParseUnitSuffix, Percent) {
    std::string math; UnitXform xf;
    EXPECT_TRUE(ParseUnitSuffix("21 [%]", math, xf));
    EXPECT_NEAR(xf.ToBase(21.0), 0.21, 1e-12);   // 21% → 0.21 mol/mol
}

TEST(ParseUnitSuffix, PPM) {
    std::string math; UnitXform xf;
    EXPECT_TRUE(ParseUnitSuffix("400 [ppm]", math, xf));
    EXPECT_NEAR(xf.ToBase(400.0), 4e-4, 1e-15);
}

// ─── Round-trip: ToBase / FromBase ──────────────────────────────────────────

TEST(ParseUnitSuffix, RoundTrip_bar) {
    std::string math; UnitXform xf;
    ASSERT_TRUE(ParseUnitSuffix("10 [bar]", math, xf));
    double base = xf.ToBase(10.0);
    EXPECT_DOUBLE_EQ(xf.FromBase(base), 10.0);
}

TEST(ParseUnitSuffix, RoundTrip_degC) {
    std::string math; UnitXform xf;
    ASSERT_TRUE(ParseUnitSuffix("100 [degC]", math, xf));
    double base = xf.ToBase(100.0);
    EXPECT_NEAR(xf.FromBase(base), 100.0, 1e-10);
}

TEST(ParseUnitSuffix, RoundTrip_degF) {
    std::string math; UnitXform xf;
    ASSERT_TRUE(ParseUnitSuffix("212 [degF]", math, xf));
    double base = xf.ToBase(212.0);
    EXPECT_NEAR(xf.FromBase(base), 212.0, 1e-10);
}

// ─── Equivalence: different units for the same quantity ──────────────────────

TEST(ParseUnitSuffix, Equivalence_1atm_101325Pa) {
    std::string m1, m2; UnitXform xf1, xf2;
    ASSERT_TRUE(ParseUnitSuffix("1 [atm]",    m1, xf1));
    ASSERT_TRUE(ParseUnitSuffix("101325 [Pa]", m2, xf2));
    EXPECT_NEAR(xf1.ToBase(1.0), xf2.ToBase(101325.0), 1e-6);
}

TEST(ParseUnitSuffix, Equivalence_1bar_100kPa) {
    std::string m1, m2; UnitXform xf1, xf2;
    ASSERT_TRUE(ParseUnitSuffix("1 [bar]",    m1, xf1));
    ASSERT_TRUE(ParseUnitSuffix("100 [kPa]",  m2, xf2));
    EXPECT_NEAR(xf1.ToBase(1.0), xf2.ToBase(100.0), 1e-9);
}

TEST(ParseUnitSuffix, Equivalence_1h_3600s) {
    std::string m1, m2; UnitXform xf1, xf2;
    ASSERT_TRUE(ParseUnitSuffix("1 [h]",    m1, xf1));
    ASSERT_TRUE(ParseUnitSuffix("3600 [s]", m2, xf2));
    EXPECT_DOUBLE_EQ(xf1.ToBase(1.0), xf2.ToBase(3600.0));
}

TEST(ParseUnitSuffix, Equivalence_0degC_273_15K) {
    std::string m1, m2; UnitXform xf1, xf2;
    ASSERT_TRUE(ParseUnitSuffix("0 [degC]",      m1, xf1));
    ASSERT_TRUE(ParseUnitSuffix("273.15 [K]",    m2, xf2));
    EXPECT_NEAR(xf1.ToBase(0.0), xf2.ToBase(273.15), 1e-10);
}

TEST(ParseUnitSuffix, Equivalence_32degF_0degC) {
    std::string m1, m2; UnitXform xf1, xf2;
    ASSERT_TRUE(ParseUnitSuffix("32 [degF]", m1, xf1));
    ASSERT_TRUE(ParseUnitSuffix("0 [degC]",  m2, xf2));
    EXPECT_NEAR(xf1.ToBase(32.0), xf2.ToBase(0.0), 1e-6);
}

// ─── Normal/standard volumetric → molar flow ────────────────────────────────
// V_m at STP (0°C, 1 atm) = 22.414e-3 m³/mol, so 1 Nm³ ≈ 44.6149 mol.

TEST(ParseUnitSuffix, MolarFlow_Nm3_per_h) {
    std::string math; UnitXform xf;
    ASSERT_TRUE(ParseUnitSuffix("2000 [Nm^3/h]", math, xf));
    // 2000 / 22.414e-3 / 3600 ≈ 24.786 mol/s
    EXPECT_NEAR(xf.ToBase(2000.0), 2000.0 / 22.414e-3 / 3600.0, 1e-9);
}

TEST(ParseUnitSuffix, MolarFlow_Nm3_per_s) {
    std::string math; UnitXform xf;
    ASSERT_TRUE(ParseUnitSuffix("1 [Nm^3/s]", math, xf));
    EXPECT_NEAR(xf.ToBase(1.0), 1.0 / 22.414e-3, 1e-9);
}

TEST(ParseUnitSuffix, MolarFlow_SCCM) {
    std::string math; UnitXform xf;
    ASSERT_TRUE(ParseUnitSuffix("100 [SCCM]", math, xf));
    EXPECT_NEAR(xf.ToBase(100.0), 100.0 * 1e-6 / 22.414e-3 / 60.0, 1e-15);
}

TEST(ParseUnitSuffix, MolarFlow_SCFM) {
    std::string math; UnitXform xf;
    ASSERT_TRUE(ParseUnitSuffix("1 [SCFM]", math, xf));
    EXPECT_NEAR(xf.ToBase(1.0), 0.028316846592 / 22.414e-3 / 60.0, 1e-12);
}

// ─── Category output ────────────────────────────────────────────────────────

TEST(ParseUnitSuffix, CategoryOut_Pressure) {
    std::string math; UnitXform xf; UnitCategory cat = UnitCategory::Dimensionless;
    ASSERT_TRUE(ParseUnitSuffix("1 [bar]", math, xf, &cat));
    EXPECT_EQ(cat, UnitCategory::Pressure);
}

TEST(ParseUnitSuffix, CategoryOut_VolumetricFlow_m3_per_h) {
    std::string math; UnitXform xf; UnitCategory cat = UnitCategory::Dimensionless;
    ASSERT_TRUE(ParseUnitSuffix("2000 [m^3/h]", math, xf, &cat));
    EXPECT_EQ(cat, UnitCategory::VolumetricFlowRate);
}

TEST(ParseUnitSuffix, CategoryOut_MolarFlow_Nm3_per_h) {
    std::string math; UnitXform xf; UnitCategory cat = UnitCategory::Dimensionless;
    ASSERT_TRUE(ParseUnitSuffix("2000 [Nm^3/h]", math, xf, &cat));
    EXPECT_EQ(cat, UnitCategory::MolarFlowRate);
}

TEST(ParseUnitSuffix, CategoryOut_PowerIsCustom) {
    std::string math; UnitXform xf; UnitCategory cat = UnitCategory::Dimensionless;
    ASSERT_TRUE(ParseUnitSuffix("5 [kW]", math, xf, &cat));
    EXPECT_EQ(cat, UnitCategory::Custom);
}

TEST(ParseUnitSuffix, CategoryOut_NullptrIsAllowed) {
    std::string math; UnitXform xf;
    EXPECT_TRUE(ParseUnitSuffix("1 [Pa]", math, xf, nullptr));
}

// ─── ParseUnitExpression — compound units ───────────────────────────────────

static UnitInfo Unit(const std::string& text) {
    UnitInfo u;
    std::string why;
    EXPECT_TRUE(ParseUnitExpression(text, u, &why)) << text << ": " << why;
    return u;
}

TEST(ParseUnitExpression, SymbolsPrefixesAndNamedUnits) {
    EXPECT_DOUBLE_EQ(Unit("s").factor, 1.0);
    EXPECT_DOUBLE_EQ(Unit("min").factor, 60.0);
    EXPECT_DOUBLE_EQ(Unit("h").factor, 3600.0);
    EXPECT_DOUBLE_EQ(Unit("ms").factor, 1e-3);
    EXPECT_DOUBLE_EQ(Unit("mm").factor, 1e-3);
    EXPECT_DOUBLE_EQ(Unit("kPa").factor, 1e3);
    EXPECT_DOUBLE_EQ(Unit("mbar").factor, 100.0);
    EXPECT_DOUBLE_EQ(Unit("hPa").factor, 100.0);
    EXPECT_DOUBLE_EQ(Unit("kg").factor, 1.0);
    EXPECT_DOUBLE_EQ(Unit("mmol").factor, 1e-3);
    EXPECT_DOUBLE_EQ(Unit("\xC2\xB5m").factor, 1e-6);
    EXPECT_EQ(Unit("bar").dimension, Unit("Pa").dimension);
    EXPECT_EQ(DimensionName(Unit("atm").dimension), "pressure");
    EXPECT_EQ(DimensionName(Unit("min").dimension), "time");
}

TEST(ParseUnitExpression, CompoundUnitsInEveryNotationTheRegistryUses) {
    const Dimension htc = Unit("W/(m^2*K)").dimension;
    EXPECT_EQ(DimensionName(htc), "heat transfer coefficient");
    // Everything after '/' is in the denominator: W/m2-K is W/(m²·K).
    EXPECT_EQ(Unit("W/m2-K").dimension, htc);
    EXPECT_EQ(Unit("W/(m2.K)").dimension, htc);
    EXPECT_EQ(Unit("W m^-2 K^-1").dimension, htc);
    EXPECT_EQ(Unit("J/(kg-K)").dimension, Unit("J/(kg\xC2\xB7K)").dimension);
    EXPECT_EQ(Unit("kg/m3").dimension, Unit("kg m**-3").dimension);
    EXPECT_EQ(Unit("m\xC2\xB2/s").dimension, Unit("m2/s").dimension);
    EXPECT_DOUBLE_EQ(Unit("cm^2/s").factor, 1e-4);
    EXPECT_DOUBLE_EQ(Unit("kJ/mol").factor, 1e3);
    EXPECT_DOUBLE_EQ(Unit("mmol/g").factor, 1.0);
    EXPECT_DOUBLE_EQ(Unit("L/min").factor, 1e-3 / 60.0);
    EXPECT_DOUBLE_EQ(Unit("1/bar").factor, 1e-5);
    EXPECT_EQ(DimensionName(Unit("1/Pa").dimension), "inverse pressure");
    EXPECT_EQ(DimensionName(Unit("Nm3/h").dimension), "molar flow");
    EXPECT_NEAR(Unit("Nm3/h").factor, 1.0 / 22.414e-3 / 3600.0, 1e-15);
    EXPECT_TRUE(Unit("mol/mol").dimension.IsDimensionless());
    EXPECT_TRUE(Unit("-").dimension.IsDimensionless());
    EXPECT_DOUBLE_EQ(Unit("%").factor, 0.01);
}

TEST(ParseUnitExpression, OffsetTemperaturesAreAbsoluteOnlyAlone) {
    EXPECT_DOUBLE_EQ(Unit("degC").offset, 273.15);
    EXPECT_DOUBLE_EQ(Unit("\xC2\xB0" "C").offset, 273.15);
    // Inside a compound unit a degree is a temperature difference.
    EXPECT_DOUBLE_EQ(Unit("J/(kg-degC)").offset, 0.0);
    EXPECT_EQ(Unit("J/(kg-degC)").dimension, Unit("J/(kg-K)").dimension);
}

TEST(ParseUnitExpression, RejectsWhatItCannotRead) {
    UnitInfo u;
    for (const char* bad : {"xyz", "megafoo", "s-1", "m^", "kg/", "(m", "m)", "m /"}) {
        std::string why;
        EXPECT_FALSE(ParseUnitExpression(bad, u, &why)) << bad;
        EXPECT_FALSE(why.empty()) << bad;
    }
}

TEST(ConvertUnit, ConvertsWithinADimensionAndRefusesAcross) {
    EXPECT_DOUBLE_EQ(ConvertUnit(2.0, "min", "s"), 120.0);
    EXPECT_DOUBLE_EQ(ConvertUnit(25.0, "degC", "K"), 298.15);
    EXPECT_NEAR(ConvertUnit(298.15, "K", "degC"), 25.0, 1e-12);
    EXPECT_DOUBLE_EQ(ConvertUnit(1.0, "1/kPa", "1/bar"), 100.0);
    EXPECT_THROW(ConvertUnit(1.0, "bar", "s"), std::invalid_argument);
    EXPECT_THROW(ConvertUnit(1.0, "flurgle", "m"), std::invalid_argument);
}

// Unit symbols come from the LLNL units library: far more units than a
// hand-kept table, read strictly, with process-engineering overrides.
TEST(ParseUnitExpression, UsesTheLlnlUnitsLibrary) {
    EXPECT_NEAR(ConvertUnit(1.0, "furlong", "m"), 201.17, 0.01);
    EXPECT_NEAR(ConvertUnit(1.0, "inHg", "Pa"), 3386.39, 0.01);
    EXPECT_NEAR(ConvertUnit(1.0, "BTU/(lb*degF)", "J/(kg*K)"), 4186.8, 0.5);
    EXPECT_DOUBLE_EQ(ConvertUnit(1.0, "GJ/t", "kJ/kg"), 1e3);
    // C, F and A are degrees and the ångström here, not coulomb, farad, ampere.
    EXPECT_DOUBLE_EQ(Unit("C").offset, 273.15);
    EXPECT_EQ(DimensionName(Unit("F").dimension), "temperature");
    EXPECT_DOUBLE_EQ(Unit("A").factor, 1e-10);
    EXPECT_NEAR(ConvertUnit(60.0, "ccm", "m3/s"), 1e-6, 1e-18);
    // Units outside length, mass, time, temperature and amount are refused.
    UnitInfo u;
    for (const char* other : {"V", "ohm", "cd", "rad", "dB", "USD"})
        EXPECT_FALSE(ParseUnitExpression(other, u)) << other;
}

TEST(ParseUnitSuffix, CategoryFollowsTheDimension) {
    std::string math; UnitXform xf; UnitCategory cat = UnitCategory::Custom;
    ASSERT_TRUE(ParseUnitSuffix("10 [W/(m^2*K)]", math, xf, &cat));
    EXPECT_EQ(cat, UnitCategory::HeatTransferCoefficient);
    ASSERT_TRUE(ParseUnitSuffix("1.8e-5 [Pa*s]", math, xf, &cat));
    EXPECT_EQ(cat, UnitCategory::Viscosity);
    ASSERT_TRUE(ParseUnitSuffix("15 [%]", math, xf, &cat));
    EXPECT_EQ(cat, UnitCategory::MoleFraction);
    ASSERT_TRUE(ParseUnitSuffix("3 [1/bar]", math, xf, &cat));
    EXPECT_EQ(cat, UnitCategory::Custom);
}

TEST(SplitQuantity, ANumberFollowedByAUnit) {
    double n = 0.0;
    std::string u;
    EXPECT_TRUE(SplitQuantity("120 s", n, u));
    EXPECT_DOUBLE_EQ(n, 120.0);
    EXPECT_EQ(u, "s");
    EXPECT_TRUE(SplitQuantity(" 2 [min] ", n, u));
    EXPECT_EQ(u, "min");
    EXPECT_TRUE(SplitQuantity("-5 degC", n, u));
    EXPECT_DOUBLE_EQ(n, -5.0);
    EXPECT_TRUE(SplitQuantity("1.5e5 Pa", n, u));
    EXPECT_DOUBLE_EQ(n, 1.5e5);
    EXPECT_TRUE(SplitQuantity("2 W/(m^2*K)", n, u));
    for (const char* not_quantity : {"120", "2*x", "2 + 3", "2 x", "feed [min]", "0x10 s", ""})
        EXPECT_FALSE(SplitQuantity(not_quantity, n, u)) << not_quantity;
}

TEST(ParseUnitSuffix, FallsBackToUnitExpressions) {
    std::string math;
    UnitXform xf;
    UnitCategory cat = UnitCategory::Dimensionless;
    ASSERT_TRUE(ParseUnitSuffix("5 [kJ/(kg-K)]", math, xf, &cat));
    EXPECT_EQ(math, "5");
    EXPECT_DOUBLE_EQ(xf.mul, 1e3);
    EXPECT_EQ(cat, UnitCategory::HeatCapacity);
    ASSERT_TRUE(ParseUnitSuffix("2 [W/(m^2*K)]", math, xf, &cat));
    EXPECT_DOUBLE_EQ(xf.ToBase(2.0), 2.0);
}

} // namespace adsim
