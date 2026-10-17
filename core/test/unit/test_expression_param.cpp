// adsim includes
#include "adsim/math/expression_param.h"
#include "adsim/math/unit_parser.h"
#include "adsim/math/user_variables.h"
#include "adsim/objects/user_data.h"

// gtest includes
#include <gtest/gtest.h>

#include <cmath>
#include <limits>
#include <stdexcept>
#include <variant>

namespace adsim {

// ─── UserVariableRegistry ───────────────────────────────────────────────────

TEST(UserVariableRegistry, AddAndFind) {
    UserVariableRegistry reg;
    reg.Add({"A", 120.0, "", "pre-exp"});
    ASSERT_NE(reg.Find("A"), nullptr);
    EXPECT_DOUBLE_EQ(reg.Find("A")->value, 120.0);
}

TEST(UserVariableRegistry, AddOverwritesDuplicate) {
    UserVariableRegistry reg;
    reg.Add({"A", 1.0});
    reg.Add({"A", 2.0});
    EXPECT_EQ(reg.GetAll().size(), 1u);
    EXPECT_DOUBLE_EQ(reg.Find("A")->value, 2.0);
}

TEST(UserVariableRegistry, Remove) {
    UserVariableRegistry reg;
    reg.Add({"A", 1.0});
    reg.Remove("A");
    EXPECT_TRUE(reg.GetAll().empty());
    EXPECT_EQ(reg.Find("A"), nullptr);
}

TEST(UserVariableRegistry, BuildVarMap) {
    UserVariableRegistry reg;
    reg.Add({"A", 5.0});
    reg.Add({"B", 10.0});
    auto map = reg.BuildVarMap();
    EXPECT_DOUBLE_EQ(map.at("A"), 5.0);
    EXPECT_DOUBLE_EQ(map.at("B"), 10.0);
}

TEST(UserVariableRegistry, IsReservedName) {
    EXPECT_TRUE(UserVariableRegistry::IsReservedName("T"));
    EXPECT_TRUE(UserVariableRegistry::IsReservedName("P"));
    EXPECT_TRUE(UserVariableRegistry::IsReservedName("c"));
    EXPECT_TRUE(UserVariableRegistry::IsReservedName("R"));
    EXPECT_TRUE(UserVariableRegistry::IsReservedName("PI"));
    EXPECT_TRUE(UserVariableRegistry::IsReservedName("exp"));
    EXPECT_FALSE(UserVariableRegistry::IsReservedName("A"));
    EXPECT_FALSE(UserVariableRegistry::IsReservedName("myVar"));
}

// ─── ExpressionParam — plain double mode ────────────────────────────────────

TEST(ExpressionParam, DefaultIsZero) {
    ExpressionParam p;
    EXPECT_DOUBLE_EQ(p.Get(), 0.0);
    EXPECT_FALSE(p.HasExpression());
}

TEST(ExpressionParam, SetAndGetValue) {
    ExpressionParam p(3.14);
    EXPECT_DOUBLE_EQ(p.Get(), 3.14);
}

TEST(ExpressionParam, EvaluateNoExpressionIsNoOp) {
    ExpressionParam p(42.0);
    EvalContext ctx;
    EXPECT_TRUE(p.Evaluate(ctx));
    EXPECT_DOUBLE_EQ(p.Get(), 42.0);
}

// ─── ExpressionParam — expression mode ──────────────────────────────────────

TEST(ExpressionParam, SimpleArithmetic) {
    ExpressionParam p;
    p.SetExpr("2 + 3 * 4");
    EXPECT_TRUE(p.Evaluate(EvalContext{}));
    EXPECT_DOUBLE_EQ(p.Get(), 14.0);
}

TEST(ExpressionParam, BuiltInConstantR) {
    ExpressionParam p;
    p.SetExpr("R");
    EXPECT_TRUE(p.Evaluate(EvalContext{}));
    EXPECT_NEAR(p.Get(), 8.31446261815324, 1e-10);
}

TEST(ExpressionParam, BuiltInConstantPI) {
    ExpressionParam p;
    p.SetExpr("PI");
    EXPECT_TRUE(p.Evaluate(EvalContext{}));
    EXPECT_NEAR(p.Get(), 3.14159265358979, 1e-10);
}

TEST(ExpressionParam, FunctionExp) {
    ExpressionParam p;
    p.SetExpr("exp(1)");
    EXPECT_TRUE(p.Evaluate(EvalContext{}));
    EXPECT_NEAR(p.Get(), std::exp(1.0), 1e-12);
}

TEST(ExpressionParam, FunctionLog) {
    ExpressionParam p;
    p.SetExpr("log(exp(1))");
    EXPECT_TRUE(p.Evaluate(EvalContext{}));
    EXPECT_NEAR(p.Get(), 1.0, 1e-12);
}

TEST(ExpressionParam, FunctionSqrt) {
    ExpressionParam p;
    p.SetExpr("sqrt(9)");
    EXPECT_TRUE(p.Evaluate(EvalContext{}));
    EXPECT_NEAR(p.Get(), 3.0, 1e-12);
}

TEST(ExpressionParam, PowerOperator) {
    ExpressionParam p;
    p.SetExpr("2^10");
    EXPECT_TRUE(p.Evaluate(EvalContext{}));
    EXPECT_NEAR(p.Get(), 1024.0, 1e-9);
}

TEST(ExpressionParam, SimStateVariableT) {
    ExpressionParam p;
    p.SetExpr("T");
    EvalContext ctx;
    ctx.T = 350.0;
    EXPECT_TRUE(p.Evaluate(ctx));
    EXPECT_DOUBLE_EQ(p.Get(), 350.0);
}

TEST(ExpressionParam, SimStateVariableP) {
    ExpressionParam p;
    p.SetExpr("P * 2");
    EvalContext ctx;
    ctx.P = 5e5;
    EXPECT_TRUE(p.Evaluate(ctx));
    EXPECT_DOUBLE_EQ(p.Get(), 1e6);
}

TEST(ExpressionParam, ArrheniusExpression) {
    // A * exp(-B / (R * T))
    ExpressionParam p;
    p.SetExpr("A * exp(-B / (R * T))");

    UserVariableRegistry reg;
    reg.Add({"A", 120.0});
    reg.Add({"B", 1e4});

    EvalContext ctx;
    ctx.registry = &reg;
    ctx.T = 300.0;

    EXPECT_TRUE(p.Evaluate(ctx));
    double expected = 120.0 * std::exp(-1e4 / (8.31446261815324 * 300.0));
    EXPECT_NEAR(p.Get(), expected, std::abs(expected) * 1e-9);
}

TEST(ExpressionParam, UserDefinedVariable) {
    ExpressionParam p;
    p.SetExpr("A * 2");

    UserVariableRegistry reg;
    reg.Add({"A", 5.0});

    EvalContext ctx;
    ctx.registry = &reg;

    EXPECT_TRUE(p.Evaluate(ctx));
    EXPECT_DOUBLE_EQ(p.Get(), 10.0);
}

TEST(ExpressionParam, InvalidExpressionReturnsError) {
    ExpressionParam p;
    p.SetExpr("2 +* 3");
    EXPECT_FALSE(p.Evaluate(EvalContext{}));
    EXPECT_FALSE(p.LastError().empty());
}

TEST(ExpressionParam, UndefinedVariableReturnsError) {
    ExpressionParam p;
    p.SetExpr("UNDEFINED_VAR");
    EXPECT_FALSE(p.Evaluate(EvalContext{}));
    EXPECT_FALSE(p.LastError().empty());
}

TEST(ExpressionParam, ClearExprRestoresPlainMode) {
    ExpressionParam p(5.0);
    p.SetExpr("999");
    p.Evaluate(EvalContext{});
    EXPECT_DOUBLE_EQ(p.Get(), 999.0);

    p.ClearExpr();
    EXPECT_FALSE(p.HasExpression());
    // value_ retains the last evaluated result after clearing
    EXPECT_DOUBLE_EQ(p.Get(), 999.0);
}

TEST(ExpressionParam, NoRegistryNullSafe) {
    ExpressionParam p;
    p.SetExpr("2 + 2");
    EvalContext ctx;
    ctx.registry = nullptr;
    EXPECT_TRUE(p.Evaluate(ctx));
    EXPECT_DOUBLE_EQ(p.Get(), 4.0);
}

// ─── ExpressionParam — [unit] suffix: pressure ──────────────────────────────

TEST(ExpressionParam, UnitSuffix_1bar_Equals_1e5Pa) {
    ExpressionParam p;
    p.SetExpr("1 [bar]");
    EXPECT_TRUE(p.Evaluate(EvalContext{}));
    EXPECT_DOUBLE_EQ(p.Get(), 1e5);
}

TEST(ExpressionParam, UnitSuffix_10bar_Equals_1MPa) {
    ExpressionParam p;
    p.SetExpr("10 [bar]");
    EXPECT_TRUE(p.Evaluate(EvalContext{}));
    EXPECT_DOUBLE_EQ(p.Get(), 1e6);
}

TEST(ExpressionParam, UnitSuffix_1MPa) {
    ExpressionParam p;
    p.SetExpr("1 [MPa]");
    EXPECT_TRUE(p.Evaluate(EvalContext{}));
    EXPECT_DOUBLE_EQ(p.Get(), 1e6);
}

TEST(ExpressionParam, UnitSuffix_1atm) {
    ExpressionParam p;
    p.SetExpr("1 [atm]");
    EXPECT_TRUE(p.Evaluate(EvalContext{}));
    EXPECT_NEAR(p.Get(), 101325.0, 1e-6);
}

TEST(ExpressionParam, UnitSuffix_100kPa_Equals_1bar) {
    ExpressionParam p1, p2;
    p1.SetExpr("100 [kPa]");
    p2.SetExpr("1 [bar]");
    EXPECT_TRUE(p1.Evaluate(EvalContext{}));
    EXPECT_TRUE(p2.Evaluate(EvalContext{}));
    EXPECT_NEAR(p1.Get(), p2.Get(), 1e-9);
}

// ─── ExpressionParam — [unit] suffix: temperature (affine) ──────────────────

TEST(ExpressionParam, UnitSuffix_0degC_Equals_273_15K) {
    ExpressionParam p;
    p.SetExpr("0 [degC]");
    EXPECT_TRUE(p.Evaluate(EvalContext{}));
    EXPECT_NEAR(p.Get(), 273.15, 1e-10);
}

TEST(ExpressionParam, UnitSuffix_100degC_Equals_373_15K) {
    ExpressionParam p;
    p.SetExpr("100 [degC]");
    EXPECT_TRUE(p.Evaluate(EvalContext{}));
    EXPECT_NEAR(p.Get(), 373.15, 1e-10);
}

TEST(ExpressionParam, UnitSuffix_UnicodeDegreeCelsius) {
    ExpressionParam p;
    p.SetExpr("25 [°C]");
    EXPECT_TRUE(p.Evaluate(EvalContext{}));
    EXPECT_NEAR(p.Get(), 298.15, 1e-10);
}

TEST(ExpressionParam, UnitSuffix_32degF_Equals_273_15K) {
    ExpressionParam p;
    p.SetExpr("32 [degF]");
    EXPECT_TRUE(p.Evaluate(EvalContext{}));
    EXPECT_NEAR(p.Get(), 273.15, 1e-6);
}

TEST(ExpressionParam, UnitSuffix_212degF_Equals_373_15K) {
    ExpressionParam p;
    p.SetExpr("212 [degF]");
    EXPECT_TRUE(p.Evaluate(EvalContext{}));
    EXPECT_NEAR(p.Get(), 373.15, 1e-6);
}

// ─── ExpressionParam — [unit] suffix: length / time ─────────────────────────

TEST(ExpressionParam, UnitSuffix_2_5cm) {
    ExpressionParam p;
    p.SetExpr("2.5 [cm]");
    EXPECT_TRUE(p.Evaluate(EvalContext{}));
    EXPECT_NEAR(p.Get(), 0.025, 1e-15);
}

TEST(ExpressionParam, UnitSuffix_500um) {
    ExpressionParam p;
    p.SetExpr("500 [um]");
    EXPECT_TRUE(p.Evaluate(EvalContext{}));
    EXPECT_NEAR(p.Get(), 5e-4, 1e-15);
}

TEST(ExpressionParam, UnitSuffix_1min_Equals_60s) {
    ExpressionParam p;
    p.SetExpr("1 [min]");
    EXPECT_TRUE(p.Evaluate(EvalContext{}));
    EXPECT_DOUBLE_EQ(p.Get(), 60.0);
}

TEST(ExpressionParam, UnitSuffix_30kJ_mol) {
    ExpressionParam p;
    p.SetExpr("30 [kJ/mol]");
    EXPECT_TRUE(p.Evaluate(EvalContext{}));
    EXPECT_DOUBLE_EQ(p.Get(), 30000.0);
}

// ─── ExpressionParam — [unit] with math expressions ─────────────────────────

TEST(ExpressionParam, UnitSuffix_MathExpression_SumThenConvert) {
    // (2 + 8) [bar] = 10 bar = 1e6 Pa
    ExpressionParam p;
    p.SetExpr("2 + 8 [bar]");
    EXPECT_TRUE(p.Evaluate(EvalContext{}));
    EXPECT_DOUBLE_EQ(p.Get(), 1e6);
}

TEST(ExpressionParam, UnitSuffix_ScientificNotation) {
    ExpressionParam p;
    p.SetExpr("1e6 [Pa]");
    EXPECT_TRUE(p.Evaluate(EvalContext{}));
    EXPECT_DOUBLE_EQ(p.Get(), 1e6);
}

TEST(ExpressionParam, UnitSuffix_MathWithPower) {
    // sqrt(4) [bar] = 2 bar = 2e5 Pa
    ExpressionParam p;
    p.SetExpr("sqrt(4) [bar]");
    EXPECT_TRUE(p.Evaluate(EvalContext{}));
    EXPECT_NEAR(p.Get(), 2e5, 1e-9);
}

TEST(ExpressionParam, UnitSuffix_MathWithConstants) {
    // PI [bar] = π × 1e5 Pa
    ExpressionParam p;
    p.SetExpr("PI [bar]");
    EXPECT_TRUE(p.Evaluate(EvalContext{}));
    EXPECT_NEAR(p.Get(), 3.14159265358979 * 1e5, 1.0);
}

TEST(ExpressionParam, UnitSuffix_NoBracket_NoSpaceChange) {
    // No unit suffix — evaluated as before
    ExpressionParam p;
    p.SetExpr("1e5");
    EXPECT_TRUE(p.Evaluate(EvalContext{}));
    EXPECT_DOUBLE_EQ(p.Get(), 1e5);
}

// ─── ExpressionParam — [unit] with simulation state variables ────────────────

TEST(ExpressionParam, UnitSuffix_SimStateScaled_bar) {
    // P [bar] where P = 10 (Pa by convention in sim) → 10 * 1e5 = 1e6
    // (This tests that state variables work alongside unit conversion)
    ExpressionParam p;
    p.SetExpr("P [bar]");
    EvalContext ctx;
    ctx.P = 10.0;
    EXPECT_TRUE(p.Evaluate(ctx));
    EXPECT_DOUBLE_EQ(p.Get(), 1e6);
}

TEST(ExpressionParam, UnitSuffix_TemperatureExpression_WithK) {
    // T [K] = T * 1 (identity — result is T in kelvin)
    ExpressionParam p;
    p.SetExpr("T [K]");
    EvalContext ctx;
    ctx.T = 350.0;
    EXPECT_TRUE(p.Evaluate(ctx));
    EXPECT_DOUBLE_EQ(p.Get(), 350.0);
}

// ─── ExpressionParam — [unit] with user-defined variables ───────────────────

TEST(ExpressionParam, UnitSuffix_UserVariable_Converted) {
    // P_hi is stored as 10 (in bar units by user intent);
    // P_hi [bar] converts it to Pa → 1e6
    ExpressionParam p;
    p.SetExpr("P_hi [bar]");

    UserVariableRegistry reg;
    reg.Add({"P_hi", 10.0});

    EvalContext ctx;
    ctx.registry = &reg;
    EXPECT_TRUE(p.Evaluate(ctx));
    EXPECT_DOUBLE_EQ(p.Get(), 1e6);
}

TEST(ExpressionParam, UnitSuffix_UserVariable_ArithmeticThenConvert) {
    // (P_lo + delta) [bar]
    ExpressionParam p;
    p.SetExpr("(P_lo + delta) [bar]");

    UserVariableRegistry reg;
    reg.Add({"P_lo",  3.0});
    reg.Add({"delta", 7.0});

    EvalContext ctx;
    ctx.registry = &reg;
    EXPECT_TRUE(p.Evaluate(ctx));
    EXPECT_DOUBLE_EQ(p.Get(), 1e6);  // (3 + 7) bar = 1e6 Pa
}

// ─── User-variable inline units vs declared category ────────────────────────
// A user variable left at the DEFAULT (Dimensionless) category — the "unspecified"
// state — must accept an inline `[unit]` and convert its value rather than
// aborting the whole EvaluateAllExpressions cascade. Regression for
// "1 [bar] breaks global vars".
//
// The adopted dimension is used only for THIS evaluation: it is deliberately NOT
// persisted onto the variable. Persisting it would overwrite the user's own
// category (the web exposes a picker) and make a later edit to a different
// dimension look like a mismatch against the auto-adopted one — see
// DimensionlessVariableCanBeReDimensioned below.
TEST(EvaluateAllExpressions, DimensionlessVariableAdoptsInlineUnit) {
    UserData ud;
    UserVariable ph;
    ph.name = "P_h";
    ph.value_expr = "1 [bar]";
    ph.unit_category = UnitCategory::Dimensionless;  // web default
    ud.variable_registry.Add(std::move(ph));

    std::string err;
    ASSERT_TRUE(ud.EvaluateAllExpressions(err)) << err;

    const UserVariable* v = ud.variable_registry.Find("P_h");
    ASSERT_NE(v, nullptr);
    EXPECT_DOUBLE_EQ(v->value, 1e5);  // 1 bar -> 1e5 Pa (SI base)
    // The category is NOT mutated: it stays at the unspecified default.
    EXPECT_EQ(v->unit_category, UnitCategory::Dimensionless);
}

// A Dimensionless variable that once evaluated with a `[bar]` unit must still
// accept a later expression in a DIFFERENT dimension. Previously the first
// evaluation persisted Pressure onto the variable, so the second (`[K]`) hit the
// explicit-mismatch path and aborted the whole cascade.
TEST(EvaluateAllExpressions, DimensionlessVariableCanBeReDimensioned) {
    UserData ud;
    UserVariable v;
    v.name = "X";
    v.value_expr = "1 [bar]";
    v.unit_category = UnitCategory::Dimensionless;
    ud.variable_registry.Add(std::move(v));

    std::string err;
    ASSERT_TRUE(ud.EvaluateAllExpressions(err)) << err;
    EXPECT_DOUBLE_EQ(ud.variable_registry.Find("X")->value, 1e5);

    // Repurpose the same variable to a temperature — must not be rejected.
    ud.variable_registry.Find("X")->value_expr = "5 [K]";
    ASSERT_TRUE(ud.EvaluateAllExpressions(err)) << err;
    EXPECT_DOUBLE_EQ(ud.variable_registry.Find("X")->value, 5.0);
    EXPECT_EQ(ud.variable_registry.Find("X")->unit_category, UnitCategory::Dimensionless);
}

TEST(EvaluateAllExpressions, ExplicitCategoryStillRejectsMismatchedUnit) {
    UserData ud;
    UserVariable t;
    t.name = "T_set";
    t.value_expr = "1 [bar]";
    t.unit_category = UnitCategory::Temperature;  // explicit, contradicts [bar]
    ud.variable_registry.Add(std::move(t));

    std::string err;
    EXPECT_FALSE(ud.EvaluateAllExpressions(err));
    EXPECT_NE(err.find("T_set"), std::string::npos) << err;
}

TEST(EvaluateAllExpressions, ExplicitMatchingCategoryConverts) {
    UserData ud;
    UserVariable p;
    p.name = "P_h";
    p.value_expr = "1 [bar]";
    p.unit_category = UnitCategory::Pressure;  // explicit, matches [bar]
    ud.variable_registry.Add(std::move(p));

    std::string err;
    ASSERT_TRUE(ud.EvaluateAllExpressions(err)) << err;
    EXPECT_DOUBLE_EQ(ud.variable_registry.Find("P_h")->value, 1e5);
    EXPECT_EQ(ud.variable_registry.Find("P_h")->unit_category, UnitCategory::Pressure);
}

TEST(ExpressionParam, UnitSuffix_Arrhenius_WithEaInkJmol) {
    // The [unit] suffix only applies to the WHOLE trailing expression.
    // The idiomatic way to use a kJ/mol activation energy:
    //   1) Evaluate "30 [kJ/mol]" separately → 30000 J/mol
    //   2) Store that in a user variable Ea
    //   3) Reference Ea in the Arrhenius expression
    ExpressionParam ea_expr;
    ea_expr.SetExpr("30 [kJ/mol]");
    EXPECT_TRUE(ea_expr.Evaluate(EvalContext{}));
    EXPECT_DOUBLE_EQ(ea_expr.Get(), 30000.0);

    UserVariableRegistry reg;
    reg.Add({"A",  1e5});
    reg.Add({"Ea", ea_expr.Get()});   // Ea = 30000 J/mol

    ExpressionParam p_ref, p_via_var;
    p_ref.SetExpr("A * exp(-30000 / (R * T))");
    p_via_var.SetExpr("A * exp(-Ea / (R * T))");

    EvalContext ctx;
    ctx.registry = &reg;
    ctx.T = 500.0;

    EXPECT_TRUE(p_ref.Evaluate(ctx));
    EXPECT_TRUE(p_via_var.Evaluate(ctx));
    EXPECT_NEAR(p_ref.Get(), p_via_var.Get(), std::abs(p_ref.Get()) * 1e-9);
}

// ─── ExpressionParam — [unit] error / edge cases ────────────────────────────

TEST(ExpressionParam, UnitSuffix_UnknownSymbol_ProducesError) {
    // [xyz] is not in the dictionary → ParseUnitSuffix returns false →
    // muParser receives "10 [xyz]" verbatim → parse error
    ExpressionParam p;
    p.SetExpr("10 [xyz]");
    EXPECT_FALSE(p.Evaluate(EvalContext{}));
    EXPECT_FALSE(p.LastError().empty());
}

TEST(ExpressionParam, UnitSuffix_WhitespaceVariants) {
    // Extra spaces around the bracket should not affect the result
    ExpressionParam p1, p2, p3;
    p1.SetExpr("10[bar]");
    p2.SetExpr("10 [bar]");
    p3.SetExpr("10  [bar]  ");
    EXPECT_TRUE(p1.Evaluate(EvalContext{}));
    EXPECT_TRUE(p2.Evaluate(EvalContext{}));
    EXPECT_TRUE(p3.Evaluate(EvalContext{}));
    EXPECT_DOUBLE_EQ(p1.Get(), p2.Get());
    EXPECT_DOUBLE_EQ(p2.Get(), p3.Get());
}

TEST(ExpressionParam, UnitSuffix_NegativeValue) {
    // Negative values should convert correctly (no offset issue for multiplicative units)
    ExpressionParam p;
    p.SetExpr("-5 [bar]");
    EXPECT_TRUE(p.Evaluate(EvalContext{}));
    EXPECT_DOUBLE_EQ(p.Get(), -5e5);
}

TEST(ExpressionParam, UnitSuffix_ZeroValueMultiplicative) {
    ExpressionParam p;
    p.SetExpr("0 [bar]");
    EXPECT_TRUE(p.Evaluate(EvalContext{}));
    EXPECT_DOUBLE_EQ(p.Get(), 0.0);
}

TEST(ExpressionParam, UnitSuffix_ZeroValueAffine_0degC) {
    // 0°C → 273.15 K (affine: zero in Celsius is NOT zero in Kelvin)
    ExpressionParam p;
    p.SetExpr("0 [degC]");
    EXPECT_TRUE(p.Evaluate(EvalContext{}));
    EXPECT_NEAR(p.Get(), 273.15, 1e-10);
}

TEST(ExpressionParam, UnitSuffix_IdentityUnit_Pa) {
    // [Pa] has mul=1, add=0 — result equals the numeric value
    ExpressionParam p;
    p.SetExpr("101325 [Pa]");
    EXPECT_TRUE(p.Evaluate(EvalContext{}));
    EXPECT_DOUBLE_EQ(p.Get(), 101325.0);
}

TEST(ExpressionParam, UnitSuffix_FloatingPointLiteral) {
    ExpressionParam p;
    p.SetExpr("2.5 [bar]");
    EXPECT_TRUE(p.Evaluate(EvalContext{}));
    EXPECT_NEAR(p.Get(), 2.5e5, 1e-9);
}

TEST(ExpressionParam, UnitSuffix_Percent) {
    ExpressionParam p;
    p.SetExpr("21 [%]");
    EXPECT_TRUE(p.Evaluate(EvalContext{}));
    EXPECT_NEAR(p.Get(), 0.21, 1e-12);
}

// ─── Re-evaluation with updated registry (regression for stale-double bug) ────
//
// The bug: an expression like "P_high" is evaluated once when the user commits
// a BC field; the result is stored in a plain double.  Later the user changes
// P_high in the Parameters tab but does NOT recommit the BC field.  The plain
// double still holds the old value.
//
// The fix: (1) per-frame sync when the widget is visible, (2) SyncAllExpressions()
// called before the solver starts.  These tests verify the core invariant that
// Evaluate() always uses the *current* registry contents, never a stale snapshot.

// Changing a plain-value variable in the registry is visible on the very next Evaluate().
TEST(ExpressionParam, ReEval_RegistryValueChange_IsVisible) {
    UserVariableRegistry reg;
    reg.Add({"P_high", 1e6});  // 1 MPa

    ExpressionParam p;
    p.SetExpr("P_high");

    EvalContext ctx;
    ctx.registry = &reg;

    EXPECT_TRUE(p.Evaluate(ctx));
    EXPECT_DOUBLE_EQ(p.Get(), 1e6);

    // User changes P_high in the Parameters tab
    reg.Find("P_high")->value = 2e6;

    EXPECT_TRUE(p.Evaluate(ctx));
    EXPECT_DOUBLE_EQ(p.Get(), 2e6);  // must reflect the new value
}

// The plain-double snapshot pattern: storing Evaluate() result in a double
// is correct only when followed by a re-evaluation before each use.
TEST(ExpressionParam, ReEval_PlainDoubleSnapshotGoesStalWithoutReEval) {
    UserVariableRegistry reg;
    reg.Add({"P_high", 1e6});

    ExpressionParam expr;
    expr.SetExpr("P_high");

    EvalContext ctx;
    ctx.registry = &reg;
    expr.Evaluate(ctx);

    // Simulate what ConstantBoundaryValue::constant_value_ holds after first commit
    double snapshot = expr.Get();
    EXPECT_DOUBLE_EQ(snapshot, 1e6);

    // User changes the variable
    reg.Find("P_high")->value = 3e6;

    // Without re-evaluation the snapshot is stale
    EXPECT_DOUBLE_EQ(snapshot, 1e6);  // stale!

    // After re-evaluation and writeback the snapshot is fresh
    expr.Evaluate(ctx);
    snapshot = expr.Get();
    EXPECT_DOUBLE_EQ(snapshot, 3e6);  // fresh
}

// Expressions that combine multiple variables reflect all of their updates.
TEST(ExpressionParam, ReEval_MultiVariableExpression_AllUpdatesVisible) {
    UserVariableRegistry reg;
    reg.Add({"P_lo", 1e5});
    reg.Add({"P_hi", 1e6});

    ExpressionParam p;
    p.SetExpr("P_hi - P_lo");

    EvalContext ctx;
    ctx.registry = &reg;

    EXPECT_TRUE(p.Evaluate(ctx));
    EXPECT_DOUBLE_EQ(p.Get(), 9e5);

    reg.Find("P_lo")->value = 2e5;
    reg.Find("P_hi")->value = 5e5;

    EXPECT_TRUE(p.Evaluate(ctx));
    EXPECT_DOUBLE_EQ(p.Get(), 3e5);
}

// Consecutive value changes: each re-evaluation tracks the latest registry state.
TEST(ExpressionParam, ReEval_ConsecutiveChangesAllTracked) {
    UserVariableRegistry reg;
    reg.Add({"k", 100.0});

    ExpressionParam p;
    p.SetExpr("k * 2");

    EvalContext ctx;
    ctx.registry = &reg;

    for (double v : {100.0, 200.0, 50.0, 0.0, -10.0}) {
        reg.Find("k")->value = v;
        EXPECT_TRUE(p.Evaluate(ctx));
        EXPECT_DOUBLE_EQ(p.Get(), v * 2.0);
    }
}

// Unit-suffix expressions also pick up registry changes on re-evaluation.
TEST(ExpressionParam, ReEval_UnitSuffix_RegistryChange) {
    UserVariableRegistry reg;
    reg.Add({"P_bar", 10.0});  // 10 bar = 1 MPa

    ExpressionParam p;
    p.SetExpr("P_bar [bar]");

    EvalContext ctx;
    ctx.registry = &reg;

    EXPECT_TRUE(p.Evaluate(ctx));
    EXPECT_DOUBLE_EQ(p.Get(), 1e6);

    reg.Find("P_bar")->value = 20.0;  // 20 bar = 2 MPa

    EXPECT_TRUE(p.Evaluate(ctx));
    EXPECT_DOUBLE_EQ(p.Get(), 2e6);
}

// Temperature (affine) unit suffix picks up registry changes correctly.
TEST(ExpressionParam, ReEval_TemperatureAffine_RegistryChange) {
    UserVariableRegistry reg;
    reg.Add({"T_amb", 25.0});  // 25°C = 298.15 K

    ExpressionParam p;
    p.SetExpr("T_amb [degC]");

    EvalContext ctx;
    ctx.registry = &reg;

    EXPECT_TRUE(p.Evaluate(ctx));
    EXPECT_NEAR(p.Get(), 298.15, 1e-10);

    reg.Find("T_amb")->value = 100.0;  // 100°C = 373.15 K

    EXPECT_TRUE(p.Evaluate(ctx));
    EXPECT_NEAR(p.Get(), 373.15, 1e-10);
}

// ─── Chain variable updates ───────────────────────────────────────────────────
//
// Variables B="A*2" and C="B*3" form a dependency chain.  When A changes,
// re-evaluating B then C in order gives the correct final value for C.

TEST(ExpressionParam, ChainDeps_BDependsOnA) {
    UserVariableRegistry reg;
    reg.Add({"A", 5.0});
    reg.Add({"B", 0.0});  // B = A * 2, maintained manually in this test

    ExpressionParam exprB;
    exprB.SetExpr("A * 2");

    EvalContext ctx;
    ctx.registry = &reg;

    exprB.Evaluate(ctx);
    reg.Find("B")->value = exprB.Get();  // propagate B
    EXPECT_DOUBLE_EQ(reg.Find("B")->value, 10.0);

    reg.Find("A")->value = 7.0;
    exprB.Evaluate(ctx);
    reg.Find("B")->value = exprB.Get();  // re-propagate B
    EXPECT_DOUBLE_EQ(reg.Find("B")->value, 14.0);
}

TEST(ExpressionParam, ChainDeps_CDependsOnB_DependsOnA) {
    UserVariableRegistry reg;
    reg.Add({"A", 3.0});
    reg.Add({"B", 0.0});

    ExpressionParam exprB, exprC;
    exprB.SetExpr("A * 2");
    exprC.SetExpr("B * 3");

    EvalContext ctx;
    ctx.registry = &reg;

    // First pass: evaluate B, write back, then evaluate C
    exprB.Evaluate(ctx);
    reg.Find("B")->value = exprB.Get();
    exprC.Evaluate(ctx);
    EXPECT_DOUBLE_EQ(exprC.Get(), 3.0 * 2.0 * 3.0);  // A=3 → B=6 → C=18

    // Change A and re-propagate
    reg.Find("A")->value = 5.0;
    exprB.Evaluate(ctx);
    reg.Find("B")->value = exprB.Get();
    exprC.Evaluate(ctx);
    EXPECT_DOUBLE_EQ(exprC.Get(), 5.0 * 2.0 * 3.0);  // A=5 → B=10 → C=30
}

// ─── ExpressionParam fields never go stale (contrast with plain-double) ───────
//
// When an ExpressionParam stores its own expression (like MTC params do), the
// solver calls Evaluate() at runtime and always sees the current registry.

TEST(ExpressionParam, ExprFieldAlwaysFresh_ContrastToDoubleSnapshot) {
    UserVariableRegistry reg;
    reg.Add({"k0", 100.0});

    // ExpressionParam stores the expression — mirrors what MTC / isotherm params do
    ExpressionParam mtc_param;
    mtc_param.SetExpr("k0");

    // Separate plain-double that mimics ConstantBoundaryValue::constant_value_
    double snapshot = 0.0;

    EvalContext ctx;
    ctx.registry = &reg;

    // Initial commit: both see k0 = 100
    mtc_param.Evaluate(ctx);
    snapshot = mtc_param.Get();
    EXPECT_DOUBLE_EQ(mtc_param.Get(), 100.0);
    EXPECT_DOUBLE_EQ(snapshot, 100.0);

    // User changes k0
    reg.Find("k0")->value = 250.0;

    // ExpressionParam is fresh on every Evaluate() call
    mtc_param.Evaluate(ctx);
    EXPECT_DOUBLE_EQ(mtc_param.Get(), 250.0);  // fresh ✓

    // Plain-double snapshot is stale until explicitly refreshed
    EXPECT_DOUBLE_EQ(snapshot, 100.0);         // stale ✗

    // Refreshing the snapshot is the fix
    snapshot = mtc_param.Get();
    EXPECT_DOUBLE_EQ(snapshot, 250.0);         // now fresh ✓
}

// ─── Registry mutation methods all propagate to Evaluate() ───────────────────

TEST(ExpressionParam, RegistryMutation_DirectFieldWrite) {
    UserVariableRegistry reg;
    reg.Add({"x", 1.0});
    ExpressionParam p;
    p.SetExpr("x");
    EvalContext ctx;
    ctx.registry = &reg;

    reg.Find("x")->value = 42.0;
    p.Evaluate(ctx);
    EXPECT_DOUBLE_EQ(p.Get(), 42.0);
}

TEST(ExpressionParam, RegistryMutation_AddOverwrite) {
    UserVariableRegistry reg;
    reg.Add({"x", 1.0});
    ExpressionParam p;
    p.SetExpr("x");
    EvalContext ctx;
    ctx.registry = &reg;

    reg.Add({"x", 99.0});  // overwrite via Add
    p.Evaluate(ctx);
    EXPECT_DOUBLE_EQ(p.Get(), 99.0);
}

TEST(ExpressionParam, RegistryMutation_BuildVarMapReflectsLatestValue) {
    UserVariableRegistry reg;
    reg.Add({"A", 1.0});
    reg.Add({"B", 2.0});

    reg.Find("A")->value = 10.0;

    auto m = reg.BuildVarMap();
    EXPECT_DOUBLE_EQ(m.at("A"), 10.0);
    EXPECT_DOUBLE_EQ(m.at("B"), 2.0);
}

// ─── ExpressionParam — inline unit category exposure ────────────────────────

TEST(ExpressionParam, InlineUnitCategory_PressureCaptured) {
    ExpressionParam p;
    p.SetExpr("1 [bar]");
    ASSERT_TRUE(p.Evaluate(EvalContext{}));
    EXPECT_TRUE(p.HasInlineUnit());
    EXPECT_EQ(p.InlineUnitCategory(), UnitCategory::Pressure);
}

TEST(ExpressionParam, InlineUnitCategory_PlainNumberHasNoUnit) {
    ExpressionParam p;
    p.SetExpr("42");
    ASSERT_TRUE(p.Evaluate(EvalContext{}));
    EXPECT_FALSE(p.HasInlineUnit());
}

TEST(ExpressionParam, InlineUnitCategory_MathExprWithoutUnitHasNoUnit) {
    ExpressionParam p;
    p.SetExpr("2 + 3");
    ASSERT_TRUE(p.Evaluate(EvalContext{}));
    EXPECT_FALSE(p.HasInlineUnit());
}

TEST(ExpressionParam, InlineUnitCategory_Nm3PerHIsMolarFlow) {
    // 2000 Nm³/h ≈ 24.78 mol/s (V_m at STP = 22.414e-3 m³/mol)
    ExpressionParam p;
    p.SetExpr("2000 [Nm^3/h]");
    ASSERT_TRUE(p.Evaluate(EvalContext{}));
    EXPECT_TRUE(p.HasInlineUnit());
    EXPECT_EQ(p.InlineUnitCategory(), UnitCategory::MolarFlowRate);
    EXPECT_NEAR(p.Get(), 2000.0 / 22.414e-3 / 3600.0, 1e-9);
}

TEST(ExpressionParam, InlineUnitCategory_M3PerHIsVolumetricFlow) {
    // 2000 m³/h converts to 0.5556 m³/s — but the category is VolumetricFlowRate,
    // NOT MolarFlowRate. The user-data validator uses this to reject the mismatch.
    ExpressionParam p;
    p.SetExpr("2000 [m^3/h]");
    ASSERT_TRUE(p.Evaluate(EvalContext{}));
    EXPECT_TRUE(p.HasInlineUnit());
    EXPECT_EQ(p.InlineUnitCategory(), UnitCategory::VolumetricFlowRate);
}

TEST(ExpressionParam, InlineUnitCategory_StaleStateClearedOnSetExpr) {
    ExpressionParam p;
    p.SetExpr("1 [bar]");
    ASSERT_TRUE(p.Evaluate(EvalContext{}));
    ASSERT_TRUE(p.HasInlineUnit());

    p.SetExpr("42");
    EXPECT_FALSE(p.HasInlineUnit());
    ASSERT_TRUE(p.Evaluate(EvalContext{}));
    EXPECT_FALSE(p.HasInlineUnit());
}

// ─── F-05: non-finite expression values ─────────────────────────────────────

TEST(ExpressionParam, UsesStateVariables) {
    auto uses = [](const char* expr) {
        ExpressionParam p;
        p.SetExpr(expr);
        return p.UsesStateVariables();
    };
    EXPECT_TRUE(uses("2*T"));
    EXPECT_TRUE(uses("k/v"));          // k undefined: still parsed
    EXPECT_TRUE(uses("1/y [bar]"));    // unit suffix stripped first
    EXPECT_TRUE(uses("P + c"));
    EXPECT_FALSE(uses("1/0"));
    EXPECT_FALSE(uses("Tin*2"));       // a name containing T is not T
    EXPECT_FALSE(uses("R*PI"));
    EXPECT_FALSE(uses("1.5e5"));       // plain literal
    EXPECT_FALSE(uses("1/("));         // unparseable
}

TEST(ExpressionParam, IsNonFiniteConstant) {
    ExpressionParam p;
    p.SetExpr("1/0");
    ASSERT_TRUE(p.Evaluate(EvalContext{}));
    EXPECT_TRUE(std::isinf(p.Get()));
    EXPECT_TRUE(p.IsNonFiniteConstant());

    p.SetExpr("inf");  // a plain literal strtod reads as infinity
    ASSERT_TRUE(p.Evaluate(EvalContext{}));
    EXPECT_TRUE(p.IsNonFiniteConstant());

    p.SetExpr("1/y");  // inf at the placeholder y = 0, but run-time state
    ASSERT_TRUE(p.Evaluate(EvalContext{}));
    EXPECT_TRUE(std::isinf(p.Get()));
    EXPECT_FALSE(p.IsNonFiniteConstant());

    p.SetExpr("1/2");
    ASSERT_TRUE(p.Evaluate(EvalContext{}));
    EXPECT_FALSE(p.IsNonFiniteConstant());
}

// Audit mode only inspects: it never re-evaluates (so it is safe on a model
// whose values were written directly, e.g. by the optimiser).
TEST(ExpressionParam, AuditModeDoesNotEvaluate) {
    ExpressionParam p(2.0);
    p.SetExpr("1/0");  // not evaluated yet: cached value still 2
    NonFiniteExpressionAudit audit;
    EvalContext ctx;
    ctx.audit = &audit;
    EXPECT_TRUE(p.Evaluate(ctx));
    EXPECT_DOUBLE_EQ(p.Get(), 2.0);
    EXPECT_TRUE(audit.entries.empty());

    ASSERT_TRUE(p.Evaluate(EvalContext{}));
    p.Evaluate(ctx);
    p.Evaluate(ctx);  // recorded once
    ASSERT_EQ(audit.entries.size(), 1u);
    EXPECT_EQ(audit.entries[0].expr, "1/0");
}

namespace {
    constexpr const char* kNonFiniteDeck = "../examples/decks/Ward2024.json";

    // Loads the Ward deck, applies `edit`, then evaluates and checks it the
    // way every run path does. Returns PerformChecks' result.
    template <typename Edit>
    bool EvaluateAndCheck(Edit&& edit, std::string& err) {
        UserData ud;
        LoadFromJson(ud, kNonFiniteDeck);
        edit(ud);
        err.clear();
        if (!ud.EvaluateAllExpressions(err)) return false;
        return ud.PerformChecks(err);
    }

    CylinderBed& Cylinder(UserData& ud) {
        auto* cyl = std::get_if<CylinderBed>(&ud.bed.shape_);
        if (!cyl) throw std::runtime_error("Ward deck is expected to be a cylinder");
        return *cyl;
    }

    MassTransferCoefficientConstant& ConstantMtc(UserData& ud) {
        for (auto& [id, ads] : ud.adsorbent_directory.GetAll()) {
            for (auto& [comp, mtm] : ads.adsorption_source_.mass_transfer_models_) {
                if (!mtm) continue;
                if (auto* m = dynamic_cast<MassTransferCoefficientConstant*>(mtm->mtc_model_.get()))
                    return *m;
            }
        }
        throw std::runtime_error("Ward deck is expected to have a constant MTC");
    }
}

TEST(NonFiniteExpressions, FiniteDeckPasses) {
    std::string err;
    EXPECT_TRUE(EvaluateAndCheck([](UserData&) {}, err)) << err;

    UserData ud;
    LoadFromJson(ud, kNonFiniteDeck);
    ASSERT_TRUE(ud.EvaluateAllExpressions(err)) << err;
    EXPECT_TRUE(ud.FindNonFiniteExpressions().empty());
}

TEST(NonFiniteExpressions, InfiniteBedDiameterFailsChecksNamingIt) {
    std::string err;
    EXPECT_FALSE(EvaluateAndCheck([](UserData& ud) { Cylinder(ud).diameter.SetExpr("1/0"); }, err));
    EXPECT_NE(err.find("Non-finite"), std::string::npos) << err;
    EXPECT_NE(err.find("Bed geometry: \"1/0\" = inf"), std::string::npos) << err;
}

TEST(NonFiniteExpressions, ListsEveryOffender) {
    UserData ud;
    LoadFromJson(ud, kNonFiniteDeck);
    Cylinder(ud).diameter.SetExpr("1/0");
    auto& layer = ud.bed.layers_.layers_.at(0);
    layer.particle_diameter_.SetExpr("0/0");
    std::string err;
    ASSERT_TRUE(ud.EvaluateAllExpressions(err)) << err;

    const auto found = ud.FindNonFiniteExpressions();
    ASSERT_EQ(found.size(), 2u);
    EXPECT_EQ(found[0], "Bed geometry: \"1/0\" = inf");
    EXPECT_EQ(found[1], "Layer '" + layer.name_ + "' particle diameter: \"0/0\" = NaN");

    EXPECT_FALSE(ud.PerformChecks(err));
    EXPECT_NE(err.find(found[0]), std::string::npos) << err;
    EXPECT_NE(err.find(found[1]), std::string::npos) << err;
}

TEST(NonFiniteExpressions, InfiniteStepPressureNamesStepAndEnd) {
    UserData ud;
    LoadFromJson(ud, kNonFiniteDeck);
    std::string expected;
    for (auto& step : ud.cycle.step_list_) {
        for (auto [side, bc] : {std::pair{"top", step.top_boundary_condition_.get()},
                                std::pair{"bottom", step.bottom_boundary_condition_.get()}}) {
            auto* pbc = dynamic_cast<PressureBoundaryCondition*>(bc);
            if (!pbc || !expected.empty()) continue;
            if (auto* c = dynamic_cast<ConstantBoundaryValue*>(pbc->pressure_.get()))
                c->constant_value_.SetExpr("1/0");
            else if (auto* l = dynamic_cast<LinearBoundaryValue*>(pbc->pressure_.get()))
                l->end_value_.SetExpr("1/0");
            else
                continue;
            expected = "Step '" + step.step_name_ + "' " + side + " pressure: \"1/0\" = inf";
        }
    }
    ASSERT_FALSE(expected.empty()) << "Ward deck is expected to have a constant/linear pressure BC";
    std::string err;
    ASSERT_TRUE(ud.EvaluateAllExpressions(err)) << err;
    const auto found = ud.FindNonFiniteExpressions();
    ASSERT_EQ(found.size(), 1u);
    EXPECT_EQ(found[0], expected);
}

// `k/v` and `1/y` are infinite at the placeholder state (y = v = 0) used outside
// the run, but they are evaluated per cell during it: never flagged.
TEST(NonFiniteExpressions, StateDependentExpressionNotFlagged) {
    for (const char* expr : {"1/y", "2/v", "1/(T-298.15)"}) {
        std::string err;
        EXPECT_TRUE(EvaluateAndCheck([&](UserData& ud) { ConstantMtc(ud).constant_mtc_.SetExpr(expr); }, err))
            << expr << ": " << err;
    }
}

TEST(NonFiniteExpressions, InfiniteUserVariableFlagged) {
    std::string err;
    EXPECT_FALSE(EvaluateAndCheck([](UserData& ud) {
        UserVariable k;
        k.name = "k";
        k.value_expr = "1/0";
        ud.variable_registry.Add(std::move(k));
        // A state-dependent MTC using k is not itself flagged, but k is (user
        // variables are resolved once, so the run would use k = inf).
        ConstantMtc(ud).constant_mtc_.SetExpr("k*T");
    }, err));
    EXPECT_NE(err.find("Variable 'k': \"1/0\" = inf"), std::string::npos) << err;
    EXPECT_EQ(err.find("mass transfer"), std::string::npos) << err;
}

// In FRACTION mode step_time_ is derived by Cycle::PerformChecks; a stale
// non-finite value must not block the check that refreshes it.
TEST(NonFiniteExpressions, DerivedStepTimeIgnoredInFractionMode) {
    UserData ud;
    LoadFromJson(ud, kNonFiniteDeck);
    auto& steps = ud.cycle.step_list_;
    ASSERT_FALSE(steps.empty());
    ud.cycle.step_time_mode_ = StepTimeMode::FRACTION;
    for (auto& s : steps) s.step_fraction_ = 1.0 / static_cast<double>(steps.size());
    steps[0].step_time_ = std::numeric_limits<double>::infinity();
    EXPECT_TRUE(ud.FindNonFiniteExpressions().empty());

    ud.cycle.total_cycle_time_.SetExpr("1/0");
    std::string err;
    ASSERT_TRUE(ud.EvaluateAllExpressions(err)) << err;
    const auto found = ud.FindNonFiniteExpressions();
    ASSERT_EQ(found.size(), 1u);
    EXPECT_EQ(found[0], "Cycle total time: \"1/0\" = inf");
}

// ─── CompileForRun / EvaluateAt (per-cell residual evaluation) ──────────────

// EvaluateAt on a compiled parameter must return exactly what Evaluate does.
TEST(CompiledExpression, MatchesEvaluateBitForBit) {
    UserVariableRegistry reg;
    reg.Add({"k0", 1.7});
    const char* exprs[] = {"4.02", "12.07 [1/s]", "k0*2", "k0*exp(-300/(R*T))*(P/1e5)^0.05*(1+0.1*y)",
                           "3 [1/min]*T/300"};
    for (const char* e : exprs) {
        ExpressionParam compiled, reference;
        compiled.SetExpr(e);
        reference.SetExpr(e);
        compiled.CompileForRun(&reg);
        for (double T : {250.0, 298.15, 371.3}) {
            EvalContext ctx{&reg, T, 1.3e5, 0.0, 0.21, 0.0};
            reference.Evaluate(ctx);
            EXPECT_EQ(compiled.EvaluateAt(ctx), reference.Get()) << e << " at T=" << T;
        }
    }
}

// Editing the text after compiling falls back to Evaluate (never a stale value).
TEST(CompiledExpression, EditedTextFallsBackToEvaluate) {
    ExpressionParam p;
    p.SetExpr("2*T");
    p.CompileForRun(nullptr);
    EvalContext ctx{nullptr, 300.0, 1e5, 0.0, 0.0, 0.0};
    EXPECT_DOUBLE_EQ(p.EvaluateAt(ctx), 600.0);
    p.SetExpr("3*T");
    EXPECT_DOUBLE_EQ(p.EvaluateAt(ctx), 900.0);
    p = 5.0;
    EXPECT_DOUBLE_EQ(p.EvaluateAt(ctx), 5.0);
}

// A copy starts uncompiled (its parser would be bound to the source's storage)
// but still evaluates correctly through the fallback.
TEST(CompiledExpression, CopyEvaluatesIndependently) {
    auto src = std::make_unique<ExpressionParam>();
    src->SetExpr("T+1");
    src->CompileForRun(nullptr);
    ExpressionParam copy = *src;
    src.reset();
    EvalContext ctx{nullptr, 10.0, 1e5, 0.0, 0.0, 0.0};
    EXPECT_DOUBLE_EQ(copy.EvaluateAt(ctx), 11.0);
}

} // namespace adsim
