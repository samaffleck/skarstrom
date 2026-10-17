// adsim includes
#include "adsim/boundary_conditions/boundary_value.h"

// gtest includes
#include <gtest/gtest.h>

#include <cmath>

namespace adsim {

namespace {

LinearBoundaryValue MakeLinear(double start, double end, double ramp) {
  LinearBoundaryValue bv;
  bv.start_value_ = start;
  bv.end_value_ = end;
  bv.ramp_time_ = ramp;
  return bv;
}

ExponentialBoundaryValue MakeExponential(double start, double end, double step_time) {
  ExponentialBoundaryValue bv;
  bv.start_value_ = start;
  bv.end_value_ = end;
  bv.step_time_ = step_time;
  return bv;
}

} // namespace

// ─── LinearBoundaryValue ────────────────────────────────────────────────────

TEST(LinearBoundaryValue, RampsLinearlyWithinRamp) {
  auto bv = MakeLinear(0.0, 0.1, 10.0);
  EXPECT_DOUBLE_EQ(bv.GetValue(0.0), 0.0);
  EXPECT_DOUBLE_EQ(bv.GetValue(5.0), 0.05);
}

TEST(LinearBoundaryValue, HoldsEndValueAtAndAfterRampEnd) {
  // Flow-profile scenario from E-05: ramp 10 s inside a 100 s step.
  auto bv = MakeLinear(0.0, 0.1, 10.0);
  EXPECT_DOUBLE_EQ(bv.GetValue(10.0), 0.1);
  EXPECT_DOUBLE_EQ(bv.GetValue(50.0), 0.1);
  EXPECT_DOUBLE_EQ(bv.GetValue(100.0), 0.1);
}

TEST(LinearBoundaryValue, ZeroRampGivesEndValue) {
  auto bv = MakeLinear(0.0, 0.1, 0.0);
  for (double t : {0.0, 1e-9, 1.0, 100.0}) {
    const double v = bv.GetValue(t);
    EXPECT_TRUE(std::isfinite(v)) << "t = " << t;
    EXPECT_DOUBLE_EQ(v, 0.1) << "t = " << t;
  }
}

TEST(LinearBoundaryValue, NegativeTimeHoldsStartValue) {
  auto bv = MakeLinear(1.0, 2.0, 10.0);
  EXPECT_DOUBLE_EQ(bv.GetValue(-1.0), 1.0);
}

TEST(LinearBoundaryValue, PressureProfileUnchangedWithinStep) {
  // Pressure BCs call InitialiseStartValue/InitialiseStepTime, so ramp == step.
  const double step_time = 300.0;
  LinearBoundaryValue bv;
  bv.end_value_ = 10e5;
  bv.InitialiseStartValue(1e5);
  bv.InitialiseStepTime(step_time);
  for (double t : {0.0, 1.0, 150.0, 299.0, 300.0}) {
    const double expected = 1e5 + (10e5 - 1e5) * (t / step_time);
    EXPECT_DOUBLE_EQ(bv.GetValue(t), expected) << "t = " << t;
  }
}

TEST(LinearBoundaryValue, RuntimeInitialiseLeavesUserExpressionsUntouched) {
  // L-05: pressure/valve parents' runtime start/ramp must not rewrite the
  // user's (saved) expressions.
  LinearBoundaryValue bv;
  bv.start_value_.SetExpr("3e5 / 2");
  bv.ramp_time_.SetExpr("2 * 5");
  bv.end_value_ = 10e5;
  bv.InitialiseStartValue(1e5);
  bv.InitialiseStepTime(300.0);
  EXPECT_EQ(bv.start_value_.GetExpr(), "3e5 / 2");
  EXPECT_EQ(bv.ramp_time_.GetExpr(), "2 * 5");
  // ...while GetValue follows the runtime start/ramp.
  EXPECT_DOUBLE_EQ(bv.GetValue(0.0), 1e5);
  EXPECT_DOUBLE_EQ(bv.GetValue(150.0), 5.5e5);
  EXPECT_DOUBLE_EQ(bv.GetValue(300.0), 10e5);
}

// ─── ExponentialBoundaryValue ───────────────────────────────────────────────

TEST(ExponentialBoundaryValue, ReachesNinetyNinePercentAtStepTime) {
  auto bv = MakeExponential(0.0, 1.0, 100.0);
  EXPECT_DOUBLE_EQ(bv.GetValue(0.0), 0.0);
  EXPECT_NEAR(bv.GetValue(100.0), 0.99, 1e-12);
  // Saturates towards (but never past) the end value after the step time.
  const double late = bv.GetValue(1000.0);
  EXPECT_GT(late, 0.99);
  EXPECT_LE(late, 1.0);
}

TEST(ExponentialBoundaryValue, ZeroTimeConstantIsFinite) {
  auto bv = MakeExponential(1e5, 4e4, 0.0);
  for (double t : {0.0, 1e-9, 1.0, 100.0}) {
    const double v = bv.GetValue(t);
    EXPECT_TRUE(std::isfinite(v)) << "t = " << t;
    EXPECT_DOUBLE_EQ(v, 4e4) << "t = " << t;
  }
}

TEST(ExponentialBoundaryValue, PressureProfileUnchangedWithinStep) {
  const double step_time = 200.0;
  ExponentialBoundaryValue bv;
  bv.end_value_ = 1e4;
  bv.InitialiseStartValue(1e5);
  bv.InitialiseStepTime(step_time);
  for (double t : {0.0, 1.0, 100.0, 200.0}) {
    const double expected = 1e5 + (1e4 - 1e5) *
      (1.0 - std::exp(-(std::log(100.0) / step_time) * t));
    EXPECT_DOUBLE_EQ(bv.GetValue(t), expected) << "t = " << t;
  }
}

TEST(ExponentialBoundaryValue, RuntimeInitialiseLeavesUserExpressionsUntouched) {
  ExponentialBoundaryValue bv;
  bv.start_value_.SetExpr("3e5 / 2");
  bv.step_time_.SetExpr("2 * 5");
  bv.end_value_ = 1e4;
  bv.InitialiseStartValue(1e5);
  bv.InitialiseStepTime(200.0);
  EXPECT_EQ(bv.start_value_.GetExpr(), "3e5 / 2");
  EXPECT_EQ(bv.step_time_.GetExpr(), "2 * 5");
  EXPECT_DOUBLE_EQ(bv.GetValue(0.0), 1e5);
  EXPECT_NEAR(bv.GetValue(200.0), 1e5 + (1e4 - 1e5) * 0.99, 1e-6);
}

} // namespace adsim
