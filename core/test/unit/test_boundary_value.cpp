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

LogLinearBoundaryValue MakeLogLinear(double start, double end, double ramp) {
  LogLinearBoundaryValue bv;
  bv.start_value_ = start;
  bv.end_value_ = end;
  bv.ramp_time_ = ramp;
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

// ─── DoubleExponentialBoundaryValue ─────────────────────────────────────────

TEST(DoubleExponentialBoundaryValue, SumsTwoDecaysWithTheGivenWeights) {
  DoubleExponentialBoundaryValue bv;
  bv.end_value_ = 3.2e3;
  bv.rate_1_ = 0.019;
  bv.rate_2_ = 0.101;
  bv.weight_ = 0.3;
  bv.InitialiseStartValue(6.0e3);
  EXPECT_DOUBLE_EQ(bv.GetValue(0.0), 6.0e3);
  for (double t : {1.0, 20.0, 57.0, 300.0}) {
    const double expected = 3.2e3 + (6.0e3 - 3.2e3) * (0.3 * std::exp(-0.019 * t) + 0.7 * std::exp(-0.101 * t));
    EXPECT_DOUBLE_EQ(bv.GetValue(t), expected) << "t = " << t;
  }
  EXPECT_GT(bv.GetValue(1000.0), 3.2e3);
  EXPECT_LT(bv.GetValue(1000.0), 3.2e3 + 1.0);
}

TEST(DoubleExponentialBoundaryValue, WeightOneIsASingleExponential) {
  DoubleExponentialBoundaryValue bv;
  bv.end_value_ = 0.0;
  bv.rate_1_ = 0.1742;
  bv.rate_2_ = 5.0;
  bv.weight_ = 1.0;
  bv.InitialiseStartValue(1.0);
  EXPECT_DOUBLE_EQ(bv.GetValue(10.0), std::exp(-1.742));
}

TEST(DoubleExponentialBoundaryValue, ZeroRatesHoldTheStartValue) {
  DoubleExponentialBoundaryValue bv;
  bv.end_value_ = 1e4;
  bv.rate_1_ = 0.0;
  bv.rate_2_ = 0.0;
  bv.weight_ = 0.5;
  bv.InitialiseStartValue(1e5);
  for (double t : {0.0, 1.0, 100.0}) EXPECT_DOUBLE_EQ(bv.GetValue(t), 1e5) << "t = " << t;
}

// ─── LogLinearBoundaryValue ─────────────────────────────────────────────────

TEST(LogLinearBoundaryValue, RisesGeometrically) {
  // 3 kPa to 300 kPa over 20 s: a factor of 10 every 10 s.
  auto bv = MakeLogLinear(3e3, 3e5, 20.0);
  EXPECT_DOUBLE_EQ(bv.GetValue(0.0), 3e3);
  EXPECT_NEAR(bv.GetValue(10.0), 3e4, 1e-6);
  EXPECT_NEAR(bv.GetValue(5.0), 3e3 * std::sqrt(10.0), 1e-6);
  EXPECT_DOUBLE_EQ(bv.GetValue(20.0), 3e5);
}

TEST(LogLinearBoundaryValue, RelativeRateOfRiseIsConstant) {
  auto bv = MakeLogLinear(3e3, 4.6e5, 20.0);
  const double rate = std::log(4.6e5 / 3e3) / 20.0;
  const double h = 1e-4;
  for (double t : {0.5, 5.0, 10.0, 19.5}) {
    const double dlnv = (std::log(bv.GetValue(t + h)) - std::log(bv.GetValue(t - h))) / (2 * h);
    EXPECT_NEAR(dlnv, rate, 1e-6) << "t = " << t;
  }
}

TEST(LogLinearBoundaryValue, FallsGeometricallyForADecreasingRamp) {
  auto bv = MakeLogLinear(1e5, 1e3, 10.0);
  EXPECT_NEAR(bv.GetValue(5.0), 1e4, 1e-6);
  EXPECT_DOUBLE_EQ(bv.GetValue(10.0), 1e3);
}

TEST(LogLinearBoundaryValue, HoldsEndValueAndStartValue) {
  auto bv = MakeLogLinear(1.0, 100.0, 10.0);
  EXPECT_DOUBLE_EQ(bv.GetValue(-1.0), 1.0);
  EXPECT_DOUBLE_EQ(bv.GetValue(10.0), 100.0);
  EXPECT_DOUBLE_EQ(bv.GetValue(50.0), 100.0);
}

TEST(LogLinearBoundaryValue, ZeroRampGivesEndValue) {
  auto bv = MakeLogLinear(1.0, 100.0, 0.0);
  for (double t : {0.0, 1e-9, 1.0}) EXPECT_DOUBLE_EQ(bv.GetValue(t), 100.0) << "t = " << t;
}

TEST(LogLinearBoundaryValue, NonPositiveValuesRampLinearly) {
  // A flow that starts at zero has no ratio: the ramp is linear instead.
  auto from_zero = MakeLogLinear(0.0, 0.1, 10.0);
  EXPECT_DOUBLE_EQ(from_zero.GetValue(5.0), 0.05);
  auto negative = MakeLogLinear(-1.0, 1.0, 10.0);
  EXPECT_DOUBLE_EQ(negative.GetValue(5.0), 0.0);
  for (double t : {0.0, 2.5, 10.0}) EXPECT_TRUE(std::isfinite(from_zero.GetValue(t)));
}

TEST(LogLinearBoundaryValue, PressureProfileSpansStepFromBedPressure) {
  const double step_time = 20.0;
  LogLinearBoundaryValue bv;
  bv.end_value_ = 4.6e5;
  bv.start_value_.SetExpr("1e5");
  bv.ramp_time_.SetExpr("2 * 5");
  bv.InitialiseStartValue(3e3);
  bv.InitialiseStepTime(step_time);
  // The user's (saved) expressions are untouched...
  EXPECT_EQ(bv.start_value_.GetExpr(), "1e5");
  EXPECT_EQ(bv.ramp_time_.GetExpr(), "2 * 5");
  // ...while GetValue follows the runtime start and ramp.
  for (double t : {0.0, 2.66, 10.0, 20.0}) {
    const double expected = 3e3 * std::pow(4.6e5 / 3e3, t / step_time);
    EXPECT_NEAR(bv.GetValue(t), expected, 1e-6 * expected) << "t = " << t;
  }
}

TEST(LogLinearBoundaryValue, FactoryAndTemperatureProfileCreateIt) {
  auto v = CreateBoundaryValue(BoundaryValueType::LOG_LINEAR);
  ASSERT_NE(v, nullptr);
  EXPECT_EQ(v->GetType(), BoundaryValueType::LOG_LINEAR);
  EXPECT_STREQ(BoundaryValueTypeList[static_cast<int>(BoundaryValueType::LOG_LINEAR)], "Log-linear");
  auto t = CreateInletTemperatureProfile(BoundaryValueType::LOG_LINEAR, 300.0, 60.0);
  ASSERT_NE(t, nullptr);
  EXPECT_EQ(t->GetType(), BoundaryValueType::LOG_LINEAR);
  EXPECT_DOUBLE_EQ(t->GetValue(30.0), 300.0);
}

} // namespace adsim
