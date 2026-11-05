// adsim includes
#include "adsim/boundary_conditions/boundary_value.h"

// gtest includes
#include <gtest/gtest.h>

#include <cereal/archives/json.hpp>
#include <cereal/types/memory.hpp>

#include <cmath>
#include <memory>
#include <regex>
#include <sstream>
#include <string>

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

// The factories the document layers use (the web contract builds a profile
// from its kind string through CreateBoundaryValue; a null here collapses the
// profile to a constant downstream).
TEST(DoubleExponentialBoundaryValue, FactoryAndTemperatureProfileCreateIt) {
  auto v = CreateBoundaryValue(BoundaryValueType::DOUBLE_EXPONENTIAL);
  ASSERT_NE(v, nullptr);
  EXPECT_EQ(v->GetType(), BoundaryValueType::DOUBLE_EXPONENTIAL);
  EXPECT_STREQ(BoundaryValueTypeList[static_cast<int>(BoundaryValueType::DOUBLE_EXPONENTIAL)], "Double exponential");
  auto t = CreateInletTemperatureProfile(BoundaryValueType::DOUBLE_EXPONENTIAL, 300.0, 60.0);
  ASSERT_NE(t, nullptr);
  EXPECT_EQ(t->GetType(), BoundaryValueType::DOUBLE_EXPONENTIAL);
  EXPECT_DOUBLE_EQ(t->GetValue(30.0), 300.0);
}

// ─── Pressure-profile ramp time (use_step_time_) ────────────────────────────

// Off: a pressure-driven parent still sets the start value, but the ramp lasts
// the profile's own ramp time and then holds the end value.
TEST(LinearBoundaryValue, OwnRampTimeOverridesTheStepTime) {
  LinearBoundaryValue bv;
  bv.end_value_ = 10e5;
  bv.ramp_time_ = 100.0;
  bv.use_step_time_ = false;
  bv.InitialiseStartValue(1e5);
  bv.InitialiseStepTime(300.0);
  EXPECT_DOUBLE_EQ(bv.GetValue(50.0), 5.5e5);
  EXPECT_DOUBLE_EQ(bv.GetValue(100.0), 10e5);
  EXPECT_DOUBLE_EQ(bv.GetValue(250.0), 10e5);
}

// Switching the flag takes effect at the next step start (no stale override).
TEST(LinearBoundaryValue, StepTimeSwitchAppliesAtTheNextStepStart) {
  LinearBoundaryValue bv;
  bv.end_value_ = 2.0;
  bv.ramp_time_ = 10.0;
  bv.InitialiseStartValue(1.0);
  bv.InitialiseStepTime(100.0);
  EXPECT_DOUBLE_EQ(bv.GetValue(10.0), 1.1);
  bv.use_step_time_ = false;
  bv.InitialiseStepTime(100.0);
  EXPECT_DOUBLE_EQ(bv.GetValue(10.0), 2.0);
  bv.use_step_time_ = true;
  bv.InitialiseStepTime(100.0);
  EXPECT_DOUBLE_EQ(bv.GetValue(10.0), 1.1);
}

TEST(ExponentialBoundaryValue, OwnRampTimeOverridesTheStepTime) {
  ExponentialBoundaryValue bv;
  bv.end_value_ = 0.0;
  bv.step_time_ = 20.0;
  bv.use_step_time_ = false;
  bv.InitialiseStartValue(1.0);
  bv.InitialiseStepTime(200.0);
  EXPECT_NEAR(bv.GetValue(20.0), 0.01, 1e-12);  // 99 % covered at its own 20 s
}

TEST(LogLinearBoundaryValue, OwnRampTimeOverridesTheStepTime) {
  LogLinearBoundaryValue bv;
  bv.end_value_ = 100.0;
  bv.ramp_time_ = 10.0;
  bv.use_step_time_ = false;
  bv.InitialiseStartValue(1.0);
  bv.InitialiseStepTime(40.0);
  EXPECT_NEAR(bv.GetValue(5.0), 10.0, 1e-12);
  EXPECT_DOUBLE_EQ(bv.GetValue(30.0), 100.0);
}

// ─── Exponential factor ─────────────────────────────────────────────────────

// v(t) = v1 + (v0 - v1) F^(-t / t_r): the change left at t_r is 1/F of the initial.
TEST(ExponentialBoundaryValue, FactorSetsTheChangeLeftAtTheRampTime) {
  for (double F : {2.0, 10.0, 100.0, 1000.0}) {
    auto bv = MakeExponential(5.0, 1.0, 40.0);
    bv.factor_ = F;
    EXPECT_NEAR(bv.GetValue(40.0), 1.0 + 4.0 / F, 1e-12) << "F = " << F;
    EXPECT_NEAR(bv.GetValue(20.0), 1.0 + 4.0 / std::sqrt(F), 1e-12) << "F = " << F;
    EXPECT_DOUBLE_EQ(bv.GetValue(0.0), 5.0) << "F = " << F;
  }
}

TEST(ExponentialBoundaryValue, DefaultFactorIsOneHundred) {
  ExponentialBoundaryValue bv;
  EXPECT_DOUBLE_EQ(bv.factor_.Get(), 100.0);
  EXPECT_TRUE(bv.use_step_time_);
}

// F <= 1 is rejected by the checks; GetValue stays finite and holds the start.
TEST(ExponentialBoundaryValue, FactorAtOrBelowOneHoldsTheStartValue) {
  for (double F : {1.0, 0.5, 0.0, -3.0}) {
    auto bv = MakeExponential(5.0, 1.0, 40.0);
    bv.factor_ = F;
    EXPECT_DOUBLE_EQ(bv.GetValue(40.0), 5.0) << "F = " << F;
    std::string err;
    EXPECT_FALSE(CheckProfileFactor(&bv, "test profile", err)) << "F = " << F;
    EXPECT_NE(err.find("greater than 1"), std::string::npos) << err;
  }
  auto ok = MakeExponential(5.0, 1.0, 40.0);
  ok.factor_ = 1.5;
  std::string err;
  EXPECT_TRUE(CheckProfileFactor(&ok, "test profile", err)) << err;
}

// ─── CheckPressureProfile ───────────────────────────────────────────────────

TEST(CheckPressureProfile, OwnRampTimeMustFitInTheStep) {
  LinearBoundaryValue lin;
  lin.ramp_time_ = 30.0;
  std::string err;
  // Over the step time: the user's ramp time is ignored.
  EXPECT_TRUE(CheckPressureProfile(&lin, 20.0, "p", err)) << err;
  lin.use_step_time_ = false;
  EXPECT_TRUE(CheckPressureProfile(&lin, 30.0, "p", err)) << err;  // equal is fine
  EXPECT_TRUE(CheckPressureProfile(&lin, 40.0, "p", err)) << err;
  EXPECT_FALSE(CheckPressureProfile(&lin, 20.0, "p", err));
  EXPECT_NE(err.find("longer than the step"), std::string::npos) << err;

  err.clear();
  lin.ramp_time_ = -1.0;
  EXPECT_FALSE(CheckPressureProfile(&lin, 20.0, "p", err));
  EXPECT_NE(err.find("zero or positive"), std::string::npos) << err;

  LogLinearBoundaryValue lg;
  lg.ramp_time_ = 30.0;
  lg.use_step_time_ = false;
  err.clear();
  EXPECT_FALSE(CheckPressureProfile(&lg, 20.0, "p", err));

  ExponentialBoundaryValue ex;
  ex.step_time_ = 30.0;
  ex.use_step_time_ = false;
  err.clear();
  EXPECT_FALSE(CheckPressureProfile(&ex, 20.0, "p", err));
  ex.step_time_ = 10.0;
  ex.factor_ = 1.0;
  err.clear();
  EXPECT_FALSE(CheckPressureProfile(&ex, 20.0, "p", err));  // the factor is checked too

  ConstantBoundaryValue c(1e5);
  DoubleExponentialBoundaryValue de;
  err.clear();
  EXPECT_TRUE(CheckPressureProfile(&c, 1.0, "p", err));
  EXPECT_TRUE(CheckPressureProfile(&de, 1.0, "p", err));
  EXPECT_TRUE(CheckPressureProfile(nullptr, 1.0, "p", err));
}

// ─── Serialisation ──────────────────────────────────────────────────────────

namespace {

std::string Save(std::unique_ptr<IBoundaryValue> v) {
  std::ostringstream os;
  {
    cereal::JSONOutputArchive ar(os);
    ar(cereal::make_nvp("profile", v));
  }
  return os.str();
}

std::unique_ptr<IBoundaryValue> Load(const std::string& json) {
  std::istringstream is(json);
  std::unique_ptr<IBoundaryValue> v;
  cereal::JSONInputArchive ar(is);
  ar(cereal::make_nvp("profile", v));
  return v;
}

// The JSON as an older engine wrote it: the class at `version`, without `fields`.
std::string AsOlderVersion(std::string json, int version, std::initializer_list<const char*> fields) {
  json = std::regex_replace(json, std::regex("\"cereal_class_version\": [0-9]+"),
                            "\"cereal_class_version\": " + std::to_string(version), std::regex_constants::format_first_only);
  for (const char* f : fields) {
    // A field and its value (a bool or an ExpressionParam object), with the comma before it.
    json = std::regex_replace(json, std::regex(",\\s*\"" + std::string(f) + "\": (true|false|\\{[^}]*\\})"), "");
  }
  return json;
}

}  // namespace

TEST(BoundaryValueSerialisation, RampTimeSwitchAndFactorRoundTrip) {
  auto lin = std::make_unique<LinearBoundaryValue>();
  lin->use_step_time_ = false;
  lin->ramp_time_ = 12.5;
  auto back = Load(Save(std::move(lin)));
  auto* l = dynamic_cast<LinearBoundaryValue*>(back.get());
  ASSERT_NE(l, nullptr);
  EXPECT_FALSE(l->use_step_time_);
  EXPECT_DOUBLE_EQ(l->ramp_time_.Get(), 12.5);

  auto ex = std::make_unique<ExponentialBoundaryValue>();
  ex->use_step_time_ = false;
  ex->factor_.SetExpr("2 * 500");
  auto back_ex = Load(Save(std::move(ex)));
  auto* e = dynamic_cast<ExponentialBoundaryValue*>(back_ex.get());
  ASSERT_NE(e, nullptr);
  EXPECT_FALSE(e->use_step_time_);
  EXPECT_EQ(e->factor_.GetExpr(), "2 * 500");

  auto lg = std::make_unique<LogLinearBoundaryValue>();
  lg->use_step_time_ = false;
  auto back_lg = Load(Save(std::move(lg)));
  auto* g = dynamic_cast<LogLinearBoundaryValue*>(back_lg.get());
  ASSERT_NE(g, nullptr);
  EXPECT_FALSE(g->use_step_time_);
}

// Decks written before the switch and the factor existed ramp over the step
// time with F = 100, as they always did.
TEST(BoundaryValueSerialisation, OlderProfilesRampOverTheStepWithFactor100) {
  {
    auto lin = std::make_unique<LinearBoundaryValue>();
    lin->use_step_time_ = false;
    const std::string old = AsOlderVersion(Save(std::move(lin)), 3, {"Use_Step_Time"});
    ASSERT_EQ(old.find("Use_Step_Time"), std::string::npos) << old;
    auto v = Load(old);
    auto* l = dynamic_cast<LinearBoundaryValue*>(v.get());
    ASSERT_NE(l, nullptr);
    EXPECT_TRUE(l->use_step_time_);
  }
  {
    auto ex = std::make_unique<ExponentialBoundaryValue>();
    ex->use_step_time_ = false;
    ex->factor_ = 7.0;
    const std::string old = AsOlderVersion(Save(std::move(ex)), 3, {"Factor", "Use_Step_Time"});
    ASSERT_EQ(old.find("Factor"), std::string::npos) << old;
    auto v = Load(old);
    auto* e = dynamic_cast<ExponentialBoundaryValue*>(v.get());
    ASSERT_NE(e, nullptr);
    EXPECT_TRUE(e->use_step_time_);
    EXPECT_DOUBLE_EQ(e->factor_.Get(), 100.0);
  }
  {
    auto lg = std::make_unique<LogLinearBoundaryValue>();
    lg->use_step_time_ = false;
    auto v = Load(AsOlderVersion(Save(std::move(lg)), 1, {"Use_Step_Time"}));
    auto* g = dynamic_cast<LogLinearBoundaryValue*>(v.get());
    ASSERT_NE(g, nullptr);
    EXPECT_TRUE(g->use_step_time_);
  }
}

// Every listed kind has a factory entry.
TEST(BoundaryValueFactory, MakesEveryListedKind) {
  const int n = static_cast<int>(sizeof(BoundaryValueTypeList) / sizeof(BoundaryValueTypeList[0]));
  for (int i = 0; i < n; ++i) {
    auto v = CreateBoundaryValue(static_cast<BoundaryValueType>(i));
    ASSERT_NE(v, nullptr) << BoundaryValueTypeList[i];
    EXPECT_EQ(static_cast<int>(v->GetType()), i) << BoundaryValueTypeList[i];
  }
}

} // namespace adsim
