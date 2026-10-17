// Tests for the pure parameterisation helper functions defined in
// core/include/adsim/solver/parameterisation.h

// STL includes
#include <cmath>
#include <map>
#include <vector>

// adsim includes
#include "adsim/solver/parameterisation.h"

// GTest
#include <gtest/gtest.h>


namespace adsim {

// =============================================================================
// InterpLinear
// =============================================================================

TEST(InterpLinear, ExactGridPoints)
{
  std::vector<double> x = {0.0, 1.0, 2.0, 3.0};
  std::vector<double> y = {10.0, 20.0, 30.0, 40.0};

  EXPECT_DOUBLE_EQ(InterpLinear(x, y, 0.0), 10.0);
  EXPECT_DOUBLE_EQ(InterpLinear(x, y, 1.0), 20.0);
  EXPECT_DOUBLE_EQ(InterpLinear(x, y, 2.0), 30.0);
  EXPECT_DOUBLE_EQ(InterpLinear(x, y, 3.0), 40.0);
}

TEST(InterpLinear, MidpointBetweenPoints)
{
  std::vector<double> x = {0.0, 2.0, 4.0};
  std::vector<double> y = {0.0, 4.0, 8.0};  // linear y = 2x

  EXPECT_NEAR(InterpLinear(x, y, 1.0), 2.0, 1e-12);
  EXPECT_NEAR(InterpLinear(x, y, 3.0), 6.0, 1e-12);
}

TEST(InterpLinear, ClampsBelow)
{
  std::vector<double> x = {1.0, 2.0, 3.0};
  std::vector<double> y = {10.0, 20.0, 30.0};

  // Query below x[0] → returns y[0]
  EXPECT_DOUBLE_EQ(InterpLinear(x, y, -5.0), 10.0);
  EXPECT_DOUBLE_EQ(InterpLinear(x, y,  0.0), 10.0);
}

TEST(InterpLinear, ClampsAbove)
{
  std::vector<double> x = {1.0, 2.0, 3.0};
  std::vector<double> y = {10.0, 20.0, 30.0};

  // Query above x.back() → returns y.back()
  EXPECT_DOUBLE_EQ(InterpLinear(x, y, 4.0), 30.0);
  EXPECT_DOUBLE_EQ(InterpLinear(x, y, 100.0), 30.0);
}

TEST(InterpLinear, EmptyVector)
{
  std::vector<double> x, y;
  EXPECT_DOUBLE_EQ(InterpLinear(x, y, 1.0), 0.0);
}

TEST(InterpLinear, SingleElement)
{
  std::vector<double> x = {5.0};
  std::vector<double> y = {42.0};

  EXPECT_DOUBLE_EQ(InterpLinear(x, y, 0.0),  42.0);
  EXPECT_DOUBLE_EQ(InterpLinear(x, y, 5.0),  42.0);
  EXPECT_DOUBLE_EQ(InterpLinear(x, y, 99.0), 42.0);
}

TEST(InterpLinear, NonUniformSpacing)
{
  // x = {0, 1, 10}, y = {0, 1, 10} → f(x) = x
  std::vector<double> x = {0.0, 1.0, 10.0};
  std::vector<double> y = {0.0, 1.0, 10.0};

  EXPECT_NEAR(InterpLinear(x, y, 0.5), 0.5, 1e-12);
  EXPECT_NEAR(InterpLinear(x, y, 5.5), 5.5, 1e-12);
}

// =============================================================================
// ComputeRMSE
// =============================================================================

TEST(ComputeRMSE, ZeroErrorWhenDataMatchesExactly)
{
  // sim and exp are identical grid → RMSE should be 0
  std::vector<double> sim_t = {0.0, 1.0, 2.0, 3.0, 4.0};
  std::vector<double> sim_v = {0.0, 1.0, 4.0, 9.0, 16.0};  // y = x²

  XYData exp;
  exp.x = {0.0, 1.0, 2.0, 3.0, 4.0};
  exp.y = {0.0, 1.0, 4.0, 9.0, 16.0};

  EXPECT_NEAR(ComputeRMSE(sim_t, sim_v, exp), 0.0, 1e-12);
}

TEST(ComputeRMSE, ConstantOffset)
{
  // sim = x, exp = x + 2  →  residuals all = -2, RMSE = 2
  std::vector<double> sim_t = {0.0, 1.0, 2.0, 3.0};
  std::vector<double> sim_v = {0.0, 1.0, 2.0, 3.0};

  XYData exp;
  exp.x = {0.0, 1.0, 2.0, 3.0};
  exp.y = {2.0, 3.0, 4.0, 5.0};  // sim + 2

  EXPECT_NEAR(ComputeRMSE(sim_t, sim_v, exp), 2.0, 1e-12);
}

TEST(ComputeRMSE, KnownMixedResiduals)
{
  // residuals: -1, +1, -1, +1 → SSE = 4, RMSE = sqrt(4/4) = 1
  std::vector<double> sim_t = {0.0, 1.0, 2.0, 3.0};
  std::vector<double> sim_v = {0.0, 2.0, 0.0, 2.0};

  XYData exp;
  exp.x = {0.0, 1.0, 2.0, 3.0};
  exp.y = {1.0, 1.0, 1.0, 1.0};

  EXPECT_NEAR(ComputeRMSE(sim_t, sim_v, exp), 1.0, 1e-12);
}

TEST(ComputeRMSE, EmptySimTimesReturns1e10)
{
  std::vector<double> sim_t, sim_v;
  XYData exp;
  exp.x = {1.0};
  exp.y = {2.0};

  EXPECT_EQ(ComputeRMSE(sim_t, sim_v, exp), 1e10);
}

TEST(ComputeRMSE, EmptyExpDataReturns1e10)
{
  std::vector<double> sim_t = {0.0, 1.0};
  std::vector<double> sim_v = {0.0, 1.0};
  XYData exp;  // empty x, y

  EXPECT_EQ(ComputeRMSE(sim_t, sim_v, exp), 1e10);
}

TEST(ComputeRMSE, TimeUnitConversion)
{
  // sim in seconds, exp x-axis in minutes (mul = 60)
  // sim: t=0→0, t=60→1  (linear)
  // exp point at x=1 min = 60 s → sim_y = 1, exp_y = 1 → RMSE = 0
  std::vector<double> sim_t = {0.0, 60.0};
  std::vector<double> sim_v = {0.0, 1.0};

  XYData exp;
  exp.x = {1.0};    // 1 minute
  exp.y = {1.0};

  EXPECT_NEAR(ComputeRMSE(sim_t, sim_v, exp, /*time_unit_mul=*/60.0), 0.0, 1e-12);
}

TEST(ComputeRMSE, InterpolatedSimValue)
{
  // sim: t=0→0, t=10→10 (linear)
  // exp at t=5 → sim interpolates to 5, exp_y = 3 → residual = 5-3 = 2, RMSE = 2
  std::vector<double> sim_t = {0.0, 10.0};
  std::vector<double> sim_v = {0.0, 10.0};

  XYData exp;
  exp.x = {5.0};
  exp.y = {3.0};

  EXPECT_NEAR(ComputeRMSE(sim_t, sim_v, exp), 2.0, 1e-12);
}

TEST(ComputeRMSE, ExperimentalPointOutsideSimRange)
{
  // exp at t=100 > sim max t=10 → clamps to sim_v.back() = 10
  // residual = 10 - 5 = 5, RMSE = 5
  std::vector<double> sim_t = {0.0, 10.0};
  std::vector<double> sim_v = {0.0, 10.0};

  XYData exp;
  exp.x = {100.0};  // outside range
  exp.y = {5.0};

  EXPECT_NEAR(ComputeRMSE(sim_t, sim_v, exp), 5.0, 1e-12);
}

TEST(ComputeRMSE, MultipleExperimentalPoints)
{
  // sim: t=0→0, t=4→4 (linear, y=t)
  // exp points at t=1,2,3 with y=1,2,3 → all residuals = 0
  std::vector<double> sim_t = {0.0, 4.0};
  std::vector<double> sim_v = {0.0, 4.0};

  XYData exp;
  exp.x = {1.0, 2.0, 3.0};
  exp.y = {1.0, 2.0, 3.0};

  EXPECT_NEAR(ComputeRMSE(sim_t, sim_v, exp), 0.0, 1e-12);
}

// =============================================================================
// ComputeKPIRMSE edge cases (no full simulation needed)
// =============================================================================

TEST(ComputeKPIRMSE, ReturnsLargePenaltyWhenNoCyclesMatch)
{
  // Build a minimal EquationLog (empty) and BoundaryLog
  // ComputeKPIRMSE should return 1e10 when no cycles match
  EquationLog eq_log;
  BoundaryLog b_log;

  XYData exp;
  exp.x = {0.0, 1.0, 2.0};  // cycle numbers
  exp.y = {0.9, 0.91, 0.92};

  // The log is empty, so sim_log.find() always fails → n = 0 → returns 1e10
  double rmse = ComputeKPIRMSE(
    eq_log, b_log, exp,
    /*cycle_offset=*/0,
    KPIVariableType::PURITY,
    "CO2",
    /*feed_endpoints=*/{ { 0, "Stream 1", BoundaryLocation::BOTTOM } },
    /*product_endpoints=*/{});

  EXPECT_EQ(rmse, 1e10);
}

TEST(ComputeKPIRMSE, ReturnsLargePenaltyWhenExpDataEmpty)
{
  EquationLog eq_log;
  BoundaryLog b_log;

  XYData exp;  // empty

  double rmse = ComputeKPIRMSE(
    eq_log, b_log, exp,
    0, KPIVariableType::PURITY, "CO2",
    {}, {});

  EXPECT_EQ(rmse, 1e10);
}

// =============================================================================
// CycleStartTime / Extract*TimeSeries edge cases (E-10)
// =============================================================================

TEST(CycleStartTime, EmptyCycleReturnsZero)
{
  std::map<int, Log> cycle_log;
  EXPECT_EQ(CycleStartTime(cycle_log), 0.0);
}

TEST(CycleStartTime, UsesFirstLoggedAbsoluteTime)
{
  // Logged times are absolute IDA times: the third cycle of a 100 s cycle
  // starts at t = 200 s, and its steps are laid out contiguously after that.
  std::map<int, Log> cycle_log;
  cycle_log[0].step_data.times = {200.0, 225.0, 250.0};
  cycle_log[1].step_data.times = {250.0, 275.0, 300.0};
  EXPECT_DOUBLE_EQ(CycleStartTime(cycle_log), 200.0);
}

TEST(CycleStartTime, SkipsStepsWithNoLoggedTimes)
{
  std::map<int, Log> cycle_log;
  cycle_log[0];                                  // nothing logged
  cycle_log[1].step_data.times = {150.0, 175.0};
  EXPECT_DOUBLE_EQ(CycleStartTime(cycle_log), 150.0);
}

TEST(ExtractTimeSeries, EmptyLogReturnsFalse)
{
  EquationLog eq_log;
  std::vector<double> t{1.0}, v{1.0};
  EXPECT_FALSE(ExtractBulkTimeSeries(eq_log, VariableType::PRESSURE,
                                     ParamSpatialLoc::Bottom, 0.0, t, v));
  EXPECT_FALSE(ExtractComponentTimeSeries(eq_log, ComponentVariableType::MOLE_FRACTION,
                                          "CO2", ParamSpatialLoc::Bottom, 0.0, t, v));
  EXPECT_TRUE(t.empty());
  EXPECT_TRUE(v.empty());
}

} // namespace adsim
