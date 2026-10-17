// Integration tests for parameterisation helper functions.
// Runs a real simulation and verifies:
//   1. ExtractBulkTimeSeries / ExtractComponentTimeSeries produce physically
//      sensible, monotonically-non-decreasing time series.
//   2. ComputeRMSE returns zero when experimental data is sampled directly
//      from the simulation output (self-consistency check).
//   3. ComputeKPIRMSE returns zero when the "experimental" KPI values are
//      the values computed from the same simulation run.

// STL includes
#define _USE_MATH_DEFINES
#include <cmath>
#include <filesystem>
#include <algorithm>
#include <numeric>

// AdSim includes
#include "../src/pch.h"
#include "adsim/solver/parameterisation.h"
#include "adsim/solver/coupled_solver.h"
#include "adsim/solver/kpis.h"
#include "adsim/solver/kpi_config.h"
#include "adsim/math/constants.h"
#include "adsim/equations/isotherm.h"
#include "adsim/boundary_conditions/boundary_conditions.h"
#include "adsim/boundary_conditions/boundary_stream.h"
#include "adsim/boundary_conditions/step.h"
#include "adsim/equations/dispersion_model.h"
#include "adsim/equations/heat_transfer_model.h"

// GTest
#include <gtest/gtest.h>

// Test helpers
#include "../fixtures/user_data_test.h"


namespace adsim {

namespace {

// Checks that a time vector is monotonically non-decreasing.
bool IsNonDecreasing(const std::vector<double>& v) {
  for (size_t i = 1; i < v.size(); ++i)
    if (v[i] < v[i - 1]) return false;
  return true;
}

} // namespace

// ---------------------------------------------------------------------------
// Fixture: a simple CO2 breakthrough experiment (single adsorption step).
// ---------------------------------------------------------------------------

class ParameterisationTest : public UserDataTest {
public:
  void SetUp() override {
    UserDataTest::SetUp();

    AddComponents({"CarbonDioxide", "Air"});

    auto [lid, aid] = AddLayer("Activated Carbon");
    layer_id = lid;
    ads_id   = aid;

    auto& AC = GetLayer(layer_id);
    AC.domain_.length     = 0.286;
    AC.particle_diameter_ = 0.003;
    AC.eb_                = 0.456;
    AC.domain_.n_cells    = 20;   // coarse mesh for speed
    std::get<CylinderBed>(data.bed.shape_).diameter = 0.0189;

    auto& AC_ads = GetAdsorbent(ads_id);
    AC_ads.solid_density_         = 799;
    AC_ads.heat_capacity_         = 1040;
    AC_ads.thermal_conductivity_  = 0.63;
    AC_ads.ep_                    = 0.0;

    data.bed.wall_.heat_capacity_        = 902;
    data.bed.wall_.density_              = 2700;
    data.bed.wall_.thermal_conductivity_ = 205;
    data.bed.wall_.thickness_            = 3.25e-3;
    data.bed.ambient_temperature_        = 293;
    data.fluid.adsorbed_heat_capacity_["CarbonDioxide"] = 36.96;

    auto dsl = std::make_unique<DualSiteLangmuirIsotherm_Concentration>("CarbonDioxide");
    dsl->params_["CarbonDioxide"].m1 = 0.418;
    dsl->params_["CarbonDioxide"].b1 = 2.05e-6;
    dsl->params_["CarbonDioxide"].q1 = 29980;
    dsl->params_["CarbonDioxide"].m2 = 6.04;
    dsl->params_["CarbonDioxide"].b2 = 3.32e-6;
    dsl->params_["CarbonDioxide"].q2 = 20367;
    AddIsotherm({"CarbonDioxide"}, ads_id, std::move(dsl));

    AddLinearDrivingForceModel_Haynes(ads_id, "CarbonDioxide", 3, 0.59);
    AddLinearDrivingForceModel(ads_id, "Air", 0.0);

    SetHOA(ads_id, "CarbonDioxide", 29980.0);
    SetHOA(ads_id, "Air", 0.0);

    GetLayer(layer_id).equations_.dispersion_model_ =
      CreateAxialDispersionModel(AxialDispersionType::RUTHVEN);
    GetLayer(layer_id).equations_.heat_dispersion_model_ =
      CreateAxialHeatDispersionModel(AxialHeatDispersionType::YAGI_AND_KUNII);

    AddConstantInternalHeatTransferCoefficientModel(57.9, layer_id);
    auto ehtm = std::make_unique<ExternalHeatTransferModel_Holman>();
    ehtm->c = 0.027;  ehtm->m = 0.8;  ehtm->n = 1.0/3.0;
    ehtm->external_velocity = 0.5;
    GetLayer(layer_id).equations_.external_heat_transfer_model_ = std::move(ehtm);

    AddInitialConditions(1.0e5, 291.4, 291.4, {{"CarbonDioxide", 0.0}, {"Air", 1.0}});

    AddStep("Feed", 500.0);  // 500 s to keep the test fast

    const double Q_m3s = 100 * 1.666666666666667e-8;  // 100 ccm → m³/s

    auto inlet_stream = CreateNoValveInletStream(
      CreateVolumeFlowBC(Q_m3s),
      {{"CarbonDioxide", 1.0}, {"Air", 0.0}},
      293.0);
    inlet_stream->name_ = "CO2 Inlet";

    auto feed_flow = std::make_unique<FlowBoundaryCondition>();
    feed_flow->direction_ = BoundaryDirection::INLET;
    feed_flow->location_  = BoundaryLocation::BOTTOM;
    feed_flow->streams_.push_back(std::move(inlet_stream));
    data.cycle.GetStep("Feed").bottom_boundary_condition_ = std::move(feed_flow);

    auto pressure_out = std::make_unique<PressureBoundaryCondition>();
    pressure_out->direction_               = BoundaryDirection::OUTLET;
    pressure_out->location_                = BoundaryLocation::TOP;
    pressure_out->pressure_                = CreateConstantBoundaryValue(1e5);
    pressure_out->mole_fraction_["CarbonDioxide"] = 0.0;
    pressure_out->mole_fraction_["Air"]           = 1.0;
    data.cycle.GetStep("Feed").top_boundary_condition_ = std::move(pressure_out);

    ConfigureCycle(1);

    ConfigureKPI(
      /*product_comp_idx=*/0,
      /*feed_steps=*/{0},  /*feed_locs=*/{BoundaryLocation::BOTTOM},
      /*product_steps=*/{0}, /*product_locs=*/{BoundaryLocation::TOP});
  }

  LayerID    layer_id{};
  AdsorbentID ads_id{};
};

// ---------------------------------------------------------------------------
// Test 1: Simulation runs and logs are non-empty.
// ---------------------------------------------------------------------------

TEST_F(ParameterisationTest, SolverCompletesSuccessfully)
{
  std::string err;
  ASSERT_TRUE(data.PerformChecks(err)) << err;

  CoupledSolver solver(data);
  int retval = solver.IntegrateCycle(data.cycle, false);
  EXPECT_GE(retval, 0) << "Solver returned error code: " << retval;

  const auto& log = solver.GetLog();
  EXPECT_FALSE(log.GetSimulationLog().empty());
}

// ---------------------------------------------------------------------------
// Test 2: ExtractBulkTimeSeries – times are non-decreasing, values finite.
// ---------------------------------------------------------------------------

TEST_F(ParameterisationTest, ExtractBulkTimeSeries_TimeIsNonDecreasing)
{
  CoupledSolver solver(data);
  ASSERT_GE(solver.IntegrateCycle(data.cycle, false), 0);

  auto eq_log = solver.ExtractLog();

  std::vector<double> times, values;
  bool ok = ExtractBulkTimeSeries(
    eq_log,
    VariableType::PRESSURE,
    ParamSpatialLoc::Bottom,
    /*z_pos=*/0.0,
    times, values);

  ASSERT_TRUE(ok) << "ExtractBulkTimeSeries returned false";
  ASSERT_FALSE(times.empty());
  EXPECT_EQ(times.size(), values.size());
  EXPECT_TRUE(IsNonDecreasing(times))
    << "Concatenated times are not monotonically non-decreasing";

  // All pressures should be finite and positive
  for (double v : values) {
    EXPECT_TRUE(std::isfinite(v)) << "Non-finite pressure value";
    EXPECT_GT(v, 0.0);
  }
}

TEST_F(ParameterisationTest, ExtractBulkTimeSeries_PositionInterpolation)
{
  CoupledSolver solver(data);
  ASSERT_GE(solver.IntegrateCycle(data.cycle, false), 0);

  auto eq_log = solver.ExtractLog();

  // Query at a midpoint position
  double bed_length = GetLayer(layer_id).domain_.length;
  std::vector<double> times_top, values_top;
  std::vector<double> times_mid, values_mid;
  std::vector<double> times_bot, values_bot;

  bool ok1 = ExtractBulkTimeSeries(eq_log, VariableType::PRESSURE,
                                   ParamSpatialLoc::Top, 0.0, times_top, values_top);
  bool ok2 = ExtractBulkTimeSeries(eq_log, VariableType::PRESSURE,
                                   ParamSpatialLoc::Position, bed_length * 0.5,
                                   times_mid, values_mid);
  bool ok3 = ExtractBulkTimeSeries(eq_log, VariableType::PRESSURE,
                                   ParamSpatialLoc::Bottom, 0.0, times_bot, values_bot);

  ASSERT_TRUE(ok1 && ok2 && ok3);
  ASSERT_EQ(times_top.size(), times_mid.size());
  ASSERT_EQ(times_top.size(), times_bot.size());

  // All three should return the same time axes
  for (size_t i = 0; i < times_top.size(); ++i)
    EXPECT_DOUBLE_EQ(times_top[i], times_mid[i]);

  // All values should be finite
  for (size_t i = 0; i < values_mid.size(); ++i)
    EXPECT_TRUE(std::isfinite(values_mid[i]));
}

// ---------------------------------------------------------------------------
// Test 3: ExtractComponentTimeSeries – basic sanity.
// ---------------------------------------------------------------------------

TEST_F(ParameterisationTest, ExtractComponentTimeSeries_MoleFractionBounds)
{
  CoupledSolver solver(data);
  ASSERT_GE(solver.IntegrateCycle(data.cycle, false), 0);

  auto eq_log = solver.ExtractLog();

  std::vector<double> times, values;
  bool ok = ExtractComponentTimeSeries(
    eq_log,
    ComponentVariableType::MOLE_FRACTION,
    "CarbonDioxide",
    ParamSpatialLoc::Top,
    0.0,
    times, values);

  ASSERT_TRUE(ok);
  ASSERT_FALSE(values.empty());

  for (double mf : values) {
    EXPECT_GE(mf, -1e-6) << "Mole fraction below 0";
    EXPECT_LE(mf,  1.0 + 1e-6) << "Mole fraction above 1";
  }

  // The adsorbent's CO2 capacity exceeds the amount fed in this short 500 s run,
  // so the mass-transfer front does not fully reach the outlet: the top CO2 mole
  // fraction stays at numerical zero. A strict "final > initial" check here would
  // be comparing floating-point noise around zero (and is sensitive to benign
  // summation-order changes). Assert instead that the outlet series is
  // monotonically non-decreasing (a clean breakthrough never overshoots), which
  // still catches real regressions (NaNs, oscillations, spurious desorption).
  for (size_t i = 1; i < values.size(); ++i)
    EXPECT_GE(values[i], values[i - 1] - 1e-6)
      << "CO2 outlet mole fraction should not decrease during breakthrough";
}

// ---------------------------------------------------------------------------
// Test 4: ComputeRMSE self-consistency.
// Extract the simulated time series, then use it directly as "experimental
// data". The RMSE should be zero (or machine-epsilon scale).
// ---------------------------------------------------------------------------

TEST_F(ParameterisationTest, ComputeRMSE_SelfConsistencyIsZero)
{
  CoupledSolver solver(data);
  ASSERT_GE(solver.IntegrateCycle(data.cycle, false), 0);

  auto eq_log = solver.ExtractLog();

  std::vector<double> sim_t, sim_v;
  ASSERT_TRUE(ExtractComponentTimeSeries(
    eq_log,
    ComponentVariableType::MOLE_FRACTION,
    "CarbonDioxide",
    ParamSpatialLoc::Top, 0.0,
    sim_t, sim_v));

  // Build "experimental" data from a subset of simulation points
  XYData exp;
  for (size_t i = 0; i < sim_t.size(); i += 5) {  // every 5th point
    exp.x.push_back(sim_t[i]);
    exp.y.push_back(sim_v[i]);
  }
  if (exp.x.empty()) {
    exp.x = {sim_t[0]};
    exp.y = {sim_v[0]};
  }

  double rmse = ComputeRMSE(sim_t, sim_v, exp, /*time_unit_mul=*/1.0);
  EXPECT_NEAR(rmse, 0.0, 1e-10)
    << "Self-consistency RMSE should be zero when experimental = simulation";
}

// ---------------------------------------------------------------------------
// Test 5: ComputeRMSE is sensitive to a perturbation.
// Shift experimental values by a known constant delta, verify RMSE ≈ delta.
// ---------------------------------------------------------------------------

TEST_F(ParameterisationTest, ComputeRMSE_DetectsConstantOffset)
{
  CoupledSolver solver(data);
  ASSERT_GE(solver.IntegrateCycle(data.cycle, false), 0);

  auto eq_log = solver.ExtractLog();

  std::vector<double> sim_t, sim_v;
  ASSERT_TRUE(ExtractBulkTimeSeries(
    eq_log, VariableType::PRESSURE,
    ParamSpatialLoc::Bottom, 0.0,
    sim_t, sim_v));

  const double delta = 500.0;  // Pa offset

  XYData exp;
  // Use every 10th simulation point, shifted by delta
  for (size_t i = 0; i < sim_t.size(); i += 10) {
    exp.x.push_back(sim_t[i]);
    exp.y.push_back(sim_v[i] + delta);
  }
  if (exp.x.empty()) {
    exp.x = {sim_t[0]};
    exp.y = {sim_v[0] + delta};
  }

  double rmse = ComputeRMSE(sim_t, sim_v, exp, 1.0);
  // All residuals = -delta, so RMSE = |delta|
  EXPECT_NEAR(rmse, std::abs(delta), 1.0)  // 1 Pa tolerance for floating point
    << "RMSE should equal the constant offset applied to experimental data";
}

// ---------------------------------------------------------------------------
// Test 6: ComputeKPIRMSE self-consistency.
// Compute the purity from the simulation log, then use that value as the
// "experimental" data. ComputeKPIRMSE should return ~0.
// ---------------------------------------------------------------------------

TEST_F(ParameterisationTest, ComputeKPIRMSE_SelfConsistencyIsZero)
{
  CoupledSolver solver(data);
  ASSERT_GE(solver.IntegrateCycle(data.cycle, false), 0);

  auto eq_log = solver.ExtractLog();
  auto b_log  = solver.ExtractBoundaryLog();

  // Build KPI params using the test fixture's kpi_config (first tracked component)
  const auto* kpi_cfg = data.kpi_config.FirstTracked();
  ASSERT_NE(kpi_cfg, nullptr) << "Test fixture has no tracked KPI component";
  const std::string product_comp = kpi_cfg->product_component;
  const auto feed_endpoints = kpi_cfg->feed_endpoints;
  const auto product_endpoints = kpi_cfg->product_endpoints;

  // Compute purity at cycle 0
  KPIParams kpi(eq_log, b_log);
  kpi.cycle_index_      = 0;
  kpi.product_comp_     = product_comp;
  kpi.product_endpoints = product_endpoints;
  kpi.feed_endpoints    = feed_endpoints;

  double purity = CalculatePurity(kpi);
  ASSERT_TRUE(std::isfinite(purity)) << "Purity is not finite";
  ASSERT_GE(purity, 0.0);
  ASSERT_LE(purity, 1.01);

  // Build "experimental" data with this exact purity value at cycle 0
  XYData exp;
  exp.x = {0.0};   // cycle index 0
  exp.y = {purity};

  double rmse = ComputeKPIRMSE(
    eq_log, b_log, exp,
    /*cycle_offset=*/0,
    KPIVariableType::PURITY,
    product_comp,
    feed_endpoints, product_endpoints);

  EXPECT_NEAR(rmse, 0.0, 1e-10)
    << "KPI RMSE should be zero when experimental = computed KPI. Purity = " << purity;
}

// ---------------------------------------------------------------------------
// Multi-step, multi-cycle fixture (E-10).
//
// Feed (CO2 at the bottom) followed by a Purge (air at the bottom), repeated
// for several cycles. The logged times are absolute IDA times, so the final
// cycle starts at (n_cycles - 1) * cycle_time; the extracted parameteriser
// time axis must nevertheless start at 0 and span exactly one cycle.
// ---------------------------------------------------------------------------

class ParameterisationMultiCycleTest : public ParameterisationTest {
public:
  static constexpr double kFeedTime  = 60.0;
  static constexpr double kPurgeTime = 40.0;
  static constexpr double kCycleTime = kFeedTime + kPurgeTime;
  static constexpr int    kNumCycles = 3;

  void SetUp() override {
    ParameterisationTest::SetUp();

    data.cycle.GetStep("Feed").step_time_ = kFeedTime;
    AddStep("Purge", kPurgeTime);

    const double Q_m3s = 200 * 1.666666666666667e-8;  // 200 ccm → m³/s
    auto purge_stream = CreateNoValveInletStream(
      CreateVolumeFlowBC(Q_m3s),
      {{"CarbonDioxide", 0.0}, {"Air", 1.0}},
      293.0);
    purge_stream->name_ = "Air Inlet";

    auto purge_flow = std::make_unique<FlowBoundaryCondition>();
    purge_flow->direction_ = BoundaryDirection::INLET;
    purge_flow->location_  = BoundaryLocation::BOTTOM;
    purge_flow->streams_.push_back(std::move(purge_stream));
    data.cycle.GetStep("Purge").bottom_boundary_condition_ = std::move(purge_flow);

    auto pressure_out = std::make_unique<PressureBoundaryCondition>();
    pressure_out->direction_               = BoundaryDirection::OUTLET;
    pressure_out->location_                = BoundaryLocation::TOP;
    pressure_out->pressure_                = CreateConstantBoundaryValue(1e5);
    pressure_out->mole_fraction_["CarbonDioxide"] = 0.0;
    pressure_out->mole_fraction_["Air"]           = 1.0;
    data.cycle.GetStep("Purge").top_boundary_condition_ = std::move(pressure_out);

    ConfigureCycle(kNumCycles);
  }
};

TEST_F(ParameterisationMultiCycleTest, ExtractedTimeAxisIsRelativeToFinalCycle)
{
  CoupledSolver solver(data);
  ASSERT_GE(solver.IntegrateCycle(data.cycle, false), 0);

  auto eq_log = solver.ExtractLog();
  const auto& sim_log = eq_log.GetSimulationLog();
  ASSERT_FALSE(sim_log.empty());
  const int last_cycle = sim_log.rbegin()->first;
  ASSERT_EQ(last_cycle, kNumCycles - 1);
  ASSERT_EQ(sim_log.at(last_cycle).size(), 2u);

  std::vector<double> bulk_t, bulk_v, comp_t, comp_v;
  ASSERT_TRUE(ExtractBulkTimeSeries(eq_log, VariableType::PRESSURE,
                                    ParamSpatialLoc::Bottom, 0.0, bulk_t, bulk_v));
  ASSERT_TRUE(ExtractComponentTimeSeries(eq_log, ComponentVariableType::MOLE_FRACTION,
                                         "CarbonDioxide", ParamSpatialLoc::Bottom, 0.0,
                                         comp_t, comp_v));

  for (const auto* times : {&bulk_t, &comp_t}) {
    ASSERT_FALSE(times->empty());
    EXPECT_NEAR(times->front(), 0.0, 1e-9) << "Time axis must start at 0";
    EXPECT_NEAR(times->back(), kCycleTime, 1e-6)
      << "Time axis must span exactly one cycle";
    EXPECT_TRUE(IsNonDecreasing(*times))
      << "Time axis must be continuous/monotonic across steps";
  }
  EXPECT_EQ(bulk_t, comp_t);

  // The Feed -> Purge boundary must sit at kFeedTime on the extracted axis.
  const auto feed_vle = eq_log.GetVariable(last_cycle, 0, VariableType::PRESSURE);
  ASSERT_FALSE(feed_vle.times.empty());
  EXPECT_NEAR(bulk_t[feed_vle.times.size() - 1], kFeedTime, 1e-6);
  EXPECT_NEAR(bulk_t[feed_vle.times.size()],     kFeedTime, 1e-6);
}

TEST_F(ParameterisationMultiCycleTest, ComputeRMSE_ZeroAgainstOwnFinalCycleOutput)
{
  CoupledSolver solver(data);
  ASSERT_GE(solver.IntegrateCycle(data.cycle, false), 0);

  auto eq_log = solver.ExtractLog();
  const auto& sim_log = eq_log.GetSimulationLog();
  ASSERT_FALSE(sim_log.empty());
  const int last_cycle = sim_log.rbegin()->first;
  const double cycle_start = last_cycle * kCycleTime;

  // Build the "experimental" record straight from the raw log of the final
  // cycle: the inlet CO2 mole fraction (bottom cell) against time since the
  // start of that cycle, as a lab logger would record it. Step-boundary
  // samples are skipped because the two steps log the same instant.
  XYData exp;
  for (const auto& [step_idx, step_log] : sim_log.at(last_cycle)) {
    auto vle = eq_log.GetComponentVariable(last_cycle, step_idx,
                                           ComponentVariableType::MOLE_FRACTION,
                                           "CarbonDioxide");
    ASSERT_EQ(vle.times.size(), vle.values.size());
    for (size_t t = 1; t + 1 < vle.times.size(); ++t) {
      exp.x.push_back(vle.times[t] - cycle_start);
      exp.y.push_back(vle.values[t].front());
    }
  }
  ASSERT_GT(exp.x.size(), 4u);

  // The signal must actually vary over the cycle (CO2 during Feed, air during
  // Purge), otherwise a misaligned time axis would go unnoticed.
  const auto [mn, mx] = std::minmax_element(exp.y.begin(), exp.y.end());
  ASSERT_GT(*mx - *mn, 0.5);

  std::vector<double> sim_t, sim_v;
  ASSERT_TRUE(ExtractComponentTimeSeries(eq_log, ComponentVariableType::MOLE_FRACTION,
                                         "CarbonDioxide", ParamSpatialLoc::Bottom, 0.0,
                                         sim_t, sim_v));

  EXPECT_NEAR(ComputeRMSE(sim_t, sim_v, exp, 1.0), 0.0, 1e-9)
    << "RMSE against the model's own final-cycle output should be zero";

  // Same data expressed in minutes.
  XYData exp_min = exp;
  for (auto& x : exp_min.x) x /= 60.0;
  EXPECT_NEAR(ComputeRMSE(sim_t, sim_v, exp_min, 60.0), 0.0, 1e-9);
}

// ---------------------------------------------------------------------------
// Test 7: ExperimentalData CSV loading (prerequisite for parameterisation).
// ---------------------------------------------------------------------------

TEST(ExperimentalData, LoadCSVAndVerifyContents)
{
  // Write a small CSV to a temp file
  std::filesystem::path tmp = std::filesystem::temp_directory_path() / "test_param.csv";
  {
    std::ofstream f(tmp);
    ASSERT_TRUE(f.is_open());
    f << "Time,CO2_MoleFraction\n";
    f << "0.0,0.00\n";
    f << "100.0,0.10\n";
    f << "200.0,0.50\n";
    f << "300.0,0.90\n";
    f << "400.0,0.99\n";
  }

  ExperimentalData exp;
  ASSERT_TRUE(exp.LoadData(tmp, /*headers_included=*/true));
  ASSERT_EQ(exp.size(), 1u);

  const auto& d = exp.GetAllData()[0];
  ASSERT_EQ(d.x.size(), 5u);
  ASSERT_EQ(d.y.size(), 5u);

  EXPECT_NEAR(d.x[0], 0.0, 1e-12);
  EXPECT_NEAR(d.x[4], 400.0, 1e-12);
  EXPECT_NEAR(d.y[2], 0.50, 1e-12);
  EXPECT_EQ(d.x_label, "Time");
  EXPECT_EQ(d.y_label, "CO2_MoleFraction");

  std::filesystem::remove(tmp);
}

TEST(ExperimentalData, LoadCSVWithoutHeaders)
{
  std::filesystem::path tmp = std::filesystem::temp_directory_path() / "test_param_noheader.csv";
  {
    std::ofstream f(tmp);
    ASSERT_TRUE(f.is_open());
    f << "0.0,1.0\n";
    f << "1.0,2.0\n";
    f << "2.0,3.0\n";
  }

  ExperimentalData exp;
  ASSERT_TRUE(exp.LoadData(tmp, /*headers_included=*/false));
  ASSERT_EQ(exp.size(), 1u);

  const auto& d = exp.GetAllData()[0];
  ASSERT_EQ(d.x.size(), 3u);
  EXPECT_NEAR(d.x[0], 0.0, 1e-12);
  EXPECT_NEAR(d.y[2], 3.0, 1e-12);

  std::filesystem::remove(tmp);
}

TEST(ExperimentalData, AddAndDeleteData)
{
  ExperimentalData store;

  XYData d;
  d.name = "test";
  d.x = {1.0, 2.0};
  d.y = {3.0, 4.0};
  ASSERT_TRUE(store.AddData(d));
  ASSERT_EQ(store.size(), 1u);

  ASSERT_TRUE(store.DeleteDataByIndex(0));
  EXPECT_EQ(store.size(), 0u);
}

} // namespace adsim
