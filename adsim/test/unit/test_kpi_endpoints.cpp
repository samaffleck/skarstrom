// test_kpi_endpoints.cpp
//
// Unit tests for the endpoint-keyed KPI maths (kpis.h). These build a BoundaryLog
// by hand (no full solve) so we can assert the exact selection behaviour:
//   * a flow boundary with TWO outlet streams — selecting only the product stream
//     must ignore the other (recycle) stream on the same face;
//   * folding a re-entry (recycle) endpoint into the product list subtracts it
//     via the signed sum (the replacement for the old "product - recycled" term);
//   * the velocity-sign -> entering/leaving conversion keyed off the endpoint's
//     bed end (Top negated, Bottom not).

#include <gtest/gtest.h>

#include <cmath>
#include <set>
#include <vector>

#include "adsim/solver/kpis.h"
#include "adsim/solver/log.h"
#include "adsim/solver/cycle.h"   // GetAllBoundaryStreams
#include "adsim/boundary_conditions/step.h"
#include "adsim/boundary_conditions/boundary_conditions.h"
#include "adsim/boundary_conditions/boundary_stream.h"
#include "adsim/boundary_conditions/boundary_value.h"

namespace adsim {
namespace {

  const std::string CO2 = "CarbonDioxide";

  // Make a flow stream whose logged StreamData is a constant molar flow over t=[0,1].
  // `total`/`co2` are the constant [mol/s] rates (velocity-signed: + = bottom->top).
  std::unique_ptr<BoundaryStream> MakeStream(const std::string& name, double total, double co2) {
    auto s = std::make_unique<BoundaryStream>();
    s->name_ = name;
    s->stream_data_.time_series_ = { 0.0, 1.0 };
    s->stream_data_.molar_flow_ = { total, total };
    s->stream_data_.component_molar_flow_[CO2] = { co2, co2 };
    return s;
  }

  std::unique_ptr<FlowBoundaryCondition> MakeFlowBC(
      BoundaryLocation loc, std::vector<std::unique_ptr<BoundaryStream>> streams) {
    auto bc = std::make_unique<FlowBoundaryCondition>();
    bc->location_ = loc;
    for (auto& s : streams) bc->streams_.push_back(std::move(s));
    return bc;
  }

  // A one-step cycle:
  //   TOP flow BC  : "Product" (leaving, +2 total / +1.6 CO2) and
  //                  "Recycle" (re-entering, -1 total / -0.5 CO2)
  //   BOTTOM flow BC: "Feed"   (entering, +4 total / +2 CO2)
  // Returns the Step (keeps ownership of the BCs) and a populated BoundaryLog.
  struct Rig {
    Step step;
    EquationLog eq_log;
    BoundaryLog b_log;
  };

  std::shared_ptr<Rig> BuildRig() {
    auto rig = std::make_shared<Rig>();

    std::vector<std::unique_ptr<BoundaryStream>> top;
    top.push_back(MakeStream("Product",  2.0,  1.6));   // outlet at top: velocity up
    top.push_back(MakeStream("Recycle", -1.0, -0.5));   // re-entry at top: velocity down
    rig->step.top_boundary_condition_ = MakeFlowBC(BoundaryLocation::TOP, std::move(top));

    std::vector<std::unique_ptr<BoundaryStream>> bottom;
    bottom.push_back(MakeStream("Feed", 4.0, 2.0));     // inlet at bottom: velocity up
    rig->step.bottom_boundary_condition_ = MakeFlowBC(BoundaryLocation::BOTTOM, std::move(bottom));

    rig->eq_log.sys_data.molar_mass[CO2] = 0.044;
    rig->b_log.sys_data.step_list_.push_back("Step");
    rig->b_log.LogEndOfStep(rig->step, /*cycle=*/0, /*step_idx=*/0);
    return rig;
  }

  KPIEndpointRef Ep(int step, const std::string& name, BoundaryLocation loc) {
    return KPIEndpointRef{ step, name, loc };
  }

} // namespace

// Selecting only the "Product" stream on a two-stream top face ignores "Recycle".
TEST(KpiEndpoints, ProductStreamSelectionIgnoresOtherStreamOnSameFace) {
  auto rig = BuildRig();
  KPIParams p(rig->eq_log, rig->b_log);
  p.cycle_index_ = 0;
  p.product_comp_ = CO2;
  p.feed_endpoints    = { Ep(0, "Feed", BoundaryLocation::BOTTOM) };
  p.product_endpoints = { Ep(0, "Product", BoundaryLocation::TOP) };

  ASSERT_TRUE(p.Check());
  // Product CO2 = 1.6 (Product stream only); total product = 2.0.
  EXPECT_NEAR(CalculateComponentMolesInProduct(p), 1.6, 1e-9);
  EXPECT_NEAR(CalculateTotalMolesInProduct(p),     2.0, 1e-9);
  // Feed CO2 = 2.0. Recovery = 1.6/2.0; purity = 1.6/2.0.
  EXPECT_NEAR(CalculateComponentMolesInFeed(p), 2.0, 1e-9);
  EXPECT_NEAR(CalculateRecovery(p), 0.8, 1e-9);
  EXPECT_NEAR(CalculatePurity(p),   0.8, 1e-9);
}

// Folding the re-entry "Recycle" endpoint into the product list subtracts it
// (signed sum) — reproducing the old "product - recycled product" behaviour.
TEST(KpiEndpoints, RecycleFoldedIntoProductSubtracts) {
  auto rig = BuildRig();
  KPIParams p(rig->eq_log, rig->b_log);
  p.cycle_index_ = 0;
  p.product_comp_ = CO2;
  p.feed_endpoints    = { Ep(0, "Feed", BoundaryLocation::BOTTOM) };
  p.product_endpoints = {
    Ep(0, "Product", BoundaryLocation::TOP),
    Ep(0, "Recycle", BoundaryLocation::TOP),
  };

  ASSERT_TRUE(p.Check());
  // Net product CO2 = 1.6 - 0.5 = 1.1; net total product = 2.0 - 1.0 = 1.0.
  EXPECT_NEAR(CalculateComponentMolesInProduct(p), 1.1, 1e-9);
  EXPECT_NEAR(CalculateTotalMolesInProduct(p),     1.0, 1e-9);
  EXPECT_NEAR(CalculateComponentMolesInFeed(p),    2.0, 1e-9);
  EXPECT_NEAR(CalculateRecovery(p), 1.1 / 2.0, 1e-9);
}

// Contract guard: every name from Step::GetEndpoints() (what the KPI config
// stores) must match a BoundaryLog key from GetAllBoundaryStreams (what the KPI
// compute looks up) — for flow streams AND pressure boundaries. A mismatch would
// silently zero that endpoint's KPI contribution.
TEST(KpiEndpoints, EndpointNamesMatchBoundaryLogKeys) {
  Step step;                       // default: pressure BC @ Top, flow "Stream 1" @ Bottom
  step.step_name_ = "Evacuation";

  std::set<std::string> endpoint_names, log_keys;
  for (const auto& e : step.GetEndpoints()) endpoint_names.insert(e.name);
  for (const auto& b : GetAllBoundaryStreams(step)) log_keys.insert(b.boundary_name);

  EXPECT_EQ(endpoint_names, log_keys);
  // The pressure endpoint must use the log-key form, not the bare canonical name.
  EXPECT_TRUE(endpoint_names.count("Evacuation Pressure Boundary Top") == 1);
}

// Check() fails when no product or no feed endpoint is selected.
TEST(KpiEndpoints, CheckRequiresFeedAndProduct) {
  auto rig = BuildRig();
  KPIParams p(rig->eq_log, rig->b_log);
  p.cycle_index_ = 0;
  p.product_comp_ = CO2;
  p.product_endpoints = { Ep(0, "Product", BoundaryLocation::TOP) };
  EXPECT_FALSE(p.Check());  // no feed endpoint
  p.feed_endpoints = { Ep(0, "Feed", BoundaryLocation::BOTTOM) };
  EXPECT_TRUE(p.Check());
}

// L-04: every ratio KPI with a zero denominator (no product, no feed, zero
// adsorbent mass / bed volume / cycle time) is "undefined": it returns exactly
// 0 — never NaN/Inf — and reports defined=false. A healthy KPI reports true.
namespace {
  // Every KPI with its definedness flag.
  struct KpiOut { double value; bool defined; };
  std::vector<std::pair<std::string, KpiOut>> AllKpis(const KPIParams& p) {
    std::vector<std::pair<std::string, KpiOut>> out;
    auto add = [&](const char* name, double (*fn)(const KPIParams&, bool*)) {
      bool d = true;
      double v = fn(p, &d);
      out.push_back({ name, { v, d } });
    };
    add("purity", &CalculatePurity);
    add("recovery", &CalculateRecovery);
    add("productivity", &CalculateProductivity);
    add("productivity_mass", &CalculateProductivityMass);
    add("productivity_volume", &CalculateProductivityVolume);
    add("production_rate", &CalculateProductionRate);
    add("specific_energy", &CalculateSpecificEnergyConsumption);
    return out;
  }
} // namespace

TEST(KpiEndpoints, ZeroDenominatorsReturnZeroAndUndefined) {
  // Product and feed endpoints both carry zero flow; zero adsorbent mass, zero
  // cycle time and an empty domain (zero bed volume).
  Step step;
  std::vector<std::unique_ptr<BoundaryStream>> top;
  top.push_back(MakeStream("Product", 0.0, 0.0));
  step.top_boundary_condition_ = MakeFlowBC(BoundaryLocation::TOP, std::move(top));
  std::vector<std::unique_ptr<BoundaryStream>> bottom;
  bottom.push_back(MakeStream("Feed", 0.0, 0.0));
  step.bottom_boundary_condition_ = MakeFlowBC(BoundaryLocation::BOTTOM, std::move(bottom));

  EquationLog eq_log;
  BoundaryLog b_log;
  eq_log.sys_data.molar_mass[CO2] = 0.044;
  eq_log.sys_data.adsorbent_mass = 0.0;
  eq_log.sys_data.total_cycle_time = 0.0;
  b_log.sys_data.step_list_.push_back("Step");
  b_log.LogEndOfStep(step, /*cycle=*/0, /*step_idx=*/0);

  KPIParams p(eq_log, b_log);
  p.cycle_index_ = 0;
  p.product_comp_ = CO2;
  p.feed_endpoints    = { Ep(0, "Feed", BoundaryLocation::BOTTOM) };
  p.product_endpoints = { Ep(0, "Product", BoundaryLocation::TOP) };
  ASSERT_TRUE(p.Check());

  for (const auto& [name, k] : AllKpis(p)) {
    EXPECT_EQ(k.value, 0.0) << name;
    EXPECT_FALSE(k.defined) << name;
  }
  // The default (no out-parameter) overloads give the same finite 0.
  EXPECT_EQ(CalculatePurity(p), 0.0);
  EXPECT_EQ(CalculateRecovery(p), 0.0);
  EXPECT_EQ(CalculateSpecificEnergyConsumption(p), 0.0);
}

// Product flows but the adsorbent mass / bed volume are zero: only the per-mass
// and per-volume productivities are undefined; the rest stay defined.
TEST(KpiEndpoints, ZeroAdsorbentMassOnlyUndefinesItsKpi) {
  auto rig = BuildRig();
  rig->eq_log.sys_data.adsorbent_mass = 0.0;
  rig->eq_log.sys_data.total_cycle_time = 10.0;
  KPIParams p(rig->eq_log, rig->b_log);
  p.cycle_index_ = 0;
  p.product_comp_ = CO2;
  p.feed_endpoints    = { Ep(0, "Feed", BoundaryLocation::BOTTOM) };
  p.product_endpoints = { Ep(0, "Product", BoundaryLocation::TOP) };
  ASSERT_TRUE(p.Check());

  for (const auto& [name, k] : AllKpis(p)) {
    EXPECT_TRUE(std::isfinite(k.value)) << name;
    const bool expect_defined = name != "productivity_mass" && name != "productivity_volume";
    EXPECT_EQ(k.defined, expect_defined) << name;
    if (!expect_defined) EXPECT_EQ(k.value, 0.0) << name;
  }
  bool d = false;
  EXPECT_NEAR(CalculateProductionRate(p, &d), 1.6 / 10.0, 1e-12);
  EXPECT_TRUE(d);
  EXPECT_NEAR(CalculateProductivity(p, &d), 1.6 * 0.044 / 10.0, 1e-12);
  EXPECT_TRUE(d);
  // No compressor work, product drawn: specific energy is a defined 0.
  EXPECT_EQ(CalculateSpecificEnergyConsumption(p, &d), 0.0);
  EXPECT_TRUE(d);

  rig->eq_log.sys_data.adsorbent_mass = 2.0;
  EXPECT_NEAR(CalculateProductivityMass(p, &d), 1.6 * 0.044 / 10.0 / 2.0, 1e-12);
  EXPECT_TRUE(d);
}

// A failed Check() is undefined too (0, defined=false).
TEST(KpiEndpoints, InvalidParamsAreUndefined) {
  auto rig = BuildRig();
  rig->eq_log.sys_data.adsorbent_mass = 1.0;
  rig->eq_log.sys_data.total_cycle_time = 1.0;
  KPIParams p(rig->eq_log, rig->b_log);
  p.cycle_index_ = 0;
  p.product_comp_ = CO2;
  p.product_endpoints = { Ep(0, "Product", BoundaryLocation::TOP) };  // no feed
  ASSERT_FALSE(p.Check());
  for (const auto& [name, k] : AllKpis(p)) {
    EXPECT_EQ(k.value, 0.0) << name;
    EXPECT_FALSE(k.defined) << name;
  }
}

} // namespace adsim
