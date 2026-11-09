// test_kpi_uncounted_inflow.cpp
//
// The uncounted-inflow diagnostic (KPIRecord::uncounted_inflow): the gross
// inflow of the tracked component through unconnected endpoints that are
// neither feed nor recycle. Built from a hand-made BoundaryLog (no solve):
//   Feed      bottom, in 2.0 mol CO2          (feed role: counted)
//   Product   top, out 1.6 mol CO2            (product role)
//   Vent      top, unconnected, backflows 0.3 mol CO2   (counted by the diagnostic)
//   Equalise  top, connected, in 0.7 mol CO2  (internal gas: excluded)

#include <gtest/gtest.h>

#include <cmath>
#include <memory>
#include <vector>

#include "adsim/boundary_conditions/boundary_conditions.h"
#include "adsim/boundary_conditions/boundary_stream.h"
#include "adsim/boundary_conditions/step.h"
#include "adsim/solver/kpi_config.h"
#include "adsim/solver/kpi_history.h"
#include "adsim/solver/kpis.h"
#include "adsim/solver/log.h"

namespace adsim {
namespace {

  const std::string CO2 = "CarbonDioxide";

  std::unique_ptr<BoundaryStream> MakeStream(const std::string& name, double total, double co2) {
    auto s = std::make_unique<BoundaryStream>();
    s->name_ = name;
    s->stream_data_.time_series_ = { 0.0, 1.0 };
    s->stream_data_.molar_flow_ = { total, total };
    s->stream_data_.component_molar_flow_[CO2] = { co2, co2 };
    return s;
  }

  struct Rig {
    Step step;
    BoundaryStream partner;  // what "Equalise" is connected to
    EquationLog eq_log;
    BoundaryLog b_log;
    KPIConfig config;
    std::vector<Step> steps;
  };

  // `vent_quadrature`: give Vent quadrature integrals (in 0.25, out 0.1 mol CO2).
  std::shared_ptr<Rig> BuildRig(bool vent_quadrature) {
    auto rig = std::make_shared<Rig>();

    auto top = std::make_unique<FlowBoundaryCondition>();
    top->location_ = BoundaryLocation::TOP;
    top->streams_.push_back(MakeStream("Product", 2.0, 1.6));    // leaving
    top->streams_.push_back(MakeStream("Vent", -0.5, -0.3));     // entering (backflow)
    top->streams_.push_back(MakeStream("Equalise", -1.0, -0.7)); // entering
    top->streams_.back()->connection_.MakeConnection(7, rig->partner);
    if (vent_quadrature) {
      EndpointIntegrals& q = top->streams_[1]->stream_data_.integrals_;
      q.valid = true;
      q.component_in[CO2] = 0.25;
      q.component_out[CO2] = 0.1;
    }
    rig->step.top_boundary_condition_ = std::move(top);

    auto bottom = std::make_unique<FlowBoundaryCondition>();
    bottom->location_ = BoundaryLocation::BOTTOM;
    bottom->streams_.push_back(MakeStream("Feed", 4.0, 2.0));
    rig->step.bottom_boundary_condition_ = std::move(bottom);

    rig->eq_log.sys_data.molar_mass[CO2] = 0.044;
    rig->eq_log.sys_data.total_cycle_time = 10.0;
    rig->b_log.sys_data.step_list_.push_back("Step");
    rig->b_log.LogEndOfStep(rig->step, 0, 0);

    rig->config.SyncToComponents({ CO2 });
    ComponentKPIConfig* c = rig->config.FindByComponent(CO2);
    c->tracked = true;
    c->SetFeed(0, "Feed", BoundaryLocation::BOTTOM, true);
    c->SetProduct(0, "Product", BoundaryLocation::TOP, true);
    rig->steps.push_back(rig->step);
    return rig;
  }

  const KPIRecord& Record(const std::vector<KPIRecord>& v) {
    EXPECT_EQ(v.size(), 1u);
    return v.front();
  }

}  // namespace

// Trapezoid fallback: the positive part of Vent's flow counts, the connected
// Equalise endpoint does not, and no other KPI moves.
TEST(KpiUncountedInflow, VentCountedConnectedEndpointExcluded) {
  auto rig = BuildRig(false);
  const auto recs = ComputeCycleKPIs(rig->config, rig->steps, rig->eq_log, rig->b_log, 0);
  const KPIRecord& r = Record(recs);
  ASSERT_TRUE(r.valid);
  EXPECT_NEAR(r.uncounted_inflow, 0.3, 1e-12);
  EXPECT_NEAR(r.uncounted_inflow_fraction, 0.3 / 2.0, 1e-12);
  // The existing KPIs ignore the diagnostic: recovery is net product over net feed.
  EXPECT_NEAR(r.recovery, 1.6 / 2.0, 1e-12);
  EXPECT_NEAR(r.purity, 1.6 / 2.0, 1e-12);
}

// With quadratures the gross inflow is read, not the net.
TEST(KpiUncountedInflow, QuadratureUsesGrossInflow) {
  auto rig = BuildRig(true);
  const KPIRecord& r = Record(ComputeCycleKPIs(rig->config, rig->steps, rig->eq_log, rig->b_log, 0));
  EXPECT_NEAR(r.uncounted_inflow, 0.25, 1e-12);
}

// An endpoint configured as feed or recycle is not uncounted.
TEST(KpiUncountedInflow, FeedAndRecycleRolesAreCounted) {
  auto rig = BuildRig(false);
  ComponentKPIConfig* c = rig->config.FindByComponent(CO2);
  c->SetFeed(0, "Vent", BoundaryLocation::TOP, true);
  EXPECT_NEAR(Record(ComputeCycleKPIs(rig->config, rig->steps, rig->eq_log, rig->b_log, 0)).uncounted_inflow, 0.0, 1e-12);
  c->SetFeed(0, "Vent", BoundaryLocation::TOP, false);
  c->SetRecycle(0, "Vent", BoundaryLocation::TOP, true);
  EXPECT_NEAR(Record(ComputeCycleKPIs(rig->config, rig->steps, rig->eq_log, rig->b_log, 0)).uncounted_inflow, 0.0, 1e-12);
}

// No external inflow: exactly zero, and the fraction is zero too.
TEST(KpiUncountedInflow, ZeroWithoutExternalInflow) {
  auto rig = BuildRig(false);
  // A top end without the Vent stream (the log still holds the other streams).
  auto top = std::make_unique<FlowBoundaryCondition>();
  top->location_ = BoundaryLocation::TOP;
  top->streams_.push_back(MakeStream("Product", 2.0, 1.6));
  rig->steps[0].top_boundary_condition_ = std::move(top);
  const KPIRecord& r = Record(ComputeCycleKPIs(rig->config, rig->steps, rig->eq_log, rig->b_log, 0));
  EXPECT_EQ(r.uncounted_inflow, 0.0);
  EXPECT_EQ(r.uncounted_inflow_fraction, 0.0);
}

}  // namespace adsim
