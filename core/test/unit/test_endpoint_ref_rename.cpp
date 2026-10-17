// test_endpoint_ref_rename.cpp
//
// Pins endpoint references across the edits E-20 covers (bug_review.md):
//   - UserData::SwapStepBoundaries ("Swap Boundary Conditions"): the step's
//     KPI refs follow their endpoints to the other bed end;
//   - UserData::RenameStep: buffer-tank sources naming the step follow it;
//   - UserData::RenameStream: KPI refs and buffer-tank sources naming the
//     stream follow it.
// Before E-20 the desktop swap was a bare std::swap, and neither rename
// touched buffer-tank sources (a stream rename did not touch KPI refs either),
// so the references silently dangled (KPIConfig::Prune dropped them; the
// buffer tank failed Cycle::PerformChecks or, in the GUI, was re-pointed at
// another step's endpoint).

#include <gtest/gtest.h>

#include "adsim/objects/user_data.h"
#include "adsim/solver/kpi_config.h"
#include "adsim/boundary_conditions/step.h"
#include "adsim/boundary_conditions/boundary_conditions.h"
#include "adsim/boundary_conditions/boundary_stream.h"
#include "adsim/boundary_conditions/buffer_tank_source.h"

namespace adsim {
namespace {

  const std::string CO2 = "CarbonDioxide";

  std::unique_ptr<FlowBoundaryCondition> Flow(BoundaryLocation loc,
                                              std::initializer_list<const char*> names) {
    auto flow = std::make_unique<FlowBoundaryCondition>();
    flow->location_ = loc;
    for (const char* n : names) {
      auto stream = std::make_unique<BoundaryStream>();
      stream->name_ = n;
      flow->streams_.push_back(std::move(stream));
    }
    return flow;
  }

  std::unique_ptr<PressureBoundaryCondition> Pressure(BoundaryLocation loc) {
    auto p = std::make_unique<PressureBoundaryCondition>();
    p->location_ = loc;
    return p;
  }

  IBoundaryStream& StreamAt(UserData& ud, size_t step, BoundaryLocation side, size_t i) {
    Step& s = ud.cycle.step_list_.at(step);
    auto& slot = side == BoundaryLocation::TOP ? s.top_boundary_condition_ : s.bottom_boundary_condition_;
    return *dynamic_cast<FlowBoundaryCondition&>(*slot).streams_.at(i);
  }

  // Steps:
  //   0 "Pressurise": top flow ["Stream 1"],          bottom pressure
  //   1 "Feed":       top flow ["Stream 1"],          bottom pressure
  //   2 "Purge":      top flow ["Purge In"] (inlet),  bottom flow ["Stream 1"]
  // KPI (CO2): feed = Feed's pressure boundary (bottom); product = Feed's
  // "Stream 1" (top) and Pressurise's "Stream 1" (top).
  // Buffer tanks on Purge's "Purge In": sourced from Feed / "Stream 1".
  // Buffer tank on Pressurise's bottom pressure BC: sourced from Feed's
  // pressure boundary.
  UserData BuildDeck() {
    UserData ud;
    for (const char* n : { "Pressurise", "Feed", "Purge" }) ud.cycle.AddStep(n);
    auto& steps = ud.cycle.step_list_;
    steps[0].top_boundary_condition_ = Flow(BoundaryLocation::TOP, { "Stream 1" });
    steps[0].bottom_boundary_condition_ = Pressure(BoundaryLocation::BOTTOM);
    steps[1].top_boundary_condition_ = Flow(BoundaryLocation::TOP, { "Stream 1" });
    steps[1].bottom_boundary_condition_ = Pressure(BoundaryLocation::BOTTOM);
    steps[2].top_boundary_condition_ = Flow(BoundaryLocation::TOP, { "Purge In" });
    steps[2].bottom_boundary_condition_ = Flow(BoundaryLocation::BOTTOM, { "Stream 1" });

    StreamAt(ud, 2, BoundaryLocation::TOP, 0).buffer_tank_.Activate("Feed", "Stream 1");
    dynamic_cast<PressureBoundaryCondition&>(*steps[0].bottom_boundary_condition_)
      .buffer_tank_.Activate("Feed", kPressureBoundaryEndpointName);

    ud.kpi_config.SyncToComponents({ CO2 });
    ComponentKPIConfig* c = ud.kpi_config.FindByComponent(CO2);
    c->tracked = true;
    c->SetFeed(1, "Feed Pressure Boundary Bottom", BoundaryLocation::BOTTOM, true);
    c->SetProduct(1, "Stream 1", BoundaryLocation::TOP, true);
    c->SetProduct(0, "Stream 1", BoundaryLocation::TOP, true);
    ud.kpi_config.Prune(ud.cycle.step_list_);
    return ud;
  }

  const ComponentKPIConfig& Cfg(const UserData& ud) {
    return *ud.kpi_config.FindByComponent(CO2);
  }

  const BufferTankSource& PurgeTank(UserData& ud) {
    return StreamAt(ud, 2, BoundaryLocation::TOP, 0).buffer_tank_;
  }

  const BufferTankSource& PressuriseTank(UserData& ud) {
    return dynamic_cast<PressureBoundaryCondition&>(*ud.cycle.step_list_[0].bottom_boundary_condition_)
      .buffer_tank_;
  }

  // Every buffer-tank source resolves (what Cycle::PerformChecks requires).
  bool TanksResolve(UserData& ud) {
    ud.cycle.ResolveBufferTankSources();
    return PurgeTank(ud).summary_ != nullptr && PressuriseTank(ud).summary_ != nullptr;
  }

  size_t RefCount(UserData& ud) {
    ud.kpi_config.Prune(ud.cycle.step_list_);
    return Cfg(ud).feed_endpoints.size() + Cfg(ud).product_endpoints.size();
  }

} // namespace

TEST(EndpointRefRename, BaselineResolves) {
  UserData ud = BuildDeck();
  EXPECT_EQ(RefCount(ud), 3u);
  EXPECT_TRUE(TanksResolve(ud));
}

// ---- Swap Boundary Conditions ------------------------------------------------

TEST(EndpointRefRename, SwapBoundariesRemapsKpiRefs) {
  UserData ud = BuildDeck();
  ud.SwapStepBoundaries(1);

  Step& feed = ud.cycle.step_list_[1];
  ASSERT_NE(dynamic_cast<PressureBoundaryCondition*>(feed.top_boundary_condition_.get()), nullptr);
  EXPECT_EQ(feed.top_boundary_condition_->location_, BoundaryLocation::TOP);
  EXPECT_EQ(feed.bottom_boundary_condition_->location_, BoundaryLocation::BOTTOM);

  // Before the Prune (the remap itself, not Prune's location refresh).
  ASSERT_EQ(Cfg(ud).feed_endpoints.size(), 1u);
  EXPECT_EQ(Cfg(ud).feed_endpoints[0].endpoint_name, "Feed Pressure Boundary Top");
  EXPECT_EQ(Cfg(ud).feed_endpoints[0].location, BoundaryLocation::TOP);
  EXPECT_TRUE(Cfg(ud).IsProduct(1, "Stream 1"));
  for (const auto& r : Cfg(ud).product_endpoints) {
    // Feed's stream moved to the bottom; Pressurise's did not move.
    EXPECT_EQ(r.location, r.step_index == 1 ? BoundaryLocation::BOTTOM : BoundaryLocation::TOP);
  }

  EXPECT_EQ(RefCount(ud), 3u);  // all survive Prune
  EXPECT_EQ(Cfg(ud).feed_endpoints[0].location, BoundaryLocation::TOP);
  EXPECT_TRUE(TanksResolve(ud));  // tank names are side-less

  // The same-location check in PerformChecks must not trip on stale location_.
  EXPECT_NE(feed.top_boundary_condition_->location_, feed.bottom_boundary_condition_->location_);
}

TEST(EndpointRefRename, SwapBoundariesTwiceRestores) {
  UserData ud = BuildDeck();
  const auto before_feed = Cfg(ud).feed_endpoints;
  ud.SwapStepBoundaries(1);
  ud.SwapStepBoundaries(1);
  ASSERT_EQ(Cfg(ud).feed_endpoints.size(), 1u);
  EXPECT_EQ(Cfg(ud).feed_endpoints[0].endpoint_name, before_feed[0].endpoint_name);
  EXPECT_EQ(Cfg(ud).feed_endpoints[0].location, before_feed[0].location);
  for (const auto& r : Cfg(ud).product_endpoints) EXPECT_EQ(r.location, BoundaryLocation::TOP);
}

TEST(EndpointRefRename, SwapBoundariesOutOfRangeIsNoOp) {
  UserData ud = BuildDeck();
  ud.SwapStepBoundaries(99);
  EXPECT_EQ(Cfg(ud).feed_endpoints[0].endpoint_name, "Feed Pressure Boundary Bottom");
  EXPECT_EQ(RefCount(ud), 3u);
}

// ---- Step rename ---------------------------------------------------------------

TEST(EndpointRefRename, RenameStepRewritesBufferTankSources) {
  UserData ud = BuildDeck();
  ud.RenameStep(1, "Adsorb");

  EXPECT_EQ(PurgeTank(ud).source_step_name_, "Adsorb");
  EXPECT_EQ(PurgeTank(ud).source_endpoint_name_, "Stream 1");
  EXPECT_EQ(PressuriseTank(ud).source_step_name_, "Adsorb");
  EXPECT_EQ(PressuriseTank(ud).source_endpoint_name_, kPressureBoundaryEndpointName);
  EXPECT_TRUE(TanksResolve(ud));
  EXPECT_EQ(RefCount(ud), 3u);  // E-11 pressure-endpoint rename still holds

  // Renaming a step no tank names leaves the tanks alone.
  ud.RenameStep(2, "Purge 2");
  EXPECT_EQ(PurgeTank(ud).source_step_name_, "Adsorb");
  EXPECT_TRUE(TanksResolve(ud));
}

// ---- Stream rename -------------------------------------------------------------

TEST(EndpointRefRename, RenameStreamRewritesKpiRefsAndBufferTanks) {
  UserData ud = BuildDeck();
  ud.RenameStream(1, BoundaryLocation::TOP, 0, "Product");

  EXPECT_EQ(StreamAt(ud, 1, BoundaryLocation::TOP, 0).name_, "Product");
  EXPECT_TRUE(Cfg(ud).IsProduct(1, "Product"));
  EXPECT_FALSE(Cfg(ud).IsProduct(1, "Stream 1"));
  EXPECT_TRUE(Cfg(ud).IsProduct(0, "Stream 1"));  // Pressurise's same-named stream untouched
  EXPECT_EQ(PurgeTank(ud).source_endpoint_name_, "Product");
  EXPECT_EQ(PressuriseTank(ud).source_endpoint_name_, kPressureBoundaryEndpointName);
  EXPECT_EQ(RefCount(ud), 3u);
  EXPECT_TRUE(TanksResolve(ud));
}

// Renaming a same-named stream on a step no ref names changes nothing else.
TEST(EndpointRefRename, RenameStreamOnOtherStepLeavesRefs) {
  UserData ud = BuildDeck();
  ud.RenameStream(2, BoundaryLocation::BOTTOM, 0, "Waste");
  EXPECT_TRUE(Cfg(ud).IsProduct(1, "Stream 1"));
  EXPECT_TRUE(Cfg(ud).IsProduct(0, "Stream 1"));
  EXPECT_EQ(PurgeTank(ud).source_endpoint_name_, "Stream 1");
  EXPECT_EQ(RefCount(ud), 3u);
  EXPECT_TRUE(TanksResolve(ud));
}

// Both ends of a step carry "Stream 1": only the refs at the renamed stream's
// end follow it; a buffer tank follows only the stream it resolves to (the
// first match, top before bottom).
TEST(EndpointRefRename, RenameStreamWithSameNameAtOtherEnd) {
  UserData ud = BuildDeck();
  Step& feed = ud.cycle.step_list_[1];
  feed.bottom_boundary_condition_ = Flow(BoundaryLocation::BOTTOM, { "Stream 1" });
  ComponentKPIConfig* c = ud.kpi_config.FindByComponent(CO2);
  c->feed_endpoints = { { 1, "Stream 1", BoundaryLocation::BOTTOM } };

  ud.RenameStream(1, BoundaryLocation::BOTTOM, 0, "Feed In");
  ASSERT_EQ(Cfg(ud).feed_endpoints.size(), 1u);
  EXPECT_EQ(Cfg(ud).feed_endpoints[0].endpoint_name, "Feed In");
  EXPECT_TRUE(Cfg(ud).IsProduct(1, "Stream 1"));          // top ref stays
  EXPECT_EQ(PurgeTank(ud).source_endpoint_name_, "Stream 1");  // resolved to top

  ud.RenameStream(1, BoundaryLocation::TOP, 0, "Product");
  EXPECT_TRUE(Cfg(ud).IsProduct(1, "Product"));
  EXPECT_EQ(Cfg(ud).feed_endpoints[0].endpoint_name, "Feed In");
  EXPECT_EQ(PurgeTank(ud).source_endpoint_name_, "Product");
  EXPECT_EQ(RefCount(ud), 3u);
}

// Two streams at the same end share a name: refs to it cannot be told apart,
// so renaming the second one leaves them on the first (which keeps the name);
// once the name is unique again the refs follow it.
TEST(EndpointRefRename, RenameStreamWithSameNameAtSameEnd) {
  UserData ud = BuildDeck();
  ud.cycle.step_list_[1].top_boundary_condition_ =
    Flow(BoundaryLocation::TOP, { "Stream 1", "Stream 1" });

  ud.RenameStream(1, BoundaryLocation::TOP, 1, "Other");
  EXPECT_TRUE(Cfg(ud).IsProduct(1, "Stream 1"));
  EXPECT_FALSE(Cfg(ud).IsProduct(1, "Other"));
  EXPECT_EQ(PurgeTank(ud).source_endpoint_name_, "Stream 1");

  ud.RenameStream(1, BoundaryLocation::TOP, 0, "Product");
  EXPECT_TRUE(Cfg(ud).IsProduct(1, "Product"));
  EXPECT_EQ(PurgeTank(ud).source_endpoint_name_, "Product");
  EXPECT_EQ(RefCount(ud), 3u);
  EXPECT_TRUE(TanksResolve(ud));
}

TEST(EndpointRefRename, RenameStreamBadAddressIsNoOp) {
  UserData ud = BuildDeck();
  ud.RenameStream(9, BoundaryLocation::TOP, 0, "X");       // no such step
  ud.RenameStream(1, BoundaryLocation::BOTTOM, 0, "X");    // pressure side
  ud.RenameStream(1, BoundaryLocation::TOP, 5, "X");       // no such stream
  ud.RenameStream(1, BoundaryLocation::TOP, 0, "Stream 1");  // unchanged name
  EXPECT_EQ(StreamAt(ud, 1, BoundaryLocation::TOP, 0).name_, "Stream 1");
  EXPECT_TRUE(Cfg(ud).IsProduct(1, "Stream 1"));
  EXPECT_EQ(PurgeTank(ud).source_endpoint_name_, "Stream 1");
}

// The KPIConfig remaps on their own (what the wrappers delegate to).
TEST(EndpointRefRename, KpiConfigOnStreamRenamedHonoursLocation) {
  KPIConfig k;
  k.SyncToComponents({ CO2 });
  ComponentKPIConfig& c = *k.FindByComponent(CO2);
  c.SetFeed(0, "S", BoundaryLocation::BOTTOM, true);
  c.SetProduct(0, "S", BoundaryLocation::TOP, true);
  c.SetProduct(1, "S", BoundaryLocation::TOP, true);

  k.OnStreamRenamed(0, "S", "T", BoundaryLocation::BOTTOM);
  EXPECT_TRUE(c.IsFeed(0, "T"));
  EXPECT_TRUE(c.IsProduct(0, "S"));

  k.OnStreamRenamed(0, "S", "U");
  EXPECT_TRUE(c.IsProduct(0, "U"));
  EXPECT_TRUE(c.IsProduct(1, "S"));  // other step untouched
}

} // namespace adsim
