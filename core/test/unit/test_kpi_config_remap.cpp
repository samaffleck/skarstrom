// test_kpi_config_remap.cpp
//
// Pins KPI endpoint refs (keyed by step index + endpoint name) across the
// structural cycle edits that go through UserData: duplicate, remove, swap
// (move left/right) and rename. Before E-11 none of these remapped the refs,
// so a ref silently resolved to whichever step now sat at the old index when
// it had a same-named endpoint ("Stream 1"), and a rename made Prune drop the
// refs to the step's pressure boundary.

#include <gtest/gtest.h>

#include "adsim/objects/user_data.h"
#include "adsim/solver/kpi_config.h"
#include "adsim/boundary_conditions/step.h"
#include "adsim/boundary_conditions/boundary_conditions.h"
#include "adsim/boundary_conditions/boundary_stream.h"

namespace adsim {
namespace {

  const std::string CO2 = "CarbonDioxide";

  // Every step: TOP flow BC with one stream named "Stream 1" (the GUI default,
  // so every step shares that endpoint name) and a BOTTOM pressure BC.
  void AddStep(UserData& ud, const std::string& name) {
    ud.cycle.AddStep(name);
    Step& s = ud.cycle.GetStep(name);

    auto flow = std::make_unique<FlowBoundaryCondition>();
    flow->location_ = BoundaryLocation::TOP;
    auto stream = std::make_unique<BoundaryStream>();
    stream->name_ = "Stream 1";
    flow->streams_.push_back(std::move(stream));
    s.top_boundary_condition_ = std::move(flow);

    auto p = std::make_unique<PressureBoundaryCondition>();
    p->location_ = BoundaryLocation::BOTTOM;
    s.bottom_boundary_condition_ = std::move(p);
  }

  // Steps [Pressurise, Feed, Blowdown]. Feed (step 1) is the feed (its bottom
  // pressure boundary) and the product (its top "Stream 1"); Blowdown's
  // "Stream 1" (step 2) is also tagged as product.
  UserData BuildDeck() {
    UserData ud;
    AddStep(ud, "Pressurise");
    AddStep(ud, "Feed");
    AddStep(ud, "Blowdown");

    ud.kpi_config.SyncToComponents({ CO2 });
    ComponentKPIConfig* c = ud.kpi_config.FindByComponent(CO2);
    c->tracked = true;
    c->SetFeed(1, "Feed Pressure Boundary Bottom", BoundaryLocation::BOTTOM, true);
    c->SetProduct(1, "Stream 1", BoundaryLocation::TOP, true);
    c->SetProduct(2, "Stream 1", BoundaryLocation::TOP, true);
    ud.kpi_config.Prune(ud.cycle.step_list_);
    return ud;
  }

  const ComponentKPIConfig& Cfg(const UserData& ud) {
    return *ud.kpi_config.FindByComponent(CO2);
  }

  // Name of the step each ref of `v` resolves to after a Prune (the refs that
  // survive are exactly what the KPI maths will read).
  std::vector<std::string> ResolvedSteps(UserData& ud, bool product) {
    ud.kpi_config.Prune(ud.cycle.step_list_);
    const auto& v = product ? Cfg(ud).product_endpoints : Cfg(ud).feed_endpoints;
    std::vector<std::string> out;
    for (const auto& r : v) out.push_back(ud.cycle.step_list_.at(r.step_index).step_name_);
    return out;
  }

  using Names = std::vector<std::string>;

} // namespace

TEST(KpiConfigRemap, BaselineResolvesToTaggedSteps) {
  UserData ud = BuildDeck();
  EXPECT_EQ(ResolvedSteps(ud, false), (Names{ "Feed" }));
  EXPECT_EQ(ResolvedSteps(ud, true), (Names{ "Feed", "Blowdown" }));
}

// Duplicating a step before the tagged ones shifts their indices +1; the copy
// carries no KPI roles.
TEST(KpiConfigRemap, DuplicateShiftsLaterRefs) {
  UserData ud = BuildDeck();
  ASSERT_EQ(ud.DuplicateStep(0), 1u);
  ASSERT_EQ(ud.cycle.step_list_[1].step_name_, "Pressurise (copy)");

  EXPECT_EQ(ResolvedSteps(ud, false), (Names{ "Feed" }));
  EXPECT_EQ(ResolvedSteps(ud, true), (Names{ "Feed", "Blowdown" }));
  EXPECT_TRUE(Cfg(ud).IsProduct(2, "Stream 1"));
  EXPECT_TRUE(Cfg(ud).IsProduct(3, "Stream 1"));
  EXPECT_FALSE(Cfg(ud).IsProduct(1, "Stream 1"));  // the copy is untagged
}

// Duplicating the tagged step itself leaves the original's refs in place and
// shifts only the refs after the inserted copy.
TEST(KpiConfigRemap, DuplicateTaggedStepKeepsOriginal) {
  UserData ud = BuildDeck();
  ASSERT_EQ(ud.DuplicateStep(1), 2u);

  EXPECT_EQ(ResolvedSteps(ud, false), (Names{ "Feed" }));
  EXPECT_EQ(ResolvedSteps(ud, true), (Names{ "Feed", "Blowdown" }));
  EXPECT_FALSE(Cfg(ud).IsProduct(2, "Stream 1"));  // "Feed (copy)"
}

TEST(KpiConfigRemap, OutOfRangeDuplicateIsNoOp) {
  UserData ud = BuildDeck();
  EXPECT_EQ(ud.DuplicateStep(99), static_cast<size_t>(-1));
  EXPECT_EQ(ResolvedSteps(ud, true), (Names{ "Feed", "Blowdown" }));
}

// Removing an untagged earlier step shifts later refs -1.
TEST(KpiConfigRemap, RemoveEarlierStepShiftsRefs) {
  UserData ud = BuildDeck();
  ud.RemoveStep(0);
  ASSERT_EQ(ud.cycle.step_list_.size(), 2u);

  EXPECT_EQ(ResolvedSteps(ud, false), (Names{ "Feed" }));
  EXPECT_EQ(ResolvedSteps(ud, true), (Names{ "Feed", "Blowdown" }));
}

// Removing the tagged step drops its refs instead of handing them to the step
// that slides into its index (Blowdown also has a "Stream 1").
TEST(KpiConfigRemap, RemoveTaggedStepDropsItsRefs) {
  UserData ud = BuildDeck();
  ud.RemoveStep(1);
  ASSERT_EQ(ud.cycle.step_list_.size(), 2u);

  EXPECT_TRUE(ResolvedSteps(ud, false).empty());
  EXPECT_EQ(ResolvedSteps(ud, true), (Names{ "Blowdown" }));
  EXPECT_EQ(Cfg(ud).product_endpoints.size(), 1u);
}

// Move left / move right (adjacent swaps): refs follow their step.
TEST(KpiConfigRemap, SwapStepsRefsFollowStep) {
  UserData ud = BuildDeck();
  ud.SwapSteps(1, 0);  // move Feed left -> [Feed, Pressurise, Blowdown]
  ASSERT_EQ(ud.cycle.step_list_[0].step_name_, "Feed");
  EXPECT_EQ(ResolvedSteps(ud, false), (Names{ "Feed" }));
  EXPECT_EQ(ResolvedSteps(ud, true), (Names{ "Feed", "Blowdown" }));
  EXPECT_TRUE(Cfg(ud).IsProduct(0, "Stream 1"));
  EXPECT_FALSE(Cfg(ud).IsProduct(1, "Stream 1"));

  ud.SwapSteps(1, 2);  // move Pressurise right -> [Feed, Blowdown, Pressurise]
  ASSERT_EQ(ud.cycle.step_list_[2].step_name_, "Pressurise");
  EXPECT_EQ(ResolvedSteps(ud, true), (Names{ "Feed", "Blowdown" }));
  EXPECT_TRUE(Cfg(ud).IsProduct(1, "Stream 1"));
  EXPECT_FALSE(Cfg(ud).IsProduct(2, "Stream 1"));
}

// Renaming a step rewrites refs to its pressure-boundary endpoint (whose name
// embeds the step name) so Prune keeps them; flow-stream refs are untouched.
TEST(KpiConfigRemap, RenameRewritesPressureEndpointRefs) {
  UserData ud = BuildDeck();
  ud.RenameStep(1, "Adsorption");
  ASSERT_EQ(ud.cycle.step_list_[1].step_name_, "Adsorption");

  EXPECT_EQ(ResolvedSteps(ud, false), (Names{ "Adsorption" }));
  ASSERT_EQ(Cfg(ud).feed_endpoints.size(), 1u);
  EXPECT_EQ(Cfg(ud).feed_endpoints[0].endpoint_name, "Adsorption Pressure Boundary Bottom");
  EXPECT_EQ(Cfg(ud).feed_endpoints[0].location, BoundaryLocation::BOTTOM);
  EXPECT_TRUE(Cfg(ud).IsProduct(1, "Stream 1"));

  // Renaming a different step leaves these refs alone.
  ud.RenameStep(0, "Repressurise");
  EXPECT_EQ(Cfg(ud).feed_endpoints[0].endpoint_name, "Adsorption Pressure Boundary Bottom");
}

// The endpoint name used by the rename remap is the one GetEndpoints exposes.
TEST(KpiConfigRemap, PressureEndpointNameMatchesGetEndpoints) {
  UserData ud = BuildDeck();
  const auto eps = ud.cycle.step_list_[1].GetEndpoints();
  bool found = false;
  for (const auto& ep : eps)
    if (ep.name == Step::PressureEndpointName("Feed", BoundaryLocation::BOTTOM)) found = true;
  EXPECT_TRUE(found);
}

// A sequence of edits (what a user does in the cycle editor) keeps every ref
// on its logical endpoint.
TEST(KpiConfigRemap, MixedEditSequence) {
  UserData ud = BuildDeck();
  ud.DuplicateStep(2);          // [P, Feed, Blowdown, Blowdown (copy)]
  ud.SwapSteps(0, 1);           // [Feed, P, Blowdown, Blowdown (copy)]
  ud.RenameStep(0, "Adsorb");   // [Adsorb, P, Blowdown, Blowdown (copy)]
  ud.RemoveStep(1);             // [Adsorb, Blowdown, Blowdown (copy)]

  EXPECT_EQ(ResolvedSteps(ud, false), (Names{ "Adsorb" }));
  EXPECT_EQ(ResolvedSteps(ud, true), (Names{ "Adsorb", "Blowdown" }));
  EXPECT_TRUE(Cfg(ud).IsFeed(0, "Adsorb Pressure Boundary Bottom"));
}

} // namespace adsim
