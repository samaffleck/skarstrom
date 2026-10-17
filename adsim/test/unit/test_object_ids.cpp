// test_object_ids.cpp
//
// Stable object IDs (Step::id_, endpoint_id_ on streams and pressure
// boundaries) and the references that use them: KPI endpoint refs,
// buffer-tank sources and parameter paths.
//
//   * Decks saved before IDs existed get IDs on load, deterministically (steps
//     1..N, endpoints 1..M in cycle order, bottom boundary before top), and
//     their name-based references are bound to those IDs.
//   * IDs survive save/load, including gaps left by deleted steps.
//   * Renaming or reordering steps and streams WITHOUT going through the
//     UserData rename helpers (as any editor or script might) no longer breaks
//     references: they resolve by ID and refresh their cached names.
//     (bug_review.md E-11 / E-20 / W-07 were this class of bug.)
//   * An old deck runs to the same numbers before and after it is re-saved
//     with IDs.

#include <filesystem>
#include <set>
#include <sstream>
#include <string>
#include <vector>

#include <gtest/gtest.h>

#include "adsim/boundary_conditions/boundary_conditions.h"
#include "adsim/boundary_conditions/boundary_stream.h"
#include "adsim/boundary_conditions/buffer_tank_source.h"
#include "adsim/boundary_conditions/step.h"
#include "adsim/objects/user_data.h"
#include "adsim/parameters/parameter_registry.h"
#include "adsim/simulation.h"
#include "adsim/solver/kpi_config.h"

namespace adsim {
namespace {

const std::filesystem::path kExamples = std::filesystem::current_path() / ".." / "examples" / "decks";
const std::string CO2 = "CarbonDioxide";

std::unique_ptr<FlowBoundaryCondition> Flow(BoundaryLocation loc, std::initializer_list<const char*> names) {
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

PressureBoundaryCondition& PressureAt(UserData& ud, size_t step, BoundaryLocation side) {
  Step& s = ud.cycle.step_list_.at(step);
  auto& slot = side == BoundaryLocation::TOP ? s.top_boundary_condition_ : s.bottom_boundary_condition_;
  return dynamic_cast<PressureBoundaryCondition&>(*slot);
}

// Steps (endpoint IDs in brackets, assigned bottom-then-top):
//   0 "Pressurise": bottom pressure [1], top flow ["Stream 1"] [2]
//   1 "Feed":       bottom pressure [3], top flow ["Stream 1"] [4]
//   2 "Purge":      bottom flow ["Stream 1"] [5], top flow ["Purge In"] [6]
// KPI (CO2): feed = Feed's pressure boundary; product = Feed's "Stream 1".
// Buffer tank on Purge's "Purge In", sourced from Feed / "Stream 1".
UserData BuildDeck() {
  UserData ud;
  for (const char* n : {"Pressurise", "Feed", "Purge"}) ud.cycle.AddStep(n);
  auto& steps = ud.cycle.step_list_;
  steps[0].top_boundary_condition_ = Flow(BoundaryLocation::TOP, {"Stream 1"});
  steps[0].bottom_boundary_condition_ = Pressure(BoundaryLocation::BOTTOM);
  steps[1].top_boundary_condition_ = Flow(BoundaryLocation::TOP, {"Stream 1"});
  steps[1].bottom_boundary_condition_ = Pressure(BoundaryLocation::BOTTOM);
  steps[2].top_boundary_condition_ = Flow(BoundaryLocation::TOP, {"Purge In"});
  steps[2].bottom_boundary_condition_ = Flow(BoundaryLocation::BOTTOM, {"Stream 1"});
  StreamAt(ud, 2, BoundaryLocation::TOP, 0).buffer_tank_.Activate("Feed", "Stream 1");

  ud.cycle.EnsureIds();
  ud.fluid.components_ = {CO2};  // the KPI config is synced to these on load
  ud.kpi_config.SyncToComponents({CO2});
  ComponentKPIConfig* c = ud.kpi_config.FindByComponent(CO2);
  c->tracked = true;
  c->SetFeed(1, "Feed Pressure Boundary Bottom", BoundaryLocation::BOTTOM, true);
  c->SetProduct(1, "Stream 1", BoundaryLocation::TOP, true);
  ud.kpi_config.ResolveIds(ud.cycle.step_list_);
  ud.cycle.ResolveBufferTankReferences();
  return ud;
}

std::string Save(UserData& ud) {
  std::string s;
  SaveToString(ud, s);
  return s;
}

}  // namespace

TEST(ObjectIds, AssignedInCycleOrderBottomBeforeTop) {
  UserData ud = BuildDeck();
  const auto& s = ud.cycle.step_list_;
  EXPECT_EQ(s[0].id_, 1u);
  EXPECT_EQ(s[1].id_, 2u);
  EXPECT_EQ(s[2].id_, 3u);
  EXPECT_EQ(PressureAt(ud, 0, BoundaryLocation::BOTTOM).endpoint_id_, 1u);
  EXPECT_EQ(StreamAt(ud, 0, BoundaryLocation::TOP, 0).endpoint_id_, 2u);
  EXPECT_EQ(PressureAt(ud, 1, BoundaryLocation::BOTTOM).endpoint_id_, 3u);
  EXPECT_EQ(StreamAt(ud, 1, BoundaryLocation::TOP, 0).endpoint_id_, 4u);
  EXPECT_EQ(StreamAt(ud, 2, BoundaryLocation::BOTTOM, 0).endpoint_id_, 5u);
  EXPECT_EQ(StreamAt(ud, 2, BoundaryLocation::TOP, 0).endpoint_id_, 6u);

  const auto* c = ud.kpi_config.FindByComponent(CO2);
  ASSERT_EQ(c->feed_endpoints.size(), 1u);
  EXPECT_EQ(c->feed_endpoints[0].endpoint_id, 3u);
  EXPECT_EQ(c->feed_endpoints[0].step_id, 2u);
  EXPECT_EQ(c->product_endpoints[0].endpoint_id, 4u);
  const auto& bt = StreamAt(ud, 2, BoundaryLocation::TOP, 0).buffer_tank_;
  EXPECT_EQ(bt.source_step_id_, 2u);
  EXPECT_EQ(bt.source_endpoint_id_, 4u);
}

TEST(ObjectIds, EnsureIdsFillsGapsAndRenumbersRepeats) {
  UserData ud = BuildDeck();
  // Remove the middle step: IDs of the others do not change.
  ud.RemoveStep(1);
  EXPECT_EQ(ud.cycle.step_list_[0].id_, 1u);
  EXPECT_EQ(ud.cycle.step_list_[1].id_, 3u);
  // A new step and a duplicate get fresh IDs; the original keeps its own.
  ud.cycle.AddStep("Blowdown");
  EXPECT_EQ(ud.cycle.step_list_.back().id_, 4u);
  const size_t dup = ud.DuplicateStep(0);
  EXPECT_EQ(ud.cycle.step_list_[0].id_, 1u);
  EXPECT_EQ(ud.cycle.step_list_[dup].id_, 5u);
  EXPECT_EQ(PressureAt(ud, 0, BoundaryLocation::BOTTOM).endpoint_id_, 1u);
  EXPECT_NE(PressureAt(ud, dup, BoundaryLocation::BOTTOM).endpoint_id_, 1u);

  std::set<size_t> ids;
  for (const auto& h : ud.cycle.Endpoints()) EXPECT_TRUE(ids.insert(h.id).second) << "repeated endpoint ID " << h.id;

  // Changing a stream's device type keeps its endpoint ID.
  auto& slot = dynamic_cast<FlowBoundaryCondition&>(*ud.cycle.step_list_[0].top_boundary_condition_).streams_[0];
  const size_t before = slot->endpoint_id_;
  ud.cycle.ReplaceStream(slot, CreateStream(BoundaryStreamType::VALVE));
  EXPECT_EQ(slot->endpoint_id_, before);
}

TEST(ObjectIds, SurviveSaveAndLoad) {
  UserData ud = BuildDeck();
  ud.RemoveStep(0);  // leaves a gap: step IDs 2, 3
  UserData back;
  LoadFromString(back, Save(ud));
  ASSERT_EQ(back.cycle.step_list_.size(), 2u);
  for (const auto& s : back.cycle.step_list_) {
    ASSERT_NE(s.bottom_boundary_condition_, nullptr) << s.step_name_;
    ASSERT_NE(s.top_boundary_condition_, nullptr) << s.step_name_;
  }
  EXPECT_EQ(back.cycle.step_list_[0].id_, 2u);
  EXPECT_EQ(back.cycle.step_list_[1].id_, 3u);
  EXPECT_EQ(PressureAt(back, 0, BoundaryLocation::BOTTOM).endpoint_id_, 3u);
  const auto* c = back.kpi_config.FindByComponent(CO2);
  ASSERT_NE(c, nullptr);
  ASSERT_EQ(c->product_endpoints.size(), 1u);
  EXPECT_EQ(c->product_endpoints[0].endpoint_id, 4u);
  EXPECT_EQ(c->product_endpoints[0].step_index, 0);
  EXPECT_EQ(StreamAt(back, 1, BoundaryLocation::TOP, 0).buffer_tank_.source_endpoint_id_, 4u);
}

// The bug class of E-11 / E-20 / W-07: edits that bypass the rename helpers.
TEST(ObjectIds, DirectRenamesAndReordersKeepReferences) {
  UserData ud = BuildDeck();
  // Rename the source step and the product stream directly, then reorder.
  ud.cycle.step_list_[1].step_name_ = "Adsorb";
  StreamAt(ud, 1, BoundaryLocation::TOP, 0).name_ = "Product";
  std::swap(ud.cycle.step_list_[0], ud.cycle.step_list_[1]);  // no KPI remap

  ud.kpi_config.Prune(ud.cycle.step_list_);
  const auto* c = ud.kpi_config.FindByComponent(CO2);
  ASSERT_EQ(c->feed_endpoints.size(), 1u) << "feed ref was dropped";
  ASSERT_EQ(c->product_endpoints.size(), 1u) << "product ref was dropped";
  EXPECT_EQ(c->feed_endpoints[0].step_index, 0);
  EXPECT_EQ(c->feed_endpoints[0].endpoint_name, "Adsorb Pressure Boundary Bottom");
  EXPECT_EQ(c->product_endpoints[0].step_index, 0);
  EXPECT_EQ(c->product_endpoints[0].endpoint_name, "Product");
  EXPECT_EQ(c->product_endpoints[0].location, BoundaryLocation::TOP);

  ud.cycle.ResolveBufferTankSources();
  auto& bt = StreamAt(ud, 2, BoundaryLocation::TOP, 0).buffer_tank_;
  EXPECT_EQ(bt.source_step_name_, "Adsorb");
  EXPECT_EQ(bt.source_endpoint_name_, "Product");
  EXPECT_EQ(bt.summary_, &StreamAt(ud, 0, BoundaryLocation::TOP, 0).cycle_summary_);
}

TEST(ObjectIds, BufferTankOnABottomPressureBoundaryBindsToThatBoundary) {
  // Two pressure boundaries on the source step share the name "Pressure
  // Boundary"; by name the top one wins, by ID the stored one does.
  UserData ud;
  for (const char* n : {"Source", "Consumer"}) ud.cycle.AddStep(n);
  ud.cycle.step_list_[0].bottom_boundary_condition_ = Pressure(BoundaryLocation::BOTTOM);
  ud.cycle.step_list_[0].top_boundary_condition_ = Pressure(BoundaryLocation::TOP);
  ud.cycle.step_list_[1].bottom_boundary_condition_ = Flow(BoundaryLocation::BOTTOM, {"In"});
  ud.cycle.step_list_[1].top_boundary_condition_ = Pressure(BoundaryLocation::TOP);
  ud.cycle.EnsureIds();
  auto& bt = StreamAt(ud, 1, BoundaryLocation::BOTTOM, 0).buffer_tank_;
  bt.Activate("Source", kPressureBoundaryEndpointName);
  ud.cycle.ResolveBufferTankSources();
  EXPECT_EQ(bt.summary_, &PressureAt(ud, 0, BoundaryLocation::TOP).cycle_summary_);
  bt.source_endpoint_id_ = PressureAt(ud, 0, BoundaryLocation::BOTTOM).endpoint_id_;
  ud.cycle.ResolveBufferTankSources();
  EXPECT_EQ(bt.summary_, &PressureAt(ud, 0, BoundaryLocation::BOTTOM).cycle_summary_);
}

TEST(ObjectIds, ParameterPathsUseIds) {
  UserData ud = BuildDeck();
  ParameterRegistry reg(ud);
  ASSERT_NE(reg.Find("step[2].time"), nullptr);
  ASSERT_NE(reg.Find("step[Feed].time"), nullptr);
  EXPECT_EQ(reg.Find("step[Feed].time"), reg.Find("step[2].time"));
  const double* feed_time = &ud.cycle.step_list_[1].step_time_.Get();
  EXPECT_EQ(reg.Find("step[2].time")->ValuePointer(), feed_time);

  ud.cycle.step_list_[1].step_name_ = "Adsorb";
  std::swap(ud.cycle.step_list_[0], ud.cycle.step_list_[1]);
  reg.Refresh();
  // Same path, same step (now at index 0), new alias.
  EXPECT_EQ(reg.Find("step[2].time")->ValuePointer(), &ud.cycle.step_list_[0].step_time_.Get());
  EXPECT_NE(reg.Find("step[Adsorb].time"), nullptr);
  EXPECT_EQ(reg.Find("step[Feed].time"), nullptr);
}

// Example decks predate IDs: they load with deterministic IDs, their KPI refs
// bind, and a re-saved copy (with IDs) runs to exactly the same numbers.
TEST(ObjectIds, OldDeckLoadsDeterministicallyAndRunsIdentically) {
  const std::string deck = (kExamples / "Perez2019_A.json").string();
  UserData a, b;
  LoadFromJson(a, deck);
  LoadFromJson(b, deck);
  ASSERT_EQ(a.cycle.step_list_.size(), 4u);
  for (size_t i = 0; i < a.cycle.step_list_.size(); ++i) {
    EXPECT_EQ(a.cycle.step_list_[i].id_, i + 1);
    EXPECT_EQ(b.cycle.step_list_[i].id_, i + 1);
  }
  const auto ea = a.cycle.Endpoints();
  for (size_t i = 0; i < ea.size(); ++i) EXPECT_EQ(ea[i].id, i + 1);
  const auto* c = a.kpi_config.FindByComponent(CO2);
  ASSERT_NE(c, nullptr);
  for (const auto* list : {&c->feed_endpoints, &c->product_endpoints})
    for (const auto& r : *list) {
      EXPECT_NE(r.endpoint_id, 0u) << r.endpoint_name;
      auto h = a.cycle.FindEndpoint(r.endpoint_id);
      ASSERT_TRUE(h.has_value());
      EXPECT_EQ(h->log_name, r.endpoint_name);
      EXPECT_EQ(static_cast<int>(h->step_index), r.step_index);
    }

  Simulation old_deck = Simulation::FromFile(deck);
  Simulation resaved = Simulation::FromJson(old_deck.ToJson());
  for (Simulation* s : {&old_deck, &resaved}) {
    s->SetParameter("cycle.number_of_cycles", 2);
    ASSERT_TRUE(s->Run().Ok());
  }
  const auto k1 = old_deck.KPIs();
  const auto k2 = resaved.KPIs();
  ASSERT_EQ(k1.size(), k2.size());
  ASSERT_FALSE(k1.empty());
  for (size_t i = 0; i < k1.size(); ++i) {
    EXPECT_EQ(k1[i].purity, k2[i].purity);
    EXPECT_EQ(k1[i].recovery, k2[i].recovery);
    EXPECT_EQ(k1[i].energy_consumption, k2[i].energy_consumption);
  }
  EXPECT_EQ(old_deck.Series("pressure").values, resaved.Series("pressure").values);
}

}  // namespace adsim
