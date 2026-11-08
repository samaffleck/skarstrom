/*
 * Skarstrom - pressure and temperature swing adsorption simulator
 * Copyright (c) 2024-2026 Sam Affleck and contributors
 *
 * Licensed under the terms in the LICENSE file at the repository root.
 */
// The structural cascades of adsim/objects/deck_edit.h: layers, steps,
// streams, kind switches with parameter carry-over, IDs and connection repair.
// Each test checks the precondition, the effect and the invariants (the deck
// checks pass afterwards); RoundTrip checks that a deck edited through every
// cascade saves, reloads and saves to the same text.

#include <algorithm>
#include <filesystem>
#include <fstream>
#include <map>
#include <sstream>
#include <string>

#include <gtest/gtest.h>

#include "adsim/equations/isotherm_kinds.h"
#include "adsim/objects/deck_builder.h"
#include "adsim/objects/deck_edit.h"
#include "adsim/objects/deck_io.h"
#include "adsim/objects/user_data.h"

#include <cereal/external/rapidjson/document.h>
#include <cereal/external/rapidjson/stringbuffer.h>
#include <cereal/external/rapidjson/writer.h>

namespace adsim {
namespace {

const std::filesystem::path kExamples = std::filesystem::current_path() / ".." / "examples" / "decks";

UserData Load(const std::string& name) {
  std::ifstream in(kExamples / (name + ".json"), std::ios::binary);
  std::ostringstream ss;
  ss << in.rdbuf();
  UserData data;
  LoadDeckFlexible(data, ss.str());
  return data;
}

// The deck checks pass (expressions evaluated first).
void ExpectChecksPass(UserData& data, const char* label) {
  std::string err;
  ASSERT_TRUE(data.EvaluateAllExpressions(err)) << label << ": " << err;
  std::string checks;
  EXPECT_TRUE(data.PerformChecks(checks)) << label << ": " << checks;
}

std::string Saved(UserData& data) {
  data.info.sim_version = 0;
  std::string out;
  SaveToString(data, out);
  return out;
}

// The saved text with every unordered map (cereal writes them as arrays of
// {"key", "value"} objects, in hash order) sorted by key, so two saves of the
// same deck compare equal whatever order the maps iterate in. cereal writes a
// polymorphic type's name and class version only where it first meets them
// and numbers polymorphic ids in that order, so those markers depend on the
// map order too and are left out (the fields themselves differ by type).
// Loading a deck gives its connections fresh pair ids, so pair ids are
// renumbered by first appearance.
namespace rj = CEREAL_RAPIDJSON_NAMESPACE;

std::string Write(const rj::Value& v) {
  rj::StringBuffer buf;
  rj::Writer<rj::StringBuffer> w(buf);
  v.Accept(w);
  return buf.GetString();
}

void Canonicalise(rj::Value& v, rj::Document::AllocatorType& a, std::map<uint64_t, uint64_t>& pairs) {
  if (v.IsObject()) {
    if (v.HasMember("Pair_ID") && v["Pair_ID"].IsUint64() && v["Pair_ID"].GetUint64() != 0) {
      const uint64_t id = v["Pair_ID"].GetUint64();
      if (!pairs.count(id)) pairs[id] = pairs.size() + 1;
      v["Pair_ID"].SetUint64(pairs[id]);
    }
    for (const char* marker : {"polymorphic_id", "polymorphic_name", "cereal_class_version"}) v.RemoveMember(marker);
    for (auto& m : v.GetObject()) Canonicalise(m.value, a, pairs);
  } else if (v.IsArray()) {
    bool map = !v.Empty();
    for (auto& e : v.GetArray()) {
      Canonicalise(e, a, pairs);
      map = map && e.IsObject() && e.MemberCount() == 2 && e.HasMember("key") && e.HasMember("value");
    }
    if (!map) return;
    std::vector<std::pair<std::string, size_t>> order;
    for (size_t i = 0; i < v.Size(); ++i) order.emplace_back(Write(v[static_cast<rj::SizeType>(i)]["key"]), i);
    std::sort(order.begin(), order.end());
    rj::Value sorted(rj::kArrayType);
    for (const auto& [k, i] : order) sorted.PushBack(rj::Value(v[static_cast<rj::SizeType>(i)], a), a);
    v.Swap(sorted);
  }
}

std::string Canonical(UserData& data) {
  rj::Document doc;
  doc.Parse(Saved(data).c_str());
  std::map<uint64_t, uint64_t> pairs;
  Canonicalise(doc, doc.GetAllocator(), pairs);
  return Write(doc);
}

double Sum(const ExpressionMap& y) {
  double s = 0.0;
  for (const auto& [c, v] : y) s += v.Get();
  return s;
}

Step& StepById(UserData& d, size_t id) { return *d.cycle.FindStepById(id); }

Connection& ConnOf(const Cycle::EndpointHandle& h) {
  return h.stream ? h.stream->connection_ : h.pressure->connection_;
}

}  // namespace

// ---- layers -------------------------------------------------------------------------------------

TEST(DeckEditCascades, LayersAddRemoveMove) {
  UserData d = Load("Perez2019_A");
  const size_t before = d.bed.layers_.layers_.size();
  ASSERT_GE(before, 1u);
  const std::string first = d.bed.layers_.layers_.front().name_;

  const LayerID a = AddLayer(d);
  ASSERT_EQ(d.bed.layers_.layers_.size(), before + 1);
  const auto& la = d.bed.layers_.Get(a);
  EXPECT_EQ(la.name_, "Layer " + std::to_string(before + 1));
  EXPECT_EQ(la.composition_.size(), 1u) << "the composition is seeded";
  // A taken name is made unique.
  const LayerID b = AddLayer(d, first);
  EXPECT_EQ(d.bed.layers_.Get(b).name_, first + " 2");
  EXPECT_NE(a, b);
  ExpectChecksPass(d, "two layers added");

  // Move b towards the feed end and back; the ends refuse.
  std::string err;
  EXPECT_FALSE(MoveLayer(d, b, +1, err)) << "already at the product end";
  EXPECT_FALSE(MoveLayer(d, b, 2, err));
  ASSERT_TRUE(MoveLayer(d, b, -1, err)) << err;
  EXPECT_EQ(d.bed.layers_.layers_[before].id_, b);
  EXPECT_EQ(d.bed.layers_.layers_[before + 1].id_, a);
  ASSERT_TRUE(MoveLayer(d, b, +1, err)) << err;
  EXPECT_EQ(d.bed.layers_.layers_.back().id_, b);

  // An energy source on layer a goes with it.
  StepById(d, d.cycle.step_list_.front().id_).energy_sources_[a].enabled_ = true;
  ASSERT_TRUE(RemoveLayer(d, a, err)) << err;
  for (const auto& s : d.cycle.step_list_) EXPECT_EQ(s.energy_sources_.count(a), 0u);
  EXPECT_FALSE(RemoveLayer(d, a, err)) << "already removed";
  EXPECT_NE(err.find("no layer"), std::string::npos);
  ASSERT_TRUE(RemoveLayer(d, b, err)) << err;
  EXPECT_EQ(d.bed.layers_.layers_.size(), before);
  ExpectChecksPass(d, "layers removed");
}

// ---- steps --------------------------------------------------------------------------------------

TEST(DeckEditCascades, AddStepSeedsCompositionAndPassesChecks) {
  UserData d = Load("Perez2019_A");
  const size_t n = d.cycle.step_list_.size();
  const size_t id = AddStep(d);
  ASSERT_EQ(d.cycle.step_list_.size(), n + 1);
  Step& s = StepById(d, id);
  EXPECT_EQ(s.step_name_, "Step " + std::to_string(n + 1));
  auto* bottom = dynamic_cast<FlowBoundaryCondition*>(s.bottom_boundary_condition_.get());
  auto* top = dynamic_cast<PressureBoundaryCondition*>(s.top_boundary_condition_.get());
  ASSERT_TRUE(bottom && top);
  ASSERT_EQ(bottom->streams_.size(), 1u);
  EXPECT_EQ(bottom->direction_, BoundaryDirection::INLET);
  EXPECT_EQ(top->direction_, BoundaryDirection::OUTLET);
  for (const ExpressionMap* y : {&bottom->streams_[0]->mole_fraction, &top->mole_fraction_}) {
    EXPECT_EQ(y->size(), d.fluid.components_.size());
    EXPECT_DOUBLE_EQ(Sum(*y), 1.0);
  }
  EXPECT_NE(bottom->streams_[0]->endpoint_id_, 0u);
  // A taken name is made unique.
  const size_t id2 = AddStep(d, s.step_name_);
  EXPECT_EQ(StepById(d, id2).step_name_, s.step_name_ + " 2");
  ExpectChecksPass(d, "steps added");
}

TEST(DeckEditCascades, RemoveStepDisconnectsPartnersAndTanks) {
  UserData d = Load("FiveStepProductPurge");
  // A connected pair across two steps.
  std::optional<Cycle::EndpointHandle> a, b;
  for (const auto& h : d.cycle.Endpoints()) {
    if (!ConnOf(h).IsConnected()) continue;
    if (!a) a = h;
    else if (ConnOf(h).GetID() == ConnOf(*a).GetID()) b = h;
  }
  ASSERT_TRUE(a && b);
  ASSERT_NE(a->step->id_, b->step->id_);
  const size_t gone = a->step->id_, partner = b->id;

  // A buffer tank on another step drawing from an endpoint of `gone`. The
  // tank needs an unconnected inlet: give a third step a new inlet stream.
  const size_t extra = AddStep(d, "Tank user");
  const size_t inlet = dynamic_cast<FlowBoundaryCondition*>(StepById(d, extra).bottom_boundary_condition_.get())
                           ->streams_[0]->endpoint_id_;
  build::SetBufferTank(d, inlet, a->id);
  ASSERT_TRUE(d.cycle.FindEndpoint(inlet)->stream->buffer_tank_.IsActive());

  std::string err;
  ASSERT_TRUE(RemoveStep(d, gone, err)) << err;
  EXPECT_EQ(d.cycle.FindStepById(gone), nullptr);
  EXPECT_FALSE(ConnOf(*d.cycle.FindEndpoint(partner)).IsConnected()) << "the partner is disconnected";
  EXPECT_FALSE(d.cycle.FindEndpoint(inlet)->stream->buffer_tank_.IsActive()) << "the tank is deactivated";
  EXPECT_FALSE(RemoveStep(d, gone, err));
  EXPECT_EQ(RepairConnections(d), 0u) << "no orphan left";
}

TEST(DeckEditCascades, DuplicateAndMoveStep) {
  UserData d = Load("Perez2019_A");
  const size_t first = d.cycle.step_list_[0].id_;
  const size_t second = d.cycle.step_list_[1].id_;
  size_t copy = 0;
  std::string err;
  ASSERT_TRUE(DuplicateStep(d, first, copy, err)) << err;
  EXPECT_EQ(d.cycle.step_list_[1].id_, copy);
  EXPECT_NE(copy, first);
  EXPECT_EQ(d.cycle.step_list_[1].step_name_, d.cycle.step_list_[0].step_name_ + " (copy)");

  EXPECT_FALSE(MoveStep(d, first, -1, err)) << "already first";
  ASSERT_TRUE(MoveStep(d, first, +1, err)) << err;
  EXPECT_EQ(d.cycle.step_list_[0].id_, copy);
  EXPECT_EQ(d.cycle.step_list_[1].id_, first);
  EXPECT_EQ(d.cycle.step_list_[2].id_, second);
  EXPECT_FALSE(DuplicateStep(d, 9999, copy, err));
  ExpectChecksPass(d, "duplicated and moved");
}

// ---- streams ------------------------------------------------------------------------------------

TEST(DeckEditCascades, AddAndRemoveStreams) {
  UserData d = Load("Perez2019_A");
  const size_t step = AddStep(d, "S");
  std::string err;
  size_t valve = 0, pump = 0, flow = 0;
  EXPECT_FALSE(AddStream(d, step, BoundaryLocation::TOP, BoundaryStreamType::VALVE, valve, err))
      << "the top is a pressure boundary";
  ASSERT_TRUE(AddStream(d, step, BoundaryLocation::BOTTOM, BoundaryStreamType::VALVE, valve, err)) << err;
  ASSERT_TRUE(AddStream(d, step, BoundaryLocation::BOTTOM, BoundaryStreamType::PUMP, pump, err)) << err;
  auto* fbc = dynamic_cast<FlowBoundaryCondition*>(StepById(d, step).bottom_boundary_condition_.get());
  ASSERT_EQ(fbc->streams_.size(), 3u);
  EXPECT_EQ(fbc->streams_[1]->name_, "Stream 2");
  EXPECT_EQ(fbc->streams_[2]->name_, "Stream 3");
  EXPECT_DOUBLE_EQ(Sum(fbc->streams_[2]->mole_fraction), 1.0);
  auto* p = dynamic_cast<PumpBoundaryStream*>(fbc->streams_[2].get());
  ASSERT_TRUE(p && p->pump_);
  EXPECT_EQ(p->pump_->GetType(), PumpType::COMPRESSOR) << "an inlet's pump is a compressor";
  EXPECT_TRUE(p->pump_->use_pump_);

  // Remove the valve stream. The next name counts on from the stream count
  // ("Stream 3", taken, so "Stream 4").
  ASSERT_TRUE(RemoveStream(d, step, BoundaryLocation::BOTTOM, 1, err)) << err;
  ASSERT_NE(valve, 0u);
  EXPECT_EQ(d.cycle.FindEndpoint(valve), std::nullopt);
  ASSERT_TRUE(AddStream(d, step, BoundaryLocation::BOTTOM, BoundaryStreamType::NONE, flow, err)) << err;
  EXPECT_EQ(fbc->streams_.back()->name_, "Stream 4");
  EXPECT_FALSE(RemoveStream(d, step, BoundaryLocation::BOTTOM, 7, err));
  ExpectChecksPass(d, "streams added and removed");
}

TEST(DeckEditCascades, RemoveConnectedStreamDisconnectsPartner) {
  UserData d = Load("FiveStepProductPurge");
  for (auto& step : d.cycle.step_list_) {
    for (auto side : {BoundaryLocation::BOTTOM, BoundaryLocation::TOP}) {
      auto& bc = side == BoundaryLocation::TOP ? step.top_boundary_condition_ : step.bottom_boundary_condition_;
      auto* fbc = dynamic_cast<FlowBoundaryCondition*>(bc.get());
      if (!fbc) continue;
      for (size_t i = 0; i < fbc->streams_.size(); ++i) {
        if (!fbc->streams_[i]->connection_.IsConnected()) continue;
        const size_t pid = fbc->streams_[i]->connection_.GetID();
        std::string err;
        ASSERT_TRUE(RemoveStream(d, step.id_, side, i, err)) << err;
        for (const auto& h : d.cycle.Endpoints())
          EXPECT_FALSE(ConnOf(h).IsConnected() && ConnOf(h).GetID() == pid) << "partner still connected";
        return;
      }
    }
  }
  FAIL() << "no connected stream in FiveStepProductPurge";
}

// ---- IDs and connections ----------------------------------------------------------------------------

TEST(DeckEditCascades, EnsureIdsNumbersPairsByFirstAppearance) {
  UserData d = Load("FiveStepProductPurge");
  EnsureIds(d);
  size_t expected = 0;
  std::vector<size_t> seen;
  for (const auto& h : d.cycle.Endpoints()) {
    const Connection& c = ConnOf(h);
    if (!c.IsConnected()) continue;
    if (std::find(seen.begin(), seen.end(), c.GetID()) == seen.end()) {
      EXPECT_EQ(c.GetID(), ++expected);
      seen.push_back(c.GetID());
    }
  }
  ASSERT_GT(expected, 0u);
  const std::string once = Saved(d);
  EnsureIds(d);
  EXPECT_EQ(Saved(d), once) << "idempotent";
  ExpectChecksPass(d, "renumbered");
}

TEST(DeckEditCascades, ConnectAndRepair) {
  UserData d = Load("Perez2019_A");
  const size_t s1 = AddStep(d, "Giver"), s2 = AddStep(d, "Taker");
  const size_t e1 = dynamic_cast<PressureBoundaryCondition*>(StepById(d, s1).top_boundary_condition_.get())->endpoint_id_;
  const size_t e2 = dynamic_cast<PressureBoundaryCondition*>(StepById(d, s2).top_boundary_condition_.get())->endpoint_id_;
  const size_t f2 = dynamic_cast<FlowBoundaryCondition*>(StepById(d, s2).bottom_boundary_condition_.get())
                        ->streams_[0]->endpoint_id_;
  std::string err;
  EXPECT_FALSE(ConnectEndpoints(d, e1, e2, err)) << "two pressure boundaries";
  EXPECT_FALSE(ConnectEndpoints(d, f2, f2, err));
  ASSERT_TRUE(ConnectEndpoints(d, e1, f2, err)) << err;
  EXPECT_FALSE(ConnectEndpoints(d, e1, f2, err)) << "already connected";
  auto h1 = *d.cycle.FindEndpoint(e1), h2 = *d.cycle.FindEndpoint(f2);
  EXPECT_EQ(ConnOf(h1).stream_, &h2.stream->stream_data_);
  EXPECT_EQ(ConnOf(h2).stream_, &h1.pressure->stream_data_);
  EXPECT_EQ(RepairConnections(d), 0u);

  // An orphan: erase the partner by hand.
  auto* fbc = dynamic_cast<FlowBoundaryCondition*>(StepById(d, s2).bottom_boundary_condition_.get());
  fbc->streams_.clear();
  EXPECT_EQ(RepairConnections(d), 1u);
  EXPECT_FALSE(ConnOf(*d.cycle.FindEndpoint(e1)).IsConnected());
}

TEST(DeckEditCascades, FillDeckGapsIsIdempotentOnExamples) {
  for (const char* name : {"Perez2019_A", "FiveStepProductPurge"}) {
    UserData d = Load(name);
    // An example deck already has everything; filling it changes no value.
    FillDeckGaps(d);
    const std::string once = Saved(d);
    FillDeckGaps(d);
    EXPECT_EQ(Saved(d), once) << name;
    // A component added without its cascade gets its entries.
    d.fluid.AddComponent("Argon");
    FillDeckGaps(d);
    EXPECT_EQ(d.ic.mole_fraction.count("Argon"), 1u) << name;
    EXPECT_EQ(d.fluid.adsorbed_heat_capacity_.count("Argon"), 1u) << name;
    for (const auto& [id, ads] : d.adsorbent_directory.GetAll())
      EXPECT_EQ(ads.adsorption_source_.mass_transfer_models_.count("Argon"), 1u) << name;
    EXPECT_NE(d.kpi_config.FindByComponent("Argon"), nullptr) << name;
  }
}

// ---- kind switches ----------------------------------------------------------------------------------

namespace {
// One adsorbent with one isotherm parameter set we can recognise.
struct IsoFixture {
  UserData d;
  AdsorbentID ads = 0;
  IsoFixture() {
    build::AddComponent(d, "CarbonDioxide");
    build::AddComponent(d, "Nitrogen");
    ads = build::AddAdsorbent(d, "A");
  }
  ExpressionParam* Field(const std::string& c, const std::string& name) {
    auto& src = d.adsorbent_directory.Get(ads).adsorption_source_;
    for (auto& [n, p] : IsothermFields(src.GetIsotherm(c).get(), c))
      if (n == name) return p;
    return nullptr;
  }
};
}  // namespace

TEST(DeckEditCascades, SetIsothermKindCarriesSharedParameters) {
  IsoFixture f;
  std::string err;
  ASSERT_TRUE(SetIsothermKind(f.d, f.ads, "CarbonDioxide", IsothermType::DUAL_SITE_LANGMUIR_P, err)) << err;
  f.Field("CarbonDioxide", "m1")->SetExpr("4.5");
  f.Field("CarbonDioxide", "m1")->Get() = 4.5;
  *f.Field("CarbonDioxide", "b1") = 2e-9;
  *f.Field("CarbonDioxide", "q2") = 31000.0;

  // The paired kind keeps every parameter.
  ASSERT_TRUE(SetIsothermKind(f.d, f.ads, "CarbonDioxide", IsothermType::EXTENDED_DUAL_SITE_LANGMUIR_P, err)) << err;
  EXPECT_DOUBLE_EQ(f.Field("CarbonDioxide", "m1")->Get(), 4.5);
  EXPECT_DOUBLE_EQ(f.Field("CarbonDioxide", "b1")->Get(), 2e-9);
  EXPECT_DOUBLE_EQ(f.Field("CarbonDioxide", "q2")->Get(), 31000.0);

  // The concentration basis keeps m and q, not b (1/Pa vs m3/mol).
  ASSERT_TRUE(SetIsothermKind(f.d, f.ads, "CarbonDioxide", IsothermType::DUAL_SITE_LANGMUIR_C, err)) << err;
  EXPECT_DOUBLE_EQ(f.Field("CarbonDioxide", "m1")->Get(), 4.5);
  EXPECT_DOUBLE_EQ(f.Field("CarbonDioxide", "q2")->Get(), 31000.0);
  EXPECT_DOUBLE_EQ(f.Field("CarbonDioxide", "b1")->Get(),
                   IsothermKind(IsothermType::DUAL_SITE_LANGMUIR_C).Find("b1")->default_value);

  // The other component is untouched; unknown inputs are refused.
  EXPECT_EQ(f.d.adsorbent_directory.Get(f.ads).adsorption_source_.IsothermMap().at("Nitrogen"), IsothermType::INERT);
  EXPECT_FALSE(SetIsothermKind(f.d, f.ads, "Argon", IsothermType::SIPS, err));
  EXPECT_FALSE(SetIsothermKind(f.d, 99, "Nitrogen", IsothermType::SIPS, err));
}

TEST(DeckEditCascades, SetMassTransferKindCarriesSharedParameters) {
  IsoFixture f;
  auto& mtm = f.d.adsorbent_directory.Get(f.ads).adsorption_source_.mass_transfer_models_;
  std::string err;
  ASSERT_TRUE(SetMassTransferKind(f.d, f.ads, "CarbonDioxide", MassTransferModelType::LDF,
                                  MassTransferCoefficientType::MACROPORE, err)) << err;
  auto* m = dynamic_cast<MassTransferCoefficientMacropore*>(mtm.at("CarbonDioxide")->mtc_model_.get());
  ASSERT_TRUE(m);
  m->tortuosity_ = 4.2;
  m->macropore_void_fraction_.SetExpr("0.3");
  ASSERT_TRUE(SetMassTransferKind(f.d, f.ads, "CarbonDioxide", MassTransferModelType::LDF,
                                  MassTransferCoefficientType::MACROPORE_KNUDSEN, err)) << err;
  auto* k = dynamic_cast<MassTransferCoefficientMacroporeKnudsen*>(mtm.at("CarbonDioxide")->mtc_model_.get());
  ASSERT_TRUE(k);
  EXPECT_DOUBLE_EQ(k->tortuosity_.Get(), 4.2);
  EXPECT_EQ(k->macropore_void_fraction_.GetExpr(), "0.3") << "the expression text carries over";

  // The rate constant does not survive a change of rate-model order.
  ASSERT_TRUE(SetMassTransferKind(f.d, f.ads, "Nitrogen", MassTransferModelType::LDF,
                                  MassTransferCoefficientType::CONSTANT, err)) << err;
  dynamic_cast<MassTransferCoefficientConstant*>(mtm.at("Nitrogen")->mtc_model_.get())->constant_mtc_ = 0.07;
  ASSERT_TRUE(SetMassTransferKind(f.d, f.ads, "Nitrogen", MassTransferModelType::LDF2,
                                  MassTransferCoefficientType::CONSTANT, err)) << err;
  EXPECT_EQ(mtm.at("Nitrogen")->GetModelType(), MassTransferModelType::LDF2);
  EXPECT_DOUBLE_EQ(dynamic_cast<MassTransferCoefficientConstant*>(mtm.at("Nitrogen")->mtc_model_.get())->constant_mtc_.Get(),
                   0.5);
}

TEST(DeckEditCascades, LayerModelKindsCarryShared) {
  UserData d = Load("Perez2019_A");
  const LayerID l = d.bed.layers_.layers_.front().id_;
  std::string err;
  ASSERT_TRUE(SetDispersionKind(d, l, AxialDispersionType::RUTHVEN, err)) << err;
  dynamic_cast<DispersionModelRuthven*>(d.bed.layers_.Get(l).equations_.dispersion_model_.get())->Dm = 2.5e-5;
  ASSERT_TRUE(SetDispersionKind(d, l, AxialDispersionType::TAYLOR_ARIS, err)) << err;
  EXPECT_DOUBLE_EQ(
      dynamic_cast<DispersionModelTaylorAris*>(d.bed.layers_.Get(l).equations_.dispersion_model_.get())->Dm.Get(), 2.5e-5);
  ASSERT_TRUE(SetDispersionKind(d, l, AxialDispersionType::RUTHVEN, err)) << err;
  ASSERT_TRUE(SetInternalHeatTransferKind(d, l, InternalHeatTransferType::LEVA_1947, err)) << err;
  ASSERT_TRUE(SetExternalHeatTransferKind(d, l, ExternalHeatTransferType::HOLMAN, err)) << err;
  ASSERT_TRUE(SetHeatDispersionKind(d, l, AxialHeatDispersionType::RUTHVEN, err)) << err;
  ASSERT_TRUE(SetParticleHeatTransferKind(d, ParticleHeatTransferType::CONSTANT, err)) << err;
  EXPECT_FALSE(SetDispersionKind(d, 999, AxialDispersionType::RUTHVEN, err));
  ExpectChecksPass(d, "layer models switched");
}

TEST(DeckEditCascades, SetBoundaryConditionKindCarriesInlet) {
  UserData d = Load("Perez2019_A");
  const size_t s = AddStep(d, "B");
  auto* fbc = dynamic_cast<FlowBoundaryCondition*>(StepById(d, s).bottom_boundary_condition_.get());
  fbc->streams_[0]->temperature_ = 310.0;
  const std::string c0 = d.fluid.components_.back();
  for (auto& [c, y] : fbc->streams_[0]->mole_fraction) y = (c == c0) ? 1.0 : 0.0;
  std::string err;
  ASSERT_TRUE(SetBoundaryConditionKind(d, s, BoundaryLocation::BOTTOM, BoundaryType::PRESSURE, err)) << err;
  auto* p = dynamic_cast<PressureBoundaryCondition*>(StepById(d, s).bottom_boundary_condition_.get());
  ASSERT_TRUE(p);
  EXPECT_EQ(p->location_, BoundaryLocation::BOTTOM);
  EXPECT_EQ(p->direction_, BoundaryDirection::INLET);
  EXPECT_DOUBLE_EQ(p->temperature_.Get(), 310.0);
  EXPECT_DOUBLE_EQ(p->mole_fraction_.at(c0).Get(), 1.0);
  EXPECT_NE(p->endpoint_id_, 0u);

  ASSERT_TRUE(SetBoundaryConditionKind(d, s, BoundaryLocation::BOTTOM, BoundaryType::CLOSED, err)) << err;
  ASSERT_TRUE(SetBoundaryConditionKind(d, s, BoundaryLocation::BOTTOM, BoundaryType::FLOW, err)) << err;
  auto* f = dynamic_cast<FlowBoundaryCondition*>(StepById(d, s).bottom_boundary_condition_.get());
  ASSERT_TRUE(f);
  ASSERT_EQ(f->streams_.size(), 1u);
  EXPECT_DOUBLE_EQ(Sum(f->streams_[0]->mole_fraction), 1.0);
  ASSERT_TRUE(SetBoundaryConditionKind(d, s, BoundaryLocation::TOP, BoundaryType::FLOW, err)) << err;
  ExpectChecksPass(d, "boundaries switched");
}

TEST(DeckEditCascades, SetStreamKindKeepsTheEndpoint) {
  UserData d = Load("Perez2019_A");
  const size_t s = AddStep(d, "St");
  auto* fbc = dynamic_cast<FlowBoundaryCondition*>(StepById(d, s).bottom_boundary_condition_.get());
  const size_t id = fbc->streams_[0]->endpoint_id_;
  fbc->streams_[0]->name_ = "Feed";
  fbc->streams_[0]->temperature_ = 305.0;
  std::string err;
  ASSERT_TRUE(SetStreamKind(d, s, BoundaryLocation::BOTTOM, 0, BoundaryStreamType::PUMP, err)) << err;
  ASSERT_EQ(fbc->streams_[0]->GetType(), BoundaryStreamType::PUMP);
  EXPECT_EQ(fbc->streams_[0]->endpoint_id_, id);
  EXPECT_EQ(fbc->streams_[0]->name_, "Feed");
  EXPECT_DOUBLE_EQ(fbc->streams_[0]->temperature_.Get(), 305.0);
  EXPECT_TRUE(dynamic_cast<PumpBoundaryStream*>(fbc->streams_[0].get())->pump_->use_pump_);
  ASSERT_TRUE(SetStreamKind(d, s, BoundaryLocation::BOTTOM, 0, BoundaryStreamType::VALVE, err)) << err;
  EXPECT_EQ(fbc->streams_[0]->endpoint_id_, id);
  ASSERT_TRUE(SetStreamKind(d, s, BoundaryLocation::BOTTOM, 0, BoundaryStreamType::NONE, err)) << err;
  EXPECT_FALSE(SetStreamKind(d, s, BoundaryLocation::BOTTOM, 3, BoundaryStreamType::NONE, err));
  ExpectChecksPass(d, "stream kinds switched");
}

TEST(DeckEditCascades, FlowAndProfileKinds) {
  UserData d = Load("Perez2019_A");
  const size_t s = AddStep(d, "P");
  auto* fbc = dynamic_cast<FlowBoundaryCondition*>(StepById(d, s).bottom_boundary_condition_.get());
  const size_t stream = fbc->streams_[0]->endpoint_id_;
  const size_t top = dynamic_cast<PressureBoundaryCondition*>(StepById(d, s).top_boundary_condition_.get())->endpoint_id_;
  std::string err;

  // Pressure: constant value -> exponential end -> linear end; use_step_time kept.
  auto* pbc = dynamic_cast<PressureBoundaryCondition*>(StepById(d, s).top_boundary_condition_.get());
  dynamic_cast<ConstantBoundaryValue*>(pbc->pressure_.get())->constant_value_ = 2.5e5;
  ASSERT_TRUE(SetProfileKind(d, top, ProfileTarget::Pressure, BoundaryValueType::EXPONENTIAL, err)) << err;
  auto* ex = dynamic_cast<ExponentialBoundaryValue*>(pbc->pressure_.get());
  ASSERT_TRUE(ex);
  EXPECT_DOUBLE_EQ(ex->end_value_.Get(), 2.5e5);
  ex->use_step_time_ = false;
  ex->step_time_ = 12.0;
  ASSERT_TRUE(SetProfileKind(d, top, ProfileTarget::Pressure, BoundaryValueType::LINEAR, err)) << err;
  auto* lin = dynamic_cast<LinearBoundaryValue*>(pbc->pressure_.get());
  ASSERT_TRUE(lin);
  EXPECT_DOUBLE_EQ(lin->end_value_.Get(), 2.5e5);
  EXPECT_DOUBLE_EQ(lin->ramp_time_.Get(), 12.0);
  EXPECT_FALSE(lin->use_step_time_);
  ASSERT_TRUE(SetProfileKind(d, top, ProfileTarget::Pressure, BoundaryValueType::CONSTANT, err)) << err;
  EXPECT_DOUBLE_EQ(dynamic_cast<ConstantBoundaryValue*>(pbc->pressure_.get())->constant_value_.Get(), 2.5e5);
  EXPECT_FALSE(SetProfileKind(d, top, ProfileTarget::Flow, BoundaryValueType::LINEAR, err));

  // Flow: a linear velocity ramp; switching to a molar flow keeps the ramp
  // kind and time, not the values.
  ASSERT_TRUE(SetProfileKind(d, stream, ProfileTarget::Flow, BoundaryValueType::LINEAR, err)) << err;
  auto* bs = dynamic_cast<BoundaryStream*>(fbc->streams_[0].get());
  auto* vl = dynamic_cast<LinearBoundaryValue*>(bs->boundary_flow_->profile_.get());
  ASSERT_TRUE(vl);
  vl->ramp_time_ = 7.0;
  vl->end_value_ = 0.3;
  ASSERT_TRUE(SetFlowKind(d, stream, FlowType::MOLAR, err)) << err;
  EXPECT_EQ(bs->boundary_flow_->GetType(), FlowType::MOLAR);
  auto* ml = dynamic_cast<LinearBoundaryValue*>(bs->boundary_flow_->profile_.get());
  ASSERT_TRUE(ml);
  EXPECT_DOUBLE_EQ(ml->ramp_time_.Get(), 7.0);
  EXPECT_NE(ml->end_value_.Get(), 0.3);
  EXPECT_FALSE(SetFlowKind(d, top, FlowType::MOLAR, err)) << "a pressure boundary has no flow";

  // Temperature: constant -> linear starts and ends at the constant.
  fbc->streams_[0]->temperature_ = 320.0;
  ASSERT_TRUE(SetProfileKind(d, stream, ProfileTarget::Temperature, BoundaryValueType::LINEAR, err)) << err;
  auto* tl = dynamic_cast<LinearBoundaryValue*>(fbc->streams_[0]->temperature_profile_.get());
  ASSERT_TRUE(tl);
  EXPECT_DOUBLE_EQ(tl->start_value_.Get(), 320.0);
  EXPECT_DOUBLE_EQ(tl->end_value_.Get(), 320.0);
  ASSERT_TRUE(SetProfileKind(d, stream, ProfileTarget::Temperature, BoundaryValueType::CONSTANT, err)) << err;
  EXPECT_EQ(fbc->streams_[0]->temperature_profile_, nullptr);
  ExpectChecksPass(d, "profiles switched");
}

TEST(DeckEditCascades, PumpKinds) {
  UserData d = Load("Perez2019_A");
  const size_t s = AddStep(d, "Pumps");
  auto* pbc = dynamic_cast<PressureBoundaryCondition*>(StepById(d, s).top_boundary_condition_.get());
  const size_t top = pbc->endpoint_id_;
  std::string err;
  ASSERT_TRUE(SetPumpKind(d, top, 2, err)) << err;  // vacuum pump
  ASSERT_TRUE(pbc->pump_ && pbc->pump_->use_pump_);
  EXPECT_EQ(pbc->pump_->GetType(), PumpType::VACUUM_PUMP);
  ASSERT_TRUE(SetPumpEfficiencyKind(d, top, PumpEfficiencyType::VACUUM_PUMP, err)) << err;
  dynamic_cast<VacuumPumpEfficiencyModel*>(pbc->pump_->efficiency_model_.get())->eta_max_ = 0.65;
  pbc->pump_->poly_index_ = 1.3;
  ASSERT_TRUE(SetPumpKind(d, top, 1, err)) << err;  // compressor
  EXPECT_EQ(pbc->pump_->GetType(), PumpType::COMPRESSOR);
  EXPECT_DOUBLE_EQ(pbc->pump_->poly_index_.Get(), 1.3);
  auto* eff = dynamic_cast<VacuumPumpEfficiencyModel*>(pbc->pump_->efficiency_model_.get());
  ASSERT_TRUE(eff) << "the efficiency model moves with the pump";
  EXPECT_DOUBLE_EQ(eff->eta_max_.Get(), 0.65);
  ASSERT_TRUE(SetPumpKind(d, top, 0, err)) << err;
  EXPECT_FALSE(pbc->pump_->use_pump_);

  size_t ps = 0;
  ASSERT_TRUE(AddStream(d, s, BoundaryLocation::BOTTOM, BoundaryStreamType::PUMP, ps, err)) << err;
  EXPECT_FALSE(SetPumpKind(d, ps, 0, err)) << "a pump stream always has a pump";
  ASSERT_TRUE(SetPumpPowerKind(d, ps, PumpPowerType::POLYNOMIAL, err)) << err;
  EXPECT_FALSE(SetPumpPowerKind(d, top, PumpPowerType::POLYNOMIAL, err));
  size_t valve = 0;
  ASSERT_TRUE(AddStream(d, s, BoundaryLocation::BOTTOM, BoundaryStreamType::VALVE, valve, err)) << err;
  EXPECT_FALSE(SetPumpKind(d, valve, 1, err)) << "a valve stream has no pump";
}

// ---- round trip ------------------------------------------------------------------------------------

TEST(DeckEditCascades, RoundTripAfterEveryCascade) {
  UserData d = Load("FiveStepProductPurge");
  std::string err;
  const LayerID layer = AddLayer(d);
  ASSERT_TRUE(MoveLayer(d, layer, -1, err)) << err;
  const size_t s = AddStep(d);
  size_t copy = 0, ep = 0;
  ASSERT_TRUE(DuplicateStep(d, s, copy, err)) << err;
  ASSERT_TRUE(MoveStep(d, copy, -1, err)) << err;
  ASSERT_TRUE(AddStream(d, s, BoundaryLocation::BOTTOM, BoundaryStreamType::VALVE, ep, err)) << err;
  ASSERT_TRUE(SetStreamKind(d, s, BoundaryLocation::BOTTOM, 1, BoundaryStreamType::PUMP, err)) << err;
  ASSERT_TRUE(SetBoundaryConditionKind(d, copy, BoundaryLocation::TOP, BoundaryType::FLOW, err)) << err;
  const size_t top = dynamic_cast<PressureBoundaryCondition*>(StepById(d, s).top_boundary_condition_.get())->endpoint_id_;
  ASSERT_TRUE(SetProfileKind(d, top, ProfileTarget::Pressure, BoundaryValueType::LOG_LINEAR, err)) << err;
  ASSERT_TRUE(SetPumpKind(d, top, 2, err)) << err;
  const AdsorbentID ads = d.adsorbent_directory.GetAll().begin()->first;
  ASSERT_TRUE(SetIsothermKind(d, ads, d.fluid.components_.front(), IsothermType::TOTH, err)) << err;
  ASSERT_TRUE(SetMassTransferKind(d, ads, d.fluid.components_.front(), MassTransferModelType::LDF,
                                  MassTransferCoefficientType::ARRHENIUS, err)) << err;
  ASSERT_TRUE(SetDispersionKind(d, layer, AxialDispersionType::RUTHVEN, err)) << err;
  EnsureIds(d);
  ExpectChecksPass(d, "edited");

  const std::string text = Saved(d);
  UserData back;
  LoadFromString(back, text);
  EXPECT_EQ(Canonical(back), Canonical(d)) << "the reloaded deck equals the edited one";
  ExpectChecksPass(back, "reloaded");
}

}  // namespace adsim
