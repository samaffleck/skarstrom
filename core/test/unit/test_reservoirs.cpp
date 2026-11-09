/*
 * Skarstrom - pressure and temperature swing adsorption simulator
 * Copyright (c) 2024-2026 Sam Affleck and contributors
 *
 * Licensed under the terms in the LICENSE file at the repository root.
 */
// Named sources and sinks (adsim/boundary_conditions/reservoir.h): the edit
// cascades, how an attached endpoint takes the reservoir's values, the
// registry and deck-view entries, the deck round trip, and that a deck run
// through a source gives the run of the same values typed on the endpoint.

#include <cmath>
#include <filesystem>
#include <fstream>
#include <sstream>
#include <string>

#include <gtest/gtest.h>

#include "adsim/objects/deck_builder.h"
#include "adsim/objects/deck_edit.h"
#include "adsim/objects/deck_io.h"
#include "adsim/objects/deck_view.h"
#include "adsim/objects/user_data.h"
#include "adsim/parameters/parameter_registry.h"
#include "adsim/simulation.h"

namespace adsim {
namespace {

const std::filesystem::path kExamples = std::filesystem::current_path() / ".." / "examples" / "decks";

std::string ReadText(const std::filesystem::path& p) {
  std::ifstream in(p, std::ios::binary);
  std::ostringstream ss;
  ss << in.rdbuf();
  return ss.str();
}

// Perez2019_A: Feed (flow inlet at the bottom, pressure outlet at the top),
// Blowdown, Evacuation, Pressurisation (pressure inlet at the bottom).
UserData LoadPerez() {
  UserData data;
  LoadDeckFlexible(data, ReadText(kExamples / "Perez2019_A.json"));
  return data;
}

Cycle::EndpointHandle EndpointOf(UserData& data, const std::string& step, BoundaryLocation side) {
  for (const auto& h : data.cycle.Endpoints())
    if (h.step->step_name_ == step && h.location == side) return h;
  ADD_FAILURE() << "no endpoint at that end of step " << step;
  return {};
}

bool HasPath(UserData& data, const std::string& path) {
  ParameterRegistry registry(data);
  return registry.Find(path) != nullptr;
}

}  // namespace

TEST(Reservoirs, AddSeedsNamesIdsAndMoleFractions) {
  UserData data = LoadPerez();
  const size_t a = AddReservoir(data, ReservoirKind::SOURCE);
  const size_t b = AddReservoir(data, ReservoirKind::SINK, "Ambient");
  const size_t c = AddReservoir(data, ReservoirKind::SINK, "Ambient");
  EXPECT_NE(a, 0u);
  EXPECT_NE(a, b);
  EXPECT_NE(b, c);
  EXPECT_EQ(data.cycle.FindReservoir(a)->name_, "Source 1");
  EXPECT_EQ(data.cycle.FindReservoir(b)->name_, "Ambient");
  EXPECT_NE(data.cycle.FindReservoir(c)->name_, "Ambient");  // names stay unique
  const Reservoir& r = *data.cycle.FindReservoir(a);
  EXPECT_EQ(r.mole_fraction_.size(), data.fluid.components_.size());
  double sum = 0.0;
  for (const auto& [name, y] : r.mole_fraction_) sum += y.Get();
  EXPECT_DOUBLE_EQ(sum, 1.0);

  std::string err;
  EXPECT_FALSE(RenameReservoir(data, a, "Ambient", err));
  EXPECT_FALSE(RenameReservoir(data, a, "", err));
  EXPECT_TRUE(RenameReservoir(data, a, "Feed", err)) << err;
  EXPECT_FALSE(RenameReservoir(data, 999, "x", err));
}

TEST(Reservoirs, AttachGivesTheEndpointTheReservoirValues) {
  UserData data = LoadPerez();
  const auto feed = EndpointOf(data, "Feed", BoundaryLocation::BOTTOM);
  const auto press = EndpointOf(data, "Pressurisation", BoundaryLocation::BOTTOM);
  ASSERT_TRUE(feed.stream);
  ASSERT_TRUE(press.pressure);

  const size_t source = AddReservoir(data, ReservoirKind::SOURCE, "Feed gas");
  Reservoir& r = *data.cycle.FindReservoir(source);
  r.temperature_ = 310.0;
  r.pressure_ = 1.7e5;
  for (auto& [name, y] : r.mole_fraction_) y = 0.0;
  r.mole_fraction_.begin()->second = 1.0;
  const std::string pure = r.mole_fraction_.begin()->first;

  std::string err;
  ASSERT_TRUE(AttachReservoir(data, feed.id, source, err)) << err;
  ASSERT_TRUE(AttachReservoir(data, press.id, source, err)) << err;
  EXPECT_DOUBLE_EQ(feed.stream->temperature_.Get(), 310.0);
  EXPECT_DOUBLE_EQ(feed.stream->mole_fraction.at(pure).Get(), 1.0);
  EXPECT_DOUBLE_EQ(press.pressure->temperature_.Get(), 310.0);
  EXPECT_DOUBLE_EQ(press.pressure->mole_fraction_.at(pure).Get(), 1.0);
  // The pressure boundary ramps to the source's pressure.
  ParameterRegistry before(data);
  bool found_target = false;
  for (const auto& p : before.All())
    if (p.Path().rfind("endpoint[" + std::to_string(press.id) + "].pressure", 0) == 0 &&
        (p.Path().size() >= 4 && p.Path().compare(p.Path().size() - 4, 4, "_end") == 0))
      found_target = true;
  EXPECT_FALSE(found_target) << "an attached endpoint must not list its pressure target";

  // A later change of the reservoir reaches both endpoints on evaluation.
  r.temperature_ = 320.0;
  data.EvaluateAllExpressions();
  EXPECT_DOUBLE_EQ(feed.stream->temperature_.Get(), 320.0);
  EXPECT_DOUBLE_EQ(press.pressure->temperature_.Get(), 320.0);

  // The registry lists the values once, on the reservoir.
  const std::string rp = "reservoir[" + std::to_string(source) + "]";
  EXPECT_TRUE(HasPath(data, rp + ".temperature"));
  EXPECT_TRUE(HasPath(data, rp + ".pressure"));
  EXPECT_TRUE(HasPath(data, rp + ".y[" + pure + "]"));
  EXPECT_TRUE(HasPath(data, "reservoir[Feed gas].temperature"));
  EXPECT_FALSE(HasPath(data, "endpoint[" + std::to_string(feed.id) + "].temperature"));
  EXPECT_FALSE(HasPath(data, "endpoint[" + std::to_string(feed.id) + "].y[" + pure + "]"));
  // The flow stays on the endpoint.
  bool has_flow = false;
  {
    ParameterRegistry registry(data);
    for (const auto& p : registry.All())
      if (p.Path().rfind("endpoint[" + std::to_string(feed.id) + "].", 0) == 0) has_flow = true;
  }
  EXPECT_TRUE(has_flow);

  // Detaching keeps the values as the endpoint's own and lists them again.
  ASSERT_TRUE(AttachReservoir(data, feed.id, 0, err)) << err;
  EXPECT_DOUBLE_EQ(feed.stream->temperature_.Get(), 320.0);
  EXPECT_TRUE(HasPath(data, "endpoint[" + std::to_string(feed.id) + "].temperature"));
}

TEST(Reservoirs, KindsMustAgreeWithTheDirection) {
  UserData data = LoadPerez();
  const auto feed = EndpointOf(data, "Feed", BoundaryLocation::BOTTOM);   // inlet
  const auto vent = EndpointOf(data, "Feed", BoundaryLocation::TOP);      // outlet
  const size_t source = AddReservoir(data, ReservoirKind::SOURCE);
  const size_t sink = AddReservoir(data, ReservoirKind::SINK);
  std::string err;
  EXPECT_FALSE(AttachReservoir(data, feed.id, sink, err));
  EXPECT_NE(err.find("inlet"), std::string::npos);
  EXPECT_FALSE(AttachReservoir(data, vent.id, source, err));
  EXPECT_FALSE(AttachReservoir(data, feed.id, 999, err));
  EXPECT_FALSE(AttachReservoir(data, 999, source, err));
  EXPECT_TRUE(AttachReservoir(data, vent.id, sink, err)) << err;
  EXPECT_TRUE(AttachReservoir(data, feed.id, source, err)) << err;
  // A sink is a pressure only: an outlet is zero-gradient in temperature and composition.
  EXPECT_TRUE(HasPath(data, "reservoir[" + std::to_string(sink) + "].pressure"));
  EXPECT_FALSE(HasPath(data, "reservoir[" + std::to_string(sink) + "].temperature"));
  EXPECT_TRUE(HasPath(data, "reservoir[" + std::to_string(source) + "].temperature"));

  // A direction changed behind the cascade is caught by the checks and
  // repaired by DetachMismatchedReservoirs.
  feed.boundary->direction_ = BoundaryDirection::OUTLET;
  vent.boundary->direction_ = BoundaryDirection::INLET;
  std::string msg;
  EXPECT_FALSE(data.cycle.PerformChecks(data.fluid, msg));
  EXPECT_NE(msg.find("attaches to a source"), std::string::npos) << msg;
  EXPECT_EQ(DetachMismatchedReservoirs(data), 2u);
  EXPECT_EQ(Cycle::ReservoirIdOf(feed), 0u);
  EXPECT_EQ(Cycle::ReservoirIdOf(vent), 0u);
}

TEST(Reservoirs, FromEndpointTakesItsValuesAndAttachesIt) {
  UserData data = LoadPerez();
  const auto feed = EndpointOf(data, "Feed", BoundaryLocation::BOTTOM);
  const auto vent = EndpointOf(data, "Feed", BoundaryLocation::TOP);
  const double t = feed.stream->temperature_.Get();
  const ExpressionMap y = feed.stream->mole_fraction;
  size_t source = 0, sink = 0;
  std::string err;
  ASSERT_TRUE(AddReservoirFromEndpoint(data, feed.id, "Feed gas", source, err)) << err;
  ASSERT_TRUE(AddReservoirFromEndpoint(data, vent.id, "", sink, err)) << err;
  const Reservoir& s = *data.cycle.FindReservoir(source);
  const Reservoir& k = *data.cycle.FindReservoir(sink);
  EXPECT_EQ(s.kind_, ReservoirKind::SOURCE);
  EXPECT_EQ(k.kind_, ReservoirKind::SINK);
  EXPECT_EQ(k.name_, "Sink 1");
  EXPECT_DOUBLE_EQ(s.temperature_.Get(), t);
  for (const auto& [name, yi] : y) EXPECT_DOUBLE_EQ(s.mole_fraction_.at(name).Get(), yi.Get());
  EXPECT_EQ(Cycle::ReservoirIdOf(feed), source);
  EXPECT_EQ(Cycle::ReservoirIdOf(vent), sink);
  // The set-flow feed has no pressure: the source takes the pressure the
  // step's outlet holds, which is the sink's.
  EXPECT_DOUBLE_EQ(s.pressure_.Get(), k.pressure_.Get());
  EXPECT_FALSE(AddReservoirFromEndpoint(data, 999, "", source, err));
}

TEST(Reservoirs, RemoveDetachesAndKeepsTheValues) {
  UserData data = LoadPerez();
  const auto feed = EndpointOf(data, "Feed", BoundaryLocation::BOTTOM);
  const size_t source = AddReservoir(data, ReservoirKind::SOURCE);
  data.cycle.FindReservoir(source)->temperature_ = 333.0;
  std::string err;
  ASSERT_TRUE(AttachReservoir(data, feed.id, source, err)) << err;
  ASSERT_TRUE(RemoveReservoir(data, source, err)) << err;
  EXPECT_TRUE(data.cycle.reservoirs_.empty());
  EXPECT_EQ(Cycle::ReservoirIdOf(feed), 0u);
  EXPECT_DOUBLE_EQ(feed.stream->temperature_.Get(), 333.0);
  EXPECT_FALSE(RemoveReservoir(data, source, err));
}

TEST(Reservoirs, ComponentsStepsAndKindSwitchesKeepThemConsistent) {
  UserData data = LoadPerez();
  const auto feed = EndpointOf(data, "Feed", BoundaryLocation::BOTTOM);
  size_t source = 0;
  std::string err;
  ASSERT_TRUE(AddReservoirFromEndpoint(data, feed.id, "Feed gas", source, err)) << err;

  // A new component gets an entry on the reservoir, a removed one loses it.
  ASSERT_TRUE(AddComponentToDeck(data, "Argon", err)) << err;
  EXPECT_EQ(data.cycle.FindReservoir(source)->mole_fraction_.count("Argon"), 1u);
  RemoveComponentFromDeck(data, "Argon");
  EXPECT_EQ(data.cycle.FindReservoir(source)->mole_fraction_.count("Argon"), 0u);

  // A copied step stays attached to the same source: two feed steps, one source.
  size_t copy = 0;
  ASSERT_TRUE(DuplicateStep(data, feed.step->id_, copy, err)) << err;
  const Step* dup = data.cycle.FindStepById(copy);
  ASSERT_NE(dup, nullptr);
  const auto dup_feed = EndpointOf(data, dup->step_name_, BoundaryLocation::BOTTOM);
  EXPECT_EQ(Cycle::ReservoirIdOf(dup_feed), source);

  // Switching the stream to a valve keeps the attachment, and the valve's
  // back pressure becomes the source's pressure.
  data.cycle.FindReservoir(source)->pressure_ = 2.5e5;
  ASSERT_TRUE(SetStreamKind(data, copy, BoundaryLocation::BOTTOM, 0, BoundaryStreamType::VALVE, err)) << err;
  const auto valve = EndpointOf(data, dup->step_name_, BoundaryLocation::BOTTOM);
  EXPECT_EQ(Cycle::ReservoirIdOf(valve), source);
  auto* vs = dynamic_cast<ValveBoundaryStream*>(valve.stream);
  ASSERT_NE(vs, nullptr);
  auto* constant = dynamic_cast<ConstantBoundaryValue*>(vs->boundary_value_.get());
  ASSERT_NE(constant, nullptr);
  EXPECT_DOUBLE_EQ(constant->constant_value_.Get(), 2.5e5);

  // A reservoir has one temperature: no ramp while attached.
  EXPECT_FALSE(SetProfileKind(data, valve.id, ProfileTarget::Temperature, BoundaryValueType::LINEAR, err));
}

TEST(Reservoirs, DeckRoundTripAndView) {
  UserData data = LoadPerez();
  const auto feed = EndpointOf(data, "Feed", BoundaryLocation::BOTTOM);
  const auto vent = EndpointOf(data, "Feed", BoundaryLocation::TOP);
  size_t source = 0, sink = 0;
  std::string err;
  ASSERT_TRUE(AddReservoirFromEndpoint(data, feed.id, "Feed gas", source, err)) << err;
  ASSERT_TRUE(AddReservoirFromEndpoint(data, vent.id, "Ambient", sink, err)) << err;

  std::string text;
  SaveToString(data, text);
  UserData back;
  LoadDeckFlexible(back, text);
  ASSERT_EQ(back.cycle.reservoirs_.size(), 2u);
  EXPECT_EQ(back.cycle.reservoirs_[0].name_, "Feed gas");
  EXPECT_EQ(back.cycle.reservoirs_[0].kind_, ReservoirKind::SOURCE);
  EXPECT_EQ(back.cycle.reservoirs_[1].kind_, ReservoirKind::SINK);
  EXPECT_EQ(Cycle::ReservoirIdOf(EndpointOf(back, "Feed", BoundaryLocation::BOTTOM)), source);
  EXPECT_EQ(Cycle::ReservoirIdOf(EndpointOf(back, "Feed", BoundaryLocation::TOP)), sink);

  const std::string view = DeckView(back);
  EXPECT_NE(view.find("\"reservoirs\":[{\"id\":" + std::to_string(source) + ",\"name\":\"Feed gas\",\"kind\":\"source\""),
            std::string::npos)
      << view.substr(0, 400);
  EXPECT_NE(view.find("\"reservoir\":" + std::to_string(source)), std::string::npos);
  EXPECT_NE(view.find("reservoir[" + std::to_string(source) + "].temperature"), std::string::npos);
  EXPECT_NE(DeckViewSchema().find("\"Reservoir\""), std::string::npos);
}

TEST(Reservoirs, BuilderAddsAndAttaches) {
  UserData data = LoadPerez();
  const auto feed = EndpointOf(data, "Feed", BoundaryLocation::BOTTOM);
  EXPECT_EQ(build::ReservoirKinds(), (std::vector<std::string>{"source", "sink"}));
  const size_t source = build::AddReservoir(data, "source", "Feed gas");
  const size_t sink = build::AddReservoir(data, "sink", "Ambient");
  EXPECT_THROW(build::AddReservoir(data, "tank", "x"), std::invalid_argument);
  build::AttachReservoir(data, feed.id, source);
  EXPECT_EQ(Cycle::ReservoirIdOf(feed), source);
  EXPECT_THROW(build::AttachReservoir(data, feed.id, sink), std::invalid_argument);
  build::AttachReservoir(data, feed.id, 0);
  EXPECT_EQ(Cycle::ReservoirIdOf(feed), 0u);
}

// The same values, typed on the endpoints or held by a source and a sink,
// give the same run.
TEST(Reservoirs, ARunThroughASourceEqualsTheRunOfItsValues) {
  auto run = [](bool with_reservoirs) {
    Simulation sim = Simulation::FromJson(ReadText(kExamples / "Perez2019_A.json"));
    if (with_reservoirs) {
      UserData& data = sim.Data();
      size_t id = 0;
      std::string err;
      for (const auto& h : data.cycle.Endpoints())
        EXPECT_TRUE(AddReservoirFromEndpoint(data, h.id, "", id, err)) << err;
      sim.RefreshParameters();
    }
    RunOptions options;
    options.cycles = 2;
    const RunResult result = sim.Run(options);
    EXPECT_TRUE(result.Ok()) << result.message;
    return sim.KPIs(-1);
  };
  const auto plain = run(false);
  const auto attached = run(true);
  ASSERT_EQ(plain.size(), attached.size());
  ASSERT_FALSE(plain.empty());
  for (size_t i = 0; i < plain.size(); ++i) {
    EXPECT_DOUBLE_EQ(plain[i].purity, attached[i].purity);
    EXPECT_DOUBLE_EQ(plain[i].recovery, attached[i].recovery);
  }
}

}  // namespace adsim
