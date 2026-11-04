/*
 * Skarstrom - pressure and temperature swing adsorption simulator
 * Copyright (c) 2024-2026 Sam Affleck and contributors
 *
 * Licensed under the terms in the LICENSE file at the repository root.
 */

// Dynamic tanks in the deck: serialisation round trip, old decks, endpoint
// IDs shared with the bed endpoints, the builder and the registry paths.

#include "../fixtures/user_data_test.h"
#include "adsim/objects/deck_builder.h"
#include "adsim/parameters/parameter_registry.h"

#include <regex>

namespace adsim {

  namespace {
    // A small deck with one tank of every thermal and stream kind in use, a
    // stream port and a pressure-boundary port.
    void BuildTankDeck(UserData& data, size_t& tank, size_t& feed_ep, size_t& pressure_ep) {
      namespace b = build;
      b::AddComponent(data, "Nitrogen");
      b::AddComponent(data, "Oxygen");
      const auto ads = b::AddAdsorbent(data, "Packing");
      const auto layer = b::AddLayer(data, "Layer");
      b::SetLayerComposition(data, layer, {{ads, 1.0}});
      data.ic.mole_fraction["Nitrogen"] = 0.79;
      data.ic.mole_fraction["Oxygen"] = 0.21;
      const size_t feed = b::AddStep(data, "Feed");
      b::SetFlowBoundary(data, feed, "bottom", "inlet");
      feed_ep = b::AddStream(data, feed, "bottom", "Feed gas", "molar_flow", "constant");
      b::SetFlowBoundary(data, feed, "top", "outlet");
      const size_t product_ep = b::AddValveStream(data, feed, "top", "To tank", "constant");
      const size_t purge = b::AddStep(data, "Purge");
      pressure_ep = b::SetPressureBoundary(data, purge, "top", "inlet", "constant");
      b::SetClosed(data, purge, "bottom");

      tank = b::AddTank(data, "Product tank");
      b::SetTankThermal(data, tank, "wall");
      b::AddTankStream(data, tank, "make_up", "Make-up");
      b::AddTankStream(data, tank, "valve", "Vent");
      b::AddTankStream(data, tank, "regulator", "Regulator");
      const size_t hold = b::AddTankStream(data, tank, "hold", "Product");
      (void)hold;
      b::SetTankPort(data, product_ep, tank);
      b::SetTankPort(data, pressure_ep, tank);
      Tank* t = data.tanks.Find(tank);
      t->volume_ = 0.0314;
      t->initial_pressure_ = 11.12e5;
      t->initial_mole_fraction_["Nitrogen"] = 0.999;
      t->initial_mole_fraction_["Oxygen"] = 0.001;
      t->h_in_ = 12.5;
    }

    std::string WithoutVersion(const std::string& json) {
      return std::regex_replace(json, std::regex("\"Version\": [0-9]+"), "\"Version\": 0");
    }
  }  // namespace

  TEST(TankDeck, RoundTripKeepsTanksStreamsAndPorts) {
    UserData data;
    size_t tank = 0, feed_ep = 0, pressure_ep = 0;
    BuildTankDeck(data, tank, feed_ep, pressure_ep);
    std::string first;
    SaveToString(data, first);

    UserData loaded;
    LoadFromString(loaded, first);
    ASSERT_EQ(loaded.tanks.Size(), 1u);
    const Tank* t = loaded.tanks.Find(tank);
    ASSERT_NE(t, nullptr);
    EXPECT_EQ(t->name_, "Product tank");
    EXPECT_EQ(t->thermal_, TankThermalMode::WALL);
    EXPECT_DOUBLE_EQ(t->volume_.Get(), 0.0314);
    EXPECT_DOUBLE_EQ(t->initial_pressure_.Get(), 11.12e5);
    EXPECT_DOUBLE_EQ(t->h_in_.Get(), 12.5);
    EXPECT_DOUBLE_EQ(t->initial_mole_fraction_.at("Oxygen").Get(), 0.001);
    ASSERT_EQ(t->streams_.size(), 4u);
    EXPECT_EQ(t->streams_[0]->GetKind(), TankStreamKind::MAKE_UP);
    EXPECT_EQ(t->streams_[1]->GetKind(), TankStreamKind::VALVE);
    EXPECT_EQ(t->streams_[2]->GetKind(), TankStreamKind::REGULATOR);
    EXPECT_EQ(t->streams_[3]->GetKind(), TankStreamKind::HOLD);
    EXPECT_EQ(t->streams_[3]->name_, "Product");

    // The ports survive, the feed stream has none.
    auto feed = loaded.cycle.FindEndpoint(feed_ep);
    auto top = loaded.cycle.FindEndpoint(pressure_ep);
    ASSERT_TRUE(feed && top);
    EXPECT_EQ(feed->stream->tank_port_.tank_id_, 0u);
    EXPECT_EQ(top->pressure->tank_port_.tank_id_, tank);
    size_t ported_streams = 0;
    for (auto& h : loaded.cycle.Endpoints())
      if (h.stream && h.stream->tank_port_.tank_id_ == tank) ++ported_streams;
    EXPECT_EQ(ported_streams, 1u);

    // Saving the loaded deck gives the same archive (but the save counter).
    std::string second;
    SaveToString(loaded, second);
    EXPECT_EQ(WithoutVersion(first), WithoutVersion(second));
    std::string error;
    EXPECT_TRUE(loaded.PerformChecks(error)) << error;
  }

  TEST(TankDeck, TankStreamsShareTheEndpointIdSpace) {
    UserData data;
    size_t tank = 0, feed_ep = 0, pressure_ep = 0;
    BuildTankDeck(data, tank, feed_ep, pressure_ep);
    std::set<size_t> ids;
    for (auto& h : data.cycle.Endpoints()) EXPECT_TRUE(ids.insert(h.id).second) << "repeated endpoint ID " << h.id;
    for (auto& s : data.tanks.Find(tank)->streams_)
      EXPECT_TRUE(ids.insert(s->endpoint_id_).second) << "tank stream repeats endpoint ID " << s->endpoint_id_;
    // A stream added to the cycle afterwards does not take a tank stream's ID.
    const size_t step = build::AddStep(data, "Extra");
    build::SetFlowBoundary(data, step, "bottom", "inlet");
    const size_t extra = build::AddStream(data, step, "bottom", "Extra gas", "molar_flow", "constant");
    EXPECT_TRUE(ids.insert(extra).second) << "new endpoint took a tank stream's ID " << extra;

    // The registry lists the tank (tank[ID]) and its streams (endpoint[ID]).
    ParameterRegistry reg(data);
    EXPECT_NE(reg.Find("tank[" + std::to_string(tank) + "].volume"), nullptr);
    EXPECT_NE(reg.Find("tank[Product tank].initial_pressure"), nullptr);
    EXPECT_NE(reg.Find("tank[" + std::to_string(tank) + "].y[Oxygen]"), nullptr);
    EXPECT_NE(reg.Find("tank[" + std::to_string(tank) + "].h_in"), nullptr);
    const size_t hold = data.tanks.Find(tank)->streams_[3]->endpoint_id_;
    Parameter* band = reg.Find("endpoint[" + std::to_string(hold) + "].band");
    ASSERT_NE(band, nullptr);
    EXPECT_DOUBLE_EQ(band->GetValue(), TankHold::kDefaultBand);
    EXPECT_NE(reg.Find("tank[Product tank].stream[Product].set_pressure"), nullptr);
  }

  TEST(TankDeck, RemovingATankClearsItsPorts) {
    UserData data;
    size_t tank = 0, feed_ep = 0, pressure_ep = 0;
    BuildTankDeck(data, tank, feed_ep, pressure_ep);
    build::RemoveTank(data, tank);
    EXPECT_TRUE(data.tanks.Empty());
    for (auto& h : data.cycle.Endpoints()) {
      const TankPort& port = h.stream ? h.stream->tank_port_ : h.pressure->tank_port_;
      EXPECT_FALSE(port.IsSet()) << h.log_name;
    }
  }

  TEST(TankDeck, ChecksRejectAPortWithAConnectionOrAnUnknownTank) {
    UserData data;
    size_t tank = 0, feed_ep = 0, pressure_ep = 0;
    BuildTankDeck(data, tank, feed_ep, pressure_ep);
    std::string error;
    // A port naming a tank that is not in the deck.
    data.cycle.FindEndpoint(feed_ep)->stream->tank_port_.tank_id_ = 99;
    EXPECT_FALSE(data.PerformChecks(error));
    EXPECT_NE(error.find("does not exist"), std::string::npos) << error;
    data.cycle.FindEndpoint(feed_ep)->stream->tank_port_.tank_id_ = 0;
    // More than one bed sharing a tank is not implemented yet.
    data.tanks.Find(tank)->beds_ = 2;
    error.clear();
    EXPECT_FALSE(data.PerformChecks(error));
    EXPECT_NE(error.find("beds"), std::string::npos) << error;
    // The builder refuses a port on a connected endpoint.
    data.tanks.Find(tank)->beds_ = 1;
    build::SetTankPort(data, pressure_ep, 0);
    build::Connect(data, feed_ep, pressure_ep);
    EXPECT_THROW(build::SetTankPort(data, pressure_ep, tank), std::invalid_argument);
  }

  // Every shipped deck predates tanks (UserData < 5, IBoundaryStream and
  // PressureBoundaryCondition < 6): it loads with no tanks and no ports.
  TEST(TankDeck, OldDecksLoadWithoutTanks) {
    for (const char* name : {"Ward2024", "Casas2012_Fig4_50_50_15bar_25C", "Perez2019_D1"}) {
      const std::string path = std::string("../examples/decks/") + name + ".json";
      if (!std::filesystem::exists(path)) continue;
      UserData data;
      LoadFromJson(data, path);
      EXPECT_TRUE(data.tanks.Empty()) << name;
      for (auto& h : data.cycle.Endpoints()) {
        const TankPort& port = h.stream ? h.stream->tank_port_ : h.pressure->tank_port_;
        EXPECT_FALSE(port.IsSet()) << name << ": " << h.log_name;
      }
      std::string error;
      EXPECT_TRUE(data.PerformChecks(error)) << name << ": " << error;
    }
    EXPECT_TRUE(std::filesystem::exists("../examples/decks/Ward2024.json"));
  }

}  // namespace adsim
