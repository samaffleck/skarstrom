// G-04: a failed load must leave the existing model untouched.
//
// cereal deserialises in place, so loading a deck that fails part-way through
// used to leave the target half-overwritten (Bed/Fluid/Cycle from the new file,
// the rest from the old one). The loaders now deserialise into a temporary and
// only move-assign on success. These tests use decks that are syntactically
// valid JSON (so rapidjson parses them up-front) but fail late in the cereal
// walk, after several top-level members have already been read.

// STL includes
#include <filesystem>
#include <fstream>
#include <set>
#include <sstream>
#include <string>

// AdSim includes
#include "adsim/objects/user_data.h"
#include "adsim/solver/cycle.h"
#include "adsim/boundary_conditions/boundary_conditions.h"
#include "adsim/boundary_conditions/boundary_stream.h"

// GTest includes
#include <gtest/gtest.h>


namespace adsim {
namespace {

  const std::filesystem::path kExamples = std::filesystem::current_path() / ".." / "examples" / "decks";

  std::string ReadFile(const std::filesystem::path& p) {
    std::ifstream ifs(p);
    std::stringstream ss;
    ss << ifs.rdbuf();
    return ss.str();
  }

  std::filesystem::path WriteTempFile(const std::string& name, const std::string& content) {
    const auto p = std::filesystem::temp_directory_path() / name;
    std::ofstream ofs(p);
    ofs << content;
    return p;
  }

  // Renames the (unique) key `key` so cereal fails to find it. Everything
  // serialised before it has already been read when the throw happens.
  std::string BreakKey(std::string json, const std::string& key) {
    const std::string quoted = "\"" + key + "\"";
    const auto pos = json.find(quoted);
    EXPECT_NE(pos, std::string::npos) << key;
    json.replace(pos, quoted.size(), "\"" + key + "_Corrupted\"");
    return json;
  }

  // Serialises without SaveToJson's sim_version bump, for before/after equality.
  std::string Snapshot(UserData& data) {
    std::ostringstream oss;
    {
      cereal::JSONOutputArchive ar(oss);
      ar(data);
    }
    return oss.str();
  }

  std::string Snapshot(Cycle& cycle) {
    std::ostringstream oss;
    {
      cereal::JSONOutputArchive ar(oss);
      ar(cycle);
    }
    return oss.str();
  }

  // Every connected endpoint's raw StreamData* must point at a stream/pressure
  // BC owned by `cycle` (i.e. not at a destroyed temporary). Returns the
  // number of connected endpoints checked.
  int ExpectConnectionsPointIntoCycle(Cycle& cycle) {
    std::set<const StreamData*> owned;
    std::vector<Connection*> conns;
    for (Step& step : cycle.step_list_) {
      for (IBoundaryCondition* bc : { step.bottom_boundary_condition_.get(), step.top_boundary_condition_.get() }) {
        if (auto* flow = dynamic_cast<FlowBoundaryCondition*>(bc)) {
          for (auto& s : flow->streams_) {
            owned.insert(&s->stream_data_);
            conns.push_back(&s->connection_);
          }
        } else if (auto* p = dynamic_cast<PressureBoundaryCondition*>(bc)) {
          owned.insert(&p->stream_data_);
          conns.push_back(&p->connection_);
        }
      }
    }
    int connected = 0;
    for (Connection* c : conns) {
      if (!c->IsConnected()) continue;
      ++connected;
      EXPECT_TRUE(owned.count(c->stream_)) << "connection points outside the loaded cycle";
    }
    return connected;
  }

} // namespace

  TEST(UserDataLoadTest, SuccessfulLoadBumpsGenerationAndKeepsConnectionsValid) {
    UserData data;
    data.load_generation_ = 41;
    LoadFromJson(data, (kExamples / "Perez2019_D1.json").string());
    EXPECT_EQ(data.load_generation_, 42u);
    EXPECT_GT(ExpectConnectionsPointIntoCycle(data.cycle), 0);

    // Second load into an already-populated object.
    LoadFromJson(data, (kExamples / "FiveStepProductPurge.json").string());
    EXPECT_EQ(data.load_generation_, 43u);
    EXPECT_GT(ExpectConnectionsPointIntoCycle(data.cycle), 0);
  }

  TEST(UserDataLoadTest, MalformedFileLeavesExistingUserDataUnchanged) {
    UserData data;
    LoadFromJson(data, (kExamples / "Perez2019_A.json").string());
    const std::string before = Snapshot(data);
    const size_t gen_before = data.load_generation_;

    // A different deck that fails after Bed/Fluid/Cycle/... have been read.
    const auto bad = WriteTempFile("adsim_g04_bad_deck.json",
      BreakKey(ReadFile(kExamples / "Perez2019_D1.json"), "System_Information"));

    EXPECT_ANY_THROW(LoadFromJson(data, bad.string()));
    EXPECT_EQ(Snapshot(data), before);
    EXPECT_EQ(data.load_generation_, gen_before);
    std::filesystem::remove(bad);
  }

  TEST(UserDataLoadTest, MissingFileThrowsAndLeavesUserDataUnchanged) {
    UserData data;
    LoadFromJson(data, (kExamples / "Perez2019_A.json").string());
    const std::string before = Snapshot(data);

    EXPECT_ANY_THROW(LoadFromJson(data, (kExamples / "does_not_exist.json").string()));
    EXPECT_EQ(Snapshot(data), before);
  }

  TEST(UserDataLoadTest, MalformedStringLeavesExistingUserDataUnchanged) {
    UserData other;
    LoadFromJson(other, (kExamples / "Perez2019_D1.json").string());
    std::string good;
    SaveToString(other, good);

    UserData data;
    LoadFromJson(data, (kExamples / "Perez2019_A.json").string());
    const std::string before = Snapshot(data);
    const size_t gen_before = data.load_generation_;

    EXPECT_ANY_THROW(LoadFromString(data, BreakKey(good, "System_Information")));
    EXPECT_EQ(Snapshot(data), before);
    EXPECT_EQ(data.load_generation_, gen_before);

    // And the intact string still loads.
    LoadFromString(data, good);
    EXPECT_EQ(data.load_generation_, gen_before + 1);
    EXPECT_GT(ExpectConnectionsPointIntoCycle(data.cycle), 0);
  }

  TEST(UserDataLoadTest, MalformedCycleFileLeavesExistingCycleUnchanged) {
    UserData data;
    LoadFromJson(data, (kExamples / "Perez2019_A.json").string());
    const std::string before = Snapshot(data.cycle);

    UserData other;
    LoadFromJson(other, (kExamples / "Perez2019_D1.json").string());
    std::ostringstream oss;
    {
      cereal::JSONOutputArchive ar(oss);
      ar(other.cycle);
    }
    // Step_List is read after Number_Of_Cycles; break a key inside the steps
    // that appears after earlier steps have been loaded.
    std::string cycle_json = oss.str();
    const auto last = cycle_json.rfind("\"Step_Name\"");
    ASSERT_NE(last, std::string::npos);
    cycle_json.replace(last, std::string("\"Step_Name\"").size(), "\"Step_Name_Corrupted\"");
    const auto bad = WriteTempFile("adsim_g04_bad_cycle.json", cycle_json);

    EXPECT_ANY_THROW(LoadCycleFromJSON(bad.string(), data.cycle));
    EXPECT_EQ(Snapshot(data.cycle), before);
    std::filesystem::remove(bad);
  }

} // namespace adsim
