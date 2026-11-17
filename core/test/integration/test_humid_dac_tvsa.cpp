// Illustrative humid direct-air-capture TVSA cycle that crosses water
// saturation: examples/decks/Illustrative_Humid_DAC_TVSA.json.
//
// ILLUSTRATIVE ONLY, not a reproduction of a paper. Five steps on the library
// adsorbent APDES-NFC (GAB water, humid Toth CO2): humid air adsorption
// (400 ppm CO2, 49 % RH at 20 degC, 4 h), evacuation to 0.2 bar, heating with
// a 95 degC medium at 0.2 bar (jacket and internal exchanger), closed cooling
// with 20 degC water, feed pressurisation. The feed, column diameter and wall
// heat transfer are those of Illustrative_Humid_DAC_Breakthrough; the step
// times, vacuum and media are round numbers for the example.
//
// The deck is built by examples/illustrative_humid_dac_tvsa.py
// (python/tests/test_examples_build.py keeps the script and the deck equal).
// This test loads it and saves it back in the form a load gives, so the
// shipped file follows the engine's serialisation like the other decks.
//
// The water-rich gas that the heating step leaves in the bed passes its
// saturation pressure at the wall temperature in the cooling step: the jacket
// cools the thin steel wall faster than the packing. The engine does not model
// condensation, and the SUPERSATURATION run-health check must say so.
//
// Measured (Linux, release, 3 cycles of the shipped deck): the warning names
// Water at the wall in step "Cooling" of cycle 3, relative saturation 2.56 at
// the top cell 32 s into the step. Per cycle, the wall measure is 2.13, 2.39
// and 2.56 (always in Cooling) and the gas measure 0.999, 0.98 and 0.96
// (always in Pressurisation, where the feed compresses the water vapour left
// in the bed): the sorbent takes the water back up as the bed cools and is
// compressed, which holds the gas just below saturation. Run to CSS (cycle 12),
// the last cycle gives 1.41 at the wall and 0.91 in the gas.
//
// Fast tier: 3 cycles, CSS off (~2 s).

#include <cmath>
#include <filesystem>

#include <gtest/gtest.h>

#include "../fixtures/atomic_write.h"
#include "adsim/simulation.h"

namespace adsim {

namespace {

const std::filesystem::path kDeck =
    std::filesystem::path("..") / "examples" / "decks" / "Illustrative_Humid_DAC_TVSA.json";

}  // namespace

TEST(HumidDacTvsa, CoolingCrossesWaterSaturationAtTheWall) {
  const std::filesystem::path deck = std::filesystem::current_path() / kDeck;
  ASSERT_TRUE(std::filesystem::exists(deck)) << deck;
  {
    Simulation shipped = Simulation::FromFile(deck.string());
    test::SaveDeckAtomically(shipped.Data(), deck);
  }

  Simulation sim = Simulation::FromFile(deck.string());
  EXPECT_TRUE(sim.Data().settings.warnings_.supersaturation);  // on in the shipped deck
  sim.SetParameter("cycle.number_of_cycles", 3);
  sim.SetParameter("css.method", 0);  // none: three cycles only
  const RunResult r = sim.Run();
  ASSERT_TRUE(r.Ok()) << r.message;
  EXPECT_EQ(r.cycles_run, 3);
  EXPECT_FALSE(r.healthy);

  // The only run-health warning is the supersaturation of water at the wall,
  // in the cooling step.
  ASSERT_EQ(r.health_warnings.size(), 1u);
  const RunHealthWarning& w = r.health_warnings.front();
  std::printf("  %s\n", w.message.c_str());
  EXPECT_EQ(w.code, "SUPERSATURATION");
  EXPECT_EQ(w.component, "Water");
  EXPECT_EQ(w.step_name, "Cooling");
  EXPECT_DOUBLE_EQ(w.threshold, 1.0);
  EXPECT_GT(w.value, 1.5);
  EXPECT_NE(w.message.find("at the wall temperature"), std::string::npos) << w.message;

  // Every cycle crosses saturation at the wall. The gas comes close (0.96-1.00:
  // the feed compresses the water vapour left by the cooling step, and the
  // sorbent takes it back up), so its measure is checked for a lower bound only.
  const RunHealthReport& health = sim.RunHealth();
  ASSERT_EQ(health.cycles.size(), 3u);
  for (const RunHealthCycle& c : health.cycles) {
    std::printf("  cycle %d: relative saturation %.3f in the gas, %.3f at the wall\n", c.cycle + 1,
                c.supersaturation_gas, c.supersaturation_wall);
    EXPECT_GT(c.supersaturation_wall, 1.5) << "cycle " << c.cycle;
    EXPECT_GT(c.supersaturation_gas, 0.9) << "cycle " << c.cycle;
  }
}

}  // namespace adsim
