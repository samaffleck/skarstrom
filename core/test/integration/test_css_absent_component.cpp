// test_css_absent_component.cpp
//
// The per-component mass-balance CSS test with a component that is in the deck
// but never fed. Under the IAST mixing rule such a component takes round-off
// loadings (its infinite-dilution loading times a round-off mole fraction), so
// its boundary flows are round-off-sized but not zero. Its relative balance
// error is then O(1) on about half the cycles, and when the test skipped only
// an inflow below 1e-30 mol, the run passed three cycles in a row only by
// chance (cycle 313 here; in a variant of the deck, never in 800 cycles). The
// test now skips a component whose inflow is negligible against the moles fed
// in the cycle (CyclicSteadyState::NegligibleThroughput), so the deck with the
// absent component stops where the same deck without it stops.
//
// The deck is Nguyen2023_Exp_VII with the light-product draw of Feed 1 removed,
// pure N2 pressurisation, a pure-component dual-site Langmuir isotherm for CO2
// and N2 under IAST, and (in one of the two runs) an inert Water at mole
// fraction 0 everywhere.
//
// Measured at the shipped 80 cells: both decks reach CSS at cycle 274. With
// the 1e-30 mol rule the deck with Water stopped at cycle 313, and its balance
// error differed from the dry deck's on 156 cycles (from cycle 19; 1.0 or up
// to 1e5 where the dry deck had 0.08 down to 1e-3). At 30 cells the absent
// component's flows stay below 1e-30 mol and the old rule passes too, so the
// test keeps the shipped grid.
//
// Slow tier (ADSIM_SLOW_TESTS): two runs of 275 cycles, side by side (about
// 10 min on two cores).

#include <algorithm>
#include <cmath>
#include <cstdio>
#include <map>
#include <string>
#include <thread>
#include <vector>

#include "adsim/simulation.h"
#include "adsim/objects/deck_builder.h"
#include "adsim/objects/deck_edit.h"
#include "adsim/parameters/parameter_registry.h"

#include <gtest/gtest.h>

namespace adsim {

  namespace {

    // The path of the y[component] of the Pressurisation top (its one stream
    // or boundary, whatever its name).
    std::string PressurisationTopY(Simulation& sim, const std::string& component) {
      const std::string prefix = "step[Pressurisation].top", suffix = ".y[" + component + "]";
      for (const auto& p : sim.Parameters().All()) {
        const std::string& a = p.Info().alias;
        if (a.rfind(prefix, 0) == 0 && a.size() >= suffix.size() &&
            a.compare(a.size() - suffix.size(), suffix.size(), suffix) == 0)
          return p.Path();
      }
      ADD_FAILURE() << "no " << prefix << "..." << suffix;
      return {};
    }

    struct Outcome {
      RunResult result;
      std::vector<CSSRecord> css;
      std::vector<std::string> messages;
    };

    Outcome RunVariant(bool with_water) {
      Simulation sim = Simulation::FromFile("../examples/decks/Nguyen2023_Exp_VII.json");
      UserData& data = sim.Data();
      std::string err;

      // Feed 1 top: drop the light-product draw (stream 1), keep the vent.
      size_t feed_id = 0, press_id = 0;
      for (const auto& step : data.cycle.step_list_) {
        if (step.step_name_ == "Feed 1") feed_id = step.id_;
        if (step.step_name_ == "Pressurisation") press_id = step.id_;
      }
      EXPECT_NE(feed_id, 0u);
      EXPECT_NE(press_id, 0u);
      EXPECT_TRUE(RemoveStream(data, feed_id, BoundaryLocation::TOP, 1, err)) << err;
      sim.RefreshParameters();
      sim.SetParameter(PressurisationTopY(sim, "CarbonDioxide"), 0.0);
      sim.SetParameter(PressurisationTopY(sim, "Nitrogen"), 1.0);

      if (with_water) {
        EXPECT_TRUE(AddComponentToDeck(data, "Water", err)) << err;
        sim.RefreshParameters();
        sim.SetParameter(PressurisationTopY(sim, "Water"), 0.0);
      }

      // The same dual-site Langmuir parameters, as pure-component isotherms
      // combined by IAST.
      const AdsorbentID ads = data.adsorbent_directory.GetAll().begin()->first;
      const std::string prefix = "adsorbent[" + std::to_string(ads) + "].isotherm[";
      const std::vector<std::string> components{"CarbonDioxide", "Nitrogen"};
      const std::vector<std::string> fields{"m1", "b1", "q1", "m2", "b2", "q2"};
      std::map<std::string, double> values;
      for (const auto& c : components)
        for (const auto& f : fields) values[c + "." + f] = sim.GetParameter(prefix + c + "]." + f);
      build::SetIsotherm(data, ads, "dual_site_langmuir_c", components);
      build::SetMixingRule(data, ads, "iast");
      sim.RefreshParameters();
      for (const auto& c : components)
        for (const auto& f : fields) sim.SetParameter(prefix + c + "]." + f, values[c + "." + f]);

      sim.SetParameter("cycle.number_of_cycles", 400);

      Outcome out;
      RunOptions run;
      run.on_message = [&](LogLevel, const std::string& m) { out.messages.push_back(m); };
      out.result = sim.Run(run);
      out.css = sim.CSSHistory();
      return out;
    }

  }  // namespace

  TEST(CSSAbsentComponent, IastDeckWithAnAbsentComponentReachesCssWithTheDryDeck) {
    // The two runs are independent: run them side by side.
    Outcome dry;
    std::thread dry_thread([&dry] { dry = RunVariant(false); });
    const Outcome wet = RunVariant(true);
    dry_thread.join();

    ASSERT_TRUE(dry.result.Ok()) << dry.result.message;
    ASSERT_TRUE(dry.result.css_reached) << "the deck without the absent component did not reach CSS";
    ASSERT_TRUE(wet.result.Ok()) << wet.result.message;
    ASSERT_TRUE(wet.result.css_reached) << "the deck with the absent component did not reach CSS";

    // Within a few cycles of the deck without it (the stop then depends only
    // on the state test, which the absent component does not change).
    EXPECT_LE(std::abs(wet.result.css_cycle - dry.result.css_cycle), 5)
        << "dry CSS at cycle " << dry.result.css_cycle << ", with the absent component at cycle "
        << wet.result.css_cycle;

    std::printf("CSS at cycle %d without the absent component, %d with it\n", dry.result.css_cycle,
                wet.result.css_cycle);

    // Cycle by cycle, the balance error is the one without the absent
    // component: it is not tested, so its round-off error does not enter.
    const size_t n = std::min(dry.css.size(), wet.css.size());
    ASSERT_GT(n, 0u);
    for (size_t c = 0; c < n; ++c)
      EXPECT_NEAR(wet.css[c].error, dry.css[c].error, 0.05 * dry.css[c].error + 1e-4)
          << "balance error at cycle " << c << ": " << wet.css[c].error << " with the absent component, "
          << dry.css[c].error << " without";

    // The decision is logged once.
    int logged = 0;
    for (const auto& m : wet.messages)
      if (m.find("component Water: negligible throughput") != std::string::npos) ++logged;
    EXPECT_EQ(logged, 1);
  }

}  // namespace adsim
