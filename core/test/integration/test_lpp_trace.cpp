// test_lpp_trace.cpp
//
// An inert trace in a light-product pressurisation deck stays at round-off.
// In Nguyen2023_Exp_VII the Feed 1 top holds the vent and the light-product
// draw, and the draw takes about twice what the bed supplies there; the vent
// back-feeds the rest. When the back-fed gas was a copy of the bed-end gas,
// the pressurisation put back about twice what left the top, and a component
// present only at round-off (water, added to the deck) grew about 1.9 times
// per cycle, to a mole fraction of 6e-5 at the top by cycle 20 and 0.95 by
// cycle 36. The vents carry N2 and the outlet node mixes the back-fed gas
// (FlowBoundaryCondition::Mixing), so the water stays at round-off.
//
// Slow tier (ADSIM_SLOW_TESTS): 20 cycles of the shipped deck.

#include <algorithm>
#include <string>

#include "adsim/simulation.h"
#include "adsim/boundary_conditions/step.h"
#include "adsim/objects/deck_edit.h"

#include <gtest/gtest.h>

namespace adsim {

  TEST(LppTrace, InertTraceStaysAtRoundOffInNguyen2023ExpVII) {
    Simulation sim = Simulation::FromFile("../examples/decks/Nguyen2023_Exp_VII.json");
    std::string err;
    ASSERT_TRUE(AddComponentToDeck(sim.Data(), "Water", err)) << err;
    sim.Data().settings.css.method_ = CSSMethod::NONE;

    double y_max = 0.0;  // water at the top of the bed, as the LPP stream logs it
    RunOptions run;
    run.cycles = 20;
    run.on_message = [](LogLevel, const std::string&) {};
    run.progress = [&](const RunProgress& p) {
      if (p.step != 0) return true;
      for (auto& e : sim.Data().cycle.step_list_[0].LoggedEndpoints()) {
        if (!e.stream || !e.stream->connection_.IsConnected()) continue;
        const auto it = e.data->mole_fraction_boundary_.find("Water");
        if (it != e.data->mole_fraction_boundary_.end())
          for (double y : it->second) y_max = std::max(y_max, y);
      }
      return true;
    };
    const RunResult result = sim.Run(run);
    ASSERT_TRUE(result.Ok()) << result.message;
    EXPECT_LE(y_max, 1e-12) << "the water trace grew to " << y_max;
  }

} // namespace adsim
