// test_css_start_independence.cpp
//
// Cyclic steady state does not depend on where a run starts. Each deck is run
// to CSS (its own stop rule, no acceleration) from:
//   * the deck as shipped (the baseline);
//   * a clean bed: the gas is the weakly adsorbed component alone, so only
//     it is adsorbed at the start (loadings start in equilibrium with the
//     initial gas);
//   * a bed saturated with the feed: the feed gas, loadings in equilibrium;
//   * a hot and a cold bed: the state after one cycle from the shipped start
//     with every temperature (fluid, solid, wall) raised or lowered by 40 K,
//     given through the warm-start API. Not by editing the deck's initial
//     temperature: with the default CoolProp gas heat capacities that also
//     sets the heat capacities for the whole run (Fluid::Initialise), so it
//     changes the model, not only the start (issue #145);
//   * a neighbouring deck's CSS: the CSS of a deck with the same cycle at
//     another operating point, given through the warm-start API (decks with
//     no such neighbour skip this start).
// A start that coincides with the shipped initial conditions is not run again.
//
// Every start must reach the baseline's CSS within c * eps, eps being the
// stop's state tolerance, in the full end-of-cycle state, the endpoint
// histories and the KPIs (css_state_compare.h), with the same bound for every
// deck and start as CSSReference (test_css_reference.cpp): c <= 100. Two stops
// approaching the CSS from different sides can be up to the sum of their own
// distances from it apart; CSSReference measured those at most 37 eps (Perez
// 2019 E), so two stops at most 74 eps apart; measured here: at most 74
// (Perez 2019 E, clean and hot against shipped). The measured c of each start is
// printed and tabulated in docs/cyclic-steady-state.md. A disagreement beyond
// the bound would mean several steady states or a stop rule that is too
// loose; it is reported, never hidden by a looser bound.
//
// The decks are those of CSSReference (one or more per family and cycle
// type). Slow tier: about 2-10 minutes per deck.

#include <cstdio>
#include <functional>
#include <map>
#include <string>

#include <gtest/gtest.h>

#include "adsim/simulation.h"
#include "css_state_compare.h"

namespace adsim {

  using namespace adsim::validation;

  namespace {

    constexpr double kBound = 100.0;     // as CSSReference
    constexpr int kMaxCycles = 1000;     // cap on any run to CSS
    constexpr double kTemperatureShift = 40.0;  // [K] hot and cold starts

    using Composition = std::map<std::string, double>;

    struct DeckSpec {
      const char* deck;
      Composition clean;      // gas of a clean bed
      Composition feed;       // feed gas
      const char* neighbour;  // deck with the same cycle at another operating point, or nullptr
    };

    const Composition kCleanN2 = {{"Nitrogen", 1.0}, {"CarbonDioxide", 0.0}};
    const Composition kFlueGas = {{"Nitrogen", 0.85}, {"CarbonDioxide", 0.15}};

    struct Result {
      CycleState state;
      std::vector<KPIRecord> kpis;
      int css_cycle = -1;
      double eps = 0.0;
    };

    std::string DeckPath(const std::string& deck) { return "../examples/decks/" + deck + ".json"; }

    // A run of `deck` to CSS, after `edit` changes its initial conditions,
    // or from `from` (warm start) when given.
    Result RunToCss(const std::string& deck, const std::function<void(UserData&)>& edit,
                    const CycleState* from = nullptr) {
      Simulation sim = Simulation::FromFile(DeckPath(deck));
      auto& data = sim.Data();
      if (data.settings.css.method_ == CSSMethod::NONE) data.settings.css.method_ = CSSMethod::MASS_BALANCE_COMPONENT;
      data.settings.acc_css.method_ = AccelerationMethod::NONE;
      if (edit) edit(data);
      RunOptions options;
      options.cycles = kMaxCycles;
      options.initial_state = from;
      options.on_message = [](LogLevel, const std::string&) {};
      const RunResult run = sim.Run(options);
      EXPECT_TRUE(run.Ok()) << deck << ": " << run.message;
      Result r;
      r.eps = data.settings.css.state_tolerance_.Get();
      if (!run.Ok()) return r;
      r.state = sim.FinalState();
      r.kpis = sim.KPIs(-1);
      r.css_cycle = run.css_reached ? run.css_cycle : -1;
      return r;
    }

    bool SameComposition(const UserData& data, const Composition& y) {
      for (const auto& [c, v] : y) {
        const auto it = data.ic.mole_fraction.find(c);
        const double ic = it == data.ic.mole_fraction.end() ? 0.0 : it->second.Get();
        if (std::abs(ic - v) > 1e-12) return false;
      }
      return true;
    }

    std::function<void(UserData&)> WithGas(const Composition& y) {
      return [y](UserData& data) {
        for (const auto& [c, v] : y) data.ic.mole_fraction[c] = v;
      };
    }

    // The shipped deck's state after one cycle, with every temperature
    // shifted by dT and its cycle count reset (the minimum cycles apply).
    CycleState ShiftedState(const std::string& deck, double dT) {
      Simulation sim = Simulation::FromFile(DeckPath(deck));
      sim.Data().settings.css.method_ = CSSMethod::NONE;
      RunOptions options;
      options.cycles = 1;
      options.on_message = [](LogLevel, const std::string&) {};
      const RunResult run = sim.Run(options);
      EXPECT_TRUE(run.Ok()) << deck << ": " << run.message;
      CycleState s = sim.FinalState();
      for (const auto& f : s.fields)
        if (f.name.find("Temperature") != std::string::npos)
          for (size_t i : f.index) s.y[i] += dT;
      s.cycles = 0;
      return s;
    }

  }  // namespace

  class CSSStartIndependence : public ::testing::TestWithParam<DeckSpec> {};

  TEST_P(CSSStartIndependence, EveryStartReachesTheSameCss) {
    const DeckSpec spec = GetParam();
    const std::string deck = spec.deck;
    Simulation shipped = Simulation::FromFile(DeckPath(deck));
    const UserData& shipped_data = shipped.Data();

    const Result base = RunToCss(deck, nullptr);
    ASSERT_GE(base.css_cycle, 0) << deck << ": the shipped start did not reach CSS";
    ASSERT_GT(base.eps, 0.0);
    std::printf("[CSSSTART] %-24s %-10s stop %4d  (state tolerance %.1e)\n", deck.c_str(), "shipped",
                base.css_cycle, base.eps);

    std::vector<std::pair<std::string, Result>> starts;
    if (!SameComposition(shipped_data, spec.clean)) starts.emplace_back("clean", RunToCss(deck, WithGas(spec.clean)));
    if (!SameComposition(shipped_data, spec.feed)) starts.emplace_back("feed", RunToCss(deck, WithGas(spec.feed)));
    const CycleState hot = ShiftedState(deck, kTemperatureShift);
    starts.emplace_back("hot", RunToCss(deck, nullptr, &hot));
    const CycleState cold = ShiftedState(deck, -kTemperatureShift);
    starts.emplace_back("cold", RunToCss(deck, nullptr, &cold));
    if (spec.neighbour) {
      const Result other = RunToCss(spec.neighbour, nullptr);
      ASSERT_GE(other.css_cycle, 0) << spec.neighbour << ": no CSS";
      CycleState from = other.state;
      from.cycles = 0;  // the minimum cycles apply again
      starts.emplace_back(std::string("neighbour"), RunToCss(deck, nullptr, &from));
    }

    for (const auto& [name, r] : starts) {
      const CssDistance d = Distance(r.state, r.kpis, base.state, base.kpis);
      std::printf("[CSSSTART] %-24s %-10s stop %4d  c_state %8.3f (%s)  c_history %8.3f  c_kpi %8.3f\n",
                  deck.c_str(), name.c_str(), r.css_cycle, d.state / base.eps, d.state_field.c_str(),
                  d.history / base.eps, d.kpi / base.eps);
      EXPECT_GE(r.css_cycle, 0) << deck << " " << name << ": no CSS within " << kMaxCycles << " cycles";
      EXPECT_LE(d.state, kBound * base.eps) << deck << " " << name << ": bed state (" << d.state_field << ")";
      EXPECT_LE(d.history, kBound * base.eps) << deck << " " << name << ": endpoint histories";
      EXPECT_LE(d.kpi, kBound * base.eps) << deck << " " << name << ": KPIs";
    }
  }

  INSTANTIATE_TEST_SUITE_P(
      ExampleDecks, CSSStartIndependence,
      ::testing::Values(DeckSpec{"Illustrative_TSA_13X", kCleanN2, kFlueGas, nullptr},
                        DeckSpec{"Perez2019_A", kCleanN2, kFlueGas, "Perez2019_B"},
                        DeckSpec{"Perez2019_E", kCleanN2, kFlueGas, "Perez2019_D1"},
                        DeckSpec{"Nguyen2023_Exp_I", kCleanN2, kFlueGas, "Nguyen2023_Exp_II"},
                        DeckSpec{"Nguyen2023_Exp_VI", kCleanN2, kFlueGas, "Nguyen2023_Exp_VII"},
                        DeckSpec{"Vysyaraju2024_PSA_exp1", {{"Helium", 1.0}, {"Nitrogen", 0.0}},
                                 {{"Helium", 0.1}, {"Nitrogen", 0.9}}, nullptr},
                        DeckSpec{"Ward2024", kCleanN2, kFlueGas, nullptr},
                        DeckSpec{"FiveStepProductPurge", kCleanN2, kFlueGas, nullptr}),
      [](const ::testing::TestParamInfo<DeckSpec>& info) { return std::string(info.param.deck); });

}  // namespace adsim
