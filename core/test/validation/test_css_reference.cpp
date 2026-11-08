// test_css_reference.cpp
//
// Cyclic-steady-state stops against a long reference, on the example decks.
//
// For each deck: a plain run to CSS (the deck's own rule), runs with Aitken
// and with Anderson acceleration, and a warm start from a perturbed state are
// compared with a reference: the plain stop continued for a fixed number of
// cycles with the CSS check off. The comparison covers the full end-of-cycle
// state (loadings, mole fractions, temperatures, pressure), the endpoint
// histories and the KPIs (css_state_compare.h).
//
// The tolerance is tied to the stop rule: a run that stops at the state
// tolerance eps must be within c * eps of the reference, with one c for every
// deck and method (not tuned per deck). The measured c of each deck and method
// is printed and tabulated in docs/cyclic-steady-state.md.
//
// Reference. A CSS tolerance 100x tighter than the default (2e-7) is below the
// noise of the state change itself (about 1e-6 at the default reltol of
// 1e-6), so such a run stops only by chance. The reference is therefore a
// fixed long run: from the plain stop, two further segments of
// max(kMinSegment, plain stop) cycles each with the CSS check off. The
// distance between the ends of the two segments (the reference drift) bounds
// how far the reference itself still is from the true CSS, and is checked to
// be well inside the bound.
//
// Warm start. The reference end state with every loading scaled by 0.8 and
// every temperature raised by 10 K, run to CSS with the deck's rule (its
// cycle count reset, so the minimum cycles apply again).
//
// Subsample. The example decks that cycle to CSS stop between cycle 26 and
// cycle 1554; a reference of three plain lengths for all of them would cost
// hours. The test runs eight decks that stop by cycle ~200, at least one from
// each family and cycle type (TSA, VPSA, PSA, five-step and the product-purge
// cycle); the decks left out, and why, are listed in docs/cyclic-steady-state.md.
// Slow tier: about 1-9 minutes per deck.

#include <cstdio>
#include <string>

#include <gtest/gtest.h>

#include "adsim/simulation.h"
#include "css_state_compare.h"

namespace adsim {

  using namespace adsim::validation;

  namespace {

    // Every distance must be within kBound times the stop's state tolerance.
    // One value for all decks and methods. With a geometric approach of rate
    // rho the error left at a stop is about eps * rho / (1 - rho), so slow
    // decks end further away. Measured: c up to 37 (Perez 2019 E, plain and
    // warm start); 100 covers it with margin.
    constexpr double kBound = 100.0;
    // The reference drift must be below this multiple of the state tolerance:
    // the reference then resolves the bound with a factor of 100 to spare. A
    // stop closer to the reference than the drift is at the reference's
    // resolution (Anderson on Perez 2019 E: c_kpi 0.5 against a drift of 0.3).
    // Measured drifts: 0.003-0.34 (state, histories and KPIs), except the KPI
    // drift of Perez 2019 E, 1.28 since its decks hold the gas in the bed
    // voids only (its recovery approaches CSS slowly; state 0.20, histories
    // 0.27). This bound is a regression guard on the reference, not a paper
    // tolerance: 1.5 covers that deck and still resolves kBound 65-fold.
    constexpr double kMaxReferenceDrift = 1.5;
    constexpr int kMinSegment = 100;   // cycles per reference segment, at least
    constexpr int kMaxCycles = 600;    // cap on any run to CSS

    struct Result {
      CycleState state;
      std::vector<KPIRecord> kpis;
      int css_cycle = -1;  // -1: no CSS (or CSS off)
    };

    std::string DeckPath(const std::string& deck) { return "../examples/decks/" + deck + ".json"; }

    // One run of `deck`: to CSS with `method` (css true) or for `cycles`
    // cycles with the CSS check off (css false), from `from` when given.
    Result RunDeck(const std::string& deck, AccelerationMethod method, bool css, int cycles,
                   const CycleState* from, double* eps = nullptr) {
      Simulation sim = Simulation::FromFile(DeckPath(deck));
      auto& data = sim.Data();
      // Decks without a rule (they run a fixed number of cycles) use the
      // default per-component rule here.
      if (data.settings.css.method_ == CSSMethod::NONE) data.settings.css.method_ = CSSMethod::MASS_BALANCE_COMPONENT;
      if (eps) *eps = data.settings.css.state_tolerance_.Get();
      if (!css) data.settings.css.method_ = CSSMethod::NONE;
      data.settings.acc_css.method_ = method;
      RunOptions options;
      options.cycles = cycles;
      options.initial_state = from;
      options.on_message = [](LogLevel, const std::string&) {};
      const RunResult run = sim.Run(options);
      EXPECT_TRUE(run.Ok()) << deck << ": " << run.message;
      Result r;
      if (!run.Ok()) return r;
      r.state = sim.FinalState();
      r.kpis = sim.KPIs(-1);
      r.css_cycle = run.css_reached ? run.css_cycle : -1;
      return r;
    }

    CycleState Perturbed(const CycleState& s) {
      CycleState p = s;
      for (const auto& f : p.fields) {
        if (f.name.find("Adsorbed Concentration") != std::string::npos)
          for (size_t i : f.index) p.y[i] *= 0.8;
        else if (f.name.find("Temperature") != std::string::npos)
          for (size_t i : f.index) p.y[i] += 10.0;
      }
      p.cycles = 0;
      return p;
    }

    void Report(const std::string& deck, const char* method, int stop, const CssDistance& d, double eps) {
      std::printf("[CSSREF] %-24s %-9s stop %4d  c_state %8.3f (%s)  c_history %8.3f  c_kpi %8.3f\n",
                  deck.c_str(), method, stop, d.state / eps, d.state_field.c_str(), d.history / eps, d.kpi / eps);
    }

    void Check(const std::string& deck, const char* method, int stop, const CssDistance& d, double eps) {
      Report(deck, method, stop, d, eps);
      EXPECT_GE(stop, 0) << deck << " " << method << ": no CSS within " << kMaxCycles << " cycles";
      EXPECT_LE(d.state, kBound * eps) << deck << " " << method << ": bed state (" << d.state_field << ")";
      EXPECT_LE(d.history, kBound * eps) << deck << " " << method << ": endpoint histories";
      EXPECT_LE(d.kpi, kBound * eps) << deck << " " << method << ": KPIs";
    }

  }  // namespace

  class CSSReference : public ::testing::TestWithParam<const char*> {};

  TEST_P(CSSReference, StopsAndWarmStartAgreeWithALongRun) {
    const std::string deck = GetParam();
    double eps = 0.0;

    const Result plain = RunDeck(deck, AccelerationMethod::NONE, true, kMaxCycles, nullptr, &eps);
    ASSERT_GE(plain.css_cycle, 0) << deck << ": the plain run did not reach CSS";
    ASSERT_GT(eps, 0.0) << deck << ": the stop rule has no state tolerance";

    // Reference: two segments from the plain stop, CSS check off.
    const int segment = std::max(kMinSegment, plain.css_cycle + 1);
    const Result half = RunDeck(deck, AccelerationMethod::NONE, false, segment, &plain.state);
    ASSERT_FALSE(half.state.Empty());
    const Result ref = RunDeck(deck, AccelerationMethod::NONE, false, segment, &half.state);
    ASSERT_FALSE(ref.state.Empty());
    const CssDistance drift = Distance(half.state, half.kpis, ref.state, ref.kpis);
    std::printf("[CSSREF] %-24s reference: 2 x %d cycles after the stop; drift c_state %.4f c_history %.4f c_kpi %.4f\n",
                deck.c_str(), segment, drift.state / eps, drift.history / eps, drift.kpi / eps);
    EXPECT_LE(drift.state, kMaxReferenceDrift * eps) << deck << ": the reference has not settled";
    EXPECT_LE(drift.history, kMaxReferenceDrift * eps) << deck << ": the reference has not settled";
    EXPECT_LE(drift.kpi, kMaxReferenceDrift * eps) << deck << ": the reference has not settled";

    Check(deck, "plain", plain.css_cycle, Distance(plain.state, plain.kpis, ref.state, ref.kpis), eps);

    const Result aitken = RunDeck(deck, AccelerationMethod::AITKEN, true, kMaxCycles, nullptr);
    Check(deck, "Aitken", aitken.css_cycle, Distance(aitken.state, aitken.kpis, ref.state, ref.kpis), eps);

    const Result anderson = RunDeck(deck, AccelerationMethod::ANDERSON, true, kMaxCycles, nullptr);
    Check(deck, "Anderson", anderson.css_cycle, Distance(anderson.state, anderson.kpis, ref.state, ref.kpis), eps);

    const CycleState perturbed = Perturbed(ref.state);
    const Result warm = RunDeck(deck, AccelerationMethod::NONE, true, kMaxCycles, &perturbed);
    Check(deck, "warm", warm.css_cycle, Distance(warm.state, warm.kpis, ref.state, ref.kpis), eps);
  }

  INSTANTIATE_TEST_SUITE_P(ExampleDecks, CSSReference,
                           ::testing::Values("Illustrative_TSA_13X", "Perez2019_A", "Perez2019_E", "Nguyen2023_Exp_I",
                                             "Nguyen2023_Exp_VI", "Vysyaraju2024_PSA_exp1", "Ward2024",
                                             "FiveStepProductPurge"),
                           [](const ::testing::TestParamInfo<const char*>& info) { return std::string(info.param); });

}  // namespace adsim
