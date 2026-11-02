// The cyclic-steady-state stopping rule (CyclicSteadyState) on synthetic
// histories, and the loading of decks saved before the rule (class version < 2).
//
// The histories mimic what a slow PSA cycle does on its way to CSS (the
// Yancy-Caballero 2020 13X cycle, case C2, at 10 cells): the per-cycle
// component balance error is a decaying oscillation that swings through zero
// (the bed alternately gains and loses CO2), while the cycle-to-cycle change of
// the bed state decays geometrically and never crosses zero. A single-cycle
// balance test with a 0.5 % band passes on the fifth cycle, long before the
// purity and recovery have settled.

#include <cmath>
#include <optional>
#include <sstream>
#include <string>
#include <vector>

#include <gtest/gtest.h>

#include <cereal/archives/json.hpp>

#include "adsim/solver/cyclic_steady_state.h"

namespace adsim {

  namespace {

    // The first cycle (0-based) at which the mass-balance rule stops: the
    // single-cycle test has passed on `k` consecutive cycles and at least
    // `min_cycles` cycles have run (CyclicSteadyState::CheckCSS on numbers).
    std::optional<int> StoppingCycle(const std::vector<double>& balance, const std::vector<double>& change,
                                     double balance_tolerance, double state_tolerance, int k,
                                     int min_cycles = 5) {
      auto passes = [&](int c) {
        return CyclicSteadyState::MassBalanceCyclePasses(balance[c], balance_tolerance, change[c],
                                                         state_tolerance);
      };
      for (int c = 0; c < static_cast<int>(balance.size()); ++c) {
        if (c >= min_cycles - 1 && CyclicSteadyState::TrailingPasses(c, k, passes) >= k) return c;
      }
      return std::nullopt;
    }

    // Synthetic slow approach to CSS. The KPI error decays as rho^c; the balance
    // error is that decay times a slow oscillation (period ~40 cycles), plus a
    // constant offset standing in for a boundary-flux defect; the state change is
    // the (non-negative) cycle-to-cycle step of the decay. `kpi_error[c]` is the
    // distance of the KPIs from their converged values.
    struct History {
      std::vector<double> balance, change, kpi_error;
    };

    History SlowApproach(int cycles, double offset = 0.0) {
      History h;
      const double rho = 0.92, pi = std::acos(-1.0);
      for (int c = 0; c < cycles; ++c) {
        const double decay = std::pow(rho, c);
        h.kpi_error.push_back(6e-3 * decay);
        h.balance.push_back(std::abs(offset + 2e-2 * decay * std::cos(2.0 * pi * c / 40.0 + 1.2)));
        h.change.push_back(c == 0 ? -1.0 : 5e-2 * decay * (1.0 - rho));
      }
      return h;
    }

  }  // namespace

  TEST(CSSStoppingRule, SingleCycleBalanceStopsByAccident) {
    // The old rule (balance only, one cycle): the balance is inside 0.5 % at the
    // minimum cycle count, so it stops while the KPIs are still far off.
    const History h = SlowApproach(300);
    const auto stop = StoppingCycle(h.balance, h.change, 0.005, 0.0, 1);
    ASSERT_TRUE(stop.has_value());
    EXPECT_LT(*stop, 30);
    EXPECT_GT(h.kpi_error[*stop], 1e-3);
  }

  TEST(CSSStoppingRule, ConsecutiveCyclesAloneDoNotFixALooseBand) {
    // The balance stays inside a loose band for many cycles in a row near a
    // zero crossing, so requiring consecutive passes is not enough on its own.
    const History h = SlowApproach(300);
    const auto stop = StoppingCycle(h.balance, h.change, 0.005, 0.0, 3);
    ASSERT_TRUE(stop.has_value());
    EXPECT_GT(h.kpi_error[*stop], 1e-3);
  }

  TEST(CSSStoppingRule, StateTestStopsCloseToCSS) {
    const History h = SlowApproach(300);
    const auto stop = StoppingCycle(h.balance, h.change, 0.005, CyclicSteadyState::kDefaultStateTolerance,
                                    CyclicSteadyState::kDefaultConsecutiveCycles);
    ASSERT_TRUE(stop.has_value());
    EXPECT_GT(*stop, 30);
    EXPECT_LT(h.kpi_error[*stop], 1e-4);
    // The state change was small on every cycle of the window.
    for (int c = *stop - CyclicSteadyState::kDefaultConsecutiveCycles + 1; c <= *stop; ++c)
      EXPECT_LE(h.change[c], CyclicSteadyState::kDefaultStateTolerance);
  }

  TEST(CSSStoppingRule, StateTestIsIndependentOfABalanceOffset) {
    // A constant balance offset inside the band (a boundary-flux defect, which
    // a much tighter balance tolerance could never get past) moves nothing.
    const History clean = SlowApproach(300), offset = SlowApproach(300, 3e-3);
    const double tol = CyclicSteadyState::kDefaultStateTolerance;
    const int k = CyclicSteadyState::kDefaultConsecutiveCycles;
    EXPECT_EQ(StoppingCycle(clean.balance, clean.change, 0.005, tol, k),
              StoppingCycle(offset.balance, offset.change, 0.005, tol, k));
    // A balance-only rule tight enough for the clean history never converges.
    EXPECT_FALSE(StoppingCycle(offset.balance, offset.change, 1e-4, 0.0, k).has_value());
  }

  TEST(CSSStoppingRule, ASingleQuietCycleDoesNotStop) {
    // One accidental quiet cycle (a dip in the state change) inside a drifting
    // stretch passes the single-cycle test but not the consecutive rule.
    std::vector<double> balance(40, 1e-3), change(40, 1e-3);
    change[0] = -1.0;
    change[20] = 1e-6;
    EXPECT_EQ(StoppingCycle(balance, change, 0.005, 5e-5, 1), 20);
    EXPECT_FALSE(StoppingCycle(balance, change, 0.005, 5e-5, 3).has_value());
  }

  TEST(CSSStoppingRule, MinimumCyclesAndMissingState) {
    // Everything passes from cycle 1: the stop is the minimum-cycle floor.
    std::vector<double> balance(20, 0.0), change(20, 0.0);
    change[0] = -1.0;  // cycle 0 has no previous cycle
    EXPECT_EQ(StoppingCycle(balance, change, 0.005, 5e-5, 3, 5), 4);
    EXPECT_EQ(StoppingCycle(balance, change, 0.005, 5e-5, 3, 1), 3);   // cycles 1, 2, 3
    EXPECT_EQ(StoppingCycle(balance, change, 0.005, 0.0, 3, 1), 2);    // state off: cycles 0, 1, 2
    EXPECT_FALSE(CyclicSteadyState::MassBalanceCyclePasses(0.0, 0.005, -1.0, 5e-5));
    EXPECT_TRUE(CyclicSteadyState::MassBalanceCyclePasses(0.0, 0.005, -1.0, 0.0));
    EXPECT_FALSE(CyclicSteadyState::MassBalanceCyclePasses(std::nan(""), 0.005, 0.0, 5e-5));
  }

  TEST(CSSStoppingRule, CyclesToLog) {
    CyclicSteadyState css;
    EXPECT_EQ(css.CyclesToLog(), 0);  // None
    css.method_ = CSSMethod::MASS_BALANCE_COMPONENT;
    EXPECT_EQ(css.CyclesToLog(), CyclicSteadyState::kDefaultConsecutiveCycles + 1);
    css.state_tolerance_ = 0.0;
    css.consecutive_cycles_ = 1;
    EXPECT_EQ(css.CyclesToLog(), 1);
    css.method_ = CSSMethod::MOLAR_PROFILE;
    EXPECT_EQ(css.CyclesToLog(), 2);
    css.consecutive_cycles_ = 0;  // counts as 1
    EXPECT_EQ(css.CyclesToLog(), 2);
  }

  namespace {

    // Serialises `css`, then rewrites the document to class version 1 and drops
    // the fields version 2 added (what a deck saved before them contains).
    CyclicSteadyState RoundTrip(const CyclicSteadyState& css, bool as_version_1) {
      std::ostringstream out;
      {
        cereal::JSONOutputArchive archive(out);
        archive(cereal::make_nvp("css", css));
      }
      std::string json = out.str();
      if (as_version_1) {
        const std::string current = "\"cereal_class_version\": 2";
        const auto pos = json.find(current);
        EXPECT_NE(pos, std::string::npos) << json;
        json.replace(pos, current.size(), "\"cereal_class_version\": 1");
        // Cut from the comma before "state_tolerance" to the end of "consecutive_cycles".
        const auto start = json.rfind(',', json.find("\"state_tolerance\""));
        const auto key = json.find("\"consecutive_cycles\"");
        EXPECT_NE(key, std::string::npos);
        const auto end = json.find_first_of("\n}", json.find(':', key));
        json.erase(start, end - start);
      }
      std::istringstream in(json);
      cereal::JSONInputArchive archive(in);
      CyclicSteadyState loaded;
      loaded.state_tolerance_ = 0.123;   // a stale value must not survive the load
      loaded.consecutive_cycles_ = 7;
      archive(cereal::make_nvp("css", loaded));
      return loaded;
    }

  }  // namespace

  TEST(CSSStoppingRule, VersionOneDecksLoadWithTheRobustRule) {
    CyclicSteadyState css;
    css.method_ = CSSMethod::MASS_BALANCE_COMPONENT;
    css.min_cycles_ = 9;
    css.mass_balance_tolerance_ = 0.002;
    const CyclicSteadyState loaded = RoundTrip(css, true);
    EXPECT_EQ(loaded.method_, CSSMethod::MASS_BALANCE_COMPONENT);
    EXPECT_EQ(loaded.min_cycles_, 9);
    EXPECT_DOUBLE_EQ(loaded.mass_balance_tolerance_.Get(), 0.002);
    EXPECT_DOUBLE_EQ(loaded.state_tolerance_.Get(), CyclicSteadyState::kDefaultStateTolerance);
    EXPECT_EQ(loaded.consecutive_cycles_, CyclicSteadyState::kDefaultConsecutiveCycles);
  }

  TEST(CSSStoppingRule, NewFieldsRoundTrip) {
    CyclicSteadyState css;
    css.method_ = CSSMethod::MASS_BALANCE_TOTAL;
    css.state_tolerance_ = 0.0;
    css.consecutive_cycles_ = 1;
    const CyclicSteadyState loaded = RoundTrip(css, false);
    EXPECT_EQ(loaded.method_, CSSMethod::MASS_BALANCE_TOTAL);
    EXPECT_DOUBLE_EQ(loaded.state_tolerance_.Get(), 0.0);
    EXPECT_EQ(loaded.consecutive_cycles_, 1);
  }

}  // namespace adsim
