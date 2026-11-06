/*
 * Skarstrom - pressure and temperature swing adsorption simulator
 * Copyright (c) 2024-2026 Sam Affleck and contributors
 *
 * Licensed under the terms in the LICENSE file at the repository root.
 */

// Anderson acceleration of the approach to cyclic steady state
// (AccelerationMethod::ANDERSON) on a synthetic cycle map.
//
// The end-of-cycle map here is exactly linear with two error modes,
//     x_{k+1} = x* + rho1 a1 e1 + rho2 a2 e2   when x_k = x* + a1 e1 + a2 e2,
// applied to the slow fields of a bundled deck's state (loadings, solid and
// wall temperatures, mole fractions; e1 and e2 are the first two sine modes
// over the cells). Aitken removes one mode per jump; Anderson combines the
// earlier cycle-to-cycle steps, so once it holds two independent steps its
// jump lands on x* (to rounding). The tests check that, that a caller which
// rejects jumps (the re-seed probe in CoupledSolver) does not corrupt the
// pairs Anderson keeps, and that the mole fractions still sum to one.

#include <cmath>
#include <string>
#include <vector>

#include <gtest/gtest.h>

#include "adsim/objects/user_data.h"
#include "adsim/solver/cyclic_steady_state.h"

namespace adsim {

  namespace {
    struct Field {
      Equation* eq;
      double base;  // x*, uniform over the cells
      double amp;   // error amplitude relative to base (mole fractions pair up)
    };

    double Mode(int k, size_t c, size_t n) {
      const double pi = std::acos(-1.0);
      return std::sin(pi * k * (static_cast<double>(c) + 0.5) / static_cast<double>(n));
    }

    // The slow fields, with a base value and an error amplitude each. The
    // first two mole fractions carry opposite errors and any others none, so
    // the composition sums to one.
    std::vector<Field> SlowFields(Equations& state) {
      std::vector<Field> fields;
      for (Equation* eq : state.GetEquations()) {
        if (eq->GetEquationType() != EquationType::Differential) continue;
        if (eq == &state.concentration || eq == &state.fluid_temperature) continue;
        double base = 1.0, amp = 0.01;
        if (eq == &state.solid_temperature || eq == &state.wall_temperature) base = 300.0;
        const size_t n_comp = state.mole_fraction.size();
        size_t index = 0;
        for (auto& [component, y] : state.mole_fraction) {
          if (eq == &y) {
            base = index == 0 ? 0.3 : index == 1 ? 0.6 : 0.1 / static_cast<double>(n_comp - 2);
            amp = index == 0 ? 0.03 / 0.3 : index == 1 ? -0.03 / 0.6 : 0.0;
          }
          ++index;
        }
        fields.push_back({eq, base, amp});
      }
      return fields;
    }

    class AndersonSynthetic : public ::testing::Test {
    protected:
      void SetUp() override {
        LoadFromJson(data, "../examples/decks/Nguyen2022_Fig4.json");
        std::string error;
        ASSERT_TRUE(data.PerformChecks(error)) << error;
        data.Initialise();
        ASSERT_GE(data.state.mole_fraction.size(), 2u);  // GUARD: a mixture
        fields = SlowFields(data.state);
        ASSERT_GE(fields.size(), 4u);
      }

      void Set(double a1, double a2) {
        for (const auto& f : fields) {
          const size_t n = f.eq->var.size();
          for (size_t c = 0; c < n; ++c)
            f.eq->var[c] = f.base * (1.0 + f.amp * (a1 * Mode(1, c, n) + a2 * Mode(2, c, n)));
        }
      }

      // The mode coefficients of the state, from the first loading field (the
      // sine modes are orthogonal on the cell centres), and the largest
      // relative misfit of any field against that reconstruction.
      void Read(double& a1, double& a2, double& misfit) const {
        const Field& f = fields.front();
        const size_t n = f.eq->var.size();
        double p1 = 0.0, p2 = 0.0, n1 = 0.0, n2 = 0.0;
        for (size_t c = 0; c < n; ++c) {
          const double d = (f.eq->var[c] / f.base - 1.0) / f.amp;
          p1 += d * Mode(1, c, n); n1 += Mode(1, c, n) * Mode(1, c, n);
          p2 += d * Mode(2, c, n); n2 += Mode(2, c, n) * Mode(2, c, n);
        }
        a1 = p1 / n1;
        a2 = p2 / n2;
        misfit = 0.0;
        for (const auto& g : fields) {
          const size_t m = g.eq->var.size();
          for (size_t c = 0; c < m; ++c) {
            const double expected = g.base * (1.0 + g.amp * (a1 * Mode(1, c, m) + a2 * Mode(2, c, m)));
            misfit = std::max(misfit, std::abs(g.eq->var[c] - expected) / g.base);
          }
        }
      }

      double WorstMoleFractionSum() const {
        double worst = 0.0;
        const size_t n = data.state.mole_fraction.begin()->second.var.size();
        for (size_t c = 0; c < n; ++c) {
          double sum = 0.0;
          for (const auto& [component, y] : data.state.mole_fraction) sum += y.var[c];
          worst = std::max(worst, std::abs(sum - 1.0));
        }
        return worst;
      }

      struct Run {
        int jumps = 0;
        int first_exact = -1;           // first cycle whose start state is x* to 1e-10
        std::vector<double> error;      // max(|a1|, |a2|) at the start of each cycle
      };

      // Runs `cycles` cycles of the linear map from (1, 1). reject(k) says
      // whether the caller rejects a jump made at the end of cycle k (it then
      // restores the natural end state, as CoupledSolver does when the re-seed
      // probe fails).
      template <class Reject>
      Run Cycle(AccelerationMethod method, double rho1, double rho2, int cycles, Reject reject) {
        AccelerateCyclicSteadyState acc;
        acc.method_ = method;
        Run run;
        double a1 = 1.0, a2 = 1.0;
        for (int k = 0; k < cycles; ++k) {
          run.error.push_back(std::max(std::abs(a1), std::abs(a2)));
          if (run.first_exact < 0 && run.error.back() < 1e-10) run.first_exact = k;
          const double e1 = rho1 * a1, e2 = rho2 * a2;  // natural end of cycle k
          Set(e1, e2);
          AccelerationNote note;
          const bool moved = acc.Accelerate(data.state, k, &note);
          if (moved && !reject(k)) {
            ++run.jumps;
            double misfit = 0.0;
            Read(a1, a2, misfit);
            EXPECT_LT(misfit, 1e-12) << "cycle " << k << ": the jump left the two-mode form";
            EXPECT_LT(WorstMoleFractionSum(), 1e-12) << "cycle " << k;
          } else {
            a1 = e1;
            a2 = e2;
          }
        }
        return run;
      }

      UserData data;
      std::vector<Field> fields;
    };
  }  // namespace

  // Two modes of rate 0.9 and 0.6. Pairs come from cycles 2, 3 and 4. At the
  // end of cycle 3 the fit has one earlier step, which cannot remove both
  // modes: it leaves more than half the residual and is skipped. At the end of
  // cycle 4 it combines two independent steps, and the jump lands on x*. Plain
  // substitution is still 0.9^9 = 0.39 of the way from x* at the start of
  // cycle 9.
  TEST_F(AndersonSynthetic, TwoModeLinearMapLandsOnTheFixedPoint) {
    const auto never = [](int) { return false; };
    const Run plain = Cycle(AccelerationMethod::NONE, 0.9, 0.6, 12, never);
    const Run anderson = Cycle(AccelerationMethod::ANDERSON, 0.9, 0.6, 12, never);
    const Run aitken = Cycle(AccelerationMethod::AITKEN, 0.9, 0.6, 12, never);

    EXPECT_EQ(plain.jumps, 0);
    EXPECT_GT(plain.error[9], 0.3);
    EXPECT_EQ(anderson.jumps, 1);
    EXPECT_EQ(anderson.first_exact, 5) << "Anderson should reach x* at the start of cycle 5";
    EXPECT_LT(anderson.error.back(), 1e-10);
    // Aitken removes one mode per window, so two modes stay away from x*.
    EXPECT_GT(aitken.error[6], 1e-6);
  }

  // A slow mode (0.97, which asks for a jump of 32 natural steps) is shortened
  // to kMaxStepMultiple steps, never dropped, and the run still converges
  // much faster than plain substitution.
  TEST_F(AndersonSynthetic, SlowModeJumpIsLimitedNotDropped) {
    const auto never = [](int) { return false; };
    const Run plain = Cycle(AccelerationMethod::NONE, 0.97, 0.5, 30, never);
    const Run anderson = Cycle(AccelerationMethod::ANDERSON, 0.97, 0.5, 30, never);
    EXPECT_GE(anderson.jumps, 3);
    EXPECT_LT(anderson.error.back(), 1e-3 * plain.error.back());
  }

  // A caller that rejects jumps keeps the natural end state. The next cycle
  // then gives no pair (its start state is unknown to Anderson), and the
  // earlier pairs stay valid: the run still reaches x*, later.
  TEST_F(AndersonSynthetic, RejectedJumpsDoNotCorruptThePairs) {
    int rejected = 0;
    const auto first_two = [&rejected](int) { return rejected++ < 2; };
    const Run anderson = Cycle(AccelerationMethod::ANDERSON, 0.9, 0.6, 16, first_two);
    EXPECT_EQ(rejected, anderson.jumps + 2);
    EXPECT_GE(anderson.first_exact, 0);
    EXPECT_LT(anderson.error.back(), 1e-10);
  }

}  // namespace adsim
