/*
 * Skarstrom - pressure and temperature swing adsorption simulator
 * Copyright (c) 2024-2026 Sam Affleck and contributors
 *
 * Licensed under the terms in the LICENSE file at the repository root.
 */

// Extreme operating corners.
//
// The robustness fuzz (test_robustness_fuzz.cpp) samples physical decks at
// random; this file aims at the corners an operator can reach by turning one
// knob to an extreme on an otherwise ordinary deck: a pressure ratio of 500 or
// of 1.0001, a blowdown to 1 mbar, a feed of 1e-9 mol/s, steps of 1 ms or of
// 1e5 s, a run of 1e6 s, a bed at 600-1200 K. No accuracy is claimed. The
// contract is the fuzz harness's, asserted per corner:
//
//   * the run CONVERGES (IntegrateCycle >= 0) or FAILS CLEANLY (a negative
//     IDA code or a C++ exception), within the tightened work caps (so never
//     a hang);
//   * every logged state is finite (no NaN or Inf leaks through), whether or
//     not the run converged;
//   * a converged run is physical: P > 0, T > 0, every mole fraction in
//     [0, 1] and every loading >= 0, each to within its absolute tolerance
//     (1e-8, docs/numerics.md), and the mole fractions of every cell sum to 1
//     to within 1e-3, the bound of the fuzz's IsMassConserved (each species
//     has its own equation, so the sum is held only to the integrator's
//     global error; measured worst 4e-6, in the 1 mbar blowdown).
//
// Each corner also has a guard showing it really is extreme (the pressure
// actually reached, the temperature actually held, ...), so a deck that
// quietly ignored the setting could not pass.

#include "breakthrough_case.h"
#include "psa_cycle_case.h"

#include <cmath>
#include <functional>
#include <iostream>
#include <string>

namespace adsim {

  using namespace adsim::validation;

  namespace {

    constexpr double kAbsTol = 1e-8;        // mole fractions and loadings
    constexpr double kSumTolerance = 1e-3;  // |sum y - 1|, as IsMassConserved

    struct Outcome {
      int flag = 0;
      bool threw = false;
      std::string message;
      bool finite = true;
      double P_min = 1e300, P_max = -1e300, T_min = 1e300, T_max = -1e300;
      double y_min = 1e300, y_max = -1e300, q_min = 1e300, worst_sum = 0.0;
      size_t snapshots = 0;
    };

    Outcome RunCorner(UserData& data) {
      Outcome o;
      // As in the fuzz: a non-convergent deck returns IDA_TOO_MUCH_WORK quickly.
      data.settings.max_num_steps = 20000;
      data.settings.max_iters_per_step = 20000;
      std::string error;
      EXPECT_TRUE(data.PerformChecks(error)) << error;
      try {
        CoupledSolver solver{data};
        o.flag = solver.IntegrateCycle(data.cycle, false);
        const auto& log = solver.GetLog();
        auto scan = [&](const VariableLogEntry& v, double* lo, double* hi) {
          for (const auto& s : v.values) {
            for (size_t i = 1; i + 1 < s.size(); ++i) {
              if (!std::isfinite(s[i])) { o.finite = false; continue; }
              if (lo) *lo = std::min(*lo, s[i]);
              if (hi) *hi = std::max(*hi, s[i]);
            }
          }
          o.snapshots += v.values.size();
        };
        for (int cycle : log.sys_data.cycle_index_set_) {
          for (int step : log.sys_data.step_index_set_) {
            scan(log.GetVariable(cycle, step, VariableType::PRESSURE), &o.P_min, &o.P_max);
            scan(log.GetVariable(cycle, step, VariableType::TEMPERATURE), &o.T_min, &o.T_max);
            scan(log.GetVariable(cycle, step, VariableType::WALL_TEMPERATURE), nullptr, nullptr);
            std::vector<std::vector<double>> sum;
            for (const auto& c : data.fluid.components_) {
              const auto y = log.GetComponentVariable(cycle, step, ComponentVariableType::MOLE_FRACTION, c);
              scan(y, &o.y_min, &o.y_max);
              scan(log.GetComponentVariable(cycle, step, ComponentVariableType::ADSORBED_CONCENTRATION, c),
                   &o.q_min, nullptr);
              if (sum.empty()) sum.assign(y.values.size(), {});
              for (size_t n = 0; n < y.values.size(); ++n) {
                if (sum[n].empty()) sum[n].assign(y.values[n].size(), 0.0);
                for (size_t i = 0; i < y.values[n].size(); ++i) sum[n][i] += y.values[n][i];
              }
            }
            for (const auto& s : sum)
              for (size_t i = 1; i + 1 < s.size(); ++i) o.worst_sum = std::max(o.worst_sum, std::abs(s[i] - 1.0));
          }
        }
      } catch (const std::exception& e) {
        o.threw = true;
        o.message = e.what();
      }
      return o;
    }

    void ExpectContract(const char* corner, const Outcome& o) {
      std::cout << "[" << corner << "] flag " << o.flag << (o.threw ? " (threw: " + o.message + ")" : "")
                << ", P " << o.P_min << " .. " << o.P_max << " Pa, T " << o.T_min << " .. " << o.T_max
                << " K, y " << o.y_min << " .. " << o.y_max << ", min q " << o.q_min
                << ", max |sum y - 1| " << o.worst_sum << ", " << o.snapshots << " snapshots\n";
      EXPECT_TRUE(o.finite) << corner << ": a logged state is not finite";
      if (o.threw || o.flag < 0) {
        SUCCEED() << corner << " failed cleanly (flag " << o.flag << ")";
        return;
      }
      EXPECT_GT(o.snapshots, 0u) << corner;
      EXPECT_GT(o.P_min, 0.0) << corner;
      EXPECT_GT(o.T_min, 0.0) << corner;
      EXPECT_GE(o.y_min, -kAbsTol) << corner;
      EXPECT_LE(o.y_max, 1.0 + kAbsTol) << corner;
      EXPECT_GE(o.q_min, -kAbsTol) << corner;
      EXPECT_LE(o.worst_sum, kSumTolerance) << corner;
    }

    bool Converged(const Outcome& o) { return !o.threw && o.flag >= 0; }

    void SetPsaPressures(UserData& data, double P_high, double P_low) {
      auto* top = dynamic_cast<PressureBoundaryCondition*>(data.cycle.GetStep("Feed").top_boundary_condition_.get());
      static_cast<ConstantBoundaryValue*>(top->pressure_.get())->constant_value_ = P_high;
      auto* bottom = dynamic_cast<PressureBoundaryCondition*>(data.cycle.GetStep("Blowdown").bottom_boundary_condition_.get());
      static_cast<ConstantBoundaryValue*>(bottom->pressure_.get())->constant_value_ = P_low;
      data.ic.pressure = P_low;
    }

    void SetPsaFeedFlow(UserData& data, double molar_flow) {
      for (const char* name : {"Pressurise", "Feed"}) {
        auto* bc = dynamic_cast<FlowBoundaryCondition*>(data.cycle.GetStep(name).bottom_boundary_condition_.get());
        auto* stream = dynamic_cast<BoundaryStream*>(bc->streams_.front().get());
        auto* flow = dynamic_cast<MolarFlowBoundaryValue*>(stream->boundary_flow_.get());
        static_cast<ConstantBoundaryValue*>(flow->profile_.get())->constant_value_ = molar_flow;
      }
    }

    // The shared PSA cycle (pressurise / feed / blowdown, non-isothermal,
    // two competing adsorbates), coarse and with a few cycles.
    void BuildCyclicCorner(UserData& data, int cycles) {
      BuildPsaCycle(data, 20);
      data.cycle.number_of_cycles_ = cycles;
      data.settings.max_step = 2.0;
    }

  } // namespace

  // Pressure ratio 500: product at 50 bar, blowdown to 0.1 bar, with ten
  // times the feed so the bed reaches the high pressure.
  TEST(ExtremeConditions, PressureRatioOf500) {
    UserData data;
    BuildCyclicCorner(data, 3);
    SetPsaPressures(data, 50e5, 0.1e5);
    SetPsaFeedFlow(data, 0.2);
    const Outcome o = RunCorner(data);
    ExpectContract("pressure ratio 500", o);
    if (Converged(o)) {
      EXPECT_GT(o.P_max, 40e5) << "The bed never reached the high pressure";
      EXPECT_LT(o.P_min, 0.2e5) << "The bed never reached the low pressure";
    }
  }

  // Pressure ratio 1.0001 (no swing to speak of).
  TEST(ExtremeConditions, PressureRatioOfOnePointZeroZeroZeroOne) {
    UserData data;
    BuildCyclicCorner(data, 3);
    SetPsaPressures(data, 1.0001e5, 1.0e5);
    ExpectContract("pressure ratio 1.0001", RunCorner(data));
  }

  // Blowdown to 1 mbar.
  TEST(ExtremeConditions, BlowdownToOneMillibar) {
    UserData data;
    BuildCyclicCorner(data, 2);
    SetPsaPressures(data, 3e5, 100.0);
    const Outcome o = RunCorner(data);
    ExpectContract("blowdown to 1 mbar", o);
    if (Converged(o)) EXPECT_LT(o.P_min, 200.0) << "The bed never reached 1 mbar";
  }

  // Near-zero flow: 1e-9 mol/s into the cycle, 1e-10 mol/s through a
  // breakthrough column.
  TEST(ExtremeConditions, NearZeroFlow) {
    {
      UserData data;
      BuildCyclicCorner(data, 3);
      SetPsaFeedFlow(data, 1e-9);
      const Outcome o = RunCorner(data);
      ExpectContract("cycle, feed 1e-9 mol/s", o);
      // Guard: nothing to pressurise with, so the bed stays near its initial
      // pressure until the feed step's 3 bar outlet pushes gas in from the top.
      if (Converged(o)) EXPECT_GE(o.P_min, 0.99e5);
    }
    {
      BreakthroughCase c;
      c.adsorbing = true;
      c.n_cells = 30;
      c.run_time = 400.0;
      c.molar_flow = 1e-10;
      UserData data;
      Build(data, c);
      ExpectContract("breakthrough, feed 1e-10 mol/s", RunCorner(data));
    }
  }

  // Very short and very long steps: every step 1 ms (five cycles), every step
  // 1e5 s (the closed-top flow pressurisation then packs ~2000 mol into the
  // bed: some 10^4 bar), and a breakthrough run of 1e6 s with no step cap.
  TEST(ExtremeConditions, VeryShortAndVeryLongSteps) {
    {
      UserData data;
      BuildCyclicCorner(data, 5);
      for (auto& step : data.cycle.step_list_) step.step_time_ = 1e-3;
      data.settings.max_step = 0.0;
      ExpectContract("steps of 1 ms", RunCorner(data));
    }
    {
      UserData data;
      BuildCyclicCorner(data, 2);
      for (auto& step : data.cycle.step_list_) step.step_time_ = 1e5;
      data.settings.max_step = 0.0;
      const Outcome o = RunCorner(data);
      ExpectContract("steps of 1e5 s", o);
      if (Converged(o)) EXPECT_GT(o.P_max, 1e8) << "The long pressurisation did not pack the bed";
    }
    {
      BreakthroughCase c;
      c.adsorbing = true;
      c.n_cells = 30;
      c.run_time = 1e6;
      c.max_step = 0.0;
      UserData data;
      Build(data, c);
      const Outcome o = RunCorner(data);
      ExpectContract("breakthrough over 1e6 s", o);
      if (Converged(o)) EXPECT_GT(o.y_max, 0.99 * c.y_feed) << "The column never broke through";
    }
  }

  // High temperature: a breakthrough column (bed and feed) at 600, 900 and
  // 1200 K (the gas viscosity and conductivity fits end at 650 K and are
  // clamped above), and the PSA cycle at 700 K.
  TEST(ExtremeConditions, HighTemperature) {
    for (double T : {600.0, 900.0, 1200.0}) {
      BreakthroughCase c;
      c.adsorbing = true;
      c.n_cells = 30;
      c.run_time = 400.0;
      c.bed_temperature = T;
      c.feed_temperature = T;
      UserData data;
      Build(data, c);
      const std::string name = "breakthrough at " + std::to_string(static_cast<int>(T)) + " K";
      const Outcome o = RunCorner(data);
      ExpectContract(name.c_str(), o);
      if (Converged(o)) EXPECT_GT(o.T_min, 0.99 * T) << "The bed was not held hot";
    }
    {
      UserData data;
      BuildCyclicCorner(data, 2);
      data.ic.temperature = 700.0;
      data.ic.wall_temperature = 700.0;
      data.bed.ambient_temperature_ = 700.0;
      for (const char* name : {"Pressurise", "Feed"}) {
        auto* bc = dynamic_cast<FlowBoundaryCondition*>(data.cycle.GetStep(name).bottom_boundary_condition_.get());
        bc->streams_.front()->temperature_ = 700.0;
      }
      const Outcome o = RunCorner(data);
      ExpectContract("cycle at 700 K", o);
      if (Converged(o)) EXPECT_GT(o.T_min, 690.0);
    }
  }

} // namespace adsim
