/*
 * Skarstrom - pressure and temperature swing adsorption simulator
 * Copyright (c) 2024-2026 Sam Affleck and contributors
 *
 * Licensed under the terms in the LICENSE file at the repository root.
 */

// Robustness guards: the solver must never silently produce a NaN.
//
// These tests drive the discretisation hard enough to exercise the robustness
// guards added alongside them (denominator floors on P/T and the wall annulus,
// and clamping of the limited mole-fraction face reconstruction). The contract
// asserted here is deliberately weak on *accuracy* and strong on *safety*:
//
//   - The run either COMPLETES (IntegrateCycle returns >= 0) or FAILS with a
//     specific negative IDA return code -- it must not hang or return a NaN.
//   - Every logged field is finite at every logged instant (no Inf/NaN leaks
//     into the state, which would mean a guarded denominator slipped through).
//   - Every reconstructed mole-fraction FACE value stays within [0,1], so a
//     high-order limiter overshoot at a steep front cannot inject a negative or
//     >1 composition into the convective flux.
//
// A steep composition front (large feed fraction, high flow, coarse-ish mesh)
// is the natural stress case for the face-value clamp; the P/T floors are
// exercised implicitly because any transient non-positive state the BDF solver
// iterates through must not poison the residual vector.

#include "breakthrough_case.h"

#include <cmath>

namespace adsim {

  using namespace adsim::validation;

  namespace {
    // Assert every logged snapshot of a cell/face variable is finite.
    void ExpectAllFinite(const VariableLogEntry& v, const char* what) {
      for (size_t n = 0; n < v.values.size(); ++n) {
        for (size_t i = 0; i < v.values[n].size(); ++i) {
          ASSERT_TRUE(std::isfinite(v.values[n][i]))
            << "Non-finite " << what << " at snapshot " << n << ", cell " << i
            << " = " << v.values[n][i];
        }
      }
    }
  } // namespace

  // A steep, fast composition front with a moderately coarse mesh: the classic
  // case where a limited (WENO / van Leer) reconstruction can overshoot the
  // [0,1] bound at the front. The face-value clamp must keep every reconstructed
  // face composition physical, and no field may become NaN.
  TEST(RobustnessGuards, SteepFrontStaysFiniteAndBounded) {
    BreakthroughCase c;
    c.adsorbing = true;      // adsorption sharpens the front further
    c.n_cells = 60;
    c.molar_flow = 0.02;     // fast: strong convective front
    c.y_feed = 0.9;          // near-pure feed against a near-pure bed -> steep step
    c.run_time = 600.0;
    c.max_step = 1.0;
    c.reltol = 1e-8;

    UserData data;
    Build(data, c);

    CoupledSolver solver{data};
    std::string error;
    ASSERT_TRUE(data.PerformChecks(error)) << error;

    const int flag = solver.IntegrateCycle(data.cycle, false);
    // Contract: complete (>= 0) OR fail with a specific code. Either way the run
    // must have terminated (this line reached) rather than hung, and whatever it
    // produced must be finite below.
    SUCCEED() << "Solver returned flag " << flag;

    const auto& log = solver.GetLog();

    // Every logged state field must be finite -- a guarded denominator that
    // slipped through would surface here as an Inf/NaN.
    ExpectAllFinite(log.GetVariable(0, 0, VariableType::PRESSURE), "pressure");
    ExpectAllFinite(log.GetVariable(0, 0, VariableType::TEMPERATURE), "temperature");
    ExpectAllFinite(log.GetVariable(0, 0, VariableType::WALL_TEMPERATURE), "wall temperature");

    // Reconstructed mole-fraction FACE values must stay physical everywhere.
    for (const char* component : {"Nitrogen", "Oxygen"}) {
      const auto y = log.GetComponentFaceVariable(
          0, 0, ComponentVariableType::MOLE_FRACTION, component);
      double worst_low = 0.0;   // how far below 0 any face value went
      double worst_high = 0.0;  // how far above 1 any face value went
      for (size_t n = 0; n < y.values.size(); ++n) {
        for (size_t i = 0; i < y.values[n].size(); ++i) {
          const double yf = y.values[n][i];
          ASSERT_TRUE(std::isfinite(yf))
            << "Non-finite face " << component << " at snapshot " << n << ", face " << i;
          worst_low = std::max(worst_low, -yf);
          worst_high = std::max(worst_high, yf - 1.0);
        }
      }
      // The clamp is exact, so the only slack allowed is floating-point noise.
      EXPECT_LT(worst_low, 1e-12)
        << "Face " << component << " went negative by " << worst_low;
      EXPECT_LT(worst_high, 1e-12)
        << "Face " << component << " exceeded 1 by " << worst_high;
    }
  }

  // A thin column wall exercises the (d_o^2 - d_i^2) floor in the wall energy
  // balance: with wall heat transfer on and a very thin wall, the annulus
  // cross-section is small but the run must still stay finite (no divide-by-zero
  // Inf/NaN from a near-degenerate wall).
  TEST(RobustnessGuards, ThinWallEnergyBalanceStaysFinite) {
    BreakthroughCase c;
    c.adsorbing = false;
    c.run_time = 2000.0;
    c.max_step = 10.0;
    c.feed_temperature = 360.0;  // a real thermal front so the wall term is active
    c.bed_temperature = 300.0;

    UserData data;
    Build(data, c);

    // Make the wall extremely thin so d_o^2 - d_i^2 is tiny, then turn on wall
    // heat transfer so the (now near-singular) denominator is actually used.
    data.bed.wall_.thickness_ = 1e-7;  // [m]
    auto& layer = data.bed.layers_.layers_.front();
    auto ihtc = std::make_unique<ConstantInternalHeatTransferModel>();
    ihtc->heat_transfer_coefficient_ = 20.0;
    layer.equations_.internal_heat_transfer_model_ = std::move(ihtc);
    auto ehtc = std::make_unique<ConstantExternalHeatTransferModel>();
    ehtc->heat_transfer_coefficient_ = 5.0;
    layer.equations_.external_heat_transfer_model_ = std::move(ehtc);

    CoupledSolver solver{data};
    std::string error;
    ASSERT_TRUE(data.PerformChecks(error)) << error;

    const int flag = solver.IntegrateCycle(data.cycle, false);
    SUCCEED() << "Solver returned flag " << flag;

    const auto& log = solver.GetLog();
    ExpectAllFinite(log.GetVariable(0, 0, VariableType::TEMPERATURE), "temperature");
    ExpectAllFinite(log.GetVariable(0, 0, VariableType::WALL_TEMPERATURE), "wall temperature");
    ExpectAllFinite(log.GetVariable(0, 0, VariableType::PRESSURE), "pressure");
  }

} // namespace adsim
