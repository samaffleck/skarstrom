// test_solver_reinitialise.cpp
//
// Regression tests for E-13: SolverData::Initialise allocated a fresh IDA
// memory block, six N_Vectors, a band matrix and a band linear solver on every
// call without releasing the previous set, and only the destructor freed the
// last one. The GUI optimiser and parameteriser construct one CoupledSolver
// and call IntegrateCycle once per objective evaluation, so every evaluation
// leaked the whole integrator.
//
//  * Re-running IntegrateCycle on the same CoupledSolver must reproduce the
//    first run bit-for-bit (state and integrator statistics).
//  * On glibc, the heap in use must not grow by roughly one integrator per
//    re-run. (Built with -fsanitize=address, LeakSanitizer reports the leak
//    directly; the heap probe is skipped there.)

#include <string>
#include <vector>

#if defined(__GLIBC__) && !defined(__SANITIZE_ADDRESS__)
#include <malloc.h>
#define ADSIM_HEAP_PROBE 1
#endif

#include "../src/pch.h"
#include "adsim/solver/coupled_solver.h"
#include "adsim/boundary_conditions/boundary_conditions.h"
#include "adsim/boundary_conditions/boundary_stream.h"
#include "adsim/boundary_conditions/step.h"

#include <gtest/gtest.h>

#include "../fixtures/user_data_test.h"

namespace adsim {

  namespace {

    struct RunResult {
      int retval = 0;
      std::vector<double> y;
      std::vector<double> yp;
      IntegrationStats stats;
    };

#ifdef ADSIM_HEAP_PROBE
    size_t HeapInUse() { return mallinfo2().uordblks; }
#endif

  } // namespace

  // Small two-step deck: a short feed step followed by a closed hold step.
  class SolverReinitialiseTest : public UserDataTest {
  protected:
    void SetUp() override {
      UserDataTest::SetUp();

      AddComponents({ "Oxygen", "Nitrogen" });

      auto [layer_id, ads_id] = AddLayer("13X");
      GetLayer(layer_id).domain_.n_cells = 10;
      AddConstantInternalHeatTransferCoefficientModel(0, layer_id);
      AddConstantExternalHeatTransferCoefficientModel(0, layer_id);

      AddInitialConditions(101325, 298, 298, { {"Oxygen", 0.0}, {"Nitrogen", 1.0} });

      AddStep("Feed", 5.0);
      Step& feed = data.cycle.GetStep("Feed");
      auto feed_flow_bc = std::make_unique<FlowBoundaryCondition>();
      feed_flow_bc->direction_ = BoundaryDirection::INLET;
      feed_flow_bc->location_ = BoundaryLocation::BOTTOM;
      feed_flow_bc->streams_.push_back(CreateNoValveInletStream(
        CreateMolarFlowBC(0.1), { {"Oxygen", 0.21}, {"Nitrogen", 0.79} }, 298));
      feed.bottom_boundary_condition_ = std::move(feed_flow_bc);

      auto feed_pressure_out = std::make_unique<PressureBoundaryCondition>();
      feed_pressure_out->direction_ = BoundaryDirection::OUTLET;
      feed_pressure_out->location_ = BoundaryLocation::TOP;
      feed_pressure_out->pressure_ = CreateConstantBoundaryValue(101325);
      feed_pressure_out->mole_fraction_["Oxygen"] = 0.0;
      feed_pressure_out->mole_fraction_["Nitrogen"] = 1.0;
      feed.top_boundary_condition_ = std::move(feed_pressure_out);

      AddStep("Hold", 5.0);
      AddClosedBC(BoundaryLocation::BOTTOM, "Hold");
      AddClosedBC(BoundaryLocation::TOP, "Hold");

      ConfigureCycle(2);
      data.settings.max_step = 1.0;
    }

    RunResult Run(CoupledSolver& solver) {
      RunResult r;
      r.retval = solver.IntegrateCycle(data.cycle, false);
      const int n = data.state.GetSystemSize();
      r.y.assign(n, 0.0);
      r.yp.assign(n, 0.0);
      data.state.PackBanded(r.y.data(), r.yp.data());
      r.stats = solver.GetStatistics();
      return r;
    }
  };

  TEST_F(SolverReinitialiseTest, RepeatedIntegrateCycleIsReproducible) {
    CoupledSolver solver{data};
    const RunResult first = Run(solver);
    ASSERT_GE(first.retval, 0) << "first run did not complete";
    ASSERT_FALSE(first.y.empty());
    ASSERT_GT(first.stats.num_steps, 0);

    for (int run = 1; run < 4; ++run) {
      const RunResult again = Run(solver);
      ASSERT_EQ(again.retval, first.retval) << "run " << run;
      ASSERT_EQ(again.y.size(), first.y.size()) << "run " << run;
      for (size_t i = 0; i < first.y.size(); ++i) {
        ASSERT_EQ(again.y[i], first.y[i]) << "run " << run << ", y[" << i << "]";
        ASSERT_EQ(again.yp[i], first.yp[i]) << "run " << run << ", yp[" << i << "]";
      }
      EXPECT_EQ(again.stats.num_steps, first.stats.num_steps) << "run " << run;
      EXPECT_EQ(again.stats.num_res_evals, first.stats.num_res_evals) << "run " << run;
      EXPECT_EQ(again.stats.num_jac_evals, first.stats.num_jac_evals) << "run " << run;
    }
  }

#ifdef ADSIM_HEAP_PROBE
  TEST_F(SolverReinitialiseTest, RepeatedIntegrateCycleDoesNotLeakIntegrator) {
    CoupledSolver solver{data};
    // Warm-up runs: populate one-off caches (fluid properties, logs, the
    // integrator itself) so later runs only reuse or replace them.
    for (int i = 0; i < 2; ++i) ASSERT_GE(solver.IntegrateCycle(data.cycle, false), 0);

    // Lower bound on what one leaked integrator would hold: the band matrix
    // alone stores sys_size * (smu + ml + 1) doubles, smu = mu + ml.
    const size_t sys_size = static_cast<size_t>(data.state.GetSystemSize());
    const size_t bw = data.state.GetEquations().size() * 2;
    const size_t per_run_leak = sys_size * (3 * bw + 1) * sizeof(double);

    constexpr int kRuns = 5;
    const size_t before = HeapInUse();
    for (int i = 0; i < kRuns; ++i) ASSERT_GE(solver.IntegrateCycle(data.cycle, false), 0);
    const size_t after = HeapInUse();

    const long long growth = static_cast<long long>(after) - static_cast<long long>(before);
    RecordProperty("heap_growth_bytes", std::to_string(growth));
    RecordProperty("one_band_matrix_bytes", std::to_string(per_run_leak));
    EXPECT_LT(growth, static_cast<long long>(per_run_leak))
        << "heap grew by " << growth << " bytes over " << kRuns
        << " re-runs; one leaked band matrix is " << per_run_leak << " bytes";
  }
#endif

} // namespace adsim
