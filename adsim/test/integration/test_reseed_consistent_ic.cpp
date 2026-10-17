// test_reseed_consistent_ic.cpp
//
// Regression tests for L-03: IDACalcIC writes the corrected (consistent)
// initial condition only into IDA's history array (phi[0], phi[1]), never into
// the y / yp vectors passed to IDAInit / IDAReInit. SolverData::
// ValidateAndCommitReseed unpacked y_ / yp_ straight after the live IDACalcIC,
// so user_data.state (and the first logged row of the accelerated cycle) held
// the uncorrected extrapolation; it also ignored the IDAReInit / IDACalcIC
// return codes.
//
// What IDA actually integrated from is recovered independently of the code
// under test: after one accepted step, IDA's interpolating polynomial passes
// through the previous point, so IDAGetDky at the start time returns the
// consistent IC IDA started from.
//
//  * After SolverData::Initialise, user_data.state (the logged t = 0 row) is
//    the consistent IC.
//  * After a committed re-seed, user_data.state is IDA's consistent IC, not the
//    extrapolated input.
//  * A live IDACalcIC failure in the commit path (forced via a test hook) takes
//    the reject branch: the pre-reseed state is restored and integration
//    continues.

#include <algorithm>
#include <cmath>
#include <string>
#include <vector>

#include "../src/pch.h"
#include "adsim/solver/coupled_solver.h"
#include "adsim/boundary_conditions/boundary_conditions.h"
#include "adsim/boundary_conditions/boundary_stream.h"
#include "adsim/boundary_conditions/step.h"

#include <gtest/gtest.h>

#include "../fixtures/user_data_test.h"

namespace adsim {

  namespace {

    struct Packed {
      std::vector<double> y;
      std::vector<double> yp;
    };

    Packed PackState(Equations& state) {
      const int n = state.GetSystemSize();
      Packed p{ std::vector<double>(n, 0.0), std::vector<double>(n, 0.0) };
      state.PackBanded(p.y.data(), p.yp.data());
      return p;
    }

    // max_i |a_i - b_i| / (|b_i| + 1e-10), and the entry where it occurs.
    struct Diff { double rel = 0.0; size_t at = 0; };
    Diff MaxRelDiff(const std::vector<double>& a, const std::vector<double>& b) {
      Diff d;
      for (size_t i = 0; i < a.size(); ++i) {
        const double r = std::abs(a[i] - b[i]) / (std::abs(b[i]) + 1e-10);
        if (r > d.rel) { d.rel = r; d.at = i; }
      }
      return d;
    }

  } // namespace

  // Feed step (flow inlet, pressure outlet) on an initially uniform bed: the
  // t = 0 state is not consistent with the feed BCs, so IDACalcIC has real
  // corrections to make.
  class ReseedConsistentICTest : public UserDataTest {
  protected:
    void SetUp() override {
      UserDataTest::SetUp();

      AddComponents({ "Oxygen", "Nitrogen" });

      auto [layer_id, ads_id] = AddLayer("13X");
      GetLayer(layer_id).domain_.n_cells = 10;
      AddConstantInternalHeatTransferCoefficientModel(0, layer_id);
      AddConstantExternalHeatTransferCoefficientModel(0, layer_id);

      AddInitialConditions(101325, 298, 298, { {"Oxygen", 0.0}, {"Nitrogen", 1.0} });

      AddStep("Feed", 20.0);
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

      ConfigureCycle(1);
      data.settings.max_step = 1.0;
    }

    // Mirrors CoupledSolver::Initialise up to (not including) SolverData.
    void PrepareFirstStep() {
      data.Initialise();
      Step& first = data.cycle.step_list_[0];
      data.state.SetCurrentStep(&first);
      data.step_start_time_ = 0.0;
      first.Initialise(data.state);
    }

    // Integrates the live solver to t_end and loads its state into data.state.
    void IntegrateTo(SolverData& sd, double t_end) {
      ASSERT_EQ(sd.SetStopTime(t_end), IDA_SUCCESS);
      double t = sd.GetTime();
      for (int i = 0; i < 100000 && t < t_end; ++i) {
        const int f = sd.OneStep(t_end, t);
        ASSERT_GE(f, 0) << "IDASolve failed at t = " << t;
        if (f == IDA_TSTOP_RETURN) break;
      }
      ASSERT_NEAR(sd.GetTime(), t_end, 1e-9);
      sd.Unpack(data);
    }

    // The y IDA actually started integrating from at t0 (where nothing has been
    // integrated yet): take one step, then evaluate IDA's interpolant back at t0.
    Packed StartingPointOfNextStep(SolverData& sd, double t0) {
      double t = t0;
      EXPECT_EQ(sd.SetStopTime(t0 + 1.0), IDA_SUCCESS);
      EXPECT_GE(sd.OneStep(t0 + 1.0, t), 0);
      EXPECT_GT(t, t0);
      EXPECT_TRUE(sd.UnpackInterpolated(data, t0));
      return PackState(data.state);
    }

    // A stand-in for the CSS extrapolation: the live state with the (differential)
    // gas mole fractions nudged, so the algebraic variables and yp no longer
    // match it. (Nudging fluid temperature or total concentration independently
    // makes even the scratch IDACalcIC fail, so the re-seed would be rejected.)
    Packed Perturbed(const Packed& p) {
      std::vector<double> id(p.y.size(), 0.0);
      data.state.SetIDBanded(id.data());
      Packed q = p;
      for (size_t i = 0; i < q.y.size(); ++i) {
        if (id[i] == 1.0
            && data.state.DescribeBandedIndex(i).find("Mole Fraction") != std::string::npos) {
          q.y[i] *= (1.0 + 1e-4 * ((i % 3) + 1));
        }
      }
      return q;
    }

    // Tolerance for "the same vector" here: the interpolant at the previous
    // point reproduces it only up to round-off.
    static constexpr double kSame = 1e-9;
  };

  TEST_F(ReseedConsistentICTest, InitialiseLeavesConsistentICInState) {
    PrepareFirstStep();
    const Packed guess = PackState(data.state);

    SolverData sd;
    std::string err;
    ASSERT_EQ(sd.Initialise(data, err), IDA_SUCCESS) << err;
    const Packed logged = PackState(data.state); // what the t = 0 row logs

    const Packed ic = StartingPointOfNextStep(sd, 0.0);
    const Diff guess_vs_ic = MaxRelDiff(guess.y, ic.y);
    ASSERT_GT(guess_vs_ic.rel, 1e-6)
        << "the fixture's initial guess is already consistent; the test is vacuous";
    const Diff d = MaxRelDiff(logged.y, ic.y);
    EXPECT_LT(d.rel, kSame)
        << "state after Initialise is not the IC IDA integrated from: "
        << data.state.DescribeBandedIndex(d.at) << " is " << logged.y[d.at]
        << ", IDA started from " << ic.y[d.at];
  }

  TEST_F(ReseedConsistentICTest, CommittedReseedLeavesConsistentICInState) {
    PrepareFirstStep();
    SolverData sd;
    std::string err;
    ASSERT_EQ(sd.Initialise(data, err), IDA_SUCCESS) << err;
    IntegrateTo(sd, 5.0);
    const double t0 = sd.GetTime();

    const Packed seed = Perturbed(PackState(data.state));
    data.state.UnpackBanded(seed.y.data(), seed.yp.data());

    std::string why = "unset";
    ASSERT_TRUE(sd.ValidateAndCommitReseed(data, 1.0, &why)) << why;
    EXPECT_TRUE(why.empty()) << why;
    EXPECT_DOUBLE_EQ(sd.GetTime(), t0);
    const Packed committed = PackState(data.state); // the accelerated cycle's first row

    // IDA_YA_YDP_INIT keeps the differential variables and corrects the
    // algebraic ones (and the derivatives).
    std::vector<double> id(seed.y.size(), 0.0);
    data.state.SetIDBanded(id.data());
    for (size_t i = 0; i < id.size(); ++i) {
      if (id[i] == 1.0) EXPECT_EQ(committed.y[i], seed.y[i]) << "differential y[" << i << "]";
    }

    const Packed ic = StartingPointOfNextStep(sd, t0);
    const Diff seed_vs_ic = MaxRelDiff(seed.y, ic.y);
    ASSERT_GT(seed_vs_ic.rel, 1e-6)
        << "IDACalcIC made no correction to the re-seed; the test is vacuous";
    const Diff d = MaxRelDiff(committed.y, ic.y);
    EXPECT_LT(d.rel, kSame)
        << "state after the re-seed commit is not the IC IDA integrated from: "
        << data.state.DescribeBandedIndex(d.at) << " is " << committed.y[d.at]
        << ", IDA started from " << ic.y[d.at] << " (re-seed input " << seed.y[d.at] << ")";

    // The live integrator carries on from the committed state.
    IntegrateTo(sd, t0 + 5.0);
  }

  TEST_F(ReseedConsistentICTest, FailedLiveCalcICFallsBackToPreReseedState) {
    PrepareFirstStep();
    SolverData sd;
    std::string err;
    ASSERT_EQ(sd.Initialise(data, err), IDA_SUCCESS) << err;
    IntegrateTo(sd, 5.0);
    const double t0 = sd.GetTime();
    const Packed pre = PackState(data.state);

    const Packed seed = Perturbed(pre);
    data.state.UnpackBanded(seed.y.data(), seed.yp.data());

    sd.FailNextReseedCommitForTesting();
    std::string why;
    EXPECT_FALSE(sd.ValidateAndCommitReseed(data, 1.0, &why));
    EXPECT_NE(why.find("IDACalcIC"), std::string::npos) << why;
    EXPECT_DOUBLE_EQ(sd.GetTime(), t0);

    const Packed after = PackState(data.state);
    ASSERT_EQ(after.y.size(), pre.y.size());
    for (size_t i = 0; i < pre.y.size(); ++i) {
      ASSERT_EQ(after.y[i], pre.y[i]) << "y[" << i << "] not restored";
      ASSERT_EQ(after.yp[i], pre.yp[i]) << "yp[" << i << "] not restored";
    }

    // The live integrator restarts from the pre-reseed state, not from the
    // (half-applied) re-seed, and the un-accelerated trajectory continues.
    const Packed start = StartingPointOfNextStep(sd, t0);
    const Diff d = MaxRelDiff(start.y, pre.y);
    EXPECT_LT(d.rel, kSame) << "integrator did not restart from the pre-reseed state at "
                            << data.state.DescribeBandedIndex(d.at);
    IntegrateTo(sd, t0 + 5.0);

    // The hook is one-shot: the next re-seed commits normally.
    const Packed seed2 = Perturbed(PackState(data.state));
    data.state.UnpackBanded(seed2.y.data(), seed2.yp.data());
    EXPECT_TRUE(sd.ValidateAndCommitReseed(data, 1.0, &why)) << why;
  }

} // namespace adsim
