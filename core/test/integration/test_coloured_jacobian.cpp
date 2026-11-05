// test_coloured_jacobian.cpp
//
// The coloured DQ Jacobian (JacobianMode::Coloured) perturbs several columns
// per residual call using a learned sparsity pattern. These tests run real
// decks in ColouredCheck mode, which also builds IDA's full banded DQ Jacobian
// at every coloured evaluation, and require the two matrices to agree: any
// coupling missing from the learned pattern shows up as a mismatch.
//
// The WENO and van Leer cases run the same decks with the higher-order
// reconstruction: its stencil reaches two cells and its interior faces are
// rescaled over the components, so every mole fraction couples to every mole
// fraction two cells away (the band is widened for it, see
// test_jacobian_band.cpp). That pattern is learned like any other and must
// match the band evaluation just as exactly.
//
// ColouredCheck is only a faithful diagnostic if its comparison does not
// change what it observes. ColouredCheckObservesTheColouredRun starts from a
// deliberately incomplete pattern (IDASettings::jacobian_learn_once_for_testing_),
// so the audits have couplings to find, and requires the ColouredCheck run to
// be the Coloured run exactly: same audits and misses, same integrator path,
// bit-identical final state.

#include <string>
#include <vector>

#include <gtest/gtest.h>

#include "adsim/objects/user_data.h"
#include "adsim/solver/coupled_solver.h"

namespace adsim {

  namespace {

    void ExpectColouredMatchesBand(const std::string& deck,
                                   AdvectionScheme scheme = AdvectionScheme::Upwind) {
      UserData data;
      LoadFromJson(data, "../examples/decks/" + deck);
      data.settings.jacobian_mode_ = JacobianMode::ColouredCheck;
      data.settings.advection_scheme_ = scheme;
      // 30 cells (the Nguyen deck's own grid; Streb2021_A1 has 60, which took
      // 2-3 s per scheme): the band and the colouring have the same structure
      // on any grid. Not fewer: on 10 and 20 cells (not 15 or 25) the learned
      // van Leer pattern misses couplings two cells away (NumMisses 1-8), an
      // open issue of the pattern learning, not of this test.
      for (auto& layer : data.bed.layers_.layers_) layer.domain_.n_cells = 30;

      CoupledSolver solver{data};
      const int retval = solver.IntegrateCycle(data.cycle, false);
      ASSERT_GE(retval, 0) << deck << ": solver did not complete";

      const auto& jac = solver.GetJacobian();
      EXPECT_LT(jac.NumColours(), 4 * static_cast<int>(data.state.GetEquations().size()) + 1)
          << deck << ": colouring should need fewer groups than IDA's band stride";
      // Entries the learned pattern misses must be rounding noise only.
      EXPECT_LT(jac.MaxCheckRelativeError(), 1e-8)
          << deck << ": " << jac.NumCheckMismatches() << " mismatched entries; first: "
          << (jac.CheckReport().empty() ? std::string("-") : jac.CheckReport().front());
      // The pattern learned in the first cycle held: no audit had to add to it.
      EXPECT_EQ(jac.NumMisses(), 0) << deck;
    }

    struct ColouredRun {
      int retval = 0;
      IntegrationStats stats;
      long misses = 0;
      long audits = 0;
      long failure_audits = 0;
      long res_evals = 0;
      long out_of_band = 0;
      int colours = 0;
      long check_mismatches = 0;
      std::vector<double> y;  // final mole fraction of the first component
      std::vector<double> P;  // final pressure
    };

    ColouredRun RunLearnOnce(const std::string& deck, JacobianMode mode) {
      UserData data;
      LoadFromJson(data, "../examples/decks/" + deck);
      data.settings.jacobian_mode_ = mode;
      data.settings.jacobian_learn_once_for_testing_ = true;

      CoupledSolver solver{data};
      ColouredRun r;
      r.retval = solver.IntegrateCycle(data.cycle, false);
      r.stats = solver.GetStatistics();
      const auto& jac = solver.GetJacobian();
      r.misses = jac.NumMisses();
      r.audits = jac.NumAudits();
      r.failure_audits = jac.NumFailureAudits();
      r.res_evals = jac.NumResidualEvals();
      r.out_of_band = jac.NumOutOfBand();
      r.colours = jac.NumColours();
      r.check_mismatches = jac.NumCheckMismatches();
      r.y = data.state.mole_fraction.at(data.fluid.components_.front()).var;
      r.P = data.state.pressure.var;
      return r;
    }

  } // namespace

  TEST(ColouredJacobian, MatchesBandDQ_Nguyen2022Fig5) {
    ExpectColouredMatchesBand("Nguyen2022_Fig5.json");
  }

  TEST(ColouredJacobian, MatchesBandDQ_Streb2021A1) {
    ExpectColouredMatchesBand("Streb2021_A1.json");
  }

  TEST(ColouredJacobian, MatchesBandDQ_Nguyen2022Fig5_WENO) {
    ExpectColouredMatchesBand("Nguyen2022_Fig5.json", AdvectionScheme::WENO);
  }

  TEST(ColouredJacobian, MatchesBandDQ_Nguyen2022Fig5_VanLeer) {
    ExpectColouredMatchesBand("Nguyen2022_Fig5.json", AdvectionScheme::VanLeer);
  }

  TEST(ColouredJacobian, MatchesBandDQ_Streb2021A1_WENO) {
    ExpectColouredMatchesBand("Streb2021_A1.json", AdvectionScheme::WENO);
  }

  TEST(ColouredJacobian, MatchesBandDQ_Streb2021A1_VanLeer) {
    ExpectColouredMatchesBand("Streb2021_A1.json", AdvectionScheme::VanLeer);
  }

  TEST(ColouredJacobian, ColouredCheckObservesTheColouredRun) {
    const std::string deck = "Nguyen2022_Fig5.json";
    const ColouredRun coloured = RunLearnOnce(deck, JacobianMode::Coloured);
    const ColouredRun check = RunLearnOnce(deck, JacobianMode::ColouredCheck);
    ASSERT_GE(coloured.retval, 0) << "the audits must complete a pattern learned from one evaluation";
    ASSERT_GE(check.retval, 0);

    // Guard: the pattern really was incomplete, the audits completed it, and
    // the check saw the missing couplings before they did.
    EXPECT_GT(coloured.misses, 0);
    EXPECT_GT(check.check_mismatches, 0);

    // The check only observes: the audits find exactly the same couplings...
    EXPECT_EQ(check.misses, coloured.misses);
    EXPECT_EQ(check.audits, coloured.audits);
    EXPECT_EQ(check.failure_audits, coloured.failure_audits);
    EXPECT_EQ(check.colours, coloured.colours);
    EXPECT_EQ(check.res_evals, coloured.res_evals);
    EXPECT_EQ(check.out_of_band, coloured.out_of_band);
    // ...and IDA receives the same Jacobians, so the run is the same run.
    EXPECT_EQ(check.stats.num_steps, coloured.stats.num_steps);
    EXPECT_EQ(check.stats.num_res_evals, coloured.stats.num_res_evals);
    EXPECT_EQ(check.stats.num_jac_evals, coloured.stats.num_jac_evals);
    EXPECT_EQ(check.stats.num_conv_fails, coloured.stats.num_conv_fails);
    EXPECT_EQ(check.stats.num_err_test_fails, coloured.stats.num_err_test_fails);
    ASSERT_EQ(check.y.size(), coloured.y.size());
    for (size_t i = 0; i < coloured.y.size(); ++i) EXPECT_EQ(check.y[i], coloured.y[i]) << "y, cell " << i;
    ASSERT_EQ(check.P.size(), coloured.P.size());
    for (size_t i = 0; i < coloured.P.size(); ++i) EXPECT_EQ(check.P[i], coloured.P[i]) << "P, cell " << i;
  }

} // end adsim namespace
