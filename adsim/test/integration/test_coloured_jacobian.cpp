// test_coloured_jacobian.cpp
//
// The coloured DQ Jacobian (JacobianMode::Coloured) perturbs several columns
// per residual call using a learned sparsity pattern. These tests run real
// decks in ColouredCheck mode, which also builds IDA's full banded DQ Jacobian
// at every coloured evaluation, and require the two matrices to agree: any
// coupling missing from the learned pattern shows up as a mismatch.

#include <string>

#include <gtest/gtest.h>

#include "adsim/objects/user_data.h"
#include "adsim/solver/coupled_solver.h"

namespace adsim {

  namespace {

    void ExpectColouredMatchesBand(const std::string& deck) {
      UserData data;
      LoadFromJson(data, "../examples/decks/" + deck);
      data.settings.jacobian_mode_ = JacobianMode::ColouredCheck;

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
    }

  } // namespace

  TEST(ColouredJacobian, MatchesBandDQ_Nguyen2022Fig5) {
    ExpectColouredMatchesBand("Nguyen_2022_Fig5.json");
  }

  TEST(ColouredJacobian, MatchesBandDQ_Streb2021A1) {
    ExpectColouredMatchesBand("streb_2021_A1.json");
  }

} // end adsim namespace
