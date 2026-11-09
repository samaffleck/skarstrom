// test_chunked_statistics.cpp
//
// The chunked solver (the one a front end runs in slices) reports the same
// integrator counters as a blocking run, through GetStatistics().

#include <string>

#include "../src/pch.h"
#include "adsim/solver/coupled_solver.h"

#include <gtest/gtest.h>

#include "../fixtures/user_data_test.h"

namespace adsim {

  class ChunkedStatisticsTest : public UserDataTest {
  protected:
    void SetUp() override {
      UserDataTest::SetUp();
      AddFeedHoldDeck(2);
    }
  };

  TEST_F(ChunkedStatisticsTest, MatchesBlockingRun) {
    std::string error;
    ASSERT_TRUE(data.PerformChecks(error)) << error;

    CoupledSolver blocking{data};
    ASSERT_GE(blocking.IntegrateCycle(data.cycle, false), 0);
    const auto b = blocking.GetStatistics();
    ASSERT_GT(b.num_steps, 0);

    CoupledSolverChuncked chunked{data};
    chunked.Start();
    ASSERT_TRUE(chunked.IsActive()) << chunked.GetStatus().GetMessage();
    int chunks = 0;
    while (chunked.IsActive() && chunks < 20000) {
      chunked.Integrate(5);
      ++chunks;
    }
    ASSERT_FALSE(chunked.IsActive());
    const auto c = chunked.GetStatistics();

    EXPECT_EQ(c.num_steps, b.num_steps);
    EXPECT_EQ(c.num_res_evals, b.num_res_evals);
    EXPECT_EQ(c.num_nonlin_iters, b.num_nonlin_iters);
    EXPECT_EQ(c.num_err_test_fails, b.num_err_test_fails);
    EXPECT_EQ(c.num_conv_fails, b.num_conv_fails);
    EXPECT_EQ(c.num_jac_evals, b.num_jac_evals);
    EXPECT_EQ(c.num_lin_setups, b.num_lin_setups);
    EXPECT_EQ(c.recovered_steps, b.recovered_steps);
  }

} // namespace adsim
