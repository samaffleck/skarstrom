// STL includes
#include <filesystem>

// ads includes
#include "adsim/math/linear_algebra.h"

// Gtest includes
#include <gtest/gtest.h>

namespace adsim {

  TEST(LinearAlgebra, TDMA_SolvesTridiagonalSystem)
  {
    constexpr int n = 5;
    Eigen::MatrixXd A = Eigen::MatrixXd::Zero(n, n);

    for (int i = 0; i < n; ++i)
    {
      A(i, i) = 2.0;
      if (i > 0)     A(i, i - 1) = 1.0;
      if (i < n - 1) A(i, i + 1) = -1.0;
    }

    Eigen::VectorXd x_expected(n);
    x_expected << 1, 2, 3, 4, 5;

    Eigen::VectorXd S = A * x_expected;

    Eigen::VectorXd x = Eigen::VectorXd::Zero(n);
    TDMA(A, S, x);

    for (int i = 0; i < n; ++i)
      EXPECT_NEAR(x(i), x_expected(i), 1e-9);
  }

  TEST(NonlinearSolver, ApproximateJacobian_MatchesAnalytical)
  {
    constexpr int n = 5;
    Eigen::VectorXd x(n);
    x << 0.5, -1.1, 2.0, -0.3, 4.2;

    auto fx = [](const Eigen::VectorXd& in, Eigen::VectorXd& out)
    {
      out = in.array().square().matrix();
    };

    Eigen::VectorXd y(n);
    fx(x, y);

    Eigen::MatrixXd J = Eigen::MatrixXd::Zero(n, n);
    ApproximateJacobianForward(fx, J, x, y);

    Eigen::MatrixXd J_expected = Eigen::MatrixXd::Zero(n, n);
    for (int i = 0; i < n; ++i) J_expected(i, i) = 2.0 * x(i);

    for (int i = 0; i < n; ++i)
    {
      for (int j = std::max(0, i - 1); j <= std::min(n - 1, i + 1); ++j)
      {
        EXPECT_NEAR(J(i, j), J_expected(i, j), 1e-5);
      }
    }
  }

  TEST(NonlinearSolver, NewtonStep_UpdatesSolutionCorrectly)
  {
    constexpr int n = 5;
    Eigen::MatrixXd J = Eigen::MatrixXd::Zero(n, n);

    for (int i = 0; i < n; ++i)
    {
      J(i, i) = 2.0;
      if (i > 0)     J(i, i - 1) = 1.0;
      if (i < n - 1) J(i, i + 1) = -1.0;
    }

    Eigen::VectorXd x_target(n);
    x_target << 1, 2, 3, 4, 5;

    Eigen::VectorXd R = -J * x_target;

    Eigen::VectorXd x = Eigen::VectorXd::Zero(n);

    NewtonStep(J, R, x, /*relaxFactor=*/1.0);

    ASSERT_TRUE(x.isApprox(x_target, 1e-12));
  }

} // end adsim namespace
