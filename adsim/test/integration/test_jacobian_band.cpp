// test_jacobian_band.cpp
//
// The Newton matrix is a SUNDIALS band matrix with half-bandwidths mu = ml =
// 2*neq (neq = number of equation types). Any residual dependency outside the
// band is not merely dropped: band DQ perturbs columns j, j + (mu+ml+1), ...
// together, so an out-of-band coupling of column j is mis-attributed to another
// column of the same group. These tests build a dense finite-difference
// dependency map of the residual (one column per residual call, at a randomised
// generic state with both flow directions) and require every nonzero to lie
// inside the band the solver actually allocates, for every advection scheme.
// WENO and van Leer faces reach two cells away, but only within the same mole
// fraction, so their reach in the cell-major layout is exactly 2*neq.

#include <algorithm>
#include <cmath>
#include <random>
#include <string>
#include <vector>

#include <gtest/gtest.h>

#include "adsim/objects/user_data.h"
#include "adsim/solver/coupled_solver.h"

namespace adsim {

  namespace {

    struct BandCheck {
      long lower_reach = 0;  // max (row - col) over residual dependencies
      long upper_reach = 0;  // max (col - row) over residual dependencies
      std::string lower_at, upper_at;
      sunindextype mu = 0;   // solver's upper half-bandwidth
      sunindextype ml = 0;   // solver's lower half-bandwidth
    };

    BandCheck MeasureBand(const std::string& deck, AdvectionScheme scheme) {
      BandCheck out;
      UserData data;
      LoadFromJson(data, "../examples/decks/" + deck);
      data.settings.advection_scheme_ = scheme;
      data.settings.jacobian_mode_ = JacobianMode::Coloured;
      std::string err;
      EXPECT_TRUE(data.PerformChecks(err)) << deck << ": " << err;

      // Same set-up as CoupledSolver::Initialise.
      data.Initialise();
      Step& first = data.cycle.step_list_[0];
      data.state.SetCurrentStep(&first);
      data.step_start_time_ = 0.0;
      first.Initialise(data.state);

      SolverData solver;
      const int flag = solver.Initialise(data, err);
      EXPECT_EQ(flag, IDA_SUCCESS) << deck << ": " << err;
      out.mu = solver.GetJacobian().UpperBandwidth();
      out.ml = solver.GetJacobian().LowerBandwidth();

      auto& state = data.state;
      const int n = state.GetSystemSize();
      const auto& eqs = state.GetEquations();
      const int neq = static_cast<int>(eqs.size());

      // Equation type of each banded entry (same loop order as PackBanded).
      size_t n_cells = 0;
      for (const auto* eq : eqs) n_cells = std::max(n_cells, eq->Size());
      std::vector<int> var;
      for (size_t c = 0; c < n_cells; ++c)
        for (int e = 0; e < neq; ++e)
          if (eqs[e]->IsActiveAt(c, &state.domain_vectors.is_void_cell)) var.push_back(e);
      EXPECT_EQ(static_cast<int>(var.size()), n);

      SUNContext ctx = nullptr;
      SUNContext_Create(SUN_COMM_NULL, &ctx);
      N_Vector y = N_VNew_Serial(n, ctx);
      N_Vector yp = N_VNew_Serial(n, ctx);
      N_Vector r0 = N_VNew_Serial(n, ctx);
      N_Vector r1 = N_VNew_Serial(n, ctx);
      double* yd = N_VGetArrayPointer(y);
      double* ypd = N_VGetArrayPointer(yp);
      const double* r0d = N_VGetArrayPointer(r0);
      const double* r1d = N_VGetArrayPointer(r1);
      state.PackBanded(yd, ypd);

      // Generic state: non-uniform mole fractions (so limiter / smoothness
      // weights have nonzero derivatives) and face velocities of random sign
      // (so both the left- and right-biased reconstructions are exercised).
      constexpr int kVelocity = 5;  // equation index of u (Equations::Initialise order)
      constexpr int kFirstMoleFraction = 6;
      const int n_comp = static_cast<int>(data.fluid.components_.size());
      std::mt19937_64 rng(0xBA9D);
      std::uniform_real_distribution<double> jitter(-1.0, 1.0);
      for (int j = 0; j < n; ++j) {
        const int e = var[j];
        if (e == kVelocity) yd[j] = (std::abs(yd[j]) + 0.05) * jitter(rng);
        else if (e >= kFirstMoleFraction && e < kFirstMoleFraction + n_comp) yd[j] = 0.5 + 0.4 * jitter(rng);
        else yd[j] *= 1.0 + 1e-2 * jitter(rng);
        ypd[j] = (std::abs(ypd[j]) + 1e-3 * (std::abs(yd[j]) + 1e-3)) * jitter(rng);
      }

      UpdateResiduals(0.0, y, yp, r0, &data);
      for (int i = 0; i < n; ++i)
        EXPECT_TRUE(std::isfinite(r0d[i])) << deck << ": non-finite residual at " << state.DescribeBandedIndex(i);

      // Dense FD dependency map over both y and yp (IDA's Jacobian is dF/dy + c_j dF/dyp).
      for (double* v : { yd, ypd }) {
        for (int j = 0; j < n; ++j) {
          const double saved = v[j];
          v[j] = saved + 1e-6 * (std::abs(saved) + 1e-6);
          UpdateResiduals(0.0, y, yp, r1, &data);
          v[j] = saved;
          for (int i = 0; i < n; ++i) {
            if (r1d[i] == r0d[i]) continue;
            if (i - j > out.lower_reach) {
              out.lower_reach = i - j;
              out.lower_at = state.DescribeBandedIndex(i) + " <- " + state.DescribeBandedIndex(j);
            }
            if (j - i > out.upper_reach) {
              out.upper_reach = j - i;
              out.upper_at = state.DescribeBandedIndex(i) + " <- " + state.DescribeBandedIndex(j);
            }
          }
        }
      }

      N_VDestroy(y);
      N_VDestroy(yp);
      N_VDestroy(r0);
      N_VDestroy(r1);
      SUNContext_Free(&ctx);
      return out;
    }

    void ExpectResidualWithinBand(const std::string& deck) {
      const BandCheck upwind = MeasureBand(deck, AdvectionScheme::Upwind);
      for (auto scheme : { AdvectionScheme::Upwind, AdvectionScheme::WENO, AdvectionScheme::VanLeer }) {
        const BandCheck b = MeasureBand(deck, scheme);
        const int s = static_cast<int>(scheme);
        EXPECT_LE(b.lower_reach, b.ml) << deck << " scheme " << s << ": row-col reach "
                                       << b.lower_reach << " (" << b.lower_at << ") exceeds ml = " << b.ml;
        EXPECT_LE(b.upper_reach, b.mu) << deck << " scheme " << s << ": col-row reach "
                                       << b.upper_reach << " (" << b.upper_at << ") exceeds mu = " << b.mu;
        if (scheme != AdvectionScheme::Upwind) {
          // The probe must actually see the two-cell stencil, otherwise the
          // check above is vacuous.
          EXPECT_GT(b.lower_reach, upwind.lower_reach) << deck << " scheme " << s;
          EXPECT_GT(b.upper_reach, upwind.upper_reach) << deck << " scheme " << s;
        }
      }
    }

  } // namespace

  TEST(JacobianBand, ResidualWithinBand_Nguyen2022Fig5) {
    ExpectResidualWithinBand("Nguyen_2022_Fig5.json");  // 3 components
  }

  TEST(JacobianBand, ResidualWithinBand_Arora2021PSA) {
    ExpectResidualWithinBand("Arora2021_LiLSX_PSA.json");  // 2 components
  }

} // end adsim namespace
