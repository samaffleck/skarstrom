// test_jacobian_dense_check.cpp
//
// Checks the Jacobian the integrator actually uses against a dense
// finite-difference Jacobian. test_coloured_jacobian.cpp compares the coloured
// evaluation with IDA's banded difference quotient, and test_jacobian_band.cpp
// measures how far the residual couplings reach, at a jittered initial state.
// Neither compares the matrix values with a reference that makes no
// assumption: the band DQ itself perturbs columns j, j + (mu + ml + 1), ...
// together, so it shares any error caused by a coupling outside the band.
//
// Here the decks run normally (Coloured mode, the default), and at states in
// the middle of each step (20, 50 and 80 % of the step time, in the last cycle
// run) the matrix IDA receives is compared, entry by entry over the whole
// n x n matrix, with a dense forward difference: every column perturbed on its
// own, no colouring, no band. Each column uses the increment the solver used
// (EvaluationView::increments), so the two agree to rounding wherever the
// matrix is right; a coupling the learned pattern misses, or one outside the
// band, shows up as an entry of order one relative to its column.
//
// The cases cover the three advection schemes, thermal runs (fluid and wall
// energy balances) and isothermal ones, Kozeny-Carman and Ergun momentum, a
// breakthrough, a PSA and a TSA cycle (states in the second cycle, after the
// pattern was learned), and each deck's own grid as well as a coarse 10-cell
// grid, including the van Leer / WENO cases on 10 and 20 cells where the
// learned pattern used to miss couplings two cells away (GitHub issue #77).
// Every checked evaluation (coloured, and the band evaluations that learn or
// audit the pattern) must match. DenseCheckDetectsAnIncompletePattern shows
// that the comparison does see a missing coupling.

#include <algorithm>
#include <chrono>
#include <cmath>
#include <cstdio>
#include <sstream>
#include <string>
#include <vector>

#include <gtest/gtest.h>
#include <sunmatrix/sunmatrix_band.h>
#include <nvector/nvector_serial.h>

#include "adsim/objects/user_data.h"
#include "adsim/solver/coupled_solver.h"
#include "adsim/solver/jacobian.h"

namespace adsim {

  namespace {

    // Entries whose difference from the dense reference exceeds this fraction
    // of their column's largest dense entry are mismatches. Agreeing entries
    // come out of the same arithmetic up to rounding (observed < 1e-10).
    constexpr double kRelTol = 1e-8;

    struct DenseCase {
      std::string deck;
      AdvectionScheme scheme = AdvectionScheme::Upwind;
      int cells = 0;            // 0: the deck's own grid
      bool isothermal = false;  // force the isothermal energy balance
      int cycles = 1;           // cycles to run; states are taken in the last one
      bool learn_once = false;  // IDASettings::jacobian_learn_once_for_testing_
      int first_coloured = 0;   // also check the first N coloured evaluations of the run
    };

    struct DenseResult {
      int retval = 0;
      int n = 0;
      int checks = 0;            // evaluations compared
      int coloured_checks = 0;   // of which coloured (the rest were band evaluations)
      long mismatches = 0;       // entries above kRelTol (coloured + band evaluations)
      long band_mismatches = 0;  // of which in band evaluations
      long outside_band = 0;     // dense nonzeros the band cannot hold
      double max_rel = 0.0;
      double max_rel_coloured = 0.0;
      double max_rel_band = 0.0;
      std::vector<std::string> report;
      std::vector<std::string> variables;  // equation types of the run
      long misses = 0;           // NumMisses() of the run: couplings the audits added
      bool ergun = false;        // every layer uses Ergun momentum
    };

    std::string SchemeName(AdvectionScheme s) {
      switch (s) {
        case AdvectionScheme::Upwind: return "upwind";
        case AdvectionScheme::WENO: return "WENO";
        case AdvectionScheme::VanLeer: return "van Leer";
      }
      return "?";
    }

    std::string Describe(const DenseCase& c) {
      std::ostringstream os;
      os << c.deck << ", " << SchemeName(c.scheme) << ", "
         << (c.cells > 0 ? std::to_string(c.cells) + " cells" : std::string("deck grid"))
         << (c.isothermal ? ", isothermal" : "");
      return os.str();
    }

    // Compares the matrix of one evaluation with the dense forward difference.
    void CompareWithDense(const BandedJacobian::EvaluationView& v, DenseResult& out) {
      UserData* udata = v.udata;
      const sunindextype n = N_VGetLength(v.yy);
      N_Vector y = N_VClone(v.yy);
      N_Vector yp = N_VClone(v.yp);
      N_Vector r0 = N_VClone(v.yy);
      N_Vector r1 = N_VClone(v.yy);
      N_VScale(1.0, v.yy, y);
      N_VScale(1.0, v.yp, yp);
      double* yd = N_VGetArrayPointer(y);
      double* ypd = N_VGetArrayPointer(yp);
      const double* r0d = N_VGetArrayPointer(r0);
      const double* r1d = N_VGetArrayPointer(r1);
      const std::string residual_error = udata->residual_error_;

      bool ok = UpdateResiduals(v.tt, y, yp, r0, udata) == 0;
      std::vector<double> column(n);
      const bool coloured = v.kind == BandedJacobian::EvaluationView::Kind::Coloured;
      for (sunindextype j = 0; ok && j < n; ++j) {
        const double inc = (*v.increments)[j];
        const double y0 = yd[j], yp0 = ypd[j];
        yd[j] += inc;
        ypd[j] += v.c_j * inc;
        ok = UpdateResiduals(v.tt, y, yp, r1, udata) == 0;
        yd[j] = y0;
        ypd[j] = yp0;
        if (!ok) break;
        double col_max = 0.0;
        for (sunindextype i = 0; i < n; ++i) {
          column[i] = (r1d[i] - r0d[i]) / inc;
          col_max = std::max(col_max, std::abs(column[i]));
        }
        if (col_max == 0.0) continue;
        for (sunindextype i = 0; i < n; ++i) {
          const bool in_band = i >= j - v.mu && i <= j + v.ml;
          if (!in_band && column[i] != 0.0) ++out.outside_band;
          const double used = in_band ? SM_ELEMENT_B(v.jac, i, j) : 0.0;
          const double rel = std::abs(used - column[i]) / col_max;
          out.max_rel = std::max(out.max_rel, rel);
          (coloured ? out.max_rel_coloured : out.max_rel_band) =
              std::max(coloured ? out.max_rel_coloured : out.max_rel_band, rel);
          if (rel <= kRelTol) continue;
          ++out.mismatches;
          if (!coloured) ++out.band_mismatches;
          if (out.report.size() < 12) {
            char buf[320];
            std::snprintf(buf, sizeof(buf), "%s t=%.4g: row %s col %s used=%.4e dense=%.4e rel=%.2e%s",
                          coloured ? "coloured" : "band", v.tt,
                          udata->state.DescribeBandedIndex(static_cast<int>(i)).c_str(),
                          udata->state.DescribeBandedIndex(static_cast<int>(j)).c_str(),
                          used, column[i], rel, in_band ? "" : " (outside band)");
            out.report.push_back(buf);
          }
        }
      }
      EXPECT_TRUE(ok) << "residual failed during the dense difference quotient: " << udata->residual_error_;
      udata->residual_error_ = residual_error;
      ++out.checks;
      if (coloured) ++out.coloured_checks;

      // Leave the state where IDA left it (its next residual call overwrites it anyway).
      if (ok) UpdateResiduals(v.tt, y, yp, r0, udata);
      udata->residual_error_ = residual_error;
      N_VDestroy(y);
      N_VDestroy(yp);
      N_VDestroy(r0);
      N_VDestroy(r1);
    }

    DenseResult RunDenseCheck(const DenseCase& c) {
      DenseResult out;
      UserData data;
      LoadFromJson(data, "../examples/decks/" + c.deck);
      data.settings.jacobian_mode_ = JacobianMode::Coloured;
      data.settings.advection_scheme_ = c.scheme;
      if (c.cells > 0)
        for (auto& layer : data.bed.layers_.layers_) layer.domain_.n_cells = c.cells;
      if (c.isothermal) data.bed.equations_.energy_balance_mode_ = EnergyBalanceMode::ISOTHERMAL;
      data.cycle.number_of_cycles_ = c.cycles;
      data.settings.jacobian_learn_once_for_testing_ = c.learn_once;
      out.ergun = !data.bed.layers_.layers_.empty();
      for (const auto& layer : data.bed.layers_.layers_)
        out.ergun = out.ergun && layer.momentum_model_ == MomentumModel::Ergun;

      // Step durations (resolved expressions) to place the checks mid-step.
      std::vector<double> durations;
      for (auto& step : data.cycle.step_list_) durations.push_back(static_cast<double>(step.step_time_));
      const int n_steps = static_cast<int>(durations.size());
      const int last_cycle = c.cycles - 1;
      const double fractions[] = { 0.2, 0.5, 0.8 };
      constexpr int kNumFractions = 3;

      // Step boundaries are recognised by the step start time changing.
      double last_start = -1.0;
      int step_counter = -1;
      int next_fraction = 0;
      int coloured_seen = 0;
      bool want_coloured = false;  // a sample point checked a band evaluation
      BandedJacobian::SetEvaluationObserverForTesting([&](const BandedJacobian::EvaluationView& v) {
        if (v.udata->step_start_time_ != last_start || step_counter < 0) {
          last_start = v.udata->step_start_time_;
          ++step_counter;
          next_fraction = 0;
        }
        const bool is_coloured = v.kind == BandedJacobian::EvaluationView::Kind::Coloured;
        if (is_coloured && coloured_seen++ < c.first_coloured) {
          CompareWithDense(v, out);
          return;
        }
        // A sample point whose evaluation was a band one (the pattern being
        // learned or audited) also checks the next coloured evaluation, so
        // every case checks the matrix IDA usually receives whatever the
        // solver's path (which evaluations fall after the sample points
        // differs between platforms).
        if (is_coloured && want_coloured) {
          want_coloured = false;
          CompareWithDense(v, out);
          return;
        }
        if (n_steps == 0 || step_counter / n_steps != last_cycle) return;
        const double duration = durations[step_counter % n_steps];
        const double t_rel = v.tt - v.udata->step_start_time_;
        if (next_fraction >= kNumFractions || t_rel < fractions[next_fraction] * duration) return;
        while (next_fraction < kNumFractions && t_rel >= fractions[next_fraction] * duration) ++next_fraction;
        CompareWithDense(v, out);
        if (!is_coloured) want_coloured = true;
      });

      CoupledSolver solver{data};
      out.retval = solver.IntegrateCycle(data.cycle, false);
      BandedJacobian::SetEvaluationObserverForTesting({});
      out.n = data.state.GetSystemSize();
      out.misses = solver.GetJacobian().NumMisses();
      const auto& eqs = data.state.GetEquations();
      const std::pair<const Equation*, const char*> named[] = {
        { &data.state.pressure, "P" }, { &data.state.fluid_temperature, "Tf" },
        { &data.state.solid_temperature, "Ts" }, { &data.state.wall_temperature, "Tw" },
        { &data.state.velocity, "u" } };
      for (const auto* eq : eqs)
        for (const auto& [ptr, name] : named)
          if (eq == ptr) out.variables.push_back(name);
      return out;
    }

    bool Has(const std::vector<std::string>& v, const std::string& s) {
      return std::find(v.begin(), v.end(), s) != v.end();
    }

    std::string Summary(const DenseCase& c, const DenseResult& r) {
      std::ostringstream os;
      os << Describe(c) << ": n = " << r.n << ", " << r.checks << " states (" << r.coloured_checks
         << " coloured), max rel error " << r.max_rel << " (coloured " << r.max_rel_coloured
         << ", band " << r.max_rel_band << "), " << r.mismatches << " mismatched entries ("
         << r.band_mismatches << " in band evaluations), " << r.outside_band
         << " dense couplings outside the band, " << r.misses << " pattern misses added by audits";
      for (const auto& line : r.report) os << "\n    " << line;
      return os.str();
    }

    // The matrix IDA receives must equal the dense difference quotient at
    // every checked state, over the whole matrix.
    void ExpectMatchesDense(const DenseCase& c) {
      const DenseResult r = RunDenseCheck(c);
      std::printf("[dense jacobian] %s\n", Summary(c, r).c_str());
      ASSERT_GE(r.retval, 0) << Describe(c) << ": solver did not complete";
      ASSERT_GT(r.coloured_checks, 0) << Describe(c) << ": no coloured evaluation was checked";
      EXPECT_EQ(r.outside_band, 0) << Summary(c, r);
      EXPECT_LT(r.max_rel, kRelTol) << Summary(c, r);
      EXPECT_EQ(r.mismatches, 0) << Summary(c, r);
    }

    struct Expect {
      bool thermal = true;  // fluid and wall energy balances are solved
      bool ergun = false;
    };

    // The matrix IDA receives must equal the dense difference quotient at
    // every checked state, over the whole matrix.
    void ExpectMatchesDense(const DenseCase& c, Expect expect) {
      const DenseResult r = RunDenseCheck(c);
      std::printf("[dense jacobian] %s\n", Summary(c, r).c_str());
      ASSERT_GE(r.retval, 0) << Describe(c) << ": solver did not complete";
      ASSERT_GT(r.coloured_checks, 0) << Describe(c) << ": no coloured evaluation was checked";
      // The case covers what it is meant to cover.
      EXPECT_TRUE(Has(r.variables, "P") && Has(r.variables, "u")) << Describe(c);
      EXPECT_EQ(Has(r.variables, "Tf"), expect.thermal) << Describe(c);
      EXPECT_EQ(Has(r.variables, "Tw"), expect.thermal) << Describe(c);
      EXPECT_EQ(r.ergun, expect.ergun) << Describe(c);

      EXPECT_EQ(r.outside_band, 0) << Summary(c, r);
      EXPECT_LT(r.max_rel, kRelTol) << Summary(c, r);
      EXPECT_EQ(r.mismatches, 0) << Summary(c, r);
      EXPECT_EQ(r.misses, 0) << Summary(c, r);
    }

    constexpr auto kUpwind = AdvectionScheme::Upwind;
    constexpr auto kWENO = AdvectionScheme::WENO;
    constexpr auto kVanLeer = AdvectionScheme::VanLeer;
    constexpr AdvectionScheme kSchemes[] = { kUpwind, kWENO, kVanLeer };

  } // namespace

  // Breakthrough, three components, Kozeny-Carman, fluid + wall energy
  // balances: every scheme on the deck's grid (30 cells), 10 and 20 cells.
  TEST(JacobianDenseCheck, Nguyen2022Fig5_AllSchemesAndGrids) {
    for (auto scheme : kSchemes)
      for (int cells : { 0, 10, 20 })
        ExpectMatchesDense({ "Nguyen2022_Fig5.json", scheme, cells }, {});
  }

  TEST(JacobianDenseCheck, Nguyen2022Fig5_Isothermal) {
    for (auto scheme : kSchemes)
      for (int cells : { 0, 10 })
        ExpectMatchesDense({ "Nguyen2022_Fig5.json", scheme, cells, true }, { .thermal = false });
  }

  // Ergun momentum (the pressure couples through the nonlinear Ergun root),
  // fixed wall temperature boundary.
  TEST(JacobianDenseCheck, Ward2022Q100ccm_Ergun) {
    for (auto scheme : kSchemes)
      for (int cells : { 0, 10 })
        ExpectMatchesDense({ "Ward2022_Q100ccm.json", scheme, cells }, { .ergun = true });
  }

  // Four-step PSA cycle; states in the second cycle (pressurisation, feed,
  // blowdown, purge), after the pattern was learned in the first.
  TEST(JacobianDenseCheck, Perez2019A_PsaCycle) {
    ExpectMatchesDense({ "Perez2019_A.json", kUpwind, 0, false, 2 }, {});
    ExpectMatchesDense({ "Perez2019_A.json", kVanLeer, 10, false, 2 }, {});
    ExpectMatchesDense({ "Perez2019_A.json", kWENO, 10, true, 2 }, { .thermal = false });
  }

  // TSA cycle with Ergun momentum and heating; second cycle.
  TEST(JacobianDenseCheck, IllustrativeTsa13X_ErgunTsaCycle) {
    for (auto scheme : kSchemes)
      ExpectMatchesDense({ "Illustrative_TSA_13X.json", scheme, 0, false, 2 }, { .ergun = true });
    ExpectMatchesDense({ "Illustrative_TSA_13X.json", kVanLeer, 10, false, 2 }, { .ergun = true });
  }

  // The coarse-grid higher-order cases of issue #77 (Streb2021_A1, van Leer
  // and WENO on 10 and 20 cells, also over three cycles).
  TEST(JacobianDenseCheck, Streb2021A1_CoarseHigherOrder) {
    for (auto scheme : { kWENO, kVanLeer })
      for (int cells : { 10, 20 })
        ExpectMatchesDense({ "Streb2021_A1.json", scheme, cells, false, 3 }, {});
  }

  // The check has teeth: a pattern learned from a single evaluation is
  // incomplete until the audits complete it, and the dense comparison must
  // see the coloured entries it misses (and nothing wrong in the band
  // evaluations or outside the band).
  TEST(JacobianDenseCheck, DenseCheckDetectsAnIncompletePattern) {
    DenseCase c{ "Nguyen2022_Fig5.json", kUpwind, 0 };
    c.learn_once = true;
    c.first_coloured = 5;  // the audits complete the pattern long before mid-step
    const DenseResult r = RunDenseCheck(c);
    std::printf("[dense jacobian, learn once] %s\n", Summary(c, r).c_str());
    ASSERT_GE(r.retval, 0);
    ASSERT_GT(r.misses, 0) << "the guard needs an incomplete pattern";
    EXPECT_GT(r.mismatches, 0) << Summary(c, r);
    EXPECT_GT(r.max_rel_coloured, 1e-3) << Summary(c, r);
    EXPECT_EQ(r.band_mismatches, 0) << Summary(c, r);
    EXPECT_EQ(r.outside_band, 0) << Summary(c, r);
  }

} // end adsim namespace
