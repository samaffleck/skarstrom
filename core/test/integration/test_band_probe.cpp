// test_band_probe.cpp
//
// The banded Jacobian matrix is sized from the coupling reach a difference-
// quotient probe measures at the initial state (IDASettings::band_from_probe_)
// instead of the 2 neq each side that covers every scheme. These tests check
// that the probed band is narrower, that it reproduces the full-band results
// exactly, and that a band found too narrow by an audit widens itself.

#include <cmath>
#include <string>
#include <vector>

#include <gtest/gtest.h>

#include "adsim/objects/user_data.h"
#include "adsim/solver/coupled_solver.h"

namespace adsim {

  namespace {

    struct BandRun {
      int retval = 0;
      sunindextype mu = 0;
      sunindextype ml = 0;
      bool widened = false;
      long out_of_band = 0;
      std::vector<double> y;  // final mole fraction of the first component
    };

    BandRun RunDeck(const std::string& deck, bool probe, int shrink) {
      UserData data;
      LoadFromJson(data, "../examples/decks/" + deck);
      data.settings.band_from_probe_ = probe;
      data.settings.band_probe_shrink_for_testing_ = shrink;
      CoupledSolver solver{data};
      BandRun r;
      r.retval = solver.IntegrateCycle(data.cycle, false);
      r.mu = solver.GetSolverData().BandUpper();
      r.ml = solver.GetSolverData().BandLower();
      r.widened = solver.GetSolverData().BandWidened();
      r.out_of_band = solver.GetJacobian().NumOutOfBand();
      r.y = data.state.mole_fraction.at(data.fluid.components_.front()).var;
      return r;
    }

  } // namespace

  TEST(BandProbe, ProbedBandIsNarrowerAndReproducesTheFullBand) {
    const BandRun full = RunDeck("Nguyen2022_Fig5.json", false, 0);
    const BandRun probed = RunDeck("Nguyen2022_Fig5.json", true, 0);
    ASSERT_GE(full.retval, 0);
    ASSERT_GE(probed.retval, 0);
    EXPECT_LT(probed.mu, full.mu);
    EXPECT_LT(probed.ml, full.ml);
    EXPECT_EQ(probed.out_of_band, 0);
    EXPECT_FALSE(probed.widened);
    // The band LU performs the same operations on the same nonzeros, so the
    // trajectory is bit for bit the one of the full band.
    ASSERT_EQ(probed.y.size(), full.y.size());
    for (size_t i = 0; i < full.y.size(); ++i) EXPECT_EQ(probed.y[i], full.y[i]) << "cell " << i;
  }

  TEST(BandProbe, TooNarrowBandWidensItselfAndStillConverges) {
    const BandRun full = RunDeck("Nguyen2022_Fig5.json", false, 0);
    const BandRun narrow = RunDeck("Nguyen2022_Fig5.json", true, 4);
    ASSERT_GE(full.retval, 0);
    ASSERT_GE(narrow.retval, 0) << "the run must survive a band that was probed too narrow";
    EXPECT_TRUE(narrow.widened) << "an audit should have found couplings outside the band";
    EXPECT_EQ(narrow.mu, full.mu);
    EXPECT_EQ(narrow.ml, full.ml);
    // The restart changes the integration path, so compare within tolerance.
    ASSERT_EQ(narrow.y.size(), full.y.size());
    for (size_t i = 0; i < full.y.size(); ++i)
      EXPECT_NEAR(narrow.y[i], full.y[i], 1e-3) << "cell " << i;
  }

} // end adsim namespace
