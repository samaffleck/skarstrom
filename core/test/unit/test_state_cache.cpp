// test_state_cache.cpp
//
// CycleStateCache (adsim/solver/state_cache.h): nearest-neighbour lookup in
// the bounds-scaled decision space, replacement of an equal key, eviction of
// the oldest entry at capacity, and the argument checks. The runs that use
// the cache are in integration/test_warm_start_cache.cpp.

#include <cmath>
#include <stdexcept>
#include <vector>

#include <gtest/gtest.h>

#include "adsim/solver/state_cache.h"

namespace adsim {

namespace {

// A state whose y[0] identifies it.
CycleState Tagged(double tag) {
  CycleState s;
  s.y = {tag};
  s.yp = {0.0};
  return s;
}

double TagOf(const CycleState* s) { return s ? s->y.at(0) : std::nan(""); }

}  // namespace

TEST(CycleStateCache, NearestIsInScaledVariablesNotTheMostRecent) {
  // Variable 0 spans 0..1000, variable 1 spans 0..1: unscaled, variable 0
  // would decide every distance.
  CycleStateCache cache({0.0, 0.0}, {1000.0, 1.0});
  EXPECT_EQ(cache.Dimension(), 2u);
  EXPECT_EQ(cache.Nearest({500.0, 0.5}), nullptr);

  cache.Add({100.0, 0.9}, Tagged(1));  // scaled (0.1, 0.9)
  cache.Add({900.0, 0.1}, Tagged(2));  // scaled (0.9, 0.1)
  cache.Add({500.0, 0.5}, Tagged(3));  // scaled (0.5, 0.5), the most recent

  double d = -1.0;
  // (200, 0.8) is scaled (0.2, 0.8): 0.14 from entry 1.
  EXPECT_EQ(TagOf(cache.Nearest({200.0, 0.8}, &d)), 1.0);
  EXPECT_NEAR(d, std::hypot(0.1, 0.1), 1e-12);
  EXPECT_EQ(TagOf(cache.Nearest({850.0, 0.0})), 2.0);
  // Raw distances would pick entry 3 for (450, 0.95): |450-500| < |450-100|.
  // Scaled: (0.45, 0.95) is 0.36 from entry 1 and 0.45 from entry 3.
  EXPECT_EQ(TagOf(cache.Nearest({450.0, 0.95})), 1.0);
  EXPECT_EQ(TagOf(cache.Nearest({500.0, 0.5}, &d)), 3.0);
  EXPECT_EQ(d, 0.0);

  const std::vector<double> key = cache.Normalise({250.0, 0.25});
  EXPECT_DOUBLE_EQ(key[0], 0.25);
  EXPECT_DOUBLE_EQ(key[1], 0.25);
}

TEST(CycleStateCache, WithoutBoundsKeysAreUsedAsGiven) {
  CycleStateCache cache;
  EXPECT_EQ(cache.Dimension(), 0u);
  cache.Add({0.2, 0.2, 0.2}, Tagged(1));
  cache.Add({0.8, 0.8, 0.8}, Tagged(2));
  EXPECT_EQ(cache.Dimension(), 3u);
  EXPECT_EQ(TagOf(cache.Nearest({0.4, 0.4, 0.4})), 1.0);
  EXPECT_EQ(TagOf(cache.Nearest({0.6, 0.6, 0.6})), 2.0);
  EXPECT_THROW(cache.Nearest({0.5, 0.5}), std::invalid_argument);
  EXPECT_THROW(cache.Add({0.5}, Tagged(3)), std::invalid_argument);
}

TEST(CycleStateCache, AnEqualKeyIsReplacedAndTiesGoToTheMostRecent) {
  CycleStateCache cache({0.0}, {1.0});
  cache.Add({0.5}, Tagged(1));
  cache.Add({0.5}, Tagged(2));
  EXPECT_EQ(cache.Size(), 1u);
  EXPECT_EQ(TagOf(cache.Nearest({0.5})), 2.0);

  cache.Add({0.3}, Tagged(3));
  cache.Add({0.7}, Tagged(4));
  EXPECT_EQ(TagOf(cache.Nearest({0.5})), 2.0);  // exact hit
  // 0.5 is equally far from 0.25 and 0.75: the newer one wins.
  CycleStateCache tie({0.0}, {1.0});
  tie.Add({0.25}, Tagged(5));
  tie.Add({0.75}, Tagged(6));
  EXPECT_EQ(TagOf(tie.Nearest({0.5})), 6.0);
}

TEST(CycleStateCache, TheOldestEntryGoesAtCapacity) {
  CycleStateCache cache({0.0}, {8.0}, 3);  // keys are exact binary fractions
  EXPECT_EQ(cache.Capacity(), 3u);
  cache.Add({1.0}, Tagged(1));
  cache.Add({2.0}, Tagged(2));
  cache.Add({3.0}, Tagged(3));
  cache.Add({1.0}, Tagged(11));  // replaces entry 1 and makes it the newest
  cache.Add({7.0}, Tagged(7));   // evicts entry 2, now the oldest
  EXPECT_EQ(cache.Size(), 3u);
  EXPECT_EQ(TagOf(cache.Nearest({2.0})), 11.0);  // 1 and 3 are both 1 away; 11 is newer
  EXPECT_EQ(TagOf(cache.Nearest({1.0})), 11.0);
  EXPECT_EQ(TagOf(cache.Nearest({7.5})), 7.0);
  cache.Clear();
  EXPECT_TRUE(cache.Empty());
}

TEST(CycleStateCache, BadArgumentsAreRejected) {
  EXPECT_THROW(CycleStateCache({0.0, 0.0}, {1.0}), std::invalid_argument);
  EXPECT_THROW(CycleStateCache({1.0}, {1.0}), std::invalid_argument);
  EXPECT_THROW(CycleStateCache({0.0}, {std::nan("")}), std::invalid_argument);
  CycleStateCache cache({0.0}, {1.0});
  EXPECT_THROW(cache.Add({0.5}, CycleState{}), std::invalid_argument);
  EXPECT_THROW(cache.Add({std::nan("")}, Tagged(1)), std::invalid_argument);
  EXPECT_THROW(cache.Add({0.5, 0.5}, Tagged(1)), std::invalid_argument);
}

}  // namespace adsim
