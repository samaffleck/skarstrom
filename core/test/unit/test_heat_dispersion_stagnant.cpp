// Stagnant (no-flow) effective conductivity of a packed bed
// (dispersion_model.h): the Zehner-Schlunder ratio k0/k_g and the
// Krupiczka-type ratio of the yagi_kunii kind. These tests pin:
//   * the limits: kappa = 1 gives k0 = k_g, eps -> 1 gives k0 -> k_g;
//   * monotonicity: k0 rises with kappa and falls with eps;
//   * the removable singularity at kappa = B (N = 0): finite, equal to the
//     series limit, continuous across the switch between the series and the
//     closed form at |N| = 0.1;
//   * a few values. These are computed from the closed form in double
//     precision (an independent Python evaluation), not literature values.

#include <cmath>

#include <gtest/gtest.h>

#include "adsim/equations/dispersion_model.h"

namespace adsim {
namespace {

double DeformationB(double eps) { return 1.25 * std::pow((1.0 - eps) / eps, 10.0 / 9.0); }

TEST(HeatDispersionStagnant, ZehnerSchlunderKappaOneGivesGasConductivity) {
  for (double eps : {0.3, 0.36, 0.4, 0.5, 0.7, 0.9})
    EXPECT_NEAR(ZehnerSchlunderConductivityRatio(eps, 1.0), 1.0, 1e-12) << "eps " << eps;
}

TEST(HeatDispersionStagnant, ZehnerSchlunderEmptyBedGivesGasConductivity) {
  EXPECT_DOUBLE_EQ(ZehnerSchlunderConductivityRatio(1.0, 50.0), 1.0);
  EXPECT_NEAR(ZehnerSchlunderConductivityRatio(0.999999, 50.0), 1.0, 1e-7);
}

TEST(HeatDispersionStagnant, ZehnerSchlunderRisesWithKappaAndFallsWithVoidage) {
  double previous = 0.0;
  for (double kappa : {1.0, 3.0, 10.0, 30.0, 100.0, 1000.0}) {
    const double k0 = ZehnerSchlunderConductivityRatio(0.4, kappa);
    EXPECT_GT(k0, previous) << "kappa " << kappa;
    previous = k0;
  }
  previous = 1e300;
  for (double eps : {0.3, 0.4, 0.5, 0.6, 0.8}) {
    const double k0 = ZehnerSchlunderConductivityRatio(eps, 100.0);
    EXPECT_LT(k0, previous) << "eps " << eps;
    previous = k0;
  }
}

TEST(HeatDispersionStagnant, ZehnerSchlunderIsFiniteAndContinuousAtKappaEqualsB) {
  const double eps = 0.4;
  const double B = DeformationB(eps);
  const double s = std::sqrt(1.0 - eps);
  // N = 0 exactly: the series limit 1 - s + 2 s ((B - 1)/3 + 1/2).
  const double limit = 1.0 - s + 2.0 * s * ((B - 1.0) / 3.0 + 0.5);
  EXPECT_NEAR(ZehnerSchlunderConductivityRatio(eps, B), limit, 1e-14);
  // Approaching from both sides converges to the limit.
  for (double d : {1e-4, 1e-8, 1e-12}) {
    EXPECT_NEAR(ZehnerSchlunderConductivityRatio(eps, B * (1.0 + d)), limit, 2.0 * d + 1e-13) << d;
    EXPECT_NEAR(ZehnerSchlunderConductivityRatio(eps, B * (1.0 - d)), limit, 2.0 * d + 1e-13) << d;
  }
  // Across the series / closed-form switch at |N| = 0.1 (kappa = B / 0.9 and
  // B / 1.1) the value does not jump.
  for (double n : {0.1, -0.1}) {
    const double kappa = B / (1.0 - n);
    const double below = ZehnerSchlunderConductivityRatio(eps, kappa * (1.0 - 1e-12));
    const double above = ZehnerSchlunderConductivityRatio(eps, kappa * (1.0 + 1e-12));
    EXPECT_NEAR(below, above, 1e-11) << "N " << n;
  }
}

TEST(HeatDispersionStagnant, ZehnerSchlunderValues) {
  EXPECT_NEAR(ZehnerSchlunderConductivityRatio(0.4, 10.0), 3.642726060371148, 1e-12);
  EXPECT_NEAR(ZehnerSchlunderConductivityRatio(0.4, 100.0), 8.886955605455523, 1e-12);
  EXPECT_NEAR(ZehnerSchlunderConductivityRatio(0.4, 1000.0), 15.467469424870831, 1e-11);
  EXPECT_NEAR(ZehnerSchlunderConductivityRatio(0.36, 50.0), 8.060452103938445, 1e-12);
}

TEST(HeatDispersionStagnant, KrupiczkaKappaOneGivesGasConductivity) {
  for (double eps : {0.3, 0.4, 0.5}) EXPECT_DOUBLE_EQ(KrupiczkaConductivityRatio(eps, 1.0), 1.0);
  EXPECT_GT(KrupiczkaConductivityRatio(0.4, 100.0), KrupiczkaConductivityRatio(0.4, 10.0));
}

}  // namespace
}  // namespace adsim
