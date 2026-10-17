// ads includes
#include "adsim/objects/porous_media.h"
#include "adsim/objects/reactor.h"
#include "adsim/solver/cost_estimator.h"

// Gtest includes
#include <gtest/gtest.h>


namespace adsim {

  // ---------- CostEstimator::CalculateAdsorbentFillCost ----------
  // adsorbent_cost_per_kg is keyed by AdsorbentID, never by LayerID.

  TEST(CostEstimatorTest, AdsorbentFillCostTwoLayersSameAdsorbent) {
    AdsorbentDirectory dir;
    AdsorbentID id_a = dir.Add("Adsorbent A");
    dir.Get(id_a).solid_density_ = 1200.0;
    AdsorbentID id_b = dir.Add("Adsorbent B"); // priced, but in no layer
    dir.Get(id_b).solid_density_ = 900.0;

    Bed bed;
    std::get<CylinderBed>(bed.shape_).diameter = 0.1;
    const double area = (constants::kPI / 4.0) * 0.1 * 0.1;

    // Both layers are filled with adsorbent A. Layer 2's LayerID coincides
    // with adsorbent B's AdsorbentID, which is what a LayerID lookup hits.
    LayerID lid1 = bed.layers_.Add("Bottom");
    auto& layer1 = bed.layers_.Get(lid1);
    layer1.domain_.length = 0.5;
    layer1.eb_ = 0.4;
    layer1.AddAdsorbent(id_a, 1.0);

    LayerID lid2 = bed.layers_.Add("Top");
    auto& layer2 = bed.layers_.Get(lid2);
    layer2.domain_.length = 1.5;
    layer2.eb_ = 0.3;
    layer2.AddAdsorbent(id_a, 1.0);
    ASSERT_EQ(lid2, id_b);

    CostEstimator ce;
    ce.adsorbent_cost_per_kg[id_a] = 10.0;
    ce.adsorbent_cost_per_kg[id_b] = 1000.0;

    const double mass1 = area * 0.5 * (1.0 - 0.4) * 1200.0;
    const double mass2 = area * 1.5 * (1.0 - 0.3) * 1200.0;
    EXPECT_NEAR(ce.CalculateAdsorbentFillCost(bed, dir), (mass1 + mass2) * 10.0, 1e-9);
  }

  TEST(CostEstimatorTest, AdsorbentFillCostMixedLayerByMassFraction) {
    AdsorbentDirectory dir;
    AdsorbentID id_a = dir.Add("Light");
    dir.Get(id_a).solid_density_ = 600.0;
    AdsorbentID id_b = dir.Add("Heavy");
    dir.Get(id_b).solid_density_ = 1400.0;

    Bed bed;
    std::get<CylinderBed>(bed.shape_).diameter = 0.1;
    const double area = (constants::kPI / 4.0) * 0.1 * 0.1;
    LayerID lid = bed.layers_.Add("Mixed");
    auto& layer = bed.layers_.Get(lid);
    layer.domain_.length = 1.0;
    layer.eb_ = 0.37;
    layer.AddAdsorbent(id_a, 0.3);
    layer.AddAdsorbent(id_b, 0.7);

    CostEstimator ce;
    ce.adsorbent_cost_per_kg[id_a] = 5.0;
    ce.adsorbent_cost_per_kg[id_b] = 40.0;

    // Mixture particle density 1 / (0.3/600 + 0.7/1400) = 1000 kg/m3, so the
    // layer holds area * 1 * 0.63 * 1000 kg, split 30/70 by mass.
    const double layer_mass = area * 1.0 * (1.0 - 0.37) * 1000.0;
    ASSERT_NEAR(GetAdsorbentMass(bed, dir), layer_mass, 1e-10);
    const double expected = layer_mass * (0.3 * 5.0 + 0.7 * 40.0);
    EXPECT_NEAR(ce.CalculateAdsorbentFillCost(bed, dir), expected, 1e-9);

    // An adsorbent with no price entry contributes nothing.
    ce.adsorbent_cost_per_kg.erase(id_b);
    EXPECT_NEAR(ce.CalculateAdsorbentFillCost(bed, dir), layer_mass * 0.3 * 5.0, 1e-9);
  }

  // ---------- CostEstimator::CalculateCost (LCOC) ----------
  // L-04: with no production LCOC is undefined -> 0 and defined=false (kpis.h
  // convention), so the optimiser does not read it as a free £0/kg.
  TEST(CostEstimatorTest, LcocWithNoProductionIsUndefined) {
    AdsorbentDirectory dir;
    Bed bed;
    EquationLog eq_log;
    BoundaryLog b_log;
    eq_log.sys_data.adsorbent_mass = 1.0;
    eq_log.sys_data.total_cycle_time = 1.0;
    KPIParams params(eq_log, b_log);  // no endpoints: nothing is produced

    CostEstimator ce;
    bool defined = true;
    const double lcoc = ce.CalculateCost(bed, dir, params, &defined);
    EXPECT_EQ(lcoc, 0.0);
    EXPECT_FALSE(defined);
  }

} // namespace adsim
