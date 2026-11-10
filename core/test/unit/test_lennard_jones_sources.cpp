// Lennard-Jones parameters: provenance, and agreement with CoolProp (values
// recorded in coolprop_reference_values.h).
#include <cmath>
#include <iostream>
#include <string>
#include <vector>

#include <gtest/gtest.h>

#include "adsim/objects/fluid.h"
#include "coolprop_reference_values.h"

namespace adsim {

  // The engine's viscosity (polynomial) is within 2 % of CoolProp's
  // dilute-gas viscosity at 300 K and 500 K for every component CoolProp has a
  // viscosity model for. (Measured worst case: 0.8 %, pentane.)
  TEST(TransportProperties, LennardJonesMatchesReference) {
    const double tol = 0.02;
    double worst = 0.0; std::string worst_name;
    int checked = 0;
    for (const std::string& name : Components::componentNames) {
      const auto it = coolprop_reference::kViscosity.find(name);
      if (it == coolprop_reference::kViscosity.end()) continue;
      for (double T : {300.0, 500.0}) {
        const double ref = T == 300.0 ? it->second.mu300 : it->second.mu500;
        ASSERT_TRUE(std::isfinite(ref) && ref > 0) << name;
        const double mu = Fluid::GetComponentViscosity(name, T);
        const double err = std::abs(mu - ref) / ref;
        if (err > worst) { worst = err; worst_name = name; }
        EXPECT_LT(err, tol) << name << " at " << T << " K";
        ++checked;
      }
    }
    EXPECT_GT(checked, 20);
    std::cout << "worst viscosity deviation " << 100 * worst << " % (" << worst_name << ")\n";
  }

  // The CoolProp-sourced pairs equal what CoolProp reports (recorded from the
  // n-Octane and Toluene fluid files).
  TEST(TransportProperties, CoolPropSourcedPairsEqualCoolProp) {
    for (const std::string name : {"Octane", "Toluene"}) {
      const auto& ref = coolprop_reference::kViscosity.at(name);
      const auto& d = Components::GetComponentData(name);
      EXPECT_NEAR(d.eps, ref.eps, 0.01) << name;
      EXPECT_NEAR(d.sigma, ref.sigma, 0.001) << name;
    }
  }

  TEST(TransportProperties, LennardJonesProvenance) {
    EXPECT_EQ(Components::GetComponentDataSource("CarbonDioxide"), LJSource::Svehla1962);
    EXPECT_EQ(Components::GetComponentData("CarbonDioxide").eps, 195.2);
    EXPECT_EQ(Components::GetComponentData("CarbonDioxide").sigma, 3.941);
    EXPECT_EQ(Components::GetComponentDataSource("Isopentane"), LJSource::IsomerStandIn);
    EXPECT_EQ(Components::GetComponentData("Isopentane").eps, Components::GetComponentData("Pentane").eps);
    EXPECT_EQ(Components::GetComponentData("Isopentane").sigma, Components::GetComponentData("Pentane").sigma);
    EXPECT_NE(Components::LennardJonesCitation("Nitrogen").find("TR R-132"), std::string::npos);
    EXPECT_NE(Components::LennardJonesCitation("Octane").find("CoolProp"), std::string::npos);
  }

} // namespace adsim
