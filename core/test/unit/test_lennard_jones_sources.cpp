// Lennard-Jones parameters: provenance, and agreement with CoolProp.
#include <cmath>
#include <string>
#include <vector>

#include <gtest/gtest.h>

#include "CoolProp.h"
#include "adsim/objects/fluid.h"

namespace adsim {

  namespace {
    // Number after `"key":` in a JSON text (cereal and CoolProp bundle clashing
    // rapidjson copies, so the two values are read from the text).
    bool JsonNumber(const std::string& js, const std::string& key, double& out) {
      const auto k = js.find("\"" + key + "\"");
      if (k == std::string::npos) return false;
      const auto c = js.find(':', k);
      if (c == std::string::npos) return false;
      try { out = std::stod(js.substr(c + 1)); } catch (...) { return false; }
      return true;
    }

    // CoolProp's own LJ pair for the viscosity model, or false if it has none.
    bool CoolPropPair(const std::string& fluid, double& eps, double& sigma) {
      try {
        const std::string js = CoolProp::get_fluid_param_string(fluid, "JSON");
        double s_m;
        if (!JsonNumber(js, "sigma_eta", s_m) || !JsonNumber(js, "epsilon_over_k", eps)) return false;
        sigma = s_m * 1e10;  // m -> Angstrom
        return true;
      } catch (...) { return false; }
    }
  }

  // The engine's viscosity (polynomial) is within 2 % of CoolProp's
  // dilute-gas viscosity at 300 K and 500 K for every component CoolProp has a
  // viscosity model for. (Measured worst case: 0.8 %, pentane.)
  TEST(TransportProperties, LennardJonesMatchesReference) {
    const double tol = 0.02;
    double worst = 0.0; std::string worst_name;
    int checked = 0;
    for (const std::string& name : Components::componentNames) {
      double e, s;
      if (!CoolPropPair(name, e, s)) continue;
      for (double T : {300.0, 500.0}) {
        double ref;
        try { ref = CoolProp::PropsSI("V", "T", T, "P", 1.0e3, name); } catch (...) { continue; }
        if (!std::isfinite(ref) || ref <= 0) continue;
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

  // The CoolProp-sourced pairs equal what CoolProp reports, so a CoolProp
  // upgrade that changes them is noticed.
  TEST(TransportProperties, CoolPropSourcedPairsEqualCoolProp) {
    const std::vector<std::pair<std::string, std::string>> cp = {
      {"Octane", "n-Octane"}, {"Toluene", "Toluene"}};
    for (const auto& [name, fluid] : cp) {
      double e, s;
      ASSERT_TRUE(CoolPropPair(fluid, e, s)) << fluid;
      const auto& d = Components::GetComponentData(name);
      EXPECT_NEAR(d.eps, e, 0.01) << name;
      EXPECT_NEAR(d.sigma, s, 0.001) << name;
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
