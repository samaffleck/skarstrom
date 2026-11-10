// STL includes
#include <iostream>
#include <cmath>
#include <algorithm>

#include <cereal/archives/json.hpp>

// ads includes
#include "adsim/objects/fluid.h"
#include "adsim/equations/binary_diffusion_model.h"

// Gtest includes
#include <gtest/gtest.h>

#include "adsim/equations/isotherm.h"
#include "adsim/objects/component_properties.h"
#include "adsim/objects/manifests.h"
#include "adsim/simulation.h"

#include <fstream>

#include <set>
#include <sstream>
#include <stdexcept>


namespace adsim {

  // Every engine component has a row in the generated table, and every row is
  // an engine component.
  TEST(FluidTests, EveryComponentHasPropertiesInTheTable) {
    ASSERT_EQ(kComponentProperties.size(), Components::componentNames.size());
    std::set<std::string_view> rows;
    for (const auto& c : kComponentProperties) rows.insert(c.name);
    for (const std::string& component : Components::componentNames) {
      ASSERT_TRUE(rows.count(component)) << component;
      const double mr = GetMolarMass(component);
      EXPECT_TRUE(std::isfinite(mr) && mr > 1e-3 && mr < 1.0) << component;
      const ComponentProperties* p = FindComponentProperties(component);
      EXPECT_GT(p->T_critical, p->T_triple) << component;
      EXPECT_GT(p->p_critical, 0.0) << component;
      EXPECT_FALSE(p->cas.empty()) << component;
    }
    EXPECT_EQ(GetMolarMass("NotAComponent"), 0.0);
    EXPECT_EQ(FindComponentProperties("NotAComponent"), nullptr);
  }

  // Spot values of the table: the molar masses of the CoolProp 6.7.0 fluid
  // files (Nitrogen.json, CarbonDioxide.json, Water.json: EOS molar_mass), and
  // water's saturation pressure from its vapour-pressure ancillary (Water.json
  // ANCILLARIES.pS): the values CoolProp's saturation_ancillary(iP, 1, iT, T)
  // gave the engine before the table (recorded from the 0.9.1.dev0 engine).
  TEST(FluidTests, TableSpotValues) {
    EXPECT_DOUBLE_EQ(GetMolarMass("Nitrogen"), 0.02801348);
    EXPECT_DOUBLE_EQ(GetMolarMass("CarbonDioxide"), 0.0440098);
    EXPECT_DOUBLE_EQ(GetMolarMass("Water"), 0.018015268);
    EXPECT_NEAR(SaturationPressure("Water", 298.15), 3169.8911918548201, 1e-9);
    EXPECT_NEAR(SaturationPressure("Water", 313.15), 7385.1909268284771, 1e-9);
    const ComponentProperties* w = FindComponentProperties("Water");
    EXPECT_EQ(w->formula, "H2O");
    EXPECT_EQ(w->cas, "7732-18-5");
    EXPECT_DOUBLE_EQ(w->T_critical, 647.096);
    EXPECT_DOUBLE_EQ(w->p_critical, 22064000.0);
    EXPECT_DOUBLE_EQ(w->T_triple, 273.16);
  }

  // The component manifest reads the table.
  TEST(FluidTests, ComponentManifestReadsTheTable) {
    const ComponentInfo co2 = DescribeComponent("CarbonDioxide");
    EXPECT_EQ(co2.formula, "CO2");
    EXPECT_EQ(co2.cas, "124-38-9");
    EXPECT_DOUBLE_EQ(co2.molar_mass, 0.0440098);
    EXPECT_NE(std::find(co2.aliases.begin(), co2.aliases.end(), "R744"), co2.aliases.end());
    EXPECT_TRUE(co2.ideal_gas_heat_capacity);
    const std::string manifest = ComponentManifest();
    EXPECT_NE(manifest.find(kComponentPropertiesSource), std::string::npos);
    EXPECT_EQ(manifest.find("\"heat_capacity_held\":true"), std::string::npos);
  }

  // Hexane was removed: it is not a component, and a deck that has it fails
  // to load with an error that names it.
  TEST(FluidTests, HexaneIsRemoved) {
    EXPECT_EQ(std::find(Components::componentNames.begin(), Components::componentNames.end(), "Hexane"),
              Components::componentNames.end());
    EXPECT_EQ(FindComponentProperties("Hexane"), nullptr);
    std::istringstream in(R"({"value0": {"cereal_class_version": 4, "Components": ["Nitrogen", "Hexane"]}})");
    Fluid f;
    try {
      cereal::JSONInputArchive ar(in);
      ar(f);
      FAIL() << "a fluid with Hexane loaded";
    } catch (const std::invalid_argument& e) {
      const std::string what = e.what();
      EXPECT_NE(what.find("Hexane"), std::string::npos) << what;
      EXPECT_NE(what.find("removed"), std::string::npos) << what;
    }
  }

  // The same through a whole deck: a shipped deck with Hexane added to its
  // component list does not load, and the error names Hexane.
  TEST(FluidTests, DeckWithHexaneFailsToLoad) {
    std::ifstream in("../examples/decks/Arora2021_LiLSX_PSA.json");
    ASSERT_TRUE(in) << "run from core/";
    std::stringstream ss;
    ss << in.rdbuf();
    std::string deck = ss.str();
    const std::string key = "\"Components\": [";
    const auto at = deck.find(key);
    ASSERT_NE(at, std::string::npos);
    deck.insert(at + key.size(), "\"Hexane\", ");
    try {
      Simulation::FromJson(deck);
      FAIL() << "a deck with Hexane loaded";
    } catch (const std::exception& e) {
      EXPECT_NE(std::string(e.what()).find("Hexane"), std::string::npos) << e.what();
    }
  }

  TEST(FluidTests, BinaryDiffusionCheck) {
    std::string CO2 = "CarbonDioxide";
    std::string N2 = "Nitrogen";
    BinaryDiffusionCoefficientModel dif_model;
    dif_model.Initialise(1, { CO2, N2 });
    dif_model.Update({ 590.0 }, { 1.0e5 });
    auto coeffs = dif_model.m_DiffusionCoefficients.at({CO2, N2});
    ASSERT_NEAR(coeffs.Dm_ij[0], 5.16e-5, 1e-5);
  }

} // end adsim namespace
