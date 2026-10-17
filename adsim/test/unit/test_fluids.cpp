// STL includes
#include <iostream>
#include <cmath>

// ads includes
#include "adsim/objects/fluid.h"
#include "adsim/equations/binary_diffusion_model.h"

// Gtest includes
#include <gtest/gtest.h>

// CoolProp includes
#include "CoolProp.h"


namespace adsim {

  TEST(FluidTests, MolarMass) {
    for (const std::string& component : Components::componentNames) {
      double mr = CoolProp::Props1SI(component, "molar_mass");
      assert(std::isfinite(mr));
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
