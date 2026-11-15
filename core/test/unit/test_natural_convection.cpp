/*
 * Skarstrom - pressure and temperature swing adsorption simulator
 * Copyright (c) 2024-2026 Sam Affleck and contributors
 *
 * Licensed under the terms in the LICENSE file at the repository root.
 */
// Wall-to-ambient natural convection (adsim/equations/ambient_air.h):
//   Churchill and Chu (1975, Int. J. Heat Mass Transfer 18, 1323), vertical plate:
//     Nu_L = {0.825 + 0.387 Ra^(1/6) / [1 + (0.492/Pr)^(9/16)]^(8/27)}^2;
//   Churchill and Chu (1975, Int. J. Heat Mass Transfer 18, 1049), horizontal cylinder:
//     Nu_D = {0.60 + 0.387 Ra^(1/6) / [1 + (0.559/Pr)^(9/16)]^(8/27)}^2;
// and the ambient air properties (Sutherland's law, White 2006). Checked here:
// the analytic limits, values worked out by hand (Python, double precision)
// from the formulas above, the air properties against Incropera et al.'s
// Table A.4, the turbulent-plate property that h does not depend on L, the
// kinds' storage (factory, deck archive, builder, registry) and the
// slender-column deck warning. The solver-level checks are in
// validation/test_natural_convection_solver.cpp.

#include <cmath>
#include <sstream>
#include <string>

#include <cereal/archives/json.hpp>
#include <gtest/gtest.h>

#include "adsim/equations/ambient_air.h"
#include "adsim/equations/heat_transfer_model.h"
#include "adsim/objects/deck_builder.h"
#include "adsim/objects/model_kinds.h"
#include "adsim/objects/reactor.h"
#include "adsim/solver/deck_warnings.h"

namespace adsim {
namespace {

namespace air = ambient_air;

// Ra -> 0: the conduction limits of the two correlations, 0.825^2 and 0.60^2.
TEST(NaturalConvection, StagnantLimits) {
  for (double pr : {0.7, 1.0, 5.0}) {
    EXPECT_DOUBLE_EQ(air::ChurchillChuVerticalPlateNu(0.0, pr), 0.825 * 0.825);
    EXPECT_DOUBLE_EQ(air::ChurchillChuHorizontalCylinderNu(0.0, pr), 0.36);
  }
}

// Values worked out by hand from the formulas in the header comment.
TEST(NaturalConvection, HandValues) {
  EXPECT_NEAR(air::ChurchillChuVerticalPlateNu(1e4, 0.71), 5.4327454632929415, 1e-12);
  EXPECT_NEAR(air::ChurchillChuVerticalPlateNu(1e8, 0.71), 61.06517223358536, 1e-10);
  EXPECT_NEAR(air::ChurchillChuHorizontalCylinderNu(1e4, 0.71), 4.373272099562677, 1e-12);
  EXPECT_NEAR(air::ChurchillChuHorizontalCylinderNu(1e8, 0.71), 56.57610488040866, 1e-10);
  // A wall at 350 K in air at 300 K: a 1 m vertical column (Ra = 3.24e9) and
  // a 0.1 m horizontal one (Ra = 3.24e6), film properties at 325 K.
  EXPECT_NEAR(air::Rayleigh(350.0, 300.0, 1.0), 3241024140.803931, 1e-3);
  EXPECT_NEAR(air::ChurchillChuCoefficient(true, 350.0, 300.0, 1.0, 0.1), 4.978324913015165, 1e-12);
  EXPECT_NEAR(air::ChurchillChuCoefficient(false, 350.0, 300.0, 1.0, 0.1), 5.716152079683699, 1e-12);
}

// The air properties against Incropera et al., Fundamentals of Heat and Mass
// Transfer, Table A.4 (air at 1 atm): within 2 % from 300 to 400 K.
TEST(NaturalConvection, AirPropertiesMatchTables) {
  struct Row { double T, mu, k, rho, pr; };
  const Row rows[] = {{300.0, 184.6e-7, 26.3e-3, 1.1614, 0.707},
                      {350.0, 208.2e-7, 30.0e-3, 0.9950, 0.700},
                      {400.0, 230.1e-7, 33.8e-3, 0.8711, 0.690}};
  for (const auto& r : rows) {
    EXPECT_NEAR(air::Viscosity(r.T) / r.mu, 1.0, 0.02) << r.T;
    EXPECT_NEAR(air::ThermalConductivity(r.T) / r.k, 1.0, 0.02) << r.T;
    EXPECT_NEAR(air::Density(r.T) / r.rho, 1.0, 0.02) << r.T;
    EXPECT_NEAR(air::Prandtl(r.T) / r.pr, 1.0, 0.02) << r.T;
  }
}

// Turbulent plate: Nu -> Ra^(1/3) at large Ra, so h does not depend on L.
TEST(NaturalConvection, TurbulentPlateCoefficientIndependentOfHeight) {
  const double h10 = air::ChurchillChuCoefficient(true, 343.0, 293.0, 10.0, 1.0);
  const double h20 = air::ChurchillChuCoefficient(true, 343.0, 293.0, 20.0, 1.0);
  EXPECT_NEAR(h10 / h20, 1.0119118262427629, 1e-12);
  // Laminar plate (Nu ~ Ra^(1/4)): h falls with the height.
  EXPECT_GT(air::ChurchillChuCoefficient(true, 303.0, 293.0, 0.05, 0.01),
            1.5 * air::ChurchillChuCoefficient(true, 303.0, 293.0, 0.5, 0.01));
}

// Monotone in the temperature difference; finite (the conduction limit) at 0.
TEST(NaturalConvection, MonotoneInTemperatureDifference) {
  double last = 0.0;
  for (double dT : {0.0, 1.0, 10.0, 50.0, 150.0}) {
    const double h = air::ChurchillChuCoefficient(true, 300.0 + dT, 300.0, 1.0, 0.1);
    EXPECT_GT(h, last) << dT;
    last = h;
  }
  EXPECT_DOUBLE_EQ(air::ChurchillChuCoefficient(true, 300.0, 300.0, 2.0, 0.1),
                   0.825 * 0.825 * air::ThermalConductivity(300.0) / 2.0);
  // A wall colder than the air loses the same way (|dT|).
  EXPECT_GT(air::ChurchillChuCoefficient(false, 280.0, 300.0, 1.0, 0.1), 0.36 * air::ThermalConductivity(290.0) / 0.1);
}

TEST(NaturalConvection, KindsAndFactory) {
  const auto& tags = build::ExternalHeatTransferKinds();
  ASSERT_EQ(tags.size(), 4u);
  EXPECT_EQ(tags[2], "churchill_chu_vertical");
  EXPECT_EQ(tags[3], "churchill_chu_horizontal");
  auto v = CreateExternalHeatTransferModel(ExternalHeatTransferType::CHURCHILL_CHU_VERTICAL);
  auto h = CreateExternalHeatTransferModel(ExternalHeatTransferType::CHURCHILL_CHU_HORIZONTAL);
  EXPECT_EQ(v->GetModelType(), ExternalHeatTransferType::CHURCHILL_CHU_VERTICAL);
  EXPECT_EQ(h->GetModelType(), ExternalHeatTransferType::CHURCHILL_CHU_HORIZONTAL);
  EXPECT_TRUE(ModelFields(v.get()).empty());
}

// The orientation survives a deck save and load.
TEST(NaturalConvection, OrientationRoundTripsThroughTheArchive) {
  for (bool vertical : {true, false}) {
    std::unique_ptr<IExternalHeatTransferModel> in = std::make_unique<ExternalHeatTransferModel_ChurchillChu>(vertical);
    std::stringstream ss;
    {
      cereal::JSONOutputArchive ar(ss);
      ar(cereal::make_nvp("model", in));
    }
    std::unique_ptr<IExternalHeatTransferModel> out;
    {
      cereal::JSONInputArchive ar(ss);
      ar(cereal::make_nvp("model", out));
    }
    ASSERT_TRUE(out);
    EXPECT_EQ(out->GetModelType(), in->GetModelType());
  }
}

// A non-isothermal bed of one layer with a 2 mm wall.
Bed SmallBed(ExternalHeatTransferType kind, double length, double diameter) {
  Bed bed;
  bed.shape_ = CylinderBed{ExpressionParam{diameter}};
  bed.wall_.thickness_ = ExpressionParam{0.002};
  bed.ambient_temperature_ = ExpressionParam{293.0};
  const auto id = bed.layers_.Add("Main");
  bed.layers_.Get(id).domain_.length = ExpressionParam{length};
  bed.layers_.Get(id).equations_.external_heat_transfer_model_ = CreateExternalHeatTransferModel(kind);
  return bed;
}

TEST(NaturalConvection, SlenderColumnWarning) {
  using T = ExternalHeatTransferType;
  // A 1 m lab column 30 mm across: a plate needs about 0.18 m.
  {
    const auto w = NaturalConvectionWarnings(SmallBed(T::CHURCHILL_CHU_VERTICAL, 1.0, 0.03));
    ASSERT_EQ(w.size(), 1u);
    EXPECT_NE(w[0].find("too slender"), std::string::npos) << w[0];
  }
  // A 5 m column 1.5 m across is a plate.
  EXPECT_TRUE(NaturalConvectionWarnings(SmallBed(T::CHURCHILL_CHU_VERTICAL, 5.0, 1.5)).empty());
  // The horizontal kind and the other kinds are not checked, nor an isothermal run.
  EXPECT_TRUE(NaturalConvectionWarnings(SmallBed(T::CHURCHILL_CHU_HORIZONTAL, 1.0, 0.03)).empty());
  EXPECT_TRUE(NaturalConvectionWarnings(SmallBed(T::CONSTANT, 1.0, 0.03)).empty());
  Bed iso = SmallBed(T::CHURCHILL_CHU_VERTICAL, 1.0, 0.03);
  iso.equations_.energy_balance_mode_ = EnergyBalanceMode::ISOTHERMAL;
  EXPECT_TRUE(NaturalConvectionWarnings(iso).empty());
}

}  // namespace
}  // namespace adsim
