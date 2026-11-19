/*
 * Skarstrom - pressure and temperature swing adsorption simulator
 * Copyright (c) 2024-2026 Sam Affleck and contributors
 *
 * Licensed under the terms in the LICENSE file at the repository root.
 */
// Wall-to-ambient natural and forced convection (adsim/equations/ambient_air.h):
//   Churchill and Chu (1975, Int. J. Heat Mass Transfer 18, 1323), vertical plate:
//     Nu_L = {0.825 + 0.387 Ra^(1/6) / [1 + (0.492/Pr)^(9/16)]^(8/27)}^2;
//   Churchill and Chu (1975, Int. J. Heat Mass Transfer 18, 1049), horizontal cylinder:
//     Nu_D = {0.60 + 0.387 Ra^(1/6) / [1 + (0.559/Pr)^(9/16)]^(8/27)}^2;
//   Yang (1985), vertical cylinder at any H/d:
//     Nu_H = {(0.36 H/d)^(1/2) + [0.150 (Ra_H / [1 + (0.492/Pr)^(9/16)]^(16/9))^(1/3)]^(1/2)}^2;
//   Churchill and Bernstein (1977, J. Heat Transfer 99, 300), cylinder in crossflow:
//     Nu_D = 0.3 + 0.62 Re^(1/2) Pr^(1/3) / [1 + (0.4/Pr)^(2/3)]^(1/4) [1 + (Re/282000)^(5/8)]^(4/5);
//   mixed convection (Churchill 1977 blend), h = (h_F^n + h_N^n)^(1/n), with h_N Yang
//     (vertical) or Churchill-Chu horizontal and h_F Churchill-Bernstein less its 0.3 k/d;
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
#include "adsim/objects/deck_edit.h"
#include "adsim/objects/model_kinds.h"
#include "adsim/objects/reactor.h"
#include "adsim/objects/user_data.h"
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
  // Laminar plate (Ra 1.2e5 to 1.2e8, Nu nearer Ra^(1/4)): h falls with the
  // height. Hand value of h(0.05 m) / h(0.5 m).
  EXPECT_NEAR(air::ChurchillChuCoefficient(true, 303.0, 293.0, 0.05, 0.01) /
                  air::ChurchillChuCoefficient(true, 303.0, 293.0, 0.5, 0.01),
              1.4928534406221199, 1e-12);
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

// Yang (1985), vertical cylinder at any H/d: Ra -> 0 gives Nu_H = 0.36 H/d,
// so h -> 0.36 k/d (the conduction limit of a horizontal cylinder).
TEST(NaturalConvection, YangStagnantLimit) {
  for (double hd : {0.5, 1.0, 33.0}) EXPECT_DOUBLE_EQ(air::YangVerticalCylinderNu(0.0, 0.71, hd), 0.36 * hd);
}

// Hand values (Python, double precision) of Yang's formula in the header.
TEST(NaturalConvection, YangHandValues) {
  EXPECT_NEAR(air::YangVerticalCylinderNu(1e4, 0.71, 1.0), 4.439389527427472, 1e-12);
  EXPECT_NEAR(air::YangVerticalCylinderNu(1e8, 0.71, 10.0), 79.07078442854775, 1e-10);
  EXPECT_NEAR(air::YangVerticalCylinderNu(1e9, 0.71, 33.3), 188.49489930411173, 1e-10);
  // A 1 m lab column 30 mm across, wall 350 K in air at 300 K: 44 % above
  // the vertical-plate (Churchill-Chu) value of 4.978 W/(m2 K).
  EXPECT_NEAR(air::YangCoefficient(350.0, 300.0, 1.0, 0.03), 7.159659645323737, 1e-12);
  EXPECT_NEAR(air::YangCoefficient(350.0, 300.0, 1.0, 1.0), 4.818319943587659, 1e-12);
}

// With H/d = 0.825^2 / 0.36 the two conduction terms are equal, so Yang and
// the Churchill-Chu plate differ only by 0.150 against 0.387^2 in the
// convection term: within 0.16 % at every Ra.
TEST(NaturalConvection, YangMatchesThePlateWhenTheConductionTermsAgree) {
  const double hd = 0.825 * 0.825 / 0.36;
  for (double ra : {1e2, 1e5, 1e8, 1e11}) {
    const double ratio = air::YangVerticalCylinderNu(ra, 0.71, hd) / air::ChurchillChuVerticalPlateNu(ra, 0.71);
    EXPECT_NEAR(ratio, 1.0, 0.0016) << ra;
  }
}

// A slender column loses more heat than a plate of its height; the gap
// closes as the column gets wider.
TEST(NaturalConvection, YangSlenderColumnAbovePlate) {
  const double plate = air::ChurchillChuCoefficient(true, 330.0, 300.0, 1.0, 0.1);
  double last = 1e9;
  for (double d : {0.01, 0.03, 0.1, 0.3}) {
    const double h = air::YangCoefficient(330.0, 300.0, 1.0, d);
    EXPECT_LT(h, last) << d;
    last = h;
  }
  EXPECT_GT(air::YangCoefficient(330.0, 300.0, 1.0, 0.03), 1.2 * plate);
}

// Churchill and Bernstein (1977): Re -> 0 leaves the constant 0.3 at any Pr.
TEST(ForcedConvection, ChurchillBernsteinStagnantLimit) {
  for (double pr : {0.5, 0.71, 7.0}) EXPECT_DOUBLE_EQ(air::ChurchillBernsteinNu(0.0, pr), 0.3);
  // A negative Re (a sign slip) is clamped to the stagnant value, not NaN.
  EXPECT_DOUBLE_EQ(air::ChurchillBernsteinNu(-5.0, 0.71), 0.3);
  // With no air flow, h = 0.3 k / d (film temperature 300 K, d = 30 mm).
  EXPECT_DOUBLE_EQ(air::ChurchillBernsteinCoefficient(300.0, 300.0, 0.0, 0.03),
                   0.3 * air::ThermalConductivity(300.0) / 0.03);
}

// Hand values (Python, double precision) of the formula in the header, at
// Pr = 0.71. A textbook check: Re = 6071, Pr = 0.7 gives Nu = 40.6, the value
// of Incropera et al.'s worked example of a cylinder in crossflow (to be
// confirmed against the book).
TEST(ForcedConvection, ChurchillBernsteinHandValues) {
  EXPECT_NEAR(air::ChurchillBernsteinNu(1e3, 0.71), 16.018791873942707, 1e-11);
  EXPECT_NEAR(air::ChurchillBernsteinNu(1e5, 0.71), 215.34609302481786, 1e-10);
  EXPECT_NEAR(air::ChurchillBernsteinNu(1e6, 0.71), 1233.7195753575413, 1e-9);
  EXPECT_NEAR(air::ChurchillBernsteinNu(6071.0, 0.7), 40.6, 0.05);
  // A lab column 30 mm across, wall 310 K in air at 300 K crossing at
  // 0.1 m/s (Re_D = 186): h = 6.18 W/(m2 K). A 0.1 m column at 350 K in
  // air at 300 K and 2 m/s (Re_D = 11 065): h = 15.9 W/(m2 K).
  EXPECT_NEAR(air::ChurchillBernsteinCoefficient(310.0, 300.0, 0.1, 0.03), 6.179707809873068, 1e-11);
  EXPECT_NEAR(air::ChurchillBernsteinCoefficient(350.0, 300.0, 2.0, 0.1), 15.90337881484203, 1e-11);
}

// h rises with the air velocity, and the wind direction does not matter.
TEST(ForcedConvection, ChurchillBernsteinMonotoneInVelocity) {
  double last = air::ChurchillBernsteinCoefficient(320.0, 300.0, 0.0, 0.05);
  for (double v : {0.01, 0.1, 0.5, 2.0, 10.0, 50.0}) {
    const double h = air::ChurchillBernsteinCoefficient(320.0, 300.0, v, 0.05);
    EXPECT_GT(h, last) << v;
    last = h;
  }
  EXPECT_DOUBLE_EQ(air::ChurchillBernsteinCoefficient(320.0, 300.0, -1.0, 0.05),
                   air::ChurchillBernsteinCoefficient(320.0, 300.0, 1.0, 0.05));
}

namespace mixed {
double ForcedTerm(double Tw, double Tinf, double v, double d) {
  return air::ChurchillBernsteinCoefficient(Tw, Tinf, v, d) - 0.3 * air::ThermalConductivity(0.5 * (Tw + Tinf)) / d;
}
}  // namespace mixed

// v = 0 gives exactly the natural-convection kind of the orientation: Yang's
// vertical cylinder, or Churchill and Chu's horizontal cylinder.
TEST(MixedConvection, NoWindIsTheNaturalKind) {
  for (double n : {1.0, 3.0, 6.0}) {
    EXPECT_DOUBLE_EQ(air::MixedConvectionCoefficient(true, 310.0, 300.0, 0.0, 1.0, 0.03, n),
                     air::YangCoefficient(310.0, 300.0, 1.0, 0.03));
    EXPECT_DOUBLE_EQ(air::MixedConvectionCoefficient(false, 310.0, 300.0, 0.0, 1.0, 0.03, n),
                     air::ChurchillChuCoefficient(false, 310.0, 300.0, 1.0, 0.03));
  }
}

// Strong wind (Re_D >= 1e4): within 1 % of churchill_bernstein. The forced
// term lacks the 0.3 of Nu >= 53.6, and the natural term adds little.
TEST(MixedConvection, StrongWindIsChurchillBernstein) {
  const double nu_air = air::Viscosity(305.0) / air::Density(305.0);
  for (double re : {1e4, 1e5}) {
    for (bool vertical : {true, false}) {
      const double v = re * nu_air / 0.03;
      const double h = air::MixedConvectionCoefficient(vertical, 310.0, 300.0, v, 1.0, 0.03, 3.0);
      EXPECT_NEAR(h / air::ChurchillBernsteinCoefficient(310.0, 300.0, v, 0.03), 1.0, 0.01) << re << vertical;
    }
  }
}

// Hand values (Python, double precision) of the blend. A lab column 30 mm
// across and 1 m tall, wall 310 K in air at 300 K crossing at 0.1 m/s:
// h_F = 5.913 (6.180 by Churchill-Bernstein less 0.3 k/d), h_N = 4.849
// (Yang), h = 6.846 W/(m2 K). (With the full Churchill-Bernstein value the
// blend would be 7.05; with the vertical-plate h_N = 3.13 and the full
// value, 6.44.)
TEST(MixedConvection, HandValues) {
  EXPECT_NEAR(mixed::ForcedTerm(310.0, 300.0, 0.1, 0.03), 5.913366110760218, 1e-11);
  EXPECT_NEAR(air::MixedConvectionCoefficient(true, 310.0, 300.0, 0.1, 1.0, 0.03, 3.0), 6.845732030430517, 1e-11);
  EXPECT_NEAR(air::MixedConvectionCoefficient(false, 310.0, 300.0, 0.1, 1.0, 0.03, 3.0), 6.8084109905233285, 1e-11);
  EXPECT_NEAR(air::MixedConvectionCoefficient(true, 350.0, 300.0, 1.0, 2.0, 0.1, 4.0), 11.036954398488174, 1e-11);
  // n = 1 is the plain sum.
  EXPECT_NEAR(air::MixedConvectionCoefficient(true, 310.0, 300.0, 0.1, 1.0, 0.03, 1.0),
              mixed::ForcedTerm(310.0, 300.0, 0.1, 0.03) + air::YangCoefficient(310.0, 300.0, 1.0, 0.03), 1e-12);
}

// h rises with the wind and with the wall-to-air difference, it is
// continuous at v = 0, and it lies between max(h_F, h_N) and h_F + h_N.
TEST(MixedConvection, MonotoneAndContinuous) {
  for (bool vertical : {true, false}) {
    double last = air::MixedConvectionCoefficient(vertical, 320.0, 300.0, 0.0, 1.0, 0.03, 3.0);
    for (double v : {1e-6, 1e-3, 0.01, 0.05, 0.1, 0.5, 2.0, 10.0}) {
      const double h = air::MixedConvectionCoefficient(vertical, 320.0, 300.0, v, 1.0, 0.03, 3.0);
      EXPECT_GT(h, last) << v;
      if (v == 1e-6) EXPECT_NEAR(h, last, 1e-3 * last);
      const double hF = mixed::ForcedTerm(320.0, 300.0, v, 0.03);
      const double hN = vertical ? air::YangCoefficient(320.0, 300.0, 1.0, 0.03)
                                 : air::ChurchillChuCoefficient(false, 320.0, 300.0, 1.0, 0.03);
      EXPECT_GE(h, std::max(hF, hN));
      EXPECT_LE(h, hF + hN);
      last = h;
    }
    last = 0.0;
    for (double dT : {1.0, 5.0, 10.0, 30.0, 80.0}) {
      const double h = air::MixedConvectionCoefficient(vertical, 300.0 + dT, 300.0, 0.1, 1.0, 0.03, 3.0);
      EXPECT_GT(h, last) << dT;
      last = h;
    }
  }
}

TEST(NaturalConvection, KindsAndFactory) {
  const auto& tags = build::ExternalHeatTransferKinds();
  ASSERT_EQ(tags.size(), 8u);
  EXPECT_EQ(tags[6], "mixed_convection_vertical");
  EXPECT_EQ(tags[7], "mixed_convection_horizontal");
  for (auto t : {ExternalHeatTransferType::MIXED_CONVECTION_VERTICAL, ExternalHeatTransferType::MIXED_CONVECTION_HORIZONTAL}) {
    auto m = CreateExternalHeatTransferModel(t);
    ASSERT_TRUE(m);
    EXPECT_EQ(m->GetModelType(), t);
    const auto f = ModelFields(m.get());
    ASSERT_EQ(f.size(), 2u);
    EXPECT_EQ(f[0].first, "external_velocity");
    EXPECT_EQ(f[1].first, "mixed_convection_exponent");
    EXPECT_DOUBLE_EQ(f[1].second->Get(), 3.0);
  }
  EXPECT_EQ(tags[1], "holman");
  EXPECT_EQ(tags[2], "churchill_chu_vertical");
  EXPECT_EQ(tags[3], "churchill_chu_horizontal");
  EXPECT_EQ(tags[4], "yang_vertical_cylinder");
  EXPECT_EQ(tags[5], "churchill_bernstein");
  auto cb = CreateExternalHeatTransferModel(ExternalHeatTransferType::CHURCHILL_BERNSTEIN);
  ASSERT_TRUE(cb);
  EXPECT_EQ(cb->GetModelType(), ExternalHeatTransferType::CHURCHILL_BERNSTEIN);
  const auto fields = ModelFields(cb.get());
  ASSERT_EQ(fields.size(), 1u);
  EXPECT_EQ(fields[0].first, "external_velocity");
  EXPECT_EQ(CreateExternalHeatTransferModel(ExternalHeatTransferType::YANG_VERTICAL_CYLINDER)->GetModelType(),
            ExternalHeatTransferType::YANG_VERTICAL_CYLINDER);
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

// The external velocity survives a deck save and load.
TEST(ForcedConvection, ChurchillBernsteinRoundTripsThroughTheArchive) {
  auto model = std::make_unique<ExternalHeatTransferModel_ChurchillBernstein>();
  model->external_velocity = ExpressionParam{1.7};
  std::unique_ptr<IExternalHeatTransferModel> in = std::move(model);
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
  auto* cb = dynamic_cast<ExternalHeatTransferModel_ChurchillBernstein*>(out.get());
  ASSERT_TRUE(cb);
  EXPECT_DOUBLE_EQ(cb->external_velocity.Get(), 1.7);
}

// The orientation, velocity and exponent survive a deck save and load.
TEST(MixedConvection, RoundTripsThroughTheArchive) {
  for (bool vertical : {true, false}) {
    auto model = std::make_unique<ExternalHeatTransferModel_MixedConvection>(vertical);
    model->external_velocity = ExpressionParam{0.2};
    model->exponent = ExpressionParam{4.0};
    std::unique_ptr<IExternalHeatTransferModel> in = std::move(model);
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
    auto* m = dynamic_cast<ExternalHeatTransferModel_MixedConvection*>(out.get());
    ASSERT_TRUE(m);
    EXPECT_EQ(m->Vertical(), vertical);
    EXPECT_DOUBLE_EQ(m->external_velocity.Get(), 0.2);
    EXPECT_DOUBLE_EQ(m->exponent.Get(), 4.0);
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
    EXPECT_NE(w[0].find("yang_vertical_cylinder"), std::string::npos) << w[0];
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

// Switching holman <-> churchill_bernstein keeps the shared external velocity.
TEST(ForcedConvection, KindSwitchKeepsTheExternalVelocity) {
  UserData d;
  d.bed = SmallBed(ExternalHeatTransferType::HOLMAN, 1.0, 0.03);
  auto& layer = d.bed.layers_.layers_.front();
  dynamic_cast<ExternalHeatTransferModel_Holman&>(*layer.equations_.external_heat_transfer_model_).external_velocity =
      ExpressionParam{2.5};
  std::string err;
  ASSERT_TRUE(SetExternalHeatTransferKind(d, layer.id_, ExternalHeatTransferType::CHURCHILL_BERNSTEIN, err)) << err;
  auto* cb = dynamic_cast<ExternalHeatTransferModel_ChurchillBernstein*>(
      d.bed.layers_.layers_.front().equations_.external_heat_transfer_model_.get());
  ASSERT_TRUE(cb);
  EXPECT_DOUBLE_EQ(cb->external_velocity.Get(), 2.5);
}

// Switching churchill_bernstein -> mixed convection keeps the external velocity.
TEST(MixedConvection, KindSwitchKeepsTheExternalVelocity) {
  UserData d;
  d.bed = SmallBed(ExternalHeatTransferType::CHURCHILL_BERNSTEIN, 1.0, 0.03);
  auto& layer = d.bed.layers_.layers_.front();
  dynamic_cast<ExternalHeatTransferModel_ChurchillBernstein&>(*layer.equations_.external_heat_transfer_model_)
      .external_velocity = ExpressionParam{0.15};
  std::string err;
  ASSERT_TRUE(SetExternalHeatTransferKind(d, layer.id_, ExternalHeatTransferType::MIXED_CONVECTION_HORIZONTAL, err))
      << err;
  auto* m = dynamic_cast<ExternalHeatTransferModel_MixedConvection*>(
      d.bed.layers_.layers_.front().equations_.external_heat_transfer_model_.get());
  ASSERT_TRUE(m);
  EXPECT_FALSE(m->Vertical());
  EXPECT_DOUBLE_EQ(m->external_velocity.Get(), 0.15);
  EXPECT_DOUBLE_EQ(m->exponent.Get(), 3.0);
}

}  // namespace
}  // namespace adsim
