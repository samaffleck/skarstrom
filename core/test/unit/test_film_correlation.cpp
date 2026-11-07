/*
 * Skarstrom - pressure and temperature swing adsorption simulator
 * Copyright (c) 2024-2026 Sam Affleck and contributors
 *
 * Licensed under the terms in the LICENSE file at the repository root.
 */
// Packed-bed film correlations (adsim/equations/film_correlation.h):
//   Wakao and Funazkri (1978, Chem. Eng. Sci. 33, 1375) / Wakao, Kaguei and
//   Funazkri (1979, Chem. Eng. Sci. 34, 325): N = 2 + 1.1 X^(1/3) Re^0.6;
//   Gunn (1978, Int. J. Heat Mass Transfer 21, 467):
//     N = (7 - 10e + 5e^2)(1 + 0.7 Re^0.2 X^(1/3)) + (1.33 - 2.4e + 1.2e^2) Re^0.7 X^(1/3);
//   Ranz and Marshall (1952, Chem. Eng. Prog. 48, 141): N = 2 + 0.6 Re^(1/2) X^(1/3);
// with X = Pr for Nu and Sc for Sh. Checked here: each source's analytic
// limit, values worked out by hand (Python, double precision) from the
// formulas above, monotonicity, the layer setting's storage (deck archive,
// older archives, registry, builder) and the single-sphere deck warning.
// The solver-level check (the Film heat model and the Haynes coefficient
// both follow the layer's choice) is FilmCorrelation.HeatAndMassFollowTheLayerChoice
// in validation/test_structured_transfer.cpp.

#include <cmath>
#include <sstream>
#include <string>

#include <cereal/archives/json.hpp>
#include <gtest/gtest.h>

#include "adsim/equations/film_correlation.h"
#include "adsim/objects/deck_builder.h"
#include "adsim/objects/porous_media.h"
#include "adsim/objects/reactor.h"
#include "adsim/simulation.h"
#include "adsim/solver/deck_warnings.h"

namespace adsim {
namespace {

constexpr FilmCorrelation kAll[] = {FilmCorrelation::WakaoFunazkri, FilmCorrelation::Gunn,
                                    FilmCorrelation::RanzMarshall};

// Re -> 0: the conduction limit of a sphere, Nu = Sh = 2 (Wakao, Ranz-Marshall
// at any Pr; Gunn at e = 1). Gunn in a bed: 7 - 10 e + 5 e^2.
TEST(FilmCorrelation, StagnantLimits) {
  for (double x : {0.7, 1.0, 2.5}) {
    EXPECT_DOUBLE_EQ(PackedFilmNumber(FilmCorrelation::WakaoFunazkri, 0.0, x, 0.4), 2.0);
    EXPECT_DOUBLE_EQ(PackedFilmNumber(FilmCorrelation::RanzMarshall, 0.0, x, 0.4), 2.0);
    EXPECT_DOUBLE_EQ(PackedFilmNumber(FilmCorrelation::Gunn, 0.0, x, 1.0), 2.0);
    EXPECT_DOUBLE_EQ(PackedFilmNumber(FilmCorrelation::Gunn, 0.0, x, 0.4), 7.0 - 4.0 + 0.8);
  }
}

// Values worked out by hand from the formulas in the header comment.
TEST(FilmCorrelation, HandValues) {
  struct Case { FilmCorrelation c; double re, x, eps, expected; };
  const Case cases[] = {
    {FilmCorrelation::WakaoFunazkri, 10, 0.7, 0.4, 5.888290494654143},
    {FilmCorrelation::WakaoFunazkri, 100, 0.7, 0.4, 17.47956327116819},
    {FilmCorrelation::WakaoFunazkri, 1000, 2.0, 0.4, 89.44520636818697},
    {FilmCorrelation::RanzMarshall, 10, 0.7, 0.4, 3.684679393450839},
    {FilmCorrelation::RanzMarshall, 100, 0.7, 0.4, 7.327424010455604},
    {FilmCorrelation::RanzMarshall, 1000, 2.0, 0.4, 25.905321137950686},
    {FilmCorrelation::Gunn, 10, 0.7, 0.4, 10.04417436609237},
    {FilmCorrelation::Gunn, 100, 0.7, 0.4, 22.267000039731258},
    {FilmCorrelation::Gunn, 1000, 2.0, 0.4, 106.28356428729988},
    {FilmCorrelation::Gunn, 10, 0.7, 1.0, 4.548634206557106},
    {FilmCorrelation::Gunn, 100, 0.7, 1.0, 8.021847839042804},
    {FilmCorrelation::Gunn, 1000, 2.0, 1.0, 29.642076605493724},
  };
  for (const auto& c : cases)
    EXPECT_NEAR(PackedFilmNumber(c.c, c.re, c.x, c.eps), c.expected, 1e-12 * c.expected)
        << FilmCorrelationKinds()[static_cast<int>(c.c)] << " Re " << c.re << " X " << c.x << " eps " << c.eps;
}

// The default is the engine's Sherwood number from before the setting existed,
// to the last bit (decks keep their results).
TEST(FilmCorrelation, DefaultIsTheFormerSherwoodNumberBitForBit) {
  for (double re : {0.0, 0.37, 12.0, 450.0})
    for (double sc : {0.2, 0.8, 1.9})
      EXPECT_EQ(PackedFilmNumber(FilmCorrelation::WakaoFunazkri, re, sc, 0.4),
                2 + 1.1 * std::pow(sc, 0.3333333) * std::pow(re, 0.6));
  EXPECT_EQ(MixedAdsorbentLayer{}.film_correlation_, FilmCorrelation::WakaoFunazkri);
}

// Increasing in Re and in Pr (Sc) for every correlation; Gunn decreasing in
// the void fraction over the packed range at fixed Re (its packed-bed
// enhancement).
TEST(FilmCorrelation, Monotone) {
  for (auto c : kAll) {
    double prev = PackedFilmNumber(c, 0.0, 0.7, 0.4);
    for (double re = 0.5; re < 1e5; re *= 2.0) {
      const double now = PackedFilmNumber(c, re, 0.7, 0.4);
      EXPECT_GT(now, prev) << FilmCorrelationKinds()[static_cast<int>(c)] << " Re " << re;
      EXPECT_GT(PackedFilmNumber(c, re, 0.9, 0.4), now);
      prev = now;
    }
  }
  for (double re : {1.0, 100.0, 1e4}) {
    double prev = PackedFilmNumber(FilmCorrelation::Gunn, re, 0.7, 0.35);
    for (double e = 0.4; e <= 1.0 + 1e-12; e += 0.05) {
      const double now = PackedFilmNumber(FilmCorrelation::Gunn, re, 0.7, e);
      EXPECT_LT(now, prev) << "Re " << re << " eps " << e;
      prev = now;
    }
  }
}

std::string SaveLayer(const MixedAdsorbentLayer& layer) {
  std::ostringstream os;
  {
    cereal::JSONOutputArchive ar(os);
    ar(cereal::make_nvp("layer", layer));
  }
  return os.str();
}

MixedAdsorbentLayer LoadLayer(const std::string& json) {
  MixedAdsorbentLayer layer;
  std::istringstream is(json);
  cereal::JSONInputArchive ar(is);
  ar(cereal::make_nvp("layer", layer));
  return layer;
}

// The choice is saved with the layer; an archive older than v6 (without it)
// loads the default.
TEST(FilmCorrelation, StoredWithTheLayerAndOlderArchivesLoadTheDefault) {
  MixedAdsorbentLayer layer;
  layer.film_correlation_ = FilmCorrelation::Gunn;
  const std::string json = SaveLayer(layer);
  EXPECT_NE(json.find("\"Film_Correlation\": 1"), std::string::npos) << json;
  EXPECT_EQ(LoadLayer(json).film_correlation_, FilmCorrelation::Gunn);

  // The same archive as version 5 wrote it: no Film_Correlation field.
  std::string v5 = json;
  const std::string version = "\"cereal_class_version\": 6";
  const auto at = v5.find(version);
  ASSERT_NE(at, std::string::npos) << json;
  v5.replace(at, version.size(), "\"cereal_class_version\": 5");
  const auto field = v5.find(",\n        \"Film_Correlation\": 1");
  ASSERT_NE(field, std::string::npos) << json;
  v5.erase(field, std::string(",\n        \"Film_Correlation\": 1").size());
  EXPECT_EQ(LoadLayer(v5).film_correlation_, FilmCorrelation::WakaoFunazkri);
}

// The registry path layer[ID].film_correlation (by name or index) and the
// builder's kind names; a structured layer does not list it.
TEST(FilmCorrelation, RegistryAndBuilder) {
  Simulation sim;
  UserData& d = sim.Data();
  build::AddComponent(d, "CarbonDioxide");
  build::AddComponent(d, "Nitrogen");
  const auto ads = build::AddAdsorbent(d, "13X");
  const auto layer = build::AddLayer(d, "Bed");
  build::SetLayerComposition(d, layer, {{ads, 1.0}});
  sim.RefreshParameters();

  const std::string path = "layer[" + std::to_string(layer) + "].film_correlation";
  ASSERT_NE(sim.Parameters().Find(path), nullptr);
  EXPECT_EQ(sim.GetParameter(path), 0.0);
  sim.SetParameterText(path, "gunn");
  EXPECT_EQ(d.bed.layers_.Get(layer).film_correlation_, FilmCorrelation::Gunn);
  sim.SetParameter(path, 2.0);
  EXPECT_EQ(d.bed.layers_.Get(layer).film_correlation_, FilmCorrelation::RanzMarshall);
  EXPECT_ANY_THROW(sim.SetParameter(path, 3.0));

  EXPECT_EQ(build::FilmCorrelationKinds(),
            (std::vector<std::string>{"wakao_funazkri", "gunn", "ranz_marshall"}));
  build::SetFilmCorrelation(d, layer, "wakao_funazkri");
  EXPECT_EQ(d.bed.layers_.Get(layer).film_correlation_, FilmCorrelation::WakaoFunazkri);
  EXPECT_THROW(build::SetFilmCorrelation(d, layer, "colburn"), std::invalid_argument);

  build::SetLayerKind(d, layer, "laminate");
  sim.RefreshParameters();
  EXPECT_EQ(sim.Parameters().Find(path), nullptr);
}

// Ranz-Marshall on a packed layer: a warning (single sphere), not an error.
TEST(FilmCorrelation, RanzMarshallSingleSphereWarning) {
  Bed bed;
  const auto a = bed.layers_.Add("Inlet");
  const auto b = bed.layers_.Add("Main");
  EXPECT_TRUE(FilmCorrelationWarnings(bed).empty());
  bed.layers_.Get(a).film_correlation_ = FilmCorrelation::Gunn;
  EXPECT_TRUE(FilmCorrelationWarnings(bed).empty());
  bed.layers_.Get(b).film_correlation_ = FilmCorrelation::RanzMarshall;
  const auto w = FilmCorrelationWarnings(bed);
  ASSERT_EQ(w.size(), 1u);
  EXPECT_EQ(w[0].rfind("Warning: layer 'Main'", 0), 0u) << w[0];
  EXPECT_NE(w[0].find("single sphere"), std::string::npos) << w[0];
  // A structured layer does not use the setting.
  bed.layers_.Get(b).kind_ = LayerKind::Laminate;
  EXPECT_TRUE(FilmCorrelationWarnings(bed).empty());
}

}  // namespace
}  // namespace adsim
