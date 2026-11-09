/*
 * Skarstrom - pressure and temperature swing adsorption simulator
 * Copyright (c) 2024-2026 Sam Affleck and contributors
 *
 * Licensed under the terms in the LICENSE file at the repository root.
 */
// The packed particle shape of a layer (MixedAdsorbentLayer::particle_shape_):
// the LDF diffusion geometry per shape, the generic coefficient's terms for a
// cylinder against closed forms written out here (k = 8 D / r_f^2 and
// V/S = r_f / 2 for an infinite cylinder of radius r_f = d_p / 2; Glueckauf's
// parabolic profile, see docs/kinetic-models.md), the storage (deck archive,
// archives older than v7 load spheres), the registry path, the deck edit and
// the deck view field.

#include <cmath>
#include <sstream>
#include <string>

#include <cereal/archives/json.hpp>
#include <gtest/gtest.h>

#include "adsim/equations/equation.h"
#include "adsim/equations/mass_transfer_model.h"
#include "adsim/math/constants.h"
#include "adsim/objects/deck_builder.h"
#include "adsim/objects/deck_edit.h"
#include "adsim/objects/deck_view.h"
#include "adsim/objects/fluid.h"
#include "adsim/objects/porous_media.h"
#include "adsim/simulation.h"

namespace adsim {
namespace {

constexpr int kCells = 5;
constexpr AdsorbentID kAds = 1;
const std::string kCO2 = "CarbonDioxide";
constexpr double kT = 298.15;   // [K]
constexpr double kDp = 2.0e-3;  // particle (cylinder: extrudate) diameter [m]
constexpr double kRf = kDp / 2.0;
constexpr double kEp = 0.35;    // macropore void fraction [-]
constexpr double kTau = 3.0;    // tortuosity [-]
constexpr double kDm = 1.6e-5;  // molecular diffusivity of the state [m2/s]
constexpr double kSh = 7.5;     // Sherwood number of the state [-]

using Generic = MassTransferCoefficientSeriesResistance;
using Geometry = MixedAdsorbentLayer::DiffusionGeometry;

TEST(ParticleShape, DiffusionGeometryPerShape) {
  MixedAdsorbentLayer layer;
  layer.particle_diameter_ = kDp;
  EXPECT_EQ(layer.particle_shape_, ParticleShape::Sphere);  // the default

  Geometry g = layer.GetDiffusionGeometry();
  EXPECT_EQ(g.shape, MixedAdsorbentLayer::Shape::Sphere);
  EXPECT_EQ(g.c, 60.0);  // on d_p (15 on r_p)
  EXPECT_EQ(g.length, kDp);
  EXPECT_EQ(g.volume_to_surface, kDp / 6.0);

  layer.particle_shape_ = ParticleShape::Cylinder;
  g = layer.GetDiffusionGeometry();
  EXPECT_EQ(g.shape, MixedAdsorbentLayer::Shape::Cylinder);
  EXPECT_EQ(g.c, 8.0);  // on r_f
  EXPECT_EQ(g.length, kRf);
  EXPECT_EQ(g.volume_to_surface, kRf / 2.0);
  EXPECT_EQ(layer.GetDiffusionLength(), kRf);
  // a_v stays the equivalent-sphere 6 (1 - eb) / d_p (the solver's value).
  EXPECT_DOUBLE_EQ(layer.GetSpecificSurfaceArea(), 6.0 * (1.0 - layer.eb_.Get()) / kDp);

  // A structured layer is a slab whatever its stored particle shape.
  layer.kind_ = LayerKind::Laminate;
  g = layer.GetDiffusionGeometry();
  EXPECT_EQ(g.shape, MixedAdsorbentLayer::Shape::Slab);
  EXPECT_EQ(g.c, 3.0);
  EXPECT_EQ(g.length, layer.GetDiffusionLength());
  EXPECT_EQ(g.volume_to_surface, g.length);
}

// A uniform state for one coefficient evaluation.
struct State {
  Equations eq;
  MixedAdsorbentLayer layer{"bed"};
  Fluid fluid;

  explicit State(ParticleShape shape) {
    eq.fluid_temperature.var.assign(kCells, kT);
    eq.solid_temperature.var.assign(kCells, kT);
    eq.pressure.var.assign(kCells, 1e5);
    eq.sherwood_number.var.assign(kCells, kSh);
    Equation yi{"Mole Fraction", CellType::CellCentre, EquationType::Algebraic, nullptr};
    yi.var.assign(kCells, 0.15);
    eq.mole_fraction.emplace(kCO2, std::move(yi));
    IndependentEquation dm{"Molecular Diffusion", CellType::CellCentre, nullptr};
    dm.var.assign(kCells, kDm);
    eq.molecular_diffusion.emplace(kCO2, std::move(dm));
    IndependentEquation k{"MTC", CellType::CellCentre, nullptr};
    k.var.assign(kCells, -1.0);
    eq.mass_transfer_coefficient[kAds].emplace(kCO2, std::move(k));
    IndependentEquation q{"q*", CellType::CellCentre, nullptr};
    q.var.assign(kCells, 0.0);
    eq.eq_ads_concentration[kAds].emplace(kCO2, std::move(q));
    layer.domain_.start_index = 1;
    layer.domain_.end_index = kCells - 2;
    layer.particle_diameter_ = kDp;
    layer.particle_shape_ = shape;
    fluid.molar_mass_[kCO2] = 0.04401;
  }

  double Update(IMassTransferCoefficientModel& m) {
    m.UpdateMassTransferCoefficient(eq, kAds, fluid, kCO2, layer, nullptr);
    const auto& k = eq.mass_transfer_coefficient.at(kAds).at(kCO2).var;
    for (int n = layer.domain_.start_index; n <= layer.domain_.end_index; ++n)
      EXPECT_EQ(k[n], k[layer.domain_.start_index]);
    return k[layer.domain_.start_index];
  }
};

Generic Model(bool film, bool macropore, bool knudsen, bool micropore, CapacityFactorMode mode,
              double fixed_k = 0.0) {
  Generic m;
  MassTransferResistances r;
  r.film = film;
  r.macropore = macropore;
  r.knudsen = knudsen;
  r.micropore = micropore;
  r.lumped = false;
  m.resistances_ = r;
  m.capacity_mode_ = static_cast<int>(mode);
  m.capacity_factor_ = fixed_k;
  m.tortuosity_ = kTau;
  m.macropore_void_fraction_ = kEp;
  m.molecular_diffusion_coefficient_ = 0.0;  // the state's D_m
  return m;
}

double Rel(double a, double b) { return std::abs(a - b) / std::abs(b); }

// Film only: 1/k = (V/S) / k_f, k_f = Sh D_m / d_p (the packed correlation on
// d_p). Sphere V/S = r_p / 3: k = 6 Sh D / d_p^2. Cylinder V/S = r_f / 2:
// k = 4 Sh D / d_p^2, so k is proportional to Sh.
TEST(ParticleShape, CylinderFilmCoefficient) {
  auto m = Model(true, false, false, false, CapacityFactorMode::NONE);
  State sphere(ParticleShape::Sphere), cyl(ParticleShape::Cylinder);
  EXPECT_LT(Rel(sphere.Update(m), 6.0 * kSh * kDm / (kDp * kDp)), 1e-14);
  const double k = cyl.Update(m);
  EXPECT_LT(Rel(k, 4.0 * kSh * kDm / (kDp * kDp)), 1e-14);
  cyl.eq.sherwood_number.var.assign(kCells, 2.0 * kSh);
  EXPECT_LT(Rel(cyl.Update(m), 2.0 * k), 1e-14);
}

// Macropore only: k = 8 eps_p D / (tau r_f^2), times 1/K with a capacity
// factor; mode none (the macropore-alone path) and fixed K (the series path).
TEST(ParticleShape, CylinderMacroporeCoefficient) {
  const double expected = 8.0 * kEp * kDm / (kTau * kRf * kRf);
  State cyl(ParticleShape::Cylinder);
  auto none = Model(false, true, false, false, CapacityFactorMode::NONE);
  EXPECT_LT(Rel(cyl.Update(none), expected), 1e-14);
  auto fixed = Model(false, true, false, false, CapacityFactorMode::FIXED, 2.5);
  EXPECT_LT(Rel(cyl.Update(fixed), expected / 2.5), 1e-14);

  State sphere(ParticleShape::Sphere);
  EXPECT_LT(Rel(sphere.Update(none), 15.0 * kEp * kDm / (kTau * 0.25 * kDp * kDp)), 1e-14);
}

// Film and macropore in series; Knudsen enters through D; the micropore term
// is on the crystal radius and does not follow the particle shape.
TEST(ParticleShape, CylinderSeriesTerms) {
  State cyl(ParticleShape::Cylinder);
  const double kf = kSh * kDm / kDp;
  const double film = kRf / (2.0 * kf);
  auto fm = Model(true, true, false, false, CapacityFactorMode::NONE);
  EXPECT_LT(Rel(cyl.Update(fm), 1.0 / (film + kTau * kRf * kRf / (8.0 * kEp * kDm))), 1e-14);

  fm.resistances_.knudsen = true;
  fm.pore_radius_ = 5e-9;
  const double R = constants::kIdealGasConstant;
  const double Dk = (2.0 / 3.0) * 5e-9 * std::sqrt(8.0 * R * kT / (constants::kPI * 0.04401));
  const double Db = 1.0 / (1.0 / Dk + 1.0 / kDm);
  EXPECT_LT(Rel(cyl.Update(fm), 1.0 / (film + kTau * kRf * kRf / (8.0 * kEp * Db))), 1e-14);

  auto micro = Model(false, false, false, true, CapacityFactorMode::NONE);
  State sphere(ParticleShape::Sphere);
  EXPECT_EQ(cyl.Update(micro), sphere.Update(micro));
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
  std::istringstream is(json);
  cereal::JSONInputArchive ar(is);
  MixedAdsorbentLayer layer;
  ar(cereal::make_nvp("layer", layer));
  return layer;
}

// Saved with the layer (v7); a v6 archive (no Particle_Shape) loads a sphere.
TEST(ParticleShape, StoredWithTheLayerAndOlderArchivesLoadSpheres) {
  MixedAdsorbentLayer layer;
  layer.particle_shape_ = ParticleShape::Cylinder;
  const std::string json = SaveLayer(layer);
  EXPECT_NE(json.find("\"Particle_Shape\": 1"), std::string::npos) << json;
  EXPECT_EQ(LoadLayer(json).particle_shape_, ParticleShape::Cylinder);

  std::string v6 = json;
  const std::string version = "\"cereal_class_version\": 7";
  const auto at = v6.find(version);
  ASSERT_NE(at, std::string::npos) << json;
  v6.replace(at, version.size(), "\"cereal_class_version\": 6");
  const std::string field = ",\n        \"Particle_Shape\": 1";
  const auto f = v6.find(field);
  ASSERT_NE(f, std::string::npos) << json;
  v6.erase(f, field.size());
  EXPECT_EQ(LoadLayer(v6).particle_shape_, ParticleShape::Sphere);
}

// The registry path layer[ID].particle_shape (by name or index), the deck
// edit and the deck view field; a structured layer lists no path and shows
// null.
TEST(ParticleShape, RegistryEditAndDeckView) {
  Simulation sim;
  UserData& d = sim.Data();
  build::AddComponent(d, "CarbonDioxide");
  build::AddComponent(d, "Nitrogen");
  const auto ads = build::AddAdsorbent(d, "13X");
  const auto layer = build::AddLayer(d, "Bed");
  build::SetLayerComposition(d, layer, {{ads, 1.0}});
  sim.RefreshParameters();

  const std::string path = "layer[" + std::to_string(layer) + "].particle_shape";
  ASSERT_NE(sim.Parameters().Find(path), nullptr);
  EXPECT_EQ(sim.GetParameter(path), 0.0);
  sim.SetParameterText(path, "cylinder");
  EXPECT_EQ(d.bed.layers_.Get(layer).particle_shape_, ParticleShape::Cylinder);
  sim.SetParameter(path, 0.0);
  EXPECT_EQ(d.bed.layers_.Get(layer).particle_shape_, ParticleShape::Sphere);
  EXPECT_ANY_THROW(sim.SetParameter(path, 2.0));

  EXPECT_EQ(ParticleShapeKinds(), (std::vector<std::string>{"sphere", "cylinder"}));
  std::string err;
  EXPECT_TRUE(SetParticleShape(d, layer, "cylinder", err)) << err;
  EXPECT_EQ(d.bed.layers_.Get(layer).particle_shape_, ParticleShape::Cylinder);
  EXPECT_FALSE(SetParticleShape(d, layer, "cube", err));
  EXPECT_NE(err.find("sphere, cylinder"), std::string::npos) << err;
  EXPECT_EQ(d.bed.layers_.Get(layer).particle_shape_, ParticleShape::Cylinder);
  EXPECT_FALSE(SetParticleShape(d, 999, "sphere", err));

  std::string view = DeckView(d);
  EXPECT_NE(view.find("\"particle_shape\":\"cylinder\""), std::string::npos);

  build::SetLayerKind(d, layer, "laminate");
  sim.RefreshParameters();
  EXPECT_EQ(sim.Parameters().Find(path), nullptr);
  view = DeckView(d);
  EXPECT_NE(view.find("\"particle_shape\":null"), std::string::npos);
}

}  // namespace
}  // namespace adsim
