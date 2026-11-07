// Structured layers (MixedAdsorbentLayer with a monolith or laminate kind):
// the geometry getters, the geometry checks, the builder (SetLayerKind,
// SetChannelShape), the parameter registry's per-kind fields, and
// serialisation (v4; older decks load as packed beds, bit for bit). The solved
// pressure drop is in validation/test_structured_layers.cpp.

#include <cmath>
#include <filesystem>
#include <set>
#include <sstream>
#include <string>

#include <cereal/archives/json.hpp>
#include <gtest/gtest.h>

#include "adsim/objects/deck_builder.h"
#include "adsim/objects/deck_io.h"
#include "adsim/objects/porous_media.h"
#include "adsim/objects/user_data.h"
#include "adsim/parameters/parameter_registry.h"
#include "adsim/simulation.h"
#include "adsim/solver/deck_warnings.h"
#include "adsim/equations/dispersion_model.h"
#include "adsim/equations/mass_transfer_model.h"

namespace adsim {
namespace {

  const std::filesystem::path kExamples = std::filesystem::current_path() / ".." / "examples" / "decks";
  constexpr double kPi = 3.14159265358979323846;

  MixedAdsorbentLayer Monolith(ChannelShape shape, double n, double w, double t_wc = 0.0) {
    MixedAdsorbentLayer l("M");
    l.kind_ = LayerKind::Monolith;
    l.structured_.channel_shape_ = shape;
    l.structured_.cell_density_ = n;
    l.structured_.wall_thickness_ = w;
    l.structured_.washcoat_thickness_ = t_wc;
    return l;
  }

  MixedAdsorbentLayer Laminate(double t, double g, double t_wc = 0.0) {
    MixedAdsorbentLayer l("L");
    l.kind_ = LayerKind::Laminate;
    l.structured_.sheet_thickness_ = t;
    l.structured_.gap_ = g;
    l.structured_.washcoat_thickness_ = t_wc;
    return l;
  }

  TEST(StructuredGeometry, DefaultsAreA400CpsiMonolithAndPacked) {
    MixedAdsorbentLayer l;
    EXPECT_EQ(l.kind_, LayerKind::Packed);
    EXPECT_FALSE(l.IsStructured());
    EXPECT_EQ(l.structured_.channel_shape_, ChannelShape::Square);
    EXPECT_DOUBLE_EQ(l.structured_.cell_density_.Get(), 6.2e5);
    EXPECT_DOUBLE_EQ(l.structured_.wall_thickness_.Get(), 1.65e-4);
    EXPECT_DOUBLE_EQ(l.structured_.sheet_thickness_.Get(), 2e-4);
    EXPECT_DOUBLE_EQ(l.structured_.gap_.Get(), 2e-4);
    EXPECT_DOUBLE_EQ(l.structured_.washcoat_thickness_.Get(), 0.0);
  }

  TEST(StructuredGeometry, SquareMonolith) {
    const double n = 6.2e5, w = 1.65e-4, t = 2e-5;
    const auto l = Monolith(ChannelShape::Square, n, w, t);
    const double p = 1.0 / std::sqrt(n), a = p - w - 2 * t, b = p - w;
    EXPECT_NEAR(l.GetOpenFraction(), a * a / (p * p), 1e-15);
    EXPECT_NEAR(l.GetHydraulicDiameter(), a, 1e-18);
    EXPECT_NEAR(l.GetSpecificSurfaceArea(), 4 * a / (p * p), 1e-9);
    EXPECT_DOUBLE_EQ(l.GetFrictionFactorReynolds(), 14.23);
    EXPECT_NEAR(l.GetPermeability(), (a * a / (p * p)) * a * a / (2 * 14.23), 1e-22);
    EXPECT_DOUBLE_EQ(l.GetInertialCoefficient(), 0.0);
    EXPECT_DOUBLE_EQ(l.GetDiffusionLength(), t);  // the washcoat
    EXPECT_NEAR(l.GetWashcoatFraction(), (b * b - a * a) / (p * p - a * a), 1e-14);
    // 400 cpsi, 6.5 mil walls: the textbook ~0.72 open fraction.
    const auto plain = Monolith(ChannelShape::Square, n, w);
    EXPECT_NEAR(plain.GetOpenFraction(), 0.76, 0.01);
    EXPECT_DOUBLE_EQ(plain.GetDiffusionLength(), w / 2);
    EXPECT_DOUBLE_EQ(plain.GetWashcoatFraction(), 0.0);
  }

  TEST(StructuredGeometry, CircularMonolith) {
    const double n = 6.2e5, w = 1.65e-4;
    const auto l = Monolith(ChannelShape::Circular, n, w);
    const double p = 1.0 / std::sqrt(n), a = p - w;
    EXPECT_NEAR(l.GetOpenFraction(), kPi * a * a / (4 * p * p), 1e-15);
    EXPECT_NEAR(l.GetHydraulicDiameter(), a, 1e-18);
    EXPECT_NEAR(l.GetSpecificSurfaceArea(), kPi * a / (p * p), 1e-9);
    EXPECT_DOUBLE_EQ(l.GetFrictionFactorReynolds(), 16.0);
  }

  TEST(StructuredGeometry, Laminate) {
    const double t = 2e-4, g = 3e-4, t_wc = 5e-5;
    const auto l = Laminate(t, g, t_wc);
    EXPECT_NEAR(l.GetOpenFraction(), (g - 2 * t_wc) / (g + t), 1e-15);
    EXPECT_NEAR(l.GetHydraulicDiameter(), 2 * (g - 2 * t_wc), 1e-18);
    EXPECT_NEAR(l.GetSpecificSurfaceArea(), 2 / (g + t), 1e-9);
    EXPECT_DOUBLE_EQ(l.GetFrictionFactorReynolds(), 24.0);
    EXPECT_DOUBLE_EQ(l.GetDiffusionLength(), t_wc);
    EXPECT_NEAR(l.GetWashcoatFraction(), 2 * t_wc / (t + 2 * t_wc), 1e-15);
    EXPECT_DOUBLE_EQ(Laminate(t, g).GetDiffusionLength(), t / 2);
  }

  TEST(StructuredGeometry, PackedKeepsItsMomentumModel) {
    MixedAdsorbentLayer l;
    l.eb_ = 0.4;
    l.particle_diameter_ = 2e-3;
    l.momentum_model_ = MomentumModel::Ergun;
    EXPECT_DOUBLE_EQ(l.GetPermeability(), GetLayerPermeability(0.4, 2e-3, MomentumModel::Ergun));
    EXPECT_DOUBLE_EQ(l.GetInertialCoefficient(), GetLayerInertialCoefficient(0.4, 2e-3, MomentumModel::Ergun));
    EXPECT_DOUBLE_EQ(l.GetOpenFraction(), 0.4);
    EXPECT_NEAR(l.GetHydraulicDiameter(), (2.0 / 3.0) * 0.4 * 2e-3 / 0.6, 1e-18);
    EXPECT_NEAR(l.GetSpecificSurfaceArea(), 6 * 0.6 / 2e-3, 1e-9);
    EXPECT_DOUBLE_EQ(l.GetDiffusionLength(), 1e-3);
    // The momentum model is ignored once the layer is structured.
    l.kind_ = LayerKind::Monolith;
    EXPECT_DOUBLE_EQ(l.GetInertialCoefficient(), 0.0);
  }

  // The open frontal area takes the bed void fraction's place: gas holdup and
  // bulk density follow the geometry.
  TEST(StructuredGeometry, OpenAreaIsTheVoidFraction) {
    AdsorbentDirectory dir;
    const AdsorbentID id = dir.Add("S");
    dir.Get(id).solid_density_ = 1200.0;
    dir.Get(id).ep_ = 0.3;
    auto l = Monolith(ChannelShape::Square, 6.2e5, 1.65e-4);
    l.composition_[id] = 1.0;
    l.eb_ = 0.37;  // unused once structured
    const double ofa = l.GetOpenFraction();
    EXPECT_DOUBLE_EQ(l.GetVoidFraction(), ofa);
    EXPECT_DOUBLE_EQ(l.GetTotalVoidFraction(dir), ofa + (1 - ofa) * 0.3);
    EXPECT_DOUBLE_EQ(l.GetBulkDensity(dir), (1 - ofa) * 1200.0);
    EXPECT_DOUBLE_EQ(l.GetAdsorbentBulkDensities(dir).at(id), (1 - ofa) * 1200.0);
    l.kind_ = LayerKind::Packed;
    EXPECT_DOUBLE_EQ(l.GetVoidFraction(), 0.37);
    EXPECT_DOUBLE_EQ(l.GetBulkDensity(dir), (1 - 0.37) * 1200.0);
  }

  TEST(StructuredGeometry, SummaryForFrontEnds) {
    const auto l = Monolith(ChannelShape::Square, 6.2e5, 1.65e-4);
    const double area = 0.25 * kPi * 0.05 * 0.05;
    const LayerGeometry s = LayerGeometrySummary(l, area);
    EXPECT_EQ(s.kind, LayerKind::Monolith);
    EXPECT_DOUBLE_EQ(s.open_fraction, l.GetOpenFraction());
    EXPECT_DOUBLE_EQ(s.hydraulic_diameter, l.GetHydraulicDiameter());
    EXPECT_DOUBLE_EQ(s.specific_surface_area, l.GetSpecificSurfaceArea());
    EXPECT_DOUBLE_EQ(s.diffusion_length, l.GetDiffusionLength());
    EXPECT_DOUBLE_EQ(s.channel_count, 6.2e5 * area);
    EXPECT_DOUBLE_EQ(s.friction_factor_reynolds, 14.23);
    EXPECT_TRUE(std::isnan(LayerGeometrySummary(Laminate(2e-4, 2e-4), area).channel_count));
  }

  TEST(StructuredGeometry, ChecksRejectGeometryWithoutAnOpenChannel) {
    EXPECT_TRUE(CheckLayerGeometry(MixedAdsorbentLayer{}).empty());  // packed: nothing to check
    EXPECT_TRUE(CheckLayerGeometry(Monolith(ChannelShape::Square, 6.2e5, 1.65e-4, 2e-5)).empty());
    EXPECT_TRUE(CheckLayerGeometry(Laminate(2e-4, 2e-4, 5e-5)).empty());

    // Walls and washcoat fill the pitch (p = 1.27 mm).
    auto e = CheckLayerGeometry(Monolith(ChannelShape::Square, 6.2e5, 1.0e-3, 2e-4));
    ASSERT_EQ(e.size(), 1u);
    EXPECT_NE(e[0].find("leaving no open channel"), std::string::npos) << e[0];
    EXPECT_NE(e[0].find("Layer \"M\" (Monolith)"), std::string::npos) << e[0];
    // Washcoat fills the gap.
    e = CheckLayerGeometry(Laminate(2e-4, 2e-4, 1e-4));
    ASSERT_EQ(e.size(), 1u);
    EXPECT_NE(e[0].find("fills the gap"), std::string::npos) << e[0];
    // Non-positive inputs.
    EXPECT_FALSE(CheckLayerGeometry(Monolith(ChannelShape::Square, 0.0, 1.65e-4)).empty());
    EXPECT_FALSE(CheckLayerGeometry(Monolith(ChannelShape::Square, 6.2e5, 0.0)).empty());
    EXPECT_FALSE(CheckLayerGeometry(Monolith(ChannelShape::Square, 6.2e5, 1.65e-4, -1e-6)).empty());
    EXPECT_FALSE(CheckLayerGeometry(Laminate(0.0, 2e-4)).empty());
    EXPECT_FALSE(CheckLayerGeometry(Laminate(2e-4, -1e-4)).empty());
  }

  // Switching kind keeps the composition, domain, name and both sets of
  // fields, and the deck stays valid; the registry shows only the fields of
  // the layer's kind.
  TEST(StructuredGeometry, BuilderSwitchesKindAndRegistryFollows) {
    Simulation sim = Simulation::FromFile((kExamples / "Perez2019_A.json").string());
    UserData& d = sim.Data();
    auto& layer = d.bed.layers_.layers_.front();
    const LayerID id = layer.id_;
    const auto composition = layer.composition_;
    const double length = layer.domain_.length.Get();
    const double eb = layer.eb_.Get();
    const std::string name = layer.name_;
    const std::string lp = "layer[" + std::to_string(id) + "].";

    auto paths = [&]() {
      std::set<std::string> s;
      sim.RefreshParameters();
      for (const auto& p : sim.Parameters().All()) s.insert(p.Info().path);
      return s;
    };
    auto packed = paths();
    EXPECT_TRUE(packed.count(lp + "bed_void_fraction"));
    EXPECT_TRUE(packed.count(lp + "momentum_model"));
    EXPECT_FALSE(packed.count(lp + "cell_density"));

    build::SetLayerKind(d, id, "monolith");
    EXPECT_EQ(layer.kind_, LayerKind::Monolith);
    EXPECT_EQ(layer.composition_, composition);
    EXPECT_DOUBLE_EQ(layer.domain_.length.Get(), length);
    EXPECT_EQ(layer.name_, name);
    std::string err;
    EXPECT_TRUE(sim.Check(err)) << err;
    auto mono = paths();
    for (const char* f : {"channel_shape", "cell_density", "wall_thickness", "washcoat_thickness"})
      EXPECT_TRUE(mono.count(lp + f)) << f;
    for (const char* f : {"bed_void_fraction", "momentum_model", "sheet_thickness", "gap"})
      EXPECT_FALSE(mono.count(lp + f)) << f;
    // Kept listed (an existing path), though no accepted closure reads it.
    EXPECT_TRUE(mono.count(lp + "particle_diameter"));
    // The channel heat closures were installed: channel conduction and a
    // constant wall coefficient.
    EXPECT_NE(dynamic_cast<HeatDispersionModelChannel*>(layer.equations_.heat_dispersion_model_.get()), nullptr);
    EXPECT_NE(dynamic_cast<ConstantInternalHeatTransferModel*>(layer.equations_.internal_heat_transfer_model_.get()),
              nullptr);
    sim.SetParameter(lp + "cell_density", 9.3e5);  // 600 cpsi
    EXPECT_DOUBLE_EQ(layer.structured_.cell_density_.Get(), 9.3e5);
    build::SetChannelShape(d, id, "circular");
    EXPECT_EQ(layer.structured_.channel_shape_, ChannelShape::Circular);
    EXPECT_THROW(build::SetChannelShape(d, id, "hexagonal"), std::invalid_argument);

    build::SetLayerKind(d, id, "laminate");
    auto lam = paths();
    for (const char* f : {"sheet_thickness", "gap", "washcoat_thickness"}) EXPECT_TRUE(lam.count(lp + f)) << f;
    for (const char* f : {"channel_shape", "cell_density", "wall_thickness", "bed_void_fraction"})
      EXPECT_FALSE(lam.count(lp + f)) << f;
    err.clear();
    EXPECT_TRUE(sim.Check(err)) << err;

    // Back to packed: the packed fields were kept, and the packed default
    // (constant) heat dispersion replaces channel conduction.
    build::SetLayerKind(d, id, "packed");
    EXPECT_DOUBLE_EQ(layer.eb_.Get(), eb);
    EXPECT_NE(dynamic_cast<HeatDispersionModelConstant*>(layer.equations_.heat_dispersion_model_.get()), nullptr);
    err.clear();
    EXPECT_TRUE(sim.Check(err)) << err;
    EXPECT_DOUBLE_EQ(layer.structured_.cell_density_.Get(), 9.3e5);  // kept for next time
    EXPECT_THROW(build::SetLayerKind(d, id, "foam"), std::invalid_argument);
    EXPECT_EQ(build::LayerKinds(), (std::vector<std::string>{"packed", "monolith", "laminate"}));
    EXPECT_EQ(build::ChannelShapeKinds(), (std::vector<std::string>{"square", "circular"}));
  }

  std::string Save(const MixedAdsorbentLayer& l) {
    std::ostringstream os;
    {
      cereal::JSONOutputArchive ar(os);
      ar(cereal::make_nvp("layer", l));
    }
    return os.str();
  }

  TEST(StructuredGeometry, SerialisationRoundTripAndOldDecksLoadPacked) {
    auto l = Laminate(2.5e-4, 3e-4, 4e-5);
    l.structured_.channel_shape_ = ChannelShape::Circular;
    l.structured_.cell_density_.SetExpr("2 * 3.1e5");
    std::istringstream is(Save(l));
    MixedAdsorbentLayer back;
    {
      cereal::JSONInputArchive ar(is);
      ar(cereal::make_nvp("layer", back));
    }
    EXPECT_EQ(back.kind_, LayerKind::Laminate);
    EXPECT_EQ(back.structured_.channel_shape_, ChannelShape::Circular);
    EXPECT_EQ(back.structured_.cell_density_.GetExpr(), "2 * 3.1e5");
    EXPECT_DOUBLE_EQ(back.structured_.sheet_thickness_.Get(), 2.5e-4);
    EXPECT_DOUBLE_EQ(back.structured_.gap_.Get(), 3e-4);
    EXPECT_DOUBLE_EQ(back.structured_.washcoat_thickness_.Get(), 4e-5);

    // A bundled deck (written before v4) loads as packed with the default
    // geometry and the same permeability as its momentum model gives.
    UserData d;
    LoadFromJson(d, (kExamples / "Perez2019_A.json").string());
    for (const auto& layer : d.bed.layers_.layers_) {
      EXPECT_EQ(layer.kind_, LayerKind::Packed);
      EXPECT_DOUBLE_EQ(layer.structured_.cell_density_.Get(), 6.2e5);
      EXPECT_EQ(layer.GetPermeability(), GetLayerPermeability(layer.eb_, layer.particle_diameter_, layer.momentum_model_));
    }
  }

  // Structured layers warn about LTE only; a packed deck gets no warning.
  TEST(StructuredGeometry, WarningForLte) {
    UserData d;
    LoadFromJson(d, (kExamples / "Perez2019_A.json").string());
    EXPECT_TRUE(StructuredLayerWarnings(d.adsorbent_directory, d.bed).empty());
    auto& layer = d.bed.layers_.layers_.front();
    layer.kind_ = LayerKind::Monolith;
    d.bed.equations_.thermal_equilibrium_mode_ = ThermalEquilibriumMode::LOCAL;
    auto w = StructuredLayerWarnings(d.adsorbent_directory, d.bed);
    ASSERT_EQ(w.size(), 1u);
    EXPECT_NE(w[0].find("LTNE) is recommended"), std::string::npos) << w[0];

    d.bed.equations_.thermal_equilibrium_mode_ = ThermalEquilibriumMode::NON_LOCAL;
    w = StructuredLayerWarnings(d.adsorbent_directory, d.bed);
    for (const auto& s : w) ADD_FAILURE() << s;
  }

  // The packed-bed heat closures, which read the particle diameter, are deck
  // errors on a structured layer now that channel closures replace them; the
  // channel conduction model is an error on a packed bed.
  TEST(StructuredGeometry, DeckChecksRefuseHeatClosuresOfTheOtherKind) {
    Simulation sim = Simulation::FromFile((kExamples / "Perez2019_A.json").string());
    UserData& d = sim.Data();
    const LayerID id = d.bed.layers_.layers_.front().id_;
    auto& layer = d.bed.layers_.Get(id);
    std::string err;
    ASSERT_TRUE(sim.Check(err)) << err;

    build::SetHeatDispersion(d, id, "channel");
    err.clear();
    EXPECT_FALSE(sim.Check(err));
    EXPECT_NE(err.find("is a packed bed but uses channel heat conduction"), std::string::npos) << err;

    build::SetLayerKind(d, id, "monolith");
    err.clear();
    EXPECT_TRUE(sim.Check(err)) << err;
    struct Case { const char* op; const char* kind; const char* says; };
    for (const Case& c : {Case{"heat", "ruthven", "Ruthven's packed-bed heat dispersion"},
                          Case{"heat", "yagi_kunii", "Yagi and Kunii's packed-bed heat dispersion"},
                          Case{"wall", "dixon", "Dixon's packed-bed wall heat transfer coefficient"},
                          Case{"wall", "leva", "Leva's packed-bed wall heat transfer coefficient"}}) {
      if (std::string(c.op) == "heat") build::SetHeatDispersion(d, id, c.kind);
      else build::SetInternalHeatTransfer(d, id, c.kind);
      err.clear();
      EXPECT_FALSE(sim.Check(err)) << c.kind;
      EXPECT_NE(err.find(std::string("is a Monolith but uses ") + c.says), std::string::npos) << err;
      build::SetHeatDispersion(d, id, "channel");
      build::SetInternalHeatTransfer(d, id, "constant");
    }
    err.clear();
    EXPECT_TRUE(sim.Check(err)) << err;

    // A packed layer switched to a structured kind with Dixon gets a constant
    // h_i; a constant one keeps its value.
    build::SetLayerKind(d, id, "packed");
    build::SetInternalHeatTransfer(d, id, "dixon");
    build::SetLayerKind(d, id, "laminate");
    EXPECT_NE(dynamic_cast<ConstantInternalHeatTransferModel*>(layer.equations_.internal_heat_transfer_model_.get()),
              nullptr);
    build::SetLayerKind(d, id, "packed");
    auto hi = std::make_unique<ConstantInternalHeatTransferModel>();
    hi->heat_transfer_coefficient_ = 42.0;
    layer.equations_.internal_heat_transfer_model_ = std::move(hi);
    build::SetLayerKind(d, id, "monolith");
    auto* kept = dynamic_cast<ConstantInternalHeatTransferModel*>(layer.equations_.internal_heat_transfer_model_.get());
    ASSERT_NE(kept, nullptr);
    EXPECT_DOUBLE_EQ(kept->heat_transfer_coefficient_.Get(), 42.0);
    EXPECT_EQ(build::HeatDispersionKinds(), (std::vector<std::string>{"constant", "ruthven", "yagi_kunii", "channel"}));

    // The channel model survives a save and load.
    std::string json;
    SaveToString(d, json);
    UserData back;
    LoadDeckFlexible(back, json);
    EXPECT_EQ(back.bed.layers_.layers_.front().equations_.heat_dispersion_model_->GetModelType(),
              AxialHeatDispersionType::CHANNEL);
  }

}  // namespace
}  // namespace adsim
