/*
 * Skarstrom - pressure and temperature swing adsorption simulator
 * Copyright (c) 2024-2026 Sam Affleck and contributors
 *
 * Licensed under the terms in the LICENSE file at the repository root.
 */
// The registry leaves of the generic mass-transfer coefficient model: one
// Boolean per resistance (set by path through the deck-edit cascades), the
// capacity-factor mode, the parameters of the enabled terms only, the slow
// site's leaves only with a slow site, the former arrhenius leaf name; and
// the deck view's kinetics fields (switches, mode, preset, equations).

#include <gtest/gtest.h>

#include <algorithm>
#include <string>

#include <cereal/external/rapidjson/document.h>

#include "adsim/objects/deck_builder.h"
#include "adsim/objects/deck_edit.h"
#include "adsim/objects/deck_view.h"
#include "adsim/objects/manifests.h"
#include "adsim/objects/model_kinds.h"
#include "adsim/parameters/parameter_registry.h"
#include "adsim/simulation.h"

namespace adsim {
namespace {

const std::string kCO2 = "CarbonDioxide";

struct Deck {
  Simulation sim;
  AdsorbentID a = 0;
  std::string kp;
  explicit Deck(const char* rate = "ldf", const char* coefficient = "constant") {
    UserData& d = sim.Data();
    build::AddComponent(d, kCO2);
    a = build::AddAdsorbent(d, "13X");
    build::SetIsotherm(d, a, "henry", {kCO2});
    build::SetRateModel(d, a, kCO2, rate, coefficient);
    sim.RefreshParameters();
    kp = "adsorbent[" + std::to_string(a) + "].kinetics[" + kCO2 + "].";
  }
  MassTransferCoefficientSeriesResistance& Fast() {
    return dynamic_cast<MassTransferCoefficientSeriesResistance&>(
        *sim.Data().adsorbent_directory.Get(a).adsorption_source_.mass_transfer_models_.at(kCO2)->mtc_model_);
  }
  MassTransferCoefficientSeriesResistance* Slow() {
    auto* m = sim.Data().adsorbent_directory.Get(a).adsorption_source_.mass_transfer_models_.at(kCO2)->SlowCoefficientModel();
    return dynamic_cast<MassTransferCoefficientSeriesResistance*>(m);
  }
  bool Has(const std::string& leaf) { return sim.Parameters().Find(kp + leaf) != nullptr; }
};

// The registry lists, for one site, exactly what the model's switches say.
void ExpectAgrees(Deck& deck, MassTransferCoefficientSeriesResistance& m, const std::string& pre) {
  const KindDescriptor& kind = ResistanceKind();
  for (const auto& r : kind.resistances) {
    const bool listed = r.parent.empty() || m.Enabled(r.parent);
    ASSERT_EQ(deck.Has(pre + "resistances." + r.tag), listed) << pre << r.tag;
    if (listed) EXPECT_EQ(deck.sim.Parameters().Get(deck.kp + pre + "resistances." + r.tag).GetValue() != 0.0, m.Enabled(r.tag));
    for (const auto& p : r.parameters) EXPECT_EQ(deck.Has(pre + p), m.Enabled(r.tag)) << pre << p;
  }
  EXPECT_EQ(deck.Has(pre + "capacity_factor_mode"), m.UsesCapacityFactor()) << pre;
  if (m.UsesCapacityFactor())
    EXPECT_EQ(static_cast<int>(deck.sim.Parameters().Get(deck.kp + pre + "capacity_factor_mode").GetValue()),
              static_cast<int>(m.CapacityMode()));
  const auto cap = m.CapacityParameters();
  for (const char* p : {"concentration_floor", "capacity_factor"}) {
    const bool shown = std::find(cap.begin(), cap.end(), p) != cap.end();
    EXPECT_EQ(deck.Has(pre + p), shown) << pre << p;
  }
}

}  // namespace

TEST(ParameterRegistry, ResistanceSwitches) {
  Deck deck;
  auto& m = deck.Fast();
  ExpectAgrees(deck, m, "");
  EXPECT_FALSE(deck.Has("tortuosity"));
  EXPECT_FALSE(deck.Has("slow_resistances.lumped")) << "no slow site on an LDF model";

  // Switches set by path go through the cascade (its rules apply).
  deck.sim.SetParameterText(deck.kp + "resistances.macropore", "true");
  EXPECT_TRUE(m.Enabled("macropore"));
  ExpectAgrees(deck, m, "");
  EXPECT_TRUE(deck.Has("tortuosity"));
  EXPECT_TRUE(deck.Has("resistances.knudsen"));
  EXPECT_TRUE(deck.Has("resistances.molecular"));
  EXPECT_TRUE(m.Enabled("molecular")) << "macropore on switches molecular on";
  EXPECT_TRUE(deck.Has("molecular_diffusion_coefficient"));
  deck.sim.SetParameter(deck.kp + "resistances.lumped", 0.0);
  EXPECT_EQ(m.PresetTag(), "macropore");
  ExpectAgrees(deck, m, "");
  // The last term, and knudsen without macropore, are refused; nothing changes.
  EXPECT_THROW(deck.sim.SetParameterText(deck.kp + "resistances.macropore", "false"), std::invalid_argument);
  EXPECT_TRUE(m.Enabled("macropore"));
  deck.sim.SetParameterText(deck.kp + "resistances.knudsen", "true");
  EXPECT_EQ(m.PresetTag(), "macropore_knudsen");
  EXPECT_TRUE(deck.Has("pore_radius"));
  // Knudsen alone: molecular off hides its leaf; the last of the two is refused off.
  deck.sim.SetParameterText(deck.kp + "resistances.molecular", "false");
  EXPECT_FALSE(m.Enabled("molecular"));
  EXPECT_FALSE(deck.Has("molecular_diffusion_coefficient"));
  EXPECT_THROW(deck.sim.SetParameterText(deck.kp + "resistances.knudsen", "false"), std::invalid_argument);
  ExpectAgrees(deck, m, "");
  deck.sim.SetParameterText(deck.kp + "resistances.molecular", "true");
  deck.sim.SetParameterText(deck.kp + "resistances.film", "true");
  deck.sim.SetParameterText(deck.kp + "resistances.macropore", "false");  // clears knudsen
  EXPECT_FALSE(m.Enabled("knudsen"));
  EXPECT_FALSE(deck.Has("resistances.knudsen"));
  EXPECT_FALSE(deck.Has("pore_radius"));
  ExpectAgrees(deck, m, "");

  // The mode by name or number, through SetCapacityFactorMode; its parameters follow.
  deck.sim.SetParameterText(deck.kp + "capacity_factor_mode", "secant");
  EXPECT_EQ(m.CapacityMode(), CapacityFactorMode::SECANT);
  EXPECT_TRUE(deck.Has("concentration_floor"));
  m.capacity_factor_ = 0.0;
  deck.sim.SetParameter(deck.kp + "capacity_factor_mode", 3.0);
  EXPECT_EQ(m.CapacityMode(), CapacityFactorMode::FIXED);
  EXPECT_EQ(m.capacity_factor_.Get(), 1.0) << "the cascade's rule: fixed with K <= 0 sets K = 1";
  ExpectAgrees(deck, m, "");
  // Without a gas-side term there is no mode leaf.
  deck.sim.SetParameterText(deck.kp + "resistances.lumped", "true");
  deck.sim.SetParameterText(deck.kp + "resistances.film", "false");
  EXPECT_FALSE(deck.Has("capacity_factor_mode"));
  EXPECT_FALSE(deck.Has("capacity_factor"));
  ExpectAgrees(deck, m, "");

  // Parameters keep their values across a switch off and on, by path too.
  deck.sim.SetParameterText(deck.kp + "resistances.micropore", "true");
  deck.sim.SetParameterText(deck.kp + "crystal_radius", "2 [um]");
  deck.sim.SetParameterText(deck.kp + "resistances.micropore", "false");
  EXPECT_FALSE(deck.Has("crystal_radius"));
  deck.sim.SetParameterText(deck.kp + "resistances.micropore", "true");
  EXPECT_DOUBLE_EQ(deck.sim.GetParameter(deck.kp + "crystal_radius"), 2e-6);
}

TEST(ParameterRegistry, ResistanceSwitchesOfTheSlowSite) {
  Deck deck("dual_rate_ldf", "macropore");
  ASSERT_NE(deck.Slow(), nullptr);
  ExpectAgrees(deck, deck.Fast(), "");
  ExpectAgrees(deck, *deck.Slow(), "slow_");
  EXPECT_TRUE(deck.Has("slow_tortuosity"));
  deck.sim.SetParameterText(deck.kp + "slow_resistances.micropore", "true");
  EXPECT_TRUE(deck.Slow()->Enabled("micropore"));
  EXPECT_FALSE(deck.Fast().Enabled("micropore"));
  EXPECT_TRUE(deck.Has("slow_crystal_radius"));
  EXPECT_FALSE(deck.Has("crystal_radius"));
  deck.sim.SetParameterText(deck.kp + "slow_capacity_factor_mode", "slope");
  EXPECT_EQ(deck.Slow()->CapacityMode(), CapacityFactorMode::SLOPE);
  EXPECT_EQ(deck.Fast().CapacityMode(), CapacityFactorMode::NONE);
  ExpectAgrees(deck, deck.Fast(), "");
  ExpectAgrees(deck, *deck.Slow(), "slow_");
  // A rate model without a slow site has no slow leaves.
  std::string err;
  ASSERT_TRUE(SetMassTransferKind(deck.sim.Data(), deck.a, kCO2, MassTransferModelType::LDF,
                                  MassTransferCoefficientType::SERIES_RESISTANCE, err))
      << err;
  deck.sim.RefreshParameters();
  for (const auto& p : deck.sim.Parameters().All()) EXPECT_EQ(p.Path().find("].slow_"), std::string::npos) << p.Path();
}

TEST(ParameterRegistry, FormerArrheniusLeafIsTheLumpedK) {
  Deck deck("dual_rate_ldf", "arrhenius");
  EXPECT_EQ(deck.Fast().PresetTag(), "arrhenius");
  deck.sim.SetParameter(deck.kp + "pre_exp_factor", 42.0);
  EXPECT_EQ(deck.Fast().k_.Get(), 42.0);
  deck.sim.SetParameter(deck.kp + "slow_pre_exp_factor", 4.0);
  EXPECT_EQ(deck.Slow()->k_.Get(), 4.0);
  EXPECT_EQ(deck.sim.GetParameter(deck.kp + "pre_exp_factor"), deck.sim.GetParameter(deck.kp + "k"));
  // The name alias works too, and no entry is listed twice.
  EXPECT_EQ(deck.sim.Parameters().Find("adsorbent[13X].kinetics[" + kCO2 + "].pre_exp_factor"),
            deck.sim.Parameters().Find(deck.kp + "k"));
  size_t ks = 0;
  for (const auto& p : deck.sim.Parameters().All()) ks += p.Path() == deck.kp + "k";
  EXPECT_EQ(ks, 1u);
  // Without the lumped term there is nothing to resolve to.
  deck.sim.SetParameterText(deck.kp + "resistances.macropore", "true");
  deck.sim.SetParameterText(deck.kp + "resistances.lumped", "false");
  EXPECT_EQ(deck.sim.Parameters().Find(deck.kp + "pre_exp_factor"), nullptr);
}

TEST(DeckView, KineticsCarryTheSwitchesPresetAndEquations) {
  Deck deck("dual_rate_ldf", "haynes");
  deck.sim.SetParameterText(deck.kp + "slow_resistances.micropore", "true");
  rapidjson::Document view;
  view.Parse(DeckView(deck.sim.Data()).c_str());
  ASSERT_FALSE(view.HasParseError());
  const auto& comps = view["adsorbents"][0]["components"];
  ASSERT_EQ(comps.Size(), 1u);
  const auto& c = comps[0];
  EXPECT_STREQ(c["rate_coefficient"].GetString(), "series_resistance");
  EXPECT_STREQ(c["preset"].GetString(), "haynes");
  EXPECT_TRUE(c["resistances"]["film"].GetBool());
  EXPECT_TRUE(c["resistances"]["macropore"].GetBool());
  EXPECT_TRUE(c["resistances"]["molecular"].GetBool());
  EXPECT_FALSE(c["resistances"]["knudsen"].GetBool());
  EXPECT_FALSE(c["resistances"]["lumped"].GetBool());
  EXPECT_STREQ(c["capacity_factor_mode"].GetString(), "none");
  // In no layer: the packed (sphere) form with the default film correlation;
  // no slab form.
  const FilmCorrelation wf = FilmCorrelation::WakaoFunazkri;
  EXPECT_EQ(std::string(c["equation"].GetString()), MassTransferEquationLatex(deck.Fast(), DiffusionShape::Sphere, &wf));
  EXPECT_NE(std::string(c["equation"].GetString()).find("Sh &= 2 + 1.1"), std::string::npos);
  EXPECT_TRUE(c["equation_slab"].IsNull());
  EXPECT_STREQ(c["slow_preset"].GetString(), "custom");
  EXPECT_TRUE(c["slow_resistances"]["micropore"].GetBool());
  EXPECT_EQ(std::string(c["slow_equation"].GetString()), MassTransferEquationLatex(*deck.Slow(), DiffusionShape::Sphere, &wf));
  EXPECT_NE(std::string(c["slow_equation"].GetString()).find("r_c^{2}"), std::string::npos);

  // Recomputed after an edit: a single-site model has no slow fields.
  std::string err;
  ASSERT_TRUE(SetMassTransferKind(deck.sim.Data(), deck.a, kCO2, MassTransferModelType::LDF,
                                  MassTransferCoefficientType::MACROPORE_CAPACITY, err))
      << err;
  view.Parse(DeckView(deck.sim.Data()).c_str());
  const auto& c2 = view["adsorbents"][0]["components"][0];
  EXPECT_STREQ(c2["preset"].GetString(), "macropore_capacity");
  EXPECT_STREQ(c2["capacity_factor_mode"].GetString(), "secant");
  EXPECT_EQ(std::string(c2["equation"].GetString()),
            R"(\begin{aligned}\frac{1}{k_i} &= \frac{K\,\tau\,r_p^{2}}{15\,\varepsilon_p\,D_{m,i}} \\ K &= \max\!\left(\frac{\rho_p\,q_i^{*}}{c_i},\,1\right)\end{aligned})");
  EXPECT_TRUE(c2["slow_equation"].IsNull());
  EXPECT_TRUE(c2["slow_resistances"].IsNull());
}

TEST(ModelKindManifest, CarriesTheResistanceDescriptor) {
  rapidjson::Document m;
  m.Parse(ModelKindManifest().c_str());
  ASSERT_FALSE(m.HasParseError());
  const auto& kinds = m["rate_coefficient"];
  ASSERT_EQ(kinds.Size(), KindDescriptors(ModelFamily::RateCoefficient).size());
  for (rapidjson::SizeType i = 0; i + 1 < kinds.Size(); ++i) EXPECT_TRUE(kinds[i]["legacy"].GetBool());
  const auto& g = kinds[kinds.Size() - 1];
  EXPECT_STREQ(g["tag"].GetString(), "series_resistance");
  EXPECT_FALSE(g["legacy"].GetBool());
  ASSERT_EQ(g["resistances"].Size(), 6u);
  EXPECT_STREQ(g["resistances"][2]["tag"].GetString(), "molecular");
  EXPECT_STREQ(g["resistances"][2]["parent"].GetString(), "macropore");
  EXPECT_STREQ(g["resistances"][3]["tag"].GetString(), "knudsen");
  EXPECT_STREQ(g["resistances"][3]["parent"].GetString(), "macropore");
  EXPECT_TRUE(g["resistances"][1]["needs_child"].GetBool());
  EXPECT_FALSE(g["resistances"][0]["needs_child"].GetBool());
  EXPECT_TRUE(g["resistances"][0]["parent"].IsNull());
  EXPECT_TRUE(g["resistances"][0]["scaled_by_capacity"].GetBool());
  EXPECT_EQ(std::string(g["resistances"][0]["fragments"]["sphere"]["term"].GetString()),
            ResistanceKind().resistances[0].sphere.term);
  EXPECT_EQ(g["capacity_modes"].Size(), 4u);
  ASSERT_EQ(g["presets"].Size(), 6u);
  EXPECT_STREQ(g["presets"][5]["tag"].GetString(), "macropore_capacity");
  EXPECT_STREQ(g["presets"][5]["capacity_mode"].GetString(), "secant");
  EXPECT_EQ(g["presets"][4]["values"]["film_resistance_factor"].GetDouble(), 1.0);
  // Every other family is listed, without the resistance fields.
  EXPECT_TRUE(m.HasMember("dispersion"));
  EXPECT_FALSE(m["dispersion"][0].HasMember("resistances"));
}

}  // namespace adsim
