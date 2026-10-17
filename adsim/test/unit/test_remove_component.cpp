/*
 * Skarstrom - pressure and temperature swing adsorption simulator
 * Copyright (c) 2024-2026 Sam Affleck and contributors
 *
 * Licensed under the terms in the LICENSE file at the repository root.
 */

// Removing a fluid component must strip it from EVERY library adsorbent (not only
// those currently in a layer), and AdsorptionSource::Initialise must purge any
// per-component entry for a component the fluid lacks. Competitive isotherms
// iterate all of params_ and look each key up in the mole-fraction map, so a stale
// key surfaces as `unordered_map::at: key not found` on the first residual.

#include <sstream>

#include "../fixtures/user_data_test.h"

namespace adsim {

  namespace {
    const std::string kCO2 = "CarbonDioxide";
    const std::string kN2 = "Nitrogen";
    const std::string kO2 = "Oxygen";

    // Serialise an adsorbent exactly as the library/deck JSON does.
    std::string ToJson(const Adsorbent& adsorbent) {
      std::ostringstream os;
      {
        cereal::JSONOutputArchive archive(os);
        archive(adsorbent);
      }
      return os.str();
    }

    // Emulates LoadAdsorbentFromJSON (an imported library adsorbent).
    void FromJson(const std::string& json, Adsorbent& adsorbent) {
      std::istringstream is(json);
      cereal::JSONInputArchive archive(is);
      archive(adsorbent);
    }

    bool Mentions(const std::string& json, const std::string& component) {
      return json.find("\"" + component + "\"") != std::string::npos;
    }
  } // namespace

  class RemoveComponentTest : public UserDataTest {
  public:
    LayerID layer_id{};
    AdsorbentID ads_a{};

    void BuildBed(const std::vector<std::string>& components) {
      AddComponents(components);
      auto ids = AddLayer("Packing");
      layer_id = ids.first;
      ads_a = ids.second;

      auto& layer = GetLayer(layer_id);
      layer.domain_.n_cells = 10;
      std::get<CylinderBed>(data.bed.shape_).diameter = 0.1;
      AddConstantInternalHeatTransferCoefficientModel(0.0, layer_id);
      AddConstantExternalHeatTransferCoefficientModel(0.0, layer_id);

      auto& a = GetAdsorbent(ads_a);
      a.solid_density_ = 1000.0;
      a.heat_capacity_ = 900.0;
    }

    // Adsorbent with an Extended Langmuir isotherm (+ LDF + HOA) on `components`.
    void InstallExtendedLangmuir(AdsorbentID id, const std::vector<std::string>& components) {
      auto iso = std::make_unique<ExtendedLangmuirIsotherm>();
      for (const auto& c : components) {
        iso->params_[c].m1 = 3.0;
        iso->params_[c].b1 = 5.0e-9;
        iso->params_[c].q1 = 18000.0;
      }
      AddIsotherm(components, id, std::move(iso));
      for (const auto& c : components) {
        AddLinearDrivingForceModel(id, c, 1.0);
        SetHOA(id, c, 20000.0);
      }
    }

    void RunClosedHold(const std::unordered_map<std::string, double>& y) {
      AddInitialConditions(1.0e5, 300.0, 300.0, y, 0.0);
      AddStep("Hold", 1.0);
      AddClosedBC(BoundaryLocation::BOTTOM, "Hold");
      AddClosedBC(BoundaryLocation::TOP, "Hold");
      ConfigureCycle(1);

      std::string error;
      ASSERT_TRUE(data.PerformChecks(error)) << error;
      try {
        CoupledSolver solver{data};
        EXPECT_GE(solver.IntegrateCycle(data.cycle, false), 0) << "Solver did not complete";
      } catch (const std::exception& e) {
        FAIL() << "Solver threw: " << e.what();
      }
    }
  };

  // Adsorbent B is in the library but NOT in a layer when N2 is removed; it is
  // assigned to the layer afterwards. Before the fix B kept N2 in its Extended
  // Langmuir params_ and the run threw on yi_map.at("Nitrogen").
  TEST_F(RemoveComponentTest, UnassignedLibraryAdsorbentIsPurgedAndRuns) {
    BuildBed({kCO2, kN2, kO2});
    const AdsorbentID ads_b = data.adsorbent_directory.Add("B");
    GetAdsorbent(ads_b).solid_density_ = 1000.0;
    InstallExtendedLangmuir(ads_b, {kCO2, kN2});

    data.fluid.RemoveComponent(kN2);
    data.bed.RemoveComponent(kN2, data.adsorbent_directory);

    // Strip from the library adsorbent immediately (not deferred to Initialise).
    EXPECT_FALSE(Mentions(ToJson(GetAdsorbent(ads_b)), kN2));

    auto& layer = GetLayer(layer_id);
    layer.composition_.clear();
    layer.composition_[ads_b] = 1.0;

    RunClosedHold({{kCO2, 0.2}, {kO2, 0.8}});
  }

  // RemoveComponent must erase from every per-component map and never insert:
  // the old code's GetIsotherm(component) lookup inserted isotherm_map_[N2] = INERT.
  TEST_F(RemoveComponentTest, LeavesNoEntryForRemovedComponent) {
    BuildBed({kCO2, kN2, kO2});
    const AdsorbentID ads_b = data.adsorbent_directory.Add("B");
    InstallExtendedLangmuir(ads_b, {kCO2, kN2});
    // Layer adsorbent A: Inert everywhere, seeded by Initialise (isotherm map,
    // Inert params, LDF, HOA for every component).
    GetAdsorbent(ads_a).adsorption_source_.Initialise(data.fluid);
    ASSERT_TRUE(Mentions(ToJson(GetAdsorbent(ads_a)), kN2));

    data.fluid.RemoveComponent(kN2);
    data.bed.RemoveComponent(kN2, data.adsorbent_directory);

    for (const auto& [id, ads] : data.adsorbent_directory.GetAll()) {
      const auto& src = ads.adsorption_source_;
      EXPECT_EQ(src.heat_of_adsorption_.count(kN2), 0u) << ads.name_;
      EXPECT_EQ(src.mass_transfer_models_.count(kN2), 0u) << ads.name_;
      // Covers isotherm_map_ and every isotherm's params_ (incl. Inert).
      EXPECT_FALSE(Mentions(ToJson(ads), kN2)) << ads.name_;
      EXPECT_TRUE(Mentions(ToJson(ads), kCO2)) << ads.name_;
    }

    // Removing a component nobody has must not insert it either.
    data.bed.RemoveComponent("Argon", data.adsorbent_directory);
    for (const auto& [id, ads] : data.adsorbent_directory.GetAll()) {
      EXPECT_FALSE(Mentions(ToJson(ads), "Argon")) << ads.name_;
    }
  }

  // An imported adsorbent (LoadAdsorbentFromJSON path) carrying params for a
  // component the fluid lacks must initialise and run; Initialise purges it.
  TEST_F(RemoveComponentTest, ImportedAdsorbentWithExtraComponentRuns) {
    BuildBed({kCO2, kO2});

    Adsorbent library("Imported");
    library.solid_density_ = 1000.0;
    library.heat_capacity_ = 900.0;
    {
      auto iso = std::make_unique<ExtendedLangmuirIsotherm>();
      for (const auto& c : {kCO2, kN2, kO2}) {
        iso->params_[c].m1 = 3.0;
        iso->params_[c].b1 = 5.0e-9;
        iso->params_[c].q1 = 18000.0;
      }
      library.adsorption_source_.AddIsotherm({kCO2, kN2, kO2}, std::move(iso));
      for (const auto& c : {kCO2, kN2, kO2}) {
        auto ldf = std::make_unique<LinearDrivingForce>();
        auto mtc = std::make_unique<MassTransferCoefficientConstant>();
        mtc->constant_mtc_ = 1.0;
        ldf->mtc_model_ = std::move(mtc);
        library.adsorption_source_.mass_transfer_models_[c] = std::move(ldf);
        library.adsorption_source_.heat_of_adsorption_[c] = 20000.0;
      }
    }
    FromJson(ToJson(library), GetAdsorbent(ads_a));
    ASSERT_TRUE(Mentions(ToJson(GetAdsorbent(ads_a)), kN2));

    GetAdsorbent(ads_a).adsorption_source_.Initialise(data.fluid);
    const auto& src = GetAdsorbent(ads_a).adsorption_source_;
    EXPECT_EQ(src.heat_of_adsorption_.count(kN2), 0u);
    EXPECT_EQ(src.mass_transfer_models_.count(kN2), 0u);
    EXPECT_FALSE(Mentions(ToJson(GetAdsorbent(ads_a)), kN2));
    EXPECT_EQ(src.heat_of_adsorption_.count(kCO2), 1u);

    RunClosedHold({{kCO2, 0.2}, {kO2, 0.8}});
  }

  // Same, but without an explicit Initialise: the solver's own initialisation
  // must purge before the first residual.
  TEST_F(RemoveComponentTest, SolverInitialisePurgesExtraComponent) {
    BuildBed({kCO2, kO2});
    InstallExtendedLangmuir(ads_a, {kCO2, kN2, kO2});
    RunClosedHold({{kCO2, 0.2}, {kO2, 0.8}});
    EXPECT_FALSE(Mentions(ToJson(GetAdsorbent(ads_a)), kN2));
  }

} // namespace adsim
