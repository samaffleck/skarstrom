// STL includes
#include <filesystem>

// ads includes
#include "adsim/equations/isotherm.h"
#include "adsim/objects/porous_media.h"
#include "adsim/objects/reactor.h"
#include "adsim/objects/domain.h"
#include "adsim/objects/fluid.h"

// Gtest includes
#include <gtest/gtest.h>


namespace adsim {
  
  const std::filesystem::path db_path = std::filesystem::current_path() / ".." / "data" / "adsorbents";

  // Blobs are saved as <file_stem>.json, where file_stem is the entry's key in
  // data/adsorbents/library.json ([a-z0-9_]+), not the display name:
  // spaces and parentheses in data file names break the asset-copy build steps
  // of front ends (see cmake/check_data_names.cmake).
  void SaveAdsorbent(const Adsorbent& adsorbent, const std::string& file_stem) {
    std::filesystem::create_directories(db_path);
    
    std::string json_filename = file_stem + ".json";
    std::filesystem::path json_path = db_path / json_filename;

    SaveAdsorbentToJSON(json_path.string(), adsorbent);
  }

  void LoadAdsorbent(const std::string& file_stem, Adsorbent& adsorbent) {
    std::string json_filename = file_stem + ".json";
    std::filesystem::path json_path = db_path / json_filename;
    
    LoadAdsorbentFromJSON(json_path.string(), adsorbent);
  }


  TEST(AdsorbentTest, MgMOF74) {
    // Ward and Pini, 2024
    // https://doi.org/10.1021/acs.iecr.4c01345
    // Table 1: CO2/N2 DSL isotherm parameters (concentration form); Sect. 2.1
    // states the extended (competitive) DSL, so the extended kind is used.
    const std::string CO2 = "CarbonDioxide";
    const std::string N2 = "Nitrogen";

    Adsorbent adsorbent;
    std::string adsorbent_name = "Mg-MOF-74";
    adsorbent.name_ = adsorbent_name;
    constexpr double eb = 0.52; // Now lives on MixedAdsorbentLayer
    adsorbent.solid_density_ = 437.0 / (1.0 - eb);
    adsorbent.heat_capacity_ = 900;
    adsorbent.thermal_conductivity_ = 0.15;

    auto dsl = std::make_unique<ExtendedDualSiteLangmuirIsotherm_Concentration>(CO2);
    dsl->params_[CO2].m1 = 6.80;    // mol/kg  (qsb,CO2)
    dsl->params_[CO2].b1 = 1.81e-7; // m3/mol  (b0,CO2)
    dsl->params_[CO2].q1 = 39300;   // J/mol   (-ΔUb,CO2)
    dsl->params_[CO2].m2 = 9.90;    // mol/kg  (qsd,CO2)
    dsl->params_[CO2].b2 = 1.06e-6; // m3/mol  (d0,CO2)
    dsl->params_[CO2].q2 = 21200;   // J/mol   (-ΔUd,CO2)
    dsl->params_[N2].m1 = 14.0;
    dsl->params_[N2].b1 = 3.45e-6;
    dsl->params_[N2].q1 = 15500;
    dsl->params_[N2].m2 = 0.0;
    dsl->params_[N2].b2 = 0.0;
    dsl->params_[N2].q2 = 0.0;
    adsorbent.adsorption_source_.AddIsotherm({CO2, N2}, std::move(dsl));
    adsorbent.adsorption_source_.heat_of_adsorption_[CO2] = 30.85e3;  // J/mol
    adsorbent.adsorption_source_.heat_of_adsorption_[N2] = 15.93e3;   // J/mol

    SaveAdsorbent(adsorbent, "mg_mof_74");
    Adsorbent adsorbent_cpy;
    LoadAdsorbent("mg_mof_74", adsorbent_cpy);
    EXPECT_EQ(adsorbent.name_, adsorbent_cpy.name_);
    EXPECT_DOUBLE_EQ(adsorbent.solid_density_, adsorbent_cpy.solid_density_);
    EXPECT_DOUBLE_EQ(adsorbent.heat_capacity_, adsorbent_cpy.heat_capacity_);
    EXPECT_DOUBLE_EQ(adsorbent.thermal_conductivity_, adsorbent_cpy.thermal_conductivity_);

    auto* iso = dynamic_cast<ExtendedDualSiteLangmuirIsotherm_Concentration*>(
      adsorbent_cpy.adsorption_source_.GetIsotherm(CO2).get());
    ASSERT_NE(iso, nullptr);
    EXPECT_DOUBLE_EQ(iso->params_.at(CO2).m1, 6.80);
    EXPECT_DOUBLE_EQ(iso->params_.at(CO2).q1, 39300);
    EXPECT_DOUBLE_EQ(iso->params_.at(N2).m1, 14.0);
    EXPECT_DOUBLE_EQ(iso->params_.at(N2).q1, 15500);
  }

  TEST(AdsorbentTest, ActivatedCarbon) {
    // Ward and Pini, 2024
    // https://doi.org/10.1021/acs.iecr.4c01345
    // Table 1: CO2/N2 DSL isotherm parameters (concentration form); Sect. 2.1
    // states the extended (competitive) DSL, so the extended kind is used.
    const std::string CO2 = "CarbonDioxide";
    const std::string N2 = "Nitrogen";

    Adsorbent adsorbent;
    std::string adsorbent_name = "Activated Carbon";
    adsorbent.name_ = adsorbent_name;
    constexpr double eb = 0.46; // Now lives on MixedAdsorbentLayer
    adsorbent.solid_density_ = 435.0 / (1.0 - eb);
    adsorbent.heat_capacity_ = 1040;
    adsorbent.thermal_conductivity_ = 0.15;

    auto dsl = std::make_unique<ExtendedDualSiteLangmuirIsotherm_Concentration>(CO2);
    dsl->params_[CO2].m1 = 0.59;    // mol/kg  (qsb,CO2)
    dsl->params_[CO2].b1 = 9.40e-6; // m3/mol  (b0,CO2)
    dsl->params_[CO2].q1 = 25600;   // J/mol   (-ΔUb,CO2)
    dsl->params_[CO2].m2 = 7.51;    // mol/kg  (qsd,CO2)
    dsl->params_[CO2].b2 = 1.04e-5; // m3/mol  (d0,CO2)
    dsl->params_[CO2].q2 = 17600;   // J/mol   (-ΔUd,CO2)
    dsl->params_[N2].m1 = 0.16;
    dsl->params_[N2].b1 = 1.81e-3;
    dsl->params_[N2].q1 = 8700;
    dsl->params_[N2].m2 = 41.3;
    dsl->params_[N2].b2 = 1.72e-12;
    dsl->params_[N2].q2 = 44900;
    adsorbent.adsorption_source_.AddIsotherm({CO2, N2}, std::move(dsl));

    SaveAdsorbent(adsorbent, "activated_carbon");
    Adsorbent adsorbent_cpy;
    LoadAdsorbent("activated_carbon", adsorbent_cpy);
    EXPECT_EQ(adsorbent.name_, adsorbent_cpy.name_);
    EXPECT_DOUBLE_EQ(adsorbent.solid_density_, adsorbent_cpy.solid_density_);
    EXPECT_DOUBLE_EQ(adsorbent.heat_capacity_, adsorbent_cpy.heat_capacity_);
    EXPECT_DOUBLE_EQ(adsorbent.thermal_conductivity_, adsorbent_cpy.thermal_conductivity_);

    auto* iso = dynamic_cast<ExtendedDualSiteLangmuirIsotherm_Concentration*>(
      adsorbent_cpy.adsorption_source_.GetIsotherm(CO2).get());
    ASSERT_NE(iso, nullptr);
    EXPECT_DOUBLE_EQ(iso->params_.at(CO2).m1, 0.59);
    EXPECT_DOUBLE_EQ(iso->params_.at(CO2).q2, 17600);
    EXPECT_DOUBLE_EQ(iso->params_.at(N2).m2, 41.3);
    EXPECT_DOUBLE_EQ(iso->params_.at(N2).q2, 44900);
  }

  // ---------- GetAdsorbentMass tests ----------

  TEST(AdsorbentMassTest, SingleLayerCylinder) {
    AdsorbentDirectory dir;
    AdsorbentID id = dir.Add("Zeolite 13X");
    dir.Get(id).solid_density_ = 1200.0;

    Bed bed;
    std::get<CylinderBed>(bed.shape_).diameter = 0.1; // area = pi/4 * 0.01
    LayerID lid = bed.layers_.Add("Layer 1");
    auto& layer = bed.layers_.Get(lid);
    layer.domain_.length = 1.0;
    layer.eb_ = 0.4;
    layer.AddAdsorbent(id, 1.0);

    // mass = area * length * (1 - eb) * solid_density
    const double area = (constants::kPI / 4.0) * 0.1 * 0.1;
    const double expected = area * 1.0 * (1.0 - 0.4) * 1200.0;

    EXPECT_NEAR(GetAdsorbentMass(bed, dir), expected, 1e-10);
  }

  TEST(AdsorbentMassTest, SingleLayerSquare) {
    AdsorbentDirectory dir;
    AdsorbentID id = dir.Add("Adsorbent A");
    dir.Get(id).solid_density_ = 800.0;

    Bed bed;
    bed.shape_ = SquareBed{0.05};
    LayerID lid = bed.layers_.Add("Layer 1");
    auto& layer = bed.layers_.Get(lid);
    layer.domain_.length = 0.5;
    layer.eb_ = 0.35;
    layer.AddAdsorbent(id, 1.0);

    const double area = 0.05 * 0.05;
    const double expected = area * 0.5 * (1.0 - 0.35) * 800.0;

    EXPECT_NEAR(GetAdsorbentMass(bed, dir), expected, 1e-10);
  }

  TEST(AdsorbentMassTest, SingleLayerRectangle) {
    AdsorbentDirectory dir;
    AdsorbentID id = dir.Add("Adsorbent B");
    dir.Get(id).solid_density_ = 1000.0;

    Bed bed;
    bed.shape_ = RectangularBed{0.04, 0.08};
    LayerID lid = bed.layers_.Add("Layer 1");
    auto& layer = bed.layers_.Get(lid);
    layer.domain_.length = 2.0;
    layer.eb_ = 0.5;
    layer.AddAdsorbent(id, 1.0);

    const double area = 0.04 * 0.08;
    const double expected = area * 2.0 * (1.0 - 0.5) * 1000.0;

    EXPECT_NEAR(GetAdsorbentMass(bed, dir), expected, 1e-10);
  }

  TEST(AdsorbentMassTest, TwoLayersDifferentProperties) {
    AdsorbentDirectory dir;
    AdsorbentID id_a = dir.Add("Adsorbent A");
    dir.Get(id_a).solid_density_ = 1200.0;
    AdsorbentID id_b = dir.Add("Adsorbent B");
    dir.Get(id_b).solid_density_ = 900.0;

    Bed bed;
    std::get<CylinderBed>(bed.shape_).diameter = 0.05;
    const double area = (constants::kPI / 4.0) * 0.05 * 0.05;

    // Layer 1: high void fraction, dense adsorbent
    LayerID lid1 = bed.layers_.Add("Bottom");
    auto& layer1 = bed.layers_.Get(lid1);
    layer1.domain_.length = 0.5;
    layer1.eb_ = 0.45;
    layer1.AddAdsorbent(id_a, 1.0);

    // Layer 2: low void fraction, lighter adsorbent
    LayerID lid2 = bed.layers_.Add("Top");
    auto& layer2 = bed.layers_.Get(lid2);
    layer2.domain_.length = 1.5;
    layer2.eb_ = 0.30;
    layer2.AddAdsorbent(id_b, 1.0);

    const double mass1 = area * 0.5 * (1.0 - 0.45) * 1200.0;
    const double mass2 = area * 1.5 * (1.0 - 0.30) * 900.0;

    EXPECT_NEAR(GetAdsorbentMass(bed, dir), mass1 + mass2, 1e-10);
  }

  TEST(AdsorbentMassTest, MixedCompositionLayer) {
    AdsorbentDirectory dir;
    AdsorbentID id_a = dir.Add("Light");
    dir.Get(id_a).solid_density_ = 600.0;
    AdsorbentID id_b = dir.Add("Heavy");
    dir.Get(id_b).solid_density_ = 1400.0;

    Bed bed;
    std::get<CylinderBed>(bed.shape_).diameter = 0.1;
    LayerID lid = bed.layers_.Add("Mixed");
    auto& layer = bed.layers_.Get(lid);
    layer.domain_.length = 1.0;
    layer.eb_ = 0.37;
    layer.AddAdsorbent(id_a, 0.3);
    layer.AddAdsorbent(id_b, 0.7);

    // Composition is by MASS: 0.3 kg of 600 kg/m3 + 0.7 kg of 1400 kg/m3 occupy
    // 0.3/600 + 0.7/1400 = 1e-3 m3, so the mixture particle density is 1000 kg/m3.
    const double area = (constants::kPI / 4.0) * 0.1 * 0.1;
    const double effective_density = 1.0 / (0.3 / 600.0 + 0.7 / 1400.0);
    EXPECT_NEAR(effective_density, 1000.0, 1e-9);
    const double expected = area * 1.0 * (1.0 - 0.37) * effective_density;

    EXPECT_NEAR(GetAdsorbentMass(bed, dir), expected, 1e-10);
  }

  // Two adsorbents of different densities mixed by mass (E-04): the mixture
  // density is the harmonic mean, the particle voidage and conductivity are
  // weighted by solid volume fraction, Cp by mass fraction, and the
  // per-adsorbent bulk densities split the layer bulk density in the mass ratio.
  TEST(MixedAdsorbentPropertiesTest, TwoDensityMassFractionMixture) {
    AdsorbentDirectory dir;
    AdsorbentID id_a = dir.Add("Light");
    auto& a = dir.Get(id_a);
    a.solid_density_ = 600.0;
    a.ep_ = 0.6;
    a.thermal_conductivity_ = 0.2;
    a.heat_capacity_ = 1000.0;
    AdsorbentID id_b = dir.Add("Heavy");
    auto& b = dir.Get(id_b);
    b.solid_density_ = 1400.0;
    b.ep_ = 0.2;
    b.thermal_conductivity_ = 1.0;
    b.heat_capacity_ = 800.0;

    MixedAdsorbentLayer layer("Mixed");
    layer.eb_ = 0.4;
    layer.AddAdsorbent(id_a, 0.3);
    layer.AddAdsorbent(id_b, 0.7);

    // rho_mix = 1 / (0.3/600 + 0.7/1400) = 1000; phi_a = 0.3*1000/600 = 0.5,
    // phi_b = 0.7*1000/1400 = 0.5.
    EXPECT_NEAR(layer.GetEffectiveSolidDensity(dir), 1000.0, 1e-9);
    const auto phi = layer.GetSolidVolumeFractions(dir);
    EXPECT_NEAR(phi.at(id_a), 0.5, 1e-12);
    EXPECT_NEAR(phi.at(id_b), 0.5, 1e-12);
    EXPECT_NEAR(layer.GetEffectiveParticleVoidFraction(dir), 0.5 * 0.6 + 0.5 * 0.2, 1e-12);
    EXPECT_NEAR(layer.GetEffectiveThermalConductivity(dir), 0.5 * 0.2 + 0.5 * 1.0, 1e-12);
    EXPECT_NEAR(layer.GetEffectiveHeatCapacity(dir), 0.3 * 1000.0 + 0.7 * 800.0, 1e-9);
    EXPECT_NEAR(layer.GetTotalVoidFraction(dir), 0.4 + 0.6 * 0.4, 1e-12);

    const double rho_b = layer.GetBulkDensity(dir);
    EXPECT_NEAR(rho_b, 0.6 * 1000.0, 1e-9);
    const auto rho_bk = layer.GetAdsorbentBulkDensities(dir);
    EXPECT_NEAR(rho_bk.at(id_a), 0.6 * 300.0, 1e-9);
    EXPECT_NEAR(rho_bk.at(id_b), 0.6 * 700.0, 1e-9);
    EXPECT_NEAR(rho_bk.at(id_a) + rho_bk.at(id_b), rho_b, 1e-9);
    EXPECT_NEAR(rho_bk.at(id_a) / rho_bk.at(id_b), 0.3 / 0.7, 1e-12);
  }

  // The per-cell decomposition used by the balances must agree with the layer.
  TEST(MixedAdsorbentPropertiesTest, DomainVectorsPerAdsorbentBulkDensity) {
    AdsorbentDirectory dir;
    AdsorbentID id_a = dir.Add("Light");
    dir.Get(id_a).solid_density_ = 600.0;
    AdsorbentID id_b = dir.Add("Heavy");
    dir.Get(id_b).solid_density_ = 1400.0;

    Bed bed;
    std::get<CylinderBed>(bed.shape_).diameter = 0.1;
    LayerID lid = bed.layers_.Add("Mixed");
    auto& layer = bed.layers_.Get(lid);
    layer.domain_.length = 1.0;
    layer.domain_.n_cells = 5;
    layer.eb_ = 0.4;
    layer.AddAdsorbent(id_a, 0.3);
    layer.AddAdsorbent(id_b, 0.7);

    Fluid fluid;
    fluid.AddComponent("N2");
    bed.Initialise(dir, fluid);
    DomainVectors dv;
    dv.Initialise(bed, dir, fluid);

    for (int i = layer.domain_.start_index; i <= layer.domain_.end_index; ++i) {
      EXPECT_NEAR(dv.bulk_density[i], 600.0, 1e-9);
      ASSERT_EQ(dv.cell_adsorbents[i].size(), 2u);
      double sum = 0.0;
      for (const auto& ca : dv.cell_adsorbents[i]) {
        const double expected = (ca.id == id_a) ? 0.6 * 300.0 : 0.6 * 700.0;
        EXPECT_NEAR(ca.bulk_density, expected, 1e-9);
        sum += ca.bulk_density;
      }
      EXPECT_NEAR(sum, dv.bulk_density[i], 1e-9);
    }
  }

  // A single adsorbent must reproduce its own properties exactly (bit-for-bit),
  // and an un-normalised composition keeps scaling the properties linearly with
  // sum_k w_k, as before the mass-fraction mixing rules.
  TEST(MixedAdsorbentPropertiesTest, SingleAdsorbentExactAndUnnormalisedScaling) {
    AdsorbentDirectory dir;
    AdsorbentID id = dir.Add("Only");
    auto& ads = dir.Get(id);
    ads.solid_density_ = 1123.456789;
    ads.ep_ = 0.3456789;
    ads.thermal_conductivity_ = 0.1234567;
    ads.heat_capacity_ = 876.54321;

    MixedAdsorbentLayer layer("Single");
    layer.eb_ = 0.3711;
    layer.AddAdsorbent(id, 1.0);

    const double rho = ads.solid_density_;
    const double ep = ads.ep_;
    const double k = ads.thermal_conductivity_;
    EXPECT_EQ(layer.GetEffectiveSolidDensity(dir), rho);
    EXPECT_EQ(layer.GetEffectiveParticleVoidFraction(dir), ep);
    EXPECT_EQ(layer.GetEffectiveThermalConductivity(dir), k);
    EXPECT_EQ(layer.GetBulkDensity(dir), (1 - 0.3711) * rho);
    EXPECT_EQ(layer.GetAdsorbentBulkDensities(dir).at(id), (1.0 - 0.3711) * 1.0 * rho);

    layer.AddAdsorbent(id, 0.5);
    EXPECT_DOUBLE_EQ(layer.GetEffectiveSolidDensity(dir), 0.5 * rho);
    EXPECT_DOUBLE_EQ(layer.GetEffectiveParticleVoidFraction(dir), 0.5 * ep);
    EXPECT_DOUBLE_EQ(layer.GetAdsorbentBulkDensities(dir).at(id), layer.GetBulkDensity(dir));
  }

  TEST(AdsorbentMassTest, EmptyBed) {
    AdsorbentDirectory dir;
    Bed bed;
    EXPECT_DOUBLE_EQ(GetAdsorbentMass(bed, dir), 0.0);
  }

  // ---------- Adsorbents from literature ----------
  // Reference: Santos et al., 2010 (Table 2 and Table 4)
  // Bin Wu, Xiangping Zhang, Yajing Xu, Di Bao, Suojiang Zhang,
  // Assessment of the energy consumption of the biogas upgrading process with pressure swing adsorption using novel adsorbents,
  // Journal of Cleaner Production,
  // Volume 101,
  // 2015,
  // Pages 251-261,
  // ISSN 0959-6526,
  // https://doi.org/10.1016/j.jclepro.2015.03.082.

  TEST(AdsorbentTest, Zeolite13X_CO2CH4_ExtendedSingleSiteLangmuir) {
    const std::string CO2 = "CarbonDioxide";
    const std::string CH4 = "Methane";

    Adsorbent adsorbent;
    adsorbent.name_ = "Zeolite 13X";
    adsorbent.solid_density_ = 1324.0;   // kg/m3 (pellet density)
    adsorbent.heat_capacity_ = 900.0;    // J/(kg.K)
    adsorbent.thermal_conductivity_ = 0.15;

    auto ssl = std::make_unique<ExtendedLangmuirIsotherm>(CO2);
    ssl->params_[CO2].m1 = 3.75e-3 * 1e3;   // kmol/kg -> mol/kg
    ssl->params_[CO2].b1 = 3.96e0 * 1e-5;   // 1/bar -> 1/Pa
    ssl->params_[CO2].q1 = 31.164 * 1e3;    // kJ/mol -> J/mol
    ssl->params_[CH4].m1 = 1.42e-2 * 1e3;   // kmol/kg -> mol/kg
    ssl->params_[CH4].b1 = 1.32e-2 * 1e-5;  // 1/bar -> 1/Pa
    ssl->params_[CH4].q1 = 20.856 * 1e3;    // kJ/mol -> J/mol
    adsorbent.adsorption_source_.AddIsotherm({CO2, CH4}, std::move(ssl));
    adsorbent.adsorption_source_.heat_of_adsorption_[CO2] = 31.164 * 1e3;
    adsorbent.adsorption_source_.heat_of_adsorption_[CH4] = 20.856 * 1e3;

    // Not a library entry (the unsourced legacy 13X blob was removed from
    // data/adsorbents/), so the round trip goes through a temporary file.
    const auto json_path = std::filesystem::temp_directory_path() / "adsim_zeolite_13x.json";
    SaveAdsorbentToJSON(json_path.string(), adsorbent);
    Adsorbent adsorbent_cpy;
    LoadAdsorbentFromJSON(json_path.string(), adsorbent_cpy);
    std::filesystem::remove(json_path);

    EXPECT_EQ(adsorbent.name_, adsorbent_cpy.name_);
    EXPECT_DOUBLE_EQ(adsorbent.solid_density_, adsorbent_cpy.solid_density_);
    EXPECT_DOUBLE_EQ(adsorbent.heat_capacity_, adsorbent_cpy.heat_capacity_);

    auto* iso = dynamic_cast<ExtendedLangmuirIsotherm*>(
      adsorbent_cpy.adsorption_source_.GetIsotherm(CO2).get());
    ASSERT_NE(iso, nullptr);
    EXPECT_DOUBLE_EQ(iso->params_.at(CO2).m1, 3.75e-3 * 1e3);
    EXPECT_DOUBLE_EQ(iso->params_.at(CO2).b1, 3.96e0 * 1e-5);
    EXPECT_DOUBLE_EQ(iso->params_.at(CO2).q1, 31.164 * 1e3);
    EXPECT_DOUBLE_EQ(iso->params_.at(CH4).m1, 1.42e-2 * 1e3);
    EXPECT_DOUBLE_EQ(iso->params_.at(CH4).b1, 1.32e-2 * 1e-5);
    EXPECT_DOUBLE_EQ(iso->params_.at(CH4).q1, 20.856 * 1e3);
    EXPECT_TRUE(adsorbent_cpy.adsorption_source_.mass_transfer_models_.empty());
  }

  TEST(AdsorbentTest, CMS3K_CO2CH4_ExtendedSingleSiteLangmuir) {
    const std::string CO2 = "CarbonDioxide";
    const std::string CH4 = "Methane";

    Adsorbent adsorbent;
    adsorbent.name_ = "CMS-3K";
    adsorbent.solid_density_ = 1060.0;   // kg/m3 (pellet density)
    adsorbent.heat_capacity_ = 880.0;    // J/(kg.K)
    adsorbent.thermal_conductivity_ = 0.15;

    auto ssl = std::make_unique<ExtendedLangmuirIsotherm>(CO2);
    ssl->params_[CO2].m1 = 2.47e-3 * 1e3;   // kmol/kg -> mol/kg
    ssl->params_[CO2].b1 = 1.13e0 * 1e-5;   // 1/bar -> 1/Pa
    ssl->params_[CO2].q1 = 38.947 * 1e3;    // kJ/mol -> J/mol
    ssl->params_[CH4].m1 = 2.13e-3 * 1e3;   // kmol/kg -> mol/kg
    ssl->params_[CH4].b1 = 3.10e-1 * 1e-5;  // 1/bar -> 1/Pa
    ssl->params_[CH4].q1 = 33.674 * 1e3;    // kJ/mol -> J/mol
    adsorbent.adsorption_source_.AddIsotherm({CO2, CH4}, std::move(ssl));
    adsorbent.adsorption_source_.heat_of_adsorption_[CO2] = 38.947 * 1e3;
    adsorbent.adsorption_source_.heat_of_adsorption_[CH4] = 33.674 * 1e3;

    SaveAdsorbent(adsorbent, "cms_3k");
    Adsorbent adsorbent_cpy;
    LoadAdsorbent("cms_3k", adsorbent_cpy);

    EXPECT_EQ(adsorbent.name_, adsorbent_cpy.name_);
    EXPECT_DOUBLE_EQ(adsorbent.solid_density_, adsorbent_cpy.solid_density_);
    EXPECT_DOUBLE_EQ(adsorbent.heat_capacity_, adsorbent_cpy.heat_capacity_);

    auto* iso = dynamic_cast<ExtendedLangmuirIsotherm*>(
      adsorbent_cpy.adsorption_source_.GetIsotherm(CO2).get());
    ASSERT_NE(iso, nullptr);
    EXPECT_DOUBLE_EQ(iso->params_.at(CO2).m1, 2.47e-3 * 1e3);
    EXPECT_DOUBLE_EQ(iso->params_.at(CO2).b1, 1.13e0 * 1e-5);
    EXPECT_DOUBLE_EQ(iso->params_.at(CO2).q1, 38.947 * 1e3);
    EXPECT_DOUBLE_EQ(iso->params_.at(CH4).m1, 2.13e-3 * 1e3);
    EXPECT_DOUBLE_EQ(iso->params_.at(CH4).b1, 3.10e-1 * 1e-5);
    EXPECT_DOUBLE_EQ(iso->params_.at(CH4).q1, 33.674 * 1e3);
    EXPECT_TRUE(adsorbent_cpy.adsorption_source_.mass_transfer_models_.empty());
  }

  TEST(AdsorbentTest, MOF508b_CO2CH4_ExtendedSingleSiteLangmuir) {
    const std::string CO2 = "CarbonDioxide";
    const std::string CH4 = "Methane";

    Adsorbent adsorbent;
    adsorbent.name_ = "MOF-508b";
    adsorbent.solid_density_ = 1492.0;   // kg/m3 (pellet density)
    adsorbent.heat_capacity_ = 760.0;    // J/(kg.K)
    adsorbent.thermal_conductivity_ = 0.15;

    auto ssl = std::make_unique<ExtendedLangmuirIsotherm>(CO2);
    ssl->params_[CO2].m1 = 1.11e-2 * 1e3;   // kmol/kg -> mol/kg
    ssl->params_[CO2].b1 = 0.12e0 * 1e-5;   // 1/bar -> 1/Pa
    ssl->params_[CO2].q1 = 18.80 * 1e3;     // kJ/mol -> J/mol
    ssl->params_[CH4].m1 = 5.25e-3 * 1e3;   // kmol/kg -> mol/kg
    ssl->params_[CH4].b1 = 1.20e-1 * 1e-5;  // 1/bar -> 1/Pa
    ssl->params_[CH4].q1 = 3.65 * 1e3;      // kJ/mol -> J/mol
    adsorbent.adsorption_source_.AddIsotherm({CO2, CH4}, std::move(ssl));
    adsorbent.adsorption_source_.heat_of_adsorption_[CO2] = 18.80 * 1e3;
    adsorbent.adsorption_source_.heat_of_adsorption_[CH4] = 3.65 * 1e3;

    SaveAdsorbent(adsorbent, "mof_508b");
    Adsorbent adsorbent_cpy;
    LoadAdsorbent("mof_508b", adsorbent_cpy);

    EXPECT_EQ(adsorbent.name_, adsorbent_cpy.name_);
    EXPECT_DOUBLE_EQ(adsorbent.solid_density_, adsorbent_cpy.solid_density_);
    EXPECT_DOUBLE_EQ(adsorbent.heat_capacity_, adsorbent_cpy.heat_capacity_);

    auto* iso = dynamic_cast<ExtendedLangmuirIsotherm*>(
      adsorbent_cpy.adsorption_source_.GetIsotherm(CO2).get());
    ASSERT_NE(iso, nullptr);
    EXPECT_DOUBLE_EQ(iso->params_.at(CO2).m1, 1.11e-2 * 1e3);
    EXPECT_DOUBLE_EQ(iso->params_.at(CO2).b1, 0.12e0 * 1e-5);
    EXPECT_DOUBLE_EQ(iso->params_.at(CO2).q1, 18.80 * 1e3);
    EXPECT_DOUBLE_EQ(iso->params_.at(CH4).m1, 5.25e-3 * 1e3);
    EXPECT_DOUBLE_EQ(iso->params_.at(CH4).b1, 1.20e-1 * 1e-5);
    EXPECT_DOUBLE_EQ(iso->params_.at(CH4).q1, 3.65 * 1e3);
    EXPECT_TRUE(adsorbent_cpy.adsorption_source_.mass_transfer_models_.empty());
  }

  // Reference: Arora, A., Hasan, M.M.F. Flexible oxygen concentrators for medical applications.
  // Sci Rep 11, 14317 (2021). https://doi.org/10.1038/s41598-021-93796-3

  TEST(AdsorbentTest, LiXZeolite) {
    // Table 2 parameters for LiX zeolite (13X base, ~25% Li exchange)
    const std::string O2 = "Oxygen";
    const std::string N2 = "Nitrogen";

    Adsorbent adsorbent;
    adsorbent.name_ = "LiX Zeolite";
    adsorbent.solid_density_ = 1200.0;  // kg/m3
    adsorbent.heat_capacity_ = 1172.0;  // J/(kg·K)
    adsorbent.ep_ = 0.33;  // particle porosity (Table 2)
    adsorbent.thermal_conductivity_ = 0.15;

    auto dsl = std::make_unique<DualSiteLangmuirIsotherm>(O2);
    dsl->params_[O2].m1 = 1.5;       // mol/kg
    dsl->params_[O2].b1 = 1e-10;     // Pa^-1
    dsl->params_[O2].q1 = 21774.6;   // J/mol  (-ΔU1)
    dsl->params_[O2].m2 = 0.6075;    // mol/kg
    dsl->params_[O2].b2 = 1e-10;     // Pa^-1
    dsl->params_[O2].q2 = 24624.1;   // J/mol  (-ΔU2)
    dsl->params_[N2].m1 = 1.5;
    dsl->params_[N2].b1 = 1e-10;
    dsl->params_[N2].q1 = 29195.0;
    dsl->params_[N2].m2 = 1.5;
    dsl->params_[N2].b2 = 1e-10;
    dsl->params_[N2].q2 = 26260.5;
    adsorbent.adsorption_source_.AddIsotherm({O2, N2}, std::move(dsl));
    adsorbent.adsorption_source_.heat_of_adsorption_[O2] = 13221.44;  // J/mol (-ΔH)
    adsorbent.adsorption_source_.heat_of_adsorption_[N2] = 23430.4;   // J/mol (-ΔH)

    SaveAdsorbent(adsorbent, "lix_arora2021");
    Adsorbent adsorbent_cpy;
    LoadAdsorbent("lix_arora2021", adsorbent_cpy);
    EXPECT_EQ(adsorbent.name_, adsorbent_cpy.name_);
    EXPECT_DOUBLE_EQ(adsorbent.solid_density_, adsorbent_cpy.solid_density_);
    EXPECT_DOUBLE_EQ(adsorbent.heat_capacity_, adsorbent_cpy.heat_capacity_);

    auto* iso = dynamic_cast<DualSiteLangmuirIsotherm*>(adsorbent_cpy.adsorption_source_.GetIsotherm(O2).get());
    ASSERT_NE(iso, nullptr);
    EXPECT_DOUBLE_EQ(iso->params_.at(O2).m2, 0.6075);
    EXPECT_DOUBLE_EQ(iso->params_.at(O2).q1, 21774.6);
    EXPECT_DOUBLE_EQ(iso->params_.at(N2).q1, 29195.0);
    EXPECT_DOUBLE_EQ(iso->params_.at(N2).q2, 26260.5);

    EXPECT_DOUBLE_EQ(adsorbent_cpy.adsorption_source_.heat_of_adsorption_.at(O2), 13221.44);
    EXPECT_DOUBLE_EQ(adsorbent_cpy.adsorption_source_.heat_of_adsorption_.at(N2), 23430.4);
  }

  TEST(AdsorbentTest, LiLSXZeolite) {
    // Table 2 parameters for LiLSX zeolite (~27% Li exchange, low silica X)
    const std::string O2 = "Oxygen";
    const std::string N2 = "Nitrogen";

    Adsorbent adsorbent;
    adsorbent.name_ = "LiLSX Zeolite";
    adsorbent.solid_density_ = 1035.0;  // kg/m3
    adsorbent.heat_capacity_ = 1172.0;  // J/(kg·K)
    adsorbent.ep_ = 0.33;  // particle porosity (Table 2)
    adsorbent.thermal_conductivity_ = 0.15;

    auto dsl = std::make_unique<DualSiteLangmuirIsotherm>(O2);
    dsl->params_[O2].m1 = 1.5;
    dsl->params_[O2].b1 = 1e-10;
    dsl->params_[O2].q1 = 19606.9;
    dsl->params_[O2].m2 = 1.5;
    dsl->params_[O2].b2 = 1e-10;
    dsl->params_[O2].q2 = 22906.2;
    dsl->params_[N2].m1 = 1.5;
    dsl->params_[N2].b1 = 1e-10;
    dsl->params_[N2].q1 = 27923.9;
    dsl->params_[N2].m2 = 0.6613;
    dsl->params_[N2].b2 = 1e-10;
    dsl->params_[N2].q2 = 27923.9;
    adsorbent.adsorption_source_.AddIsotherm({O2, N2}, std::move(dsl));
    adsorbent.adsorption_source_.heat_of_adsorption_[O2] = 14661.44;  // J/mol (-ΔH)
    adsorbent.adsorption_source_.heat_of_adsorption_[N2] = 25451.65;  // J/mol (-ΔH)

    auto ldf_o2 = std::make_unique<LinearDrivingForce>();
    auto mtc_o2 = std::make_unique<MassTransferCoefficientConstant>();
    mtc_o2->constant_mtc_ = 104.07;  // s^-1
    ldf_o2->mtc_model_ = std::move(mtc_o2);
    adsorbent.adsorption_source_.mass_transfer_models_[O2] = std::move(ldf_o2);

    auto ldf_n2 = std::make_unique<LinearDrivingForce>();
    auto mtc_n2 = std::make_unique<MassTransferCoefficientConstant>();
    mtc_n2->constant_mtc_ = 42.84;  // s^-1
    ldf_n2->mtc_model_ = std::move(mtc_n2);
    adsorbent.adsorption_source_.mass_transfer_models_[N2] = std::move(ldf_n2);

    SaveAdsorbent(adsorbent, "lilsx_arora2021");
    Adsorbent adsorbent_cpy;
    LoadAdsorbent("lilsx_arora2021", adsorbent_cpy);
    EXPECT_EQ(adsorbent.name_, adsorbent_cpy.name_);
    EXPECT_DOUBLE_EQ(adsorbent.solid_density_, adsorbent_cpy.solid_density_);
    EXPECT_DOUBLE_EQ(adsorbent.heat_capacity_, adsorbent_cpy.heat_capacity_);

    auto* iso = dynamic_cast<DualSiteLangmuirIsotherm*>(adsorbent_cpy.adsorption_source_.GetIsotherm(O2).get());
    ASSERT_NE(iso, nullptr);
    EXPECT_DOUBLE_EQ(iso->params_.at(O2).q1, 19606.9);
    EXPECT_DOUBLE_EQ(iso->params_.at(O2).q2, 22906.2);
    EXPECT_DOUBLE_EQ(iso->params_.at(N2).m2, 0.6613);
    EXPECT_DOUBLE_EQ(iso->params_.at(N2).q1, 27923.9);

    auto* mtm_o2 = dynamic_cast<MassTransferCoefficientConstant*>(
      adsorbent_cpy.adsorption_source_.mass_transfer_models_.at(O2)->mtc_model_.get());
    ASSERT_NE(mtm_o2, nullptr);
    EXPECT_DOUBLE_EQ(mtm_o2->constant_mtc_, 104.07);

    auto* mtm_n2 = dynamic_cast<MassTransferCoefficientConstant*>(
      adsorbent_cpy.adsorption_source_.mass_transfer_models_.at(N2)->mtc_model_.get());
    ASSERT_NE(mtm_n2, nullptr);
    EXPECT_DOUBLE_EQ(mtm_n2->constant_mtc_, 42.84);
  }

  TEST(AdsorbentTest, Zeolite5A) {
    // Table 2 parameters for 5A zeolite
    const std::string O2 = "Oxygen";
    const std::string N2 = "Nitrogen";

    Adsorbent adsorbent;
    adsorbent.name_ = "Zeolite 5A";
    adsorbent.solid_density_ = 1160.0;  // kg/m3
    adsorbent.heat_capacity_ = 1338.9;  // J/(kg·K)
    adsorbent.ep_ = 0.65;  // particle porosity (Table 2)
    adsorbent.thermal_conductivity_ = 0.15;

    auto dsl = std::make_unique<DualSiteLangmuirIsotherm>(O2);
    dsl->params_[O2].m1 = 1.5;
    dsl->params_[O2].b1 = 1e-10;
    dsl->params_[O2].q1 = 22312.1;
    dsl->params_[O2].m2 = 1.5;
    dsl->params_[O2].b2 = 1e-10;
    dsl->params_[O2].q2 = 22312.1;
    dsl->params_[N2].m1 = 1.5;
    dsl->params_[N2].b1 = 1e-10;
    dsl->params_[N2].q1 = 26140.8;
    dsl->params_[N2].m2 = 1.5;
    dsl->params_[N2].b2 = 1e-10;
    dsl->params_[N2].q2 = 23284.3;
    adsorbent.adsorption_source_.AddIsotherm({O2, N2}, std::move(dsl));
    adsorbent.adsorption_source_.heat_of_adsorption_[O2] = 13221.44;  // J/mol (-ΔH)
    adsorbent.adsorption_source_.heat_of_adsorption_[N2] = 22886.48;  // J/mol (-ΔH)

    auto ldf_o2 = std::make_unique<LinearDrivingForce>();
    auto mtc_o2 = std::make_unique<MassTransferCoefficientConstant>();
    mtc_o2->constant_mtc_ = 12.07;  // s^-1
    ldf_o2->mtc_model_ = std::move(mtc_o2);
    adsorbent.adsorption_source_.mass_transfer_models_[O2] = std::move(ldf_o2);

    auto ldf_n2 = std::make_unique<LinearDrivingForce>();
    auto mtc_n2 = std::make_unique<MassTransferCoefficientConstant>();
    mtc_n2->constant_mtc_ = 4.02;  // s^-1
    ldf_n2->mtc_model_ = std::move(mtc_n2);
    adsorbent.adsorption_source_.mass_transfer_models_[N2] = std::move(ldf_n2);

    SaveAdsorbent(adsorbent, "zeolite_5a_arora2021");
    Adsorbent adsorbent_cpy;
    LoadAdsorbent("zeolite_5a_arora2021", adsorbent_cpy);
    EXPECT_EQ(adsorbent.name_, adsorbent_cpy.name_);
    EXPECT_DOUBLE_EQ(adsorbent.solid_density_, adsorbent_cpy.solid_density_);
    EXPECT_DOUBLE_EQ(adsorbent.heat_capacity_, adsorbent_cpy.heat_capacity_);

    auto* iso = dynamic_cast<DualSiteLangmuirIsotherm*>(adsorbent_cpy.adsorption_source_.GetIsotherm(O2).get());
    ASSERT_NE(iso, nullptr);
    EXPECT_DOUBLE_EQ(iso->params_.at(O2).q1, 22312.1);
    EXPECT_DOUBLE_EQ(iso->params_.at(O2).q2, 22312.1);
    EXPECT_DOUBLE_EQ(iso->params_.at(N2).q1, 26140.8);
    EXPECT_DOUBLE_EQ(iso->params_.at(N2).q2, 23284.3);

    auto* mtm_o2 = dynamic_cast<MassTransferCoefficientConstant*>(
      adsorbent_cpy.adsorption_source_.mass_transfer_models_.at(O2)->mtc_model_.get());
    ASSERT_NE(mtm_o2, nullptr);
    EXPECT_DOUBLE_EQ(mtm_o2->constant_mtc_, 12.07);

    auto* mtm_n2 = dynamic_cast<MassTransferCoefficientConstant*>(
      adsorbent_cpy.adsorption_source_.mass_transfer_models_.at(N2)->mtc_model_.get());
    ASSERT_NE(mtm_n2, nullptr);
    EXPECT_DOUBLE_EQ(mtm_n2->constant_mtc_, 4.02);
  }

}
