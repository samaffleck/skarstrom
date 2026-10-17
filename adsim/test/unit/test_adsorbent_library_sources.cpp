// Builds (and round-trips) the curated library blobs added for roadmap item 10
// from the values printed in their open-access sources. Running these tests
// re-saves the blobs into data/adsorbents/, exactly like the older
// generators in test_adsorbents.cpp; data/adsorbents/library.json carries
// the per-parameter provenance and scripts/docs/gen_adsorbent_library.py checks the
// catalogue against these blobs.
//
// Rules: every number below is transcribed from the cited table/equation (unit
// conversions only: kPa^-1 -> Pa^-1 is x1e-3, bar^-1 -> Pa^-1 is x1e-5,
// MPa^-1 -> Pa^-1 is x1e-6, mmol/g == mol/kg, kJ/mol -> J/mol is x1e3). The only
// computed values are (a) the Toth b0 re-referencing for Low 2023 (exact algebra,
// see ReferenceB0) and (b) heat capacities evaluated from the authors' fitted
// cp(T) correlations at 298.15 K. Quantities a source does not state keep the
// engine default and are marked "not stated in source" in library.json.

#include <cmath>
#include <filesystem>
#include <memory>
#include <string>
#include <vector>

#include <gtest/gtest.h>

#include "adsim/equations/isotherm.h"
#include "adsim/math/constants.h"
#include "adsim/objects/porous_media.h"

namespace adsim {
namespace {

  const std::filesystem::path kLibDir = std::filesystem::current_path() / ".." / "data" / "adsorbents";

  // Save to data/adsorbents/<key>.json and load it back. `key` is the
  // entry's library key (library.json), never the display name: spaces and
  // parentheses in data file names break the asset-copy build steps of front
  // ends (see cmake/check_data_names.cmake).
  void SaveAndReload(const Adsorbent& a, const std::string& key, Adsorbent& out) {
    std::filesystem::create_directories(kLibDir);
    const auto path = kLibDir / (key + ".json");
    SaveAdsorbentToJSON(path.string(), a);
    LoadAdsorbentFromJSON(path.string(), out);
  }

  void ExpectSameScalars(const Adsorbent& a, const Adsorbent& b) {
    EXPECT_EQ(a.name_, b.name_);
    EXPECT_DOUBLE_EQ(a.solid_density_.Get(), b.solid_density_.Get());
    EXPECT_DOUBLE_EQ(a.heat_capacity_.Get(), b.heat_capacity_.Get());
    EXPECT_DOUBLE_EQ(a.thermal_conductivity_.Get(), b.thermal_conductivity_.Get());
    EXPECT_DOUBLE_EQ(a.ep_.Get(), b.ep_.Get());
  }

  // Low et al. 2023/2024 write b(T) = b0 exp(-dH/(R T)) with an absolute
  // pre-exponential b0. The engine's TOTH_2 writes b(T) = b0' exp(dH0/(R T0)
  // (T0/T - 1)), i.e. b0' is b at T0. Equating both forms for every T gives
  // b0' = b0 exp(dH0/(R T0)) with dH0 = -dH > 0 (exact; no fitting).
  double ReferenceB0(double b0_per_Pa, double dH0, double T0) {
    return b0_per_Pa * std::exp(dH0 / (constants::kIdealGasConstant * T0));
  }

  const std::string CO2 = "CarbonDioxide";
  const std::string N2 = "Nitrogen";
  const std::string H2 = "Hydrogen";
  const std::string CH4 = "Methane";
  const std::string O2 = "Oxygen";

}  // namespace

// ---------------------------------------------------------------------------
// Low, Danaci, Azzan, Woodward, Petit (2023) J. Chem. Eng. Data 68, 3499.
// https://doi.org/10.1021/acs.jced.3c00401  (CC BY 4.0)
// CO2: td-Toth, Table 7 (eqs 8-11). N2: SSL, Table 8 (eq 15). Densities and
// porosity: Table 5. cp(T) = a T^-2 + b T - c (eq 7), Table 6.
// ---------------------------------------------------------------------------
TEST(AdsorbentLibrarySources, LewatitVPOC1065_Low2023) {
  Adsorbent a;
  a.name_ = "Lewatit VP OC 1065 (Low 2023)";
  a.solid_density_ = 744.0;  // Table 5, rho_particle 0.744 g/cm3
  a.ep_ = 0.338;             // Table 5, eps_particle
  {
    const double T = 298.15, ca = -3.23e4, cb = 0.00227, cc = -0.994;  // Table 6 (J K/g, J/(g K2), J/(g K))
    a.heat_capacity_ = 1000.0 * (ca / (T * T) + cb * T - cc);           // eq 7 at 298.15 K -> J/(kg K)
  }
  // thermal conductivity: not stated in the source -> engine default (0.5).

  auto toth = std::make_unique<Toth2Isotherm>(CO2);
  auto& t = toth->params_[CO2];
  t.q_sat = 3.00;                                   // q_inf,0 [mmol/g]
  t.x = 0.7;                                        // chi
  t.T_0 = 298.15;                                   // T0 [K]
  t.dH0 = 70.1e3;                                   // -dH [J/mol]
  t.b0 = ReferenceB0(2.0e-11 * 1e-3, 70.1e3, 298.15);  // b0 = 2.0e-11 kPa^-1
  t.w0 = 0.431;                                     // tau0
  t.alpha = 1.4e-10;                                // alpha
  a.adsorption_source_.AddIsotherm({CO2}, std::move(toth));

  auto ssl = std::make_unique<LangmuirIsotherm>(N2);
  ssl->params_[N2].m1 = 0.315;            // q_inf [mmol/g]
  ssl->params_[N2].b1 = 5.06e-8 * 1e-3;   // b0 = 5.06e-8 kPa^-1
  ssl->params_[N2].q1 = 20.30e3;          // -dH [J/mol]
  a.adsorption_source_.AddIsotherm({N2}, std::move(ssl));

  a.adsorption_source_.heat_of_adsorption_[CO2] = 70.1e3;   // Table 7 (td-Toth -dH)
  a.adsorption_source_.heat_of_adsorption_[N2] = 20.30e3;   // Table 8 (SSL -dH)

  Adsorbent b;
  SaveAndReload(a, "lewatit_vp_oc_1065_low2023", b);
  ExpectSameScalars(a, b);
  // At T0 the engine affinity equals the paper's b0 exp(-dH/(R T0)).
  const double bT0 = 2.0e-14 * std::exp(70.1e3 / (constants::kIdealGasConstant * 298.15));
  auto* iso = dynamic_cast<Toth2Isotherm*>(b.adsorption_source_.GetIsotherm(CO2).get());
  ASSERT_NE(iso, nullptr);
  EXPECT_NEAR(iso->params_.at(CO2).b0.Get(), bT0, 1e-12 * bT0);
  EXPECT_NEAR(a.heat_capacity_.Get(), 1307.44, 0.01);
}

TEST(AdsorbentLibrarySources, PuroliteA110_Low2023) {
  Adsorbent a;
  a.name_ = "Purolite A110 (Low 2023)";
  a.solid_density_ = 663.0;  // Table 5, rho_particle 0.663 g/cm3
  a.ep_ = 0.422;             // Table 5, eps_particle
  {
    const double T = 298.15, ca = -5.86e4, cb = 0.000994, cc = -1.69;  // Table 6
    a.heat_capacity_ = 1000.0 * (ca / (T * T) + cb * T - cc);          // eq 7 at 298.15 K
  }

  auto toth = std::make_unique<Toth2Isotherm>(CO2);
  auto& t = toth->params_[CO2];
  t.q_sat = 3.74;
  t.x = 0.3;
  t.T_0 = 298.15;
  t.dH0 = 70.1e3;
  t.b0 = ReferenceB0(3.4e-11 * 1e-3, 70.1e3, 298.15);  // b0 = 3.4e-11 kPa^-1
  t.w0 = 0.401;
  t.alpha = 1.43e-4;
  a.adsorption_source_.AddIsotherm({CO2}, std::move(toth));

  auto ssl = std::make_unique<LangmuirIsotherm>(N2);
  ssl->params_[N2].m1 = 0.0419;
  ssl->params_[N2].b1 = 75.6e-8 * 1e-3;   // b0 = 75.6e-8 kPa^-1
  ssl->params_[N2].q1 = 18.65e3;
  a.adsorption_source_.AddIsotherm({N2}, std::move(ssl));

  a.adsorption_source_.heat_of_adsorption_[CO2] = 70.1e3;
  a.adsorption_source_.heat_of_adsorption_[N2] = 18.65e3;

  Adsorbent b;
  SaveAndReload(a, "purolite_a110_low2023", b);
  ExpectSameScalars(a, b);
  EXPECT_NEAR(a.heat_capacity_.Get(), 1327.14, 0.01);
}

// ---------------------------------------------------------------------------
// Grossmann, Saenz-Cavazos, Ferru, Williams, Mazzotti (2025) Ind. Eng. Chem.
// Res. 64, 7165. https://doi.org/10.1021/acs.iecr.4c04581  (CC BY 4.0)
// Dry CO2 temperature-dependent Toth, Table 4 (eqs 8-11 are the engine's
// TOTH_2 form exactly, b0 referenced to T0 = 308 K).
// ---------------------------------------------------------------------------
TEST(AdsorbentLibrarySources, TriamineGraftedAlumina_Grossmann2025) {
  Adsorbent a;
  a.name_ = "Triamine-Grafted Alumina (Grossmann 2025)";
  // density, cp, conductivity, particle voidage: not stated in the source
  // (characterised in its ref 34) -> engine defaults.

  auto toth = std::make_unique<Toth2Isotherm>(CO2);
  auto& t = toth->params_[CO2];
  t.q_sat = 1.17;          // q_s,0 [mmol/g]
  t.x = 0.49;              // chi
  t.T_0 = 308.0;           // T0 [K]
  t.b0 = 2.96e3 * 1e-3;    // b0 = 2.96e3 kPa^-1
  t.dH0 = 77.8e3;          // dH = -77.8 kJ/mol
  t.w0 = 0.25;             // t0
  t.alpha = 0.19;          // alpha
  a.adsorption_source_.AddIsotherm({CO2}, std::move(toth));
  a.adsorption_source_.heat_of_adsorption_[CO2] = 77.8e3;  // Table 4

  Adsorbent b;
  SaveAndReload(a, "triamine_alumina_grossmann2025", b);
  ExpectSameScalars(a, b);
}

// ---------------------------------------------------------------------------
// Krishnamurthy, Blom, Ferrari, Brandani (2020) Adsorption 26, 711.
// https://doi.org/10.1007/s10450-019-00162-x  (CC BY 4.0)
// CO2/N2 dual-site Langmuir (eq 1, b0 in kPa^-1), Table 1; saturation
// capacities shared "for thermodynamic consistency" -> extended (competitive)
// DSL. Particle porosity: Table 2.
// ---------------------------------------------------------------------------
TEST(AdsorbentLibrarySources, CPO27Ni_Krishnamurthy2020) {
  Adsorbent a;
  a.name_ = "CPO-27-Ni (Krishnamurthy 2020)";
  a.ep_ = 0.383;  // Table 2, porosity (mercury intrusion)

  auto dsl = std::make_unique<ExtendedDualSiteLangmuirIsotherm>(CO2);
  auto& c = dsl->params_[CO2];
  c.m1 = 7.02;  c.b1 = 1.5e-8 * 1e-3;  c.q1 = 39.9e3;
  c.m2 = 6.95;  c.b2 = 9e-8 * 1e-3;    c.q2 = 23.6e3;
  auto& n = dsl->params_[N2];
  n.m1 = 7.02;  n.b1 = 8.3e-7 * 1e-3;  n.q1 = 17.3e3;
  n.m2 = 6.95;  n.b2 = 8.3e-7 * 1e-3;  n.q2 = 17.3e3;
  a.adsorption_source_.AddIsotherm({CO2, N2}, std::move(dsl));
  // Section 4: "the heat of adsorption of CO2 and N2 are found to be -39.7
  // and -17.3 kJ/mol".
  a.adsorption_source_.heat_of_adsorption_[CO2] = 39.7e3;
  a.adsorption_source_.heat_of_adsorption_[N2] = 17.3e3;

  Adsorbent b;
  SaveAndReload(a, "cpo_27_ni_krishnamurthy2020", b);
  ExpectSameScalars(a, b);
  // Table 1 lists the Henry constant at 303 K: 80.5 mol/(kg bar) for CO2.
  const double R = constants::kIdealGasConstant, T = 303.0;
  const double henry = (7.02 * 1.5e-8 * std::exp(39.9e3 / (R * T)) + 6.95 * 9e-8 * std::exp(23.6e3 / (R * T))) * 100.0;
  EXPECT_NEAR(henry, 80.5, 1.0);
}

// ---------------------------------------------------------------------------
// Azzan, Danaci, Petit, Pini (2023) J. Chem. Eng. Data 68, 3512.
// https://doi.org/10.1021/acs.jced.3c00504  (CC BY 4.0)
// Single-site Langmuir for absolute adsorption (eq 12, p in bar), Table 8;
// n_s,b is shared by CO2/N2/H2 of a material "for thermodynamic consistency"
// -> extended (competitive) SSL. For H-Y the published N2 parameters do not
// reproduce the paper's own tabulated N2 data (SI Table S3: ~10x higher
// loadings), so N2 is left out of that entry rather than "corrected".
// ---------------------------------------------------------------------------
namespace {
  struct SslRow { const char* comp; double b0_1e7_per_bar; double dU_kJ; };
  void BuildAzzan(const char* name, const std::string& key, double ns, const std::vector<SslRow>& rows) {
    Adsorbent a;
    a.name_ = name;
    // particle density, cp, conductivity, voidage, particle size: not stated
    // (powder samples) -> engine defaults.
    auto ssl = std::make_unique<ExtendedLangmuirIsotherm>(rows[0].comp);
    for (const auto& r : rows) {
      auto& p = ssl->params_[r.comp];
      p.m1 = ns;
      p.b1 = r.b0_1e7_per_bar * 1e-7 * 1e-5;  // (x1e-7) bar^-1 -> Pa^-1
      p.q1 = r.dU_kJ * 1e3;                   // -dU_b [J/mol]
    }
    std::vector<std::string> comps;
    for (const auto& r : rows) comps.emplace_back(r.comp);
    a.adsorption_source_.AddIsotherm(comps, std::move(ssl));
    for (const auto& r : rows) a.adsorption_source_.heat_of_adsorption_[r.comp] = r.dU_kJ * 1e3;
    Adsorbent b;
    SaveAndReload(a, key, b);
    ExpectSameScalars(a, b);
    for (const auto& r : rows) {
      auto* iso = dynamic_cast<ExtendedLangmuirIsotherm*>(b.adsorption_source_.GetIsotherm(r.comp).get());
      ASSERT_NE(iso, nullptr) << r.comp;
      EXPECT_DOUBLE_EQ(iso->params_.at(r.comp).m1.Get(), ns);
    }
  }
}  // namespace

TEST(AdsorbentLibrarySources, ZeoliteHY_Azzan2023) {
  // Table 8 also lists N2 (b0 = 193.42e-7 bar^-1, -dU_b = 11.35 kJ/mol); omitted, see above.
  BuildAzzan("Zeolite H-Y (Azzan 2023)", "zeolite_h_y_azzan2023", 5.86, {{"CarbonDioxide", 370.44, 23.31}, {"Hydrogen", 1406.86, 8.96}});
}
TEST(AdsorbentLibrarySources, ZeoliteNaY_Azzan2023) {
  BuildAzzan("Zeolite Na-Y (Azzan 2023)", "zeolite_na_y_azzan2023", 7.03,
             {{"CarbonDioxide", 102.10, 33.21}, {"Nitrogen", 1347.59, 13.44}, {"Hydrogen", 1035.24, 10.79}});
}
TEST(AdsorbentLibrarySources, ZeoliteNaTMAY_Azzan2023) {
  BuildAzzan("Zeolite NaTMA-Y (Azzan 2023)", "zeolite_natma_y_azzan2023", 6.91,
             {{"CarbonDioxide", 88.59, 31.98}, {"Nitrogen", 2078.48, 12.00}, {"Hydrogen", 744.94, 11.16}});
}

// ---------------------------------------------------------------------------
// Low, Danaci, Azzan, Jiayi, Yong, Itskou, Petit (2024) Energy Fuels 38, 11947.
// https://doi.org/10.1021/acs.energyfuels.4c01368  (CC BY 4.0)
// CO2 DSL (eqs 7-9, P in bar), N2 SSL, Table 1. Fitted separately (no shared
// saturation capacity) -> non-competitive DSL + SSL. Pellet densities and
// porosity: Table 5. cp = a (T+273.15)^-2 + b (T+273.15) - c (eq 6, T in C).
// ---------------------------------------------------------------------------
TEST(AdsorbentLibrarySources, TIFSIX3Ni_Low2024) {
  Adsorbent a;
  a.name_ = "TIFSIX-3-Ni (Low 2024)";
  a.solid_density_ = 810.0;  // Table 5, rho_pellet 0.81 g/cm3 (binderless pellet)
  a.ep_ = 0.36;              // Table 5, eps_pellet
  {
    const double Tk = 25.0 + 273.15, ca = 1.077e5, cb = 5.65e-3, cc = 1.773;  // Section 3.1 (eq 6)
    a.heat_capacity_ = 1000.0 * (ca / (Tk * Tk) + cb * Tk - cc);               // at 25 C -> J/(kg K)
  }

  auto dsl = std::make_unique<DualSiteLangmuirIsotherm>(CO2);
  auto& c = dsl->params_[CO2];
  c.m1 = 0.38737;  c.b1 = 3.4932e-5 * 1e-5;  c.q1 = 29.256e3;
  c.m2 = 1.8730;   c.b2 = 1.2200e-7 * 1e-5;  c.q2 = 59.798e3;
  a.adsorption_source_.AddIsotherm({CO2}, std::move(dsl));

  auto ssl = std::make_unique<LangmuirIsotherm>(N2);
  ssl->params_[N2].m1 = 1.7934;
  ssl->params_[N2].b1 = 1.9277e-5 * 1e-5;
  ssl->params_[N2].q1 = 21.175e3;
  a.adsorption_source_.AddIsotherm({N2}, std::move(ssl));
  // CO2: only site energies + a graphical isosteric heat (Fig. 5) -> unset.
  a.adsorption_source_.heat_of_adsorption_[N2] = 21.175e3;  // SSL -dU, Table 1

  Adsorbent b;
  SaveAndReload(a, "tifsix_3_ni_low2024", b);
  ExpectSameScalars(a, b);
  EXPECT_NEAR(a.heat_capacity_.Get(), 1123.11, 0.01);
  // Section 3.2.1: 1.11 mmol/g at 0.4 mbar and 25 C.
  const double R = constants::kIdealGasConstant, T = 298.15, P = 40.0;
  const double b1 = 3.4932e-10 * std::exp(29.256e3 / (R * T)), b2 = 1.22e-12 * std::exp(59.798e3 / (R * T));
  EXPECT_NEAR(0.38737 * b1 * P / (1 + b1 * P) + 1.873 * b2 * P / (1 + b2 * P), 1.11, 0.02);
}

// ---------------------------------------------------------------------------
// Zhang, Liu, Zhou, Li, Wei, Song (2021) ACS Omega 6, 6914.
// https://doi.org/10.1021/acsomega.0c06148  (CC BY-NC-ND 4.0)
// Sips fits at 308 K (eq 2: q = qm (bp)^(1/n) / (1 + (bp)^(1/n)), b in MPa^-1),
// Table 1. The source fits each temperature separately, so the entry is the
// 308 K fit (no temperature dependence). Isosteric heats: Section 2.2.
// ---------------------------------------------------------------------------
TEST(AdsorbentLibrarySources, CMS3KT_Zhang2021) {
  Adsorbent a;
  a.name_ = "CMS-3KT (Zhang 2021)";
  auto sips = std::make_unique<SipsIsotherm>(CH4);
  auto set = [&](const std::string& comp, double qm, double b_perMPa, double n) {
    auto& p = sips->params_[comp];
    p.q_sat = qm;
    p.b = b_perMPa * 1e-6;
    p.v = n;
  };
  set(CH4, 3.73, 1.26, 1.37);
  set(O2, 2.25, 1.74, 1.03);
  set(N2, 1.63, 1.59, 1.01);
  a.adsorption_source_.AddIsotherm({CH4, O2, N2}, std::move(sips));
  a.adsorption_source_.heat_of_adsorption_[CH4] = 20.6e3;
  a.adsorption_source_.heat_of_adsorption_[O2] = 18.1e3;
  a.adsorption_source_.heat_of_adsorption_[N2] = 19.1e3;

  Adsorbent b;
  SaveAndReload(a, "cms_3kt_zhang2021", b);
  ExpectSameScalars(a, b);
}

}  // namespace adsim
