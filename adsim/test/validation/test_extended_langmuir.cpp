/*
 * Skarstrom - pressure and temperature swing adsorption simulator
 * Copyright (c) 2024-2026 Sam Affleck and contributors
 *
 * Licensed under the terms in the LICENSE file at the repository root.
 */

// Thermodynamic consistency of the EXTENDED (competitive) Langmuir isotherm.
//
// The single- and dual-site consistency tests in test_thermodynamic_consistency
// exercise the pure-component isotherm forms. The extended Langmuir is the
// multi-component generalisation the solver actually uses for competitive
// adsorption:
//
//     q_i = m1_i b_i(T) P y_i / ( 1 + sum_j b_j(T) P y_j ),
//     b_j(T) = b1_j exp( q1_j / (R T) ).
//
// (ExtendedLangmuirIsotherm::Update in isotherm.cpp.) This test drives the REAL
// isotherm through a closed, spatially-uniform bed and checks the equilibrium
// loadings the solver seeds against a closed-form reference written straight from
// the algebra above -- the reference never calls the solver's evaluation. Because
// the bed starts uniform it starts at equilibrium (Equations seeds q = q* from
// the gas initial condition), so the loadings at t = 0 ARE the competitive
// equilibrium values; no time integration is needed for the comparison.
//
// Three consistency properties are asserted, each a defining feature of the
// competitive form and none satisfiable by an independent single-component fit:
//
//   1. LOADINGS MATCH THE COMPETITIVE FORMULA at a genuinely competitive state,
//      where the shared denominator suppresses each loading well below its
//      single-component value.
//
//   2. SINGLE-COMPONENT REDUCTION: with only one adsorbate present (y_j = 0 for
//      the others) the extended form collapses exactly to the single-site
//      Langmuir q = m1 b P y / (1 + b P y). We check this against the same
//      IsothermFunction::Langmuir the single-component tests use.
//
//   3. COMPETITION SUPPRESSES LOADING: adding a second adsorbate at fixed partial
//      pressure of the first must LOWER the first's loading (shared sites), and
//      the selectivity q_1/q_2 must equal (m1_1 b_1 y_1)/(m1_2 b_2 y_2),
//      independent of the total pressure -- the competitive Langmuir signature.

#include "validation_helpers.h"

namespace adsim {

  using namespace adsim::validation;

  namespace {

    // Competitive Langmuir loading [mol/kg] for one component, written from the
    // isotherm algebra. `params` maps component -> (m1, b1, q1); `y` maps
    // component -> mole fraction. p_i = P y_i.
    struct ELParams { double m1, b1, q1; };
    double ExtendedLangmuir(const std::string& comp,
                            const std::unordered_map<std::string, ELParams>& params,
                            const std::unordered_map<std::string, double>& y,
                            double P, double T) {
      double denom = 1.0;
      for (const auto& [c, p] : params) {
        const double b = p.b1 * std::exp(p.q1 / (kR * T));
        denom += b * P * y.at(c);
      }
      const auto& pc = params.at(comp);
      const double b = pc.b1 * std::exp(pc.q1 / (kR * T));
      return pc.m1 * b * P * y.at(comp) / denom;
    }

  } // namespace

  // Closed, uniform, isothermal bed carrying an ExtendedLangmuir isotherm on two
  // adsorbates. Reads the seeded equilibrium loadings straight from the log.
  class ExtendedLangmuirTest : public UserDataTest {
  public:
    LayerID layer_id{};
    AdsorbentID ads_id{};

    static constexpr double kT = 300.0;      // [K]
    static constexpr double kSolidDensity = 1000.0;

    std::unordered_map<std::string, ELParams> params;

    void BuildBed() {
      AddComponents({"Nitrogen", "Oxygen"});
      auto ids = AddLayer("Packing");
      layer_id = ids.first;
      ads_id = ids.second;

      auto& layer = GetLayer(layer_id);
      layer.domain_.length = 1.0;
      layer.domain_.n_cells = 10;
      layer.eb_ = 0.4;
      layer.particle_diameter_ = 0.002;
      std::get<CylinderBed>(data.bed.shape_).diameter = 0.1;

      auto& adsorbent = GetAdsorbent(ads_id);
      adsorbent.ep_ = 0.0;
      adsorbent.solid_density_ = kSolidDensity;
      adsorbent.heat_capacity_ = 900.0;
      adsorbent.thermal_conductivity_ = 0.0;

      AddConstantInternalHeatTransferCoefficientModel(0.0, layer_id);
      AddConstantExternalHeatTransferCoefficientModel(0.0, layer_id);
    }

    // Install an ExtendedLangmuir isotherm on both adsorbates with fast LDF, so
    // the seeded loading equals the competitive equilibrium value.
    void InstallExtendedLangmuir(const std::unordered_map<std::string, ELParams>& p) {
      params = p;
      auto& adsorbent = GetAdsorbent(ads_id);
      auto iso = std::make_unique<ExtendedLangmuirIsotherm>();
      for (const auto& [comp, pr] : p) {
        iso->params_[comp].m1 = pr.m1;
        iso->params_[comp].b1 = pr.b1;
        iso->params_[comp].q1 = pr.q1;
      }
      std::vector<std::string> comps;
      for (const auto& [comp, pr] : p) comps.push_back(comp);
      adsorbent.adsorption_source_.AddIsotherm(comps, std::move(iso));
      for (const auto& comp : comps) {
        AddLinearDrivingForceModel(ads_id, comp, 1.0);
        SetHOA(ads_id, comp, 0.0);
      }
    }

    // Run a short sealed hold and return the seeded equilibrium loadings of the
    // requested components (mean over the interior at t = 0), from ONE solver run.
    std::unordered_map<std::string, double> SeededLoadings(
        double P, const std::unordered_map<std::string, double>& y,
        const std::vector<std::string>& components) {
      AddInitialConditions(P, kT, kT, y, 0.0);
      AddStep("Hold", 1.0);
      AddClosedBC(BoundaryLocation::BOTTOM, "Hold");
      AddClosedBC(BoundaryLocation::TOP, "Hold");
      ConfigureCycle(1);
      data.settings.reltol = 1e-10;
      data.settings.max_step = 0.5;

      CoupledSolver solver{data};
      std::string error;
      EXPECT_TRUE(data.PerformChecks(error)) << error;
      EXPECT_GE(solver.IntegrateCycle(data.cycle, false), 0) << "Solver did not complete";

      const auto& log = solver.GetLog();
      const auto& dom = log.GetDomain();
      std::unordered_map<std::string, double> out;
      for (const auto& component : components) {
        const auto q = log.GetComponentVariable(0, 0,
            ComponentVariableType::ADSORBED_CONCENTRATION, component);
        EXPECT_FALSE(q.values.empty());
        out[component] = q.values.empty() ? 0.0 : InteriorMean(q.values.front(), dom.volume);
      }
      return out;
    }
  };

  // ===========================================================================
  // 1. Loadings match the competitive formula
  // ===========================================================================

  TEST_F(ExtendedLangmuirTest, CompetitiveLoadingsMatchClosedForm) {
    BuildBed();
    // Both adsorbates well up their isotherms so the shared denominator bites.
    InstallExtendedLangmuir({
      {"Nitrogen", {3.0, 5.0e-9, 18000.0}},
      {"Oxygen",   {4.0, 3.0e-9, 22000.0}},
    });

    const double P = 5.0e5;  // [Pa]
    const std::unordered_map<std::string, double> y{{"Nitrogen", 0.6}, {"Oxygen", 0.4}};

    const auto q_sim = SeededLoadings(P, y, {"Nitrogen", "Oxygen"});
    const double qN_sim = q_sim.at("Nitrogen");
    const double qO_sim = q_sim.at("Oxygen");
    const double qN_ref = ExtendedLangmuir("Nitrogen", params, y, P, kT);
    const double qO_ref = ExtendedLangmuir("Oxygen", params, y, P, kT);

    // GUARD: the state must be genuinely competitive -- the shared denominator
    // must suppress the loading well below the single-component value, or the
    // extended form is indistinguishable from an isolated Langmuir.
    const double bN = params["Nitrogen"].b1 * std::exp(params["Nitrogen"].q1 / (kR * kT));
    const double qN_isolated = params["Nitrogen"].m1 * bN * P * y.at("Nitrogen") /
                               (1.0 + bN * P * y.at("Nitrogen"));
    ASSERT_LT(qN_ref, 0.8 * qN_isolated)
      << "Competition is weak (competitive " << qN_ref << " vs isolated " << qN_isolated
      << " mol/kg); the shared-denominator term is not exercised";

    EXPECT_NEAR(qN_sim, qN_ref, 1e-4 * qN_ref)
      << "Extended-Langmuir N2 loading disagrees with the competitive formula "
      << "(simulated " << qN_sim << ", closed form " << qN_ref << " mol/kg)";

    EXPECT_NEAR(qO_sim, qO_ref, 1e-4 * qO_ref)
      << "Extended-Langmuir O2 loading disagrees with the competitive formula";
  }

  // ===========================================================================
  // 2. Single-component reduction to single-site Langmuir
  // ===========================================================================

  TEST_F(ExtendedLangmuirTest, ReducesToSingleSiteLangmuirWithOneAdsorbate) {
    BuildBed();
    const ELParams pN{3.5, 4.0e-9, 20000.0};
    InstallExtendedLangmuir({
      {"Nitrogen", pN},
      {"Oxygen",   {4.0, 3.0e-9, 22000.0}},
    });

    const double P = 3.0e5;
    // Only Nitrogen present: Oxygen mole fraction zero.
    const std::unordered_map<std::string, double> y{{"Nitrogen", 1.0}, {"Oxygen", 0.0}};

    const double qN_sim = SeededLoadings(P, y, {"Nitrogen"}).at("Nitrogen");

    // Reference: the single-site Langmuir the pure-component tests use.
    LangmuirIsothermParameters lp;
    lp.m1 = pN.m1; lp.b1 = pN.b1; lp.q1 = pN.q1;
    double qN_langmuir = 0.0;
    IsothermFunction::Langmuir(qN_langmuir, kT, P, /*yi=*/1.0, lp);

    ASSERT_GT(qN_langmuir, 0.1) << "Loading negligible; reduction is not exercised";
    EXPECT_NEAR(qN_sim, qN_langmuir, 1e-4 * qN_langmuir)
      << "Extended Langmuir with a single adsorbate does not reduce to single-site "
      << "Langmuir (extended " << qN_sim << ", single-site " << qN_langmuir << " mol/kg)";
  }

  // ===========================================================================
  // 3. Competition lowers loading; selectivity is pressure-independent
  // ===========================================================================

  TEST_F(ExtendedLangmuirTest, AddingCompetitorLowersLoadingAndFixesSelectivity) {
    BuildBed();
    InstallExtendedLangmuir({
      {"Nitrogen", {3.0, 5.0e-9, 18000.0}},
      {"Oxygen",   {4.0, 3.0e-9, 22000.0}},
    });

    const double P = 5.0e5;

    // Alone: only Nitrogen. Its partial pressure is 0.5 * P.
    const std::unordered_map<std::string, double> y_alone{{"Nitrogen", 0.5}, {"Oxygen", 0.0}};
    const double qN_alone = ExtendedLangmuir("Nitrogen", params, y_alone, P, kT);

    // With a competitor: SAME Nitrogen partial pressure (y_N2 = 0.5) but Oxygen
    // now fills the rest, at the same total pressure.
    const std::unordered_map<std::string, double> y_comp{{"Nitrogen", 0.5}, {"Oxygen", 0.5}};
    const double qN_comp = ExtendedLangmuir("Nitrogen", params, y_comp, P, kT);

    // Adding the competitor at fixed N2 partial pressure must lower N2 loading.
    EXPECT_LT(qN_comp, qN_alone)
      << "Adding a competitor did not reduce the loading at fixed partial pressure";

    // Selectivity q_N2/q_O2 = (m_N2 b_N2 y_N2)/(m_O2 b_O2 y_O2), independent of P.
    const double bN = params["Nitrogen"].b1 * std::exp(params["Nitrogen"].q1 / (kR * kT));
    const double bO = params["Oxygen"].b1 * std::exp(params["Oxygen"].q1 / (kR * kT));
    const double sel_expected =
      (params["Nitrogen"].m1 * bN * 0.5) / (params["Oxygen"].m1 * bO * 0.5);
    for (double P_test : {1.0e5, 5.0e5, 20.0e5}) {
      const double qN = ExtendedLangmuir("Nitrogen", params, y_comp, P_test, kT);
      const double qO = ExtendedLangmuir("Oxygen", params, y_comp, P_test, kT);
      EXPECT_NEAR(qN / qO, sel_expected, 1e-9 * sel_expected)
        << "Competitive selectivity is not pressure-independent at P = " << P_test;
    }
  }

} // namespace adsim
