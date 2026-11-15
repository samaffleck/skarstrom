/*
 * Skarstrom - pressure and temperature swing adsorption simulator
 * Copyright (c) 2024-2026 Sam Affleck and contributors
 *
 * Licensed under the terms in the LICENSE file at the repository root.
 */

// The heat of adsorption from the isotherm (HeatOfAdsorptionMode::Isotherm).
//
// With the mode set, the energy balance takes the isosteric heat
// q_st = R T^2 (d ln p / dT)_q of the component's isotherm at each cell's
// state instead of the stored constant. Checks:
//
//   * IsostericHeatAtState: Langmuir with b = b0 exp(Q/RT) gives Q at every
//     partial pressure, also at y = 0 and at a slightly negative y (the dilute
//     floor); a dual-site Langmuir gives a value between its two site energies
//     that falls as the strong site fills; an inert isotherm gives the fallback.
//   * Adiabatic breakthrough on a Henry isotherm k exp(Q/RT): the isotherm mode
//     with a wrong stored constant reproduces the constant mode with the right
//     one (Q), so the per-cell heat really replaces the constant.
//   * Dual-site Langmuir in the same column: the per-cell heats lie between the
//     two site energies and are lower where the bed is loaded than where it is
//     clean, and the plateau is lower than with the strong-site energy as a
//     constant.
//   * The mode is a registry parameter, and it survives a save and a load.
//
// Fast tier: 40 cells, one feed step.

#include "breakthrough_case.h"

#include "adsim/equations/isotherm_kinds.h"
#include "adsim/simulation.h"
#include "../fixtures/atomic_write.h"

#include <filesystem>

namespace adsim {

  using namespace adsim::validation;

  namespace {

    constexpr double kT0 = 300.0;     // [K] feed and initial temperature
    constexpr double kYFeed = 0.2;    // [-] adsorbing component (Oxygen) in the feed
    constexpr double kHenryK = 6.0e-11;  // [mol/(kg-Pa)]
    constexpr double kHenryQ = 38000.0;  // [J/mol]

    // Dual-site Langmuir: a small strong site and a large weak site.
    constexpr double kDslM1 = 0.5, kDslB1 = 2.0e-12, kDslQ1 = 45000.0;
    constexpr double kDslM2 = 3.0, kDslB2 = 1.0e-10, kDslQ2 = 25000.0;

    enum class Iso { Henry, DualSite };

    struct Column {
      UserData data;
      AdsorbentID ads_id{};

      // The Henry bed takes about 11 000 s to saturate and the dual-site bed
      // about 1200 s, so the dual-site run stops while the front is inside.
      Column(Iso kind, HeatOfAdsorptionMode mode, double stored_heat) {
        BreakthroughCase c;
        c.n_cells = 40;
        c.adsorbing = false;  // the isotherm is installed below
        c.molar_flow = 0.01;
        c.y_feed = kYFeed;
        c.feed_temperature = kT0;
        c.bed_temperature = kT0;
        c.run_time = kind == Iso::Henry ? 3000.0 : 500.0;
        c.max_step = 2.0;
        c.reltol = 1e-8;
        Build(data, c);

        auto& layer = data.bed.layers_.layers_.front();
        ads_id = layer.composition_.begin()->first;
        auto& source = data.adsorbent_directory.Get(ads_id).adsorption_source_;
        if (kind == Iso::Henry) {
          auto iso = std::make_unique<HenryIsotherm>();
          iso->params_["Oxygen"].k = kHenryK;
          iso->params_["Oxygen"].Q = kHenryQ;
          source.AddIsotherm({"Oxygen"}, std::move(iso));
        } else {
          auto iso = std::make_unique<DualSiteLangmuirIsotherm>();
          auto& p = iso->params_["Oxygen"];
          p.m1 = kDslM1; p.b1 = kDslB1; p.q1 = kDslQ1;
          p.m2 = kDslM2; p.b2 = kDslB2; p.q2 = kDslQ2;
          source.AddIsotherm({"Oxygen"}, std::move(iso));
        }
        auto ldf = std::make_unique<LinearDrivingForce>();
        auto mtc = std::make_unique<MassTransferCoefficientConstant>();
        mtc->constant_mtc_ = 2.0;
        ldf->mtc_model_ = ToResistances(std::move(mtc));
        source.mass_transfer_models_["Oxygen"] = std::move(ldf);
        source.heat_of_adsorption_["Oxygen"] = stored_heat;
        source.heat_of_adsorption_mode_["Oxygen"] = mode;
      }
    };

    struct ColumnRun {
      std::vector<double> T_mid;   // mid-plane temperature at each logged time
      double T_max = 0.0;          // highest temperature in the bed over the run
      std::vector<double> dH;      // per-cell heat of Oxygen at the end of the run
      std::vector<double> q;       // per-cell Oxygen loading at the end of the run
      int first = 0, last = 0;     // bed cells
    };

    ColumnRun Integrate(Column& col) {
      ColumnRun r;
      CoupledSolver solver{col.data};
      std::string error;
      EXPECT_TRUE(col.data.PerformChecks(error)) << error;
      EXPECT_GE(solver.IntegrateCycle(col.data.cycle, false), 0) << "solver did not complete";
      const auto& log = solver.GetLog();
      const auto T = log.GetVariable(0, 0, VariableType::TEMPERATURE);
      if (T.values.empty()) return r;
      const size_t mid = T.values.front().size() / 2;
      for (const auto& row : T.values) r.T_mid.push_back(row[mid]);
      r.first = col.data.bed.GetFirstBedCellIndex();
      r.last = col.data.bed.GetLastBedCellIndex();
      for (const auto& row : T.values)
        for (int c = r.first; c <= r.last; ++c) r.T_max = std::max(r.T_max, row[c]);
      const auto& state = col.data.state;
      r.dH = state.domain_vectors.dH_ads.at(col.ads_id).at("Oxygen");
      r.q = state.ads_concentration.at(col.ads_id).at("Oxygen").var;
      return r;
    }

    double Peak(const std::vector<double>& v) { return *std::max_element(v.begin(), v.end()); }

  }  // namespace

  TEST(IsothermHeatOfAdsorption, LangmuirGivesItsEnergyAtEveryPressure) {
    LangmuirIsotherm iso;
    auto& p = iso.params_["CarbonDioxide"];
    p.m1 = 4.0; p.b1 = 1.0e-9; p.q1 = 30000.0;
    for (double P : {1.0e3, 1.0e5, 1.0e6})
      for (double T : {280.0, 320.0, 380.0})
        EXPECT_NEAR(IsostericHeatAtState(&iso, "CarbonDioxide", T, P, 0.3, 0.0, -1.0), 30000.0, 30000.0 * 1e-5)
            << "P = " << P << ", T = " << T;
    // Absent component and a slightly negative mole fraction: the dilute
    // (Henry-limit) heat at the pressure floor, finite, not the fallback.
    EXPECT_NEAR(IsostericHeatAtState(&iso, "CarbonDioxide", 300.0, 1.0e5, 0.0, 0.0, -1.0), 30000.0, 0.5);
    EXPECT_NEAR(IsostericHeatAtState(&iso, "CarbonDioxide", 300.0, 1.0e5, -1e-9, 0.0, -1.0), 30000.0, 0.5);
  }

  TEST(IsothermHeatOfAdsorption, DualSiteLangmuirFallsBetweenItsSiteEnergies) {
    DualSiteLangmuirIsotherm iso;
    auto& p = iso.params_["Oxygen"];
    p.m1 = kDslM1; p.b1 = kDslB1; p.q1 = kDslQ1;
    p.m2 = kDslM2; p.b2 = kDslB2; p.q2 = kDslQ2;
    double previous = kDslQ1 + 1.0;
    for (double y : {0.0, 1e-4, 1e-3, 1e-2, 0.1, 1.0}) {
      const double q_st = IsostericHeatAtState(&iso, "Oxygen", kT0, 1.0e6, y, 0.0, -1.0);
      EXPECT_GT(q_st, kDslQ2) << "y = " << y;
      EXPECT_LT(q_st, kDslQ1) << "y = " << y;
      EXPECT_LE(q_st, previous + 1e-6) << "the strong site fills first, so q_st falls with loading";
      previous = q_st;
    }
    // Matches IsostericHeat away from the floor.
    EXPECT_DOUBLE_EQ(IsostericHeatAtState(&iso, "Oxygen", kT0, 1.0e5, 0.5, 0.0, -1.0),
                     IsostericHeat(&iso, "Oxygen", kT0, 0.5e5, 1.0, 0.0));
  }

  TEST(IsothermHeatOfAdsorption, InertIsothermGivesTheFallback) {
    InertIsotherm iso;
    EXPECT_EQ(IsostericHeatAtState(&iso, "Oxygen", kT0, 1.0e5, 0.5, 0.0, 1234.0), 1234.0);
  }

  TEST(IsothermHeatOfAdsorption, IsothermModeReplacesTheConstantInTheEnergyBalance) {
    // Reference: the constant mode with the right heat.
    Column ref_col(Iso::Henry, HeatOfAdsorptionMode::Constant, kHenryQ);
    const ColumnRun ref = Integrate(ref_col);
    // The isotherm mode with a wrong stored constant (only the fallback).
    Column iso_col(Iso::Henry, HeatOfAdsorptionMode::Isotherm, 1000.0);
    const ColumnRun iso = Integrate(iso_col);
    // Guard: the wrong constant on its own gives almost no heating.
    Column wrong_col(Iso::Henry, HeatOfAdsorptionMode::Constant, 1000.0);
    const ColumnRun wrong = Integrate(wrong_col);
    ASSERT_FALSE(ref.T_mid.empty());
    ASSERT_FALSE(iso.T_mid.empty());
    ASSERT_FALSE(wrong.T_mid.empty());

    const double rise_ref = ref.T_max - kT0;
    const double rise_iso = iso.T_max - kT0;
    const double rise_wrong = wrong.T_max - kT0;
    std::printf("  Henry peak rise: constant Q %.4f K, isotherm %.4f K, wrong constant %.4f K\n", rise_ref,
                rise_iso, rise_wrong);
    ASSERT_GT(rise_ref, 3.0) << "no appreciable adiabatic heating";
    EXPECT_LT(rise_wrong, 0.2 * rise_ref);
    EXPECT_NEAR(rise_iso, rise_ref, 1e-3 * rise_ref);
    EXPECT_NEAR(Peak(iso.T_mid), Peak(ref.T_mid), 1e-3 * rise_ref);

    // Every bed cell carries the isotherm's heat: Q to the finite-difference
    // error (relative steps 1e-4, so a few 1e-7 relative).
    for (int c = iso.first; c <= iso.last; ++c) EXPECT_NEAR(iso.dH[c], kHenryQ, 1e-5 * kHenryQ) << "cell " << c;
  }

  TEST(IsothermHeatOfAdsorption, DualSiteHeatFollowsTheLoading) {
    Column iso_col(Iso::DualSite, HeatOfAdsorptionMode::Isotherm, 0.0);
    const ColumnRun iso = Integrate(iso_col);
    Column strong_col(Iso::DualSite, HeatOfAdsorptionMode::Constant, kDslQ1);
    const ColumnRun strong = Integrate(strong_col);
    ASSERT_FALSE(iso.T_mid.empty());
    ASSERT_FALSE(strong.T_mid.empty());

    // At the end of the feed the inlet is loaded and the outlet nearly clean.
    ASSERT_GT(iso.q[iso.first], 10.0 * std::max(iso.q[iso.last], 1e-12));
    for (int c = iso.first; c <= iso.last; ++c) {
      EXPECT_GT(iso.dH[c], kDslQ2) << "cell " << c;
      EXPECT_LT(iso.dH[c], kDslQ1) << "cell " << c;
    }
    EXPECT_LT(iso.dH[iso.first], iso.dH[iso.last] - 100.0)
        << "the loaded inlet has a lower isosteric heat than the clean outlet";
    std::printf("  Dual-site q_st: inlet %.0f J/mol, outlet %.0f J/mol\n", iso.dH[iso.first], iso.dH[iso.last]);

    // A heat below the strong-site energy releases less heat.
    std::printf("  Dual-site peak rise: isotherm %.3f K, strong-site constant %.3f K\n", iso.T_max - kT0,
                strong.T_max - kT0);
    EXPECT_GT(iso.T_max - kT0, 0.5);
    EXPECT_LT(iso.T_max, strong.T_max - 0.1);
  }

  TEST(IsothermHeatOfAdsorption, ModeIsARegistryParameterThatSurvivesSaveAndLoad) {
    Column col(Iso::Henry, HeatOfAdsorptionMode::Constant, kHenryQ);
    std::string error;
    ASSERT_TRUE(col.data.PerformChecks(error)) << error;
    const auto dir = std::filesystem::temp_directory_path() / "adsim_isotherm_heat_of_adsorption";
    std::filesystem::create_directories(dir);
    const auto path = dir / "deck.json";
    test::SaveDeckAtomically(col.data, path);

    Simulation sim = Simulation::FromFile(path.string());
    const std::string leaf =
        "adsorbent[" + std::to_string(col.ads_id) + "].heat_of_adsorption_mode[Oxygen]";
    ASSERT_NE(sim.Parameters().Find(leaf), nullptr) << leaf;
    EXPECT_EQ(sim.GetParameter(leaf), 0.0);  // constant by default
    sim.SetParameterText(leaf, "isotherm");
    EXPECT_EQ(sim.GetParameter(leaf), 1.0);

    Simulation loaded = Simulation::FromJson(sim.ToJson());
    const auto& source = loaded.Data().adsorbent_directory.Get(col.ads_id).adsorption_source_;
    EXPECT_EQ(source.HeatOfAdsorptionModeOf("Oxygen"), HeatOfAdsorptionMode::Isotherm);
    EXPECT_EQ(source.HeatOfAdsorptionModeOf("Nitrogen"), HeatOfAdsorptionMode::Constant);
    EXPECT_EQ(loaded.GetParameter(leaf), 1.0);

    // An archive from before the mode existed (AdsorptionSource version 1)
    // loads with every component in the constant mode.
    std::string old = sim.ToJson();
    const auto at = old.find("\"Adsorption_Source\"");
    ASSERT_NE(at, std::string::npos);
    const std::string v2 = "\"cereal_class_version\": 2,";
    const auto version = old.find(v2, at);
    ASSERT_NE(version, std::string::npos);
    ASSERT_LT(version, old.find("\"Mass_Transfer_Models\"", at)) << "the version must be the source's own";
    old.replace(version, v2.size(), "\"cereal_class_version\": 1,");
    Simulation from_old = Simulation::FromJson(old);
    EXPECT_EQ(from_old.Data().adsorbent_directory.Get(col.ads_id).adsorption_source_.HeatOfAdsorptionModeOf("Oxygen"),
              HeatOfAdsorptionMode::Constant);
    EXPECT_EQ(from_old.GetParameter(leaf), 0.0);
  }

}  // namespace adsim
