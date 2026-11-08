/*
 * Skarstrom - pressure and temperature swing adsorption simulator
 * Copyright (c) 2024-2026 Sam Affleck and contributors
 *
 * Licensed under the terms in the LICENSE file at the repository root.
 */

// Per-adsorbent series of mixed layers.
//
// A layer holding several solids logs each solid's own loading q_k [mol per kg
// of solid k] beside the aggregate sum_k w_k q_k. Checked on solved runs
// through the Simulation API (Series, Profile, SeriesAdsorbents, SaveLogs):
//
//  (a) a 50/50 layer of A and an identical copy A': each solid's loading,
//      equilibrium loading and mass-transfer coefficient equal the aggregate;
//  (b) A mixed with an inert solid (w_A = 0.6): A's loading equals that of a
//      pure-A layer at A's bulk density in the blend (particle density
//      w_A rho, heat capacity cp / w_A so the solid heat capacity per bed
//      volume is the blend's), the inert holds 0, and the aggregate is the
//      mass-weighted sum w_A q_A;
//  (c) a single-solid deck logs no per-adsorbent series; its per-adsorbent
//      series are read off the aggregate and equal it.
//
// Plus the series manifest, which lists the per-adsorbent variables.

#define _USE_MATH_DEFINES
#include <cmath>

#include "../src/pch.h"
#include "adsim/boundary_conditions/boundary_conditions.h"
#include "adsim/boundary_conditions/boundary_stream.h"
#include "adsim/boundary_conditions/step.h"
#include "adsim/equations/dispersion_model.h"
#include "adsim/equations/heat_transfer_model.h"
#include "adsim/equations/isotherm.h"
#include "adsim/objects/manifests.h"
#include "adsim/simulation.h"

#include <cereal/external/rapidjson/document.h>

#include <algorithm>
#include <filesystem>
#include <set>
#include <string>
#include <vector>

#include <gtest/gtest.h>

#include "../fixtures/user_data_test.h"

namespace adsim {

  namespace {

    double MaxAbs(const std::vector<double>& v) {
      double m = 0.0;
      for (double x : v) m = std::max(m, std::abs(x));
      return m;
    }

    double MaxAbsDiff(const std::vector<double>& a, const std::vector<double>& b, double scale = 1.0) {
      EXPECT_EQ(a.size(), b.size());
      double m = 0.0;
      for (size_t i = 0; i < std::min(a.size(), b.size()); ++i) m = std::max(m, std::abs(a[i] - scale * b[i]));
      return m;
    }

    // The last row (end of the step) of a profile.
    std::vector<double> LastRow(const AxialProfile& p) {
      if (p.rows == 0) return {};
      return std::vector<double>(p.values.end() - static_cast<std::ptrdiff_t>(p.cols), p.values.end());
    }

  }  // namespace

  class PerAdsorbentSeriesTest : public UserDataTest {
  public:
    const std::string CO2 = "CarbonDioxide";
    const std::string N2 = "Nitrogen";

    static constexpr double kRho = 1100.0;  // particle density of A [kg/m3]
    static constexpr double kCp = 900.0;    // [J/(kg.K)]

    enum class Packing { Single, Copy, Inert, PureScaled };
    static constexpr double kInertW = 0.6;  // w_A in the A + inert blend

    // A: CO2 on a dual-site Langmuir isotherm with constant-MTC LDF.
    void ConfigureA(AdsorbentID id, double rho, double cp) {
      auto& ads = GetAdsorbent(id);
      ads.solid_density_ = rho;
      ads.heat_capacity_ = cp;
      ads.thermal_conductivity_ = 0.5;
      ads.ep_ = 0.35;
      auto dsl = std::make_unique<DualSiteLangmuirIsotherm_Concentration>(CO2);
      dsl->params_[CO2].m1 = 2.5;  dsl->params_[CO2].b1 = 4.0e-7; dsl->params_[CO2].q1 = 30000;
      dsl->params_[CO2].m2 = 3.0;  dsl->params_[CO2].b2 = 1.0e-7; dsl->params_[CO2].q2 = 25000;
      AddIsotherm({CO2}, id, std::move(dsl));
      AddLinearDrivingForceModel(id, CO2, 0.2);
      AddLinearDrivingForceModel(id, N2, 0.0);
      SetHOA(id, CO2, 30000.0);
      SetHOA(id, N2, 0.0);
    }

    // An inert solid with A's density, voidage, heat capacity and conductivity.
    void ConfigureInert(AdsorbentID id) {
      auto& ads = GetAdsorbent(id);
      ads.solid_density_ = kRho;
      ads.heat_capacity_ = kCp;
      ads.thermal_conductivity_ = 0.5;
      ads.ep_ = 0.35;
      for (const auto& c : {CO2, N2}) ads.adsorption_source_.SetIsothermType(c, IsothermType::INERT);
    }

    // A one-layer breakthrough column (20 cells, 100 s feed).
    UserData Build(Packing packing) {
      data = UserData{};
      UserDataTest::SetUp();
      AddComponents({CO2, N2});
      auto [layer_id, ads_id] = AddLayer("Sorbent A");
      auto& layer = GetLayer(layer_id);
      layer.domain_.length = 0.3;
      layer.particle_diameter_ = 0.0015;
      layer.eb_ = 0.4;
      layer.domain_.n_cells = 20;
      std::get<CylinderBed>(data.bed.shape_).diameter = 0.02;
      if (packing == Packing::PureScaled)
        ConfigureA(ads_id, kInertW * kRho, kCp / kInertW);
      else
        ConfigureA(ads_id, kRho, kCp);

      if (packing == Packing::Copy) {
        const auto copy = data.adsorbent_directory.Add("Sorbent A copy");
        ConfigureA(copy, kRho, kCp);
        layer.composition_.clear();
        layer.composition_[ads_id] = 0.5;
        layer.composition_[copy] = 0.5;
      } else if (packing == Packing::Inert) {
        const auto inert = data.adsorbent_directory.Add("Glass beads");
        ConfigureInert(inert);
        layer.composition_.clear();
        layer.composition_[ads_id] = kInertW;
        layer.composition_[inert] = 1.0 - kInertW;
      }

      layer.equations_.dispersion_model_ = CreateAxialDispersionModel(AxialDispersionType::RUTHVEN);
      layer.equations_.heat_dispersion_model_ = CreateAxialHeatDispersionModel(AxialHeatDispersionType::CONSTANT);
      AddConstantInternalHeatTransferCoefficientModel(8.0, layer_id);
      AddConstantExternalHeatTransferCoefficientModel(10.0, layer_id);
      data.bed.ambient_temperature_ = 293.0;
      AddInitialConditions(1.0e5, 293.0, 293.0, {{CO2, 0.0}, {N2, 1.0}});

      AddStep("Feed", 100.0);
      auto inlet = CreateNoValveInletStream(CreateVelocityBC(0.05), {{CO2, 0.2}, {N2, 0.8}}, 298.0);
      inlet->name_ = "Feed inlet";
      auto in_bc = std::make_unique<FlowBoundaryCondition>();
      in_bc->direction_ = BoundaryDirection::INLET;
      in_bc->location_ = BoundaryLocation::BOTTOM;
      in_bc->streams_.push_back(std::move(inlet));
      data.cycle.GetStep("Feed").bottom_boundary_condition_ = std::move(in_bc);
      auto out_bc = std::make_unique<PressureBoundaryCondition>();
      out_bc->pressure_ = CreateConstantBoundaryValue(1.0e5);
      out_bc->direction_ = BoundaryDirection::OUTLET;
      out_bc->location_ = BoundaryLocation::TOP;
      out_bc->mole_fraction_[CO2] = 0.2;
      out_bc->mole_fraction_[N2] = 0.8;
      data.cycle.GetStep("Feed").top_boundary_condition_ = std::move(out_bc);
      ConfigureCycle(1);
      return std::move(data);
    }

    void Run(Simulation& sim, Packing packing) {
      sim.SetData(Build(packing));
      const RunResult r = sim.Run();
      ASSERT_TRUE(r.Ok()) << r.message;
    }
  };

  // (a) Each of two identical solids holds what the aggregate says.
  TEST_F(PerAdsorbentSeriesTest, IdenticalCopiesEachEqualTheAggregate) {
    Simulation sim;
    ASSERT_NO_FATAL_FAILURE(Run(sim, Packing::Copy));
    const auto ads = sim.SeriesAdsorbents();
    ASSERT_EQ(ads.size(), 2u);
    EXPECT_TRUE(ads[0].logged);
    EXPECT_EQ(ads[0].label, "Sorbent A");
    EXPECT_EQ(ads[1].label, "Sorbent A copy");
    EXPECT_TRUE(sim.Log().LogsPerAdsorbent());

    for (const std::string var : {"adsorbed_concentration", "equilibrium_adsorbed_concentration",
                                  "mass_transfer_coefficient"}) {
      const AxialProfile agg = sim.Profile(var, CO2);
      const double scale = std::max(1e-12, MaxAbs(agg.values));
      for (const auto& a : ads) {
        const AxialProfile p = sim.Profile(var, CO2, -1, 0, a.label);
        EXPECT_EQ(p.adsorbent, a.label);
        ASSERT_EQ(p.values.size(), agg.values.size());
        EXPECT_LT(MaxAbsDiff(p.values, agg.values) / scale, 1e-9) << var << " of " << a.label;
      }
      // Series at the bed ends and inside, by label and by deck ID.
      for (auto where : {SeriesLocation::Bottom, SeriesLocation::Top, SeriesLocation::Position}) {
        const TimeSeries s_agg = sim.Series(var, CO2, where, 0.1);
        const TimeSeries s_a = sim.Series(var, CO2, where, 0.1, -1, std::to_string(ads[1].id));
        EXPECT_EQ(s_a.adsorbent, "Sorbent A copy");
        EXPECT_LT(MaxAbsDiff(s_a.values, s_agg.values) / scale, 1e-9) << var;
      }
    }
    EXPECT_GT(MaxAbs(sim.Profile("adsorbed_concentration", CO2, -1, 0, "Sorbent A").values), 0.1);

    // The log names, and the CSV files SaveLogs writes.
    const VariableLogEntry e =
        sim.Log().GetAdsorbentComponentVariable(0, 0, ComponentVariableType::ADSORBED_CONCENTRATION, CO2, ads[0].id);
    EXPECT_EQ(e.variable_name, "CarbonDioxide Adsorbed Concentration (Sorbent A)");
    EXPECT_EQ(e.units, "[mol/kg]");
    const auto dir = std::filesystem::current_path() / "Results" / "PerAdsorbentSeries" / "copy";
    std::filesystem::remove_all(dir);
    sim.SaveLogs(dir.string());
    EXPECT_TRUE(std::filesystem::exists(dir / "CarbonDioxide_Adsorbed_Concentration.csv"));
    EXPECT_TRUE(std::filesystem::exists(dir / "CarbonDioxide_Adsorbed_Concentration__Sorbent_A_.csv"));
    EXPECT_TRUE(std::filesystem::exists(dir / "CarbonDioxide_Adsorbed_Concentration__Sorbent_A_copy_.csv"));
    EXPECT_TRUE(std::filesystem::exists(dir / "Nitrogen_Equilibrium_Adsorbed_Concentration__Sorbent_A_copy_.csv"));
    EXPECT_TRUE(std::filesystem::exists(dir / "CarbonDioxide_Mass_Transfer_Coefficient__Sorbent_A_.csv"));
  }

  // (b) A diluted with an inert solid holds what pure A at its bulk density
  // holds; the inert holds nothing; the aggregate is w_A q_A.
  TEST_F(PerAdsorbentSeriesTest, InertDiluentHoldsNothingAndTheSorbentMatchesPureA) {
    Simulation blend;
    ASSERT_NO_FATAL_FAILURE(Run(blend, Packing::Inert));
    const auto ads = blend.SeriesAdsorbents();
    ASSERT_EQ(ads.size(), 2u);
    EXPECT_EQ(ads[1].label, "Glass beads");

    const AxialProfile q_a = blend.Profile("adsorbed_concentration", CO2, -1, 0, "Sorbent A");
    const AxialProfile q_i = blend.Profile("adsorbed_concentration", CO2, -1, 0, "Glass beads");
    const AxialProfile q_agg = blend.Profile("adsorbed_concentration", CO2);
    const double scale = MaxAbs(q_a.values);
    ASSERT_GT(scale, 0.1);
    EXPECT_EQ(MaxAbs(q_i.values), 0.0);
    EXPECT_EQ(MaxAbs(blend.Profile("equilibrium_adsorbed_concentration", CO2, -1, 0, "Glass beads").values), 0.0);
    // Same run, same cells: the aggregate is w_A q_A + w_I * 0.
    EXPECT_LT(MaxAbsDiff(q_agg.values, q_a.values, kInertW) / scale, 1e-12);

    Simulation pure;
    ASSERT_NO_FATAL_FAILURE(Run(pure, Packing::PureScaled));
    EXPECT_FALSE(pure.Log().LogsPerAdsorbent());
    const AxialProfile q_pure = pure.Profile("adsorbed_concentration", CO2);
    // Different runs (other unknowns, other time steps): compare the
    // end-of-step profiles to the solver tolerance.
    const auto last_a = LastRow(q_a);
    const auto last_pure = LastRow(q_pure);
    ASSERT_FALSE(last_a.empty());
    EXPECT_LT(MaxAbsDiff(last_a, last_pure) / scale, 2e-3);
    const auto qe_a = LastRow(blend.Profile("equilibrium_adsorbed_concentration", CO2, -1, 0, "Sorbent A"));
    const auto qe_pure = LastRow(pure.Profile("equilibrium_adsorbed_concentration", CO2));
    EXPECT_LT(MaxAbsDiff(qe_a, qe_pure) / scale, 2e-3);
  }

  // (c) A single-solid deck logs no per-adsorbent series; asking for one
  // reads the aggregate, which is that solid's loading.
  TEST_F(PerAdsorbentSeriesTest, SingleSolidLogsOnlyTheAggregate) {
    Simulation sim;
    ASSERT_NO_FATAL_FAILURE(Run(sim, Packing::Single));
    EXPECT_FALSE(sim.Log().LogsPerAdsorbent());
    for (const auto& [cycle, steps] : sim.Log().GetSimulationLog())
      for (const auto& [step, log] : steps) EXPECT_TRUE(log.adsorbent_variables.empty());
    const auto ads = sim.SeriesAdsorbents();
    ASSERT_EQ(ads.size(), 1u);
    EXPECT_FALSE(ads[0].logged);

    for (const std::string var : {"adsorbed_concentration", "equilibrium_adsorbed_concentration",
                                  "mass_transfer_coefficient"}) {
      const AxialProfile agg = sim.Profile(var, CO2);
      const AxialProfile one = sim.Profile(var, CO2, -1, 0, "Sorbent A");
      EXPECT_EQ(one.values, agg.values) << var;
      const TimeSeries s_agg = sim.Series(var, CO2, SeriesLocation::Bottom);
      const TimeSeries s_one = sim.Series(var, CO2, SeriesLocation::Bottom, 0.0, -1, "Sorbent A");
      EXPECT_EQ(s_one.values, s_agg.values) << var;
    }

    const auto dir = std::filesystem::current_path() / "Results" / "PerAdsorbentSeries" / "single";
    std::filesystem::remove_all(dir);
    sim.SaveLogs(dir.string());
    EXPECT_TRUE(std::filesystem::exists(dir / "CarbonDioxide_Adsorbed_Concentration.csv"));
    EXPECT_FALSE(std::filesystem::exists(dir / "CarbonDioxide_Adsorbed_Concentration__Sorbent_A_.csv"));

    // Arguments that name no per-adsorbent series are refused.
    EXPECT_THROW(sim.Series("pressure", "", SeriesLocation::Top, 0.0, -1, "Sorbent A"), std::invalid_argument);
    EXPECT_THROW(sim.Series("mole_fraction", CO2, SeriesLocation::Top, 0.0, -1, "Sorbent A"), std::invalid_argument);
    EXPECT_THROW(sim.Profile("adsorbed_concentration", CO2, -1, 0, "No such solid"), std::invalid_argument);
    EXPECT_THROW(sim.Profile("adsorbed_concentration", CO2, -1, 0, "999"), std::invalid_argument);
  }

  // The series manifest describes BulkVariableNames() and
  // ComponentVariableNames(), and marks the per-adsorbent ones.
  TEST(SeriesManifest, ListsEverySeriesAndThePerAdsorbentOnes) {
    rapidjson::Document doc;
    doc.Parse(SeriesManifest().c_str());
    ASSERT_FALSE(doc.HasParseError());
    ASSERT_TRUE(doc.HasMember("variables"));
    const auto& vars = doc["variables"];
    std::vector<std::string> keys, per_ads;
    for (const auto& v : vars.GetArray()) {
      keys.emplace_back(v["key"].GetString());
      if (v["per_adsorbent"].GetBool()) {
        per_ads.emplace_back(v["key"].GetString());
        EXPECT_STREQ(v["kind"].GetString(), "component");
      }
    }
    std::vector<std::string> expected = BulkVariableNames();
    for (const auto& c : ComponentVariableNames()) expected.push_back(c);
    EXPECT_EQ(keys, expected);
    EXPECT_EQ(per_ads, PerAdsorbentVariableNames());
    EXPECT_EQ(per_ads.size(), 4u);
    for (const auto& v : vars.GetArray()) {
      if (std::string(v["key"].GetString()) != "adsorbed_concentration") continue;
      EXPECT_STREQ(v["name"].GetString(), "{component} Adsorbed Concentration");
      EXPECT_STREQ(v["adsorbent_name"].GetString(), "{component} Adsorbed Concentration ({adsorbent})");
      EXPECT_STREQ(v["unit"].GetString(), "mol/kg");
      EXPECT_TRUE(v["porous_only"].GetBool());
    }
  }

}  // namespace adsim
