// Illustrative TSA example deck (examples/decks/Illustrative_TSA_13X.json).
//
// A four-step, indirectly heated temperature swing adsorption cycle for CO2
// capture from a 15 % CO2 / 85 % N2 flue gas: adsorption, heating with an
// electric energy source (CO2 product leaves the feed end), cooling with a
// negative energy source (column closed) and feed pressurisation.
//
// ILLUSTRATIVE ONLY. This is not a reproduction of any published experiment or
// simulation. The adsorbent is the library entry "Zeolite 13X (Perez 2019)"
// (data/adsorbents/, Estupinan Perez et al., Sep. Purif. Technol. 224
// (2019) 553-563), loaded unchanged, and the column geometry, wall and heat
// transfer coefficients are those of the Perez2019 decks. The operating
// conditions (feed velocity, step times, heater and cooler duties) are round
// numbers chosen for the docs tutorial (docs/tutorials/tsa-regeneration.md), and the
// heating step takes the isotherm above its 298-348 K fit range.
//
// The test builds the deck, saves it as an example (CLAUDE.md: re-saved decks
// are committed), runs it and checks only that it solves, conserves mass and
// behaves like a TSA cycle (CO2 is enriched in the heating-step product).
//
// The deck uses the Ergun momentum model, the default for new layers, so it is
// what a new user gets. It doubles as the regression test for Ergun at abrupt
// step switches (feed -> closed / vent, closed -> feed): before the boundary
// faces of pressure and closed conditions imposed the explicit Ergun root
// (ImposeBoundaryFaceVelocity in boundary_conditions.cpp) every Heating step,
// and some Adsorption steps, began with an IDA corrector failure rescued only
// by a restart. Both momentum models must now reach CSS with no step retry.

#define _USE_MATH_DEFINES
#include <cmath>
#include <filesystem>

#include "../src/pch.h"
#include "adsim/solver/coupled_solver.h"
#include "adsim/solver/kpis.h"
#include "adsim/boundary_conditions/boundary_conditions.h"
#include "adsim/boundary_conditions/boundary_stream.h"
#include "adsim/boundary_conditions/step.h"
#include "adsim/equations/dispersion_model.h"
#include "adsim/equations/heat_transfer_model.h"
#include "adsim/objects/wall.h"

#include <gtest/gtest.h>

#include "../fixtures/user_data_test.h"

namespace adsim {

  namespace {
    // Operating point of the illustrative cycle (round numbers, see header).
    constexpr double kFeedVelocity = 0.05;   // [m/s] superficial, adsorption step
    constexpr double kFeedTemperature = 298; // [K]
    constexpr double kPressure = 1.0e5;      // [Pa]
    constexpr double kAdsorptionTime = 600;  // [s]
    constexpr double kHeatingTime = 900;     // [s]
    constexpr double kCoolingTime = 900;     // [s]
    constexpr double kPressurisationTime = 30; // [s]
    constexpr double kHeaterPower = 30.0;    // [W] electric heating of the layer
    constexpr double kCoolerPower = -12.0;   // [W] heat removal while closed
    constexpr int kMaxCycles = 60;
  }

  class IllustrativeTSATest : public UserDataTest {
  protected:
    // Builds the deck with `momentum` and runs it to CSS. Saves the deck as
    // the example when `save_deck` is set.
    void Run(MomentumModel momentum, bool save_deck);
  };

  TEST_F(IllustrativeTSATest, IllustrativeTSA13X) {
    Run(MomentumModel::Ergun, true);
  }

  // The same cycle with the legacy Darcy / Kozeny-Carman pressure drop (not
  // saved): both models must run it cleanly, and at these velocities the
  // inertial term is small, so the KPIs barely move (tsa-regeneration.mdx).
  TEST_F(IllustrativeTSATest, IllustrativeTSA13X_KozenyCarman) {
    Run(MomentumModel::KarmenCozeny, false);
  }

  void IllustrativeTSATest::Run(MomentumModel momentum, bool save_deck) {
    data.info.title = "Illustrative TSA - CO2/N2 on 13X";
    data.info.description =
      "ILLUSTRATIVE, not a validated reproduction. Four-step indirectly heated TSA "
      "(adsorption, electric heating with CO2 product at the feed end, cooling, feed "
      "pressurisation) for a 15 % CO2 / 85 % N2 flue gas. Adsorbent: library entry "
      "Zeolite 13X (Perez 2019); column, wall and heat transfer from the Perez2019 decks. "
      "Operating conditions are round numbers for the docs tutorial; heating goes above "
      "the isotherm's 298-348 K fit range.";
    data.info.author = "Skarstrom docs tutorial";

    const std::string CO2 = "CarbonDioxide";
    const std::string N2 = "Nitrogen";
    AddComponents({ CO2, N2 });

    // Adsorbent: the curated library entry, loaded exactly as the library
    // browser does (cereal blob, then fill any missing components).
    const AdsorbentID ads_id = data.adsorbent_directory.Add("Zeolite 13X (Perez 2019)");
    auto& ads = data.adsorbent_directory.Get(ads_id);
    LoadAdsorbentFromJSON("../data/adsorbents/zeolite_13x_perez2019.json", ads);
    ads.adsorption_source_.Initialise(data.fluid);

    // Column: Perez2019 geometry (0.30 m x 3.2 cm, eps_b = 0.41, d_p = 1.5 mm).
    const LayerID layer_id = data.bed.layers_.Add("13X");
    auto& layer = data.bed.layers_.Get(layer_id);
    layer.domain_.length = 0.3;
    layer.domain_.n_cells = 30;
    layer.eb_ = 0.41;
    layer.particle_diameter_ = 1.5e-3;
    layer.composition_[ads_id] = 1.0;
    // Ergun (the default for a new layer) for the deck; the Perez2019 decks
    // themselves use Darcy (Kozeny-Carman), which the second test runs.
    layer.momentum_model_ = momentum;
    std::get<CylinderBed>(data.bed.shape_).diameter = 0.032;

    data.bed.wall_.ApplyMaterial(WallMaterial::SS316L);
    data.bed.wall_.thickness_ = 3.1e-3;
    data.bed.ambient_temperature_ = 298;
    AddConstantInternalHeatTransferCoefficientModel(8.6, layer_id);
    AddConstantExternalHeatTransferCoefficientModel(10, layer_id);
    layer.equations_.dispersion_model_ = CreateAxialDispersionModel(AxialDispersionType::RUTHVEN);
    auto kz = std::make_unique<HeatDispersionModelConstant>();
    kz->dispersion_coefficient_ = 0.0903; // [W/(m K)]
    layer.equations_.heat_dispersion_model_ = std::move(kz);
    data.fluid.adsorbed_heat_capacity_[CO2] = 30.7;
    data.fluid.adsorbed_heat_capacity_[N2] = 30.7;

    AddInitialConditions(kPressure, kFeedTemperature, kFeedTemperature, { {CO2, 0.0}, {N2, 1.0} });

    const std::unordered_map<std::string, double> feed = { {CO2, 0.15}, {N2, 0.85} };

    // 1. Adsorption: flue gas in at the bottom, N2-rich raffinate out of the top.
    AddStep("Adsorption", kAdsorptionTime);
    {
      std::vector<std::unique_ptr<IBoundaryStream>> streams;
      auto s = CreateNoValveInletStream(CreateVelocityBC(kFeedVelocity), feed, kFeedTemperature);
      s->name_ = "Flue gas";
      streams.push_back(std::move(s));
      AddFlowBC(std::move(streams), BoundaryLocation::BOTTOM, BoundaryDirection::INLET, "Adsorption");
    }
    AddPressureOutletBC(CreateConstantBoundaryValue(kPressure), BoundaryLocation::TOP, "Adsorption",
                        { {CO2, 0.0}, {N2, 1.0} });

    // 2. Heating: top closed, CO2 leaves the feed end at 1 bar while an electric
    //    heater puts kHeaterPower into the layer.
    AddStep("Heating", kHeatingTime);
    AddClosedBC(BoundaryLocation::TOP, "Heating");
    AddPressureOutletBC(CreateConstantBoundaryValue(kPressure), BoundaryLocation::BOTTOM, "Heating",
                        { {CO2, 1.0}, {N2, 0.0} });

    // 3. Cooling: both ends closed, kCoolerPower removed from the layer.
    AddStep("Cooling", kCoolingTime);
    AddClosedBC(BoundaryLocation::TOP, "Cooling");
    AddClosedBC(BoundaryLocation::BOTTOM, "Cooling");

    // 4. Pressurisation: feed back to 1 bar through the bottom, top closed.
    AddStep("Pressurisation", kPressurisationTime);
    AddClosedBC(BoundaryLocation::TOP, "Pressurisation");
    AddPressureInletBC(CreateLinearBoundaryValue(kPressure, kPressurisationTime), feed, kFeedTemperature,
                       BoundaryLocation::BOTTOM, "Pressurisation");

    EnergySource heater;
    heater.enabled_ = true;
    heater.power_ = kHeaterPower;
    data.cycle.GetStep("Heating").energy_sources_[layer_id] = heater;
    EnergySource cooler;
    cooler.enabled_ = true;
    cooler.power_ = kCoolerPower;
    data.cycle.GetStep("Cooling").energy_sources_[layer_id] = cooler;

    ConfigureCycle(kMaxCycles);
    data.settings.css.method_ = CSSMethod::MASS_BALANCE_COMPONENT;
    data.settings.css.mass_balance_tolerance_ = 0.005;
    data.settings.css.min_cycles_ = 5;
    // CO2 KPIs: fed in adsorption + pressurisation (bottom), product = heating (bottom).
    ConfigureKPI(0, {0, 3}, {BoundaryLocation::BOTTOM, BoundaryLocation::BOTTOM},
                 {1}, {BoundaryLocation::BOTTOM});

    std::string error;
    ASSERT_TRUE(data.PerformChecks(error)) << error;
    if (save_deck)
      SaveToJson(data, (std::filesystem::current_path() / "../examples/decks/Illustrative_TSA_13X.json").string());

    CoupledSolver solver{data};
    // Count step retries (each one is a recoverable IDA failure, logged at
    // MEDIUM before the step is re-attempted).
    int step_retries = 0;
    solver.SetLogLevel(LogLevel::MEDIUM);
    solver.SetMessageCallback([&step_retries](LogLevel, const std::string& msg) {
      if (msg.find("IDA failure") != std::string::npos) {
        ++step_retries;
        std::cout << "[TSA]" << msg << std::endl;
      }
    });
    const int retval = solver.IntegrateCycle(data.cycle, false);
    ASSERT_GE(retval, 0) << "Solver did not complete";
    EXPECT_TRUE(solver.ReachedCSS()) << "CSS not reached in " << kMaxCycles << " cycles";
    EXPECT_EQ(step_retries, 0) << "step switches needed solver retries";

    const EquationLog& eq_log = solver.GetLog();
    EXPECT_TRUE(IsMassConserved(eq_log));
    EXPECT_TRUE(IsMolesConserved(eq_log));
    EXPECT_TRUE(IsComponentMolesConserved(eq_log));

    const auto& cycles = eq_log.sys_data.cycle_index_set_;
    ASSERT_FALSE(cycles.empty());
    KPIParams kpi(eq_log, solver.GetBoundaryLog());
    const auto* cfg = data.kpi_config.FirstTracked();
    ASSERT_NE(cfg, nullptr);
    kpi.cycle_index_ = *cycles.rbegin();
    kpi.product_comp_ = cfg->product_component;
    kpi.feed_endpoints = cfg->feed_endpoints;
    kpi.product_endpoints = cfg->product_endpoints;
    const double purity = CalculatePurity(kpi);
    const double recovery = CalculateRecovery(kpi);
    const double production_rate = CalculateProductionRate(kpi);
    const double cycle_time = kAdsorptionTime + kHeatingTime + kCoolingTime + kPressurisationTime;
    const double heat_per_kg = kHeaterPower * kHeatingTime
                             / (production_rate * cycle_time * data.fluid.molar_mass_.at(CO2));

    double t_max = 0.0;
    double t_min = 1e9;
    for (int step = 0; step < 4; ++step) {
      const auto entry = eq_log.GetVariable(*cycles.rbegin(), step, VariableType::TEMPERATURE);
      for (const auto& profile : entry.values)
        for (size_t i = 1; i + 1 < profile.size(); ++i) {
          t_max = std::max(t_max, profile[i]);
          t_min = std::min(t_min, profile[i]);
        }
    }
    std::cout << "[TSA] " << MomentumModelList[static_cast<int>(momentum)]
              << " | cycles run = " << (*cycles.rbegin() + 1)
              << " | CO2 purity = " << purity * 100 << " % | CO2 recovery = " << recovery * 100
              << " % | production rate = " << production_rate << " mol/s"
              << " | heater energy per kg CO2 = " << heat_per_kg / 1e6 << " MJ/kg"
              << " | bed T range in last cycle = [" << t_min << ", " << t_max << "] K" << std::endl;

    // TSA behaviour, not a validation: the heating-step product is enriched in CO2
    // well above the 15 % feed, recovery is a fraction, and the bed is heated.
    EXPECT_GT(purity, 0.5);
    EXPECT_GT(recovery, 0.0);
    EXPECT_LE(recovery, 1.0);
    EXPECT_GT(t_max, kFeedTemperature + 30.0);
  }

} // namespace adsim
