// adsim includes
#include "adsim/objects/reactor_equations.h"
#include "adsim/objects/porous_media.h"
#include "adsim/objects/reactor.h"
#include "adsim/objects/domain.h"
#include "adsim/objects/fluid.h"
#include "adsim/boundary_conditions/step.h"
#include "adsim/objects/user_data.h"
#include "adsim/equations/adsorbent_equations.h"
#include "adsim/solver/coupled_solver.h"

// STL includes
#include <algorithm>
#include <cmath>
#include <tuple>
#include <vector>

// Gtest includes
#include <gtest/gtest.h>


namespace adsim {

  TEST(EnergySourceTest, DefaultState) {
    EnergySource src;
    EXPECT_FALSE(src.enabled_);
    EXPECT_DOUBLE_EQ(src.power_, 0.0);
  }

  // Helper to create a simple bed with one layer
  struct SimpleTestBed {
    AdsorbentDirectory dir;
    Bed bed;
    Fluid fluid;
    LayerID layer_id;

    SimpleTestBed(double diameter = 0.1, double length = 1.0, int n_cells = 10, double eb = 0.4) {
      AdsorbentID ads_id = dir.Add("TestAdsorbent");
      dir.Get(ads_id).solid_density_ = 1000.0;

      std::get<CylinderBed>(bed.shape_).diameter = diameter;
      layer_id = bed.layers_.Add("Layer 1");
      auto& layer = bed.layers_.Get(layer_id);
      layer.momentum_model_ = MomentumModel::KarmenCozeny;  // fixture predates Ergun: keep its Darcy / Kozeny-Carman pressure drop
      layer.domain_.length = length;
      layer.domain_.n_cells = n_cells;
      layer.eb_ = eb;
      layer.particle_diameter_ = 0.002;
      layer.AddAdsorbent(ads_id, 1.0);
      dir.Get(ads_id).ep_ = 0.3;

      fluid.AddComponent("N2");
      bed.Initialise(dir, fluid);
    }
  };

  // Helper to create a bed with two layers
  struct TwoLayerTestBed {
    AdsorbentDirectory dir;
    Bed bed;
    Fluid fluid;
    LayerID layer_id_1;
    LayerID layer_id_2;

    TwoLayerTestBed() {
      AdsorbentID ads_id = dir.Add("TestAdsorbent");
      dir.Get(ads_id).solid_density_ = 1000.0;

      std::get<CylinderBed>(bed.shape_).diameter = 0.1;

      layer_id_1 = bed.layers_.Add("Layer 1");
      auto& layer1 = bed.layers_.Get(layer_id_1);
      layer1.momentum_model_ = MomentumModel::KarmenCozeny;  // fixture predates Ergun: keep its Darcy / Kozeny-Carman pressure drop
      layer1.domain_.length = 1.0;
      layer1.domain_.n_cells = 10;
      layer1.eb_ = 0.4;
      layer1.particle_diameter_ = 0.002;
      layer1.AddAdsorbent(ads_id, 1.0);

      layer_id_2 = bed.layers_.Add("Layer 2");
      auto& layer2 = bed.layers_.Get(layer_id_2);
      layer2.momentum_model_ = MomentumModel::KarmenCozeny;  // fixture predates Ergun: keep its Darcy / Kozeny-Carman pressure drop
      layer2.domain_.length = 0.5;
      layer2.domain_.n_cells = 5;
      layer2.eb_ = 0.4;
      layer2.particle_diameter_ = 0.002;
      layer2.AddAdsorbent(ads_id, 1.0);

      dir.Get(ads_id).ep_ = 0.3;

      fluid.AddComponent("N2");
      bed.Initialise(dir, fluid);
    }
  };

  TEST(EnergySourceTest, PopulateWithEmptyMap) {
    SimpleTestBed setup;
    DomainVectors dv;
    dv.Initialise(setup.bed, setup.dir, setup.fluid);

    std::map<size_t, EnergySource> empty_sources;
    dv.PopulateEnergySources(empty_sources, setup.bed);

    for (double val : dv.energy_source_power) {
      EXPECT_DOUBLE_EQ(val, 0.0);
    }
  }

  TEST(EnergySourceTest, PopulateSingleLayerWithPower) {
    SimpleTestBed setup(0.1, 1.0, 10, 0.4);
    DomainVectors dv;
    dv.Initialise(setup.bed, setup.dir, setup.fluid);

    EnergySource src;
    src.enabled_ = true;
    src.power_ = 100.0; // 100 W

    std::map<size_t, EnergySource> sources;
    sources[setup.layer_id] = src;
    dv.PopulateEnergySources(sources, setup.bed);

    // Calculate expected volumetric power density
    double area = setup.bed.GetArea();
    auto& layer = setup.bed.layers_.Get(setup.layer_id);
    double cell_volume = layer.domain_.GetCellVolume(area);
    double V_layer = cell_volume * layer.domain_.n_cells;
    double expected_q = 100.0 / V_layer;

    // Interior cells should have the source
    for (int i = layer.domain_.start_index; i <= layer.domain_.end_index; ++i) {
      EXPECT_DOUBLE_EQ(dv.energy_source_power[i], expected_q);
    }
  }

  TEST(EnergySourceTest, PopulateTwoLayersOneActive) {
    TwoLayerTestBed setup;
    DomainVectors dv;
    dv.Initialise(setup.bed, setup.dir, setup.fluid);

    EnergySource src;
    src.enabled_ = true;
    src.power_ = 50.0;

    std::map<size_t, EnergySource> sources;
    sources[setup.layer_id_2] = src; // Only layer 2 has power

    dv.PopulateEnergySources(sources, setup.bed);

    // Layer 1 cells should be zero
    auto& layer1 = setup.bed.layers_.Get(setup.layer_id_1);
    for (int i = layer1.domain_.start_index; i <= layer1.domain_.end_index; ++i) {
      EXPECT_DOUBLE_EQ(dv.energy_source_power[i], 0.0);
    }

    // Layer 2 cells should have the source
    auto& layer2 = setup.bed.layers_.Get(setup.layer_id_2);
    double area = setup.bed.GetArea();
    double cell_volume = layer2.domain_.GetCellVolume(area);
    double V_layer = cell_volume * layer2.domain_.n_cells;
    double expected_q = 50.0 / V_layer;

    for (int i = layer2.domain_.start_index; i <= layer2.domain_.end_index; ++i) {
      EXPECT_DOUBLE_EQ(dv.energy_source_power[i], expected_q);
    }
  }

  TEST(EnergySourceTest, DisabledSourceIsZero) {
    SimpleTestBed setup;
    DomainVectors dv;
    dv.Initialise(setup.bed, setup.dir, setup.fluid);

    EnergySource src;
    src.enabled_ = false;
    src.power_ = 100.0; // Non-zero power but disabled

    std::map<size_t, EnergySource> sources;
    sources[setup.layer_id] = src;
    dv.PopulateEnergySources(sources, setup.bed);

    for (double val : dv.energy_source_power) {
      EXPECT_DOUBLE_EQ(val, 0.0);
    }
  }

  TEST(EnergySourceTest, NegativePowerForCooling) {
    SimpleTestBed setup;
    DomainVectors dv;
    dv.Initialise(setup.bed, setup.dir, setup.fluid);

    EnergySource src;
    src.enabled_ = true;
    src.power_ = -200.0; // Cooling

    std::map<size_t, EnergySource> sources;
    sources[setup.layer_id] = src;
    dv.PopulateEnergySources(sources, setup.bed);

    auto& layer = setup.bed.layers_.Get(setup.layer_id);
    double area = setup.bed.GetArea();
    double cell_volume = layer.domain_.GetCellVolume(area);
    double V_layer = cell_volume * layer.domain_.n_cells;
    double expected_q = -200.0 / V_layer;

    for (int i = layer.domain_.start_index; i <= layer.domain_.end_index; ++i) {
      EXPECT_DOUBLE_EQ(dv.energy_source_power[i], expected_q);
    }
  }

  TEST(EnergySourceTest, StepOwnsEnergySources) {
    Step step;
    EXPECT_TRUE(step.energy_sources_.empty());

    EnergySource src;
    src.enabled_ = true;
    src.power_ = 100.0;
    step.energy_sources_[0] = src;

    EXPECT_EQ(step.energy_sources_.size(), 1);
    EXPECT_TRUE(step.energy_sources_.at(0).enabled_);
    EXPECT_DOUBLE_EQ(step.energy_sources_.at(0).power_, 100.0);
  }

  // ---------------------------------------------------------------------------
  // Logged adsorption energy sources (L-12). The "<comp> Component Energy
  // Source" and "Energy Source" logs are diagnostics only; these tests pin them
  // to the term the energy residual actually uses, and the total to the sum of
  // the same evaluation's components.
  // ---------------------------------------------------------------------------

  namespace {

    // Ward 2024: two components (CO2, N2), INTERNAL_ENERGY convention with a
    // constant adsorbed-phase heat capacity Cpa = 30.7 J/(mol.K), no heaters.
    constexpr const char* kWardDeck = "../examples/decks/Ward2024.json";

    // Initialise the state the way CoupledSolver does before the first step.
    void InitialiseDeck(UserData& data, HeatOfAdsorptionConvention convention,
                        ThermalEquilibriumMode mode) {
      LoadFromJson(data, kWardDeck);
      data.fluid.heat_of_adsorption_convention_ = convention;
      data.bed.equations_.thermal_equilibrium_mode_ = mode;
      data.Initialise();
      auto& first_step = data.cycle.step_list_.front();
      first_step.Initialise(data.state);
      data.state.domain_vectors.PopulateEnergySources(first_step.energy_sources_, data.bed);
    }

    // Residual of the energy balance that carries the adsorption heat: the
    // bed (LTE) balance, or the solid balance under LTNE.
    std::vector<double> AdsorptionEnergyResidual(UserData& data) {
      if (data.bed.equations_.thermal_equilibrium_mode_ == ThermalEquilibriumMode::NON_LOCAL) {
        UpdateSolidTemperature(&data);
        return data.state.solid_temperature.res_var;
      }
      UpdateTemperature(&data);
      return data.state.fluid_temperature.res_var;
    }

  } // namespace

  class ComponentEnergySourceLog
    : public ::testing::TestWithParam<std::tuple<HeatOfAdsorptionConvention, ThermalEquilibriumMode>> {};

  // The logged per-component energy sources must sum to the adsorption source
  // term of the energy residual at the same state, under either heat-of-
  // adsorption convention and either thermal model. The residual term is
  // isolated by evaluating the residual with and without the uptake rates
  // dq/dt, which enter that residual only through the adsorption source.
  TEST_P(ComponentEnergySourceLog, MatchesResidualSourceTerm) {
    const auto [convention, mode] = GetParam();
    UserData data;
    InitialiseDeck(data, convention, mode);
    auto& state = data.state;
    ASSERT_EQ(data.fluid.heat_of_adsorption_convention_, convention);

    // Off-equilibrium state: distinct fluid and solid temperatures and
    // non-zero uptake rates of both signs, so the check sees which heat and
    // which temperature the logged source uses.
    const size_t N = state.fluid_temperature.var.size();
    for (size_t i = 0; i < N; ++i) {
      state.fluid_temperature.var[i] = 300.0 + 2.0 * i;
      state.solid_temperature.var[i] = 330.0 - 1.5 * i;
    }
    int k = 0;
    for (auto& [ads_id, comp_map] : state.ads_concentration) {
      for (auto& [comp, q] : comp_map) {
        const double sign = (k++ % 2 == 0) ? 1.0 : -0.5;
        for (size_t i = 0; i < N; ++i) q.ddt_var[i] = sign * 1e-3 * (1.0 + 0.1 * i);
      }
    }

    state.UpdateIndependentEquations(&data);
    state.UpdateDiagnosticEquations(&data);
    std::vector<double> logged(N, 0.0);
    for (const auto& [comp, Sei] : state.component_energy_source)
      for (size_t i = 0; i < N; ++i) logged[i] += Sei.var[i];

    const std::vector<double> res_with = AdsorptionEnergyResidual(data);
    for (auto& [ads_id, comp_map] : state.ads_concentration)
      for (auto& [comp, q] : comp_map)
        std::fill(q.ddt_var.begin(), q.ddt_var.end(), 0.0);
    const std::vector<double> res_without = AdsorptionEnergyResidual(data);

    const auto& dV = state.domain_vectors.volume;
    int checked = 0;
    for (size_t i = 1; i + 1 < N; ++i) {
      if (state.domain_vectors.IsVoidCell(i)) continue;
      const double expected = (res_without[i] - res_with[i]) / dV[i];  // [W/m3]
      ASSERT_GT(std::abs(expected), 1.0) << "cell " << i << ": test state has no adsorption heat";
      const double tol = 1e-12 * (std::abs(res_with[i]) + std::abs(res_without[i])) / dV[i]
                       + 1e-10 * std::abs(expected);
      EXPECT_NEAR(logged[i], expected, tol) << "cell " << i;
      ++checked;
    }
    EXPECT_GT(checked, 0);
  }

  INSTANTIATE_TEST_SUITE_P(
    ConventionAndThermalMode, ComponentEnergySourceLog,
    ::testing::Combine(
      ::testing::Values(HeatOfAdsorptionConvention::ISOSTERIC, HeatOfAdsorptionConvention::INTERNAL_ENERGY),
      ::testing::Values(ThermalEquilibriumMode::LOCAL, ThermalEquilibriumMode::NON_LOCAL)),
    [](const auto& info) {
      const bool isosteric = std::get<0>(info.param) == HeatOfAdsorptionConvention::ISOSTERIC;
      const bool lte = std::get<1>(info.param) == ThermalEquilibriumMode::LOCAL;
      return std::string(isosteric ? "Isosteric" : "InternalEnergy") + (lte ? "_LTE" : "_LTNE");
    });

  // At every logged time, the logged totals must equal the sum of the same
  // evaluation's logged per-component sources (plus the heater power for the
  // energy total), i.e. they must not lag one evaluation behind.
  TEST(EnergySourceLog, AggregatesEqualSumOfComponentsAtEveryLoggedTime) {
    UserData data;
    LoadFromJson(data, kWardDeck);
    ASSERT_EQ(data.fluid.heat_of_adsorption_convention_, HeatOfAdsorptionConvention::INTERNAL_ENERGY);
    data.cycle.number_of_cycles_ = 1;
    for (const auto& step : data.cycle.step_list_)
      ASSERT_TRUE(step.energy_sources_.empty()) << "heater power is not logged per step";

    CoupledSolver solver{ data };
    ASSERT_GE(solver.IntegrateCycle(data.cycle, false), 0) << "solver did not complete";

    const auto& sim_log = solver.GetLog().GetSimulationLog();
    ASSERT_FALSE(sim_log.empty());
    size_t rows = 0;
    for (const auto& [cycle_index, steps] : sim_log) {
      for (const auto& [step_index, log] : steps) {
        const auto& energy = log.variables.at(VariableType::ENERGY_SOURCE).values;
        const auto& mass = log.variables.at(VariableType::MASS_SOURCE).values;
        ASSERT_EQ(energy.size(), log.step_data.times.size());
        for (size_t t = 0; t < energy.size(); ++t) {
          const size_t N = energy[t].size();
          std::vector<double> e_sum(N, 0.0), e_abs(N, 0.0), m_sum(N, 0.0), m_abs(N, 0.0);
          for (const auto& [comp, vars] : log.component_variables) {
            const auto& Se = vars.at(ComponentVariableType::COMPONENT_ENERGY_SOURCE).values[t];
            const auto& Sm = vars.at(ComponentVariableType::COMPONENT_MASS_SOURCE).values[t];
            for (size_t i = 0; i < N; ++i) {
              e_sum[i] += Se[i]; e_abs[i] += std::abs(Se[i]);
              m_sum[i] += Sm[i]; m_abs[i] += std::abs(Sm[i]);
            }
          }
          for (size_t i = 0; i < N; ++i) {
            EXPECT_NEAR(energy[t][i], e_sum[i], 1e-12 * e_abs[i] + 1e-300)
              << "cycle " << cycle_index << " step " << step_index << " row " << t << " cell " << i;
            EXPECT_NEAR(mass[t][i], m_sum[i], 1e-12 * m_abs[i] + 1e-300)
              << "cycle " << cycle_index << " step " << step_index << " row " << t << " cell " << i;
          }
          ++rows;
        }
      }
    }
    EXPECT_GT(rows, 0u);
  }

} // end adsim namespace
