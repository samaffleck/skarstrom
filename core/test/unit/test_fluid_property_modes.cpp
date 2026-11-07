// Fluid property switches that let a deck reproduce models with fixed gas
// properties (e.g. Yancy-Caballero 2020 / PEESEgroup PSA, Haghpanah 2013):
//   * gas_heat_capacity_mode_ CONSTANT: user molar C_p,i [J/(mol.K)] instead of
//     CoolProp at the initial conditions;
//   * viscosity_mode_ CONSTANT: one mixture viscosity [Pa.s] in every cell;
//   * HeatOfAdsorptionConvention::ISOSTERIC_RT: the stored heat is -dU and the
//     heat released per mole adsorbed is -dU + R*T at the local temperature.
// Also: the defaults are unchanged, invalid constants fail the deck checks,
// and the new fields survive a save -> load (older decks load with defaults).

#include <cmath>
#include <sstream>

#include <cereal/archives/json.hpp>
#include <cereal/external/rapidjson/document.h>
#include <cereal/external/rapidjson/stringbuffer.h>
#include <cereal/external/rapidjson/writer.h>
#include <map>

#include "../fixtures/user_data_test.h"
#include "adsim/equations/adsorbent_equations.h"
#include "adsim/equations/isotherm.h"
#include "adsim/math/constants.h"
#include "adsim/objects/deck_edit.h"
#include "adsim/solver/initial_conditions.h"

namespace adsim {

  namespace {
    const std::string CO2 = "CarbonDioxide";
    const std::string N2 = "Nitrogen";
  }

  class FluidPropertyModesTest : public UserDataTest {
  protected:
    // Inert 6-cell bed of CO2/N2 15/85 at 1 bar, 313.15 K.
    void Build() {
      AddComponents({CO2, N2});
      auto [layer_id, ads_id] = AddLayer("Layer1");
      GetLayer(layer_id).domain_.n_cells = 6;
      auto& ads = GetAdsorbent(ads_id);
      ads.solid_density_ = 1130.0;
      ads.heat_capacity_ = 1070.0;
      ads.thermal_conductivity_ = 0.09;
      auto henry = std::make_unique<HenryIsotherm>();
      for (const auto& c : {CO2, N2}) {
        henry->params_[c].k = 0.0;
        henry->params_[c].Q = 0.0;
      }
      AddIsotherm({CO2, N2}, ads_id, std::move(henry));
      for (const auto& c : {CO2, N2}) {
        AddLinearDrivingForceModel(ads_id, c, 1.0);
        SetHOA(ads_id, c, 0.0);
      }
      AddConstantInternalHeatTransferCoefficientModel(0.0, layer_id);
      AddConstantExternalHeatTransferCoefficientModel(0.0, layer_id);
      std::get<CylinderBed>(data.bed.shape_).diameter = 0.1;
      AddInitialConditions(1.0e5, 313.15, 313.15, {{CO2, 0.15}, {N2, 0.85}});
    }
  };

  TEST_F(FluidPropertyModesTest, DefaultsAreIdealGasAndPolynomial) {
    Build();
    EXPECT_EQ(data.fluid.gas_heat_capacity_mode_, GasHeatCapacityMode::IDEAL_GAS);
    EXPECT_EQ(data.fluid.viscosity_mode_, ViscosityMode::POLYNOMIAL);
    EXPECT_LT(data.fluid.ConstantViscosity(), 0.0);
    data.Initialise();
    // Ideal-gas N2 at 313 K is about 29.2 J/(mol.K), CO2 about 37.6.
    EXPECT_NEAR(data.fluid.MolarHeatCapacity(N2, 313.15), 29.2, 0.3);
    EXPECT_NEAR(data.fluid.MolarHeatCapacity(CO2, 313.15), 38.0, 1.0);
  }

  // Constant 30.7 J/(mol.K) for both components: the molar mixture heat
  // capacity is 30.7 at every composition, so rho_g C_p,g = C * 30.7 in every
  // cell, whatever the mole fractions.
  TEST_F(FluidPropertyModesTest, ConstantGasHeatCapacityIsTheUsersMolarValue) {
    Build();
    data.fluid.gas_heat_capacity_mode_ = GasHeatCapacityMode::CONSTANT;
    data.fluid.gas_heat_capacity_[CO2] = 30.7;
    data.fluid.gas_heat_capacity_[N2] = 30.7;
    data.Initialise();
    EXPECT_NEAR(data.fluid.MolarHeatCapacity(CO2, 313.15), 30.7, 1e-12);
    EXPECT_NEAR(data.fluid.MolarHeatCapacity(N2, 313.15), 30.7, 1e-12);
    EXPECT_NEAR(data.fluid.heat_capacity_.at(CO2), 30.7 / data.fluid.molar_mass_.at(CO2), 1e-9);

    data.state.mole_fraction.at(CO2).var[2] = 0.9;
    data.state.mole_fraction.at(N2).var[2] = 0.1;
    UpdateHeatCapacity(&data);
    const auto& Cp = data.state.heat_capacity.var;
    for (size_t n = 1; n + 1 < Cp.size(); ++n) {
      const double y_co2 = data.state.mole_fraction.at(CO2).var[n];
      const double y_n2 = data.state.mole_fraction.at(N2).var[n];
      const double M = y_co2 * data.fluid.molar_mass_.at(CO2) + y_n2 * data.fluid.molar_mass_.at(N2);
      EXPECT_NEAR(Cp[n] * M, 30.7, 1e-9) << "cell " << n;
    }

    // Different constants mix on the molar basis: 0.15 * 37 + 0.85 * 29 = 30.2.
    data.fluid.gas_heat_capacity_[CO2] = 37.0;
    data.fluid.gas_heat_capacity_[N2] = 29.0;
    data.Initialise();
    const double cp = data.fluid.MixtureHeatCapacity({{CO2, 0.15}, {N2, 0.85}}, 313.15);
    const double M = 0.15 * data.fluid.molar_mass_.at(CO2) + 0.85 * data.fluid.molar_mass_.at(N2);
    EXPECT_NEAR(cp * M, 30.2, 1e-9);
  }

  TEST_F(FluidPropertyModesTest, ConstantViscosityFillsEveryCell) {
    Build();
    data.fluid.viscosity_mode_ = ViscosityMode::CONSTANT;
    data.fluid.viscosity_ = 1.72e-5;
    data.Initialise();
    EXPECT_DOUBLE_EQ(data.fluid.ConstantViscosity(), 1.72e-5);
    for (auto& T : data.state.fluid_temperature.var) T = 400.0;  // temperature does not matter
    UpdateViscosity(&data);
    for (double mu : data.state.viscosity.var) EXPECT_DOUBLE_EQ(mu, 1.72e-5);

    // The polynomial mode gives the temperature-dependent mixture value.
    data.fluid.viscosity_mode_ = ViscosityMode::POLYNOMIAL;
    UpdateViscosity(&data);
    const double expected = 0.15 * Fluid::GetComponentViscosity(CO2, 400.0) + 0.85 * Fluid::GetComponentViscosity(N2, 400.0);
    for (size_t n = 0; n < data.state.viscosity.var.size(); ++n)
      EXPECT_NEAR(data.state.viscosity.var[n], expected, 1e-15) << "cell " << n;
  }

  TEST_F(FluidPropertyModesTest, NonPositiveConstantsFailTheChecks) {
    Build();
    std::string err;
    EXPECT_TRUE(data.fluid.PerformChecks(err)) << err;

    data.fluid.gas_heat_capacity_mode_ = GasHeatCapacityMode::CONSTANT;
    data.fluid.gas_heat_capacity_[N2] = 0.0;
    err.clear();
    EXPECT_FALSE(data.fluid.PerformChecks(err));
    EXPECT_NE(err.find("Nitrogen"), std::string::npos) << err;
    data.fluid.gas_heat_capacity_[N2] = 29.1;

    data.fluid.viscosity_mode_ = ViscosityMode::CONSTANT;
    data.fluid.viscosity_ = -1.0;
    err.clear();
    EXPECT_FALSE(data.fluid.PerformChecks(err));
    data.fluid.viscosity_ = 1.72e-5;
    err.clear();
    EXPECT_TRUE(data.fluid.PerformChecks(err)) << err;
  }

  // A component added through the deck cascade gets an editable constant at
  // once (front ends list the map's entries); one added any other way gets it
  // when the fluid is initialised. Removing the component removes the entry.
  TEST(FluidPropertyModesDeck, ComponentsGetAnEditableConstantHeatCapacity) {
    UserData data;
    std::string err;
    ASSERT_TRUE(AddComponentToDeck(data, CO2, err)) << err;
    ASSERT_EQ(data.fluid.gas_heat_capacity_.count(CO2), 1u);
    EXPECT_DOUBLE_EQ(data.fluid.gas_heat_capacity_.at(CO2).Get(), Fluid::kDefaultGasHeatCapacity);
    RemoveComponentFromDeck(data, CO2);
    EXPECT_EQ(data.fluid.gas_heat_capacity_.count(CO2), 0u);

    Fluid f;
    f.AddComponent(N2);
    EXPECT_EQ(f.gas_heat_capacity_.count(N2), 0u);
    InitialConditions ic;
    ic.pressure = 1e5;
    ic.temperature = 300.0;
    f.Initialise(ic);
    EXPECT_DOUBLE_EQ(f.gas_heat_capacity_.at(N2).Get(), Fluid::kDefaultGasHeatCapacity);
  }

  // The fields round-trip through the deck format; a deck saved before they
  // existed (Fluid class version 3) loads with the defaults and seeded values.
  TEST(FluidPropertyModesSerialisation, RoundTripAndOlderDecks) {
    Fluid f;
    f.AddComponent(CO2);
    f.AddComponent(N2);
    f.gas_heat_capacity_mode_ = GasHeatCapacityMode::CONSTANT;
    f.gas_heat_capacity_[CO2] = 30.7;
    f.gas_heat_capacity_[N2] = 30.7;
    f.viscosity_mode_ = ViscosityMode::CONSTANT;
    f.viscosity_ = 1.72e-5;
    f.heat_of_adsorption_convention_ = HeatOfAdsorptionConvention::ISOSTERIC_RT;
    std::stringstream ss;
    {
      cereal::JSONOutputArchive ar(ss);
      ar(cereal::make_nvp("fluid", f));
    }
    const std::string json = ss.str();
    Fluid g;
    {
      std::istringstream is(json);
      cereal::JSONInputArchive ar(is);
      ar(cereal::make_nvp("fluid", g));
    }
    EXPECT_EQ(g.gas_heat_capacity_mode_, GasHeatCapacityMode::CONSTANT);
    EXPECT_DOUBLE_EQ(g.gas_heat_capacity_.at(CO2).Get(), 30.7);
    EXPECT_EQ(g.viscosity_mode_, ViscosityMode::CONSTANT);
    EXPECT_DOUBLE_EQ(g.viscosity_.Get(), 1.72e-5);
    EXPECT_EQ(g.heat_of_adsorption_convention_, HeatOfAdsorptionConvention::ISOSTERIC_RT);

    // A version-3 deck: the same JSON with the version set to 3 and the four
    // new fields removed.
    rapidjson::Document doc;
    doc.Parse(json.c_str());
    ASSERT_FALSE(doc.HasParseError());
    auto& fl = doc["fluid"];
    ASSERT_TRUE(fl.HasMember("cereal_class_version")) << json;
    EXPECT_EQ(fl["cereal_class_version"].GetUint(), 4u);
    fl["cereal_class_version"].SetUint(3);
    for (const char* key : {"Gas_Heat_Capacity_Mode", "Gas_Heat_Capacity", "Viscosity_Mode", "Viscosity"})
      ASSERT_TRUE(fl.RemoveMember(key)) << key;
    rapidjson::StringBuffer buf;
    rapidjson::Writer<rapidjson::StringBuffer> writer(buf);
    doc.Accept(writer);
    const std::string old_json = buf.GetString();
    Fluid h;
    {
      std::istringstream is(old_json);
      cereal::JSONInputArchive ar(is);
      ar(cereal::make_nvp("fluid", h));
    }
    EXPECT_EQ(h.gas_heat_capacity_mode_, GasHeatCapacityMode::IDEAL_GAS);
    EXPECT_EQ(h.viscosity_mode_, ViscosityMode::POLYNOMIAL);
    EXPECT_DOUBLE_EQ(h.gas_heat_capacity_.at(N2).Get(), Fluid::kDefaultGasHeatCapacity);
    EXPECT_EQ(h.heat_of_adsorption_convention_, HeatOfAdsorptionConvention::ISOSTERIC_RT);
  }

  // ISOSTERIC_RT: the heat per mole adsorbed is the stored -dU plus R*T at the
  // adsorbed-phase temperature. Compared cell by cell with the ISOSTERIC
  // source of the same state: S_rt / S_iso = (dH + R*T) / dH.
  TEST(FluidPropertyModesHeat, IsostericRTAddsRTAtTheLocalTemperature) {
    auto source_at = [](HeatOfAdsorptionConvention convention, std::vector<double>& T_out) {
      UserData data;
      LoadFromJson(data, "../examples/decks/Ward2024.json");
      data.fluid.heat_of_adsorption_convention_ = convention;
      data.Initialise();
      auto& first_step = data.cycle.step_list_.front();
      first_step.Initialise(data.state);
      auto& state = data.state;
      const size_t N = state.fluid_temperature.var.size();
      for (size_t i = 0; i < N; ++i) state.fluid_temperature.var[i] = 300.0 + 5.0 * i;
      for (auto& [ads_id, comp_map] : state.ads_concentration)
        for (auto& [comp, q] : comp_map)
          for (size_t i = 0; i < N; ++i) q.ddt_var[i] = 1e-3 * (1.0 + 0.1 * i);
      state.UpdateIndependentEquations(&data);
      state.UpdateDiagnosticEquations(&data);
      T_out = state.fluid_temperature.var;
      std::map<std::string, std::vector<double>> out;
      for (const auto& [comp, Sei] : state.component_energy_source) out[comp] = Sei.var;
      std::map<std::string, std::vector<double>> dH;
      for (const auto& [ads_id, m] : state.domain_vectors.dH_ads)
        for (const auto& [comp, v] : m) dH[comp] = v;
      return std::make_pair(out, dH);
    };
    std::vector<double> T;
    const auto [iso, dH] = source_at(HeatOfAdsorptionConvention::ISOSTERIC, T);
    const auto [rt, dH2] = source_at(HeatOfAdsorptionConvention::ISOSTERIC_RT, T);
    const double R = constants::kIdealGasConstant;
    int checked = 0;
    for (const auto& [comp, s_iso] : iso) {
      const auto& s_rt = rt.at(comp);
      for (size_t i = 1; i + 1 < s_iso.size(); ++i) {
        const double h = dH.at(comp)[i];
        if (std::abs(s_iso[i]) < 1e-9 || h == 0.0) continue;
        EXPECT_NEAR(s_rt[i] / s_iso[i], (h + R * T[i]) / h, 1e-12) << comp << " cell " << i;
        ++checked;
      }
    }
    EXPECT_GT(checked, 0);
  }

} // end adsim namespace
