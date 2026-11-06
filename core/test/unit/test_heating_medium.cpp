// test_heating_medium.cpp
//
// A step's heating medium (a jacket on the wall's outer face). These tests pin:
//   * a disabled medium (the default) changes nothing;
//   * an enabled one replaces the outer-face coefficient with h_m and the
//     ambient temperature with T_m(t), constant or ramped, in the wall energy
//     residual, through exactly the outer-face term
//     4 d_o h_m (T_w - T_m) / (d_o^2 - d_i^2);
//   * the registry lists step[ID].medium.* for every step, and the values
//     round-trip through a saved deck (and a copied step);
//   * PerformChecks rejects a negative h_m and a non-positive T_m.

#include <cmath>
#include <string>

#include <gtest/gtest.h>

#include "../src/pch.h"
#include "adsim/boundary_conditions/step.h"
#include "adsim/equations/adsorbent_equations.h"
#include "adsim/objects/user_data.h"
#include "adsim/parameters/parameter_registry.h"

#include "../fixtures/user_data_test.h"

namespace adsim {

  namespace {

    constexpr double kTw0 = 320.0;   // [K] initial wall temperature
    constexpr double kAmbient = 300.0;

    // A closed, non-adsorbing single-layer bed whose layer has no outer-face
    // coefficient, so the ambient term of the wall residual is zero and only
    // the medium's term can appear.
    class HeatingMediumTest : public UserDataTest {
    protected:
      LayerID layer_id{};

      void SetUp() override {
        UserDataTest::SetUp();
        AddComponents({"Nitrogen", "Oxygen"});
        auto ids = AddLayer("Packing");
        layer_id = ids.first;
        auto& layer = GetLayer(layer_id);
        layer.domain_.length = 0.5;
        layer.domain_.n_cells = 10;
        layer.eb_ = 0.4;
        layer.particle_diameter_ = 0.002;
        std::get<CylinderBed>(data.bed.shape_).diameter = 0.05;
        GetAdsorbent(ids.second).solid_density_ = 1000.0;
        GetAdsorbent(ids.second).ep_ = 0.0;
        AddConstantInternalHeatTransferCoefficientModel(50.0, layer_id);
        AddConstantExternalHeatTransferCoefficientModel(0.0, layer_id);
        data.bed.ambient_temperature_ = kAmbient;
        AddInitialConditions(1e5, 310.0, kTw0, {{"Nitrogen", 1.0}, {"Oxygen", 0.0}});
        AddStep("Heat", 100.0);
        AddClosedBC(BoundaryLocation::BOTTOM, "Heat");
        AddClosedBC(BoundaryLocation::TOP, "Heat");
        ConfigureCycle(1);
      }

      Step& HeatStep() { return data.cycle.GetStep("Heat"); }

      // Wall residual at the initial state, `t_rel` seconds into the step.
      std::vector<double> WallResidual(double t_rel) {
        data.Initialise();
        HeatStep().Initialise(data.state);
        data.step_start_time_ = 10.0;
        data.residual_time_ = 10.0 + t_rel;
        UpdateExternalHeatTransferCoefficient(&data);
        UpdateInternalHeatTransferCoefficient(&data);
        UpdateWallTemperature(&data);
        return data.state.wall_temperature.res_var;
      }

      // The outer-face term the medium must add to cell i.
      double OuterFaceTerm(int i, double h, double T_medium) const {
        const auto& dv = data.state.domain_vectors;
        const double d_o = dv.de[i], d_i = dv.di[i];
        return 4.0 * d_o * h * (data.state.wall_temperature.var[i] - T_medium) / (d_o * d_o - d_i * d_i);
      }
    };

  } // namespace

  TEST(HeatingMedium, DisabledByDefault) {
    Step step;
    EXPECT_FALSE(step.heating_medium_.enabled_);
    EXPECT_EQ(step.ActiveHeatingMedium(), nullptr);
    step.heating_medium_.enabled_ = true;
    EXPECT_EQ(step.ActiveHeatingMedium(), &step.heating_medium_);
  }

  TEST(HeatingMedium, TemperatureFollowsTheProfile) {
    HeatingMedium m;
    m.temperature_ = 350.0;
    EXPECT_DOUBLE_EQ(m.Temperature(0.0), 350.0);
    EXPECT_DOUBLE_EQ(m.Temperature(1e3), 350.0);
    m.temperature_profile_ = CreateInletTemperatureProfile(BoundaryValueType::LINEAR, 300.0, 100.0);
    static_cast<LinearBoundaryValue*>(m.temperature_profile_.get())->end_value_ = 400.0;
    EXPECT_DOUBLE_EQ(m.Temperature(0.0), 300.0);
    EXPECT_DOUBLE_EQ(m.Temperature(25.0), 325.0);
    EXPECT_DOUBLE_EQ(m.Temperature(500.0), 400.0);  // holds the end value after the ramp

    // Copies are deep: the profile is cloned, not shared.
    HeatingMedium copy = m;
    static_cast<LinearBoundaryValue*>(m.temperature_profile_.get())->end_value_ = 500.0;
    EXPECT_DOUBLE_EQ(copy.Temperature(500.0), 400.0);
  }

  TEST_F(HeatingMediumTest, DisabledMediumKeepsTheAmbientExchange) {
    const auto base = WallResidual(25.0);
    // The layer's outer coefficient is zero: the ambient term vanishes and
    // the logged coefficient is the layer's.
    for (int i = 1; i <= 10; ++i) EXPECT_EQ(data.state.external_heat_transfer_coefficient.var[i], 0.0);

    // Same deck with a non-zero layer coefficient: the residual moves by the
    // ambient term exactly (the medium path is not taken).
    AddConstantExternalHeatTransferCoefficientModel(20.0, layer_id);
    const auto amb = WallResidual(25.0);
    for (int i = 1; i <= 10; ++i)
      EXPECT_NEAR(amb[i] - base[i], OuterFaceTerm(i, 20.0, kAmbient), 1e-9 * std::abs(amb[i])) << i;
  }

  TEST_F(HeatingMediumTest, ConstantMediumReplacesAmbientAndCoefficient) {
    const auto base = WallResidual(25.0);

    auto& medium = HeatStep().heating_medium_;
    medium.enabled_ = true;
    medium.temperature_ = 400.0;
    medium.heat_transfer_coefficient_ = 250.0;
    // A layer coefficient that the medium must override.
    AddConstantExternalHeatTransferCoefficientModel(20.0, layer_id);
    const auto with = WallResidual(25.0);

    const auto& ho = data.state.external_heat_transfer_coefficient.var;
    for (int i = 0; i < static_cast<int>(ho.size()); ++i) EXPECT_DOUBLE_EQ(ho[i], 250.0) << i;
    for (int i = 1; i <= 10; ++i) {
      const double expected = OuterFaceTerm(i, 250.0, 400.0);
      EXPECT_LT(expected, 0.0);  // a hotter medium heats the wall
      EXPECT_NEAR(with[i] - base[i], expected, 1e-9 * std::abs(expected)) << i;
    }
  }

  TEST_F(HeatingMediumTest, RampedMediumUsesTimeSinceStepStart) {
    const auto base = WallResidual(25.0);
    auto& medium = HeatStep().heating_medium_;
    medium.enabled_ = true;
    medium.heat_transfer_coefficient_ = 100.0;
    medium.temperature_profile_ = CreateInletTemperatureProfile(BoundaryValueType::LINEAR, 300.0, 100.0);
    static_cast<LinearBoundaryValue*>(medium.temperature_profile_.get())->end_value_ = 400.0;

    const auto at25 = WallResidual(25.0);  // T_m = 325 K
    for (int i = 1; i <= 10; ++i)
      EXPECT_NEAR(at25[i] - base[i], OuterFaceTerm(i, 100.0, 325.0), 1e-9 * std::abs(at25[i])) << i;
    const auto at75 = WallResidual(75.0);  // T_m = 375 K
    for (int i = 1; i <= 10; ++i)
      EXPECT_NEAR(at75[i] - base[i], OuterFaceTerm(i, 100.0, 375.0), 1e-9 * std::abs(at75[i])) << i;
  }

  TEST_F(HeatingMediumTest, ChecksRejectInvalidValues) {
    std::string error;
    ASSERT_TRUE(data.PerformChecks(error)) << error;

    auto& medium = HeatStep().heating_medium_;
    medium.heat_transfer_coefficient_ = -1.0;  // ignored while disabled
    EXPECT_TRUE(data.PerformChecks(error)) << error;

    medium.enabled_ = true;
    error.clear();
    EXPECT_FALSE(data.PerformChecks(error));
    EXPECT_NE(error.find("heating medium of step 'Heat'"), std::string::npos) << error;

    medium.heat_transfer_coefficient_ = 100.0;
    medium.temperature_ = 0.0;
    error.clear();
    EXPECT_FALSE(data.PerformChecks(error));
    EXPECT_NE(error.find("above absolute zero"), std::string::npos) << error;

    medium.temperature_ = 350.0;
    error.clear();
    EXPECT_TRUE(data.PerformChecks(error)) << error;
  }

  TEST_F(HeatingMediumTest, RegistryParametersRoundTrip) {
    data.cycle.EnsureIds();
    const std::string sp = "step[" + std::to_string(HeatStep().id_) + "].medium";
    {
      ParameterRegistry reg(data);
      // Listed for a step without a medium, so a sweep can switch one on.
      ASSERT_NE(reg.Find(sp + ".enabled"), nullptr);
      ASSERT_NE(reg.Find(sp + ".temperature"), nullptr);
      ASSERT_NE(reg.Find(sp + ".heat_transfer_coefficient"), nullptr);
      EXPECT_EQ(reg.Find(sp + ".temperature"), reg.Find("step[Heat].medium.temperature"));  // alias
      EXPECT_EQ(reg.Get(sp + ".enabled").Info().kind, ParameterKind::Boolean);
      EXPECT_TRUE(reg.Get(sp + ".temperature").Info().tunable);
      EXPECT_EQ(reg.Get(sp + ".temperature").Info().unit, "K");
      EXPECT_TRUE(reg.Get(sp + ".heat_transfer_coefficient").Info().tunable);

      reg.Get(sp + ".enabled").SetText("true");
      reg.Get(sp + ".temperature").SetValue(373.15);
      reg.Get(sp + ".heat_transfer_coefficient").SetText("0.25 [kW/m2-K]");
      data.EvaluateAllExpressions();
    }
    EXPECT_TRUE(HeatStep().heating_medium_.enabled_);
    EXPECT_DOUBLE_EQ(HeatStep().heating_medium_.temperature_, 373.15);
    EXPECT_NEAR(HeatStep().heating_medium_.heat_transfer_coefficient_, 250.0, 1e-9);

    std::string json;
    SaveToString(data, json);
    UserData loaded;
    LoadFromString(loaded, json);
    loaded.EvaluateAllExpressions();
    const auto& m = loaded.cycle.GetStep("Heat").heating_medium_;
    EXPECT_TRUE(m.enabled_);
    EXPECT_DOUBLE_EQ(m.temperature_, 373.15);
    EXPECT_NEAR(m.heat_transfer_coefficient_, 250.0, 1e-9);
    ParameterRegistry reg2(loaded);
    EXPECT_DOUBLE_EQ(reg2.Get(sp + ".temperature").GetValue(), 373.15);
    EXPECT_EQ(reg2.Get(sp + ".enabled").GetValue(), 1.0);

    // A ramp lists its start, end and time instead of the constant.
    auto& lm = loaded.cycle.GetStep("Heat").heating_medium_;
    lm.temperature_profile_ = CreateInletTemperatureProfile(BoundaryValueType::LINEAR, 300.0, 60.0);
    ParameterRegistry reg3(loaded);
    EXPECT_EQ(reg3.Find(sp + ".temperature"), nullptr);
    ASSERT_NE(reg3.Find(sp + ".temperature_linear_end"), nullptr);
    reg3.Get(sp + ".temperature_linear_end").SetValue(420.0);
    std::string json2;
    SaveToString(loaded, json2);
    UserData loaded2;
    LoadFromString(loaded2, json2);
    loaded2.EvaluateAllExpressions();
    EXPECT_DOUBLE_EQ(loaded2.cycle.GetStep("Heat").heating_medium_.Temperature(60.0), 420.0);
    EXPECT_DOUBLE_EQ(loaded2.cycle.GetStep("Heat").heating_medium_.Temperature(0.0), 300.0);

    // A duplicated step carries its medium.
    const size_t dup = loaded2.cycle.DuplicateStep(0);
    ASSERT_LT(dup, loaded2.cycle.step_list_.size());
    EXPECT_TRUE(loaded2.cycle.step_list_[dup].heating_medium_.enabled_);
    EXPECT_DOUBLE_EQ(loaded2.cycle.step_list_[dup].heating_medium_.Temperature(60.0), 420.0);
  }

} // namespace adsim
