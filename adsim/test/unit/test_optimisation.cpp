#include <gtest/gtest.h>

#include <cmath>
#include <limits>
#include <set>
#include <string>
#include <vector>

#include "../fixtures/user_data_test.h"
#include "adsim/equations/heat_transfer_model.h"
#include "adsim/equations/isotherm.h"
#include "adsim/boundary_conditions/boundary_conditions.h"
#include "adsim/boundary_conditions/boundary_value.h"
#include "adsim/solver/optimisation.h"

namespace adsim {

// Minimise a shifted 2D quadratic: f(x, y) = (x - 1.2)^2 + (y + 0.7)^2.
// Verifies the dlib-backed GlobalOptimizer finds the optimum, respects the
// evaluation budget, and writes values back through OptimizationParameters.
TEST(GlobalOptimizerTest, FindsQuadraticMinimum) {
  double x = 0.0;
  double y = 0.0;

  OptimizationParameters params;
  params.Add("x", x, -5.0, 5.0);
  params.Add("y", y, -5.0, 5.0);

  auto objective = [&]() {
    return (x - 1.2) * (x - 1.2) + (y + 0.7) * (y + 0.7);
  };

  OptimizerSettings settings;
  settings.n_init_samples = 5;
  settings.n_iterations = 45;

  GlobalOptimizer optimizer(params, objective, settings);
  std::vector<double> result = optimizer.OptimizeParameters();

  ASSERT_EQ(result.size(), 2u);
  EXPECT_NEAR(result[0], 1.2, 1e-2);
  EXPECT_NEAR(result[1], -0.7, 1e-2);
  EXPECT_EQ(optimizer.GetIterationCount(),
            settings.n_init_samples + settings.n_iterations);
}

// The optimizer must only ever evaluate within the parameter bounds.
TEST(GlobalOptimizerTest, RespectsBounds) {
  double x = 0.5;

  OptimizationParameters params;
  params.Add("x", x, 0.1, 2.0);

  bool out_of_bounds = false;
  auto objective = [&]() {
    if (x < 0.1 || x > 2.0) out_of_bounds = true;
    return std::abs(x - 1.0);
  };

  OptimizerSettings settings;
  settings.n_init_samples = 3;
  settings.n_iterations = 12;

  GlobalOptimizer optimizer(params, objective, settings);
  std::vector<double> result = optimizer.OptimizeParameters();

  EXPECT_FALSE(out_of_bounds);
  EXPECT_GE(result[0], 0.1);
  EXPECT_LE(result[0], 2.0);
  EXPECT_NEAR(result[0], 1.0, 5e-2);
}

// L-04: a candidate whose objective or constraint KPI is undefined (non-finite
// in the cost map; the KPI itself returns 0 with defined=false) must get the
// large FINITE penalty in both minimise and maximise modes, and must never be
// cheaper than a genuinely scored candidate.
TEST(OptimizationCostFunctionTest, UndefinedKPIGetsFinitePenalty) {
  const std::string CO2 = "CarbonDioxide";
  const double nan = std::numeric_limits<double>::quiet_NaN();
  const double inf = std::numeric_limits<double>::infinity();

  for (ObjectiveDirection dir : { ObjectiveDirection::MINIMIZE, ObjectiveDirection::MAXIMIZE }) {
    OptimizationCostFunction fn;
    fn.AddObjective(CO2, KPIType::ENERGY_CONSUMPTION, 1.0, dir);

    // A scored candidate (energy 5e5 J/kg).
    PerComponentKPIValues good;
    good[CO2][KPIType::ENERGY_CONSUMPTION] = 5e5;
    const double good_cost = fn.ComputeCost(good);
    EXPECT_DOUBLE_EQ(good_cost, dir == ObjectiveDirection::MINIMIZE ? 5e5 : -5e5);

    for (double undefined : { nan, inf, -inf }) {
      PerComponentKPIValues bad;
      bad[CO2][KPIType::ENERGY_CONSUMPTION] = undefined;
      const double c = fn.ComputeCost(bad);
      EXPECT_TRUE(std::isfinite(c));
      EXPECT_DOUBLE_EQ(c, kOptimisationPenaltyCost);
      EXPECT_GT(c, good_cost);
    }
  }

  // Undefined constraint KPI (objective defined): penalty, in both senses.
  for (ConstraintSense sense : { ConstraintSense::GREATER_OR_EQUAL, ConstraintSense::LESS_OR_EQUAL }) {
    OptimizationCostFunction fn;
    fn.AddObjective(CO2, KPIType::RECOVERY, 1.0, ObjectiveDirection::MAXIMIZE);
    fn.AddConstraint(CO2, KPIType::PURITY, 0.9, 10.0, sense);
    PerComponentKPIValues kpis;
    kpis[CO2][KPIType::RECOVERY] = 0.8;
    kpis[CO2][KPIType::PURITY] = nan;
    const double c = fn.ComputeCost(kpis);
    EXPECT_TRUE(std::isfinite(c));
    EXPECT_DOUBLE_EQ(c, kOptimisationPenaltyCost);
  }
}

// GlobalOptimizer must hand dlib only finite values: an objective that returns
// NaN/Inf over part of the domain is settled as the penalty, and the optimiser
// still converges on the well-defined region.
TEST(GlobalOptimizerTest, NonFiniteObjectiveIsPenalised) {
  double x = 0.0;
  OptimizationParameters params;
  params.Add("x", x, -2.0, 2.0);

  bool saw_undefined = false;
  auto objective = [&]() {
    if (x < -1.0) {
      saw_undefined = true;
      return std::numeric_limits<double>::quiet_NaN();  // "no product" region
    }
    return (x - 0.5) * (x - 0.5);
  };

  OptimizerSettings settings;
  settings.n_init_samples = 5;
  settings.n_iterations = 25;
  GlobalOptimizer optimizer(params, objective, settings);
  std::vector<double> result = optimizer.OptimizeParameters();

  ASSERT_EQ(result.size(), 1u);
  EXPECT_TRUE(saw_undefined);  // dlib probes the bounds, so the region is hit
  EXPECT_NEAR(result[0], 0.5, 5e-2);
}

// AvailableParameters must expose every tunable coefficient of the Arrhenius
// Sips isotherm and the Leva internal-HTC correlation as decision variables
// (G-08 companion: neither model was enumerated, so e.g. the Casas 2012 fit
// parameters could not be optimised / parameterised), and selecting one must
// write through to the model's ExpressionParam.
class AvailableParametersTest : public UserDataTest {
protected:
  static OptimizationParameter* Find(const AvailableParameters& ap,
                                     const std::string& name,
                                     const std::string& category,
                                     const std::string& subcategory) {
    for (const auto& p : ap.GetParameters())
      if (p.name == name && p.category == category && p.subcategory == subcategory)
        return const_cast<OptimizationParameter*>(&p);  // SetValue is non-const
    return nullptr;
  }
};

TEST_F(AvailableParametersTest, EnumeratesArrheniusSipsAndLevaParameters) {
  const std::string CO2 = "CarbonDioxide";
  const std::string H2 = "Hydrogen";
  AddComponents({ CO2, H2 });
  auto [layer_id, ads_id] = AddLayer("AP3-60");

  auto sips = std::make_unique<ArrheniusSipsIsotherm>();
  ArrheniusSipsIsotherm* sips_ptr = sips.get();
  sips->params_[CO2].c_beta = 0.0928714;
  sips->params_[H2].k_B = 9954.73;
  AddIsotherm({ CO2, H2 }, ads_id, std::move(sips));
  AddLinearDrivingForceModel(ads_id, CO2, 0.15);
  AddLinearDrivingForceModel(ads_id, H2, 1.0);

  auto leva = std::make_unique<InternalHeatTransferModel_Leva>();
  InternalHeatTransferModel_Leva* leva_ptr = leva.get();
  leva->eta1 = 41.13;
  GetLayer(layer_id).equations_.internal_heat_transfer_model_ = std::move(leva);

  AvailableParameters ap;
  ap.UpdateParameterList(data);

  // Isotherm: all 8 coefficients, for each component, with the documented units
  // and bounds that bracket the current value.
  const std::vector<std::pair<std::string, std::string>> iso_params = {
    { "n_inf_a", "mol/kg" }, { "n_inf_b", "J/mol" }, { "k_A", "1/Pa" }, { "k_B", "J/mol" },
    { "c_ref", "-" }, { "c_amp", "-" }, { "c_beta", "1/K" }, { "c_Tref", "K" },
  };
  for (const std::string& comp : { CO2, H2 }) {
    const std::string sub = "AP3-60, AP3-60, " + comp;
    for (const auto& [name, units] : iso_params) {
      const auto* p = Find(ap, name, "isotherm", sub);
      ASSERT_NE(p, nullptr) << "missing isotherm parameter " << name << " for " << comp;
      EXPECT_EQ(p->units, units) << name;
      EXPECT_LE(p->lower_bound, p->default_value) << name << " (" << comp << ")";
      EXPECT_GE(p->upper_bound, p->default_value) << name << " (" << comp << ")";
    }
  }
  auto* c_beta = Find(ap, "c_beta", "isotherm", "AP3-60, AP3-60, " + CO2);
  ASSERT_NE(c_beta, nullptr);
  EXPECT_DOUBLE_EQ(c_beta->default_value, 0.0928714);
  c_beta->SetValue(0.05);
  EXPECT_DOUBLE_EQ(sips_ptr->params_.at(CO2).c_beta.Get(), 0.05);

  const auto* k_B = Find(ap, "k_B", "isotherm", "AP3-60, AP3-60, " + H2);
  ASSERT_NE(k_B, nullptr);
  EXPECT_DOUBLE_EQ(k_B->default_value, 9954.73);

  // Leva internal HTC: eta1, eta2, exp_coeff on the layer.
  auto* eta1 = Find(ap, "Leva Re coefficient (eta1)", "layer", "AP3-60");
  auto* eta2 = Find(ap, "Leva Re exponent (eta2)", "layer", "AP3-60");
  auto* expc = Find(ap, "Leva wall-decay coefficient (exp_coeff)", "layer", "AP3-60");
  ASSERT_NE(eta1, nullptr);
  ASSERT_NE(eta2, nullptr);
  ASSERT_NE(expc, nullptr);
  EXPECT_DOUBLE_EQ(eta1->default_value, 41.13);
  EXPECT_DOUBLE_EQ(eta2->default_value, 0.90);
  EXPECT_DOUBLE_EQ(expc->default_value, 6.0);
  for (const auto* p : { eta1, eta2, expc }) {
    EXPECT_LE(p->lower_bound, p->default_value) << p->name;
    EXPECT_GE(p->upper_bound, p->default_value) << p->name;
  }
  eta2->SetValue(0.32);
  EXPECT_DOUBLE_EQ(leva_ptr->eta2.Get(), 0.32);
}

// F-02: the LDF2 rate dq/dt = k (q* - q)|q* - q| makes k a 1/(s mol/kg) =
// kg/(mol s) coefficient, so its decision variables (constant k and the
// Arrhenius pre-factor) must carry that unit -- and LDF2-appropriate bounds --
// while first-order LDF components keep 1/s.
TEST_F(AvailableParametersTest, LDF2RateCoefficientUsesKgPerMolSecond) {
  const std::string CO2 = "CarbonDioxide";
  const std::string N2 = "Nitrogen";
  const std::string H2 = "Hydrogen";
  AddComponents({ CO2, N2, H2 });
  auto ads_id = AddLayer("Zeolite").second;
  AddIsotherm({ CO2, N2, H2 }, ads_id, std::make_unique<LangmuirIsotherm>());

  AddLinearDrivingForceModel(ads_id, CO2, 0.15);              // LDF, constant k
  AddLinearDrivingForceModelSecondOrder(ads_id, N2, 0.5);     // LDF2, constant k
  {                                                           // LDF2, Arrhenius k
    auto ldf2 = std::make_unique<LinearDrivingForceSecondOrder>();
    auto arr = std::make_unique<MassTransferCoefficientArrhenius>();
    arr->pre_exp_factor_ = 20.0;
    arr->activation_energy_ = 1.5e4;
    ldf2->mtc_model_ = std::move(arr);
    GetAdsorbent(ads_id).adsorption_source_.mass_transfer_models_[H2] = std::move(ldf2);
  }

  AvailableParameters ap;
  ap.UpdateParameterList(data);

  const std::string sub = "Zeolite, Zeolite, ";
  const auto* k_ldf = Find(ap, "mass transfer coefficient", "adsorbent", sub + CO2);
  ASSERT_NE(k_ldf, nullptr);
  EXPECT_EQ(k_ldf->units, "1/s");
  EXPECT_DOUBLE_EQ(k_ldf->lower_bound, 1e-4);
  EXPECT_DOUBLE_EQ(k_ldf->upper_bound, 100.0);

  auto* k_ldf2 = Find(ap, "mass transfer coefficient", "adsorbent", sub + N2);
  ASSERT_NE(k_ldf2, nullptr);
  EXPECT_EQ(k_ldf2->units, "kg/(mol-s)");
  EXPECT_DOUBLE_EQ(k_ldf2->lower_bound, 1e-5);
  EXPECT_DOUBLE_EQ(k_ldf2->upper_bound, 1e3);
  EXPECT_DOUBLE_EQ(k_ldf2->default_value, 0.5);

  const auto* k0_ldf2 = Find(ap, "MTC pre-exponential factor", "adsorbent", sub + H2);
  ASSERT_NE(k0_ldf2, nullptr);
  EXPECT_EQ(k0_ldf2->units, "kg/(mol-s)");
  EXPECT_DOUBLE_EQ(k0_ldf2->lower_bound, 1e-3);
  EXPECT_DOUBLE_EQ(k0_ldf2->upper_bound, 1e5);
  const auto* ea_ldf2 = Find(ap, "MTC activation energy", "adsorbent", sub + H2);
  ASSERT_NE(ea_ldf2, nullptr);
  EXPECT_EQ(ea_ldf2->units, "J/mol");  // unaffected by the rate model

  // Still writes through to the model's ExpressionParam.
  k_ldf2->SetValue(2.0);
  auto* mt = GetAdsorbent(ads_id).adsorption_source_.mass_transfer_models_.at(N2).get();
  auto* c = dynamic_cast<MassTransferCoefficientConstant*>(mt->mtc_model_.get());
  ASSERT_NE(c, nullptr);
  EXPECT_DOUBLE_EQ(c->constant_mtc_.Get(), 2.0);
}

// ---------------------------------------------------------------------------
// G-03: parameter keys. The GUI used to cache AvailableParameters (raw
// references into the model) across frames and only rebuild them when a
// structure hash changed; the hash missed most type changes, so the cached
// references dangled. Callers now keep OptimizationParameter::key and
// re-resolve it with Find after re-enumerating. These tests pin that a key
// resolves either to the parameter backed by the *live* object or to nothing.
// ---------------------------------------------------------------------------
class ParameterKeyTest : public UserDataTest {
protected:
  const std::string CO2 = "CarbonDioxide";
  const std::string N2 = "Nitrogen";
  AdsorbentID ads_id{};
  LangmuirIsotherm* langmuir = nullptr;

  void SetUp() override {
    UserDataTest::SetUp();
    AddComponents({ CO2, N2 });
    ads_id = AddLayer("Zeolite").second;
    auto iso = std::make_unique<LangmuirIsotherm>();
    iso->params_[CO2].m1 = 3.0;
    iso->params_[N2].m1 = 1.0;
    langmuir = iso.get();
    AddIsotherm({ CO2, N2 }, ads_id, std::move(iso));
    AddLinearDrivingForceModel(ads_id, CO2, 0.15);
    AddLinearDrivingForceModel(ads_id, N2, 1.0);

    // Pressure BCs at both ends: each offers a "pressure (constant)" entry in
    // the same "cycle > Feed" subcategory.
    AddStep("Feed", 100.0);
    AddPressureInletBC(CreateConstantBoundaryValue(2e5), { { CO2, 0.15 }, { N2, 0.85 } }, 300.0,
                       BoundaryLocation::BOTTOM, "Feed");
    AddPressureOutletBC(CreateConstantBoundaryValue(1e5), BoundaryLocation::TOP, "Feed");
  }

  Step& Feed() { return data.cycle.GetStep("Feed"); }

  PressureBoundaryCondition* PressureBC(BoundaryLocation loc) {
    auto& bc = loc == BoundaryLocation::TOP ? Feed().top_boundary_condition_
                                            : Feed().bottom_boundary_condition_;
    return dynamic_cast<PressureBoundaryCondition*>(bc.get());
  }

  static ConstantBoundaryValue* Constant(PressureBoundaryCondition* bc) {
    return dynamic_cast<ConstantBoundaryValue*>(bc->pressure_.get());
  }

  // Key of the entry whose reference is exactly `target`.
  static std::string KeyOf(const AvailableParameters& ap, const double& target) {
    for (const auto& p : ap)
      if (&p.value_ref == &target) return p.key;
    return {};
  }
};

TEST_F(ParameterKeyTest, KeysAreUniqueAndDistinguishBoundaryEnds) {
  AvailableParameters ap;
  ap.UpdateParameterList(data);
  ASSERT_GT(ap.Size(), 0u);

  std::set<std::string> keys;
  for (const auto& p : ap) {
    EXPECT_FALSE(p.key.empty()) << p.name;
    EXPECT_TRUE(keys.insert(p.key).second) << "duplicate key " << p.key;
    // Every key resolves to its own entry.
    EXPECT_EQ(ap.Find(p.key), &p) << p.key;
  }

  const std::string bottom = KeyOf(ap, Constant(PressureBC(BoundaryLocation::BOTTOM))->constant_value_.Get());
  const std::string top = KeyOf(ap, Constant(PressureBC(BoundaryLocation::TOP))->constant_value_.Get());
  ASSERT_FALSE(bottom.empty());
  ASSERT_FALSE(top.empty());
  EXPECT_NE(bottom, top);
  EXPECT_EQ(ap.Find(bottom)->name, "pressure (constant)");
  EXPECT_EQ(ap.Find(top)->name, "pressure (constant)");

  EXPECT_EQ(ap.Find("no such parameter"), nullptr);
}

TEST_F(ParameterKeyTest, KeysAreStableAcrossReEnumeration) {
  AvailableParameters first;
  first.UpdateParameterList(data);
  AvailableParameters second;
  second.UpdateParameterList(data);
  ASSERT_EQ(first.Size(), second.Size());
  for (const auto& p : first) {
    const auto* q = second.Find(p.key);
    ASSERT_NE(q, nullptr) << p.key;
    EXPECT_EQ(&q->value_ref, &p.value_ref) << p.key;
  }
}

TEST_F(ParameterKeyTest, FindReResolvesAfterStructuralEdits) {
  AvailableParameters ap;
  ap.UpdateParameterList(data);
  const std::string bottom_p = KeyOf(ap, Constant(PressureBC(BoundaryLocation::BOTTOM))->constant_value_.Get());
  const std::string top_p = KeyOf(ap, Constant(PressureBC(BoundaryLocation::TOP))->constant_value_.Get());
  const std::string co2_m1 = KeyOf(ap, langmuir->params_.at(CO2).m1.Get());
  const std::string inlet_t = KeyOf(ap, PressureBC(BoundaryLocation::BOTTOM)->temperature_.Get());
  for (const auto* k : { &bottom_p, &top_p, &co2_m1, &inlet_t }) ASSERT_FALSE(k->empty());

  // 1. Profile type: bottom pressure Constant -> Linear destroys the constant.
  PressureBC(BoundaryLocation::BOTTOM)->pressure_ = CreateLinearBoundaryValue(2e5, 3e5, 10.0);
  ap.UpdateParameterList(data);
  EXPECT_EQ(ap.Find(bottom_p), nullptr);
  ASSERT_NE(ap.Find(top_p), nullptr);
  EXPECT_EQ(&ap.Find(top_p)->value_ref, &Constant(PressureBC(BoundaryLocation::TOP))->constant_value_.Get());
  // The BC itself survived, so its inlet temperature still resolves to it.
  ASSERT_NE(ap.Find(inlet_t), nullptr);
  EXPECT_EQ(&ap.Find(inlet_t)->value_ref, &PressureBC(BoundaryLocation::BOTTOM)->temperature_.Get());

  // 2. BC type: top Pressure -> Flow destroys the pressure BC and its profile.
  std::vector<std::unique_ptr<IBoundaryStream>> streams;
  streams.push_back(CreateNoValveOutletStream(CreateVelocityBC(0.1)));
  AddFlowBC(std::move(streams), BoundaryLocation::TOP, BoundaryDirection::OUTLET, "Feed");
  ap.UpdateParameterList(data);
  EXPECT_EQ(ap.Find(top_p), nullptr);

  // 3. Isotherm replaced by a new object of the same type: same key, new
  //    backing object (never the destroyed one).
  auto fresh = std::make_unique<LangmuirIsotherm>();
  fresh->params_[CO2].m1 = 5.0;
  fresh->params_[N2].m1 = 2.0;
  LangmuirIsotherm* fresh_ptr = fresh.get();
  AddIsotherm({ CO2, N2 }, ads_id, std::move(fresh));  // destroys `langmuir`
  langmuir = nullptr;
  ap.UpdateParameterList(data);
  auto* m1 = ap.Find(co2_m1);
  ASSERT_NE(m1, nullptr);
  EXPECT_EQ(&m1->value_ref, &fresh_ptr->params_.at(CO2).m1.Get());
  EXPECT_EQ(m1->expr_param, &fresh_ptr->params_.at(CO2).m1);
  EXPECT_DOUBLE_EQ(m1->default_value, 5.0);
  m1->SetValue(4.2);
  EXPECT_DOUBLE_EQ(fresh_ptr->params_.at(CO2).m1.Get(), 4.2);

  // 4. Isotherm type: Langmuir -> Sips has no "m1" any more.
  GetAdsorbent(ads_id).adsorption_source_.SetIsothermType(CO2, IsothermType::SIPS);
  ap.UpdateParameterList(data);
  EXPECT_EQ(ap.Find(co2_m1), nullptr);

  // SelectByKeys keeps only what still resolves, in the order given.
  auto selected = ap.SelectByKeys({ top_p, inlet_t, bottom_p, co2_m1, "no such parameter" });
  ASSERT_EQ(selected.Size(), 1u);
  EXPECT_EQ(selected[0].key, inlet_t);
  EXPECT_EQ(&selected[0].value_ref, &PressureBC(BoundaryLocation::BOTTOM)->temperature_.Get());
}

TEST_F(ParameterKeyTest, DuplicateStreamNamesGetDistinctKeys) {
  auto a = CreateNoValveInletStream(CreateVelocityBC(0.1), { { CO2, 0.15 }, { N2, 0.85 } }, 300.0);
  auto b = CreateNoValveInletStream(CreateVelocityBC(0.2), { { CO2, 0.5 }, { N2, 0.5 } }, 310.0);
  a->name_ = "Inlet";
  b->name_ = "Inlet";
  BoundaryStream* a_ptr = a.get();
  BoundaryStream* b_ptr = b.get();
  std::vector<std::unique_ptr<IBoundaryStream>> streams;
  streams.push_back(std::move(a));
  streams.push_back(std::move(b));
  AddFlowBC(std::move(streams), BoundaryLocation::BOTTOM, BoundaryDirection::INLET, "Feed");

  AvailableParameters ap;
  ap.UpdateParameterList(data);
  const std::string ka = KeyOf(ap, a_ptr->temperature_.Get());
  const std::string kb = KeyOf(ap, b_ptr->temperature_.Get());
  ASSERT_FALSE(ka.empty());
  ASSERT_FALSE(kb.empty());
  EXPECT_NE(ka, kb);
  EXPECT_EQ(&ap.Find(ka)->value_ref, &a_ptr->temperature_.Get());
  EXPECT_EQ(&ap.Find(kb)->value_ref, &b_ptr->temperature_.Get());
}

}  // namespace adsim
