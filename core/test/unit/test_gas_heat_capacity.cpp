// Gas-mixture heat capacity: the per-kilogram pure-component C_p,i of Fluid
// must be mixed on the MOLAR basis,
//   C_p,g = sum_i y_i M_i C_p,i / sum_i y_i M_i   [J/(kg.K)],
// not as sum_i y_i C_p,i (which is right only when every M_i is equal). Checks
// Fluid::MixtureHeatCapacity against hand calculations, the solver's cell
// values (UpdateHeatCapacity) against Fluid::MixtureHeatCapacity, and that a
// single-component gas keeps its pure-component C_p.

#include "../fixtures/user_data_test.h"
#include "coolprop_reference_values.h"
#include "adsim/equations/adsorbent_equations.h"
#include "adsim/equations/isotherm.h"

namespace adsim {

  namespace {
    const std::string H2 = "Hydrogen";
    const std::string CO2 = "CarbonDioxide";
    const std::string N2 = "Nitrogen";
    const std::string O2 = "Oxygen";
    // Temperature of the hand fluids (whose held values do not depend on it)
    // and of the solver tests' initial state [K].
    constexpr double kT = 298.15;

    // A fluid with fixed, round property values so the expected numbers can be
    // worked by hand (independent of CoolProp).
    Fluid HandFluid() {
      Fluid f;
      f.molar_mass_ = {{H2, 2.01588e-3}, {CO2, 44.0098e-3}, {N2, 28.0134e-3}};  // [kg/mol]
      f.heat_capacity_ = {{H2, 14300.0}, {CO2, 850.0}, {N2, 1040.0}};          // [J/(kg.K)]
      return f;
    }
  }

  // 50/50 H2/CO2, worked by hand:
  //   molar Cp_i = M_i Cp_i:  H2 2.01588e-3 * 14300 = 28.827084, CO2 44.0098e-3 * 850 = 37.40833 J/(mol.K)
  //   molar mixture Cp = 0.5 * 28.827084 + 0.5 * 37.40833 = 33.117707 J/(mol.K)
  //   M = 0.5 * 2.01588e-3 + 0.5 * 44.0098e-3 = 23.01284e-3 kg/mol
  //   C_p,g = 33.117707 / 23.01284e-3 = 1439.0969 J/(kg.K)
  // The old mole-fraction mixing of the per-kg values, 0.5 * (14300 + 850) =
  // 7575 J/(kg.K), is 174.3 J/(mol.K): more than five times too large.
  TEST(GasHeatCapacityTest, HydrogenCarbonDioxideFiftyFiftyByHand) {
    const Fluid f = HandFluid();
    const std::unordered_map<std::string, double> y{{H2, 0.5}, {CO2, 0.5}};
    const double cp = f.MixtureHeatCapacity(y, kT);
    EXPECT_NEAR(cp, 1439.0969, 1e-3);
    EXPECT_NEAR(cp * 23.01284e-3, 33.117707, 1e-6);  // molar mixture Cp [J/(mol.K)]
    EXPECT_NEAR(f.MolarHeatCapacity(H2, kT), 28.827084, 1e-6);
    EXPECT_NEAR(f.MolarHeatCapacity(CO2, kT), 37.40833, 1e-6);
  }

  // 15/85 CO2/N2, worked by hand:
  //   molar Cp_i: CO2 37.40833, N2 28.0134e-3 * 1040 = 29.133936 J/(mol.K)
  //   molar mixture Cp = 0.15 * 37.40833 + 0.85 * 29.133936 = 30.375095 J/(mol.K)
  //   M = 0.15 * 44.0098e-3 + 0.85 * 28.0134e-3 = 30.41286e-3 kg/mol
  //   C_p,g = 30.375095 / 30.41286e-3 = 998.7583 J/(kg.K)
  // (the old rule gave 0.15 * 850 + 0.85 * 1040 = 1011.5, 1.3 % high).
  TEST(GasHeatCapacityTest, CarbonDioxideNitrogenByHand) {
    const Fluid f = HandFluid();
    const double cp = f.MixtureHeatCapacity({{CO2, 0.15}, {N2, 0.85}}, kT);
    EXPECT_NEAR(cp, 998.7583, 1e-3);
  }

  // A single component, or a mixture with all but one mole fraction zero, keeps
  // the pure-component per-kg heat capacity.
  TEST(GasHeatCapacityTest, SingleComponentIsPureComponentValue) {
    const Fluid f = HandFluid();
    EXPECT_DOUBLE_EQ(f.MixtureHeatCapacity({{N2, 1.0}}, kT), 1040.0);
    EXPECT_DOUBLE_EQ(f.MixtureHeatCapacity({{N2, 1.0}, {CO2, 0.0}, {H2, 0.0}}, kT), 1040.0);
    EXPECT_DOUBLE_EQ(f.MixtureHeatCapacity({{H2, 1.0}}, kT), 14300.0);
    // No gas at all: 0, not NaN.
    EXPECT_EQ(f.MixtureHeatCapacity({{N2, 0.0}}, kT), 0.0);
  }

  // The GAS_PHASE adsorbed-phase heat capacity is the component's own molar
  // gas-phase value C_p,i M_i (a per-component quantity, not the mixture's).
  TEST(GasHeatCapacityTest, GasPhaseAdsorbedHeatCapacityIsComponentMolarValue) {
    Fluid f = HandFluid();
    f.adsorbed_heat_capacity_[CO2] = 12.0;
    f.adsorbed_heat_capacity_mode_ = AdsorbedHeatCapacityMode::CONSTANT;
    EXPECT_DOUBLE_EQ(f.AdsorbedPhaseHeatCapacity(CO2, kT), 12.0);
    f.adsorbed_heat_capacity_mode_ = AdsorbedHeatCapacityMode::GAS_PHASE;
    EXPECT_DOUBLE_EQ(f.AdsorbedPhaseHeatCapacity(CO2, kT), 44.0098e-3 * 850.0);
  }

  // The solver's per-cell C_p,g (UpdateHeatCapacity, used by the energy
  // balances, the convective enthalpy terms, the thermal dispersion and the
  // heat-transfer correlations) uses the same rule, with the engine's properties.
  class GasHeatCapacitySolverTest : public UserDataTest {
  protected:
    void BuildInertBed(const std::vector<std::string>& components,
                       const std::unordered_map<std::string, double>& y0) {
      AddComponents(components);
      auto [layer_id, ads_id] = AddLayer("Layer1");
      GetLayer(layer_id).domain_.n_cells = 6;

      auto& ads = GetAdsorbent(ads_id);
      ads.solid_density_ = 1000.0;
      ads.heat_capacity_ = 900.0;
      ads.thermal_conductivity_ = 0.5;

      auto henry = std::make_unique<HenryIsotherm>();
      for (const auto& c : components) {
        henry->params_[c].k = 0.0;  // inert
        henry->params_[c].Q = 0.0;
      }
      AddIsotherm(components, ads_id, std::move(henry));
      for (const auto& c : components) {
        AddLinearDrivingForceModel(ads_id, c, 1.0);
        SetHOA(ads_id, c, 0.0);
      }
      AddConstantInternalHeatTransferCoefficientModel(10.0, layer_id);
      AddConstantExternalHeatTransferCoefficientModel(10.0, layer_id);
      std::get<CylinderBed>(data.bed.shape_).diameter = 0.05;
      AddInitialConditions(1.0e5, 298.15, 298.15, y0);
      data.Initialise();
    }
  };

  TEST_F(GasHeatCapacitySolverTest, HydrogenCarbonDioxideCellsUseMolarMixing) {
    BuildInertBed({H2, CO2}, {{H2, 0.5}, {CO2, 0.5}});
    const auto& fluid = data.fluid;
    const double expected = fluid.MixtureHeatCapacity({{H2, 0.5}, {CO2, 0.5}}, kT);
    const double M = 0.5 * fluid.molar_mass_.at(H2) + 0.5 * fluid.molar_mass_.at(CO2);

    // At 298.15 K: H2 about 28.8 and CO2 about 37.2 J/(mol.K),
    // so the mixture is about 33 J/(mol.K) (the old rule gave about 174).
    EXPECT_NEAR(expected * M, 33.0, 0.5);

    const auto& Cp = data.state.heat_capacity.var;
    const auto& rho = data.state.fluid_density.var;
    const auto& C = data.state.concentration.var;
    ASSERT_GT(Cp.size(), 2u);
    for (size_t n = 0; n < Cp.size(); ++n) {
      EXPECT_NEAR(Cp[n], expected, 1e-9 * expected) << "cell " << n;
      // rho_g C_p,g is the molar density times the molar mixture Cp.
      const double molar = 0.5 * fluid.MolarHeatCapacity(H2, kT) + 0.5 * fluid.MolarHeatCapacity(CO2, kT);
      EXPECT_NEAR(rho[n] * Cp[n], C[n] * molar, 1e-9 * C[n] * molar) << "cell " << n;
    }

    // A different composition in one cell: recomputed with the same rule.
    data.state.mole_fraction.at(H2).var[2] = 0.9;
    data.state.mole_fraction.at(CO2).var[2] = 0.1;
    UpdateHeatCapacity(&data);
    EXPECT_NEAR(Cp[2], fluid.MixtureHeatCapacity({{H2, 0.9}, {CO2, 0.1}}, kT), 1e-9 * Cp[2]);
  }

  TEST_F(GasHeatCapacitySolverTest, SingleComponentCellsArePureComponentValue) {
    BuildInertBed({N2}, {{N2, 1.0}});
    const double cp_pure = data.fluid.heat_capacity_.at(N2);
    for (double cp : data.state.heat_capacity.var) EXPECT_DOUBLE_EQ(cp, cp_pure);
  }

  TEST_F(GasHeatCapacitySolverTest, ZeroMoleFractionComponentDoesNotChangeCp) {
    // The N2/O2 decks of the energy tests run pure N2 with O2 at zero.
    BuildInertBed({O2, N2}, {{O2, 0.0}, {N2, 1.0}});
    const double cp_pure = data.fluid.heat_capacity_.at(N2);
    for (double cp : data.state.heat_capacity.var) EXPECT_DOUBLE_EQ(cp, cp_pure);
  }

  // ---------------------------------------------------------------------------
  // The gas in the bed is a vapour, so a component that is condensed in its
  // pure state at (P, T) (water at 1 bar and 298 K, where the liquid has about
  // 4181 J/(kg.K)) has its vapour properties: the ideal-gas C_p0(T) of its
  // polynomial, and viscosity and conductivity fitted to dilute-gas values.

  namespace {
    const std::vector<std::string> kCondensables = {"Water", "Ethanol", "Toluene"};
  }

  // Fluid::Initialise caches the vapour value in the IDEAL_GAS mode, and leaves
  // the user-defined (CONSTANT) mode exactly as it was.
  TEST(GasHeatCapacityTest, InitialiseCachesVapourValueAndKeepsConstantMode) {
    InitialConditions ic;
    ic.pressure = 1.0e5;
    ic.temperature = 298.15;

    Fluid f;
    f.AddComponent("Water");
    f.AddComponent(N2);
    f.Initialise(ic);
    EXPECT_NEAR(f.heat_capacity_.at("Water"), 1865.0, 3.0);
    EXPECT_NEAR(f.MolarHeatCapacity("Water", kT), 33.59, 0.05);
    // The ideal-gas polynomial, within 0.2 % of CoolProp's real-gas value at 1
    // bar (recorded: the value earlier engines held).
    EXPECT_NE(f.HeatCapacityPolynomial(N2), nullptr);
    EXPECT_NEAR(f.heat_capacity_.at(N2), coolprop_reference::kNitrogenHeatCapacity, 2e-3 * f.heat_capacity_.at(N2));
    // Molar values a reader can check: ethanol 65.6 and toluene 104 J/(mol.K)
    // at 298 K (ideal gas).
    EXPECT_NEAR(Fluid::EvalIdealGasHeatCapacity(*Fluid::IdealGasCoeffs("Ethanol"), kT), 65.6, 1.5);
    EXPECT_NEAR(Fluid::EvalIdealGasHeatCapacity(*Fluid::IdealGasCoeffs("Toluene"), kT), 104.0, 2.5);

    f.gas_heat_capacity_mode_ = GasHeatCapacityMode::CONSTANT;
    f.gas_heat_capacity_["Water"] = 75.3;  // a user's value (here the liquid's) is kept as given
    f.Initialise(ic);
    EXPECT_DOUBLE_EQ(f.heat_capacity_.at("Water"), 75.3 / f.molar_mass_.at("Water"));
    EXPECT_DOUBLE_EQ(f.heat_capacity_.at(N2), Fluid::kDefaultGasHeatCapacity / f.molar_mass_.at(N2));
  }

  // Viscosity and conductivity come from polynomial fits of CoolProp's
  // DILUTE-GAS values, so condensables already get their vapour values (water
  // vapour at 298 K: about 9.9e-6 Pa.s and 0.0186 W/(m.K); the liquid has
  // 8.9e-4 Pa.s and 0.61 W/(m.K)).
  TEST(GasHeatCapacityTest, CondensableViscosityAndConductivityAreVapourValues) {
    EXPECT_NEAR(Fluid::GetComponentViscosity("Water", 298.15), 9.9e-6, 0.5e-6);
    EXPECT_NEAR(Fluid::GetComponentThermalConductivity("Water", 298.15), 0.0186, 0.001);
    for (const auto& c : kCondensables) {
      ASSERT_NE(Fluid::ViscosityCoeffs(c), nullptr) << c;
      ASSERT_NE(Fluid::ThermalConductivityCoeffs(c), nullptr) << c;
      const double mu = Fluid::GetComponentViscosity(c, 298.15);
      const double k = Fluid::GetComponentThermalConductivity(c, 298.15);
      EXPECT_GT(mu, 4e-6) << c;
      EXPECT_LT(mu, 1.5e-5) << c;  // liquids: 3e-4 to 1e-3 Pa.s
      EXPECT_GT(k, 0.005) << c;
      EXPECT_LT(k, 0.03) << c;     // liquids: 0.1 to 0.6 W/(m.K)
    }
  }

} // end adsim namespace
