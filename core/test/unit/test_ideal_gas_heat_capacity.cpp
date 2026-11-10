// Ideal-gas heat capacities C_p0(T) of the gas phase (ideal_gas_heat_capacity.h,
// NASA 7-coefficient polynomials of McBride, Gordon and Reno 1993, NASA
// TM-4513). These tests pin:
//   * the polynomials agree with CoolProp's ideal-gas C_p0 (the reference
//     equations of state; values recorded in coolprop_reference_values.h)
//     within a stated bound over 200-800 K;
//   * every engine component has a polynomial;
//   * the solver's cell heat capacity follows each cell's temperature;
//   * the enthalpy is the integral of C_p0 (Kirchhoff heats of reaction);
//   * the registry name of the mode, and its former name.

#include <cmath>
#include <cstdio>
#include <set>
#include <string>
#include <vector>

#include <gtest/gtest.h>

#include "../fixtures/user_data_test.h"
#include "coolprop_reference_values.h"
#include "adsim/equations/adsorbent_equations.h"
#include "adsim/equations/isotherm.h"
#include "adsim/objects/fluid.h"
#include "adsim/objects/reaction.h"
#include "adsim/parameters/parameter_registry.h"

namespace adsim {
namespace {

// The comparison window: the PSA/TSA range and beyond. Below a fit's t_min (or
// below the fluid's triple point, where CoolProp's ideal-gas part is outside
// its equation of state) the comparison starts higher.
constexpr double kCompareTmin = 200.0;
constexpr double kCompareTmax = 800.0;

// Bounds on |C_p0,poly / C_p0,CoolProp - 1| over the window: the largest
// measured error (cyclohexane, -1.8 % near 380 K; benzene and toluene -1.4 %)
// rounded up, and a tighter bound for the light and permanent gases of PSA.
constexpr double kMaxRelativeError = 0.02;
constexpr double kMaxRelativeErrorLightGases = 0.01;
const std::set<std::string> kLightGases = {"Air",    "Argon",    "CarbonDioxide", "CarbonMonoxide", "Helium",
                                           "Hydrogen", "Krypton", "Methane",       "Neon",           "Nitrogen",
                                           "Oxygen", "Xenon"};


}  // namespace

TEST(IdealGasHeatCapacity, EveryComponentHasAPolynomial) {
  for (const auto& c : Components::componentNames) EXPECT_NE(Fluid::IdealGasCoeffs(c), nullptr) << c;
  for (const auto& [c, p] : kIdealGasHeatCapacity) {
    EXPECT_NE(std::find(Components::componentNames.begin(), Components::componentNames.end(), c),
              Components::componentNames.end()) << c << " is not an engine component";
    EXPECT_LT(p.t_min, p.t_max) << c;
    EXPECT_DOUBLE_EQ(p.t_max, 1000.0) << c;  // one set, to 1000 K
    EXPECT_LE(p.t_min, 300.0) << c;  // covers room temperature
  }
}

TEST(IdealGasHeatCapacity, ProvenanceIsData) {
  EXPECT_FALSE(kIdealGasHeatCapacityReference.empty());
  for (const auto& [c, p] : kIdealGasHeatCapacity) {
    EXPECT_FALSE(p.species.empty()) << c;
    EXPECT_FALSE(p.source.empty()) << c;
  }
  EXPECT_EQ(Fluid::IdealGasCoeffs("Butane")->species, "C4H10,n-butane");
  EXPECT_EQ(Fluid::IdealGasCoeffs("Butane")->source, "L 6/90");
}

// The table of docs/equation-of-state.md comes from this test's output.
TEST(IdealGasHeatCapacity, AgreesWithCoolPropIdealGasHeatCapacity) {
  std::printf("| Component | Fit range [K] | Compared [K] | Max error | at [K] | C_p0(298.15 K) [J/(mol K)] |\n");
  for (const auto& c : Components::componentNames) {
    const NasaPolynomial* p = Fluid::IdealGasCoeffs(c);
    ASSERT_NE(p, nullptr) << c;
    // From 200 K, or from the fit's lower limit or the triple point if higher.
    double t_lo = std::max(kCompareTmin, p->t_min);
    const ComponentProperties* props = FindComponentProperties(c);
    ASSERT_NE(props, nullptr) << c;
    t_lo = std::max(t_lo, std::ceil(props->T_triple));
    const auto& ref_table = coolprop_reference::kIdealGasHeatCapacity.at(c);
    ASSERT_EQ(ref_table.t_start, t_lo) << c;  // the values were recorded on this grid
    double worst = 0.0, at = t_lo;
    size_t i = 0;
    for (double T = t_lo; T <= kCompareTmax + 1e-9; T += coolprop_reference::kCp0Step, ++i) {
      ASSERT_LT(i, ref_table.cp0.size()) << c;
      const double ref = ref_table.cp0[i];
      ASSERT_TRUE(std::isfinite(ref) && ref > 0.0) << c << " at " << T;
      const double err = Fluid::EvalIdealGasHeatCapacity(*p, T) / ref - 1.0;
      if (std::abs(err) > std::abs(worst)) { worst = err; at = T; }
    }
    const double bound = kLightGases.count(c) ? kMaxRelativeErrorLightGases : kMaxRelativeError;
    EXPECT_LE(std::abs(worst), bound) << c << ": " << 100 * worst << " % at " << at << " K";
    std::printf("| %s | %g-%g | %g-%g | %+.2f %% | %g | %.2f |\n", c.c_str(), p->t_min, p->t_max, t_lo,
                kCompareTmax, 100 * worst, at, Fluid::EvalIdealGasHeatCapacity(*p, 298.15));
  }
}

// Outside its range a polynomial is held at the nearer limit, and the
// enthalpy continues with that slope: dH/dT = C_p0 everywhere.
TEST(IdealGasHeatCapacity, ClampingAndEnthalpyAreConsistent) {
  const NasaPolynomial& so2 = *Fluid::IdealGasCoeffs("SulfurDioxide");  // t_min = 300 K
  ASSERT_DOUBLE_EQ(so2.t_min, 300.0);
  EXPECT_DOUBLE_EQ(Fluid::EvalIdealGasHeatCapacity(so2, 250.0), Fluid::EvalIdealGasHeatCapacity(so2, 300.0));
  const NasaPolynomial& co2 = *Fluid::IdealGasCoeffs("CarbonDioxide");
  EXPECT_DOUBLE_EQ(Fluid::EvalIdealGasHeatCapacity(co2, 1500.0), Fluid::EvalIdealGasHeatCapacity(co2, 1000.0));
  for (const auto* p : {&so2, &co2}) {
    for (double T : {150.0, 250.0, 298.15, 450.0, 999.0, 1001.0, 1500.0}) {
      const double h = 1e-3;
      const double dHdT = (Fluid::EvalIdealGasEnthalpy(*p, T + h) - Fluid::EvalIdealGasEnthalpy(*p, T - h)) / (2 * h);
      EXPECT_NEAR(dHdT, Fluid::EvalIdealGasHeatCapacity(*p, T), 2e-3 * dHdT) << T;
    }
  }
  // TM-4513's datum: H0(298.15 K) is the enthalpy of formation (CO2: -393.5 kJ/mol).
  EXPECT_NEAR(Fluid::EvalIdealGasEnthalpy(co2, 298.15), -393.51e3, 0.1e3);
  // N2 heats by integral C_p0 dT = 2.92 kJ/mol from 298.15 to 398.15 K
  // (NIST-JANAF H - H298 at 400 K: 2.971 kJ/mol; at 398.15: about 2.917).
  const NasaPolynomial& n2 = *Fluid::IdealGasCoeffs("Nitrogen");
  EXPECT_NEAR(Fluid::EvalIdealGasEnthalpy(n2, 398.15) - Fluid::EvalIdealGasEnthalpy(n2, 298.15), 2917.0, 10.0);
}

// Air, absent from TM-4513, is the N2/O2/Ar mixture of CoolProp's air model.
TEST(IdealGasHeatCapacity, AirIsTheMixtureOfItsConstituents) {
  const auto& air = *Fluid::IdealGasCoeffs("Air");
  for (double T : {200.0, 300.0, 500.0, 800.0}) {
    const double mix = 0.7812 * Fluid::EvalIdealGasHeatCapacity(*Fluid::IdealGasCoeffs("Nitrogen"), T) +
                       0.2096 * Fluid::EvalIdealGasHeatCapacity(*Fluid::IdealGasCoeffs("Oxygen"), T) +
                       0.0092 * Fluid::EvalIdealGasHeatCapacity(*Fluid::IdealGasCoeffs("Argon"), T);
    EXPECT_NEAR(Fluid::EvalIdealGasHeatCapacity(air, T), mix, 1e-9 * mix) << T;
  }
}

TEST(IdealGasHeatCapacity, FluidUsesThePolynomialOrTheConstantValue) {
  InitialConditions ic;
  ic.pressure = 1.0e5;
  ic.temperature = 320.0;
  Fluid f;
  f.AddComponent("Nitrogen");
  f.AddComponent("Water");
  f.Initialise(ic);
  ASSERT_NE(f.HeatCapacityPolynomial("Nitrogen"), nullptr);
  ASSERT_NE(f.HeatCapacityPolynomial("Water"), nullptr);

  // N2 follows the temperature.
  const auto& n2 = *Fluid::IdealGasCoeffs("Nitrogen");
  EXPECT_DOUBLE_EQ(f.MolarHeatCapacity("Nitrogen", 450.0), Fluid::EvalIdealGasHeatCapacity(n2, 450.0));
  EXPECT_GT(f.MolarHeatCapacity("Nitrogen", 450.0), f.MolarHeatCapacity("Nitrogen", 300.0));

  // The cached per-kilogram value of a polynomial component is its value at
  // the initial temperature (for display).
  EXPECT_DOUBLE_EQ(f.heat_capacity_.at("Nitrogen"),
                   Fluid::EvalIdealGasHeatCapacity(n2, 320.0) / f.molar_mass_.at("Nitrogen"));

  // CONSTANT mode: no polynomial; the constant at every temperature.
  f.gas_heat_capacity_mode_ = GasHeatCapacityMode::CONSTANT;
  f.Initialise(ic);
  EXPECT_EQ(f.HeatCapacityPolynomial("Nitrogen"), nullptr);
  EXPECT_DOUBLE_EQ(f.MolarHeatCapacity("Nitrogen", 450.0), Fluid::kDefaultGasHeatCapacity);
  EXPECT_DOUBLE_EQ(f.MolarEnthalpyChange("Nitrogen", 300.0, 310.0), 10.0 * Fluid::kDefaultGasHeatCapacity);

  // GAS_PHASE adsorbed heat capacity follows the adsorbed-phase temperature.
  f.gas_heat_capacity_mode_ = GasHeatCapacityMode::IDEAL_GAS;
  f.adsorbed_heat_capacity_mode_ = AdsorbedHeatCapacityMode::GAS_PHASE;
  f.Initialise(ic);
  EXPECT_DOUBLE_EQ(f.AdsorbedPhaseHeatCapacity("Nitrogen", 400.0), Fluid::EvalIdealGasHeatCapacity(n2, 400.0));
  EXPECT_EQ(f.AdsorbedPhaseHeatCapacityPolynomial("Nitrogen"), &n2);
}

// Kirchhoff heat of reaction: dH_ref + sum_i nu_i integral C_p,i dT, exactly.
TEST(IdealGasHeatCapacity, KirchhoffHeatOfReactionIntegratesTheHeatCapacity) {
  InitialConditions ic;
  ic.pressure = 1.0e5;
  ic.temperature = 298.15;
  Fluid f;
  for (const char* c : {"Methane", "Water", "CarbonMonoxide", "Hydrogen"}) f.AddComponent(c);
  f.Initialise(ic);
  Reaction r;
  r.nu_ = {{"Methane", -1.0}, {"Water", -1.0}, {"CarbonMonoxide", 1.0}, {"Hydrogen", 3.0}};
  r.heat_of_reaction_ = 206.1e3;
  r.heat_reference_temperature_ = 298.15;
  r.heat_mode_ = HeatOfReactionMode::Kirchhoff;
  double expected = 206.1e3;
  for (const auto& [c, nu] : r.nu_) {
    const auto& p = *Fluid::IdealGasCoeffs(c);
    expected += nu * (Fluid::EvalIdealGasEnthalpy(p, 800.0) - Fluid::EvalIdealGasEnthalpy(p, 298.15));
  }
  EXPECT_NEAR(r.HeatOfReaction(800.0, f), expected, 1e-6);
  // The polynomials' own enthalpies of formation give the reforming heat at
  // 298.15 K that the reaction states (206.1 kJ/mol), so the datum is right.
  double dh298 = 0.0;
  for (const auto& [c, nu] : r.nu_) dh298 += nu * Fluid::EvalIdealGasEnthalpy(*Fluid::IdealGasCoeffs(c), 298.15);
  EXPECT_NEAR(dh298, 206.1e3, 0.5e3);
  EXPECT_GT(r.HeatOfReaction(800.0, f), 206.1e3);  // dCp > 0 for reforming
  EXPECT_NEAR(r.HeatOfReaction(298.15, f), 206.1e3, 1e-6);
}

TEST(IdealGasHeatCapacity, RegistryNamesTheModeAndAcceptsItsFormerName) {
  UserData data;
  data.fluid.AddComponent("Nitrogen");
  data.fluid.gas_heat_capacity_mode_ = GasHeatCapacityMode::CONSTANT;
  ParameterRegistry reg(data);
  Parameter& p = reg.Get("fluid.gas_heat_capacity_mode");
  EXPECT_EQ(p.Info().choices, (std::vector<std::string>{"ideal_gas", "constant"}));
  p.SetText("ideal_gas");
  EXPECT_EQ(data.fluid.gas_heat_capacity_mode_, GasHeatCapacityMode::IDEAL_GAS);
  p.SetText("constant");
  p.SetText("coolprop");  // the former name of ideal_gas
  EXPECT_EQ(data.fluid.gas_heat_capacity_mode_, GasHeatCapacityMode::IDEAL_GAS);
  EXPECT_THROW(p.SetText("refprop"), std::invalid_argument);
}

// The solver's cell heat capacity is the molar mixture of each component's
// C_p0 at that cell's gas temperature.
class IdealGasHeatCapacitySolverTest : public UserDataTest {};

TEST_F(IdealGasHeatCapacitySolverTest, CellHeatCapacityFollowsTheCellTemperature) {
  const std::vector<std::string> comps = {"CarbonDioxide", "Nitrogen"};
  AddComponents(comps);
  auto [layer_id, ads_id] = AddLayer("Layer1");
  GetLayer(layer_id).domain_.n_cells = 6;
  auto henry = std::make_unique<HenryIsotherm>();
  for (const auto& c : comps) {
    henry->params_[c].k = 0.0;
    henry->params_[c].Q = 0.0;
  }
  AddIsotherm(comps, ads_id, std::move(henry));
  for (const auto& c : comps) {
    AddLinearDrivingForceModel(ads_id, c, 1.0);
    SetHOA(ads_id, c, 0.0);
  }
  AddConstantInternalHeatTransferCoefficientModel(10.0, layer_id);
  AddConstantExternalHeatTransferCoefficientModel(10.0, layer_id);
  std::get<CylinderBed>(data.bed.shape_).diameter = 0.05;
  AddInitialConditions(1.0e5, 298.15, 298.15, {{"CarbonDioxide", 0.15}, {"Nitrogen", 0.85}});
  data.Initialise();

  auto& T = data.state.fluid_temperature.var;
  for (size_t n = 0; n < T.size(); ++n) T[n] = 280.0 + 40.0 * n;  // 280 K ... 520 K
  UpdateHeatCapacity(&data);
  const auto& Cp = data.state.heat_capacity.var;
  const std::unordered_map<std::string, double> y{{"CarbonDioxide", 0.15}, {"Nitrogen", 0.85}};
  for (size_t n = 1; n + 1 < Cp.size(); ++n)
    EXPECT_NEAR(Cp[n], data.fluid.MixtureHeatCapacity(y, T[n]), 1e-12 * Cp[n]) << "cell " << n;
  EXPECT_GT(Cp[Cp.size() - 2], 1.05 * Cp[1]);  // CO2/N2 C_p rises by > 5 % from 320 to 480 K
}

}  // namespace adsim
