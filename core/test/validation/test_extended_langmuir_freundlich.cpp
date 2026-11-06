/*
 * Skarstrom - pressure and temperature swing adsorption simulator
 * Copyright (c) 2024-2026 Sam Affleck and contributors
 *
 * Licensed under the terms in the LICENSE file at the repository root.
 */

// The EXTENDED LANGMUIR-FREUNDLICH isotherm with temperature-dependent
// parameters (Yang, "Gas Separation by Adsorption Processes", 1987; the form
// of toPSAil's calcIsothermExtLangFreu.m):
//
//     q_i  = q_s,i z_i / (1 + sum_j z_j),   z_i = b_i (P y_i / P_ref)^{n_i},
//     q_s,i = k1 + k2 T,   b_i = k3 exp(k4 / T),   1 / n_i = k5 + k6 / T.
//
// The reference below is written from that algebra and never calls the
// solver's evaluation. Parameters: H2 and CO on activated carbon as in
// toPSAil's case_study_4.0 (k3 in bar^-n, so P_ref = 1e5 Pa).
//
//   1. Competitive loadings seeded in a closed uniform bed (the solver's
//      residual path) match the closed form at the case's feed state and at a
//      second temperature, to the seeding tolerance. Guards: the shared
//      denominator suppresses the CO loading well below its single-component
//      value, and the T-dependence of q_s (k2) and n (k6), which no existing
//      kind can represent, moves the loading by far more than the tolerance.
//   2. The pure-component evaluation (EvalIsotherm, used by the Python
//      evaluate_isotherm and the front ends) and the IAST pure-isotherm
//      mapping (a Sips term with affinity b^{1/n} / P_ref and exponent n)
//      match the single-component closed form to round-off.
//   3. Trace components keep their loading: the power is continued below
//      y0 = 1e-5 as in the competitive Sips kinds (docs/isotherms.md), so the
//      pure-component form is exact down to y0 and exactly 0 at y = 0. The
//      competitive path is checked in test_competitive_sips_trace_loading.cpp.

#include "validation_helpers.h"

#include "adsim/equations/iast.h"
#include "adsim/equations/isotherm_kinds.h"

namespace adsim {

  using namespace adsim::validation;

  namespace {

    struct Elf { double k1, k2, k3, k4, k5, k6; };
    constexpr double kPref = 1.0e5;

    // case_study_4.0, 1.3_adsorbent_properties.xlsm (H2, CO).
    const Elf kH2{4.31, -1.06e-2, 2.4821120157907722e-3, 458.2, 0.986, 43.03};
    const Elf kCO{5.05, -9.05e-3, 1.1221317542561065e-3, 1617.0, 0.5245, 256.5};

    double Z(const Elf& p, double T, double P_i) {
      const double n = 1.0 / (p.k5 + p.k6 / T);
      return p.k3 * std::exp(p.k4 / T) * std::pow(P_i / kPref, n);
    }
    double Qs(const Elf& p, double T) { return p.k1 + p.k2 * T; }

    // Competitive loadings of (H2, CO) at total pressure P and CO fraction y_co.
    std::pair<double, double> Competitive(double T, double P, double y_co) {
      const double zH = Z(kH2, T, P * (1.0 - y_co));
      const double zC = Z(kCO, T, P * y_co);
      const double d = 1.0 + zH + zC;
      return {Qs(kH2, T) * zH / d, Qs(kCO, T) * zC / d};
    }

    double Pure(const Elf& p, double T, double P) {
      const double z = Z(p, T, P);
      return Qs(p, T) * z / (1.0 + z);
    }

    void SetParams(ExtendedLangmuirFreundlichIsothermParameters& to, const Elf& p) {
      to.k1 = p.k1; to.k2 = p.k2; to.k3 = p.k3; to.k4 = p.k4; to.k5 = p.k5; to.k6 = p.k6;
      to.P_ref = kPref;
    }

    std::unique_ptr<ExtendedLangmuirFreundlichIsotherm> MakeIsotherm() {
      auto iso = std::make_unique<ExtendedLangmuirFreundlichIsotherm>();
      SetParams(iso->params_["Hydrogen"], kH2);
      SetParams(iso->params_["CarbonMonoxide"], kCO);
      return iso;
    }

  } // namespace

  class ExtendedLangmuirFreundlichTest : public UserDataTest {
  public:
    // Closed, uniform, isothermal bed; returns the loadings the solver seeds at
    // t = 0 (q = q* of the initial gas), mean over the interior.
    std::pair<double, double> SeededLoadings(double T, double P, double y_co) {
      AddComponents({"Hydrogen", "CarbonMonoxide"});
      auto [layer_id, ads_id] = AddLayer("Carbon");
      auto& layer = GetLayer(layer_id);
      layer.domain_.length = 1.0;
      layer.domain_.n_cells = 5;
      layer.eb_ = 0.4;
      layer.particle_diameter_ = 0.002;
      std::get<CylinderBed>(data.bed.shape_).diameter = 0.05;
      auto& adsorbent = GetAdsorbent(ads_id);
      adsorbent.ep_ = 0.0;
      adsorbent.solid_density_ = 1000.0;
      adsorbent.heat_capacity_ = 900.0;
      adsorbent.thermal_conductivity_ = 0.0;
      AddConstantInternalHeatTransferCoefficientModel(0.0, layer_id);
      AddConstantExternalHeatTransferCoefficientModel(0.0, layer_id);
      adsorbent.adsorption_source_.AddIsotherm({"Hydrogen", "CarbonMonoxide"}, MakeIsotherm());
      for (const char* c : {"Hydrogen", "CarbonMonoxide"}) {
        AddLinearDrivingForceModel(ads_id, c, 1.0);
        SetHOA(ads_id, c, 0.0);
      }
      data.bed.equations_.energy_balance_mode_ = EnergyBalanceMode::ISOTHERMAL;
      AddInitialConditions(P, T, T, {{"Hydrogen", 1.0 - y_co}, {"CarbonMonoxide", y_co}}, 0.0);
      AddStep("Hold", 1.0);
      AddClosedBC(BoundaryLocation::BOTTOM, "Hold");
      AddClosedBC(BoundaryLocation::TOP, "Hold");
      ConfigureCycle(1);
      data.settings.reltol = 1e-10;

      CoupledSolver solver{data};
      std::string error;
      EXPECT_TRUE(data.PerformChecks(error)) << error;
      EXPECT_GE(solver.IntegrateCycle(data.cycle, false), 0) << "Solver did not complete";
      const auto& log = solver.GetLog();
      auto seeded = [&](const char* c) {
        const auto q = log.GetComponentVariable(0, 0, ComponentVariableType::ADSORBED_CONCENTRATION, c);
        EXPECT_FALSE(q.values.empty());
        return q.values.empty() ? 0.0 : InteriorMean(q.values.front(), log.GetDomain().volume);
      };
      return {seeded("Hydrogen"), seeded("CarbonMonoxide")};
    }
  };

  TEST_F(ExtendedLangmuirFreundlichTest, CompetitiveLoadingsMatchClosedFormAt303K) {
    const double T = 303.0, P = 11.14575e5, y = 0.3;  // case_study_4.0 feed
    const auto [qH, qC] = SeededLoadings(T, P, y);
    const auto [qH_ref, qC_ref] = Competitive(T, P, y);
    std::cout << "[ELF] 303 K: H2 " << qH << " (ref " << qH_ref << "), CO " << qC << " (ref "
              << qC_ref << ") mol/kg\n";
    // GUARD: genuinely competitive: CO takes a third of the H2 loading H2 would
    // have alone at its partial pressure (closed form 0.66 of it).
    ASSERT_LT(qH_ref, 0.8 * Pure(kH2, T, P * (1.0 - y)));
    EXPECT_NEAR(qH, qH_ref, 1e-6 * qH_ref);
    EXPECT_NEAR(qC, qC_ref, 1e-6 * qC_ref);
  }

  TEST_F(ExtendedLangmuirFreundlichTest, CompetitiveLoadingsMatchClosedFormAt343K) {
    const double T = 343.0, P = 5.0e5, y = 0.3;
    const auto [qH, qC] = SeededLoadings(T, P, y);
    const auto [qH_ref, qC_ref] = Competitive(T, P, y);
    std::cout << "[ELF] 343 K: H2 " << qH << " (ref " << qH_ref << "), CO " << qC << " (ref "
              << qC_ref << ") mol/kg\n";
    // GUARD: with q_s and 1/n frozen at their 303 K values (k2 = k6 = 0, the
    // only case extended_sips_2 can represent) the CO loading at 343 K is off
    // by more than 1 %, so the temperature dependence is what is being matched.
    const double T0 = 303.0;
    const Elf co_frozen{kCO.k1 + kCO.k2 * T0, 0.0, kCO.k3, kCO.k4, kCO.k5 + kCO.k6 / T0, 0.0};
    const Elf h2_frozen{kH2.k1 + kH2.k2 * T0, 0.0, kH2.k3, kH2.k4, kH2.k5 + kH2.k6 / T0, 0.0};
    const double zH = Z(h2_frozen, T, P * (1 - y)), zC = Z(co_frozen, T, P * y);
    const double qC_frozen = Qs(co_frozen, T) * zC / (1 + zH + zC);
    ASSERT_GT(std::abs(qC_frozen / qC_ref - 1.0), 0.01);
    EXPECT_NEAR(qH, qH_ref, 1e-6 * qH_ref);
    EXPECT_NEAR(qC, qC_ref, 1e-6 * qC_ref);
  }

  TEST(ExtendedLangmuirFreundlich, PureEvaluationAndIastMappingMatchClosedForm) {
    auto iso = MakeIsotherm();
    IsothermType kind{};
    ASSERT_TRUE(KindFromString("extended_langmuir_freundlich", kind));
    EXPECT_EQ(kind, IsothermType::EXTENDED_LANGMUIR_FREUNDLICH);
    EXPECT_TRUE(IsCompetitiveIsothermKind(kind));
    EXPECT_EQ(IsothermFields(iso.get(), "CarbonMonoxide").size(), 7u);
    double worst = 0.0, worst_iast = 0.0;
    for (const auto* p : {&kH2, &kCO}) {
      const char* c = (p == &kH2) ? "Hydrogen" : "CarbonMonoxide";
      for (double T : {273.15, 303.0, 373.15}) {
        const PureIsotherm pure = PreparePureIsotherm(iso.get(), c, T);
        for (double P : {1.0e3, 1.0e5, 2.0e6}) {
          const double ref = Pure(*p, T, P);
          worst = std::max(worst, std::abs(EvalIsotherm(iso.get(), c, T, P, 1.0) / ref - 1.0));
          worst_iast = std::max(worst_iast, std::abs(PureLoading(pure, P) / ref - 1.0));
        }
      }
    }
    std::cout << "[ELF] pure: max rel error " << worst << ", IAST mapping " << worst_iast << "\n";
    // Both evaluate the published power at y = 1 (round-off only).
    EXPECT_LT(worst, 1e-12);
    EXPECT_LT(worst_iast, 1e-12);
  }

  TEST(ExtendedLangmuirFreundlich, TraceComponentKeepsItsLoading) {
    auto iso = MakeIsotherm();
    // Pure-component form of 11 bar down to y0 = 1e-5 against the published one.
    // GUARD: the former step y^2 / (y^2 + 1e-10) moved it by 1e-4 at y = 1e-3,
    // 1 % at 1e-4 and 20 % at 2e-5, far outside the round-off tolerance.
    const double T = 303.0, P = 11.0e5;
    for (double y : {1e-2, 1e-3, 1e-4, 2e-5, 1e-5}) {
      const double q = EvalIsotherm(iso.get(), "CarbonMonoxide", T, P, y);
      const double ref = Pure(kCO, T, P * y);
      std::cout << "[ELF] y_CO = " << y << ": q " << q << " ref " << ref << "\n";
      EXPECT_NEAR(q, ref, 1e-12 * ref) << "y = " << y;
    }
    // ... exactly zero loading at zero mole fraction, increasing just above it.
    EXPECT_EQ(EvalIsotherm(iso.get(), "CarbonMonoxide", T, P, 0.0), 0.0);
    EXPECT_GT(EvalIsotherm(iso.get(), "CarbonMonoxide", T, P, 1e-7), 0.0);
  }

} // end adsim namespace
