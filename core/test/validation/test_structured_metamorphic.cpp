/*
 * Skarstrom - pressure and temperature swing adsorption simulator
 * Copyright (c) 2024-2026 Sam Affleck and contributors
 *
 * Licensed under the terms in the LICENSE file at the repository root.
 */

// Metamorphic relation between a packed layer and a structured one: the same
// sorbent mass, gas holdup, transfer rates and flow resistance must give the
// same breakthrough, whatever holds the sorbent.
//
// A packed layer (Kozeny-Carman, bed void eb) and a laminate or square
// monolith are matched on everything the balances read:
//   - gas holdup and sorbent mass: OFA = eb (and the same solid, ep = 0), so
//     the bulk density (1 - OFA) rho_s is the same;
//   - axial dispersion and axial heat conduction: the `constant` models with
//     the same values (Taylor-Aris and Ruthven would bring in D_h and d_p);
//   - uptake: the LDF with the same constant k (the slab and sphere forms of
//     the macropore model differ by design);
//   - flow: the channel geometry is chosen so that its laminar permeability,
//     K = OFA D_h^2 / (2 fRe), equals the packed bed's Kozeny-Carman
//     K = eb^3 d_p^2 / (180 (1 - eb)^2). Both are then Darcy flow with the
//     same K, so the pressure profile is the same too. (With unmatched K
//     the pressure drop, and through it the gas density and velocity, would
//     differ: a legitimate difference of the models, not of the sorbent.)
//   - heat: isothermal (a non-adsorbing heat balance) in the first relation;
//     in the second, LTE with a heat of adsorption and a bed-to-wall
//     coefficient, which treats both kinds alike (under LTNE they differ on
//     purpose: the wall flux of a structured layer goes into the solid, and
//     the gas-solid film area is the channel wall area).
// The two runs then solve the same equations; the outlet curves must agree to
// the integration tolerance. The runs take their own time steps, so one curve
// is read at the other's logged times by linear interpolation, whose error
// (dt^2 / 8) y'' dominates the comparison: on the isothermal front it was
// 1e-4 of y_feed with a 0.5 s step cap and 1.6e-6 with 0.05 s (falling with
// dt^2, and not with the integrator tolerance), hence the 0.05 s cap and the
// 1e-5 bound on the pointwise difference. The mean breakthrough time, an
// integral that the interpolation does not enter, must agree to 1e-6.

#include "breakthrough_case.h"

#include <algorithm>
#include <cmath>
#include <string>
#include <vector>

namespace adsim {

  using namespace adsim::validation;

  namespace {

    constexpr double kPackedDp = 0.002;  // breakthrough_case.h

    double KozenyCarmanK(double eb, double dp) {
      return eb * eb * eb * dp * dp / (180.0 * (1.0 - eb) * (1.0 - eb));
    }

    // Turn the case's packed layer into a structured one with OFA = eb and
    // the same permeability.
    void MakeStructured(MixedAdsorbentLayer& layer, LayerKind kind) {
      const double eb = layer.eb_.Get();
      const double K = KozenyCarmanK(eb, kPackedDp);
      layer.kind_ = kind;
      layer.structured_.washcoat_thickness_ = 0.0;
      if (kind == LayerKind::Laminate) {
        // OFA = g / (g + t), D_h = 2 g, fRe = 24: K = OFA (2 g)^2 / 48.
        const double g = std::sqrt(12.0 * K / eb);
        layer.structured_.gap_ = g;
        layer.structured_.sheet_thickness_ = g * (1.0 - eb) / eb;
      } else {
        // Square: OFA = a^2 / p^2, D_h = a, fRe = 14.23: K = OFA a^2 / 28.46.
        layer.structured_.channel_shape_ = ChannelShape::Square;
        const double a = std::sqrt(2.0 * 14.23 * K / eb);
        const double p = a / std::sqrt(eb);
        layer.structured_.cell_density_ = 1.0 / (p * p);
        layer.structured_.wall_thickness_ = p - a;
      }
    }

    struct Curves {
      std::vector<double> t, y, dp_bed;
    };

    Curves Run(const BreakthroughCase& c, LayerKind kind, bool thermal) {
      UserData data;
      Build(data, c);
      auto& layer = data.bed.layers_.layers_.front();
      auto disp = std::make_unique<DispersionModelConstant>();
      disp->dispersion_coefficient_ = 2.0e-5;  // [m2/s]: a real spread, the same in both
      layer.equations_.dispersion_model_ = std::move(disp);
      if (thermal) {
        auto& ads = data.adsorbent_directory.Get(layer.composition_.begin()->first);
        auto iso = std::make_unique<LangmuirIsotherm>();
        iso->params_["Oxygen"].m1 = 3.0;      // [mol/kg]
        iso->params_["Oxygen"].b1 = 1.0e-9;   // [1/Pa] at the reference
        iso->params_["Oxygen"].q1 = 20000.0;  // [J/mol]
        ads.adsorption_source_.AddIsotherm({"Oxygen"}, std::move(iso));
        ads.adsorption_source_.heat_of_adsorption_["Oxygen"] = 20000.0;
        ads.thermal_conductivity_ = 0.5;
        auto hi = std::make_unique<ConstantInternalHeatTransferModel>();
        hi->heat_transfer_coefficient_ = 30.0;
        layer.equations_.internal_heat_transfer_model_ = std::move(hi);
        auto ho = std::make_unique<ConstantExternalHeatTransferModel>();
        ho->heat_transfer_coefficient_ = 10.0;
        layer.equations_.external_heat_transfer_model_ = std::move(ho);
        auto lam = std::make_unique<HeatDispersionModelConstant>();
        lam->dispersion_coefficient_ = 0.2;  // [W/(m K)]
        layer.equations_.heat_dispersion_model_ = std::move(lam);
        data.bed.ambient_temperature_ = c.bed_temperature;
      } else {
        data.bed.equations_.energy_balance_mode_ = EnergyBalanceMode::ISOTHERMAL;
        // A Henry constant that holds the front back (retention factor
        // rho_b k R T / eb ~ 5), so there is a front to compare.
        auto& ads = data.adsorbent_directory.Get(layer.composition_.begin()->first);
        auto henry = std::make_unique<HenryIsotherm>();
        henry->params_["Oxygen"].k = 1.35e-6;  // [mol/(kg Pa)]
        henry->params_["Oxygen"].Q = 0.0;
        ads.adsorption_source_.AddIsotherm({"Oxygen"}, std::move(henry));
      }
      if (kind != LayerKind::Packed) MakeStructured(layer, kind);

      std::string err;
      EXPECT_TRUE(data.EvaluateAllExpressions(err) && data.PerformChecks(err)) << err;
      CoupledSolver solver{data};
      EXPECT_GE(solver.IntegrateCycle(data.cycle, false), 0) << "Solver did not complete";
      Curves out;
      auto [t, y] = OutletSeries(solver.GetLog(), false);
      out.t = t;
      out.y = y;
      const auto P = solver.GetLog().GetVariable(0, 0, VariableType::PRESSURE);
      for (const auto& row : P.values) out.dp_bed.push_back(row[1] - row[row.size() - 2]);
      return out;
    }

    // Largest difference of b from a, b read at a's times.
    double WorstGap(const std::vector<double>& ta, const std::vector<double>& a,
                    const std::vector<double>& tb, const std::vector<double>& b) {
      double worst = 0.0;
      for (size_t n = 0; n < ta.size(); ++n) {
        if (ta[n] < tb.front() || ta[n] > tb.back()) continue;
        worst = std::max(worst, std::abs(a[n] - InterpolateSeries(tb, b, ta[n])));
      }
      return worst;
    }

    void Compare(const BreakthroughCase& c, bool thermal, double tolerance) {
      const Curves packed = Run(c, LayerKind::Packed, thermal);
      ASSERT_GT(packed.y.size(), 10u);
      // The breakthrough must happen inside the run, or nothing is compared.
      ASSERT_LT(InterpolateSeries(packed.t, packed.y, packed.t.front() + 0.1 * c.run_time), 0.1 * c.y_feed);
      ASSERT_GT(packed.y.back(), 0.9 * c.y_feed);
      for (LayerKind kind : {LayerKind::Laminate, LayerKind::Monolith}) {
        const Curves s = Run(c, kind, thermal);
        const double gap = WorstGap(packed.t, packed.y, s.t, s.y) / c.y_feed;
        const double gap_p = WorstGap(packed.t, packed.dp_bed, s.t, s.dp_bed) /
                             std::max(1e-12, std::abs(packed.dp_bed.back()));
        std::cout << "[packed vs " << LayerKindList[static_cast<int>(kind)] << (thermal ? ", LTE" : ", isothermal")
                  << "] worst outlet difference " << gap << " of y_feed; worst bed pressure-drop difference " << gap_p
                  << " of " << packed.dp_bed.back() << " Pa\n";
        EXPECT_LT(gap, tolerance) << LayerKindList[static_cast<int>(kind)];
        // The first moment (mean breakthrough time) of each run on its own
        // grid: t1 = int (1 - y / y_feed) dt.
        auto t1 = [&](const Curves& r) {
          std::vector<double> f(r.y.size());
          for (size_t n = 0; n < f.size(); ++n) f[n] = 1.0 - r.y[n] / c.y_feed;
          return TimeIntegral(r.t, f);
        };
        const double t1p = t1(packed), t1s = t1(s);
        std::cout << "    first moment " << t1s << " s against " << t1p << " s\n";
        EXPECT_NEAR(t1s, t1p, 1e-6 * t1p) << LayerKindList[static_cast<int>(kind)];
        EXPECT_LT(gap_p, 1e-4) << LayerKindList[static_cast<int>(kind)];
      }
    }

  }  // namespace

  // The geometry matching itself: same OFA, same permeability, same bulk
  // density as the packed layer.
  TEST(StructuredMetamorphic, MatchedLayersShareVoidMassAndPermeability) {
    MixedAdsorbentLayer packed;
    packed.eb_ = kCaseVoidFraction;
    packed.particle_diameter_ = kPackedDp;
    packed.momentum_model_ = MomentumModel::KarmenCozeny;
    for (LayerKind kind : {LayerKind::Laminate, LayerKind::Monolith}) {
      MixedAdsorbentLayer s;
      s.eb_ = kCaseVoidFraction;
      MakeStructured(s, kind);
      EXPECT_NEAR(s.GetVoidFraction(), packed.GetVoidFraction(), 1e-14);
      EXPECT_NEAR(s.GetPermeability(), packed.GetPermeability(), 1e-12 * packed.GetPermeability());
      EXPECT_EQ(s.GetInertialCoefficient(), 0.0);
      EXPECT_EQ(packed.GetInertialCoefficient(), 0.0);
    }
  }

  // Isothermal linear breakthrough: packed, laminate and monolith agree.
  TEST(StructuredMetamorphic, SameSorbentAndRatesGiveTheSameIsothermalBreakthrough) {
    BreakthroughCase c;
    c.adsorbing = true;
    c.n_cells = 30;
    c.length = 0.3;
    c.molar_flow = 0.01;
    c.run_time = 100.0;
    c.reltol = 1e-9; c.max_step = 0.05;
    Compare(c, false, 1e-5);
  }

  // Non-isothermal (LTE) Langmuir breakthrough with a heat of adsorption and
  // wall exchange: the relation still holds.
  TEST(StructuredMetamorphic, SameSorbentAndRatesGiveTheSameLteBreakthrough) {
    BreakthroughCase c;
    c.adsorbing = true;
    c.n_cells = 30;
    c.length = 0.3;
    c.molar_flow = 0.01;
    c.run_time = 300.0;
    c.reltol = 1e-9; c.max_step = 0.05;
    Compare(c, true, 1e-5);
  }

}  // namespace adsim
