/*
 * Skarstrom - pressure and temperature swing adsorption simulator
 * Copyright (c) 2024-2026 Sam Affleck and contributors
 *
 * Licensed under the terms in the LICENSE file at the repository root.
 */

// Analytical validation of a Langmuir DESORPTION front: the equilibrium-theory
// simple (spreading) wave, by the method of characteristics.
//
// test_equilibrium_shock pins the favourable direction: a clean bed fed with
// adsorbate, where the Langmuir front self-sharpens into a shock whose speed
// is set by the isotherm SECANT. The reverse step -- a saturated bed purged
// with clean carrier -- is unfavourable for a Langmuir isotherm, and the
// front spreads into a simple wave in which every concentration travels at
// its own speed, set by the isotherm SLOPE at that concentration. That makes
// the whole outlet curve a probe of q*'(c) across the isotherm.
//
// Theory (de Vault 1943; Rhee, Aris & Amundson 1970; Ruthven 1984, ch. 8).
// Trace adsorbate in an inert carrier at constant pressure and temperature,
// interstitial velocity u constant, local equilibrium, no dispersion:
//
//     dc/dt + u dc/dz + (rho_b/eps) dq*/dt = 0.
//
// Each concentration travels along a characteristic at
//
//     w(c) = u / (1 + (rho_b/eps) f'(c)),   f(c) = q*(c),
//
// and since every characteristic of the wave starts at the inlet at t = 0, the
// outlet sees concentration c at t(c) = (L/u) (1 + (rho_b/eps) f'(c)). For the
// Langmuir isotherm q* = m K p / (1 + K p), p = c R T,
//
//     (rho_b/eps) f'(c) = M / (1 + K p)^2,    M = (rho_b/eps) m K R T,
//
// so the outlet partial pressure on the wave is
//
//     K p(t) = sqrt( M / (u t / L - 1) ) - 1,   t(c0) <= t <= t(0),
//
// with the plateau p = p0 before t(c0) = (L/u)(1 + M / (1 + K p0)^2) and the
// clean carrier after t(0) = (L/u)(1 + M).
//
// The solver resolves the wave with first-order upwind numerical diffusion
// D = u dx / 2 and a finite LDF rate. Both smooth the two corners of the fan
// and shift the levels inside it by O(dx) and O(1/k t); the fan itself is
// exact only in the limit. So the test asserts that the outlet time of three
// levels inside the fan (c/c0 = 0.6, 0.5, 0.4) converges at FIRST order in dx
// to the characteristic time t(c), and that the Richardson-extrapolated time
// (2 t_fine - t_coarse) agrees with t(c) to 0.1 %. The extrapolated residual
// is the sum of the neglected terms:
//   * the LDF lag, O(1/(k t)) ~ 6e-4 at k = 20 1/s;
//   * the sorption effect (the desorbed gas speeds up the carrier), O(y0):
//     measured as 1.3e-3 of t(c) at y0 = 1e-3 and 1.3e-4 at y0 = 1e-4, so
//     y0 = 1e-4 here;
//   * the solver's absolute tolerance on mole fractions (1e-8, i.e. 1e-4 of
//     y0), which moves a level by |dt/dphi| * 1e-4, about 2e-4 of t(c) at
//     c/c0 = 0.4. Lower levels, where the wave is flatter, are limited by
//     this tolerance rather than by the grid, so the test stops at 0.4.
//
// GUARDS: the isotherm is strongly nonlinear across the step (K p0 = 3, so
// f'(c) varies 16-fold over the wave), and the wave spreads over a time
// (t(0) - t(c0) ~ 450 s) far wider than the numerical smear, so a shock or a
// linear-isotherm front (one outlet time for every level) cannot pass.

#include "breakthrough_case.h"

#include <algorithm>

namespace adsim {

  using namespace adsim::validation;

  namespace {

    constexpr double kT = 298.0;
    constexpr double kY0 = 1.0e-4;        // initial adsorbate level (trace)
    constexpr double kKp0 = 3.0;          // K p0: well up the isotherm's knee
    constexpr double kM = 19.0;           // (rho_b/eps) m K R T: R(0) = 20
    constexpr double kLdf = 20.0;         // [1/s], fast
    constexpr double kRunTime = 700.0;    // [s]

    double LangmuirK() { return kKp0 / (kCasePressure * kY0); }  // [1/Pa]
    double Capacity() {                                          // m [mol/kg]
      const double rho_b = (1.0 - kCaseVoidFraction) * kCaseSolidDensity;
      return kM / ((rho_b / kCaseVoidFraction) * LangmuirK() * kR * kT);
    }

    // Characteristic outlet time of the level c/c0 = phi.
    double CharacteristicTime(double phi, double u, double L) {
      const double Kp = kKp0 * phi;
      return (L / u) * (1.0 + kM / ((1.0 + Kp) * (1.0 + Kp)));
    }

    struct WaveRun {
      std::vector<double> t, phi;  // outlet y / y0
      double u = 0.0, dP_rel = 0.0;
    };

    WaveRun RunDesorption(int n_cells) {
      BreakthroughCase c;
      c.n_cells = n_cells;
      c.adsorbing = false;   // the Langmuir isotherm is installed below
      c.y_feed = 0.0;        // purge with clean carrier
      c.run_time = kRunTime;
      c.max_step = 0.5;
      c.reltol = 1e-9;
      c.bed_temperature = kT;
      c.feed_temperature = kT;

      UserData data;
      Build(data, c);
      data.settings.max_log_points_per_step_ = 700;
      data.bed.equations_.energy_balance_mode_ = EnergyBalanceMode::ISOTHERMAL;
      // The bed starts saturated: the solver seeds q at equilibrium with the
      // initial gas.
      data.ic.mole_fraction["Nitrogen"] = 1.0 - kY0;
      data.ic.mole_fraction["Oxygen"] = kY0;

      auto& layer = data.bed.layers_.layers_.front();
      auto& adsorbent = data.adsorbent_directory.Get(layer.composition_.begin()->first);
      auto iso = std::make_unique<LangmuirIsotherm>();
      iso->params_["Oxygen"].m1 = Capacity();
      iso->params_["Oxygen"].b1 = LangmuirK();
      iso->params_["Oxygen"].q1 = 0.0;
      adsorbent.adsorption_source_.AddIsotherm({"Oxygen"}, std::move(iso));
      auto ldf = std::make_unique<LinearDrivingForce>();
      auto mtc = std::make_unique<MassTransferCoefficientConstant>();
      mtc->constant_mtc_ = kLdf;
      ldf->mtc_model_ = std::move(mtc);
      adsorbent.adsorption_source_.mass_transfer_models_["Oxygen"] = std::move(ldf);
      adsorbent.adsorption_source_.heat_of_adsorption_["Oxygen"] = 0.0;

      CoupledSolver solver{data};
      std::string error;
      EXPECT_TRUE(data.PerformChecks(error)) << error;
      EXPECT_GE(solver.IntegrateCycle(data.cycle, false), 0)
        << "Solver did not complete (n_cells = " << n_cells << ")";
      const auto& log = solver.GetLog();
      const auto& dom = log.GetDomain();

      WaveRun run;
      const auto P = log.GetVariable(0, 0, VariableType::PRESSURE);
      const auto& Pf = P.values.back();
      run.dP_rel = std::abs(Pf[1] - Pf[Pf.size() - 2]) / kCasePressure;
      const auto vel = log.GetVariable(0, 0, VariableType::VELOCITY);
      const auto& uf = vel.values.back();
      double u_sum = 0.0;
      int u_cnt = 0;
      for (size_t i = 1; i + 1 < uf.size(); ++i) { u_sum += uf[i]; ++u_cnt; }
      run.u = (u_sum / u_cnt) / dom.et[1];
      auto [times, y] = OutletSeries(log, c.feed_from_top);
      run.t = times;
      for (double yi : y) run.phi.push_back(yi / kY0);
      return run;
    }

    // First time the (falling) outlet curve crosses the level phi, by linear
    // interpolation between samples.
    double CrossingTime(const WaveRun& run, double phi) {
      for (size_t i = 1; i < run.t.size(); ++i) {
        if (run.phi[i - 1] >= phi && run.phi[i] < phi) {
          const double s = (run.phi[i - 1] - phi) / (run.phi[i - 1] - run.phi[i]);
          return run.t[i - 1] + s * (run.t[i] - run.t[i - 1]);
        }
      }
      return -1.0;
    }

  } // namespace

  TEST(DesorptionSimpleWave, LangmuirDesorptionFollowsCharacteristics) {
    const std::vector<double> levels{0.6, 0.5, 0.4};
    std::vector<WaveRun> runs;
    for (int n : {100, 200}) runs.push_back(RunDesorption(n));

    for (const auto& run : runs) {
      ASSERT_FALSE(run.t.empty());
      ASSERT_LT(run.dP_rel, 2e-3) << "Pressure drop too large for a constant velocity";
      ASSERT_NEAR(run.phi.front(), 1.0, 1e-6) << "Bed did not start saturated";
      ASSERT_LT(run.phi.back(), 1e-3) << "Bed not fully desorbed";
    }
    ASSERT_NEAR(runs[0].u, runs[1].u, 1e-3 * runs[0].u);
    const double u = runs[1].u;
    const double L = 1.0;

    // GUARD: a wide, strongly nonlinear wave.
    const double t_head = CharacteristicTime(1.0, u, L), t_tail = CharacteristicTime(0.0, u, L);
    ASSERT_GT(t_tail - t_head, 5.0 * t_head) << "The wave does not spread";
    ASSERT_LT(t_tail, kRunTime);

    for (double phi : levels) {
      const double t_ref = CharacteristicTime(phi, u, L);
      const double t_c = CrossingTime(runs[0], phi), t_f = CrossingTime(runs[1], phi);
      ASSERT_GT(t_c, 0.0);
      ASSERT_GT(t_f, 0.0);
      const double e_c = t_c - t_ref, e_f = t_f - t_ref;
      const double order = std::log(std::abs(e_c / e_f)) / std::log(2.0);
      const double t_extrap = 2.0 * t_f - t_c;
      std::cout << "[DesorptionWave] c/c0=" << phi << " t(c)=" << t_ref << " s; n=100: " << t_c
                << " n=200: " << t_f << " extrapolated " << t_extrap << " (rel "
                << (t_extrap - t_ref) / t_ref << "), order " << order << "\n";
      // The error is the O(dx) numerical diffusion ...
      EXPECT_NEAR(order, 1.0, 0.25) << "Level " << phi << " converges at order " << order;
      EXPECT_LT(std::abs(e_f), 0.02 * t_ref) << "Level " << phi << " is off the characteristic";
      // ... and the grid-converged wave is the characteristic solution.
      EXPECT_NEAR(t_extrap, t_ref, 1e-3 * t_ref)
        << "Level c/c0 = " << phi << " reaches the outlet at " << t_extrap
        << " s (extrapolated), the characteristic gives " << t_ref << " s";
    }
  }

} // namespace adsim
