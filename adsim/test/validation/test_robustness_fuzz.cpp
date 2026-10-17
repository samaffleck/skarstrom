/*
 * Skarstrom - pressure and temperature swing adsorption simulator
 * Copyright (c) 2024-2026 Sam Affleck and contributors
 *
 * Licensed under the terms in the LICENSE file at the repository root.
 */

// Robustness fuzz harness.
//
// Generates many randomised-but-physical single-layer PSA decks, runs each to
// completion, and asserts every run lands in exactly one of two acceptable
// buckets:
//
//   CONVERGED   - solver returned success AND the result is finite with Sy = 1
//                 (mass conserved). Mole-balance closure is tracked as a
//                 reported statistic.
//   CLEAN_FAIL  - solver returned a specific negative IDA code, or threw a
//                 specific C++ exception. Either way it stopped deliberately.
//
// The one outcome the dossier must rule out is DIRTY: the solver reports
// success but the state is non-finite (a silently-propagated NaN). That is
// exactly what the denominator/step-retry guards prevent, so `dirty == 0` is
// the headline guarantee this test defends.
//
// "Never a hang" is underwritten by IDA's own work caps (max_num_steps /
// max_iters_per_step, ida_settings.h): a pathological deck returns
// IDA_TOO_MUCH_WORK rather than spinning forever. We tighten those caps per
// trial so every trial terminates quickly.
//
// Determinism: the RNG is seeded from a fixed constant plus the trial index, so
// the suite is fully reproducible (no wall-clock seeding). Set the environment
// variable ADSIM_FUZZ_TRIALS to scale the count up for an offline dossier run
// (thousands); the default is CI-friendly.

#define _USE_MATH_DEFINES
#include <cmath>
#include <cstdlib>
#include <exception>
#include <iostream>
#include <random>
#include <string>

#include "../src/pch.h"
#include "adsim/solver/coupled_solver.h"
#include "adsim/math/constants.h"
#include "adsim/equations/isotherm.h"
#include "adsim/equations/dispersion_model.h"
#include "adsim/equations/heat_transfer_model.h"
#include "adsim/boundary_conditions/boundary_conditions.h"
#include "adsim/boundary_conditions/boundary_stream.h"
#include "adsim/boundary_conditions/step.h"

#include <gtest/gtest.h>

#include "../fixtures/user_data_test.h"

namespace adsim {

  namespace {

    // Small deterministic RNG wrapper. mt19937_64 is reproducible across
    // platforms for a given seed, which the fixed seeding below relies on.
    struct Rng {
      std::mt19937_64 gen;
      explicit Rng(uint64_t seed) : gen(seed) {}
      double uni(double lo, double hi) {
        return std::uniform_real_distribution<double>(lo, hi)(gen);
      }
      // Log-uniform: for quantities that span orders of magnitude (pressure,
      // mass-transfer coefficient, flow) so each decade is sampled evenly.
      double logu(double lo, double hi) {
        return std::exp(uni(std::log(lo), std::log(hi)));
      }
      int uni_int(int lo, int hi) {
        return std::uniform_int_distribution<int>(lo, hi)(gen);
      }
      bool chance(double p) { return uni(0.0, 1.0) < p; }
    };

    int TrialCount() {
      if (const char* env = std::getenv("ADSIM_FUZZ_TRIALS")) {
        const int n = std::atoi(env);
        if (n > 0) return n;
      }
      return 96;  // CI-friendly default
    }

  } // namespace

  // Each trial builds a fresh deck on the fixture's `data` member, so the
  // fixture builder helpers (AddLayer, AddIsotherm, ...) and the conservation
  // checks (IsMassConserved, IsMolesConserved) can be reused directly.
  TEST_F(UserDataTest, RobustnessFuzz_PhysicalDecksConvergeOrFailCleanly) {
    const int trials = TrialCount();

    int converged = 0, clean_fail = 0, rejected = 0, dirty = 0, threw = 0;
    int mole_ok = 0, mole_checked = 0;

    for (int t = 0; t < trials; ++t) {
      // Fresh deck. Re-apply the wall material defaults SetUp() would set.
      data = UserData{};
      data.bed.wall_.heat_capacity_ = 537;
      data.bed.wall_.density_ = 7800;
      data.bed.wall_.thermal_conductivity_ = 16.2;

      // Fixed-constant + index seeding => reproducible, no wall-clock.
      Rng rng(0x5CA97507B0DULL + static_cast<uint64_t>(t) * 0x9E3779B97F4A7C15ULL);

      // ~40% of decks are drawn from an adversarial tier: sub-millimetre
      // particles, near-zero wall thickness, sharp isotherms, very high flow,
      // coarse or very fine grids, tiny step ceilings. These are still
      // *physical* (a real, if extreme, column) but deliberately push the
      // solver toward the failure boundary so the "never NaN, never hang"
      // guarantee is tested under genuine stress, not just easy cases.
      const bool hard = rng.chance(0.40);

      AddComponents({ "Oxygen", "Nitrogen" });
      auto [layer_id, ads_id] = AddLayer("13X");
      auto& layer = GetLayer(layer_id);

      // Geometry
      layer.domain_.length = rng.uni(0.2, 2.0);            // [m]
      layer.eb_ = rng.uni(0.34, 0.50);                     // interstitial void
      layer.particle_diameter_ = hard ? rng.logu(1e-4, 0.5e-3)
                                       : rng.logu(0.5e-3, 5e-3);   // [m]
      layer.domain_.n_cells = hard ? rng.uni_int(8, 80) : rng.uni_int(15, 40);
      std::get<CylinderBed>(data.bed.shape_).diameter = rng.uni(0.02, 0.30); // [m]

      auto& ads = GetAdsorbent(ads_id);
      ads.solid_density_ = rng.uni(700.0, 1300.0);         // [kg/m3]
      ads.ep_ = rng.uni(0.30, 0.50);                       // macropore void

      // Isotherm on the adsorbing species (Oxygen); Nitrogen stays inert.
      auto iso = std::make_unique<LangmuirIsotherm>();
      iso->params_["Oxygen"].m1 = hard ? rng.uni(5.0, 12.0) : rng.uni(0.5, 5.0);
      iso->params_["Oxygen"].b1 = hard ? rng.logu(1e-9, 1e-8)
                                       : rng.logu(1e-11, 1e-9);   // affinity pre-exp
      iso->params_["Oxygen"].q1 = rng.uni(1.0e4, 4.0e4);   // heat term [J/mol]
      AddIsotherm({ "Oxygen" }, ads_id, std::move(iso));

      const double mtc_lo = hard ? 5.0 : 0.05, mtc_hi = hard ? 50.0 : 5.0;
      AddLinearDrivingForceModel(ads_id, "Oxygen", rng.logu(mtc_lo, mtc_hi));   // [1/s]
      AddLinearDrivingForceModel(ads_id, "Nitrogen", rng.logu(mtc_lo, mtc_hi));

      // Heat transfer: adiabatic half the time, otherwise a finite coefficient.
      AddConstantInternalHeatTransferCoefficientModel(rng.chance(0.5) ? 0.0 : rng.uni(1.0, 50.0), layer_id);
      AddConstantExternalHeatTransferCoefficientModel(rng.chance(0.5) ? 0.0 : rng.uni(1.0, 50.0), layer_id);

      auto disp = std::make_unique<DispersionModelConstant>();
      disp->dispersion_coefficient_ = rng.chance(0.5) ? 0.0 : rng.logu(1e-6, 1e-3);
      layer.equations_.dispersion_model_ = std::move(disp);
      auto hdisp = std::make_unique<HeatDispersionModelConstant>();
      hdisp->dispersion_coefficient_ = 0.0;
      layer.equations_.heat_dispersion_model_ = std::move(hdisp);

      // Wall thickness: nominal down to 1 mm; adversarial down to 20 microns
      // to drive the thin-wall (d_o^2 - d_i^2) denominator floor hard.
      data.bed.wall_.thickness_ = hard ? rng.logu(2e-5, 1e-3)
                                        : rng.logu(1e-3, 2e-2);    // [m]
      data.bed.ambient_temperature_ = rng.uni(283.0, 330.0);

      const double P0 = hard ? rng.logu(5e5, 2e6) : rng.logu(5e4, 5e5);  // [Pa]
      const double T0 = rng.uni(283.0, 330.0);             // [K]
      const double yO2_0 = rng.uni(0.05, 0.95);
      AddInitialConditions(P0, T0, T0,
        { {"Oxygen", yO2_0}, {"Nitrogen", 1.0 - yO2_0} }, 0.0);

      AddStep("Feed", hard ? rng.uni(300.0, 600.0) : rng.uni(50.0, 300.0));
      Step& step = data.cycle.GetStep("Feed");

      const double flow = hard ? rng.logu(2.0, 20.0) : rng.logu(0.02, 2.0);  // [mol/s]
      const double yO2_feed = rng.uni(0.05, 0.95);
      const double Tfeed = rng.uni(283.0, 330.0);
      auto inlet = CreateNoValveInletStream(CreateMolarFlowBC(flow),
        { {"Oxygen", yO2_feed}, {"Nitrogen", 1.0 - yO2_feed} }, Tfeed);
      auto feed_bc = std::make_unique<FlowBoundaryCondition>();
      feed_bc->direction_ = BoundaryDirection::INLET;
      feed_bc->location_ = BoundaryLocation::BOTTOM;
      feed_bc->streams_.push_back(std::move(inlet));
      step.bottom_boundary_condition_ = std::move(feed_bc);

      auto outlet = std::make_unique<PressureBoundaryCondition>();
      outlet->direction_ = BoundaryDirection::OUTLET;
      outlet->location_ = BoundaryLocation::TOP;
      outlet->pressure_ = CreateConstantBoundaryValue(P0);
      outlet->mole_fraction_["Oxygen"] = yO2_0;
      outlet->mole_fraction_["Nitrogen"] = 1.0 - yO2_0;
      step.top_boundary_condition_ = std::move(outlet);

      ConfigureCycle(1);
      data.settings.max_step = hard ? rng.uni(0.05, 0.5) : rng.uni(0.5, 3.0);
      // Tighten the work caps so a non-convergent deck returns TOO_MUCH_WORK
      // quickly instead of grinding to the default 100k-step ceiling.
      data.settings.max_num_steps = 20000;
      data.settings.max_iters_per_step = 20000;

      // Reject decks the input validator rejects: that is an invalid *input*,
      // not a solver robustness failure, so it does not count against us.
      std::string err;
      if (!data.PerformChecks(err)) { ++rejected; continue; }

      // A thrown exception is an acceptable *specific* failure as long as it
      // does not crash the process, so treat it as a clean failure.
      int retval = 0;
      bool caught = false;
      try {
        CoupledSolver solver{data};
        retval = solver.IntegrateCycle(data.cycle, false);

        const auto& log = solver.GetLog();
        const bool finite_and_sum1 = IsMassConserved(log);

        if (retval >= 0) {
          if (!finite_and_sum1) {
            ++dirty;
            ADD_FAILURE() << "Trial " << t << " (seed-index " << t
                          << "): solver reported success (retval=" << retval
                          << ") but the state is non-finite / Sy != 1 "
                             "(silent NaN got through the guards).";
          } else {
            ++converged;
            ++mole_checked;
            if (IsMolesConserved(log)) ++mole_ok;
          }
        } else {
          ++clean_fail;
        }
      } catch (const std::exception& e) {
        caught = true;
        (void)e;
      } catch (...) {
        caught = true;
      }
      if (caught) { ++clean_fail; ++threw; }
    }

    // Dossier summary line.
    std::cout << "\n[Robustness fuzz] trials=" << trials
              << "  converged=" << converged
              << "  clean_fail=" << clean_fail << " (of which threw=" << threw << ")"
              << "  rejected_input=" << rejected
              << "  DIRTY=" << dirty << "\n"
              << "  mole-balance closure on converged runs: " << mole_ok
              << "/" << mole_checked
              << (mole_checked ? "  (" + std::to_string(100.0 * mole_ok / mole_checked) + "%)" : "")
              << "\n" << std::flush;

    // --- Guarantees the dossier rests on ---

    // 1. No run may report success while returning NaN/garbage. This is the
    //    headline property the solver guards exist to enforce.
    EXPECT_EQ(dirty, 0)
      << "At least one deck converged to a non-finite state (silent NaN).";

    // 2. The harness must actually exercise real solves, not reject or fail
    //    everything (which would make (1) vacuously true and hide a regression
    //    that broke the solver outright). With ~40% of decks in the adversarial
    //    tier a meaningful minority may fail cleanly, so this bound is
    //    deliberately loose - it only guards against a total collapse.
    EXPECT_GT(converged, trials / 10)
      << "Too few decks converged (" << converged << "/" << trials
      << "); solver or harness regressed.";

    // 3. Mass (finite + Sy=1) holds on every converged run by construction of
    //    (1). Global mole balance should additionally close on the large
    //    majority of converged runs; a big drop signals a flux-accounting
    //    regression.
    if (mole_checked > 0) {
      EXPECT_GE(static_cast<double>(mole_ok) / mole_checked, 0.80)
        << "Mole-balance closure rate on converged runs regressed: "
        << mole_ok << "/" << mole_checked;
    }
  }

} // namespace adsim
