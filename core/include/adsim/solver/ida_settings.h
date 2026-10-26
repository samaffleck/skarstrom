/*
 * Skarstrom - pressure and temperature swing adsorption simulator
 * Copyright (c) 2024-2026 Sam Affleck and contributors
 *
 * Licensed under the terms in the LICENSE file at the repository root.
 */
#pragma once

// STL includes
#include <vector>
#include <cstdio>

// SUNDIALS includes
#include <ida/ida.h>                  // IDA core
#include <sundials/sundials_types.h>  // realtype
#include <nvector/nvector_serial.h>   // or your N_Vector of choice

// Adsim includes
#include "adsim/solver/cyclic_steady_state.h"
#include "adsim/solver/warning_settings.h"
#include "adsim/solver/jacobian.h"
#include "adsim/math/expression_param.h"


namespace adsim {

  // Convective face-reconstruction scheme for the species balance. Upwind is
  // first-order (the historical default); WENO and van Leer are higher-order
  // reconstructions of the interior face value. Selecting a scheme changes only
  // how the face mole fraction is interpolated, not the conserved cell-centred
  // state, so total species mass is unaffected.
  enum class AdvectionScheme { Upwind, WENO, VanLeer };

  class IDASettings {
  public:
    long int max_num_steps = 100000;
    // Cap on OneStep iterations within a single IntegrateStep call. Guards
    // against the outer loop spinning forever when IDA shrinks the step size
    // aggressively (e.g. near a singular Jacobian). 0 disables the cap.
    long int max_iters_per_step = 50000;
    ExpressionParam init_step{1e-10};
    ExpressionParam min_step{1e-15};
    // Maximum internal step [s]; 0 (the default) means no cap, leaving the step
    // size to IDA's error control. The logged trajectory no longer depends on
    // it (long steps are filled in from IDA's interpolant).
    ExpressionParam max_step{0.0};
    ExpressionParam reltol{1e-6};
    int max_order = 5; // BDF [1->5]
    int max_nonlin_iters = 15; // per step
    int max_conv_fails = 200;
    int max_err_test_fails = 40;
    ExpressionParam nonlin_conv_coef{0.1}; // 0 => default; otherwise (0,1]

    // Convective face-reconstruction scheme for the species balance.
    AdvectionScheme advection_scheme_ = AdvectionScheme::Upwind;

    int max_cycles_to_keep_ = 10;

    // Cap on the number of full spatial-profile snapshots recorded per step in
    // the EquationLog (the largest output by far — every internal solver step
    // was previously logged, producing ~1000 points/step and GB-scale CSVs).
    // Snapshots are spaced uniformly in time (step_time / max_log_points_per_step_);
    // the first and last point of every step are always kept, so the cyclic
    // steady-state check and acceleration (which read only those two points) are
    // unaffected and results are unchanged. The per-step outlet/boundary log that
    // feeds KPIs is never decimated. <= 0 disables the cap (log every step).
    int max_log_points_per_step_ = 200;

    // Jacobian construction (transient, not serialised): benchmarking and
    // verification switches only.
    JacobianMode jacobian_mode_ = JacobianMode::Coloured;
    bool record_jacobian_pattern_ = false;

    // Reduced DAE unknowns (transient, not serialised; off by default, see
    // solver_performance.md item 9 in samaffleck/skarstrom-gui). With local
    // thermal equilibrium the solid temperature is the algebraic constraint
    // Ts = Tf: aliasing it to Tf removes one unknown per cell. The gas
    // concentration C = P / (R T) has no time derivative in any residual:
    // deriving it from P and T removes another. Both shrink the band (2 neq
    // each side), the LU and the colour count. The equations solved are the
    // same, but the error test loses the removed entries, so the integration
    // path (step sizes, results in the last digits) changes.
    bool alias_lte_solid_temperature_ = false;
    bool derived_concentration_ = false;
    // Size the banded matrix from the coupling reach measured by a difference-
    // quotient probe at the initial state (BandedJacobian::ProbeBandwidth)
    // instead of the 2 neq each side that covers every scheme (about 45% less
    // LU work; results are byte-identical). The coloured Jacobian's band
    // audits keep scanning the full 2 neq reach; if one finds a coupling
    // outside the matrix band, the solver rebuilds the linear solver with the
    // full band and restarts from the current point (SolverData::WidenBand).
    bool band_from_probe_ = true;
    // Testing only: narrows the probed band by this many diagonals on each side
    // so that an audit finds couplings outside it and WidenBand runs.
    int band_probe_shrink_for_testing_ = 0;

    // Solver-settings experiments (transient, not serialised; see item 10 of
    // solver_performance.md in samaffleck/skarstrom-gui).
    // IDASetSuppressAlg: leave the algebraic variables (P, u, and C / Ts when
    // they are unknowns) out of the local error test.
    bool suppress_alg_ = false;
    // IDASetDeltaCjLSetup: IDA refreshes the Jacobian when cj / cj_old leaves
    // [(1 - dcj) / (1 + dcj), its inverse]; SUNDIALS' default is 0.25. < 0
    // leaves the default.
    double delta_cj_lsetup_ = -1.0;
    // Tag the pressure as a differential variable (its residual carries dP/dt)
    // instead of algebraic: with suppress_alg_ it then stays in the error test
    // while the ghost cells and the face velocity leave it.
    bool pressure_differential_ = false;

    CyclicSteadyState css{};
    AccelerateCyclicSteadyState acc_css{};
    WarningSettings warnings_{};  // advisory warnings a run reports

    void EvaluateExpressions(const EvalContext& ctx) {
      init_step.Evaluate(ctx);
      min_step.Evaluate(ctx);
      max_step.Evaluate(ctx);
      reltol.Evaluate(ctx);
      nonlin_conv_coef.Evaluate(ctx);
      css.EvaluateExpressions(ctx);
    }

  private:
    friend class ::cereal::access;
    template<class Archive>
    void serialize(Archive& ar, std::uint32_t const version) {
      ar(
        cereal::make_nvp("Maximum_Number_Of_Steps", max_num_steps),
        cereal::make_nvp("Initial_Step_Size", init_step),
        cereal::make_nvp("Minimum_Step_Size", min_step),
        cereal::make_nvp("Maximum_Step_Size", max_step),
        cereal::make_nvp("Relative_Tolerance", reltol),
        cereal::make_nvp("Maximum_Backwards_Difference_Formula_Order", max_order),
        cereal::make_nvp("Maximum_Nonlinear_Itterations_Per_Step", max_nonlin_iters),
        cereal::make_nvp("Maximum_Nonlinear_Convergence_Failures", max_conv_fails),
        cereal::make_nvp("Maximum_Error_Test_Failures", max_err_test_fails),
        cereal::make_nvp("Cyclic_Steady_State", css),
        cereal::make_nvp("max_cycles_to_keep", max_cycles_to_keep_)
      );
      if (version >= 1) {
        ar(cereal::make_nvp("Maximum_Iterations_Per_Step", max_iters_per_step));
      }
      if (version >= 2) {
        ar(cereal::make_nvp("Advection_Scheme", advection_scheme_));
      }
      if (version >= 3) {
        ar(cereal::make_nvp("Accelerate_Cyclic_Steady_State", acc_css));
      }
      if (version >= 4) {
        ar(cereal::make_nvp("Maximum_Log_Points_Per_Step", max_log_points_per_step_));
      }
      if (version >= 6) {
        ar(cereal::make_nvp("Warnings", warnings_));
      }
      // v5: the default max step changed from 2 s to no cap. Decks saved before
      // v5 that still carry the old default ("2") move to the new one; any value
      // or expression the user chose is kept.
      if constexpr (Archive::is_loading::value) {
        if (version < 5 && max_step.GetExpr() == "2") max_step = 0.0;
      }
    }
  };

  inline int ida_accumulate(int current, int flag, const char* what) {
    if (flag != IDA_SUCCESS) {
      std::fprintf(stderr, "[IDA] %s failed with flag = %d\n", what, flag);
      if (current == IDA_SUCCESS) return flag;
    }
    return current;
  }

  inline int ApplyIDASettings(void* ida_mem, const IDASettings& s) {
    int rc = IDA_SUCCESS;

    if (s.max_num_steps > 0)
      rc = ida_accumulate(rc, IDASetMaxNumSteps(ida_mem, s.max_num_steps), "IDASetMaxNumSteps");
    if (s.init_step.Get() > 0.0)
      rc = ida_accumulate(rc, IDASetInitStep(ida_mem, s.init_step.Get()), "IDASetInitStep");
    if (s.min_step.Get() > 0.0)
      rc = ida_accumulate(rc, IDASetMinStep(ida_mem, s.min_step.Get()), "IDASetMinStep");
    if (s.max_step.Get() > 0.0)
      rc = ida_accumulate(rc, IDASetMaxStep(ida_mem, s.max_step.Get()), "IDASetMaxStep");

    if (s.max_order > 0)
      rc = ida_accumulate(rc, IDASetMaxOrd(ida_mem, s.max_order), "IDASetMaxOrd");
    if (s.max_nonlin_iters > 0)
      rc = ida_accumulate(rc, IDASetMaxNonlinIters(ida_mem, s.max_nonlin_iters), "IDASetMaxNonlinIters");
    if (s.max_conv_fails > 0)
      rc = ida_accumulate(rc, IDASetMaxConvFails(ida_mem, s.max_conv_fails), "IDASetMaxConvFails");
    if (s.max_err_test_fails > 0)
      rc = ida_accumulate(rc, IDASetMaxErrTestFails(ida_mem, s.max_err_test_fails), "IDASetMaxErrTestFails");
    if (s.nonlin_conv_coef.Get() > 0.0)
      rc = ida_accumulate(rc, IDASetNonlinConvCoef(ida_mem, s.nonlin_conv_coef.Get()), "IDASetNonlinConvCoef");
    if (s.suppress_alg_)
      rc = ida_accumulate(rc, IDASetSuppressAlg(ida_mem, SUNTRUE), "IDASetSuppressAlg");
    if (s.delta_cj_lsetup_ >= 0.0 && s.delta_cj_lsetup_ < 1.0)
      rc = ida_accumulate(rc, IDASetDeltaCjLSetup(ida_mem, s.delta_cj_lsetup_), "IDASetDeltaCjLSetup");

    // The step-size ratio bounds (IDASetEtaMax/EtaMin/EtaLow/...) are left at
    // the SUNDIALS defaults. Note SUNDIALS does not reject out-of-range eta
    // values: it silently substitutes its default and still returns
    // IDA_SUCCESS (ida_io.c), so a return-code check cannot catch a bad one.
    // A former IDASetEtaMax(0.9) (must be > 1) was replaced by 2 this way, and
    // the accompanying IDASetEtaMin(0.1) never binds: after a successful step
    // IDA's proposed ratio is always >= ~0.58, above the default floor of 0.5.

    // First failing flag (IDA_SUCCESS if every call succeeded). Callers must
    // check it: a rejected setting otherwise leaves IDA on its default.
    return rc;
  }

} // end adsim namespace

CEREAL_CLASS_VERSION(adsim::IDASettings, 6);
