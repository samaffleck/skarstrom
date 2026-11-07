/*
 * Skarstrom - pressure and temperature swing adsorption simulator
 * Copyright (c) 2024-2026 Sam Affleck and contributors
 *
 * Licensed under the terms in the LICENSE file at the repository root.
 */
#include "../pch.h"
#include "adsim/solver/cycle_sensitivity.h"
#include "adsim/objects/user_data.h"
#include "adsim/solver/coupled_solver.h"
#include "adsim/solver/cycle.h"

// The engine links IDAS (a superset of IDA with the same entry points), so
// the sensitivity functions act on the integrator SolverData created.
#include <idas/idas.h>
#include <nvector/nvector_serial.h>

#include <chrono>
#include <cmath>
#include <limits>
#include <memory>

namespace adsim {

  namespace {

    // Carries the sensitivities through every IDAReInit of the run: reads
    // them (and their derivatives) at the restart time from IDA's history
    // before IDAReInit discards it, zeroes the algebraic derivatives as
    // RestartFromConsistentIC does for the state, and restarts them with
    // IDASensReInit, so the IDACalcIC that follows makes the algebraic
    // sensitivities consistent with the new step's boundary conditions.
    class SensitivityCarrier : public SolverData::ReInitObserver {
    public:
      SensitivityCarrier(N_Vector* yS, N_Vector* ypS, int ns, int ism, const sunrealtype* id)
          : yS_(yS), ypS_(ypS), ns_(ns), ism_(ism), id_(id) {}

      // The next re-initialisation starts the sensitivities from yS / ypS as
      // they are (the initial directions), not from IDA's history.
      void UseCurrentValuesNext() { pending_ = true; }

      void BeforeReInit(void* ida_mem, double t) override {
        Accumulate(ida_mem);
        // A re-initialisation with no step since the last one is the
        // solver's fallback after a failed consistent-IC solve: the
        // sensitivities were not made consistent and cannot be carried on.
        long nst = 0;
        if (!pending_ && IDAGetNumSteps(ida_mem, &nst) == IDA_SUCCESS && nst == 0 && reinits_ > 0)
          failed_ = IDA_ILL_INPUT;
        if (!pending_ && failed_ == 0) {
          int f = IDAGetSensDky(ida_mem, t, 0, yS_);
          if (f == IDA_SUCCESS) f = IDAGetSensDky(ida_mem, t, 1, ypS_);
          if (f != IDA_SUCCESS) failed_ = f;
        }
        for (int k = 0; k < ns_; ++k) {
          sunrealtype* yp = N_VGetArrayPointer(ypS_[k]);
          const sunindextype n = N_VGetLength(ypS_[k]);
          for (sunindextype i = 0; i < n; ++i)
            if (id_[i] == 0.0) yp[i] = 0.0;
        }
      }

      int AfterReInit(void* ida_mem, double) override {
        pending_ = false;
        if (failed_ != 0) return failed_;
        ++reinits_;
        return IDASensReInit(ida_mem, ism_, yS_, ypS_);
      }

      // Adds IDAS's sensitivity counters (reset by IDASensReInit) to the totals.
      void Accumulate(void* ida_mem) {
        long v = 0;
        if (IDAGetSensNumResEvals(ida_mem, &v) == IDA_SUCCESS) res_evals_ += v;
        if (IDAGetSensNumNonlinSolvIters(ida_mem, &v) == IDA_SUCCESS) nonlin_iters_ += v;
        if (IDAGetSensNumNonlinSolvConvFails(ida_mem, &v) == IDA_SUCCESS) conv_fails_ += v;
      }

      bool Failed() const { return failed_ != 0; }

      long res_evals_ = 0;
      long nonlin_iters_ = 0;
      long conv_fails_ = 0;
      int reinits_ = 0;

    private:
      N_Vector* yS_;
      N_Vector* ypS_;
      int ns_;
      int ism_;
      const sunrealtype* id_;
      bool pending_ = false;
      int failed_ = 0;
    };

    // The sensitivity residual F_y s + F_y' s' by a difference quotient of
    // the engine's residual along (s, s'). The increment is scaled to the
    // state, eps^(1/3) (centred) or eps^(1/2) (forward) relative to
    // |y| + abstol in RMS. IDAS's internal quotient instead perturbs by
    // about reltol * |y|, whose rounding error is of order one in IDA's
    // weighted norm, and IDACalcIC's line search then fails at the first
    // consistent-IC solve.
    // The quotient uses its own scratch vectors, not the three IDAS passes:
    // IDACalcIC with the simultaneous corrector passes its trial state
    // (ynew, which is IDAS's tempv2) as y and tempv2 again as the second
    // scratch vector, so writing that scratch would overwrite y.
    struct SensResidualContext {
      const sunrealtype* abstol = nullptr;
      bool centred = true;
      N_Vector ytmp = nullptr, yptmp = nullptr, rtmp = nullptr;
    };
    thread_local const SensResidualContext* g_sens_context = nullptr;

    int SensResidualDQ(int ns, sunrealtype t, N_Vector y, N_Vector yp, N_Vector res, N_Vector* yS,
                       N_Vector* ypS, N_Vector* resS, void* user_data, N_Vector, N_Vector, N_Vector) {
      const SensResidualContext& ctx = *g_sens_context;
      N_Vector tmp1 = ctx.ytmp, tmp2 = ctx.yptmp, tmp3 = ctx.rtmp;
      const sunindextype n = N_VGetLength(y);
      const sunrealtype* yd = N_VGetArrayPointer(y);
      const double eps = std::numeric_limits<double>::epsilon();
      const double rel = ctx.centred ? std::cbrt(eps) : std::sqrt(eps);
      for (int k = 0; k < ns; ++k) {
        const sunrealtype* s = N_VGetArrayPointer(yS[k]);
        double sum = 0.0;
        for (sunindextype i = 0; i < n; ++i) {
          const double q = s[i] / (std::abs(yd[i]) + ctx.abstol[i]);
          sum += q * q;
        }
        const double norm = std::sqrt(sum / static_cast<double>(n));
        if (norm == 0.0) {
          // s = 0: only F_y' s' remains; perturb y' alone, by a unit-scaled step.
          N_VScale(1.0, y, tmp1);
          N_VLinearSum(1.0, yp, rel, ypS[k], tmp2);
          int f = UpdateResiduals(t, tmp1, tmp2, resS[k], user_data);
          if (f != 0) return f;
          N_VLinearSum(1.0 / rel, resS[k], -1.0 / rel, res, resS[k]);
          continue;
        }
        const double sigma = rel / norm;
        N_VLinearSum(1.0, y, sigma, yS[k], tmp1);
        N_VLinearSum(1.0, yp, sigma, ypS[k], tmp2);
        int f = UpdateResiduals(t, tmp1, tmp2, resS[k], user_data);
        if (f != 0) return f;
        if (ctx.centred) {
          N_VLinearSum(1.0, y, -sigma, yS[k], tmp1);
          N_VLinearSum(1.0, yp, -sigma, ypS[k], tmp2);
          f = UpdateResiduals(t, tmp1, tmp2, tmp3, user_data);
          if (f != 0) return f;
          N_VLinearSum(0.5 / sigma, resS[k], -0.5 / sigma, tmp3, resS[k]);
        } else {
          N_VLinearSum(1.0 / sigma, resS[k], -1.0 / sigma, res, resS[k]);
        }
      }
      return 0;
    }

    struct VectorArray {
      std::vector<N_Vector> v;
      ~VectorArray() { for (N_Vector x : v) if (x) N_VDestroy(x); }
      N_Vector* data() { return v.empty() ? nullptr : v.data(); }
    };


  } // namespace

  CycleSensitivityResult CycleJacobianVector(UserData& data,
                                             const std::vector<std::vector<double>>& directions,
                                             const CycleSensitivityOptions& options) {
    CycleSensitivityResult out;
    auto fail = [&](const std::string& why) {
      out.ok = false;
      out.message = why;
      return out;
    };

    auto& cycle = data.cycle;
    if (cycle.step_list_.empty()) return fail("the cycle has no steps");
    if (!GetConnectedEndpoints(cycle.step_list_).empty())
      return fail("connected streams are not supported by the prototype");

    const double saved_reltol = data.settings.reltol;
    if (options.reltol > 0.0) data.settings.reltol = options.reltol;
    struct RestoreReltol {
      UserData& d; double r;
      ~RestoreReltol() { d.settings.reltol = r; }
    } restore{ data, saved_reltol };

    data.TakeResidualError();
    std::string error;
    if (!data.PerformChecks(error)) return fail("checks failed: " + error);
    for (auto& step : cycle.step_list_) step.ClearCycleSummary();

    // As CoupledSolver::Initialise.
    data.Initialise();
    data.profile_.Reset();
    auto& state = data.state;
    Step& first = cycle.step_list_[0];
    state.SetCurrentStep(&first);
    data.step_start_time_ = 0.0;
    first.Initialise(state);

    SolverData solver;
    if (solver.Initialise(data, error) != IDA_SUCCESS) return fail("initialisation failed: " + error);
    for (auto& step : cycle.step_list_) step.ClearOutletConditions();
    cycle.ResolveBufferTankSources();

    void* ida = solver.IdaMemory();
    N_Vector y = solver.StateVector();
    N_Vector yp = solver.DerivativeVector();
    const sunrealtype* id = N_VGetArrayPointer(solver.IdVector());
    const sunrealtype* abstol = N_VGetArrayPointer(solver.AbsTolVector());
    const sunindextype n = N_VGetLength(y);
    for (sunindextype i = 0; i < n; ++i) {
      if (id[i] != 1.0) continue;
      out.differential_index.push_back(static_cast<size_t>(i));
      out.description.push_back(state.DescribeBandedIndex(static_cast<size_t>(i)));
      out.abstol.push_back(abstol[i]);
    }
    const size_t nd = out.differential_index.size();
    if (!options.shift.empty() && options.shift.size() != nd)
      return fail("shift has " + std::to_string(options.shift.size()) + " entries, expected "
                  + std::to_string(nd));
    const int ns = options.sensitivities ? static_cast<int>(directions.size()) : 0;
    for (int k = 0; k < ns; ++k)
      if (directions[k].size() != nd)
        return fail("direction " + std::to_string(k) + " has " + std::to_string(directions[k].size())
                    + " entries, expected " + std::to_string(nd));

    const long max_iters = data.settings.max_iters_per_step;
    // One cycle from the integrator's current point, as CoupledSolver::
    // IntegrateStep runs its steps (without logging or retries). Returns an
    // error message, empty on success.
    auto run_cycle = [&]() -> std::string {
      int step_index = 0;
      for (auto& step : cycle.step_list_) {
        step.Initialise(state);
        state.domain_vectors.PopulateEnergySources(step.energy_sources_, data.bed);
        const double t_start = solver.GetTime();
        data.step_start_time_ = t_start;
        const double t_final = t_start + step.step_time_;
        int flag = solver.SetStopTime(t_final);
        if (flag != IDA_SUCCESS) return "IDASetStopTime failed";
        solver.OnStepStart(0);
        solver.RestartForNewStep();
        double t_curr = t_start;
        long iters = 0;
        while (true) {
          if (solver.BandNeedsWidening() && !solver.WidenBand(data))
            return "band widening failed in step " + std::to_string(step_index);
          flag = solver.OneStep(t_final, t_curr);
          if (flag < 0) {
            std::string why = "step \"" + step.step_name_ + "\" failed at t=" + std::to_string(t_curr)
                              + " with IDA code " + std::to_string(flag);
            const std::string w = data.TakeResidualError();
            if (!w.empty()) why += " (" + w + ")";
            return why;
          }
          if (flag == IDA_TSTOP_RETURN || t_curr >= t_final) break;
          if (max_iters > 0 && ++iters >= max_iters) return "too many internal steps";
        }
        ++step_index;
      }
      return {};
    };

    // Plain cycles before the one measured: its start is then a cycle
    // boundary of a running cycle rather than the deck's bed at rest.
    for (int c = 0; c < options.warmup_cycles; ++c) {
      const std::string why = run_cycle();
      if (!why.empty()) return fail("warm-up cycle " + std::to_string(c + 1) + ": " + why);
    }
    const auto measure_start = std::chrono::steady_clock::now();
    const long residual_calls_start = data.profile_.residual_calls;
    const IntegrationStats stats_start = solver.GetStatistics();

    // The cycle's initial state: the current point plus the shift, with the
    // first step's boundary conditions set up again for it (as a warm start
    // does), made consistent again below.
    sunrealtype* yd = N_VGetArrayPointer(y);
    data.step_start_time_ = solver.GetTime();
    for (size_t k = 0; k < options.shift.size(); ++k) yd[out.differential_index[k]] += options.shift[k];
    state.UnpackBanded(yd, N_VGetArrayPointer(yp));
    first.Initialise(state);

    // Sensitivities: s(0) = v on the differential components, 0 elsewhere;
    // s'(0) = 0 (IDACalcIC computes the differential derivatives and the
    // algebraic values).
    VectorArray yS, ypS, atolS, scratch;
    std::vector<sunrealtype> pbar(static_cast<size_t>(ns), 1.0), p(static_cast<size_t>(ns), 0.0);
    const int ism = options.corrector == CycleSensitivityOptions::Corrector::Staggered
                        ? IDA_STAGGERED : IDA_SIMULTANEOUS;
    std::unique_ptr<SensitivityCarrier> carrier;
    SensResidualContext sens_context;
    struct ClearContext { ~ClearContext() { g_sens_context = nullptr; } } clear_context;
    if (ns > 0) {
      for (int k = 0; k < ns; ++k) {
        yS.v.push_back(N_VClone(y));
        ypS.v.push_back(N_VClone(y));
        atolS.v.push_back(N_VClone(y));
        N_VConst(0.0, yS.v.back());
        N_VConst(0.0, ypS.v.back());
        sunrealtype* s = N_VGetArrayPointer(yS.v.back());
        for (size_t j = 0; j < nd; ++j) s[out.differential_index[j]] = directions[k][j];
        // s is y's derivative along v, so it is measured on y's scale times
        // the size of v: abstol_i + reltol |y_i|, times v's RMS norm relative
        // to |y| + abstol. With the state's own abstol alone, a component
        // whose sensitivity starts at zero (the algebraic pressures, abstol
        // 1e-6 Pa) gets a weight the IC Newton cannot meet.
        const sunrealtype* yd0 = N_VGetArrayPointer(y);
        double sum = 0.0;
        for (size_t j = 0; j < nd; ++j) {
          const size_t i = out.differential_index[j];
          const double q = directions[k][j] / (std::abs(yd0[i]) + abstol[i]);
          sum += q * q;
        }
        const double vnorm = nd > 0 ? std::sqrt(sum / static_cast<double>(nd)) : 0.0;
        sunrealtype* a = N_VGetArrayPointer(atolS.v.back());
        for (sunindextype i = 0; i < n; ++i)
          a[i] = (abstol[i] + data.settings.reltol * std::abs(yd0[i])) * (vnorm > 0.0 ? vnorm : 1.0);
      }
      sens_context.abstol = abstol;
      for (N_Vector* v : { &sens_context.ytmp, &sens_context.yptmp, &sens_context.rtmp }) {
        scratch.v.push_back(N_VClone(y));
        *v = scratch.v.back();
      }
      sens_context.centred = options.dq == CycleSensitivityOptions::DifferenceQuotient::Centred;
      g_sens_context = &sens_context;
      int f = IDASensInit(ida, ns, ism, options.idas_internal_dq ? nullptr : SensResidualDQ, yS.data(),
                          ypS.data());
      // The sensitivity residual is linear in s, so tolerances only steer the
      // corrector's convergence test; error control stays off so the step
      // sequence is the base run's.
      const double reltol = data.settings.reltol;
      const double rtolS = options.sens_reltol > 0.0 ? options.sens_reltol : reltol;
      for (int k = 0; k < ns; ++k) N_VScale(rtolS / reltol, atolS.v[k], atolS.v[k]);
      if (f == IDA_SUCCESS) f = IDASensSVtolerances(ida, rtolS, atolS.data());
      if (f == IDA_SUCCESS) f = IDASetSensErrCon(ida, SUNFALSE);
      // IDAS's internal DQ perturbs a parameter alongside y and y'; these
      // dummy parameters do not enter the residual, so it computes F_y s + F_y' s'.
      if (f == IDA_SUCCESS && options.idas_internal_dq) {
        f = IDASetSensParams(ida, p.data(), pbar.data(), nullptr);
        if (f == IDA_SUCCESS)
          f = IDASetSensDQMethod(ida, options.dq == CycleSensitivityOptions::DifferenceQuotient::Centred
                                          ? IDA_CENTERED : IDA_FORWARD, 0.0);
      }
      if (f != IDA_SUCCESS) return fail("IDAS sensitivity set-up failed with code " + std::to_string(f));
      carrier = std::make_unique<SensitivityCarrier>(yS.data(), ypS.data(), ns, ism, id);
      carrier->UseCurrentValuesNext();
      solver.SetReInitObserver(carrier.get());
    }
    if (!solver.RestartFromConsistentIC())
      return fail("consistent initial condition at the shifted state failed" + [&] {
        const std::string w = data.TakeResidualError();
        return w.empty() ? std::string{} : ": " + w;
      }());
    {
      sunrealtype* yd2 = N_VGetArrayPointer(y);
      for (size_t k = 0; k < nd; ++k) out.x.push_back(yd2[out.differential_index[k]]);
    }

    {
      const std::string why = run_cycle();
      if (!why.empty()) return fail(why);
      if (carrier && carrier->Failed())
        return fail("the sensitivities could not be restarted at a step switch "
                    "(the consistent initial condition with sensitivities failed)");
    }

    // End of the cycle.
    {
      const sunrealtype* yd2 = N_VGetArrayPointer(y);
      for (size_t k = 0; k < nd; ++k) out.g.push_back(yd2[out.differential_index[k]]);
    }
    if (ns > 0) {
      sunrealtype tret = 0.0;
      if (IDAGetSens(ida, &tret, yS.data()) != IDA_SUCCESS) return fail("IDAGetSens failed");
      for (int k = 0; k < ns; ++k) {
        const sunrealtype* s = N_VGetArrayPointer(yS.v[k]);
        std::vector<double> col(nd);
        for (size_t j = 0; j < nd; ++j) col[j] = s[out.differential_index[j]];
        out.jv.push_back(std::move(col));
      }
      carrier->Accumulate(ida);
      out.sens_res_evals = carrier->res_evals_;
      out.sens_nonlin_iters = carrier->nonlin_iters_;
      out.sens_conv_fails = carrier->conv_fails_;
      out.reinits = carrier->reinits_;
      solver.SetReInitObserver(nullptr);
    }

    const IntegrationStats st = solver.GetStatistics();
    out.steps = st.num_steps - stats_start.num_steps;
    out.nonlin_iters = st.num_nonlin_iters - stats_start.num_nonlin_iters;
    out.jac_evals = st.num_jac_evals - stats_start.num_jac_evals;
    out.err_test_fails = st.num_err_test_fails - stats_start.num_err_test_fails;
    out.conv_fails = st.num_conv_fails - stats_start.num_conv_fails;
    out.residual_calls = data.profile_.residual_calls - residual_calls_start;
    out.wall_s = std::chrono::duration<double>(std::chrono::steady_clock::now() - measure_start).count();
    out.ok = true;
    return out;
  }

} // namespace adsim
