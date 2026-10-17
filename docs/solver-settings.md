# Solver Settings

Once the fluid, adsorbents, column and [cycle](cycles.md) are defined, the deck can be run. Skarstrom integrates the coupled [governing equations](equations.md) in time as a system of differential-algebraic equations (DAEs) using a variable-order, variable-step implicit method (SUNDIALS IDA), taking the bed through every step of every cycle.

## ## Running, stopping and monitoring

- **Run.** `sim.run()` (Python), `skarstrom run deck.json` (CLI) or `adsim::Simulation::Run` (C++) validates the deck, evaluates every expression and integrates all cycles. If validation fails the run returns status `failed` with the reason in `message`; nothing is integrated.
- **Progress and cancel.** `sim.run(progress=f)` calls `f(info)` after every step with the cycle, step and completed fraction; return `False` to cancel. Ctrl-C cancels a Python run cleanly.
- **Status.** One of `completed`, `reached_css` (reached [CSS](cyclic-steady-state.md) before the cycle limit), `failed` or `cancelled`.
- **Solver log.** `sim.run(log_level="medium", on_message=f)` passes the solver's messages to `f(level, text)`; `skarstrom run --log medium` prints them. *Low* reports essential milestones, *medium* adds per-cycle and per-step progress and CSS checks, *high* adds diagnostics. The run result carries the runtime and solver statistics (steps taken, residual evaluations, non-linear iterations, error-test and convergence failures).

Solver settings

The solver parameters are stored in the deck and are parameters of the registry (`settings.*`; `skarstrom params deck.json --all` lists them with their values). The defaults are robust for typical PSA/TSA problems; adjust them only when a run struggles to converge or you need tighter accuracy. The [Numerics](numerics.md) page explains what each one does inside the integrator.

| Setting | Parameter | Default | Meaning |
|---------|-----------|---------|---------|
| **Maximum Time Step** | `settings.max_step` | 0 (no limit) | Largest step the integrator may take; 0 leaves the step size to IDA's error control. Lowering it can stabilise stiff, fast-transient steps. Results are logged from IDA's interpolant, so the logged trajectory does not depend on it. Decks saved with the old default of 2 s load with no limit. |
| **Minimum Time Step** | `settings.min_step` | 1×10⁻¹⁵ s | Smallest permitted step. |
| **Initial Step Size** | `settings.init_step` | 1×10⁻¹⁰ s | First step attempted. |
| **Relative Tolerance** | `settings.reltol` | 1×10⁻⁶ | Relative accuracy target per variable. Tighten for more accuracy, loosen to speed up a hard run. |
| **Maximum BDF Order** | `settings.max_order` | 5 | Highest order of the backward-differentiation formula (1–5). |
| **Maximum Nonlinear Iterations Per Step** | `settings.max_nonlin_iters` | 15 | Newton iterations allowed per step before the step is retried smaller. |
| **Nonlinear Convergence Coefficient** | `settings.nonlin_conv_coef` | 0.1 | Safety factor on the Newton convergence test (0–1). |
| **Maximum Convergence Failures** | `settings.max_conv_fails` | 200 | Non-linear-solver failures tolerated before giving up. |
| **Maximum Error Test Failures** | `settings.max_err_test_fails` | 40 | Local error-test failures tolerated before giving up. |
| **Maximum Number of Steps** | `settings.max_num_steps` | 100000 | Cap on total integrator steps. |
| **Maximum Number of Cycles to Log** | `settings.max_cycles_to_keep` | 10 | How many recent cycles of full profile data are retained for plotting. |

> **Absolute tolerances** are set internally on a per-variable basis (10⁻⁸ for concentrations, mole fractions and loadings, 10⁻⁶ for pressure and velocity, 10⁻⁴ K for temperatures), so only the single relative tolerance is a setting. See [Numerics](numerics.md#time-integration-sundials-ida).

## Troubleshooting convergence

- **The run fails early with a message about consistent initial conditions.** The solver could not reconcile the [initial conditions](cycles.md#initial-conditions) with the first step's [boundary conditions](boundary-conditions.md). Make the initial pressure, temperature and composition physically consistent with what the first step imposes.
- **The run stalls or reports many convergence/error-test failures.** Reduce `settings.max_step` and/or loosen `settings.reltol` slightly; check for very stiff inputs (e.g. an extremely large mass-transfer coefficient or a near-instant pressure ramp). Increasing the layer's number of cells gradually rather than in one large jump also helps.
- **Results look under-resolved.** Increase the layer's number of cells and re-run; confirm the results stop changing.

Once a run completes, retain enough cycles (`settings.max_cycles_to_keep`) to inspect the approach to [cyclic steady state](cyclic-steady-state.md) and to compute [performance metrics](performance-metrics.md).
