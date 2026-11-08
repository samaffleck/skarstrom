# Cyclic Steady State

A PSA/TSA process does not reach a steady state in the usual sense — instead it settles into a **cyclic steady state (CSS)**, where every cycle is identical to the one before. Performance metrics such as purity and recovery are only meaningful once CSS is reached, because the first cycles still carry the memory of the (arbitrary) [initial condition](cycles.md#initial-conditions).

Skarstrom can detect CSS automatically and stop the run early, so you do not have to guess how many cycles are needed. Configure it with the deck's **CSS convergence method** (parameter `css.method`, an integer) and the settings below (`css.min_cycles`, the tolerances). The [**Number of Cycles**](cycles.md#number-of-cycles) still acts as an upper limit; CSS detection simply stops the run sooner once the criterion is met, and the solver status becomes **Cyclic Steady State**.

## Convergence methods

### None

No CSS test — the solver runs the full **Number of Cycles**. Use this when you want a fixed number of cycles regardless.

### Molar Profile

Compares the **total molar-concentration profile** over all cells at the end of the **first step** of the current cycle with the same profile one cycle earlier. Convergence is declared when

$$
\lVert \mathbf{c}_{\text{curr}} - \mathbf{c}_{\text{prev}} \rVert \thickspace \le\thickspace \varepsilon_{abs} + \varepsilon_{rel}\thinspace \lVert \mathbf{c}_{\text{curr}} \rVert
$$

where the norm is taken over all cells. You set the **Molar profile absolute tolerance** $\varepsilon_{abs}$ (default 1×10⁻⁴) and the **Molar profile relative tolerance** $\varepsilon_{rel}$ (default 1×10⁻⁴). This is a direct, profile-based test of whether the bed state has stopped changing from cycle to cycle.

### Mass Balance (Total)

Checks that, over a whole cycle, the moles entering the bed equal the moles leaving:

$$
\frac{\lvert n_{\text{in}} - n_{\text{out}} \rvert}{n_{\text{in}}} \thickspace \le\thickspace \varepsilon_{mb}
$$

summed over all steps and all faces, **and** that the bed itself has stopped changing (the [state test](#state-test) below). You set the **Mass balance tolerance** $\varepsilon_{mb}$ (default 0.005, i.e. 0.5%).

### Mass Balance (Per-Component)

The same as above but with the balance evaluated **per component**, using the worst-offending component:

$$
\max_i \frac{\lvert n_{\text{in},i} - n_{\text{out},i} \rvert}{n_{\text{in},i}} \thickspace \le\thickspace \varepsilon_{mb}
$$

Components that are never fed are skipped. Together with the state test this is the safest choice when component-level purity and recovery are the outputs you care about.

### State test

On its own, a closing balance is not a safe sign of CSS. On the way to CSS the bed alternately gains and loses a component over many cycles, so the balance error swings through zero, and a cycle near a crossing passes a 0.5 % band by accident. In a five-step CO₂/N₂ cycle on 13X (Yancy-Caballero et al., 2020) the CO₂ balance was inside 0.5 % from the fifth cycle on, while purity and recovery were still moving by up to 0.6 points; they settled only after 80 to 90 cycles. A much tighter balance tolerance is not the answer either: the balance carries small systematic errors (how the boundary fluxes are reconstructed and integrated over the logged points) that a tight band would never get past.

The mass-balance methods therefore also require the **end-of-cycle bed state** to have stopped changing. At the end of the last step of each cycle the solver compares the slow fields with the previous cycle's:

$$
\Delta_{\text{state}} = \max_f \frac{\operatorname{RMS}\left(x_f^{\,c} - x_f^{\,c-1}\right)}{s_f} \thickspace \le\thickspace \varepsilon_{\text{state}}
$$

over the adsorbed loading of each component (scaled by the RMS total loading $s_f$, so that a trace component's round-off cannot dominate; $s_f$ is at least the solver's absolute tolerance on loadings, 10⁻⁸ mol/kg, so the round-off loadings of an inert packing do not dominate either), the gas mole fractions ($s_f = 1$) and the fluid, solid and wall temperatures (scaled by their RMS). The loadings hold almost all of a component's inventory, so this measures the drift that the balance can hide. You set the **State change tolerance** $\varepsilon_{\text{state}}$ (`css.state_tolerance`, default 2×10⁻⁵). Set it to 0 to switch the state test off. The state test does not depend on the balance closing exactly, so a constant offset in the balance does not delay it.

### Consecutive cycles

Whatever the method, its test must pass on **Consecutive cycles** cycles in a row (`css.consecutive_cycles`, default 3) before CSS is declared, so that one quiet cycle cannot stop the run. The solver keeps at least as many cycles in its log as the check reads (one more than this number, or this number for the balance-only test), even when **Maximum Number of Cycles to Log** is smaller.

### Minimum cycles

Whatever the method, no convergence is declared before **Minimum cycles** cycles (default 5) have run. Raise it if an early cycle can satisfy the test by accident; lower it to let a converged run stop sooner.

### What the defaults deliver

Stopping cycle (0-based) and the largest distance of the tracked purity and recovery from their converged values (the last cycle of a run of 150 to 1000 cycles with the CSS check off), for decks that use **Mass Balance (Per-Component)**; *old rule* is the single-cycle balance test used before the state test existed (`state_tolerance` 0, one cycle). The Yancy-Caballero cycle is the five-step 13X case of that paper rebuilt in Skarstrom, compared at cycle 150. `Nguyen2023_Exp_X` is still drifting after 400 cycles and stops under neither rule. These runs used the example decks as they were up to engine 0.4 (upwind advection on 30 cells, Nguyen2023 on 100); the current decks run WENO on grid-converged grids, and the Perez2019 and Nguyen2023 decks start from a feed-loaded bed (the feed gas with the loadings in equilibrium with it), which reaches the same CSS in fewer cycles (Perez2019 A-F after 49, 57, 65, 94, 206 and 517 cycles); the stopping cycles of each paper deck are given in its validation README.

| Deck | Old rule | Error | Default rule | Error |
|---|---|---|---|---|
| Yancy-Caballero 2020, C1 (10 cells) | 4 | 1.3×10⁻² | 59 | 1.3×10⁻⁶ |
| Yancy-Caballero 2020, C2 (10 cells) | 4 | 4.1×10⁻³ | 100 | 1.4×10⁻⁵ |
| `Perez2019_A` | 43 | 1.3×10⁻³ | 69 | 3.1×10⁻⁵ |
| `Perez2019_B` | 54 | 4.8×10⁻⁴ | 80 | 2.1×10⁻⁵ |
| `Perez2019_C1` | 70 | 2.4×10⁻⁴ | 96 | 1.6×10⁻⁵ |
| `Perez2019_D1` | 154 | 5.9×10⁻⁴ | 189 | 3.8×10⁻⁵ |
| `Perez2019_E` | 351 | 6.2×10⁻⁴ | 425 | 4.7×10⁻⁵ |
| `Perez2019_F` | 673 | 1.8×10⁻³ | 839 | 1.3×10⁻⁴ |
| `Nguyen2023_Exp_I` | 23 | 5.5×10⁻⁴ | 44 | 4.2×10⁻⁶ |
| `Nguyen2023_Exp_V` | 38 | 2.7×10⁻⁴ | 61 | 4.3×10⁻⁶ |
| `Ward2024` | 54 | 4.1×10⁻⁴ | 152 | 8.3×10⁻⁶ |
| `Illustrative_TSA_13X` | 8 | 1.1×10⁻³ | 15 | 2.3×10⁻⁶ |
| `Ferre2026_NaYF` (no KPI configured) | 29 | — | 123 | — |
| `Vysyaraju2024_PSA_exp1` (no KPI configured) | 59 | — | 96 | — |

The error left at the stop grows with how slowly the deck approaches CSS: if each cycle closes a fraction $1-\rho$ of the remaining gap, what is left when the per-cycle change reaches the tolerance is proportional to $1/(1-\rho)$. `Perez2019_F` closes only about 1.5 % of its gap per cycle, which is why it ends furthest away; lower `css.state_tolerance` for such a deck. Keep it well above the noise of the state change itself, which is set by the integrator tolerance: about 10⁻⁶ at the default `settings.reltol` of 10⁻⁶ and up to 1.3×10⁻⁵ at 10⁻⁴ (both measured on the Yancy cycle). A tolerance below that floor stops the run only by chance. The verification suite checks this relation on a cycle with $\rho = 0.966$: the error left at the default stop matched the geometric prediction to within 2 % ([Numerics](numerics.md#convergence-of-the-results)).

### Estimated distance to CSS

The stopping rule bounds the size of the last cycle-to-cycle change, not the distance still to go. With a geometric approach at rate $\rho$, a run whose last change is $\Delta$ is still about

$$
d \approx \Delta\,\frac{\rho}{1-\rho}
$$

from the cyclic steady state, so a deck with $\rho$ near 1 stops several tolerances away. Every CSS check therefore reports an estimate of $d$, in the norm of the [state test](#state-test), and $d/\varepsilon_{\text{state}}$. It is a diagnostic only: it never changes when a run stops, and it changes no result.

The estimate uses the state changes $\Delta$ of the last few cycles:

- **Window.** The longest run of strictly decreasing, positive changes that ends at the cycle checked, at most 5 changes long. The cycle that starts from an [accelerated](#accelerating-the-approach-to-css) state, and every cycle before it, are left out: a jump is not a step of the natural cycle. Cycles no longer in the log are left out too.
- **Rate.** $\rho = (\Delta_\text{last}/\Delta_\text{first})^{1/(n-1)}$ over the $n$ changes of the window (the geometric mean of their ratios), clipped to $[0, 0.999]$.
- **Not estimated** when the window has fewer than 3 changes. The note says why: *history too short* (fewer than 3 changes logged), *too few cycles since the last acceleration jump*, or *not a monotone decay* (for example just after a jump, or once the changes have reached the integrator's noise floor, about 10⁻⁶ at the default `settings.reltol`, where they wander instead of decaying).

$d/\varepsilon_{\text{state}}$ is reported for the mass-balance methods with the state test on. For **Molar Profile** the distance is still given in the state-change norm, without the ratio.

Where it appears:

- **Solver log.** Each CSS check line ends with, for example, `est. distance to CSS 7.69e-04 (38.5 x tolerance, rate 0.976)` (`Perez2019_E` at its stop). At the stop, a run whose estimate exceeds **100 tolerances** (`CyclicSteadyState::kDistanceWarningRatio`) also logs a warning (always shown). Lower `css.state_tolerance` for such a deck if the state itself matters, keeping it above the noise floor.
- **Python.** `sim.css_history()` records have `distance_estimated`, `convergence_rate`, `distance`, `distance_ratio` (NaN when not estimated) and `distance_note`.
- **CLI.** The `css` entries of `skarstrom run --json` carry the same fields (`null` when not estimated), `css.csv` has the numbers, and the text summary prints the estimate at the stop.

On the decks of the [long-reference verification](#verification-against-a-long-reference), the estimate at the plain stop against the measured distance of that stop from the reference (both in multiples of $\varepsilon$, state):

| Deck | CSS stop | Rate $\rho$ | Estimated $d/\varepsilon$ | Measured $c$ state |
|---|---|---|---|---|
| `Illustrative_TSA_13X` | 26 | 0.57 | 0.30 | 0.27 |
| `Perez2019_A` | 48 | 0.86 | 4.0 | 4.0 |
| `Perez2019_E` | 205 | 0.976 | 38 | 37 |
| `Nguyen2023_Exp_I` | 36 | 0.80 | 2.1 | 2.1 |
| `Nguyen2023_Exp_VI` | 79 | 0.90 | 7.5 | 7.4 |
| `Vysyaraju2024_PSA_exp1` | 72 | — | not estimated (not a monotone decay) | 5.7 |
| `Ward2024` | 153 | 0.959 | 21 | 21 |
| `FiveStepProductPurge` (per-component rule) | 94 | 0.92 | 9.5 | 9.2 |

On these decks the estimate matches the measured distance to within a few per cent, though the reference also counts the pressure. The slow-tier test `CSSDistanceDecks.Perez2019EReportsTheRemainingDistance` checks the `Perez2019_E` figure.

### Decks saved before the state test

The state test and the consecutive-cycles rule were added after Skarstrom 0.3.0 (CSS settings class version 2). A deck saved before them loads with the defaults above, so a run that used **Mass Balance (Total)** or **Mass Balance (Per-Component)** now runs on until the bed state has settled: later than before, and with different (converged) purity and recovery. A **Molar Profile** run stops at least two cycles later than before, because its test must now hold on three cycles in a row. To get the old single-cycle balance test back, set `css.state_tolerance` to 0 and `css.consecutive_cycles` to 1.

## Accelerating the approach to CSS

The approach to CSS is often a slow geometric decay. Setting **CSS acceleration method** to *Aitken (Graves-Morris)* makes the solver extrapolate the end-of-cycle state of the slow fields (loadings, solid and wall temperatures, mole fractions) from three successive cycles towards the fixed point, check the jump in a throwaway integrator, and only then apply it. Pair it with a CSS method, otherwise the run still does every cycle. It never makes a run fail: a jump that does not integrate is skipped, and so is a window of three cycles that is not a geometric decay (the log says *CSS acceleration skipped* and why). In the [Perez 2019 tutorial](perez-2019-vpsa.md) it brought CSS from cycle 94 to 48. The algorithm is described on the [Numerics](numerics.md#cyclic-steady-state-detection-and-acceleration) page.

*Anderson* moves the same slow fields along a combination of up to four earlier cycle-to-cycle steps instead of one, so it can remove several slowly decaying modes in one jump where Aitken removes one. It jumps at most every other cycle, goes through the same throwaway-integrator check, and skips a jump whose fit would remove less than half of the remaining cycle-to-cycle change. On the shipped decks, with their own CSS settings (cycle of the CSS stop; errors of CO2 purity and recovery at the stop against a long plain run with the CSS check off, cycle 499 for Perez 2019 A and 999 for Perez 2019 E):

| Deck | Method | CSS stop | Purity error | Recovery error |
|---|---|---|---|---|
| Perez 2019 A | none | 48 | 3.2×10⁻⁷ | 2.4×10⁻⁵ |
| | Aitken | 32 | 5.0×10⁻⁷ | 1.4×10⁻⁵ |
| | Anderson | 35 | 7.9×10⁻⁷ | 3.4×10⁻⁵ |
| Perez 2019 E | none | 205 | 7.1×10⁻⁵ | 6.9×10⁻⁶ |
| | Aitken | 65 | 3.9×10⁻⁵ | 8.8×10⁻⁶ |
| | Anderson | 45 | 1.1×10⁻⁶ | 1.5×10⁻⁵ |
| Nguyen 2023 X | none | over 800 | | |
| | Aitken | 548 | | |
| | Anderson | 190 | | |

The long reference runs themselves wander by about 10⁻⁷ (A) and 10⁻⁶ (E) in purity and 2×10⁻⁶ (A) and 2×10⁻⁵ (E) in recovery over their last few hundred cycles. Anderson gains most where the approach to CSS is slow (Perez 2019 E, Nguyen 2023 X); on a deck that already converges in a few dozen cycles (Perez 2019 A) it does about as well as Aitken.

### Newton–Krylov

!!! note "New in 0.9"
    The *Newton-Krylov* acceleration method (`settings.acceleration_method = "newton_krylov"`) and the `newton_*` and `floquet_multiplier` fields of `sim.css_history()`.

Aitken and Anderson extrapolate from cycles the run has already made. *Newton–Krylov* instead solves for the cyclic steady state directly. With x the bed state at the start of a cycle and G(x) the state one cycle later, the CSS is the root of

F(x) = G(x) − x = 0.

Newton's method solves (I − J) d = F(x) for a step d, with J = ∂G/∂x the cycle map's Jacobian, and moves to x + α d. J is never formed. The linear system is solved by GMRES, which only needs products J·v, and each product is a finite difference of the cycle map:

J·v ≈ (G(x + h v) − G(x)) / h,

one plain cycle per Krylov vector. The eigenvalues of I − J cluster near 1, except for the few slow Floquet modes that make the approach to CSS slow, so GMRES needs a handful of vectors per Newton step. The extra cycles run in a scratch copy of the simulation, started from the state the run holds at the end of a cycle. The run itself, its log and its KPIs only see the result: one jump, through the same projection and throwaway-integrator check as Aitken and Anderson.

**Why finite differences.** The [cycle-map Jacobian by forward sensitivities](#cycle-map-jacobian-by-forward-sensitivities-prototype) does not work across the step switches of a PSA cycle. A finite difference only ever runs plain cycles, so it does. Its increment is the crux: G carries the integrator's error, about `settings.reltol`, so h is scaled to that noise. The perturbation's RMS size, in the scaled state, is √reltol·(1 + RMS of the scaled state). That balances the difference's truncation error (of order h) against the noise divided by h, and gives J·v to about √reltol (10⁻³ at the default reltol). Newton is then inexact, and the linear solve need not be more accurate than that. GMRES stops at a residual η·‖F‖, with the forcing term η of Eisenstat and Walker (η = 0.9 (‖F_k‖/‖F_{k−1}‖)², between 10⁻³ and 0.5, starting at 0.1).

**The unknowns** are every differential field of the bed at its cells (loadings, mole fractions, total gas concentration, fluid, solid and wall temperatures), scaled as in the [state test](#state-test). The residual ‖F‖ is the state test's norm: the largest, over the fields, of the RMS scaled change. Unlike Aitken and Anderson, Newton–Krylov also moves the gas concentration and the fluid temperature. Holding them at the run's values made a solve converge in the scratch copy while the run did not: on Perez 2019 A, every following cycle moved by 2×10⁻³ again. The endpoint histories (connected streams, buffer tanks) are held during a solve and updated by the run's own cycles between solves.

**A solve** starts after the fourth cycle, from the end state of that cycle. Each Newton step is damped by the rule of Glover, Papathanasiou and Pini: when some entries would go negative, α is the smallest x_i/(2(x_i − x_i,new)) over them (half the way to zero), otherwise 1. A step is accepted when it reduces ‖F‖; otherwise it is halved once, and if that does not help either, the solve ends. The run then jumps to G of the best iterate. A solve stops at ‖F‖ ≤ 10·reltol, ten times the noise floor of the state norm, or after eight Newton steps. After a converged solve the run continues with plain cycles, and the stopping rule ends it a few cycles later. After a solve that did not converge, Anderson runs for five cycles. The next solve then runs its cycles at reltol/100 and stops at 10 times that. Without a cycle map (a C++ `CoupledSolver` run not started through `Simulation`), the method behaves as Anderson.

On the shipped decks, with their own grid, start and CSS rule (cycle equivalents are the logged cycles plus the solve's extra cycles; c is the distance from a long plain reference in multiples of the state tolerance ε = 2×10⁻⁵, as in the [reference check](#verification-against-a-long-reference)):

| Deck | Method | CSS stop | Cycle equivalents | Newton steps (Krylov vectors) | c state | c KPIs |
|---|---|---|---|---|---|---|
| Perez 2019 A | none | 48 | 48 | | 4.4 | 1.2 |
| | Anderson | 36 | 36 | | | |
| | Newton–Krylov | 8 | 35 | 3 (7, 14, 2) | 0.43 | 0.38 |
| Perez 2019 E | none | 195 | 195 | | 36 | 3.4 |
| | Anderson | 47 | 47 | | | |
| | Newton–Krylov | 8 | 48 | 6 (3, 5, 6, 7, 10, 2) | 3.7 | 3.4 |
| Nguyen 2023 X | none | over 800 | | | | |
| | Anderson | 190 | 190 | | | |
| | Newton–Krylov | 14 | 218 | 8 + 8, in two solves | — | — |

The CSS stop is the cycle at which the stopping rule fired (counted from 1). The Anderson stops are those of `CSSAndersonDecks` (`core/test/integration/test_css_anderson_decks.cpp`) on the current decks; the tables above were measured before the Perez 2019 decks held gas in the bed voids only.
The reference of Perez 2019 A drifts by 0.002 ε, that of Perez 2019 E by 0.21 ε. Nguyen 2023 X has no feasible plain reference: its slowest Floquet multiplier is about 0.998, so a plain run still drifts after a thousand cycles. Instead, 50 plain cycles continued from the Newton–Krylov end state move it by 0.47 ε (Perez 2019 A: 0.43 ε, Perez 2019 E: 3.5 ε). The first solve on Nguyen 2023 X reaches the iteration limit at ‖F‖ = 1.5×10⁻⁴ (Krylov vectors 3, 2, 4, 5, 6, 10, 13, 1). The second, after five Anderson cycles, runs at reltol 10⁻⁸ and takes ‖F‖ from 1.6×10⁻⁴ to 3.8×10⁻⁷ (Krylov vectors 24, 4, 29, 9, 30, 6, 23, 15), where that tolerance's noise floor stops it.

**What the measurements show.**

- **The finite-difference increment.** At the deck's own reltol (10⁻⁶), J·v to about 10⁻³ is enough on both Perez decks: the first solve converges, at a dominant multiplier of about 0.9 on Perez 2019 E. It is not enough on Nguyen 2023 X. A slow mode of multiplier μ puts an eigenvalue 1 − μ = 0.002 in I − J, as small as the error of J·v. The Newton direction along that mode is then wrong, and the steps stop reducing ‖F‖. With every solve at the deck's reltol, Nguyen 2023 X took 299 cycle equivalents, against 190 for Anderson. Running the solve's cycles at reltol/100 from the start costs time on the fast decks (Perez 2019 A, before its gas-holdup change: 39 s instead of 28 s, for one cycle equivalent more; Perez 2019 E: 51 instead of 65 cycle equivalents but the same 22 s; Nguyen 2023 X: 152 cycle equivalents, 56 s instead of 81 s). Hence the rule above: tighten only after a solve fails.
- **Krylov vectors per Newton step.** Two to fourteen on the Perez decks, about as many as there are slow modes, as the Floquet argument predicts. The last steps of a solve need more: near the noise floor GMRES has to resolve the slow mode itself. The slow mode of Nguyen 2023 X needs up to 30, the limit per step.
- **Cost from cold.** 35 to 48 cycle equivalents on the Perez decks and 218 on Nguyen 2023 X: about Anderson's cost on the two Perez decks, more on Nguyen 2023 X. The run stops at cycle 8 to 14, but the solves' extra cycles cost about as much as the cycles Anderson saves. With every solve kept at the deck's reltol, Nguyen 2023 X took 299 cycle equivalents and ended 7.6 ε (in 50 plain cycles) from its fixed point.
- **The main gain is the end state.** A solve stops at a residual r, which leaves about r/(1 − μ) to go along the slowest mode. On Perez 2019 A that is 0.4 ε from the reference, against 4.4 ε for the plain stop; on Perez 2019 E 3.7 ε, against 36 ε. The remaining distance shrinks with the solve's tolerance, that is with `settings.reltol`.
- **Below the noise floor.** For a deeper CSS, tighten `settings.reltol`: the solve's tolerance follows it. On Perez 2019 A (20 cells) at reltol 10⁻¹², one solve takes ‖F‖ from 5.0×10⁻³ through 5.0×10⁻⁴, 3.2×10⁻⁵, 1.7×10⁻⁶ and 1.8×10⁻⁸ to 6.0×10⁻⁹ (contraction factors 0.10, 0.064, 0.053, 0.011 down to the floor: superlinear, as inexact Newton should be), in 57 extra cycles. Below about 10⁻⁸ the integrator's absolute tolerances (10⁻⁸ mol/kg on the loadings, 10⁻⁴ K on the temperatures), which reltol does not change, set the floor: the solve stalls at 1.2×10⁻⁸ at reltol 10⁻¹⁰.
- **Warm starts.** Started from the CSS of a neighbouring design (feed step 5 % longer, coarse Perez 2019 A), the first solve needs two Newton steps (‖F‖ 1.8×10⁻⁴ → 2.4×10⁻⁵ → 7.9×10⁻⁶) and 15 extra cycles. This is the situation of the [warm-started optimiser](#warm-started-optimisation-and-studies) and of a gradient computed by perturbing a design.

**When to use it.** Use Newton–Krylov when the state at CSS matters: comparing designs that differ by little, finite-difference gradients, or a slow deck whose plain stop ends far from its CSS. From cold it costs about as much as Anderson, or somewhat more on a slow deck, and it ends closer to the CSS, much closer with a tighter reltol. For a quick look at a deck's KPIs, Anderson is as good. Each solve is reported in the solver log at level *medium* (*Newton-Krylov: |F| … → … in N Newton step(s) …*, then the residual history) and in `sim.css_history()`, at the cycle after which it ran:

```python
import skarstrom

sim = skarstrom.examples.Perez2019_A.load()
sim.settings.acceleration_method = "newton_krylov"
result = sim.run()
for rec in sim.css_history():
    if rec.newton_iterations >= 0:  # a solve ran after this cycle
        print(rec.cycle, rec.newton_iterations, rec.newton_cycles, rec.newton_krylov_vectors,
              rec.newton_residuals, rec.floquet_multiplier, rec.newton_outcome)
```

`newton_cycles` counts the solve's extra cycles (cycle equivalents). `newton_residuals` is ‖F‖ at each Newton iterate, the start first, in the norm of `state_change`. `floquet_multiplier` is the dominant Floquet multiplier of the cycle map: the largest |1 − θ| over the eigenvalues θ of GMRES's Hessenberg matrix, the Ritz values of I − J. Only the directions the Krylov basis explored count, so a solve that needed few vectors can miss the slowest mode (0.30 on the first solve of Nguyen 2023 X, 0.99 on the second). The C++ API is `AccelerationMethod::NEWTON_KRYLOV`, `AccelerateCyclicSteadyState::NewtonRecords()` and the generic solver `SolveFixedPointNewtonKrylov` in `adsim/solver/newton_krylov.h`.

## Restarting from a cycle state and warm starts

A completed run keeps its **final state**, the state at the end of its last cycle. It holds everything the next cycle depends on:

- the DAE unknowns: concentrations, temperatures, mole fractions, loadings, pressure and velocity, in every cell;
- every endpoint's logged history. A [connected](connecting-beds.md) inlet replays its partner's logged outlet, and when the partner runs later in the cycle that log comes from the previous cycle;
- the buffer-tank summaries a tank inlet reads from its source step.

So a run restarted from the state reproduces the cycle the original run would have done next. The match is to within the integrator tolerance (about 1e-5 relative on the example decks). Without the histories the restart misses by 20–25 % on a deck whose equalisation or buffer tank reads a later step.

```python
import skarstrom as sk

sim = sk.examples.Perez2019_A.load()
sim.run()                              # reaches CSS after 49 cycles
x = sim.final_state()                  # a skarstrom.CycleState

other = sk.examples.Perez2019_A.load()
other.run(initial_state=x)             # this run starts from x: CSS again after 4 cycles
other.set_initial_state(x)             # or: every later run starts from x

other.run_cycle(x)                     # exactly one cycle from x
gx = other.final_state()               # the cycle map applied to x
gx.differential                        # the differential unknowns as one vector
gx.field("Mole Fraction", "CarbonDioxide")
```

- `run_cycle()` without a state continues from the last run's final state, so repeated calls step through the cycles one by one. It leaves the deck's number of cycles unchanged.
- A warm-started run counts the cycles behind the state towards the **minimum cycles** floor. The consecutive-cycles rule still applies, and the first cycle has no earlier cycle in the log to compare with. With the defaults (3 consecutive cycles), a run from a state at CSS therefore stops after 4 cycles. On the example decks: Perez2019_A 49 → 4 cycles, Nguyen2023_Exp_VI 80 → 4, Illustrative_TSA_13X 16 → 4, with the same purity and recovery to 1e-5.
- After an edit, a warm start from the previous CSS saves less, because successive cycles still have to carry the bed to the new CSS (see [Warm starts after an edit](#warm-starts-after-an-edit)).
- A state that does not fit the deck is refused with the reason (`check_state`; `set_initial_state` raises `ValueError`). That happens when the components, layers or void regions differ, or the steps or endpoints differ. Different cell counts and step times are not reasons: the run adapts the state (see below). A different state layout, such as another energy-balance mode, fails the run with the reason.
- The state lives for the session only. It is never written to the deck file. It is cleared by the next run and by loading a deck; the warm-start state stays until `clear_initial_state()` or a deck load.

In C++ the same is `Simulation::FinalState`, `RunOptions::initial_state`, `SetInitialState`, `RunCycle` and `CheckState`, with the state in `adsim/solver/cycle_state.h`. The chunked run (`Start` / `Advance`) always starts from the initial conditions.

### Moving a state to another grid

A state taken on one grid can start a run on another grid of the same deck, with different cells per layer. The run transfers the state first:

- Inside each layer, every per-cell field is remapped by overlap-weighted cell averages on the layer's uniform cells. The remap keeps each field's layer average and never creates new extremes.
- The gas is remapped as component concentrations, and the fluid temperature weighted by the gas concentration (which keeps the gas's sensible heat). The composition and the pressure are then rebuilt from them, so the pressure is the plain cell average of the old pressure and has no new extremes.
- So every component's gas and adsorbed inventory in every layer is kept to rounding: fine to coarse, coarse to fine, and between uneven counts.
- The velocity is interpolated. It is only a starting guess: the solver's consistent-initial-condition step fixes it and the other algebraic variables, as it does for any warm start.
- The connection histories and buffer-tank summaries belong to steps, not cells, so they are carried over unchanged.
- Refining by a whole-number ratio and coarsening back gives the original state.

So a grid-refinement study can start each level from the previous level's CSS. On a two-layer test deck, the refined grid (12 + 8 cells) reached CSS in 4 cycles from the coarse grid's (6 + 4) CSS, against 22 from cold. The example decks are in the table below.

| Deck (cells) | Cold | Warm from the coarse CSS |
|---|---|---|
| Illustrative_TSA_13X (20 → 40) | 16 | 6 |
| Perez2019_A (20 → 40) | 49 | 14 |
| Nguyen2023_Exp_VI (20 → 40) | 78 | 43 |

The purity and recovery agree with the cold run to 5e-5.

In C++ the transfer is `RemapCycleState` in `adsim/solver/cycle_state.h`, also available as `Simulation::RemapState`.

### Changed step times

When a step's time has changed since the state was taken, the run stretches that step's logged histories onto the new duration. The shape stays the same in step-relative time, so a connected inlet replays the same profile over the new step. The total moles in its buffer-tank summaries are scaled with the duration, so the average flow stays the same. A step whose old or new time is zero has its histories cleared, and its consumers use their first-cycle values. The first cycle then corrects whatever the stretch got wrong.

`sim.prepare_state(x)` returns the state as the run would start from it, together with a note of what was done. The run puts the same note in its messages (*Initial state: ...*). In C++ this is `Simulation::PrepareState` / `PrepareCycleState`.

### Warm starts after an edit

A study changes the deck between runs. A warm start from the previous run's CSS reaches the new CSS with the same purity and recovery as a cold start of the edited deck, in fewer cycles. These are the example decks as shipped, each warm-started from the CSS of the unedited deck:

| Deck (cells, base CSS cycles) | Edit | Cold | Warm | Purity / recovery difference |
|---|---|---|---|---|
| Illustrative_TSA_13X (30 cells, 16) | Adsorption step +10 % | 15 | 9 | 1e-7 / 7e-7 |
| Illustrative_TSA_13X | feed CO₂ 0.15 → 0.20 | 15 | 10 | 1e-7 / 9e-7 |
| Illustrative_TSA_13X | cells × 1.5 (30 → 45) | 16 | 5 | 4e-7 / 3e-6 |
| Perez2019_A (80 cells, 49) | Feed step +10 % | 45 | 31 | 2e-6 / 7e-5 |
| Perez2019_A | feed CO₂ 0.15 → 0.20 | 50 | 45 | 2e-9 / 3e-5 |
| Perez2019_A | cells × 1.5 (80 → 120) | 49 | 6 | 6e-7 / 7e-6 |
| Nguyen2023_Exp_VI (80 cells, 80) | Feed 2 step +10 % | 68 | 49 | 1e-6 / 2e-5 |
| Nguyen2023_Exp_VI | connected Feed 1 + Pressurisation 50 → 55 s (histories stretched) | 76 | 45 | 1e-6 / 1e-5 |
| Nguyen2023_Exp_VI | feed CO₂ 0.15 → 0.20 | 52 | 48 | 4e-6 / 4e-5 |
| Nguyen2023_Exp_VI | cells × 1.5 (80 → 120) | 81 | 17 | 3e-6 / 1e-5 |

The saving is largest when the edit leaves the bed close to its old CSS, as with a grid refinement. An edit that moves the CSS a long way still leaves many cycles of slow approach. The largest differences come from where the stopping rule halts each run, not from the start.

### Warm-started optimisation and studies

An optimisation or a sensitivity study runs the same deck many times at nearby points, so each run can start from the CSS of an earlier one. The engine keeps those states in a `CycleStateCache` (`adsim/solver/state_cache.h`), keyed by the decision vector, and starts each new run from the state of the **nearest** point evaluated so far. The distance is Euclidean in the variables scaled to [0, 1] by their bounds; a sensitivity study uses the runs' probability coordinates. The nearest point, not the previous one: a global optimiser jumps across the whole box at first, so the previous candidate is often far away.

- **On by default** in `skarstrom.optimisation.optimise` / `OptimizeSimulation` and in sensitivity studies; `warm_start=False` (CLI `skarstrom sensitivity --no-warm-start`) runs every point cold.
- **The reported optimum is cold.** After the search, the optimum is re-run from the deck's initial conditions; that is the reported cost and KPIs, next to the search's own value (`search_cost`).
- **Only CSS states are cached.** A run that stops at the cycle limit ends in a state that depends on where it started, so it is not reused. A deck without a CSS check, or a fixed-cycle study that never reaches CSS, runs every point cold, exactly as before.
- **A state that does not fit is skipped.** The start goes through the warm-start path above, so changed step times and cell counts are fine. A cached state of another structure (components, layers, steps or endpoints) is refused, the run is cold and the evaluation's `note` says why. A warm-started run that fails is repeated cold.
- **Warm runs stop at the same CSS within the bound, not bit for bit.** With several worker processes in a study, each keeps its own states, so which run starts from which depends on scheduling.

Measured on three example decks as shipped (deck CSS settings, ε = `css.state_tolerance` = 2×10⁻⁵). The optimisation varies two step times over ±20 % of the deck's value, maximising purity + recovery of CO₂, with 20 evaluations (5 initial samples). The study is a Morris screening of the same two step times at ±10 %, 4 trajectories (12 runs) after the nominal run. Cycles and wall times are means per run, the first (cold) run included; wall times are from one Linux run with three jobs sharing four cores, so read them to about 15 %.

| Deck (cold CSS cycles) | Task | Cycles per run, cold → warm | Wall time per run, cold → warm | Largest KPI difference warm vs cold | Optimum: cold re-run vs search |
|---|---|---|---|---|---|
| Illustrative_TSA_13X (27) | optimisation | 18.2 → 9.9 | 6.1 → 3.3 s | 0.4 ε | 0.01 ε |
| | screening | 28.2 → 15.1 | 7.9 → 4.9 s | 0.15 ε | — |
| Perez2019_A (49) | optimisation | 51.5 → 23.9 | 55 → 26 s | 4.3 ε | 2.5 ε |
| | screening | 49.2 → 23.3 | 59 → 26 s | 4.1 ε | — |
| Nguyen2023_Exp_I (37) | optimisation | 37.8 → 17.0 | 38 → 15 s | 0.8 ε | 0.7 ε |
| | screening | 36.9 → 19.2 | 29 → 20 s | 0.7 ε | — |

The KPI difference is the largest difference in CO₂ purity or recovery between a warm and a cold run at the same point (for the optimisation, at the points both searches evaluated: the initial samples and any later coincidences). The last column compares the optimum's cold re-run with the value the warm-started search found for it.

- **A warm-started run needs about half the cycles of a cold one**: 1.8–2.2× fewer cycles per run, and 1.4–2.5× less wall time per run (the spread is mostly the shared machine). A whole 20-evaluation optimisation of `Perez2019_A` took 10 min instead of 19.
- **Everything is within the bound** used throughout this page ($c \le 100$): the largest difference is 4.3 ε, on `Perez2019_A`, the slowest-converging deck, as in [the verification](#verification-against-a-long-reference).
- **The searches end at the same optimum.** On `Illustrative_TSA_13X` both find the same corner of the box. On `Nguyen2023_Exp_I` the warm search's optimum is 0.4 s of feed time from the cold search's and its cold cost is 0.3 ε lower. On `Perez2019_A` the objective is flat near the optimum: the warm search stops 2.3 s of feed time from the cold search, with a cold cost 5.9 ε higher (1.2×10⁻⁴ in purity + recovery). Two searches whose costs differ by a few ε take different paths through the box; that is a property of the search, not an error of the warm start.
- **The screening's μ\* agree to 1 %** on every deck and output, so the ranking is the same. A run's noise from the stopping rule is a few ε at most; keep the effects being measured well above it (tighten `css.state_tolerance` for the study, see [Sensitivity Analysis](sensitivity.md#rules-that-keep-the-numbers-honest)).

## Verification against a long reference

`core/test/validation/test_css_reference.cpp` (slow tier) checks, deck by deck, that a run stopped at CSS has reached the cyclic steady state of a much longer run. It compares the end-of-cycle state, not only the KPIs. For each deck it makes:

- **plain**: a run to CSS with the deck's own rule;
- **Aitken** and **Anderson**: the same run with each acceleration method;
- **warm**: a warm start from the reference end state with every loading scaled by 0.8 and every temperature raised by 10 K, run to CSS from a cycle count of 0;
- **reference**: the plain stop continued with the CSS check off, for two segments of max(100, plain stop + 1) cycles each.

A CSS tolerance 100 times tighter than the default (2×10⁻⁷) would sit below the noise of the state change itself, about 10⁻⁶ at the default `settings.reltol`, and would stop only by chance. That is why the reference is a fixed long run instead. The distance between the ends of its two segments (the *drift*) shows how far the reference itself still is from the true CSS.

Each run is compared with the reference in three ways:

- **State**: the largest, over the fields, of the RMS difference divided by the scale of the [state test](#state-test). That is the RMS total loading for the loadings, 1 for the mole fractions, and the RMS value for pressure and the fluid, solid and wall temperatures.
- **Histories**: each endpoint's logged run, through the moles it passed (relative to the largest endpoint flow) and their average composition.
- **KPIs**: the purity and recovery of every tracked component, as absolute differences.

Each distance is divided by the stop's state tolerance $\varepsilon$ (2×10⁻⁵ on every deck) to give $c$: a stop at tolerance $\varepsilon$ must be within $c\thinspace\varepsilon$ of the reference. The test holds every deck and method to one bound, $c \le 100$, and requires a reference drift of at most $1\cdot\varepsilon$. It prints the measured $c$ of each:

| Deck | Run | CSS stop | $c$ state | $c$ histories | $c$ KPIs |
|---|---|---|---|---|---|
| `Illustrative_TSA_13X` | plain | 26 | 0.27 | 0.03 | 0.05 |
| | Aitken | 20 | 0.03 | 0.003 | 0.002 |
| | Anderson | 18 | 0.005 | 0.003 | 0.003 |
| | warm | 18 | 0.26 | 0.03 | 0.05 |
| `Perez2019_A` | plain | 48 | 4.0 | 1.3 | 1.2 |
| | Aitken | 32 | 1.4 | 0.46 | 0.71 |
| | Anderson | 35 | 4.7 | 1.0 | 1.7 |
| | warm | 59 | 4.1 | 1.3 | 1.3 |
| `Perez2019_E` | plain | 205 | 37 | 9.1 | 3.6 |
| | Aitken | 65 | 21 | 5.0 | 2.0 |
| | Anderson | 45 | 1.4 | 0.26 | 0.52 |
| | warm | 273 | 36 | 8.6 | 3.4 |
| `Nguyen2023_Exp_I` | plain | 36 | 2.1 | 0.38 | 0.25 |
| | Aitken | 22 | 2.0 | 0.36 | 0.24 |
| | Anderson | 18 | 1.1 | 0.13 | 0.43 |
| | warm | 41 | 2.2 | 0.41 | 0.20 |
| `Nguyen2023_Exp_VI` | plain | 79 | 7.4 | 2.5 | 0.14 |
| | Aitken | 35 | 3.2 | 1.2 | 0.15 |
| | Anderson | 33 | 3.3 | 1.1 | 0.08 |
| | warm | 81 | 7.0 | 2.4 | 0.40 |
| `Vysyaraju2024_PSA_exp1` (no KPI configured) | plain | 72 | 5.7 | 1.2 | — |
| | Aitken | 72 | 7.2 | 1.2 | — |
| | Anderson | 59 | 1.3 | 0.43 | — |
| | warm | 72 | 5.8 | 1.0 | — |
| `Ward2024` | plain | 153 | 21 | 1.4 | 0.42 |
| | Aitken | 125 | 16 | 1.1 | 0.64 |
| | Anderson | 102 | 14 | 0.77 | 0.56 |
| | warm | 150 | 21 | 1.4 | 0.41 |
| `FiveStepProductPurge` (no KPI configured; run with the per-component rule) | plain | 94 | 9.2 | 1.4 | — |
| | Aitken | 65 | 4.6 | 0.66 | — |
| | Anderson | 48 | 2.4 | 0.37 | — |
| | warm | 87 | 9.3 | 1.4 | — |

Reference drifts, as multiples of $\varepsilon$ (state / histories / KPIs):

| Deck | Drift |
|---|---|
| `Illustrative_TSA_13X` | 0.004 / 0.005 / 0.004 |
| `Perez2019_A` | 0.002 / 0.024 / 0.027 |
| `Perez2019_E` | 0.23 / 0.09 / 0.34 |
| `Nguyen2023_Exp_I` | 0.02 / 0.006 / 0.001 |
| `Nguyen2023_Exp_VI` | 0.04 / 0.05 / 0.02 |
| `Vysyaraju2024_PSA_exp1` | 0.17 / 0.17 / — |
| `Ward2024` | 0.05 / 0.005 / 0.005 |
| `FiveStepProductPurge` | 0.009 / 0.18 / — |

A $c$ below the drift is at the reference's resolution. That applies to the smallest TSA figures and to Anderson's KPIs on `Perez2019_E`.

What the table shows:

- **Every stop is at the reference's CSS, within the bound.** None is near $c = 100$. The furthest are the plain and warm stops of the slowly converging decks: `Perez2019_E` ($c$ 37) and `Ward2024` (21). With a geometric approach of rate $\rho$, the error left at a stop is about $\varepsilon\rho/(1-\rho)$ (see [What the defaults deliver](#what-the-defaults-deliver)), so a slow deck stops further away. Lower `css.state_tolerance` for such a deck when the state itself matters.
- **The bed state is the strictest measure.** The KPIs settle before the loadings do: on `Perez2019_E` the plain stop's loadings are 37 ε from the reference, its KPIs 3.6 ε.
- **Anderson ends closest on the slow decks.** On `Perez2019_E` it stops at cycle 45 within 1.4 ε, where the plain run needs 205 cycles and ends 37 ε away. Aitken stops sooner than the plain run on most decks and ends as close or closer, except on `Vysyaraju2024_PSA_exp1` (same stop, slightly further).
- **A warm start from a perturbed state reaches the same CSS as the plain run**, with practically the same $c$.

**Decks left out.** The test runs eight of the example decks. The others are left out:

- The breakthrough decks (Casas 2012, Nguyen 2022, Streb 2021, Ward 2022) run one step once and have no cyclic steady state.
- `Perez2019_F`, `Nguyen2023_Exp_VII` to `Exp_X` and the two Arora decks stop after 300 to 1550 cycles. A reference of three times that would take 10 to 40 minutes per deck.
- `Ferre2026_NaYF` takes about 1.2 s per cycle, so 15 minutes.
- `Perez2019_B` to `D1` and `Nguyen2023_Exp_II` to `V` are the same cycles as the decks above at other operating points.

## Independence from the starting state

`core/test/validation/test_css_start_independence.cpp` (slow tier) checks that the cyclic steady state does not depend on where a run starts. It uses the decks of the [reference check](#verification-against-a-long-reference). Each is run to CSS with its own rule from several starts:

- **shipped**: the deck as shipped.
- **clean**: a gas of the weakly adsorbed component alone. Loadings always start in equilibrium with the initial gas, so only that component is adsorbed. This is N₂ for the CO₂/N₂ decks and He for `Vysyaraju2024_PSA_exp1`.
- **feed**: the feed gas, so the bed starts saturated with the feed.
- **hot** and **cold**: the shipped deck's state after one cycle, with every temperature (fluid, solid, wall) raised or lowered by 40 K, given through the warm-start API.
- **hot deck** and **cold deck**: the deck's own initial temperatures (gas and solid, wall) edited by ±40 K.
- **neighbour**: the CSS of a deck with the same cycle at another operating point, also given through the warm-start API. The pairs are `Perez2019_A` from `Perez2019_B`, `Perez2019_E` from `Perez2019_D1`, `Nguyen2023_Exp_I` from `Exp_II` and `Exp_VI` from `Exp_VII`. The other decks have no such neighbour.

A start that coincides with the shipped initial conditions is not run again. That is the feed start of the Perez and Nguyen decks and the clean start of the others.

The hot-deck and cold-deck starts edit the deck's initial temperature. Earlier engines evaluated the gas heat capacities once at the initial-condition pressure and temperature and held them for the run, so editing the initial temperature also changed the model: `Ward2024` started 40 K colder that way converged to a different CSS (recovery 0.91746 against 0.91770; issue [#145](https://github.com/samaffleck/skarstrom/issues/145)). The heat capacities are now the ideal-gas $C_{p0}(T)$ at the local temperature ([Equation of State](equation-of-state.md#gas-heat-capacity)), so the initial temperature is only a start, and the hot and cold deck starts land where the hot and cold warm starts do. The initial composition sets no model parameter, so the clean and feed starts edit it directly.

Every start is compared with the shipped start's CSS, using the [same distances](#verification-against-a-long-reference) and the same bound, $c \le 100$. Two stops that approach the CSS from different sides can be up to the sum of their own distances from it apart.

| Deck | Start | CSS stop | $c$ state | $c$ histories | $c$ KPIs |
|---|---|---|---|---|---|
| `Illustrative_TSA_13X` | shipped (clean) | 26 | | | |
| | feed | 17 | 0.50 | 0.05 | 0.09 |
| | hot | 24 | 0.05 | 0.006 | 0.02 |
| | cold | 25 | 0.10 | 0.009 | 0.006 |
| | hot deck | 26 | 0.04 | 0.006 | 0.01 |
| | cold deck | 26 | 0.03 | 0.01 | 0.08 |
| `Perez2019_A` | shipped (feed) | 48 | | | |
| | clean | 68 | 8.1 | 2.6 | 2.3 |
| | hot | 62 | 8.5 | 2.8 | 2.5 |
| | cold | 64 | 0.28 | 0.06 | 0.20 |
| | hot deck | 65 | 8.5 | 2.8 | 2.5 |
| | cold deck | 68 | 0.04 | 0.05 | 0.20 |
| | neighbour (B) | 34 | 8.9 | 2.9 | 2.7 |
| `Perez2019_E` | shipped (feed) | 205 | | | |
| | clean | 689 | 73 | 18 | 7.0 |
| | hot | 209 | 74 | 18 | 7.1 |
| | cold | 213 | 0.41 | 0.21 | 1.1 |
| | hot deck | 291 | 74 | 18 | 7.0 |
| | cold deck | 219 | 0.56 | 0.14 | 0.36 |
| | neighbour (D1) | 174 | 0.14 | 0.20 | 0.97 |
| `Nguyen2023_Exp_I` | shipped (feed) | 36 | | | |
| | clean | 44 | 4.4 | 0.82 | 0.48 |
| | hot | 45 | 4.5 | 0.83 | 0.50 |
| | cold | 45 | 0.42 | 0.06 | 0.09 |
| | hot deck | 46 | 4.5 | 0.84 | 0.51 |
| | cold deck | 47 | 0.14 | 0.02 | 0.03 |
| | neighbour (II) | 29 | 4.3 | 0.79 | 0.48 |
| `Nguyen2023_Exp_VI` | shipped (feed) | 79 | | | |
| | clean | 103 | 14 | 4.8 | 0.53 |
| | hot | 86 | 15 | 5.0 | 0.51 |
| | cold | 86 | 0.06 | 0.04 | 0.10 |
| | hot deck | 91 | 15 | 5.0 | 0.54 |
| | cold deck | 88 | 0.22 | 0.07 | 0.12 |
| | neighbour (VII) | 78 | 15 | 5.0 | 0.54 |
| `Vysyaraju2024_PSA_exp1` (no KPI configured) | shipped (clean) | 70 | | | |
| | feed | 75 | 14 | 3.0 | — |
| | hot | 84 | 0.91 | 0.26 | — |
| | cold | 81 | 13 | 2.5 | — |
| | hot deck | 84 | 1.0 | 0.44 | — |
| | cold deck | 79 | 14 | 2.8 | — |
| `Ward2024` | shipped (clean) | 153 | | | |
| | feed | 115 | 43 | 2.7 | 0.93 |
| | hot | 184 | 0.13 | 0.01 | 0.06 |
| | cold | 167 | 41 | 2.7 | 0.77 |
| | hot deck | 185 | 0.47 | 0.03 | 0.04 |
| | cold deck | 167 | 42 | 2.7 | 0.76 |
| `FiveStepProductPurge` (no KPI configured; per-component rule) | shipped (clean) | 94 | | | |
| | feed | 76 | 0.19 | 0.20 | — |
| | hot | 105 | 0.08 | 0.14 | — |
| | cold | 90 | 18 | 2.6 | — |
| | hot deck | 106 | 0.52 | 0.22 | — |
| | cold deck | 88 | 18 | 2.7 | — |

What the table shows:

- **Every start reaches the same steady state.** No deck has a second CSS.
- **The starts fall into two groups, each about as far from the shipped start as the plain stop is from the [long reference](#verification-against-a-long-reference).** Stops that approach the CSS from the same side as the shipped start end close to it. Stops from the other side end about twice the reference distance away. On `Perez2019_E` the clean, hot, hot-deck and shipped starts are 73–74 ε apart; the plain stop is 37 ε from the reference. The cold, cold-deck and neighbour starts are within 0.6 ε. On `Ward2024` the feed, cold and cold-deck starts are 41–43 ε away; the plain stop is 21 ε from the reference.
- **Editing the deck's initial temperature changes only the start.** Each deck's hot-deck and cold-deck starts end on the same side of the CSS as its hot and cold warm starts, at distances from the shipped start that agree with theirs to about 1 ε.
- **In the KPIs the largest difference between starts is 7 ε** (`Perez2019_E`), about 1.4×10⁻⁴ in purity or recovery.
- **The neighbour start helps most on the slow deck.** `Perez2019_E` from `Perez2019_D1`'s CSS stops at cycle 174 instead of 205, within 0.14 ε of the shipped start's CSS.

## Grid sequencing

Does it pay to run a coarse grid to CSS first and warm-start the fine grid from it ([moving a state to another grid](#moving-a-state-to-another-grid))? `core/test/validation/test_grid_sequencing.cpp` (slow tier) measures the time to the fine-grid CSS three ways on six of the [reference decks](#verification-against-a-long-reference), each with its own stop rule and no acceleration:

- **cold**: the deck's own grid (N cells per layer) from its initial conditions;
- **2-level**: N/2 cells to CSS, then N cells from that state;
- **3-level**: N/4 → N/2 → N cells, each level from the last one's CSS.

The table gives the cycles run at each level (cells: cycles) and the wall time of each path. The times are from one Linux run with two decks running at a time, so read them to about 10 %. The test asserts no times, only that every level reaches CSS, that the final state agrees with the cold run's within the [same bound](#verification-against-a-long-reference) ($c \le 100$), and that the warm-started fine level needs fewer cycles than the cold run.

| Deck | Path | Cycles per level | Total cycles | Time [s] | Speed-up | $c$ state / KPIs vs cold |
|---|---|---|---|---|---|---|
| `Illustrative_TSA_13X` | cold | 30: 27 | 27 | 4.5 | | |
| | 2-level | 15: 27, 30: 13 | 40 | 4.1 | 1.1 | 0.49 / 0.08 |
| | 3-level | 8: 25, 15: 13, 30: 13 | 51 | 4.1 | 1.1 | 0.49 / 0.08 |
| `Perez2019_A` | cold | 80: 49 | 49 | 29.2 | | |
| | 2-level | 40: 49, 80: 9 | 58 | 18.4 | 1.6 | 2.7 / 1.4 |
| | 3-level | 20: 49, 40: 12, 80: 9 | 70 | 14.0 | 2.1 | 3.1 / 2.2 |
| `Perez2019_E` | cold | 80: 206 | 206 | 56.5 | | |
| | 2-level | 40: 204, 80: 36 | 240 | 38.6 | 1.5 | 14 / 3.7 |
| | 3-level | 20: 195, 40: 60, 80: 38 | 293 | 34.0 | 1.7 | 13 / 3.5 |
| `Nguyen2023_Exp_I` | cold | 80: 37 | 37 | 20.1 | | |
| | 2-level | 40: 36, 80: 10 | 46 | 12.2 | 1.6 | 0.51 / 0.20 |
| | 3-level | 20: 37, 40: 12, 80: 10 | 59 | 10.7 | 1.9 | 0.54 / 0.22 |
| `Vysyaraju2024_PSA_exp1` (no KPI configured) | cold | 80: 73 | 73 | 32.7 | | |
| | 2-level | 40: 73, 80: 36 | 109 | 31.7 | 1.0 | 12 / — |
| | 3-level | 20: 72, 40: 38, 80: 35 | 145 | 31.5 | 1.0 | 13 / — |
| `Ward2024` | cold | 100: 154 | 154 | 84.1 | | |
| | 2-level | 50: 150, 100: 48 | 198 | 58.9 | 1.4 | 40 / 0.68 |
| | 3-level | 25: 153, 50: 55, 100: 53 | 261 | 56.2 | 1.5 | 42 / 0.72 |

What the table shows:

- **A cycle costs about in proportion to the cells**: halving the grid makes a cycle 2.2-2.5 times cheaper (`Perez2019_A`: 0.60, 0.27 and 0.12 s per cycle on 80, 40 and 20 cells).
- **A coarse grid needs as many cycles to CSS as the fine one.** The approach to CSS is set by the cycle, not by the grid.
- **The warm-started fine level needs 17-49 % of the cold cycles**: 17-18 % on `Perez2019_A` and `E`, 27-31 % on `Nguyen2023_Exp_I` and `Ward2024`, about half on `Illustrative_TSA_13X` and `Vysyaraju2024_PSA_exp1`. On the TSA deck the fine level cannot go much below its 13 cycles: three consecutive passing cycles and the settling of the remapped fronts take most of them. The PSA deck still needs 36 cycles on the fine grid from the coarse CSS; it is the one deck where the warm start saves no time.
- **So sequencing saves 30-50 % of the time where the fine level needs less than about a third of the cold cycles**, and nothing where it needs half. A third level helps a little more (it makes the first CSS cheaper still), at most another 0.5 in speed-up; it never made a path slower.
- **Every path reaches the same CSS.** The largest distance from the cold run is 42 ε in the state (`Ward2024`), the same as the [cold and feed starts](#independence-from-the-starting-state) of that deck, which approach from the same side; in the KPIs at most 3.7 ε.

Recommendation:

- **For one run on one grid**, sequencing pays on decks that take more than about 40 cycles to CSS: expect 1.4-2× on the slower VSA decks. On a deck that converges in a few tens of cycles it saves nothing worth the extra step. It never changes the answer beyond the stop rule's own spread.
- **Use two or three levels, each half the last**, and keep each level's own stop rule. Going below about 20 cells on these decks gains little more, since the coarsest level costs little already.
- **A grid-refinement study should always warm-start each level from the previous level's CSS** (`RemapCycleState`). The study needs every level's CSS anyway, so the coarse levels cost nothing extra and each finer level needs only 17-50 % of its cold cycles. The last level of a study is the most expensive one, so this saves most of the study's time on the decks that need many cycles. Tighten `css.state_tolerance` as the study does; the agreement above scales with it.
- Sequencing is a two-line loop over `RunOptions::initial_state` (`sim.run(initial_state=coarse.final_state())` in Python, after setting the fine cell counts), so there is no separate grid-sequence option.

## Cycle-map Jacobian by forward sensitivities (prototype)

A Newton–Krylov solver for the cyclic steady state solves G(x) − x = 0, where x is the differential state at the start of a cycle and G(x) the state one cycle later. Each Krylov vector needs a product Φ_x·v of the cycle map's Jacobian with a direction v. One way to get it is to integrate the cycle once with IDAS forward sensitivities whose initial values are the directions. The step sequence and the matrix factorisations are then the base run's, and each direction adds only the sensitivity corrector solves. Glover, Papathanasiou & Pini (2026) imply about 0.1 cycle of extra work per direction. The engine has a prototype of this (`CycleJacobianVector`, `adsim/solver/cycle_sensitivity.h`, internal, C++ only), checked against a centred finite difference (G(x + hv) − G(x − hv))/2h of the plain cycle at reltol 10⁻¹⁰. The difference's own error is estimated from the steps h and 2h. Measured in `core/test/integration/test_cycle_sensitivity.cpp` (`CycleSensitivity.MeasureCostAndAccuracy`, slow tier).

**How it works.** The prototype runs the cycle's steps as the solver does, without logging, retries or the CSS check. It restarts the sensitivities at every step switch together with the integrator: they are read at the end-of-step point before `IDAReInit` drops IDA's history, then `IDASensReInit`, and the same `IDACalcIC` that makes the state consistent with the new step's boundary conditions also makes the algebraic sensitivities consistent. Sensitivity error control is off, so the step sequence is the base run's whenever the sensitivity corrector converges. The engine links IDAS rather than IDA for this. IDAS has the same entry points, and the end state and integrator counters of six example decks are bit-identical. Three things had to be set differently from IDAS's defaults:

- **The sensitivity residual F_y·s + F_y′·s′** is a difference quotient of the engine's residual whose increment is scaled to the state: ε^(1/3) relative to |y| + abstol for the centred quotient, ε^(1/2) for the forward one. IDAS's internal quotient perturbs by about reltol·|y|, and its rounding error is of order one in IDA's weighted norm.
- **The quotient uses its own scratch vectors.** With the simultaneous corrector, IDAS's `IDACalcIC` passes its trial state as y and the same vector again as scratch to the sensitivity residual function, so a residual function that writes its scratch overwrites y. This is a SUNDIALS issue (`idas_ic.c`, `IDAfnorm`).
- **The sensitivity tolerances** are the state's tolerances times the size of the direction. With the state's absolute tolerances alone, a component whose sensitivity starts at zero (the pressures, abstol 10⁻⁶ Pa) gets a weight no corrector can meet.

**Where it works.** On the single-step Casas 2012 breakthrough deck, started after one warm-up cycle so the measured cycle begins from a flowing bed, a smooth temperature bump gives a Φ_x·v that agrees with the finite difference to the difference's own error. On 10 cells the agreement is 7.5×10⁻⁴ against an error estimate of 7.1×10⁻⁴ (the fast-tier test). On the deck's grid (160 cells, 1440 unknowns) it is 2–3×10⁻³ against 2–4×10⁻³, but only at reltol 10⁻¹⁰. At the deck's reltol (10⁻⁶) the sensitivity corrector fails right after the cycle start for most directions.

**Cost.** Even where it converges, the cost is far from 0.1 cycle per direction. On Casas 2012 at its own grid (plain cycle 0.08 s, 1567 residual calls, 66 steps), the only runs that completed at the deck's tolerance were the 16-direction runs with the centred quotient:

| Corrector | Directions | Wall time / plain cycle | Extra wall per direction [cycles] | Extra residual calls per direction [cycles] | Steps | Sensitivity convergence failures |
|---|---|---|---|---|---|---|
| staggered | 16 | 52 | 3.2 | 2.4 | 247 | 186 |
| simultaneous | 16 | 43 | 2.6 | 1.9 | 250 | (not counted separately) |

With 1 and 4 directions, and with the forward quotient, the runs failed in the first step. Every sensitivity convergence failure makes IDA cut the step, so the sensitivities do not follow the base run's step sequence (247 steps against 66) and the saving the method relies on is lost.

**Why it fails.** The residual is not differentiable everywhere. The upwind direction switches with the sign of the velocity, the flux limiters (WENO, van Leer) and clamps have corners, and at a step switch (a valve opening, a flow reversing) the velocities pass through zero. Where a direction straddles such a corner, F_y·s from a difference quotient is not linear in s. The Newton iteration for the sensitivities then cycles instead of converging: on Perez 2019 A, after the first step cut the error by 10⁴, the iterates of the consistent-IC solve repeated with period 3 at a weighted norm of about 1, at every tolerance tried. The matrix itself is right: at every point checked it agrees with a centred difference of the residual to 10⁻¹⁴ (upwind scheme). The sensitivity IC system is well conditioned (condition number 9 after scaling). Started from the dense solution of that system, IDAS still fails. With WENO, even the finite difference of the cycle is noisy: h and 2h disagree by 1.5×10⁻² for a random direction on Perez 2019 A, against 2.5×10⁻⁵ with upwind. On Perez 2019 A (four steps) no sensitivity run completed a cycle, with any corrector, quotient or tolerance tried. A temperature bump fails at the Evacuation switch, a random direction at the cycle's first consistent-IC solve. From the deck's bed at rest the velocities sit exactly on the upwind switch at the first step's start, which is why the prototype has warm-up cycles.

**What it means for a Newton–Krylov CSS solver.** IDAS forward sensitivities with a difference-quotient residual do not give Φ_x·v on the PSA decks: they fail across the step switches. Where they do run, they cost several cycles per direction, more than a finite difference of the cycle (one cycle per direction, two if centred). The options are:

1. A **finite-difference J·v**, (G(x + hv) − G(x))/h, one plain cycle per Krylov vector. It needs no change to the integrator and works across switches, but needs a step h well above the integrator's noise (G at reltol 10⁻⁸ or tighter).
2. A sensitivity residual that is **linear in s by construction**: the product of an assembled F_y and F_y′ (coloured, as the solver's Jacobian already is) with s. Its corrector converges like the matrix it uses, at the cost of one extra coloured evaluation per step whatever the number of directions. That pays off only with many directions per cycle.
3. **Smoothing the corners** (a smooth upwind switch and limiter) so that sensitivities exist, which changes the model's results.

The prototype stays internal; it covers a single bed without connections or buffer tanks.

## Practical guidance

- Set a generous **Number of Cycles** (an upper bound) and let a CSS method stop the run when converged.
- The **Molar Profile** method is a quick check; use **Mass Balance (Per-Component)** whenever you report purity or recovery. The gas profiles and the total balance can settle long before the adsorbed inventory of one component does: in `Perez2019_E` the total balance closes within 0.5 % after 91 cycles, the CO₂ balance only after 449, and the CO₂ purity rises from 80 % to 87 % in between.
- For a study that compares small differences between runs (a [sensitivity study](sensitivity.md), an optimisation), tighten `css.state_tolerance` (the error left at the stop scales roughly with it), keeping it above the noise floor described above.
- `sim.css_history()` reports, for every logged cycle, the balance error and its tolerance (`error`, `tolerance`), the state change and its tolerance (`state_change`, `state_tolerance`), and how many consecutive cycles have passed (`passing_cycles` of `required_cycles`). The solver log at level *medium* prints the same after every cycle.
- If a run reaches the cycle limit without converging, either increase the limit or loosen the tolerance — and check that the cycle actually *can* close (e.g. that feed, product and purge flows are balanced).
- Keep enough cycles logged (**Maximum Number of Cycles to Log** in the [solver settings](solver-settings.md)) to plot the approach to CSS and confirm the last cycles overlay.
