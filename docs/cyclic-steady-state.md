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

over the adsorbed loading of each component (scaled by the RMS total loading $s_f$, so that a trace component's round-off cannot dominate), the gas mole fractions ($s_f = 1$) and the fluid, solid and wall temperatures (scaled by their RMS). The loadings hold almost all of a component's inventory, so this measures the drift that the balance can hide. You set the **State change tolerance** $\varepsilon_{\text{state}}$ (`css.state_tolerance`, default 2×10⁻⁵). Set it to 0 to switch the state test off. The state test does not depend on the balance closing exactly, so a constant offset in the balance does not delay it.

### Consecutive cycles

Whatever the method, its test must pass on **Consecutive cycles** cycles in a row (`css.consecutive_cycles`, default 3) before CSS is declared, so that one quiet cycle cannot stop the run. The solver keeps at least as many cycles in its log as the check reads (one more than this number, or this number for the balance-only test), even when **Maximum Number of Cycles to Log** is smaller.

### Minimum cycles

Whatever the method, no convergence is declared before **Minimum cycles** cycles (default 5) have run. Raise it if an early cycle can satisfy the test by accident; lower it to let a converged run stop sooner.

### What the defaults deliver

Stopping cycle (0-based) and the largest distance of the tracked purity and recovery from their converged values (the last cycle of a run of 150 to 1000 cycles with the CSS check off), for decks that use **Mass Balance (Per-Component)**; *old rule* is the single-cycle balance test used before the state test existed (`state_tolerance` 0, one cycle). The Yancy-Caballero cycle is the five-step 13X case of that paper rebuilt in Skarstrom, compared at cycle 150. `Nguyen2023_Exp_X` is still drifting after 400 cycles and stops under neither rule.

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

The error left at the stop grows with how slowly the deck approaches CSS: if each cycle closes a fraction $1-\rho$ of the remaining gap, what is left when the per-cycle change reaches the tolerance is proportional to $1/(1-\rho)$. `Perez2019_F` closes only about 1.5 % of its gap per cycle, which is why it ends furthest away; lower `css.state_tolerance` for such a deck. Keep it well above the noise of the state change itself, which is set by the integrator tolerance: about 10⁻⁶ at the default `settings.reltol` of 10⁻⁶ and up to 1.3×10⁻⁵ at 10⁻⁴ (both measured on the Yancy cycle). A tolerance below that floor stops the run only by chance. The verification suite checks this relation on a cycle with $\rho = 0.966$: the error left at the default stop matched the geometric prediction to within 8 % ([Numerics](numerics.md#convergence-of-the-results)).

### Decks saved before the state test

The state test and the consecutive-cycles rule were added after Skarstrom 0.3.0 (CSS settings class version 2). A deck saved before them loads with the defaults above, so a run that used **Mass Balance (Total)** or **Mass Balance (Per-Component)** now runs on until the bed state has settled: later than before, and with different (converged) purity and recovery. A **Molar Profile** run stops at least two cycles later than before, because its test must now hold on three cycles in a row. To get the old single-cycle balance test back, set `css.state_tolerance` to 0 and `css.consecutive_cycles` to 1.

## Accelerating the approach to CSS

The approach to CSS is often a slow geometric decay. Setting **CSS acceleration method** to *Aitken (Graves-Morris)* makes the solver extrapolate the end-of-cycle state of the slow fields (loadings, solid and wall temperatures, mole fractions) from three successive cycles towards the fixed point, check the jump in a throwaway integrator, and only then apply it. Pair it with a CSS method, otherwise the run still does every cycle. It never makes a run fail: a jump that does not integrate is skipped. In the [Perez 2019 tutorial](perez-2019-vpsa.md) it brought CSS from cycle 190 to 182. The algorithm is described on the [Numerics](numerics.md#cyclic-steady-state-detection-and-acceleration) page.

## Practical guidance

- Set a generous **Number of Cycles** (an upper bound) and let a CSS method stop the run when converged.
- The **Molar Profile** method is a quick check; use **Mass Balance (Per-Component)** whenever you report purity or recovery. The gas profiles and the total balance can settle long before the adsorbed inventory of one component does: in `Perez2019_E` the total balance closes within 0.5 % after 91 cycles, the CO₂ balance only after 449, and the CO₂ purity rises from 80 % to 87 % in between.
- For a study that compares small differences between runs (a [sensitivity study](sensitivity.md), an optimisation), tighten `css.state_tolerance` (the error left at the stop scales roughly with it), keeping it above the noise floor described above.
- `sim.css_history()` reports, for every logged cycle, the balance error and its tolerance (`error`, `tolerance`), the state change and its tolerance (`state_change`, `state_tolerance`), and how many consecutive cycles have passed (`passing_cycles` of `required_cycles`). The solver log at level *medium* prints the same after every cycle.
- If a run reaches the cycle limit without converging, either increase the limit or loosen the tolerance — and check that the cycle actually *can* close (e.g. that feed, product and purge flows are balanced).
- Keep enough cycles logged (**Maximum Number of Cycles to Log** in the [solver settings](solver-settings.md)) to plot the approach to CSS and confirm the last cycles overlay.
