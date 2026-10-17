# Cyclic Steady State

A PSA/TSA process does not reach a steady state in the usual sense — instead it settles into a **cyclic steady state (CSS)**, where every cycle is identical to the one before. Performance metrics such as purity and recovery are only meaningful once CSS is reached, because the first cycles still carry the memory of the (arbitrary) [initial condition](cycles.md#initial-conditions).

Skarstrom can detect CSS automatically and stop the run early, so you do not have to guess how many cycles are needed. Configure it with the deck's **CSS convergence method** (parameter `css.method`, an integer) and the settings below (`css.min_cycles`, the tolerances). The [**Number of Cycles**](cycles.md#number-of-cycles) still acts as an upper limit; CSS detection simply stops the run sooner once the criterion is met, and the solver status becomes **Cyclic Steady State**.

## Convergence methods

### None

No CSS test — the solver runs the full **Number of Cycles**. Use this when you want a fixed number of cycles regardless.

### Molar Profile

Compares the **total molar-concentration profile** over all cells at the end of the **first step** of the current cycle with the same profile one cycle earlier. Convergence is declared when

$$
\lVert \mathbf{c}_{\text{curr}} - \mathbf{c}_{\text{prev}} \rVert \;\le\; \varepsilon_{abs} + \varepsilon_{rel}\,\lVert \mathbf{c}_{\text{curr}} \rVert
$$

where the norm is taken over all cells. You set the **Molar profile absolute tolerance** $\varepsilon_{abs}$ (default 1×10⁻⁴) and the **Molar profile relative tolerance** $\varepsilon_{rel}$ (default 1×10⁻⁴). This is a direct, profile-based test of whether the bed state has stopped changing from cycle to cycle.

### Mass Balance (Total)

Checks that, over a whole cycle, the moles entering the bed equal the moles leaving:

$$
\frac{\lvert n_{\text{in}} - n_{\text{out}} \rvert}{n_{\text{in}}} \;\le\; \varepsilon_{mb}
$$

summed over all steps and all faces. You set the **Mass balance tolerance** $\varepsilon_{mb}$ (default 0.005, i.e. 0.5%). At CSS there is no net accumulation, so the closing mass balance is a robust convergence indicator.

### Mass Balance (Per-Component)

The same as above but evaluated **per component**, using the worst-offending component:

$$
\max_i \frac{\lvert n_{\text{in},i} - n_{\text{out},i} \rvert}{n_{\text{in},i}} \;\le\; \varepsilon_{mb}
$$

Components that are never fed are skipped. This is the strictest test and the safest choice when component-level purity and recovery are the outputs you care about.

### Minimum cycles

Whatever the method, no convergence is declared before **Minimum cycles** cycles (default 5) have run. Raise it if an early cycle can satisfy the test by accident; lower it to let a converged run stop sooner.

## Accelerating the approach to CSS

The approach to CSS is often a slow geometric decay. Setting **CSS acceleration method** to *Aitken (Graves-Morris)* makes the solver extrapolate the end-of-cycle state of the slow fields (loadings, solid and wall temperatures, mole fractions) from three successive cycles towards the fixed point, check the jump in a throwaway integrator, and only then apply it. Pair it with a CSS method, otherwise the run still does every cycle. It never makes a run fail: a jump that does not integrate is skipped. In the [Perez 2019 tutorial](perez-2019-vpsa.md) it brought CSS from cycle 182 to 176. The algorithm is described on the [Numerics](numerics.md#cyclic-steady-state-detection-and-acceleration) page.

## Practical guidance

- Set a generous **Number of Cycles** (an upper bound) and let a CSS method stop the run when converged.
- The **Molar Profile** method is a good general default; switch to **Mass Balance (Per-Component)** when you need confidence in component-level KPIs.
- If a run reaches the cycle limit without converging, either increase the limit or loosen the tolerance — and check that the cycle actually *can* close (e.g. that feed, product and purge flows are balanced).
- Keep enough cycles logged (**Maximum Number of Cycles to Log** in the [solver settings](solver-settings.md)) to plot the approach to CSS and confirm the last cycles overlay.
