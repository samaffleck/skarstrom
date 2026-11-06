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

- **shipped**: the deck as shipped;
- **clean**: a gas of the weakly adsorbed component alone. Loadings always start in equilibrium with the initial gas, so only that component is adsorbed. This is N₂ for the CO₂/N₂ decks and He for `Vysyaraju2024_PSA_exp1`;
- **feed**: the feed gas, so the bed starts saturated with the feed;
- **hot** and **cold**: the shipped bed and wall temperatures raised and lowered by 40 K;
- **neighbour**: the CSS of a deck with the same cycle at another operating point, given through the warm-start API. The pairs are `Perez2019_A` from `Perez2019_B`, `Perez2019_E` from `Perez2019_D1`, `Nguyen2023_Exp_I` from `Exp_II` and `Exp_VI` from `Exp_VII`. The other decks have no such neighbour.

A start that coincides with the shipped initial conditions is not run again. That is the feed start of the Perez and Nguyen decks and the clean start of the others.

Every start is compared with the shipped start's CSS, using the [same distances](#verification-against-a-long-reference) and the same bound, $c \le 100$. Two stops that approach the CSS from different sides can be up to the sum of their own distances from it apart.

| Deck | Start | CSS stop | $c$ state | $c$ histories | $c$ KPIs |
|---|---|---|---|---|---|
| `Illustrative_TSA_13X` | shipped (clean) | 26 | | | |
| | feed | 17 | 0.50 | 0.05 | 0.09 |
| | hot | 26 | 8.8 | 2.2 | 2.5 |
| | cold | 26 | 8.9 | 2.6 | 2.6 |
| `Perez2019_A` | shipped (feed) | 48 | | | |
| | clean | 68 | 8.0 | 2.6 | 2.4 |
| | hot | 65 | 14 | 3.2 | 0.79 |
| | cold | 68 | 6.9 | 0.55 | 1.6 |
| | neighbour (B) | 34 | 8.9 | 2.8 | 2.7 |
| `Perez2019_E` | shipped (feed) | 205 | | | |
| | clean | 689 | 73 | 18 | 7.0 |
| | hot | 291 | 74 | 18 | 7.1 |
| | cold | 219 | 2.4 | 0.34 | 0.76 |
| | neighbour (D1) | 174 | 0.14 | 0.10 | 0.01 |
| `Nguyen2023_Exp_I` | shipped (feed) | 36 | | | |
| | clean | 44 | 4.4 | 0.81 | 0.45 |
| | hot | 46 | 6.8 | 1.4 | 1.0 |
| | cold | 47 | 2.9 | 0.54 | 1.3 |
| | neighbour (II) | 29 | 4.3 | 0.78 | 0.44 |
| `Nguyen2023_Exp_VI` | shipped (feed) | 79 | | | |
| | clean | 103 | 14 | 4.8 | 0.46 |
| | hot | 91 | 16 | 5.2 | 0.34 |
| | cold | 88 | 3.0 | 0.29 | 0.27 |
| | neighbour (VII) | 78 | 15 | 5.0 | 0.56 |
| `Vysyaraju2024_PSA_exp1` (no KPI configured) | shipped (clean) | 72 | | | |
| | feed | 76 | 12 | 2.6 | — |
| | hot | 85 | 7.2 | 0.31 | — |
| | cold | 83 | 7.8 | 1.6 | — |
| `FiveStepProductPurge` (no KPI configured; per-component rule) | shipped (clean) | 94 | | | |
| | feed | 76 | 0.20 | 0.24 | — |
| | hot | 105 | 3.5 | 0.83 | — |
| | cold | 88 | 16 | 2.6 | — |
WARD_ROWS
What the table shows:

- **Every start reaches the same steady state.** No deck has a second CSS.
- **The distances between starts follow the approach rate, as in the reference check.** The decks that converge slowly end further apart. The furthest are on `Perez2019_E`: the clean start (stopped at cycle 689) and the hot start are 73 to 74 ε from the shipped start in the bed state. Their KPIs are 7 ε apart, about 1.4×10⁻⁴ in purity or recovery.
- **The neighbour start helps most on the slow deck.** `Perez2019_E` from `Perez2019_D1`'s CSS stops at cycle 174 instead of 205, within 0.14 ε of the shipped start's CSS.

WARD_TEXT
## Practical guidance

- Set a generous **Number of Cycles** (an upper bound) and let a CSS method stop the run when converged.
- The **Molar Profile** method is a quick check; use **Mass Balance (Per-Component)** whenever you report purity or recovery. The gas profiles and the total balance can settle long before the adsorbed inventory of one component does: in `Perez2019_E` the total balance closes within 0.5 % after 91 cycles, the CO₂ balance only after 449, and the CO₂ purity rises from 80 % to 87 % in between.
- For a study that compares small differences between runs (a [sensitivity study](sensitivity.md), an optimisation), tighten `css.state_tolerance` (the error left at the stop scales roughly with it), keeping it above the noise floor described above.
- `sim.css_history()` reports, for every logged cycle, the balance error and its tolerance (`error`, `tolerance`), the state change and its tolerance (`state_change`, `state_tolerance`), and how many consecutive cycles have passed (`passing_cycles` of `required_cycles`). The solver log at level *medium* prints the same after every cycle.
- If a run reaches the cycle limit without converging, either increase the limit or loosen the tolerance — and check that the cycle actually *can* close (e.g. that feed, product and purge flows are balanced).
- Keep enough cycles logged (**Maximum Number of Cycles to Log** in the [solver settings](solver-settings.md)) to plot the approach to CSS and confirm the last cycles overlay.
