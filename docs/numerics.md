## Numerics

This page describes how the solver turns the [governing equations](equations.md) into numbers: the spatial discretisation, the time integrator and its settings, how the Jacobian is built, and how cyclic steady state is detected and accelerated. Everything here is taken from the code: `core/src/equations/adsorbent_equations.cpp` (residuals), `core/src/equations/equation.cpp` (unknowns and layout), `core/src/solver/coupled_solver.cpp` and `core/include/adsim/solver/ida_settings.h` (integration), `core/src/solver/jacobian.cpp` (Jacobian) and `core/src/solver/cyclic_steady_state.cpp` (CSS).

---

### Spatial discretisation: cell-centred finite volumes

The column is divided into cells along its axis. Each layer contributes its **Number of cells** of equal width $\Delta z = L/N$; optional bottom and top [void regions](getting-started.md) add one lumped cell each; and one **ghost cell** at each end carries the boundary condition. Layers of different cell widths can sit next to each other.

The cell-centred unknowns, in every cell, are:

| Unknown | Residual it carries |
|---|---|
| total concentration $C$ | ideal gas law, $C - P/(RT) = 0$ ([Equation of State](equation-of-state.md)) |
| pressure $P$ | overall mass balance ([Mass Balance](mass-balance.md)) |
| gas temperature $T$ | bed energy balance ([Energy Balance](energy-balance.md)) |
| solid temperature $T_s$ | solid energy balance in LTNE; $T_s = T$ in LTE |
| wall temperature $T_w$ | [wall energy balance](wall-energy-balance.md) |
| mole fraction $y_i$ of every component | component mass balance |
| adsorbed concentration $q_{i,k}$ of every component on every adsorbent $k$ | the [kinetic rate model](kinetic-models.md), only in the cells of layers that contain adsorbent $k$ |

In an [isothermal run](energy-balance.md#isothermal-mode) the three temperatures are removed from the unknowns, and void cells carry no solid temperature or loadings.

**Where results are reported.** A ghost-cell value is not a value anywhere in the bed (an inlet temperature ghost holds $2T_{in} - T_1$, so that the face between it and the first cell sits at $T_{in}$), so the results never report one. Axial profiles of cell quantities are given at the cell centres, and of face quantities (velocity, flow rates) at the faces. The bed ends report the physical boundary value: the composition, temperature and pressure an inlet imposes, or the bed-end face value at an outlet, which is what the boundary-stream logs record ([Results](python.md#results), `adsim/solver/bed_series.h`).

**Staggered face velocities.** The superficial velocity $u$ is not stored at cell centres but on the **faces** between neighbouring cells, including the faces between the ghost cells and the bed. Each face carries the momentum balance as its residual, so the flux through every face is computed from the pressures on its two sides without interpolation. The finite-volume balances then use these face velocities directly: for example the convective term of the component balance in cell $n$ is the difference of the face fluxes

$$
F_{n+1/2} - F_{n-1/2},\qquad F_{n+1/2} = \frac{T_n}{P_n}\thinspace \frac{y_{i,n+1/2}\thinspace P_{n+1/2}\thinspace u_{n+1/2}}{T_{n+1/2}}
$$

with $P$ and $T$ at a face taken as the mean of the two cells. This staggering avoids the checkerboard pressure modes a collocated grid can develop, and makes the scheme conservative: what leaves one cell through a face enters its neighbour.

**Dispersion and conduction** fluxes are central differences. At a face between cells of different width or different dispersion coefficient, the face coefficient is the cell-width-weighted harmonic mean of the two cell values, so a sharp change of packing is handled consistently.

**Boundary faces (Danckwerts).** No dispersive or conductive flux crosses the two boundary faces of the bed: the flux through a boundary face is the convective one alone, evaluated with the boundary value. At an **inlet** the component flux is therefore exactly $u\thinspace C\thinspace y_{i,in}$ and the enthalpy flux $\rho\thinspace c_p\thinspace u\thinspace T_{in}$, the amounts the stream logs, the [performance metrics](performance-metrics.md) and the [cyclic-steady-state](cyclic-steady-state.md) mass balance record; this is the finite-volume form of the Danckwerts inlet condition $u\thinspace C\thinspace y_{in} = u\thinspace C\thinspace y(0^+) - \varepsilon_b D_{ax} C\thinspace \partial y/\partial z\rvert_{0^+}$, and the first cell's composition can differ from the feed's. At an **outlet** or a **closed** end the ghost cell mirrors the first cell (zero gradient), so there was no flux to drop. The species ghost cell holds the boundary composition itself; the temperature ghost cell is set so that the mean of the ghost and the first cell, the face temperature, is the boundary temperature, and the energy balance takes that face temperature as the upstream value of an inlet face.

### Momentum: Darcy and Ergun as face residuals

For a face between cells $n$ and $n+1$ the two half-cells act in series, each with its own viscous resistance $\mu\thinspace \Delta z/2K$ and, for Ergun, inertial resistance $\beta\rho\thinspace \Delta z/2$. Darcy faces (Kozeny–Carman or Blake–Kozeny, no inertial term) use the explicit residual $u - \Delta P / A = 0$. **Ergun faces are solved implicitly**: instead of imposing the root of the quadratic $\Delta P = A u + B\lvert u\rvert u$, which behaves like a square root near zero flow, each face carries the residual

$$
r_u = u - \frac{\Delta P - B\thinspace \lvert u\rvert\thinspace u}{A} = 0
$$

which is smooth through flow reversal, linear in the pressures and has a slope of at least 1 in $u$. Its solution is exactly the Ergun root. The boundary face of a **pressure** or **closed** boundary condition, which fixes that face's pressure drop, instead imposes the explicit Ergun root of the boundary pressure drop: when a step switch makes that face's velocity jump, the implicit residual's slope in $u$ changes by an order of magnitude and IDA's corrector (whose Jacobian predates the jump) fails to converge, while the explicit root is linear in $u$ and depends only on the first-cell pressure, which is continuous across the switch. Flow boundaries, which fix the face velocity, keep the implicit residual. With this split Ergun decks cross step switches as cleanly as Darcy ones. The [Momentum Balance](momentum-balance.md) page gives the details, including void regions and flow boundaries.

### Advection schemes

The mole fraction at a face, $y_{i,n+1/2}$, has to be reconstructed from the cell values. Three schemes are implemented; they change only this face value, never the conserved cell values, so total species mass is conserved whichever is used:

- **Upwind** (the default): the face takes the value of the upstream cell. First order, monotone and robust; it smears sharp fronts over a few cells, so use enough cells.
- **WENO** (`AdvectionScheme::WENO`): the third-order WENO scheme of Jiang and Shu (two stencils, $k = 2$) used by Haghpanah et al. (2013) and Yancy-Caballero et al. (2020). For flow in the $+z$ direction it blends the central interpolation $\tfrac12(y_n+y_{n+1})$ and the upwind extrapolation $\tfrac32 y_n - \tfrac12 y_{n-1}$ with linear weights 2/3 and 1/3, $y_{n+1/2} = \dfrac{\alpha_0\thinspace\tfrac12(y_n+y_{n+1}) + \alpha_1\thinspace(\tfrac32 y_n - \tfrac12 y_{n-1})}{\alpha_0+\alpha_1}$, $\alpha_0 = \dfrac{2/3}{(\beta_0+\varepsilon)^2}$, $\alpha_1 = \dfrac{1/3}{(\beta_1+\varepsilon)^2}$, with smoothness indicators $\beta_0 = (y_{n+1}-y_n)^2$, $\beta_1 = (y_n-y_{n-1})^2$ and $\varepsilon = 10^{-12}$ (the reference codes write $(\Delta y + 10^{-10})^4$, which is the same away from flat profiles). With the linear weights this is the third-order upwind-biased reconstruction $-\tfrac16 y_{n-1} + \tfrac56 y_n + \tfrac13 y_{n+1}$; near a front the weight shifts to the smoother stencil. On a non-uniform grid (layers of different cell size) both candidates are the linear interpolation or extrapolation in $z$, and $\beta$ is the slope times the cell width.
- **van Leer** (`AdvectionScheme::VanLeer`): a limited second-order reconstruction $y_{n+1/2} = y_n + \tfrac12\thinspace \phi(r)\thinspace (y_{n+1}-y_n)$ with $r = (y_n - y_{n-1})/(y_{n+1}-y_n)$ and the van Leer limiter $\phi(r) = (r + \lvert r\rvert)/(1 + \lvert r\rvert)$. The ratio is evaluated as $r = ab/(b^2 + 10^{-24})$ with $a = y_n - y_{n-1}$, $b = y_{n+1} - y_n$, which is the same for $y$ and $1-y$ and goes smoothly to zero on a flat profile.

**Boundaries.** The ghost cell holds a value *at* the bed-end face, not at a cell centre: the inlet composition (Dirichlet), or a copy of the end cell at an outlet. For the higher-order schemes the two bed-end faces therefore take the upwind value, exactly as Upwind does, so the convective flux through a bed end is the stream's own composition times its flow and equals the stream integral the KPIs and the cycle mass balance are built from. In WENO the first interior face next to an inlet extrapolates from that face value, half a cell away: $2y_1 - y_0$ with $\beta_1 = \left(2(y_1-y_0)\right)^2$, as in the reference implementations; van Leer doubles the upwind slope in $r$ there in the same way. A void region is a well-mixed tank, so both schemes treat it in the same way (its faces are upwind and its value sits on the face of the next cell). The reconstructed face value is clamped to [0, 1] so a limiter overshoot at a steep front cannot inject a negative or super-unity composition. Every face is single-valued, so each species is conserved whatever the scheme.

**Sum of the face mole fractions.** Every species and the total have their own balance, so nothing in the equations pulls the cell sum $\sum_i y_i$ back to one; the species fluxes must add up to the flux the total balance carries. Upwind faces sum to the upwind cell's sum, so a departure of the cell sum (at the level of the solver tolerance) is simply carried along and out of the bed. The nonlinear reconstructions, built component by component, do not keep that sum, and worse, linearised about a steep front they weight the downwind cell, so a departure of the cell sum grows instead of being carried: before this was handled, van Leer drove the cell sum to 90 in a CO2/N2 breakthrough on 13X. The interior faces of WENO and van Leer are therefore rescaled so that their mole fractions sum to the upwind cell's sum (blended through zero flow like the upwind value). The sum then behaves exactly as under upwinding while each component keeps its high-order profile; each face stays single-valued, so each species is still conserved. This couples every mole fraction across the stencil, which widens the Jacobian band (below).

**Smoothed upwinding at zero flow.** A hard switch of the upwind cell when $u$ changes sign makes the residual discontinuous, and Newton iterations chatter when a face sits at zero flow (a closed bed end, a bed that has just reached the boundary pressure). The face value therefore blends the two candidates with the weight $w = \tfrac12\left(1 + \tanh(u/u_s)\right)$, $u_s = 10^{-6}$ m/s. For $\lvert u\rvert > 20\thinspace u_s$ this is exactly the hard upwind choice. Flow boundaries blend their inlet composition and temperature in the same way as the boundary flow passes through zero.

The scheme is a solver setting stored in the deck (`"Advection_Scheme"` in `Solver_Settings`: 0 = Upwind, 1 = WENO, 2 = van Leer) and in the C++ API (`IDASettings::advection_scheme_`). It is also the parameter `settings.advection_scheme` (an integer) in the Python parameter registry.

### Time integration: SUNDIALS IDA

After discretisation the model is a differential-algebraic system $F(t, \mathbf{y}, \dot{\mathbf{y}}) = 0$, integrated with **IDA** from SUNDIALS, a variable-order, variable-step **backward differentiation formula (BDF)** solver. At each step IDA predicts the new state from a polynomial through the previous points, then corrects it with a Newton iteration on the implicit BDF equations.

**Settings** (Solver → Solver Settings; `IDASettings`):

| Setting | Default | Meaning |
|---|---|---|
| Max BDF order | 5 | highest order IDA may use (1–5) |
| Relative tolerance | 10⁻⁶ | one scalar for every unknown |
| Initial step | 10⁻¹⁰ s | first step attempted, and after every restart |
| Min step | 10⁻¹⁵ s | IDA fails rather than go below it |
| Max step | 0 = no limit | cap on the internal step size |
| Max nonlinear iterations / step | 15 | Newton iterations per step attempt |
| Nonlinear convergence coef | 0.1 | safety factor on the Newton convergence test |
| Max convergence failures | 200 | Newton failures allowed in one step |
| Max error-test failures | 40 | local-error failures allowed in one step |
| Max number of steps | 100 000 | IDA's limit between two output times |

**Absolute tolerances** are a vector, one value per kind of unknown, fixed in the code: $10^{-8}$ for the total concentration [mol m⁻³], the mole fractions and the loadings [mol kg⁻¹]; $10^{-6}$ for the pressure [Pa] and the face velocities [m s⁻¹]; $10^{-4}$ K for the three temperatures. IDA accepts a step when the weighted norm of its local error estimate, with weights $1/(\text{rtol}\thinspace \lvert y\rvert + \text{atol})$, is below 1.

**Differential and algebraic unknowns.** IDA is told which unknowns are algebraic (pressure, the face velocities and every ghost-cell value, which only appear through boundary conditions) and which are differential. This is what `IDACalcIC` needs.

**Consistent initial conditions.** The initial conditions you enter set the differential unknowns but rarely satisfy the algebraic equations (the boundary conditions, the momentum balance, the ideal gas law). Before the first step, `IDACalcIC` in `IDA_YA_YDP_INIT` mode solves for the algebraic unknowns and the time derivatives of the differential ones so that $F = 0$ holds, and the corrected state is what gets logged at $t = 0$. If it fails, the run stops with *Initial condition calculation failed*, which almost always means the initial conditions are far from what the first step imposes (see [Running the Solver](solver-settings.md#troubleshooting-convergence)).

**Steps of a cycle.** Each cycle step sets IDA's stop time to the end of the step and integrates to it one internal step at a time (`IDA_ONE_STEP`), so a step boundary is never stepped over. The new step's boundary conditions take effect when the step starts. Results are written at up to **200 points per step** (`max_log_points_per_step_`), evenly spaced in time and filled in from IDA's interpolating polynomial, plus the first and last point of every step; the boundary flows that feed the KPIs are logged at every accepted step. Only the last **Max cycles to keep** cycles (default 10) are retained for plotting.

**Recovery from failures.** If IDA fails on a recoverable error (corrector or error-test failures, a residual or linear-solver error, too much work), the step is retried up to four times. Each retry restarts IDA from the last accepted point with a fresh history, zeroed algebraic derivatives and a new `IDACalcIC`, which clears the large derivative estimates an abrupt boundary switch can leave in the history. The log reports each retry; after the fourth the run fails with the IDA error and the unknown with the largest local error.

### The linear solver and the Jacobian

Each Newton iteration solves a linear system with the iteration matrix $J = \partial F/\partial \mathbf{y} + c_j\thinspace \partial F/\partial \dot{\mathbf{y}}$. The unknowns are ordered **cell by cell** (all unknowns of cell 0, then cell 1, and so on), so every equation couples only to its own cell and its neighbours and $J$ is **banded**. Cross-variable couplings reach one cell either side, so a half-bandwidth of twice the number of unknowns per cell covers the upwind scheme. The WENO and van Leer stencils reach two cells, and their interior face mole fractions are rescaled over the components, so a mole fraction there reads every mole fraction two cells away; for these schemes the half-bandwidth is twice the number of unknowns per cell plus the number of components minus one. SUNDIALS' **banded direct LU solver** (`SUNLinSol_Band`) factorises it.

The Jacobian itself is approximated by difference quotients: perturb a group of unknowns, re-evaluate the residuals, and read off the columns. How the groups are chosen is the **Jacobian mode** (`JacobianMode`, a solver implementation detail that is not saved in decks; the `bench_solver` tool selects it with `--jacobian ida|band|coloured`):

| Mode | What it does |
|---|---|
| **IdaInternal** | IDA's own banded difference-quotient Jacobian: columns a full bandwidth apart are perturbed together, so each Jacobian costs one residual evaluation per column of the band (upper + lower half-bandwidth + 1). |
| **BandDQ** | Skarstrom's own replica of IDA's banded algorithm, with exactly the same increments and column groups, so every entry is bit-identical to IdaInternal. Used as the reference for the coloured mode. |
| **Coloured** (default) | Difference quotients over a **learned sparsity pattern**. Columns whose rows do not overlap are perturbed together, which needs far fewer residual evaluations per Jacobian than the band stride. |
| **ColouredCheck** | Coloured, plus a full band evaluation at every Jacobian for comparison; a diagnostic that reports any mismatch. |

How the Coloured mode stays safe:

- The pattern is learned from band evaluations: at three randomised "generic" states near the current one at the start of each step of the first cycle (so couplings that vanish only at special states such as $u = 0$ or $y = 0$ still show up), and at the first ten real evaluations of every step in the first cycle, since each step has its own boundary conditions.
- Interior columns share a translation-invariant pattern keyed by variable and cell offset; the three cells at each end keep their own exact patterns, because boundary conditions overwrite those residuals.
- Every 50th Jacobian is a full band evaluation that audits the pattern; any coupling found outside it is added and the columns are recoloured, so a missed dependency corrects itself.
- The Jacobian requested after a Newton convergence failure is always a full band evaluation, so in rare regimes (a boundary whose flow reverses, for example) the solver is never less robust than IDA's own Jacobian.

The benchmark output of `bench_solver` reports the number of colours, audits and pattern misses for each deck.

### Cyclic steady state: detection and acceleration

After the last step of every cycle the solver applies the deck's [cyclic steady state](cyclic-steady-state.md) rule (`CyclicSteadyState::CheckCSS`); if it holds, the run stops and the status becomes *Cyclic steady state reached at cycle n*. Each method has a single-cycle test (below), and the rule holds when that test has passed on **Consecutive cycles** cycles in a row (default 3) and at least **Minimum cycles** cycles have run. The check reads only the equation log, so the logs keep at least the cycles it needs.

| Method | Single-cycle test |
|---|---|
| None | never; the run does every cycle |
| Molar Profile | $\lVert \mathbf{C}^{\text{curr}} - \mathbf{C}^{\text{prev}}\rVert_2 \le \varepsilon_{abs} + \varepsilon_{rel}\lVert \mathbf{C}^{\text{curr}}\rVert_2$, where $\mathbf{C}$ is the total-concentration profile over all cells at the end of the **first step** of the current and the previous cycle |
| Mass Balance (Total) | $\lvert n_{in} - n_{out}\rvert / n_{in} \le \varepsilon_{mb}$ over the whole cycle, summing every face of every step by the sign of its net flow, **and** $\Delta_{\text{state}} \le \varepsilon_{\text{state}}$ |
| Mass Balance (Per-Component) | the balance for each component that is fed, taking the worst, **and** $\Delta_{\text{state}} \le \varepsilon_{\text{state}}$ |

$\Delta_{\text{state}}$ (`CyclicSteadyState::EndOfCycleStateChange`) compares the last logged row of the cycle's last step with the previous cycle's: the largest of $\operatorname{RMS}(x^{c} - x^{c-1}) / s$ over the adsorbed loading of each component ($s$ the RMS total loading), the mole fractions ($s = 1$) and the fluid, solid and wall temperatures ($s$ their RMS). The boundary flux in the balance is reconstructed from the logged face values and integrated over the logged points, so the balance carries a small systematic error of its own; the state test does not, and it does not cross zero the way the balance does on its way to CSS. The [Cyclic Steady State](cyclic-steady-state.md#state-test) page gives the reasoning and the measured stopping cycles.

**Aitken (Graves–Morris) acceleration.** The approach to CSS is often a slow geometric decay: each cycle closes a fixed fraction of the remaining gap. With **CSS acceleration method** set to *Aitken (Graves-Morris)*, the solver extrapolates the end-of-cycle state towards the fixed point (`AccelerateCyclicSteadyState::Accelerate`):

1. From the third cycle on, it collects three successive natural end-of-cycle states $x_0, x_1, x_2$.
2. It computes one scalar coefficient shared by all fields,

$$
\mu = \frac{\sum_f w_f\thinspace r_{2,f}\cdot d_f}{\sum_f w_f\thinspace d_f\cdot d_f},\qquad r_1 = x_1 - x_0,\thickspace r_2 = x_2 - x_1,\thickspace d = r_2 - r_1,\thickspace w_f = \frac{1}{\overline{x_f^2}}
$$

   where the weights stop fields of large magnitude (temperatures) from dominating fields of small magnitude (mole fractions), and moves every field to $x^{*} = x_2 - \mu\thinspace r_2$. Sharing one coefficient keeps the moved fields consistent with one another.
3. Only the slow differential fields are moved: the adsorbed loadings, the solid and wall temperatures and the mole fractions. The total concentration and gas temperature are left alone (they fix the pressure through the gas law and relax within the next step), and the algebraic unknowns are recomputed.
4. It skips the jump if the steps are not shrinking, if there is no curvature to fit, or if $\lvert\mu\rvert > 5$ (a decay rate above about 0.83, too close to 1 to extrapolate safely). Negative values are clipped to zero and the mole fractions renormalised to sum to one.
5. The extrapolated state is first integrated for a short probe (the smaller of 2 s and a tenth of the first step) in a throwaway copy of the integrator. Only if that succeeds is it committed to the live run, with a fresh `IDACalcIC`; otherwise it is discarded and the log says *CSS acceleration skipped*.
6. The three samples are then discarded, so the next extrapolation uses three fresh cycles. Jumps are therefore at least three cycles apart and the bed relaxes between them.

The gain depends on the deck. For the [Perez 2019 VPSA tutorial](perez-2019-vpsa.md), acceleration moved CSS from cycle 190 to cycle 182, with the same purity and recovery.

### The chunked solver for cooperative front ends

The C++ API (and the Python package) can run a whole cycle loop in one call (`CoupledSolver::IntegrateCycle`). A single-threaded host such as a browser worker cannot do that: a long blocking call could never receive a *Cancel* message. For such hosts the engine provides `CoupledSolverChuncked`, which performs the same steps, retries, CSS checks and acceleration with the same integrator code, but in slices: each call to `Integrate(max_time_ms)` advances the current step until the wall-clock budget is used or the step ends, then returns. A WebAssembly front end can call it with a short budget (e.g. 50 ms) in a loop, yielding between slices so that progress, log messages and cancel requests get through. Results are the same as a native run up to floating-point rounding differences between compilers.
