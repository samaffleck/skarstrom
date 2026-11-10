# Boundary Conditions

Every [step](cycles.md) specifies what happens at each end of the column through two boundary conditions — one at the **Top** ($z = L$) and one at the **Bottom** ($z = 0$). The boundary condition determines whether gas enters, leaves, or is sealed at that face, and with what driving force.

Every step has a **Top** and a **Bottom** boundary. The first choice for each is its **boundary condition type**.

---

## Boundary condition types

### Pressure

Imposes a specified **pressure** at the face. Use it for steps controlled by a pressure level — pressurisation, blowdown, or a constant-pressure product outlet.

- As an **inlet**, you also specify the inlet **Temperature** and **Composition**; gas enters at the imposed pressure.
- As an **outlet**, the temperature and composition are taken from the bed (zero-gradient); gas leaves at the imposed pressure. If the bed falls below the imposed pressure, gas flows back in through the outlet with the bed-end temperature and composition, and the run warns at the end ([Backflow through an outlet](#backflow-through-an-outlet)).

The pressure follows a **Pressure Profile** over the step (see [Profiles](#time-varying-profiles) below), so you can ramp the pressure up or down during the step. A **pump** (compressor or vacuum pump) can be attached to compute compression power and heating — see [Pumps](#pumps-compressors-and-vacuum-pumps).

### Flow

Imposes a **molar flow** at the face, specified through one or more **streams**. Use it for feed steps and for flow-controlled purge or product steps.

A Flow boundary contains a list of streams (add them with **Add Stream**); the total flow is the sum of the streams. Each stream has:

- a **Boundary Device Type** — **None** (flow specified directly), **Valve** (pressure-driven), or **Pump** (power-driven);
- for a *None* device, a **Flow Type** and a flow profile;
- as an inlet, its own **Temperature** and **Composition**.

Multiple streams on one face let you, for example, mix a fresh feed with a recycle stream.

### Closed

A sealed, no-flow end. Nothing enters or leaves; pressure, temperature and composition have zero gradient at the face. Use it for the closed end of a blowdown or equalisation step. A Closed boundary has no further settings.

---

## Flow direction

Each boundary has a **Flow direction**:

- **Inlet** — gas flows into the reactor. You specify what enters (temperature, composition, and — for a Flow boundary — the flow rate).
- **Outlet** — gas flows out of the reactor. The leaving gas takes the bed's local temperature and composition.

A step must have at least one non-inlet end (an outlet or a closed end); the two ends cannot both be inlets.

## Specifying inlet flow

For a Flow boundary with a *None* device, the **Flow Type** sets how the rate is expressed. All four are equivalent and are converted internally using the boundary-face area and the ideal-gas density:

| Flow Type | Quantity |
|-----------|----------|
| **Velocity** | Superficial velocity [m s⁻¹] |
| **Volumetric Flow Rate** | Volumetric flow [m³ s⁻¹] |
| **Molar Flow Rate** | Molar flow [mol s⁻¹] |
| **Mass Flow Rate** | Mass flow [kg s⁻¹] |

## Inlet composition and temperature

At an inlet you set the **Temperature** and the **Composition** (a mole fraction for each component). The mole fractions of an inlet must sum to 1. Outlets ignore these fields because the leaving gas is whatever is in the bed. When the flow through an outlet reverses, the gas that flows back into the bed carries the bed-end composition and temperature (zero gradient). An outlet's composition is used in one case: when its backflow feeds another stream that leaves through the same bed end (see [Backflow through an outlet](#backflow-through-an-outlet)).

### Inlet temperature profile

The inlet temperature of a stream (or of a Pressure inlet) can change during the step, for example to ramp a TSA feed or regeneration gas from cold to hot without splitting the step. Pick a **Temperature Profile**:

- **Constant** — the default: the single **Temperature** above, held for the whole step. Decks saved before temperature profiles existed load as Constant and are unchanged.
- **Linear** — from the **Start Temperature** $T_0$ to the **Final Temperature** $T_1$ over the **Ramp Time** $t_r$, then held at $T_1$ (it never extrapolates past the end of the ramp):

$$
T_{\text{in}}(t) = T_0 + \left(T_1 - T_0\right)\frac{\min(t,\thinspace t_r)}{t_r}
$$

- **Exponential** — a first-order approach from $T_0$ to $T_1$ whose remaining change has shrunk by the **Exponential Factor** $F$ after the ramp time $t_r$ (the default $F = 100$ has covered 99 % of the change by then):

$$
T_{\text{in}}(t) = T_1 + \left(T_0 - T_1\right)F^{-t/t_r} = T_1 + \left(T_0 - T_1\right)e^{-\ln(F)\thinspace t/t_r}
$$

- **Log-linear** — a geometric ramp from $T_0$ to $T_1$ over the ramp time, then held at $T_1$ (see [time-varying profiles](#time-varying-profiles)).

$t$ is the time since the start of the step, so the ramp restarts every time the step runs. A ramp time of zero is an immediate step to $T_1$, and a ramp longer than the step simply ends part-way. Switching from Constant starts both temperatures at the current value with a ramp time equal to the step time.

The profile replaces the typed-in temperature only. When the inlet is fed from a [buffer tank](#drawing-gas-from-another-step-buffer-tank) or a connected stream, the gas keeps that source's temperature as before, and the profile acts as the first-cycle assumption. In an [isothermal](energy-balance.md#isothermal-mode) run inlet temperatures have no effect.

### Heating and cooling the column

A step can also heat or cool the column through a [heating medium](wall-energy-balance.md#heating-medium-jacket): a fluid at a set temperature, constant or ramped with the same profiles as an inlet temperature.

- **Jacket.** The medium flows over the wall's outer face with heat transfer coefficient $h_m$, and replaces the ambient exchange while the step runs.
- **Internal heat exchanger** (optional). The same medium also runs through tubes or plates in the packing and heats the bed directly: $a_x h_x (T_m - T)$ per bed volume, with $a_x$ the exchanger's area per bed volume (see [internal heat exchanger](wall-energy-balance.md#internal-heat-exchanger)).
- Steps without a medium keep the ambient exchange.
- The heat the medium puts in and takes out is logged per step and counted in the thermal-energy and cooling-duty [KPIs](performance-metrics.md).
- The values are parameters of the step: `step[ID].medium.enabled`, `.temperature`, `.heat_transfer_coefficient`, `.exchanger_area_per_volume` and `.exchanger_heat_transfer_coefficient`.
- This is different from an [energy source](energy-balance.md#energy-source-heating-cooling), which puts a fixed power into a layer whatever its temperature. A medium's heat flow depends on how far the column is from $T_m$, so it stops once the column reaches $T_m$.

---

## Time-varying profiles

Flow rates, pressures and [inlet temperatures](#inlet-temperature-profile) can all vary during a step. The **Profile Type** (for flow) or **Pressure-Time Profile** (for pressure) offers:

- **Constant** — a single value held for the whole step.
- **Linear** — ramps linearly from a **Start Value** to a **Final Value** over a **Ramp Time**.
- **Exponential** — approaches the **Final Value** $v_1$ from the **Start Value** $v_0$ with a first-order response: after the ramp time $t_r$ the change left is $1/F$ of the initial one, where $F$ is the **Exponential Factor** (default 100, i.e. 99 % of the change covered by $t_r$; it must be greater than 1). The curve keeps approaching $v_1$ after $t_r$:

$$
v(t) = v_1 + \left(v_0 - v_1\right)F^{-t/t_r}
$$

- **Log-linear** — a geometric ramp from the **Start Value** $v_0$ to the **Final Value** $v_1$ over the ramp time $t_r$, then held at $v_1$. The value is multiplied by the same factor every second, so $\ln v$ rises linearly:

$$
v(t) = v_0\left(\frac{v_1}{v_0}\right)^{\min(t,\thinspace t_r)/t_r},
\qquad
\frac{d\ln v}{dt} = \frac{\ln(v_1/v_0)}{t_r}
$$

Ramps are the natural way to model gradual pressurisation or blowdown.

**Pressure ramps** (a Pressure boundary or a valve) start from the bed pressure at the start of the step, so they have no Start Value. By default (**Use Step Time** on) the ramp lasts the whole step. Turn it off to give the ramp its own **Ramp Time**, for example to pressurise in the first 10 s of a 30 s step and hold the final pressure for the rest; the ramp time must be no longer than the step (the run is refused otherwise). Decks saved before this option existed ramp over the step time, as before.

### Choosing a pressurisation ramp

In upward flow the pressure gradient at the bed inlet scales roughly with $d\ln P/dt$, the boundary pressure's rate of rise relative to the pressure itself (in the viscous regime $-\partial P/\partial z \propto \mu R T N / (K P)$, and the molar flux $N$ needed to fill the bed scales with $dP/dt$). A **linear** ramp from vacuum has its largest relative rise in the first seconds: from 3 kPa to 4.6 bar in 20 s it starts at 7.6 s⁻¹ and ends at 0.05 s⁻¹, which is where a spike in the approach to fluidisation appears. An **exponential** ramp front-loads the rise even more (35 s⁻¹ at the start). A **log-linear** ramp holds $d\ln P/dt$ at $\ln(P_1/P_0)/t_r$ (0.25 s⁻¹ here) for the whole step, so the gradient is spread evenly. The trade-off is that the bed stays at low pressure for longer (only $\sqrt{P_0 P_1}$ halfway through) and the pressure rises fastest in absolute terms at the end of the step, so check the approach to fluidisation there too, and the cycle's performance.

A geometric ramp cannot start from or pass through zero, so when the start or final value is zero or negative the log-linear profile ramps linearly instead (for example a flow ramp that starts from zero).

- **Double exponential** — a two-term decay towards the **Final Value**: a share $w$ of the initial excess decays at **Rate 1** and the rest at **Rate 2** (both in s⁻¹, entered as fitted, not as a 99 % time):

$$
v(t) = v_1 + (v_0 - v_1)\left[w\thinspace e^{-r_1 t} + (1 - w)\thinspace e^{-r_2 t}\right]
$$

This is the form Estupiñan Perez et al. (2019) fitted to the evacuation pressure of their VSA rig (a fast term for the column gas, a slow one for the pump and piping, $r_1 = 0.019$ and $r_2 = 0.101$ s⁻¹), where a single exponential did not describe the measured curve. For a pressure the start is the bed pressure at the start of the step; a rate of zero freezes its term at the start value.

> All numeric fields accept [expressions and units](global-variables.md), e.g. `2 [bar]` or a reference to a global variable.

### Starting a flow inlet into a bed at low pressure

A constant flow inlet switched on at the start of a step into a bed at low pressure (a feed or pressurisation step straight after evacuation or blowdown) is a step change at the boundary: the first cell sees the full molar flux at once while the gas there is thin, so the pressure and velocity at the inlet jump. The integrator often fails its first attempts at such a step, and the engine recovers by retrying (below). Ramping the flow from zero over about 1 s (a **Linear** profile with Start Value 0 and Ramp Time 1 s, or `sk.linear(v, start=0, time=1)` in Python) removes the jump and the failed attempts, and changes the amount fed by only half a second's worth of flow. A pressure boundary is smoothed the same way with a short ramp instead of a constant value.

### Recovered steps in the log

When the integrator (IDA) fails a step (repeated error-test or convergence failures, a linear-solver failure or too much work), the engine retries it up to four times, from a smaller step or, at a step-size floor, from a consistent restart at the last accepted point. A failure the retry gets past is a **recovered step**:

- it does not stop the run or change the results beyond the integration tolerance;
- the log shows it at level *medium* as `Step "<name>" recovered step: IDA failure (code <n>) at t=... ; retry k/4 ...`, followed by IDA's own message;
- it counts in the run statistics, `recovered_steps` (*New in 0.9.* `result.stats["recovered_steps"]` in Python, `stats.recovered_steps` in the CLI `--json` summary), and in the [run-health](run-health.md) `STEP_RETRIES` measure, which reads the same count.

Only a step whose retries run out is a failure: it is logged at level *low* with IDA's last message and the run stops. IDA's messages never go to stderr; at level *high* the log shows every one of them (`IDA: ...`).

---

## Valves

Set a stream's **Boundary Device Type** to **Valve** to make its flow **pressure-driven** rather than prescribed. You enter a **Valve Flow Coefficient (Cv)** and a pressure profile; the flow through the valve is driven by the pressure difference between the valve's set pressure and the bed. The valve is linear in the pressure difference (`Valve::GetMolarFlowRate` in `core/src/objects/valve.cpp`):

$$
\dot n = C_v\thinspace \frac{P_{\text{up}} - P_{\text{down}}}{10^{5}\ \text{Pa}}
$$

so $C_v$ is in mol s⁻¹ bar⁻¹: a valve with $C_v = 10$ passes 10 mol/s for a 1 bar difference. Gas flows from the higher to the lower pressure, and the flow reverses if the pressures cross.

The bed side of the law is the pressure at the bed end (the boundary face, the mean of the ghost cell and the end cell). The solver setting `settings.valve_bed_pressure` (`"Valve_Bed_Pressure"` in `Solver_Settings`; `IDASettings::valve_bed_pressure_`) can instead read the pressure of the end cell, the last cell at that end of the bed, as models that treat the column as a series of well-mixed cells do (toPSAil attaches its valves to the last cell). A connected valve then reads its partner's end cell as well. The two readings differ by the pressure drop over half a cell and agree as the grid is refined; use `end_cell` only to compare with such a model.

Valves are the mechanism behind [pressure equalisation](pressure-equalization.md): put a valve on each bed and connect them, and gas flows from the high-pressure bed to the low-pressure one at a rate set by Cv.

### Check valve

!!! note "New in 0.9"
    The check valve (`endpoint[ID].valve.check_valve`) and the outlet backflow warning (`settings.warn_outlet_backflow`).

Turn on a valve stream's **Check Valve** (`endpoint[ID].valve.check_valve`; Python `sk.ValveStream(..., check_valve=True)` or `endpoint.valve.check_valve = True`) to stop the flow from reversing. An outlet then passes gas only out of the bed, and an inlet only into it. While the valve is shut, the stream carries no flow, so the end acts as a closed end for that stream. The closure is smooth: the valve law sees the pressure difference $\Delta P$ in the stream's direction replaced by

$$
s(\Delta P) =
\begin{cases}
0 & \Delta P \le 0 \\
\Delta P\left(2 x - x^2\right),\ x = \Delta P / \delta & 0 < \Delta P < \delta \\
\Delta P & \Delta P \ge \delta
\end{cases}
$$

with $\delta = 10$ Pa (`Valve::kCheckValveBand`). Beyond 10 Pa this is the plain valve law. A hard cut-off would put a kink in the boundary residual where the flow decays through zero, and the solver's Newton iteration would then fail to converge there.

The check valve is off by default, so decks keep their results. It applies only to an **unconnected** stream. A [connected](pressure-equalization.md) valve ignores it, because flow both ways through a connection (pressure equalisation, light-product pressurisation) is expected.

### Opening ramp

!!! note "New in 0.9"
    The valve opening ramp (`endpoint[ID].valve.opening_start`, `endpoint[ID].valve.opening_ramp_time`).

A real valve takes time to open. A valve stream that opens at once onto a bed at a very different pressure passes its largest flow in the first instant, which can lift the bed (a large [approach to fluidisation](wall-design.md)) at the end it feeds. The **opening ramp** lets the valve open from a part opening to fully open over a set time from the start of the step. The flow of the valve law is multiplied by the opening fraction

$$
f(t) =
egin{cases}
f_0 & t \le 0 \
f_0 + (1 - f_0)\left(3 x^2 - 2 x^3
ight),\ x = t / t_r & 0 < t < t_r \
1 & t \ge t_r
\end{cases}
$$

with $f_0$ the **Opening Start** (`endpoint[ID].valve.opening_start`, in (0, 1]) and $t_r$ the **Opening Ramp Time** (`endpoint[ID].valve.opening_ramp_time`, s). The cubic has a zero slope at both ends of the ramp, so the flow has no kink there (a kink would stall the solver's Newton iteration, as for the check valve). The effective coefficient is $f(t)\thinspace C_v$: the ramp starts at $f_0 C_v$ and reaches $C_v$ at $t_r$.

The defaults, $f_0 = 1$ and $t_r = 0$, give a valve that is fully open for the whole step, so decks saved before this option existed keep their results. The ramp applies to connected valves as well. For a valve between two beds (pressure equalisation, light-product pressurisation), give both ends the same ramp: each end is one side of the same valve.

In Python:

```python
import skarstrom as sk

valve = sk.ValveStream("Feed header", pressure=1e5, cv=0.006, opening_start=0.05, opening_ramp_time=20)
print(valve.opening_start, valve.opening_ramp_time)
```

On a loaded deck, set `endpoint.valve.opening_start` and `endpoint.valve.opening_ramp_time`, or the registry paths above. Both accept [expressions](global-variables.md), for example a global variable shared by the two ends of a connected valve.

### Backflow through an outlet

An outlet is zero-gradient in composition and temperature, whatever the direction of the flow. When a valve outlet's bed falls below the downstream pressure, or a pressure outlet's set pressure rises above the bed, gas flows back into the bed through the outlet. This gas has the bed-end composition and temperature, not the outlet's `y` and temperature, and the [KPIs](performance-metrics.md) do not count it as feed.

If the same step also draws gas from that end, for example a light-product draw to another step, the backflow can feed that draw instead of the bed. The streams that leave through that end then carry a mix: the bed's net outflow at the bed-end composition, plus the back-fed gas at the back-feeding outlet's own composition (its `y`, normalised). So every component is conserved at the bed end, and the draw does not carry a copy of the bed's gas that the bed never gave. An outlet whose `y` is all zero has no composition of its own, and its back-fed gas takes the bed-end composition. The light-product decks give their vents the feed's inert (N2). The mix is on the composition only; the temperature of the leaving streams is the bed-end temperature.

At the end of a run, the solver warns when gas flows back in through an outlet that has no connection to another step (**Outlet backflow**, `settings.warn_outlet_backflow`). The check is on the last logged cycle. It warns when the gross moles in through one such outlet in one step exceed 0.1 % of the cycle feed (`kOutletBackflowFraction`; the feed is the gross moles in through every unconnected inlet). The warning names the step, the endpoint, the moles and the fraction of the feed. For a valve, it suggests the check valve. Connected endpoints are not checked.

## Pumps: compressors and vacuum pumps

A **Pump** device (or a pump attached to a Pressure boundary) models mechanical compression:

- an **inlet** uses a **compressor** (with a **Suction Pressure**);
- an **outlet** uses a **vacuum pump** (with a **Discharge Pressure**).

The pump settings include an **Efficiency Model** (constant **Isentropic Efficiency**, a vacuum-pump curve, or a polynomial), a **Polytropic Constant**, and an **Energy Model** (**Isothermal** or **Adiabatic**) that determines the outlet-gas temperature rise from compression. The efficiency is evaluated at the pump's suction pressure $P_1$ (in Pa); the polynomial is the full fifth-degree one:

$$
\eta(P_1) = a + b\,P_1 + c\,P_1^{2} + d\,P_1^{3} + e\,P_1^{4} + f\,P_1^{5}
$$

There are two distinct uses:

1. **Power-controlled flow** — a *Pump* device on a Flow stream sets the molar flow from a specified power (constant or a polynomial in pressure). This affects the simulation.
2. **Power for post-processing** — ticking **Calculate pump power** on a Pressure boundary or a *None* stream computes compression power and heating for reporting, without changing the flow the boundary imposes. The recorded power feeds the [energy-consumption KPI](performance-metrics.md).

### Ultimate pressure of a vacuum pump

!!! note "New in 0.9"
    `endpoint[N].pump.ultimate_pressure`.

A real vacuum pump cannot pull below its **ultimate (blank-off) pressure** $p_{ult}$. A *Pump* device on an outlet (a vacuum pump stream) takes one as `endpoint[N].pump.ultimate_pressure` [Pa] (alias `step[<step>].<end>[<stream>].pump.ultimate_pressure`). The flow that the power gives is multiplied by the effective-speed factor of a rotary pump,

$$
f = 1 - \frac{p_{ult}}{P_1},
$$

with $P_1$ the bed-end (suction) pressure. Within 2 % of the ultimate pressure ($f < 0.02$) the factor is rounded to $2u^2/d - u^3/d^2$ ($u = 1 - p_{ult}/P_1$, $d = 0.02$), so it falls to zero with a zero slope: the residual stays smooth for the integrator. At and below $p_{ult}$ the pump delivers no gas, so the bed pressure approaches $p_{ult}$ from above and is held there. The shaft power is unchanged (the pump keeps running at blank-off), so the energy KPIs still count it. The default, 0, turns the factor off, and a deck saved before 0.9 loads with it off. An inlet (compressor) stream has no ultimate pressure.

```python
import skarstrom as sk
from skarstrom import ureg as u

sim = sk.examples.Perez2019_A.load()
evac = sim.steps["Evacuation"]
evac.bottom.set_streams("outlet", [sk.PumpStream("Vacuum pump", power=10,
                                                 pump=sk.Pump("vacuum_pump", efficiency_model="vacuum_pump",
                                                              discharge_pressure=1e5, ultimate_pressure=10 * u.kPa))])
sim.run(cycles=1)
end = next(r for r in sim.step_end_pressures() if r.step_name == "Evacuation")
print(end.bottom)   # >= 1e4 Pa
```

To bound the pressure an evacuation ends at in a study or an optimisation, use the [step-end pressure](performance-metrics.md#step-end-pressures) outputs.

---

## Drawing gas from another step (buffer tank)

An inlet can be fed from the gas *produced by another step* instead of from values you type in. Tick **Use buffer tank from another step** under **Inlet Source**, choose the source step and endpoint, and set the percentage of that step's produced moles to draw. This is the basis of averaged [light reflux and recycle](light-reflux-recycle.md).

---

Next, learn how to couple two beds together in [Connecting Beds](connecting-beds.md).
