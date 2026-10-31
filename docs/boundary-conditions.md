# Boundary Conditions

Every [step](cycles.md) specifies what happens at each end of the column through two boundary conditions — one at the **Top** ($z = L$) and one at the **Bottom** ($z = 0$). The boundary condition determines whether gas enters, leaves, or is sealed at that face, and with what driving force.

Every step has a **Top** and a **Bottom** boundary. The first choice for each is its **boundary condition type**.

---

## Boundary condition types

### Pressure

Imposes a specified **pressure** at the face. Use it for steps controlled by a pressure level — pressurisation, blowdown, or a constant-pressure product outlet.

- As an **inlet**, you also specify the inlet **Temperature** and **Composition**; gas enters at the imposed pressure.
- As an **outlet**, the temperature and composition are taken from the bed (zero-gradient); gas leaves at the imposed pressure.

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

At an inlet you set the **Temperature** and the **Composition** (a mole fraction for each component). The mole fractions of an inlet must sum to 1. Outlets ignore these fields because the leaving gas is whatever is in the bed.

### Inlet temperature profile

The inlet temperature of a stream (or of a Pressure inlet) can change during the step, for example to ramp a TSA feed or regeneration gas from cold to hot without splitting the step. Pick a **Temperature Profile**:

- **Constant** — the default: the single **Temperature** above, held for the whole step. Decks saved before temperature profiles existed load as Constant and are unchanged.
- **Linear** — from the **Start Temperature** $T_0$ to the **Final Temperature** $T_1$ over the **Ramp Time** $t_r$, then held at $T_1$ (it never extrapolates past the end of the ramp):

$$
T_{\text{in}}(t) = T_0 + \left(T_1 - T_0\right)\frac{\min(t,\thinspace t_r)}{t_r}
$$

- **Exponential** — a first-order approach from $T_0$ to $T_1$ that has covered 99 % of the change after the ramp time $t_r$:

$$
T_{\text{in}}(t) = T_1 + \left(T_0 - T_1\right)e^{-\ln(100)\thinspace t/t_r}
$$

- **Log-linear** — a geometric ramp from $T_0$ to $T_1$ over the ramp time, then held at $T_1$ (see [time-varying profiles](#time-varying-profiles)).

$t$ is the time since the start of the step, so the ramp restarts every time the step runs. A ramp time of zero is an immediate step to $T_1$, and a ramp longer than the step simply ends part-way. Switching from Constant starts both temperatures at the current value with a ramp time equal to the step time.

The profile replaces the typed-in temperature only. When the inlet is fed from a [buffer tank](#drawing-gas-from-another-step-buffer-tank) or a connected stream, the gas keeps that source's temperature as before, and the profile acts as the first-cycle assumption. In an [isothermal](energy-balance.md#isothermal-mode) run inlet temperatures have no effect.

---

## Time-varying profiles

Flow rates, pressures and [inlet temperatures](#inlet-temperature-profile) can all vary during a step. The **Profile Type** (for flow) or **Pressure-Time Profile** (for pressure) offers:

- **Constant** — a single value held for the whole step.
- **Linear** — ramps linearly from a **Start Value** to a **Final Value** over a **Ramp Time**.
- **Exponential** — approaches the **Final Value** from the **Start Value** with an exponential (first-order) response over the ramp time.
- **Log-linear** — a geometric ramp from the **Start Value** $v_0$ to the **Final Value** $v_1$ over the ramp time $t_r$, then held at $v_1$. The value is multiplied by the same factor every second, so $\ln v$ rises linearly:

$$
v(t) = v_0\left(\frac{v_1}{v_0}\right)^{\min(t,\thinspace t_r)/t_r},
\qquad
\frac{d\ln v}{dt} = \frac{\ln(v_1/v_0)}{t_r}
$$

Ramps are the natural way to model gradual pressurisation or blowdown. The ramp time is clamped to the step duration.

### Choosing a pressurisation ramp

In upward flow the pressure gradient at the bed inlet scales roughly with $d\ln P/dt$, the boundary pressure's rate of rise relative to the pressure itself (in the viscous regime $-\partial P/\partial z \propto \mu R T N / (K P)$, and the molar flux $N$ needed to fill the bed scales with $dP/dt$). A **linear** ramp from vacuum has its largest relative rise in the first seconds: from 3 kPa to 4.6 bar in 20 s it starts at 7.6 s⁻¹ and ends at 0.05 s⁻¹, which is where a spike in the approach to fluidisation appears. An **exponential** ramp front-loads the rise even more (35 s⁻¹ at the start). A **log-linear** ramp holds $d\ln P/dt$ at $\ln(P_1/P_0)/t_r$ (0.25 s⁻¹ here) for the whole step, so the gradient is spread evenly. The trade-off is that the bed stays at low pressure for longer (only $\sqrt{P_0 P_1}$ halfway through) and the pressure rises fastest in absolute terms at the end of the step, so check the approach to fluidisation there too, and the cycle's performance.

A geometric ramp cannot start from or pass through zero, so when the start or final value is zero or negative the log-linear profile ramps linearly instead (for example a flow ramp that starts from zero).

> All numeric fields accept [expressions and units](global-variables.md), e.g. `2 [bar]` or a reference to a global variable.

---

## Valves

Set a stream's **Boundary Device Type** to **Valve** to make its flow **pressure-driven** rather than prescribed. You enter a **Valve Flow Coefficient (Cv)** and a pressure profile; the flow through the valve is driven by the pressure difference between the valve's set pressure and the bed. The valve is linear in the pressure difference (`Valve::GetMolarFlowRate` in `core/src/objects/valve.cpp`):

$$
\dot n = C_v\thinspace \frac{P_{\text{up}} - P_{\text{down}}}{10^{5}\ \text{Pa}}
$$

so $C_v$ is in mol s⁻¹ bar⁻¹: a valve with $C_v = 10$ passes 10 mol/s for a 1 bar difference. Gas flows from the higher to the lower pressure, and the flow reverses if the pressures cross.

Valves are the mechanism behind [pressure equalisation](pressure-equalization.md): put a valve on each bed and connect them, and gas flows from the high-pressure bed to the low-pressure one at a rate set by Cv.

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

---

## Drawing gas from another step (buffer tank)

An inlet can be fed from the gas *produced by another step* instead of from values you type in. Tick **Use buffer tank from another step** under **Inlet Source**, choose the source step and endpoint, and set the percentage of that step's produced moles to draw. This is the basis of averaged [light reflux and recycle](light-reflux-recycle.md).

---

Next, learn how to couple two beds together in [Connecting Beds](connecting-beds.md).
