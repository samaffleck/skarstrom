# Matching Another Model

This page is for anyone reproducing a result from another simulator or from the model in a paper, and for coding agents doing the same. It lists the conventions Skarstrom uses where codes commonly differ, with the formula the engine evaluates and the parameter you set. Python names (`sim.adsorbents[0].solid_density`) and registry paths (`adsorbent[1].solid_density`) are given side by side; see [Python and the Command Line](python.md) for how to set them.

All values are SI unless a unit is given (`2 * u.min`, `"2 [min]"`).

## Densities and voidages

### Bulk density and `solid_density`

`solid_density` is the **particle (pellet) density** $\rho_s$: the mass of an adsorbent particle divided by its volume, pores included. It is **not** the skeletal (crystal) density and not the bed's bulk density. The bulk density, the mass of adsorbent per m³ of bed, is

$$
\rho_b = (1 - \epsilon_b)\thinspace \rho_s
$$

with $\epsilon_b$ the layer's inter-particle voidage, `layer[ID].bed_void_fraction`. $\rho_b$ multiplies every loading rate in the [mass](mass-balance.md) and [energy](energy-balance.md) balances, so the loading $q$ (mol kg⁻¹) is per kg of particle, binder included.

Converting from what a paper gives:

| The source gives | Set `solid_density` to |
|---|---|
| particle / pellet / apparent density $\rho_p$ | $\rho_p$ |
| bulk (packing) density $\rho_b$ | $\rho_b / (1 - \epsilon_b)$ |
| skeletal / crystal / true density $\rho_{cr}$ | $(1 - \epsilon_p)\thinspace \rho_{cr}$ |

Example: a column packed with 13X at a bulk density of 660 kg m⁻³ and $\epsilon_b = 0.41$:

```python
sim.layers[0].bed_void_fraction = 0.41
sim.adsorbents["13X"].solid_density = 660 / (1 - 0.41)    # 1118.6 kg/m3 of particle
```

A layer of several adsorbents (`composition` mass fractions $w_k$) uses the mixture particle density $\rho_s = 1 / \sum_k (w_k / \rho_{s,k})$, and each adsorbent's own bulk density is $\rho_{b,k} = (1-\epsilon_b)\thinspace \rho_s\thinspace w_k$.

The mass fractions of a layer must sum to 1 (within $10^{-6}$) and none may be negative; otherwise the deck checks refuse the deck (`CheckLayerComposition` in `porous_media.h`, which `sim.check()` and every run use). The fractions are never renormalised. A solid at mass fraction 0 is allowed but takes no part in the run, and the run-start warnings say so. A solid whose components all use the Inert isotherm is a valid member of a blend (an inert diluent or packing). A solid that is part of a blend (a layer that also holds another solid) cannot be removed from the deck until it has been taken out of that layer's composition. The blend edits of `deck_edit.h` (`AddSolidToLayer`, `SetLayerMassFraction`, `RemoveSolidFromLayer`) keep the fractions summing to 1 by rescaling the other solids in proportion.

### Total voidage

The gas held in the bed (the accumulation terms of the mass and energy balances) uses the **total voidage**

$$
\epsilon_t = \epsilon_b + (1 - \epsilon_b)\thinspace \epsilon_p
$$

with $\epsilon_p$ the adsorbent's `particle_void_fraction` (`adsorbent[ID].particle_void_fraction`; for a mixed layer the solid-volume-weighted average). The gas in the macropores is assumed to have the composition of the gas between the particles: there is no separate pore-gas balance.

Three voidages enter in different places:

| Parameter | Symbol | Used in |
|---|---|---|
| `layer[ID].bed_void_fraction` | $\epsilon_b$ | $\rho_b$, $\epsilon_t$, the axial dispersion term ($\epsilon_b D_{ax}$), the momentum balance |
| `adsorbent[ID].particle_void_fraction` | $\epsilon_p$ | $\epsilon_t$ only |
| `adsorbent[ID].kinetics[C].macropore_void_fraction` | $\epsilon_p$ in the rate | the macropore rate coefficients only (see [Kinetic Models](kinetic-models.md)) |

The last two are **separate parameters**; changing one does not change the other. Set both when the source model has one particle porosity used in both places.

**When the other code has no particle porosity** (its gas accumulation term uses $\epsilon_b$ alone, or it writes $\epsilon$ for both), set `particle_void_fraction` to zero, so that $\epsilon_t = \epsilon_b$, and keep the macropore voidage of the rate model:

```python
ads = sim.adsorbents["13X"]
ads.particle_void_fraction = 0.0                   # eps_t = eps_b
ads.kinetics["CO2"].macropore_void_fraction = 0.35  # still used by the macropore k
```

Check the other model's density at the same time: if its solid term is $(1-\epsilon)\rho_s$ with $\epsilon$ its only voidage, its $\rho_s$ is a particle density in the sense above.

**Where in the code:** `MixedAdsorbentLayer::GetBulkDensity`, `GetTotalVoidFraction`, `GetEffectiveSolidDensity` and `GetEffectiveParticleVoidFraction` (`core/include/adsim/objects/porous_media.h`, `core/src/objects/porous_media.cpp`).

## Heat of adsorption sign

`adsorbent[ID].heat_of_adsorption[C]` (Python `sim.adsorbents[0].heat_of_adsorption["CO2"]`) is entered as a **positive number for exothermic adsorption**: the heat released per mole adsorbed, in J mol⁻¹. Its registry bounds are 0 to 10⁶. The energy balance adds the source

$$
S = \sum_k \sum_i \Delta H^{\mathrm{eff}}_{k,i}\thinspace \rho_{b,k}\thinspace \frac{\partial q_{k,i}}{\partial t} \quad [\mathrm{W\,m^{-3}}]
$$

on the heating side of the balance, so a loading bed ($\partial q/\partial t > 0$) warms and a desorbing bed cools. How $\Delta H^{\mathrm{eff}}$ follows from the stored value is set by `fluid.heat_of_adsorption_convention`:

| `fluid.heat_of_adsorption_convention` | The stored value is | $\Delta H^{\mathrm{eff}}$ |
|---|---|---|
| `"isosteric"` (default) | the isosteric heat $q_{st} = -\Delta H_{ads}$ | stored value |
| `"internal_energy"` | the internal-energy form of Haghpanah et al. (2013) | stored value $- C_{pa,i}\thinspace T$ |
| `"isosteric_rt"` | $-\Delta U_{ads}$ | stored value $+ R\thinspace T$ |

$T$ is the adsorbed-phase temperature (the bed temperature in LTE, the solid temperature in LTNE) and $C_{pa,i}$ the adsorbed-phase molar heat capacity (`fluid.adsorbed_heat_capacity_mode` and `fluid.adsorbed_heat_capacity[C]`; in the default constant mode it is 0 unless set, and then `internal_energy` gives the same result as `isosteric`).

**Worked sign example.** A paper reports the adsorption enthalpy of CO₂ on 13X as $\Delta H = -37$ kJ mol⁻¹ (negative: exothermic). Enter its magnitude:

```python
sim.adsorbents["13X"].heat_of_adsorption["CO2"] = 37e3     # J/mol, positive
sim.fluid.heat_of_adsorption_convention = "isosteric"       # the default
```

If the bed takes up CO₂ at $\partial q/\partial t = 10^{-3}$ mol kg⁻¹ s⁻¹ with $\rho_b = 660$ kg m⁻³, the source is $37\,000 \times 660 \times 10^{-3} = 24.4$ kW m⁻³ of heating. Had the paper instead reported $\Delta U = -35$ kJ mol⁻¹ in an energy balance written with $-(\Delta U - RT)$, enter `35e3` with `"isosteric_rt"`: at 300 K, $\Delta H^{\mathrm{eff}} = 35\,000 + 8.314 \times 300 = 37\,494$ J mol⁻¹. Entering a negative number would make adsorption cool the bed.

The energy balance always uses this stored value; it never derives a heat of adsorption from the temperature dependence of the isotherm. If the isotherm has its own energy parameter (for example the exponent of a temperature-dependent affinity), set it separately and consistently, in the sign convention of that isotherm's formula on the [Isotherms](isotherms.md) page.

**Where in the code:** `AdsorptionHeat` and the energy residuals in `core/src/equations/adsorbent_equations.cpp`; `HeatOfAdsorptionConvention` in `core/include/adsim/objects/fluid.h`.

## Composition used by the rate models

There is one gas phase per cell. Every rate model evaluates the equilibrium loading $q_i^*$ from the isotherm at the cell's **bulk gas** mole fractions $y$ and pressure $P$ and the adsorbed-phase (solid) temperature; there is no separate macropore or surface concentration.

The **capacity-factor** macropore rate (`set_rate_model(C, "ldf", "macropore_capacity", ...)`)

$$
k_i = \frac{15\thinspace \epsilon_p D_{m,i}}{\tau\thinspace r_p^2}\thinspace \frac{c_i}{\rho_s\thinspace q_i^*},
\qquad c_i = \frac{y_i P}{R T}
$$

uses:

- $c_i$: the cell's bulk gas concentration of component $i$, from that cell's $y_i$, $P$ and gas temperature $T$, floored at `kinetics[C].concentration_floor` (default 10⁻³ mol m⁻³);
- $q_i^*$: the isotherm at that same cell, with $y_i$ replaced by its floored value and every other component at its actual value (under the IAST mixing rule, the IAST mixture loading at that state);
- $\rho_s$: the adsorbent's `solid_density` (the particle density above);
- $\epsilon_p$: `kinetics[C].macropore_void_fraction` (not `particle_void_fraction`);
- the ratio $c_i/(\rho_s q_i^*)$ capped at 1.

So the rate is not evaluated at the feed composition or at a pore-gas composition, but at the local bulk gas of each cell at each instant.

**Where in the code:** `MassTransferCoefficientMacroporeCapacity::UpdateMassTransferCoefficient` (`core/src/equations/mass_transfer_model.cpp`).

## Boundary-value profiles

Each flow, inlet temperature and pressure at a step end can follow a profile in time. $t$ is the time since the start of the step, $v_0$ the start value, $v_1$ the end value, $t_r$ the ramp time.

| Kind (`sk.<kind>`) | Formula | After $t_r$ |
|---|---|---|
| `constant(v)` | $v(t) = v$ | — |
| `linear(end, start=, time=)` | $v_0 + (v_1 - v_0)\thinspace \min(t, t_r)/t_r$ | held at $v_1$ |
| `exponential(end, start=, time=, factor=)` | $v_1 + (v_0 - v_1)\thinspace F^{-t/t_r} = v_1 + (v_0 - v_1)\thinspace e^{-\ln(F)\thinspace t/t_r}$ | keeps approaching $v_1$ |
| `log_linear(end, start=, time=)` | $v_0\thinspace (v_1/v_0)^{\min(t, t_r)/t_r}$ (linear if $v_0 \le 0$ or $v_1 \le 0$) | held at $v_1$ |
| `double_exponential(end, rate_1, rate_2, weight, start=)` | $v_1 + (v_0 - v_1)\left[w\thinspace e^{-r_1 t} + (1-w)\thinspace e^{-r_2 t}\right]$ | keeps approaching $v_1$ |

$F$ is the **exponential factor** (`..._exponential_factor`, default 100, must exceed 1): after $t_r$ the change left is $1/F$ of the initial one, so $F = 100$ has covered 99 %. A ramp time of zero is an immediate step to $v_1$. To match a model written as $v_1 + (v_0 - v_1)e^{-t/\tau}$, keep any $t_r$ and set $F = e^{t_r/\tau}$; for example $\tau = 5$ s with $t_r = 20$ s gives $F = e^4 = 54.6$.

**Flows and inlet temperatures** use the start value and ramp time you give (`endpoint[ID].velocity_linear_start`, `..._linear_end`, `..._linear_time`, and so on).

**Pressures** (a pressure boundary, or a valve's set pressure) differ in two ways:

- the start value is the **bed pressure at the bed-end face at the start of the step** (the mean of the end cell and the ghost cell); there is no start parameter;
- the ramp time is the **step time** unless `endpoint[ID].pressure_<kind>_use_step_time` is false, in which case it is `endpoint[ID].pressure_<kind>_time`, which must not exceed the step time. For the exponential this means $F$ is reached at the end of the step by default.

```python
step.top.set_pressure("outlet", sk.exponential(1e4, factor=50))            # 1e4 Pa approached, 98 % by the step end
step.top.set_pressure("outlet", sk.exponential(1e4, time=10, factor=100))  # 99 % after 10 s
sim["endpoint[3].pressure_exponential_factor"] = 50                         # the same through the registry
```

**Where in the code:** `core/src/boundary_conditions/boundary_value.cpp` (the formulas); `PressureBoundaryCondition::Initialise` (`core/src/boundary_conditions/boundary_conditions.cpp`) and `core/src/boundary_conditions/boundary_stream.cpp` (the pressure start value and ramp time). [Boundary Conditions](boundary-conditions.md#time-varying-profiles) describes each profile in more detail.

## Velocity, volumetric and molar feeds

A flow boundary's `velocity` is the **superficial** velocity (flow over the empty column cross-section). Velocity and volumetric flow are converted to a molar flow with the ideal-gas concentration **at the bed-end face**, at its current pressure and temperature:

$$
\dot n = u\thinspace A\thinspace \frac{P_f}{R\thinspace T_f}, \qquad P_f = \tfrac12 (P_{\text{ghost}} + P_{\text{end cell}}), \quad T_f = \tfrac12 (T_{\text{ghost}} + T_{\text{end cell}})
$$

So a velocity or volumetric feed is an *actual* volumetric flow at the bed inlet, not at standard conditions, and its molar flow changes as the bed pressure changes during the step. A mass flow is divided by the face molar mass. If the other model specifies a molar flow, or a volumetric flow at standard conditions, use `sk.molar_flow` (convert a standard flow to mol s⁻¹ first) rather than `sk.velocity`.

**Where in the code:** `core/src/boundary_conditions/flow_boundary_value.cpp`.

## Where profile values sit

The bed is cut into finite volumes; each layer has `layer[ID].cells` equal cells. $z$ is measured from the bottom of the column ($z = 0$) to the top ($z = L$).

- **Cell variables** (pressure, temperature, mole fraction, loading) are cell averages, located at the **cell centres**: `sim.profile(...).z`, the same as `sim.positions`. The first centre is half a cell above the bottom, not at $z = 0$.
- **Face variables** (velocity, volumetric and molar flow, approach to fluidisation) are at the cell faces, from 0 to $L$ (`sim.face_positions`).
- **The bed ends**: use `sim.series(var, location="bottom")` or `location="top"`. For pressure, temperature and mole fraction this is the step endpoint's logged value: at an inlet the value the boundary imposes (feed composition and temperature, set pressure), at an outlet the bed-end face value. A closed end, and every other variable, gives the face value.
- **Any position**: `sim.series(var, z=0.15)` interpolates linearly between cell centres, and between the end cell and the bed-end value above in the last half-cell at each end.

```python
p = sim.profile("mole_fraction", "CarbonDioxide", step=0)
p.z[0]         # 0.5 * L / cells for a single layer: not the inlet
y_out = sim.series("mole_fraction", "CarbonDioxide", location="top")   # outlet value
```

When comparing with a code that stores nodes at the bed ends (a finite-difference or collocation grid with a node at $z = 0$ and $z = L$), compare its end nodes with `location=`, and its interior nodes with `series(z=...)` at the same $z$, not with `profile(...).values[:, 0]`.

**Where in the code:** `Simulation::Positions`, `FacePositions` and `Series` (`core/include/adsim/simulation.h`); `core/src/solver/bed_series.cpp`.

## Checklist: settings that most often differ

Go through these before concluding that two models disagree.

| Setting | Skarstrom default | Parameter |
|---|---|---|
| Advection scheme for the species balance | WENO (third order) | `settings.advection_scheme`: `"upwind"`, `"weno"`, `"van_leer"` |
| Number of cells | per layer | `layer[ID].cells` |
| Gas heat capacity | CoolProp, per component, evaluated once at the initial-condition pressure and temperature and held for the run | `fluid.gas_heat_capacity_mode` (`"coolprop"`, `"constant"`), `fluid.gas_heat_capacity[C]` in J mol⁻¹ K⁻¹ |
| Gas viscosity | mole-fraction average of pure-component $\mu_i(T)$ polynomials | `fluid.viscosity_mode` (`"polynomial"`, `"constant"`), `fluid.viscosity` in Pa s |
| Heat-of-adsorption form | isosteric | `fluid.heat_of_adsorption_convention` (see [above](#heat-of-adsorption-sign)) |
| Adsorbed-phase heat capacity | constant, 0 J mol⁻¹ K⁻¹ | `fluid.adsorbed_heat_capacity_mode`, `fluid.adsorbed_heat_capacity[C]` |
| Feed specification | as built | `sk.velocity` (actual flow at the bed face) vs `sk.molar_flow` (see [above](#velocity-volumetric-and-molar-feeds)) |
| Particle porosity in the gas hold-up | $\epsilon_t = \epsilon_b + (1-\epsilon_b)\epsilon_p$ | `adsorbent[ID].particle_void_fraction` (0 for $\epsilon_t = \epsilon_b$) |
| Density basis | particle density | `adsorbent[ID].solid_density` |
| Energy balance | non-isothermal, LTE | `energy_balance.energy_balance_mode`, `bed.thermal_equilibrium_mode` (`"lte"`, `"ltne"`) |
| Valve bed-side pressure | the bed-end face | `settings.valve_bed_pressure` (`"face"`, `"end_cell"`) |

Notes on the rows:

- **Advection scheme and cells.** Upwind smears fronts and needs several hundred cells for a recovery within 0.1 percentage point; WENO reaches the same on tens of cells ([Numerics](numerics.md#advection-schemes)). A reference code with first-order upwinding on 10-30 cells differs from Skarstrom's WENO default by its numerical dispersion: set `"upwind"` and the same cell count to reproduce it, or refine both. The scheme sets the convective face temperature of the energy balance as well as the face mole fractions.
- **Gas heat capacity.** The mixture value is mixed on a molar basis, $C_{p,g} = \sum_i y_i M_i C_{p,i} / \sum_i y_i M_i$ (J kg⁻¹ K⁻¹), in every cell. Many papers use a constant molar value (e.g. 30.7 J mol⁻¹ K⁻¹ for every component): set `"constant"` and that value for each component.
- **Viscosity.** A constant `fluid.viscosity` is used in every cell at every temperature and composition.
- **Dispersion.** The dispersive flux is $\epsilon_b D_{ax}\thinspace C\thinspace \partial y/\partial z$; a model that multiplies $D_{ax}$ by $\epsilon_t$, or not at all, needs its $D_{ax}$ rescaled. No dispersive flux crosses the bed-end faces (Danckwerts inlet; [Mass Balance](mass-balance.md#inlet-and-outlet-conditions)).

```python
sim.settings.advection_scheme = "upwind"
sim.layers[0].cells = 30
sim.fluid.gas_heat_capacity_mode = "constant"
for c in ("CO2", "N2"):
    sim.fluid.gas_heat_capacity[c] = 30.7          # J/(mol K)
sim.fluid.viscosity_mode = "constant"
sim.fluid.viscosity = 1.72e-5                      # Pa s
```

Finally, re-run on twice the cells: a disagreement that changes with the grid is a numerical one.
