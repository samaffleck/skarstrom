## Energy Balance

The deck's **energy balance mode** (parameter `energy_balance.energy_balance_mode`) chooses whether temperatures are solved at all:

- **Non-isothermal** — the default: the energy balance is solved. One of the two thermal models below is solved, together with the [wall energy balance](wall-energy-balance.md).
- **Isothermal** — the temperature is fixed; see [Isothermal mode](#isothermal-mode).

When the energy balance is solved, the two thermal models are selected with the **thermal equilibrium model** setting.

## Local Thermal Equilibrium (LTE)

This is the default. It assumes heat transfer between the gas and solid is fast enough that the two share a single temperature $T$, so one lumped energy equation is solved per cell.

$$
\begin{aligned}
&\Big(\epsilon_t \rho_g C_{p,g} + \rho_b C_{p,s} + \sum_k \rho_{b,k} \sum_{i=1}^{N} C_{pa,i}\thinspace q_{k,i}\Big) \frac{\partial T}{\partial t}
+ u \rho_g C_{p,g} \frac{\partial T}{\partial z}
- \frac{\partial}{\partial z}\left(\lambda \frac{\partial T}{\partial z}\right)
+ \frac{4 h_i}{d_i} (T - T_w)
- \epsilon_t \frac{\partial P}{\partial t} \\
&\qquad - \sum_k \rho_{b,k}\sum_{i=1}^{N}\Delta H_{k,i}^{\mathrm{eff}}\thinspace\frac{\partial q_{k,i}}{\partial t}
- Q_e
- Q_r
- a_x h_x \left(T_m - T\right)
= 0
\end{aligned}
$$

| Symbol          | Description                                         | Unit           |
|-----------------|-----------------------------------------------------|----------------|
| $C_{p,g}$       | Gas heat capacity, ideal-gas $C_{p0}(T)$ at the local temperature ([Equation of State](equation-of-state.md#gas-heat-capacity)) | J kg⁻¹ K⁻¹ |
| $C_{p,s}$       | Heat capacity of the layer's solid (mass-fraction average over the solids, [Mixed Layers](mixed-layers.md#mixing-rules)) | J kg⁻¹ K⁻¹     |
| $C_{pa,i}$      | Adsorbed-phase molar heat capacity of component i ([below](#adsorbed-phase-heat-capacity)) | J mol⁻¹ K⁻¹ |
| $T$             | Gas/Adsorbent temperature                           | K              |
| $T_w$           | Wall temperature                                    | K              |
| $P$             | Pressure                                             | Pa             |
| $u$             | Superficial velocity                                | m s⁻¹          |
| $\epsilon_t$    | Total void fraction (bed + intra-particle)          | -              |
| $\rho_g$        | Gas density                                         | kg m⁻³         |
| $\rho_b$        | Bulk density of the layer's solids, $\sum_k \rho_{b,k}$ | kg m⁻³         |
| $\rho_{b,k}$    | Bulk density of solid $k$ in the layer (mass of solid $k$ per bed volume) | kg m⁻³ |
| $\lambda$       | Effective axial thermal conductivity                | W m⁻¹ K⁻¹      |
| $h_i$           | Internal (bed-to-wall) heat transfer coefficient     | W m⁻² K⁻¹      |
| $d_i$           | Column internal diameter (the hydraulic diameter for a non-circular column) | m              |
| $q_{k,i}$       | Adsorbed loading of component i on solid $k$         | mol kg⁻¹       |
| $\Delta H_{k,i}^{\mathrm{eff}}$ | Heat released per mole of component i adsorbed on solid $k$ (positive for exothermic adsorption) | J mol⁻¹ |
| $Q_e$           | Step [energy source](#energy-source-heating-cooling) power per unit bed volume | W m⁻³ |
| $Q_r$           | [Heat of reaction](#heat-of-reaction) released per unit bed volume | W m⁻³ |
| $a_x h_x$, $T_m$ | Internal heat exchanger conductance per bed volume and medium temperature ([Wall Energy Balance](wall-energy-balance.md#internal-heat-exchanger)); 0 outside a step with an exchanger and in a void region | W m⁻³ K⁻¹, K |

The sums run over the solids $k$ of the layer ([Mixed Layers](mixed-layers.md)); a layer of one solid has one term, with $\rho_{b,k} = \rho_b$. The adsorption term is the heat released by adsorption: $\Delta H_{k,i}^{\mathrm{eff}}$ follows from the per-component heat of adsorption stored per adsorbent (parameter `adsorbent[<id>].heat_of_adsorption[<component>]`, entered as a **positive** number for exothermic adsorption) according to the [heat of adsorption convention](#heat-of-adsorption-convention); with the default isosteric convention it is the stored value itself. So a loading bed ($\partial q_{k,i}/\partial t > 0$) heats up and a desorbing bed cools. $q_{k,i}$ is the total loading of the solid, including the slow site of a dual-rate model. The convective term uses the face temperature reconstructed with the selected advection scheme, as the mole fractions are, and no heat is conducted through the two bed ends; [Numerics](numerics.md) gives the discrete form. [Matching Another Model](model-conventions.md#heat-of-adsorption-sign) has a worked sign example. The wall-coupling term $\tfrac{4 h_i}{d_i}(T - T_w)$ links this equation to the [wall energy balance](wall-energy-balance.md).

The effective axial thermal conductivity $\lambda$ and the internal coefficient $h_i$ are either constants or computed from the [transport correlations](transport-correlations.md).

## Local Thermal Non-Equilibrium (LTNE)

When gas–solid heat transfer is not fast (for example large particles or fast cycles), select **Local thermal non-equilibrium** (LTNE). The solver then carries **separate fluid and solid temperature equations** coupled by a gas–particle heat transfer term. This requires a **Particle heat transfer model** and coefficient, described in the [transport correlations](transport-correlations.md#particle-heat-transfer-ltne-only). With LTNE, the solid temperature is available as its own plotted variable.

The gas temperature $T$ follows

$$
\begin{aligned}
&\epsilon_t \rho_g C_{p,g} \frac{\partial T}{\partial t}
+ u \rho_g C_{p,g} \frac{\partial T}{\partial z}
- \frac{\partial}{\partial z}\left(\lambda \frac{\partial T}{\partial z}\right)
+ g\thinspace\frac{4 h_i}{d_i} (T - T_w)
+ h_p^{vol} (T - T_s) \\
&\qquad - \epsilon_t \frac{\partial P}{\partial t}
- Q_e
- \epsilon_b\thinspace a_x h_x \left(T_m - T\right)
= 0
\end{aligned}
$$

and the solid temperature $T_s$ (in every bed cell; a void region has no solid)

$$
\begin{aligned}
&\Big(\rho_b C_{p,s} + \sum_k \rho_{b,k} \sum_{i=1}^{N} C_{pa,i}\thinspace q_{k,i}\Big) \frac{\partial T_s}{\partial t}
- \frac{\partial}{\partial z}\left(\lambda_s \frac{\partial T_s}{\partial z}\right)
- h_p^{vol} (T - T_s)
+ (1 - g)\thinspace\frac{4 h_i}{d_i} (T_s - T_w) \\
&\qquad - \sum_k \rho_{b,k}\sum_{i=1}^{N}\Delta H_{k,i}^{\mathrm{eff}}\thinspace\frac{\partial q_{k,i}}{\partial t}
- Q_r
- (1 - \epsilon_b)\thinspace a_x h_x \left(T_m - T_s\right)
= 0
\end{aligned}
$$

| Symbol | Description | Unit |
|---|---|---|
| $T_s$ | Solid (adsorbent) temperature | K |
| $h_p^{vol}$ | Volumetric gas–solid heat transfer coefficient, $h_p a_v$ ([Transport Correlations](transport-correlations.md#particle-heat-transfer-ltne-only)) | W m⁻³ K⁻¹ |
| $g$ | Gas share of the bed-to-wall exchange: $\epsilon_b$ in a packed layer, 0 in a structured layer ([Wall Energy Balance](wall-energy-balance.md)) | – |
| $\epsilon_b$ | Bed void fraction (the open frontal area $\phi$ in a structured layer) | – |
| $\lambda_s$ | Axial conductivity of the solid: $(1-\phi)\thinspace k_s$ in a structured layer, 0 in a packed layer ([below](#solid-conduction-in-structured-layers)) | W m⁻¹ K⁻¹ |

The other symbols are those of the LTE equation. How the terms are shared between the two equations:

- **Gas**: the whole gas hold-up $\epsilon_t \rho_g C_{p,g}$ (the gas in the particle pores included), the convection, the pressure work $\epsilon_t\thinspace\partial P/\partial t$ and the step energy source $Q_e$: an [energy source](#energy-source-heating-cooling) heats the gas only.
- **Solid**: the solid and adsorbed-phase heat capacity, the heat of adsorption, evaluated at $T_s$ ([below](#heat-of-adsorption-convention)), and the [heat of reaction](#heat-of-reaction), released on the catalyst.
- **Shared**: the bed-to-wall exchange, the gas taking the share $g$ and the solid $1 - g$; the internal heat exchanger, the gas taking the share $\epsilon_b$ and the solid $1 - \epsilon_b$. In a structured layer the exchanger's gas share is the open frontal area $\phi$, although the gas takes no part of the wall exchange there.

The gas equation's $\lambda$ is the layer's axial thermal conductivity model. In a packed layer the solid has no conduction term, so the whole of $\lambda$ sits in the gas equation; with Yagi and Kunii that includes their stagnant term, which describes conduction through the packing. In a structured layer the channel conduction model splits the conduction between the gas and the solid (below).

### Solid conduction in structured layers

In a monolith or laminate the channel walls are continuous along the column, so under LTNE the solid energy balance carries its own axial conduction, the term $-\partial_z(\lambda_s\thinspace\partial_z T_s)$ above, with

$$
\lambda_s = (1-\phi)\thinspace k_s,
$$

$\phi$ the open frontal area and $k_s$ the layer's effective solid conductivity (the solid-volume average over the solids of its composition, [Mixed Layers](mixed-layers.md#mixing-rules)). The face between two cells conducts as their two half cells in series, $G = 1/\left(\Delta z_j/(2\lambda_{s,j}) + \Delta z_{j+1}/(2\lambda_{s,j+1})\right)$, the harmonic mean of $\lambda_s$ on a uniform grid. A packed bed's particles touch only at points, so $\lambda_s = 0$ there and in a void region: no heat crosses a face that touches a packed layer or a void, nor the two bed ends, and a bed without a structured layer has no such term. The gas equation's conductivity is then the gas's alone, $\phi\thinspace k_g$ (the [channel conduction](transport-correlations.md#channel-conduction-structured-layers) model). Under LTE that model gives the single temperature the sum $\phi\thinspace k_g + (1-\phi)\thinspace k_s$ instead.

The bed-to-wall exchange of a structured layer goes into the solid ($g = 0$): the gas reaches the column wall only through the channel walls ([Wall Energy Balance](wall-energy-balance.md)). LTE is allowed with a structured layer, but the solver warns, since it takes the exchange between the gas and the thin channel walls as instantaneous; LTNE is recommended.

## Heat of adsorption convention

The stored value is positive for exothermic adsorption; [Matching Another Model](model-conventions.md#heat-of-adsorption-sign) has a worked sign example. The value stored as `adsorbent[<id>].heat_of_adsorption[<component>]` is read according to the fluid's **heat of adsorption convention** (`fluid.heat_of_adsorption_convention`). The convention sets the heat released per mole adsorbed, $\Delta H_{k,i}^{\mathrm{eff}}$, that multiplies $\rho_{b,k}\thinspace \partial q_{k,i} / \partial t$ in the source term; $\Delta H_i$ below is the stored value of the adsorbent, positive for exothermic adsorption:

| Convention | Stored value | $\Delta H_i^{\mathrm{eff}}$ | Used by |
|---|---|---|---|
| `isosteric` (default) | isosteric heat $q_{st} = -\Delta H_{ads}$ (positive) | $\Delta H_i$ | Casas & Mazzotti, Ruthven |
| `internal_energy` | internal-energy form | $\Delta H_i - C_{pa,i}\thinspace T$ | Haghpanah 2013, Ward & Pini |
| `isosteric_rt` | $-\Delta U_i$, the internal energy of adsorption | $\Delta H_i + R\thinspace T$ | Yancy-Caballero 2020 (PEESEgroup PSA) |

$T$ is the temperature of the adsorbed phase: the bed temperature in LTE, the solid temperature in LTNE. With `isosteric_rt`, the isosteric heat follows the ideal-gas relation $q_{st} = -\Delta U + RT$ at the local temperature, where a constant `isosteric` value fixes it at one temperature. For CO₂ on 13X at 313–363 K the two differ by $R\thinspace\Delta T \approx 0.4$ kJ/mol, about 1 %.

## Adsorbed-phase heat capacity

The **adsorbed phase heat capacity** setting of the energy balance controls how the adsorbed molecules' heat capacity is treated:

- **Constant** — enter a molar heat capacity per component.
- **Gas phase** — each component's adsorbed-phase heat capacity is taken equal to its own pure-gas molar heat capacity, $C_{p,i}\thinspace M_i$, at the temperature of the adsorbed phase (see [Equation of State](equation-of-state.md#gas-heat-capacity)). The $C_{pa,i} T$ term of the internal-energy convention uses that local value too.

## Energy source (heating / cooling)

Any step can apply a heat duty to a layer — the basis of **temperature swing (TSA)** operation. Each step has a list of energy sources, one per layer it heats or cools, each with a **power** [W]; positive values heat the layer, negative values cool it. The duty is spread uniformly over the layer, $Q_e = \text{power}/V_{\text{layer}}$, and added as a source term to the energy balance for the duration of that step. Under LTNE it enters the gas equation only ([above](#local-thermal-non-equilibrium-ltne)).

A TSA feed can also be heated or cooled through its inlet: the inlet temperature of a stream or pressure boundary can follow a linear or exponential ramp during a step (see [Boundary Conditions](boundary-conditions.md#inlet-temperature-profile)).

## Heat of reaction

A [reaction](reactions.md) on a solid releases its heat into the energy balance as the volumetric source

$$
Q_r = -\sum_r \Delta H_r(T_{cat})\thinspace \rho_{b,k}\thinspace \eta_r\thinspace r_r
$$

summed over the reactions in the cell, each with the bulk density $\rho_{b,k}$ of the solid $k$ that carries it.

| Symbol | Description | Unit |
|---|---|---|
| $\Delta H_r$ | Heat of reaction per mole of extent, **negative for an exothermic reaction** | J mol⁻¹ |
| $r_r$ | Intrinsic rate per kg of the solid that carries the reaction | mol kg⁻¹ s⁻¹ |
| $\eta_r$ | Effectiveness factor (a constant, or from the Thiele modulus; [Reactions](reactions.md)) | – |
| $T_{cat}$ | Catalyst temperature: $T$ under LTE, $T_s$ under LTNE | K |

The sign convention is the opposite of the heat of adsorption's: an exothermic reaction has $\Delta H_r < 0$ and gives $Q_r > 0$, while an exothermic adsorption has a positive stored heat of adsorption. The reaction's **heat mode** sets how $\Delta H_r$ depends on temperature:

- **Kirchhoff** (the default): $\Delta H_r(T) = \Delta H_{ref} + \sum_j \nu_j \int_{T_{ref}}^{T} C_{p,j}\thinspace dT$, with the molar gas heat capacities of the energy balance (for constant heat capacities, $\Delta C_p\thinspace(T - T_{ref})$). $T_{ref}$ is the reference temperature of the reaction, 298.15 K by default.
- **Constant**: $\Delta H_r = \Delta H_{ref}$.

Under LTE $Q_r$ enters the single energy equation; under LTNE it is released on the catalyst, in the solid equation. $Q_r$ is 0 in a void region and in an isothermal run.

## Isothermal mode

Many validation cases and short-cut studies assume an isothermal column. Selecting **Isothermal** removes the fluid, solid and wall temperatures from the set of unknowns altogether and holds them at

$$
T = T_s = T_w = T_{\text{iso}}, \qquad \frac{\partial T}{\partial t} = 0
$$

for the whole run. $T_{\text{iso}}$ is the initial-condition temperature by default; untick **Use the initial-condition temperature** to enter a different value.

Because no energy equation is solved:

- the heat of adsorption, the heat of reaction, the bed–wall and wall–ambient heat transfer, the wall temperature boundary condition, a step's heating medium and internal heat exchanger, and any step **Energy Sources** have no effect (the logged energy sources and medium duties are zero);
- inlet stream temperatures (and inlet temperature ramps) are ignored — every gas entering or leaving the column, and every stream-side temperature used by valves and pumps, is at $T_{\text{iso}}$; compressors do not heat the gas they deliver;
- compressor and vacuum-pump **power** is still computed, using $T_{\text{iso}}$ as the suction temperature;
- every logged temperature (gas, solid and wall) equals $T_{\text{iso}}$.

Gas properties (density, viscosity, diffusivities), the isotherms and the reaction rates are all evaluated at $T_{\text{iso}}$. Dropping the three temperature unknowns per cell makes the system smaller (and its banded Jacobian narrower), so isothermal runs are faster than an equivalent non-isothermal run with the heat-transfer coefficients and heats of adsorption set to zero, and give the same answer as that limit.

Existing decks keep solving the energy balance; the isothermal switch is stored in the deck only from this version on.
