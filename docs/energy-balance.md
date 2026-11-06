## Energy Balance

The deck's **energy balance mode** (parameter `energy_balance.energy_balance_mode`) chooses whether temperatures are solved at all:

- **Non-isothermal** — the default: the energy balance is solved. One of the two thermal models below is solved, together with the [wall energy balance](wall-energy-balance.md).
- **Isothermal** — the temperature is fixed; see [Isothermal mode](#isothermal-mode).

When the energy balance is solved, the two thermal models are selected with the **thermal equilibrium model** setting.

## Local Thermal Equilibrium (LTE)

This is the default. It assumes heat transfer between the gas and solid is fast enough that the two share a single temperature $T$, so one lumped energy equation is solved per cell.

$$
\left(\epsilon_t \rho_g C_{p,g} + \rho_b C_{p,s} + \rho_b \sum_{i=1}^{N} C_{pa,i}\thinspace q_i\right) \frac{\partial T}{\partial t}
+ u \rho_g C_{p,g} \frac{\partial T}{\partial z}
- \lambda \frac{\partial^2 T}{\partial z^2}
+ \frac{4 h_i}{d_i} (T - T_w)
- \epsilon_t \frac{\partial P}{\partial t}
- \rho_b\sum_{i=1}^{N}\Delta H_i^{\mathrm{eff}}\thinspace\frac{\partial q_i}{\partial t}
- Q_e
= 0
$$

| Symbol          | Description                                         | Unit           |
|-----------------|-----------------------------------------------------|----------------|
| $C_{p,g}$       | Gas heat capacity                                   | J kg⁻¹ K⁻¹     |
| $C_{p,s}$       | Adsorbent heat capacity                             | J kg⁻¹ K⁻¹     |
| $C_{pa,i}$      | Adsorbed-phase molar heat capacity of component i ([below](#adsorbed-phase-heat-capacity)) | J mol⁻¹ K⁻¹ |
| $T$             | Gas/Adsorbent temperature                           | K              |
| $T_w$           | Wall temperature                                    | K              |
| $P$             | Pressure                                             | Pa             |
| $u$             | Superficial velocity                                | m s⁻¹          |
| $\epsilon_t$    | Total void fraction (bed + intra-particle)          | -              |
| $\rho_g$        | Gas density                                         | kg m⁻³         |
| $\rho_b$        | Adsorbent bulk density                              | kg m⁻³         |
| $\lambda$       | Effective axial thermal conductivity                | W m⁻¹ K⁻¹      |
| $h_i$           | Internal heat transfer coefficient (gas–wall)       | W m⁻² K⁻¹      |
| $d_i$           | Wall internal diameter                              | m              |
| $q_i$           | Adsorbed loading of component i                     | mol kg⁻¹       |
| $\Delta H_i^{\mathrm{eff}}$ | Heat released per mole of component i adsorbed (positive for exothermic adsorption) | J mol⁻¹ |
| $Q_e$           | Step [energy source](#energy-source-heating-cooling) power per unit bed volume | W m⁻³ |

The adsorption term is the heat released by adsorption: $\Delta H_i^{\mathrm{eff}}$ follows from the per-component heat of adsorption stored per adsorbent (parameter `adsorbent[<id>].heat_of_adsorption[<component>]`, entered as a **positive** number for exothermic adsorption) according to the [heat of adsorption convention](#heat-of-adsorption-convention); with the default isosteric convention it is the stored value itself. So a loading bed ($\partial q_i/\partial t > 0$) heats up and a desorbing bed cools. In a layer of several adsorbents, $\rho_b C_{pa,i} q_i$ and $\rho_b \Delta H_i^{\mathrm{eff}} \partial q_i/\partial t$ are summed over the adsorbents, each with its own bulk density and loading. The convective term uses the upwind face temperature ([Numerics](numerics.md#advection-schemes)), and no heat is conducted through the two bed ends. [Matching Another Model](model-conventions.md#heat-of-adsorption-sign) has a worked sign example. The wall-coupling term $\tfrac{4 h_i}{d_i}(T - T_w)$ links this equation to the [wall energy balance](wall-energy-balance.md).

The effective axial thermal conductivity $\lambda$ and the internal coefficient $h_i$ are either constants or computed from the [transport correlations](transport-correlations.md).

## Local Thermal Non-Equilibrium (LTNE)

When gas–solid heat transfer is not fast (for example large particles or fast cycles), select **Local thermal non-equilibrium** (LTNE). The solver then carries **separate fluid and solid temperature equations** coupled by a gas–particle heat transfer term. This requires a **Particle heat transfer model** and coefficient, described in the [transport correlations](transport-correlations.md#particle-heat-transfer-ltne-only). With LTNE, the solid temperature is available as its own plotted variable.

## Heat of adsorption convention

The stored value is positive for exothermic adsorption; [Matching Another Model](model-conventions.md#heat-of-adsorption-sign) has a worked sign example. The value stored as `adsorbent[<id>].heat_of_adsorption[<component>]` is read according to the fluid's **heat of adsorption convention** (`fluid.heat_of_adsorption_convention`). The convention sets the heat released per mole adsorbed, $\Delta H_i^{\mathrm{eff}}$, that multiplies $\rho_b\thinspace \partial q_i / \partial t$ in the source term; $\Delta H_i$ below is the stored value, positive for exothermic adsorption:

| Convention | Stored value | $\Delta H_i^{\mathrm{eff}}$ | Used by |
|---|---|---|---|
| `isosteric` (default) | isosteric heat $q_{st} = -\Delta H_{ads}$ (positive) | $\Delta H_i$ | Casas & Mazzotti, Ruthven |
| `internal_energy` | internal-energy form | $\Delta H_i - C_{pa,i}\thinspace T$ | Haghpanah 2013, Ward & Pini |
| `isosteric_rt` | $-\Delta U_i$, the internal energy of adsorption | $\Delta H_i + R\thinspace T$ | Yancy-Caballero 2020 (PEESEgroup PSA) |

$T$ is the temperature of the adsorbed phase: the bed temperature in LTE, the solid temperature in LTNE. With `isosteric_rt`, the isosteric heat follows the ideal-gas relation $q_{st} = -\Delta U + RT$ at the local temperature, where a constant `isosteric` value fixes it at one temperature. For CO₂ on 13X at 313–363 K the two differ by $R\thinspace\Delta T \approx 0.4$ kJ/mol, about 1 %.

## Adsorbed-phase heat capacity

The **adsorbed phase heat capacity** setting of the energy balance controls how the adsorbed molecules' heat capacity is treated:

- **Constant** — enter a molar heat capacity per component.
- **Gas phase** — each component's adsorbed-phase heat capacity is taken equal to its own pure-gas molar heat capacity, $C_{p,i}\thinspace M_i$ (see [Equation of State](equation-of-state.md#other-gas-properties)).

## Energy source (heating / cooling)

Any step can apply a heat duty to a layer — the basis of **temperature swing (TSA)** operation. Each step has a list of energy sources, one per layer it heats or cools, each with a **power** [W]; positive values heat the layer, negative values cool it. The duty is distributed over the layer and added as a source term to the energy balance for the duration of that step.

A TSA feed can also be heated or cooled through its inlet: the inlet temperature of a stream or pressure boundary can follow a linear or exponential ramp during a step (see [Boundary Conditions](boundary-conditions.md#inlet-temperature-profile)).

## Isothermal mode

Many validation cases and short-cut studies assume an isothermal column. Selecting **Isothermal** removes the fluid, solid and wall temperatures from the set of unknowns altogether and holds them at

$$
T = T_s = T_w = T_{\text{iso}}, \qquad \frac{\partial T}{\partial t} = 0
$$

for the whole run. $T_{\text{iso}}$ is the initial-condition temperature by default; untick **Use the initial-condition temperature** to enter a different value.

Because no energy equation is solved:

- the heat of adsorption, the gas–wall and wall–ambient heat transfer, the wall temperature boundary condition and any step **Energy Sources** have no effect (the logged energy sources are zero);
- inlet stream temperatures (and inlet temperature ramps) are ignored — every gas entering or leaving the column, and every stream-side temperature used by valves and pumps, is at $T_{\text{iso}}$; compressors do not heat the gas they deliver;
- compressor and vacuum-pump **power** is still computed, using $T_{\text{iso}}$ as the suction temperature;
- every logged temperature (gas, solid and wall) equals $T_{\text{iso}}$.

Gas properties (density, viscosity, diffusivities) and the isotherms are all evaluated at $T_{\text{iso}}$. Dropping the three temperature unknowns per cell makes the system smaller (and its banded Jacobian narrower), so isothermal runs are faster than an equivalent non-isothermal run with the heat-transfer coefficients and heats of adsorption set to zero, and give the same answer as that limit.

Existing decks keep solving the energy balance; the isothermal switch is stored in the deck only from this version on.
