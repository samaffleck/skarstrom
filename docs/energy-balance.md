## Energy Balance

The deck's **energy balance mode** (parameter `energy_balance.energy_balance_mode`) chooses whether temperatures are solved at all:

- **Non-isothermal (solve energy balance)** — the default. One of the two thermal models below is solved, together with the [wall energy balance](wall-energy-balance.md).
- **Isothermal (fixed temperature)** — see [Isothermal mode](#isothermal-mode).

When the energy balance is solved, the two thermal models are selected with the **thermal equilibrium model** setting.

## Local Thermal Equilibrium (LTE)

This is the default. It assumes heat transfer between the gas and solid is fast enough that the two share a single temperature $T$, so one lumped energy equation is solved per cell.

$$
\left(\epsilon_t C_{p,g} \rho_g + C_{p,s} \rho_b\right) \frac{\partial T}{\partial t} - 
\lambda \frac{\partial^2 T}{\partial z^2} - 
\epsilon_t \frac{\partial P}{\partial t} + 
u C_{p,g} \rho_g \frac{\partial T}{\partial z} + 
\frac{4 h_i}{d_i} (T - T_w) -
\frac{RT}{P}\rho_b\sum_{i=1}^{N}(-\Delta H_i)\frac{\partial{q_i}}{\partial{t}}
= 0
$$

| Symbol          | Description                                         | Unit           |
|-----------------|-----------------------------------------------------|----------------|
| $C_{p,g}$       | Gas heat capacity                                   | J kg⁻¹ K⁻¹     |
| $C_{p,s}$       | Adsorbent heat capacity                             | J kg⁻¹ K⁻¹     |
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
| $q_i$           | Component i mass source                             | mol kg⁻¹ s⁻¹   |
| $\Delta H_i$    | Component i heat of adsorption                      | J mol⁻¹        |

The final term is the heat released (or absorbed) by adsorption, using the per-component isosteric heat of adsorption $\Delta H_i$ stored per adsorbent and component (parameter `adsorbent[<id>].heat_of_adsorption[<component>]`). The wall-coupling term $\tfrac{4 h_i}{d_i}(T - T_w)$ links this equation to the [wall energy balance](wall-energy-balance.md).

The effective axial thermal conductivity $\lambda$ and the internal coefficient $h_i$ are either constants or computed from the [transport correlations](transport-correlations.md).

## Local Thermal Non-Equilibrium (LTNE)

When gas–solid heat transfer is not fast (for example large particles or fast cycles), select **Local thermal non-equilibrium (LTNE)**. The solver then carries **separate fluid and solid temperature equations** coupled by a gas–particle heat transfer term. This requires a **Particle heat transfer model** and coefficient, described in the [transport correlations](transport-correlations.md#particle-heat-transfer-ltne-only). With LTNE, the solid temperature is available as its own plotted variable.

## Adsorbed-phase heat capacity

The **adsorbed phase heat capacity** setting of the energy balance controls how the adsorbed molecules' heat capacity is treated:

- **Constant** — enter a molar heat capacity per component.
- **Gas phase** — the adsorbed-phase heat capacity is taken equal to the gas-phase value.

## Energy source (heating / cooling)

Any step can apply a heat duty to a layer — the basis of **temperature swing (TSA)** operation. Each step has a list of energy sources, one per layer it heats or cools, each with a **power** [W]; positive values heat the layer, negative values cool it. The duty is distributed over the layer and added as a source term to the energy balance for the duration of that step.

A TSA feed can also be heated or cooled through its inlet: the inlet temperature of a stream or pressure boundary can follow a linear or exponential ramp during a step (see [Boundary Conditions](boundary-conditions.md#inlet-temperature-profile)).

## Isothermal mode

Many validation cases and short-cut studies assume an isothermal column. Selecting **Isothermal (fixed temperature)** removes the fluid, solid and wall temperatures from the set of unknowns altogether and holds them at

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
