## Equation of State and Gas Properties

The gas phase is treated as an **ideal-gas mixture**. This page lists exactly how each gas property in the model is obtained, following `core/src/objects/fluid.cpp`, `core/include/adsim/objects/fluid.h` and the property updates in `core/src/equations/adsorbent_equations.cpp`.

### Ideal gas law

In every cell the total molar concentration $C$, the pressure $P$ and the gas temperature $T$ are tied by

$$
C = \frac{P}{R\thinspace T}
$$

This is one of the residuals the solver satisfies at every time step (the *Concentration* equation), and the same relation is used wherever a concentration is needed at a cell face or a boundary, with the face pressure and temperature taken as the mean of the two neighbouring cells. Component concentrations follow from the mole fractions, $c_i = y_i\thinspace C$. There is no compressibility factor, so the model is meant for pressures and temperatures where the gas is close to ideal; at high pressure, such as the 15 bar of the Casas (2012) decks, check that the ideal-gas error is acceptable for your gas.

| Symbol  | Description                      | Unit        |
|---------|----------------------------------|-------------|
| $C$     | Total gas molar concentration    | mol m⁻³     |
| $c_i$   | Molar concentration of component $i$ | mol m⁻³ |
| $y_i$   | Mole fraction of component $i$   | –           |
| $P$     | Pressure                         | Pa          |
| $T$     | Gas temperature                  | K           |
| $R$     | Gas constant, 8.314462618 J mol⁻¹ K⁻¹ | J mol⁻¹ K⁻¹ |

### Where CoolProp is used

Components are chosen by their [CoolProp](https://coolprop.org/) fluid names in the deck's fluid package. CoolProp is called only when a run starts (`Fluid::Initialise`), never inside the time integration:

- **Molar mass** $M_i$ of each component, from CoolProp.
- **Heat capacity** $C_{p,i}$ [J kg⁻¹ K⁻¹] of each pure component, evaluated **once**, at the pressure and temperature of the **initial conditions**, and then held constant for the whole run. A TSA cycle that heats the gas by 100 K therefore uses the gas heat capacity of the initial temperature throughout.
- The heat capacity is always the **gas-phase** value. A component that would be condensed in its pure state at the initial pressure and temperature (water, ethanol, hexane or toluene at 1 bar and 298 K) takes its **ideal-gas** heat capacity $C_{p0,i}(T)$ instead of CoolProp's liquid value: about 1865 J kg⁻¹ K⁻¹ (33.6 J mol⁻¹ K⁻¹) for water vapour at 298 K, not the liquid's 4181. In the bed such a component is a vapour at a partial pressure no higher than its saturation pressure, where it is close to ideal, and $C_{p0}$ depends on temperature only, so it exists at every state. Components in a gas or supercritical state, which includes the permanent gases (N₂, O₂, CH₄, H₂, Ar, He, and CO₂ below its saturation pressure, about 64 bar at 298 K), keep the real-gas $C_{p,i}(P, T)$.

If CoolProp cannot evaluate a property, the error is printed to the console and the value is set to 0.

### Constant gas heat capacity and viscosity

Models in the literature often use fixed gas properties; for example, Yancy-Caballero et al. (2020) and Haghpanah et al. (2013) use $C_{p,g} = 30.7$ J mol⁻¹ K⁻¹ and $\mu = 1.72\times10^{-5}$ Pa s. Two fluid switches reproduce them:

- **Gas heat capacity mode** (`fluid.gas_heat_capacity_mode`): `coolprop` (the default, as above) or `constant`. With `constant`, the molar heat capacity of each component is the value given in `fluid.gas_heat_capacity[<component>]` [J mol⁻¹ K⁻¹], and $C_{p,i}$ is that value divided by $M_i$. The mixing rule below is unchanged, so equal values give that molar heat capacity at every composition. A new component starts at 29.1 J mol⁻¹ K⁻¹ (7R/2).
- **Viscosity mode** (`fluid.viscosity_mode`): `polynomial` (the default, the mole-fraction average below) or `constant`. With `constant`, every cell has the mixture viscosity `fluid.viscosity` [Pa s], whatever its temperature and composition.

The *Gas phase* adsorbed-phase heat capacity takes the same $C_{p,i}\thinspace M_i$, so with a constant gas heat capacity it is the constant too.

### Mixture properties

All mixture properties use **linear mole-fraction mixing**, evaluated in every cell at every residual evaluation:

$$
M = \sum_i y_i\thinspace M_i,\qquad
\rho = \sum_i c_i\thinspace M_i = C\thinspace M,\qquad
C_{p,g} = \frac{\sum_i y_i\thinspace M_i\thinspace C_{p,i}}{\sum_i y_i\thinspace M_i}
$$

$$
\mu = \sum_i y_i\thinspace \mu_i(T),\qquad
k_g = \sum_i y_i\thinspace k_i(T)
$$

The heat capacity is mixed on the molar basis: mole fractions average the pure-component **molar** heat capacities $M_i C_{p,i}$ (exact for an ideal-gas mixture), and dividing by the mixture molar mass $M$ gives the per-kilogram $C_{p,g}$ [J kg⁻¹ K⁻¹] that the energy balance multiplies by $\rho$. The volumetric heat capacity is therefore $\rho\thinspace C_{p,g} = C \sum_i y_i\thinspace M_i\thinspace C_{p,i}$. Every user of the gas heat capacity (the energy balances, the convective enthalpy terms, the thermal dispersion and heat-transfer correlations and the logged heat capacity) takes this one value. For 50/50 H₂/CO₂ it is about 33 J mol⁻¹ K⁻¹, where averaging the per-kilogram $C_{p,i}$ by mole fraction would give about 174.

The linear rules for viscosity and conductivity are simpler than the Wilke or Wassiljewa mixing rules; they are adequate for the dilute, chemically similar mixtures typical of PSA, and are what the model uses.

### Viscosity and thermal conductivity polynomials

Calling CoolProp in every cell at every residual evaluation would make a run minutes slower, so the pure-component viscosity $\mu_i$ and conductivity $k_i$ are degree-4 polynomials in temperature with coefficients fitted offline and compiled into the program (`viscosity_coefficients.h`, `thermal_conductivity_coefficients.h`):

$$
\mu_i(T) = \sum_{n=0}^{4} a_{n,i}\thinspace T_c^{\thinspace n},\qquad
k_i(T) = \sum_{n=0}^{4} b_{n,i}\thinspace T_c^{\thinspace n},\qquad
T_c = \min\negthinspace \big(\max(T,\thinspace 150\ \text{K}),\thinspace 650\ \text{K}\big)
$$

- The temperature is **clamped to 150–650 K** before evaluation: outside that window the property is held at its value at the nearest end, never extrapolated.
- The polynomials are least-squares fits of CoolProp's dilute-gas values, sampled at low pressure so that condensable components (water, alcohols, hydrocarbons) get their vapour values, not the liquid's (water at 298 K: 9.7 × 10⁻⁶ Pa s and 0.0184 W m⁻¹ K⁻¹). For the components CoolProp has no viscosity model for (CarbonMonoxide, Ethylene, Fluorine, Krypton, Neon, NitrousOxide, SulfurDioxide, Xenon), the viscosity fit is to Chapman–Enskog kinetic theory with the component's Lennard-Jones parameters; for the components without a conductivity model, the conductivity fit is to the Eucken relation $k = \mu\thinspace (c_p + 1.25\thinspace R/M)$. The headers state a maximum fit error of 1 % (viscosity) and 1.8 % (conductivity) over the window.
- Pressure dependence is neglected; gas viscosity changes by less than about 1 % over 0.1–10 bar.
- A floor of 10⁻⁶ Pa s (viscosity) and 10⁻⁴ W m⁻¹ K⁻¹ (conductivity) guards against a fit undershooting at the window edges. A component with no coefficients at all falls back to 1 × 10⁻⁵ Pa s and 0.024 W m⁻¹ K⁻¹.

### Other gas properties

- **Molecular diffusivities** come from the Chapman–Enskog binary-diffusion model with the Lennard-Jones energy and length of each component (Poling, Prausnitz and O'Connell, *The Properties of Gases and Liquids*, 5th ed.). They feed the dispersion and film mass-transfer correlations; see [Transport Correlations](transport-correlations.md).
- **Adsorbed-phase heat capacity** is either a constant per component or, with the *Gas phase* option of the energy balance, $C_{p,i}\thinspace M_i$ from the same initial-condition heat capacity; see [Energy Balance](energy-balance.md#adsorbed-phase-heat-capacity).
- **Isothermal runs** evaluate all of the above at the fixed temperature $T_{\text{iso}}$; see [Isothermal mode](energy-balance.md#isothermal-mode).
