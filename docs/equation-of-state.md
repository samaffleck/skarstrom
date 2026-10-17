## Equation of State and Gas Properties

The gas phase is treated as an **ideal-gas mixture**. This page lists exactly how each gas property in the model is obtained, following `core/src/objects/fluid.cpp`, `core/include/adsim/objects/fluid.h` and the property updates in `core/src/equations/adsorbent_equations.cpp`.

### Ideal gas law

In every cell the total molar concentration $C$, the pressure $P$ and the gas temperature $T$ are tied by

$$
C = \frac{P}{R\,T}
$$

This is one of the residuals the solver satisfies at every time step (the *Concentration* equation), and the same relation is used wherever a concentration is needed at a cell face or a boundary, with the face pressure and temperature taken as the mean of the two neighbouring cells. Component concentrations follow from the mole fractions, $c_i = y_i\,C$. There is no compressibility factor, so the model is meant for pressures and temperatures where the gas is close to ideal; at high pressure, such as the 15 bar of the Casas (2012) decks, check that the ideal-gas error is acceptable for your gas.

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

If CoolProp cannot evaluate a property, the error is printed to the console and the value is set to 0.

### Mixture properties

All mixture properties use **linear mole-fraction mixing**, evaluated in every cell at every residual evaluation:

$$
M = \sum_i y_i\,M_i,\qquad
\rho = \sum_i c_i\,M_i = C\,M,\qquad
C_{p,g} = \sum_i y_i\,C_{p,i}
$$

$$
\mu = \sum_i y_i\,\mu_i(T),\qquad
k_g = \sum_i y_i\,k_i(T)
$$

Mixing the heat capacity by mole fraction, as the code does, is exact for a molar heat capacity; with the per-kilogram $C_{p,i}$ that the model stores it is an approximation for components of very different molar mass. The linear rules for viscosity and conductivity are likewise simpler than the Wilke or Wassiljewa mixing rules; they are adequate for the dilute, chemically similar mixtures typical of PSA, and are what the model uses.

### Viscosity and thermal conductivity polynomials

Calling CoolProp in every cell at every residual evaluation would make a run minutes slower, so the pure-component viscosity $\mu_i$ and conductivity $k_i$ are degree-4 polynomials in temperature with coefficients fitted offline and compiled into the program (`viscosity_coefficients.h`, `thermal_conductivity_coefficients.h`):

$$
\mu_i(T) = \sum_{n=0}^{4} a_{n,i}\,T_c^{\,n},\qquad
k_i(T) = \sum_{n=0}^{4} b_{n,i}\,T_c^{\,n},\qquad
T_c = \min\!\big(\max(T,\,150\ \text{K}),\,650\ \text{K}\big)
$$

- The temperature is **clamped to 150–650 K** before evaluation: outside that window the property is held at its value at the nearest end, never extrapolated.
- The polynomials are least-squares fits of CoolProp's dilute-gas values. For the components CoolProp has no viscosity model for (CarbonMonoxide, Ethylene, Fluorine, Krypton, Neon, NitrousOxide, SulfurDioxide, Xenon), the viscosity fit is to Chapman–Enskog kinetic theory with the component's Lennard-Jones parameters; for the components without a conductivity model, the conductivity fit is to the Eucken relation $k = \mu\,(c_p + 1.25\,R/M)$. The headers state a maximum fit error of 1 % (viscosity) and 1.8 % (conductivity) over the window.
- Pressure dependence is neglected; gas viscosity changes by less than about 1 % over 0.1–10 bar.
- A floor of 10⁻⁶ Pa s (viscosity) and 10⁻⁴ W m⁻¹ K⁻¹ (conductivity) guards against a fit undershooting at the window edges. A component with no coefficients at all falls back to 1 × 10⁻⁵ Pa s and 0.024 W m⁻¹ K⁻¹.

### Other gas properties

- **Molecular diffusivities** come from the Chapman–Enskog binary-diffusion model with the Lennard-Jones energy and length of each component (Poling, Prausnitz and O'Connell, *The Properties of Gases and Liquids*, 5th ed.). They feed the dispersion and film mass-transfer correlations; see [Transport Correlations](transport-correlations.md).
- **Adsorbed-phase heat capacity** is either a constant per component or, with the *Gas phase* option of the energy balance, $C_{p,i}\,M_i$ from the same initial-condition heat capacity; see [Energy Balance](energy-balance.md#adsorbed-phase-heat-capacity).
- **Isothermal runs** evaluate all of the above at the fixed temperature $T_{\text{iso}}$; see [Isothermal mode](energy-balance.md#isothermal-mode).
