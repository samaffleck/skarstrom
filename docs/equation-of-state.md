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

### Where the component constants come from

Components are chosen by their [CoolProp](https://coolprop.org/) fluid names in the deck's fluid package. The engine does not call CoolProp. The fixed constants of each component are in a table, `core/include/adsim/objects/component_properties.h`. The script `scripts/gen_component_properties.py` writes this table from the CoolProp fluid files:

- **Molar mass** $M_i$ of each component.
- **Formula, CAS number and aliases** (the component manifest).
- **Critical point, triple point and vapour-pressure equation** (the saturation pressure of the GAB isotherms; see [Isotherms](isotherms.md)).

The heat capacity, viscosity and thermal conductivity come from polynomials; see below.

!!! note "New in 0.10"
    The table replaces the CoolProp calls of earlier engines. The values are the same.

### Gas heat capacity

The heat capacity of each pure component is its **ideal-gas** heat capacity $C_{p0,i}(T)$, evaluated at the local gas temperature in every cell at every residual evaluation, from the NASA 7-coefficient polynomials:

$$
\frac{C_{p0,i}}{R} = a_0 + a_1 T + a_2 T^2 + a_3 T^3 + a_4 T^4,\qquad
\frac{H_{0,i}}{R\thinspace T} = a_0 + \frac{a_1 T}{2} + \frac{a_2 T^2}{3} + \frac{a_3 T^3}{4} + \frac{a_4 T^4}{5} + \frac{a_5}{T}
$$

with one coefficient set per component, valid to 1000 K: the low-temperature set of TM-4513 (adsorption processes stay well below 1000 K, so the report's 1000–6000 K set is not carried). The coefficients are those of McBride, Gordon and Reno (1993), *Coefficients for Calculating Thermodynamic and Transport Properties of Individual Species*, NASA TM-4513, copied unchanged with each species' data source (`core/include/adsim/objects/ideal_gas_heat_capacity.h`). The reference and each species' TM-4513 label and data-source code are available through the API (`skarstrom.ideal_gas_heat_capacity_reference()` and `skarstrom.ideal_gas_heat_capacity_source(component)`; C++ `kIdealGasHeatCapacityReference`, `NasaPolynomial::species` and `source`). Air, which TM-4513 does not list, is the ideal-gas mixture of its N₂, O₂ and Ar entries at the composition of CoolProp's air model (Lemmon et al. 2000: 0.7812, 0.2096, 0.0092 by mole); an ideal-gas mixture mixes the coefficients linearly, so this is exact. Outside a polynomial's range (below 200 K for most species, 298 or 300 K for a few, and above 1000 K) the value at the nearer limit is used.

Why ideal-gas and temperature-dependent:

- **The initial temperature is only a start.** Earlier engines evaluated CoolProp's heat capacity once, at the initial-condition pressure and temperature, and held it for the run. The initial temperature then also set a model parameter, so two runs differing only in their initial temperature reached different cyclic steady states (for Ward 2024 with a 40 K colder start, recovery differed by 2.4×10⁻⁴). A TSA cycle that heats the gas by 100 K also kept the heat capacity of the initial temperature throughout.
- **Gas-phase values for condensables.** A vapour such as water, ethanol or toluene takes its vapour heat capacity whatever its pure state at the bed's pressure and temperature (water at 298 K: 33.6 J mol⁻¹ K⁻¹, not the liquid's 75.3).
- **Ideal gas.** The equation of state is the ideal-gas law, and the real-gas correction to $C_p$ of the permanent gases is small at PSA pressures (at 1 bar and 298 K: CO₂ +0.8 %, CH₄ +0.3 %, N₂ +0.2 %).

The energy balance is written in temperature form, $\rho C_{p,g}(\partial T/\partial t + u\thinspace\partial T/\partial z)$ with the local $C_{p,g}(T)$, which is exact for an ideal gas whose enthalpy is $h = \int C_{p0}\thinspace dT$: no extra term is needed. The open-system energy test (`ClosedBedEnergyTest.OpenSystemEnergyBalanceConvergesUnderGridRefinement`) checks the closure with this enthalpy.

**Agreement with CoolProp.** Each polynomial is compared with CoolProp's ideal-gas $C_{p0}$ (the ideal-gas part of the reference equation of state) from 200 K, or from the fit's lower limit or the fluid's triple point if higher, to 800 K (`test_ideal_gas_heat_capacity.cpp`, which prints this table). The test fails if an error exceeds 1 % for the light and permanent gases (Air, Ar, CO₂, CO, He, H₂, Kr, CH₄, Ne, N₂, O₂, Xe) or 2 % for any other component.

| Component | Polynomial range [K] | Compared [K] | Largest error | at [K] | $C_{p0}$(298.15 K) [J mol⁻¹ K⁻¹] |
|---|---|---|---|---|---|
| Air | 200–1000 | 200–800 | −0.08 % | 675 | 29.10 |
| Ammonia | 200–1000 | 200–800 | +0.57 % | 650 | 35.63 |
| Argon | 200–1000 | 200–800 | 0.00 % | — | 20.79 |
| Benzene | 200–1000 | 279–800 | −1.35 % | 424 | 81.92 |
| Butane | 200–1000 | 200–800 | −0.74 % | 200 | 98.66 |
| CarbonDioxide | 200–1000 | 217–800 | −0.14 % | 382 | 37.14 |
| CarbonMonoxide | 200–1000 | 200–800 | +0.13 % | 420 | 29.14 |
| Cyclohexane | 200–1000 | 280–800 | −1.81 % | 380 | 105.34 |
| Ethane | 200–1000 | 200–800 | −0.49 % | 200 | 52.50 |
| Ethanol | 200–1000 | 200–800 | −0.53 % | 415 | 65.31 |
| Ethylene | 200–1000 | 200–800 | −0.78 % | 200 | 42.89 |
| Fluorine | 200–1000 | 200–800 | +0.24 % | 655 | 31.30 |
| Helium | 200–1000 | 200–800 | 0.00 % | — | 20.79 |
| Hydrogen | 200–1000 | 200–800 | +0.91 % | 200 | 28.84 |
| HydrogenSulfide | 300–1000 | 300–800 | −0.29 % | 800 | 34.21 |
| Isobutane | 200–1000 | 200–800 | −0.65 % | 200 | 96.64 |
| Isopentane | 298–1000 | 298–800 | −0.15 % | 493 | 118.85 |
| Krypton | 200–1000 | 200–800 | 0.00 % | — | 20.79 |
| Methane | 200–1000 | 200–800 | +0.81 % | 745 | 35.69 |
| Methanol | 200–1000 | 200–800 | −0.38 % | 425 | 44.04 |
| Neon | 200–1000 | 200–800 | 0.00 % | — | 20.79 |
| Nitrogen | 200–1000 | 200–800 | −0.11 % | 685 | 29.12 |
| NitrousOxide | 200–1000 | 200–800 | +0.28 % | 210 | 38.63 |
| Octane | 200–1000 | 217–800 | −1.28 % | 217 | 187.78 |
| Oxygen | 200–1000 | 200–800 | +0.06 % | 365 | 29.38 |
| Pentane | 298–1000 | 298–800 | +0.50 % | 798 | 119.95 |
| Propane | 200–1000 | 200–800 | +0.69 % | 665 | 73.59 |
| Propylene | 200–1000 | 200–800 | +0.35 % | 665 | 64.43 |
| SulfurDioxide | 300–1000 | 300–800 | +0.10 % | 365 | 39.94 |
| SulfurHexafluoride | 300–1000 | 300–800 | −0.17 % | 355 | 97.42 |
| Toluene | 200–1000 | 200–800 | −1.33 % | 200 | 103.28 |
| Water | 200–1000 | 274–800 | +0.04 % | 744 | 33.59 |
| Xenon | 200–1000 | 200–800 | 0.00 % | — | 20.79 |

The larger differences for the heavier organics are differences between the two published correlations (TM-4513's fit and the ideal-gas part of each reference equation of state), not fitting noise; below its triple point (279 K) cyclohexane's two correlations differ by up to 8 %, which is why the comparison starts there.

**Every component has a polynomial.** Earlier engines also had Hexane, which has no TM-4513 entry. Its heat capacity was CoolProp's value at the initial conditions, held for the run.

!!! note "Changed in 0.10"
    Hexane is removed. A deck that has Hexane does not load. The error names the component.

**Kirchhoff heats of reaction** integrate the same heat capacities exactly: $\Delta H_r(T) = \Delta H_r(T_\mathrm{ref}) + \sum_i \nu_i \int_{T_\mathrm{ref}}^{T} C_{p,i}\thinspace dT$.

### Constant gas heat capacity and viscosity

Models in the literature often use fixed gas properties; for example, Yancy-Caballero et al. (2020) and Haghpanah et al. (2013) use $C_{p,g} = 30.7$ J mol⁻¹ K⁻¹ and $\mu = 1.72\times10^{-5}$ Pa s. Two fluid switches reproduce them:

- **Gas heat capacity mode** (`fluid.gas_heat_capacity_mode`): `ideal_gas` (the default, as above; `coolprop`, its former name, is still accepted) or `constant`. With `constant`, the molar heat capacity of each component is the value given in `fluid.gas_heat_capacity[<component>]` [J mol⁻¹ K⁻¹], and $C_{p,i}$ is that value divided by $M_i$. The mixing rule below is unchanged, so equal values give that molar heat capacity at every composition. A new component starts at its ideal-gas heat capacity at 25 °C (the polynomial above at 298.15 K, rounded to 0.01 J mol⁻¹ K⁻¹), for example 37.14 for CO₂ and 29.12 for N₂. You can change the value.
- **Viscosity mode** (`fluid.viscosity_mode`): `polynomial` (the default, the mole-fraction average below) or `constant`. With `constant`, every cell has the mixture viscosity `fluid.viscosity` [Pa s], whatever its temperature and composition.

The *Gas phase* adsorbed-phase heat capacity takes the same $C_{p,i}\thinspace M_i$, at the temperature of the adsorbed phase (the gas temperature in LTE, the solid's in LTNE), so with a constant gas heat capacity it is the constant too.

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

### Binary diffusion coefficients (Chapman–Enskog)

The binary molecular diffusion coefficient of every pair of components is evaluated in every cell at every residual evaluation, at the cell's gas temperature and pressure, from the Chapman–Enskog relation in the form of Poling, Prausnitz and O'Connell (*The Properties of Gases and Liquids*, 5th ed., 2001, Chapter 11) (`core/src/equations/binary_diffusion_model.cpp`):

$$
D_{ij} = \frac{2.66\times10^{-7}\thinspace T^{3/2}}{\left(P/10^{5}\right)\thinspace M_{ij}^{1/2}\thinspace \sigma_{ij}^{2}\thinspace \Omega_D},
\qquad
M_{ij} = \frac{2}{1/M_i + 1/M_j}
$$

with $D_{ij}$ in m² s⁻¹, $T$ in K, $P$ in Pa (so $P/10^5$ is in bar), the molar masses $M$ in g mol⁻¹ and $\sigma_{ij}$ in Å. The collision integral is the correlation of Neufeld, Janzen and Aziz (1972) in the reduced temperature $T^* = T/\varepsilon_{ij}$:

$$
\Omega_D = \frac{A}{(T^*)^{B}} + \frac{C}{\exp(D\thinspace T^*)} + \frac{E}{\exp(F\thinspace T^*)} + \frac{G}{\exp(H\thinspace T^*)}
$$

$A = 1.06036$, $B = 0.15610$, $C = 0.19300$, $D = 0.47635$, $E = 1.03587$, $F = 1.52996$, $G = 1.76474$, $H = 3.89411$.

The pair parameters follow the Lorentz–Berthelot combining rules from each component's Lennard-Jones length $\sigma_i$ [Å] and energy $\varepsilon_i$ (as $\varepsilon/k_B$, K):

$$
\sigma_{ij} = \tfrac12\left(\sigma_i + \sigma_j\right), \qquad \varepsilon_{ij} = \sqrt{\varepsilon_i\thinspace \varepsilon_j}
$$

The Lennard-Jones parameters are those of Svehla (NASA TR R-132, 1962, Table I(a)); octane and toluene use CoolProp's viscosity-fit pairs, and isopentane uses the n-pentane row.

| Symbol | Description | Unit |
|---|---|---|
| $D_{ij}$ | Binary diffusion coefficient of the pair $i$, $j$ | m² s⁻¹ |
| $M_i$ | Molar mass of component $i$ | g mol⁻¹ |
| $\sigma_i$, $\sigma_{ij}$ | Lennard-Jones length of component $i$, of the pair | Å |
| $\varepsilon_i$, $\varepsilon_{ij}$ | Lennard-Jones energy over $k_B$ of component $i$, of the pair | K |
| $T^*$ | Reduced temperature $T/\varepsilon_{ij}$ | - |
| $\Omega_D$ | Diffusion collision integral | - |

### Molecular diffusion in a mixture

Each component $i$ has one effective molecular diffusivity $D_{i,m}$ in every cell (`UpdateMolecularDiffusionCoefficient` and `BlanksLaw` in `core/src/equations/adsorbent_equations.cpp`). It depends on the number of components $N$:

- **One component.** $D_{i,m} = D_{ii}$, the self-diffusion coefficient: the relation above with $\sigma_{ii} = \sigma_i$, $\varepsilon_{ii} = \varepsilon_i$ and $M_{ii} = M_i$.
- **Two components.** Both components take the binary value, $D_{1,m} = D_{2,m} = D_{12}$.
- **Three or more components.** Blanc's law, with a regularisation in the denominator:

$$
D_{i,m} = \frac{\left\lvert \sum_{j \ne i} y_j \right\rvert}{\sqrt{\epsilon + \left(\sum_{j \ne i} y_j / D_{ij}\right)^{2}}},
\qquad \epsilon = 10^{-12}\ \text{s}^2\thinspace\text{m}^{-4}
$$

The $y_j$ are the raw mole fractions of the cell, the solver's unknowns, not normalised to a sum of one. Away from the regularisation this is the usual Blanc form $D_{i,m} = (1 - y_i)\big/\sum_{j\ne i} (y_j/D_{ij})$, the harmonic mean of the binary coefficients weighted by the other components' mole fractions. Blanc's law describes component $i$ diffusing through a stagnant mixture of the others (see Poling, Prausnitz and O'Connell, Chapter 11). It is not a full Stefan–Maxwell treatment: the diffusion fluxes of the components are not coupled to each other.

Limits:

- **Component $i$ in trace amounts** ($y_i \to 0$, $\sum_{j\ne i} y_j \to 1$): $D_{i,m}$ tends to the mole-fraction-weighted harmonic mean of $D_{ij}$, a finite value.
- **Component $i$ nearly pure** ($\sum_{j \ne i} y_j \to 0$): Blanc's ratio itself stays finite, but the regularisation takes over once $\sum_{j\ne i} y_j / D_{ij}$ falls to about $\sqrt{\epsilon} = 10^{-6}$ s m⁻² (about $10^{-11}$ in mole fraction for $D_{ij} \approx 10^{-5}$ m² s⁻¹). Below that, $D_{i,m} \approx 10^{6}\sum_{j\ne i} y_j$ m² s⁻¹, which goes to zero with the other components' mole fractions.

| Symbol | Description | Unit |
|---|---|---|
| $D_{i,m}$ | Effective molecular diffusivity of component $i$ in the mixture | m² s⁻¹ |
| $y_j$ | Mole fraction of component $j$ (unnormalised) | - |
| $\epsilon$ | Regularisation of the Blanc denominator | s² m⁻⁴ |

**Mixture diffusivity in the Schmidt number.** The Schmidt number uses the mole-fraction-weighted mean of the component diffusivities, with a floor:

$$
\mathrm{Sc} = \frac{\mu}{\rho\thinspace D_{mix}}, \qquad D_{mix} = \max\left(\sum_i y_i\thinspace D_{i,m},\ 10^{-9}\ \text{m}^2\thinspace\text{s}^{-1}\right)
$$

**Where each diffusivity is used.**

| Model term | Diffusivity |
|---|---|
| Axial dispersion, Ruthven and Taylor–Aris ([Transport Correlations](transport-correlations.md)) | the **user-entered constant** $D_m$ of the dispersion model; the computed $D_{i,m}$ is not used |
| Film mass-transfer coefficient, $k_f = \mathrm{Sh}\thinspace D_{i,m}/d$ (the film term of the kinetic models) | $D_{i,m}$ |
| Haynes (1973) mass-transfer coefficient (film and macropore terms) | $D_{i,m}$ |
| Macropore molecular-diffusion term (Macropore, Macropore + Knudsen, Macropore capacity factor, Mass-transfer resistances) | the model's constant molecular diffusivity when it is set above 0; otherwise $D_{i,m}$ |
| Schmidt number, and through it the Sherwood correlations | $D_{mix}$ above |
| Thiele modulus of a reaction, $D_{eff} = (\varepsilon_p/\tau)\thinspace D_{i,m}$ of the key reactant ([Reactions](reactions.md)) | $D_{i,m}$ |

See [Kinetic Models](kinetic-models.md) for the mass-transfer coefficients.

### Other gas properties

- **Molecular diffusivities** come from the Chapman–Enskog binary-diffusion model and the mixture rules above.
- **Adsorbed-phase heat capacity** is either a constant per component or, with the *Gas phase* option of the energy balance, $C_{p,i}\thinspace M_i$ from the same gas heat capacity at the adsorbed phase's temperature; see [Energy Balance](energy-balance.md#adsorbed-phase-heat-capacity).
- **Isothermal runs** evaluate all of the above at the fixed temperature $T_{\text{iso}}$; see [Isothermal mode](energy-balance.md#isothermal-mode).
