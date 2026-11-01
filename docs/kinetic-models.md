## Kinetic Rate Models

The **isotherm** sets the equilibrium loading $q_i^*$; the **rate model** governs how fast the actual adsorbed loading $q_i$ moves towards that equilibrium. Together they form the adsorption source term that appears in the [mass](mass-balance.md) and [energy](energy-balance.md) balances.

Each component of each adsorbent has a **rate model**. The rate model exposes a **Mass Transfer Coefficient Type** which determines how the coefficient $k$ is obtained — either entered directly or computed from a diffusion correlation.

---

## Rate Models

### Linear Driving Force (1st Order)

The first-order linear driving force (LDF) model assumes the uptake rate is proportional to the difference between the equilibrium and current loading. This is the standard choice for the vast majority of PSA/TSA simulations.

$$
\frac{\partial q_i}{\partial t} = k_i \left(q_i^* - q_i\right)
$$

### Linear Driving Force (2nd Order)

A quadratic driving force, occasionally used for strongly non-linear uptake.

$$
\frac{\partial q_i}{\partial t} = k_i \left(q_i^* - q_i\right)^2
$$

| Symbol | Description | Unit |
|--------|-------------|------|
| $k_i$ | Mass transfer coefficient | s⁻¹ |
| $q_i^*$ | Equilibrium adsorbed loading (from the isotherm) | mol kg⁻¹ |
| $q_i$ | Current adsorbed loading | mol kg⁻¹ |

---

## Mass Transfer Coefficient

The coefficient $k_i$ in the LDF models can be supplied in several ways. The rate model's **mass transfer coefficient type** selects one.

### Constant

The coefficient is entered directly (as a number or an [expression](isotherms.md)).

$$
k_i = \text{constant}
$$

Input: **Mass transfer coefficient** [s⁻¹].

### Macropore Diffusion

The coefficient is derived from macropore molecular diffusion only:

$$
k_i = \frac{60 \thinspace D_{m,\text{eff}}}{d_p^2}, \qquad D_{m,\text{eff}} = D_m \frac{\epsilon_p}{\tau}
$$

Inputs: **Tortuosity** $\tau$, **Macropore Void Fraction** $\epsilon_p$, optionally a constant **Molecular diffusion coefficient** $D_m$ (see below). The particle diameter $d_p$ comes from the layer.

**Molecular diffusivity.** By default (molecular diffusion coefficient 0) $D_m$ is the binary-diffusion value the engine computes from the local temperature and pressure (Chapman–Enskog), which scales as $1/P$: at an evacuation pressure of 3–10 kPa it is 10–33 times its 1 bar value, and so is the macropore coefficient. A value above 0 [m² s⁻¹] is used as a constant $D_m$ instead, in every cell and at every pressure. Haghpanah et al. (2013) and Estupiñan Perez et al. (2019, Table 1) use a constant $D_m = 1.6 \times 10^{-5}$ m² s⁻¹. The option applies to the three macropore kinds (macropore, with the capacity factor, and macropore + Knudsen); the registry leaf is `adsorbent[...].kinetics[<component>].molecular_diffusion_coefficient` and the Python keyword `molecular_diffusion_coefficient=` of `set_rate_model`.

### Macropore Diffusion with the capacity factor

The macropore coefficient scaled by the ratio of the gas-phase to the adsorbed-phase capacity, the form used by the Haghpanah et al. (2013) VSA model and the papers built on it (e.g. Estupiñan Perez et al. 2019, SI Eq. S4):

$$
k_i = \frac{15 \thinspace \epsilon_p D_{m,i}}{\tau \thinspace r_p^2}\thinspace \frac{c_i}{\rho_p \thinspace q_i^*}
$$

with $c_i$ the gas-phase concentration [mol m⁻³], $q_i^*$ the equilibrium loading [mol kg⁻¹] and $\rho_p$ the particle density, so the ratio is dimensionless. For a strongly adsorbed component (CO₂ on 13X at 0.15 bar: $c/\rho_p q^* \approx 2 \times 10^{-3}$) the coefficient is orders of magnitude below the plain macropore value. As $c_i \to 0$ both $c_i$ and $q_i^*$ vanish and their ratio tends to the Henry limit $1/(\rho_p K_H)$; evaluating the state's own $c_i / q_i^*$ there is a $0/0$ that drifts with the solver's error. The engine therefore evaluates the ratio at $\max(c_i, c_{\text{floor}})$, re-evaluating the isotherm at that concentration (every other component at its actual value), which is finite and smooth; the ratio is capped at 1.

Inputs: **Tortuosity** $\tau$, **Macropore Void Fraction** $\epsilon_p$, **Concentration floor** $c_{\text{floor}}$ [mol m⁻³] (default 10⁻³, about 2.5 Pa of partial pressure at room temperature, inside the Henry regime of the usual isotherms), optionally a constant **Molecular diffusion coefficient** $D_m$ (as for the plain macropore coefficient; the papers above use $1.6 \times 10^{-5}$ m² s⁻¹).

### Macropore + Knudsen Diffusion

As above, but the effective diffusivity combines Knudsen and molecular contributions in series:

$$
\frac{1}{D_{\text{eff}}} = \frac{1}{D_{k,\text{eff}}} + \frac{1}{D_{m,\text{eff}}}, \qquad
D_{k,\text{eff}} = D_k \frac{\epsilon_p}{\tau}, \qquad
D_{m,\text{eff}} = D_m \frac{\epsilon_p}{\tau}
$$

with the Knudsen diffusivity evaluated from kinetic theory using the pore radius $r_p$ and molar mass $M$:

$$
D_k = \frac{2}{3}\thinspace r_p \sqrt{\frac{8 R T}{\pi M}}, \qquad k_i = \frac{60 \thinspace D_{\text{eff}}}{d_p^2}
$$

Inputs: **Pore radius** $r_p$, **Tortuosity** $\tau$, **Macropore Void Fraction** $\epsilon_p$, optionally a constant **Molecular diffusion coefficient** $D_m$ in the molecular term (0, the default, uses the value from $T$ and $P$).

### Haynes (1973)

Combines an external gas-film resistance and an internal macropore resistance in series (Haynes & Sarma, 1973):

$$
k_i = \left(\frac{r_p}{3 k_f} + \frac{\tau \thinspace r_p^2}{15 \thinspace \epsilon_p D_m}\right)^{-1}, \qquad k_f = \frac{\mathrm{Sh}\thinspace D_m}{d_p}
$$

where $r_p = d_p/2$ and the film coefficient $k_f$ uses the Sherwood number $\mathrm{Sh}$. Inputs: **Tortuosity** $\tau$, **Macropore Void Fraction** $\epsilon_p$. $D_m$ is always the value from $T$ and $P$: the Sherwood number comes from the Schmidt number of that same $D_m$, so a constant $D_m$ is not offered here.

### Arrhenius

An activation-energy form for temperature-dependent kinetics:

$$
k_i = k_0 \thinspace \exp\negthinspace \left(-\frac{E_a}{R T}\right)
$$

Inputs: **Pre-exponential factor** $k_0$ [s⁻¹], **Activation energy** $E_a$ [J mol⁻¹].

---

Once loading and its rate are defined, see the [transport correlations](transport-correlations.md) for the dispersion and heat-transfer closures, and the [mass balance](mass-balance.md) for how the adsorption source enters the conservation equations.
