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
k_i = \frac{60 \, D_{m,\text{eff}}}{d_p^2}, \qquad D_{m,\text{eff}} = D_m \frac{\epsilon_p}{\tau}
$$

Inputs: **Tortuosity** $\tau$, **Macropore Void Fraction** $\epsilon_p$. The molecular diffusivity $D_m$ is computed from the gas properties, and the particle diameter $d_p$ comes from the layer.

### Macropore + Knudsen Diffusion

As above, but the effective diffusivity combines Knudsen and molecular contributions in series:

$$
\frac{1}{D_{\text{eff}}} = \frac{1}{D_{k,\text{eff}}} + \frac{1}{D_{m,\text{eff}}}, \qquad
D_{k,\text{eff}} = D_k \frac{\epsilon_p}{\tau}, \qquad
D_{m,\text{eff}} = D_m \frac{\epsilon_p}{\tau}
$$

with the Knudsen diffusivity evaluated from kinetic theory using the pore radius $r_p$ and molar mass $M$:

$$
D_k = \frac{2}{3}\, r_p \sqrt{\frac{8 R T}{\pi M}}, \qquad k_i = \frac{60 \, D_{\text{eff}}}{d_p^2}
$$

Inputs: **Pore radius** $r_p$, **Tortuosity** $\tau$, **Macropore Void Fraction** $\epsilon_p$.

### Haynes (1973)

Combines an external gas-film resistance and an internal macropore resistance in series (Haynes & Sarma, 1973):

$$
k_i = \left(\frac{r_p}{3 k_f} + \frac{\tau \, r_p^2}{15 \, \epsilon_p D_m}\right)^{-1}, \qquad k_f = \frac{\mathrm{Sh}\, D_m}{d_p}
$$

where $r_p = d_p/2$ and the film coefficient $k_f$ uses the Sherwood number $\mathrm{Sh}$. Inputs: **Tortuosity** $\tau$, **Macropore Void Fraction** $\epsilon_p$.

### Arrhenius

An activation-energy form for temperature-dependent kinetics:

$$
k_i = k_0 \, \exp\!\left(-\frac{E_a}{R T}\right)
$$

Inputs: **Pre-exponential factor** $k_0$ [s⁻¹], **Activation energy** $E_a$ [J mol⁻¹].

---

Once loading and its rate are defined, see the [transport correlations](transport-correlations.md) for the dispersion and heat-transfer closures, and the [mass balance](mass-balance.md) for how the adsorption source enters the conservation equations.
