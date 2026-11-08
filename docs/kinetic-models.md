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


### Dual-Rate Linear Driving Force

Two first-order LDF sites in parallel that share the component's equilibrium loading (tag `dual_rate_ldf`). A fraction $\varphi$ of $q_i^*$ belongs to a fast site, the rest to a slow site, and each site relaxes towards its own share:

$$
\frac{\partial q_{i,1}}{\partial t} = k_{i,1} \left(\varphi\, q_i^* - q_{i,1}\right), \qquad
\frac{\partial q_{i,2}}{\partial t} = k_{i,2} \left((1 - \varphi)\, q_i^* - q_{i,2}\right), \qquad
q_i = q_{i,1} + q_{i,2}
$$

The component's loading $q_i$ is the sum of the two sites. It is what the mass and energy balances, the KPIs and the logged *Adsorbed Concentration* use, and the heat of adsorption applies to its total rate. The slow-site loading $q_{i,2}$ is a state variable of its own, logged as *Slow Adsorbed Concentration* (`slow_adsorbed_concentration` in `sim.series` and `sim.profile`). It starts in equilibrium with its share, $(1 - \varphi)\, q_i^*$, of the initial loading.

With a step in $q_i^*$ from a clean sorbent, the loading follows two exponentials:

$$
\frac{q_i(t)}{q_i^*} = \varphi \left(1 - e^{-k_{i,1} t}\right) + (1 - \varphi) \left(1 - e^{-k_{i,2} t}\right)
$$

This model describes uptake curves with a fast initial rise followed by a slow tail, as in the direct-air-capture breakthrough experiments on amine sorbents at 400 ppm CO₂ that Stampi-Bombelli et al. (2024) fit with two rates ([doi:10.1021/acs.iecr.4c01392](https://doi.org/10.1021/acs.iecr.4c01392)).

| Parameter (registry leaf) | Description | Unit | Default |
|---------------------------|-------------|------|---------|
| `split` | Fast-site fraction $\varphi$ of $q_i^*$, in (0, 1] | – | 1 |
| coefficient parameters (e.g. `k`) | Fast-site coefficient $k_{i,1}$ | as the coefficient type | |
| `slow_<name>` (e.g. `slow_k`) | Slow-site coefficient $k_{i,2}$ | as the coefficient type | the fast site's values |

The slow site has its own coefficient model, of the **same mass transfer coefficient type** as the fast site, with its own values, because published fits report both rate constants. With the constant type, the registry paths are `adsorbent[A].kinetics[C].k` ($k_{i,1}$), `adsorbent[A].kinetics[C].slow_k` ($k_{i,2}$) and `adsorbent[A].kinetics[C].split` ($\varphi$). Changing the coefficient type changes both sites. A new dual-rate model starts with the slow site a copy of the fast one, so it behaves like the LDF until you change `split` or the slow coefficient.

With $\varphi = 1$ the slow site holds nothing. The model is then exactly the first-order LDF with $k = k_{i,1}$, and the solver adds no unknown for it. The slow-site unknown exists only for components whose $\varphi$ is below 1. A warm start from a state taken with another $\varphi$, or with another rate model, gives a new slow site its share $(1 - \varphi)$ of the total loading. A slow site the deck no longer has is dropped, and the total loading keeps the inventory. $\varphi$ is one number per run; it may come from a variable but not from the local state.

```python
ads.set_rate_model("CO2", "dual_rate_ldf", "constant", k=0.05, slow_k=0.002, split=0.6)
```

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

**Geometry.** The factor is Glueckauf's (1955) linear-driving-force approximation for the shape of the solid, $k = c\thinspace D_\text{eff}/\ell^2$. For packed spheres, $c = 60$ on $\ell = d_p$, which is 15 on the radius. For a structured layer ([Momentum Balance](momentum-balance.md#structured-layers-laminar-flow-in-channels)) the solid is a slab: $c = 3$ on its half-thickness $\ell$, which is half the monolith wall, half the laminate sheet, or the washcoat thickness when there is a washcoat. The same applies to the capacity-factor and Knudsen variants below.

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

where $r_p = d_p/2$ and the film coefficient $k_f$ uses the Sherwood number $\mathrm{Sh}$ of the layer's film correlation (default Wakao and Funazkri; [Transport Correlations](transport-correlations.md#packed-bed-film-correlations)). For a slab (structured layer) the same series form uses the slab's $V/S = \ell$ and Glueckauf's factor 3: $k_i = \left(\ell/k_f + \tau\thinspace \ell^2/(3\thinspace \epsilon_p D_m)\right)^{-1}$, with $k_f = \mathrm{Sh}\thinspace D_m/D_h$ and the channel Sh ([Transport Correlations](transport-correlations.md#film-coefficients-in-structured-layers)). Inputs: **Tortuosity** $\tau$, **Macropore Void Fraction** $\epsilon_p$. $D_m$ is always the value from $T$ and $P$: the Sherwood number comes from the Schmidt number of that same $D_m$, so a constant $D_m$ is not offered here.

### Arrhenius

An activation-energy form for temperature-dependent kinetics:

$$
k_i = k_0 \thinspace \exp\negthinspace \left(-\frac{E_a}{R T}\right)
$$

Inputs: **Pre-exponential factor** $k_0$ [s⁻¹], **Activation energy** $E_a$ [J mol⁻¹].

---

Once loading and its rate are defined, see the [transport correlations](transport-correlations.md) for the dispersion and heat-transfer closures, and the [mass balance](mass-balance.md) for how the adsorption source enters the conservation equations.
