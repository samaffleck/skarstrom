## Transport Correlations

Beyond the isotherm and rate model, the bed model needs closures for **axial dispersion**, **effective axial thermal conductivity**, and the **wall heat-transfer coefficients**. Each can be set to a constant value or computed from a built-in empirical correlation.

These are configured per adsorbent layer (its dispersion and heat-transfer models) and, for the external coefficient, alongside the wall settings. Choosing a correlation instead of a constant makes the coefficient respond automatically to local velocity, temperature and gas properties.

---

## Axial Mass Dispersion

Sets the axial dispersion coefficient $D_{ax}$ in the [component mass balance](mass-balance.md).

### Constant

$$
D_{ax} = \text{constant}
$$

Input: **Dispersion coefficient** [m² s⁻¹].

### Ruthven (1984)

Combines molecular diffusion and mechanical (eddy) dispersion:

$$
D_{ax} = 0.7\, D_m + 0.5\, \frac{|u|}{\epsilon_b}\, d_p
$$

| Symbol | Description | Unit |
|--------|-------------|------|
| $D_m$ | Molecular diffusivity | m² s⁻¹ |
| $u$ | Superficial velocity | m s⁻¹ |
| $\epsilon_b$ | Bed void fraction | - |
| $d_p$ | Particle diameter | m |

---

## Axial Thermal Conductivity

Sets the effective axial conductivity $\lambda$ in the [energy balance](energy-balance.md).

### Constant

$$
\lambda = \text{constant}
$$

### Ruthven (1984)

The thermal analogue of the mechanical dispersion term:

$$
\lambda = 0.5\, \frac{|u|}{\epsilon_b}\, d_p\, \rho_g\, C_{p,g}
$$

### Yagi and Kunii (1957)

A stagnant-plus-flow correlation that accounts for the solid/gas conductivity ratio:

$$
\lambda = \lambda_{z}^{0} + 0.75\, k_g\, \mathrm{Pr}\, \mathrm{Re}
$$

$$
\lambda_{z}^{0} = k_g \left(\frac{k_s}{k_g}\right)^{n}, \qquad
n = 0.28 - 0.757 \log_{10}(\epsilon_b) - 0.057 \log_{10}\!\left(\frac{k_s}{k_g}\right)
$$

| Symbol | Description | Unit |
|--------|-------------|------|
| $\lambda_z^0$ | Stagnant (no-flow) conductivity | W m⁻¹ K⁻¹ |
| $k_g, k_s$ | Gas and solid thermal conductivity | W m⁻¹ K⁻¹ |
| $\mathrm{Pr}, \mathrm{Re}$ | Prandtl and Reynolds numbers | - |

---

## Internal Heat Transfer Coefficient

The gas–wall coefficient $h_i$ in the [energy](energy-balance.md) and [wall energy](wall-energy-balance.md) balances.

### Constant

$$
h_i = \text{constant}
$$

Input: **Heat transfer coefficient** [W m⁻² K⁻¹]. This is the most common choice for fitting experimental breakthrough data.

### Dixon (1996)

An apparent wall coefficient derived from a two-dimensional pseudo-homogeneous model, folding the near-wall resistance and the effective radial bed conductivity into a single 1-D coefficient. First the wall film and stagnant contributions,

$$
h_w = \frac{k_g}{d_p}\left[2\epsilon_b + \frac{1-\epsilon_b}{\tfrac{1}{3}(k_g/k_s) + \theta_w}\right] + \frac{k_g}{d_p}\left(0.0835\,\mathrm{Re}^{0.91}\right), \qquad \theta_w = 0.0024\left(\frac{D_c}{d_p}\right)^{1.58}
$$

then the effective radial conductivity,

$$
k_e = k_g\left[\epsilon_b + \frac{1-\epsilon_b}{\tfrac{2}{3}(k_g/k_s) + 0.22\,\epsilon_b^2}\right] + \frac{|u|\,\rho_g C_{p,g}\, d_p}{\mathrm{Pe}_H}, \qquad \mathrm{Pe}_H = 8.65\left[1 + 19.4\left(\frac{d_p}{D_c}\right)^2\right]
$$

and finally the apparent internal coefficient via the wall Biot number $\mathrm{Bi} = h_w D_c / (2 k_e)$,

$$
\frac{1}{h_i} = \frac{1}{h_w} + \frac{D_c}{6 k_e}\cdot\frac{\mathrm{Bi} + 3}{\mathrm{Bi} + 4}
$$

where $D_c$ is the column (internal) diameter. The film term is valid for roughly $10 < \mathrm{Re} < 1200$.

---

## External Heat Transfer Coefficient

The wall–ambient coefficient $h_o$ in the [wall energy balance](wall-energy-balance.md). It is set with the wall settings.

### Constant

$$
h_o = \text{constant}
$$

### Holman (2008)

Forced convection over the outside of the column from a cross-flow of air:

$$
h_o = \frac{\mathrm{Nu}\, k_{air}}{d_e}, \qquad \mathrm{Nu} = c\,\mathrm{Re}_{air}^{\,m}\,\mathrm{Pr}_{air}^{\,n}, \qquad \mathrm{Re}_{air} = \frac{d_e\, v_{ext}\, \rho_{air}}{\mu_{air}}
$$

Input: **External air flow velocity** $v_{ext}$. The air properties ($k_{air}$, $\rho_{air}$, $\mu_{air}$, $c_{p,air}$) and the correlation constants $c$, $m$, $n$ are built in; $d_e$ is the external column diameter.

---

## Particle Heat Transfer (LTNE only)

When the [energy balance](energy-balance.md) uses **local thermal non-equilibrium (LTNE)**, a gas–particle coefficient is needed. The **Constant** model converts an entered coefficient into a volumetric value using the specific surface area of the packing:

$$
h_p^{vol} = h_p\, a_v, \qquad a_v = \frac{6(1-\epsilon_b)}{d_p}
$$

Input: **Heat transfer coefficient** $h_p$ [W m⁻² K⁻¹]. Under local thermal equilibrium (LTE, the default) this model is not used.
