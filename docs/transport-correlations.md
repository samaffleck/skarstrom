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
D_{ax} = 0.7\thinspace D_m + 0.5\thinspace \frac{|u|}{\epsilon_b}\thinspace d_p
$$

| Symbol | Description | Unit |
|--------|-------------|------|
| $D_m$ | Molecular diffusivity | m² s⁻¹ |
| $u$ | Superficial velocity | m s⁻¹ |
| $\epsilon_b$ | Bed void fraction | - |
| $d_p$ | Particle diameter | m |

---

### Taylor–Aris (channels)

For a structured layer (monolith or laminate; see [Momentum Balance](momentum-balance.md#structured-layers-laminar-flow-in-channels)), laminar flow in straight channels spreads a solute by the velocity profile across the channel, on top of molecular diffusion (Taylor 1953; Aris 1956):

$$
D_{ax} = D_m + \frac{u_{ch}^2\thinspace D_h^2}{\kappa\thinspace D_m}, \qquad u_{ch} = \frac{|u|}{\phi}
$$

with $\phi$ the open frontal area and $D_h$ the hydraulic diameter. $\kappa = 192$ for a circular tube, which is Taylor's $a^2 U^2/(48 D)$ on the radius $a = D_h/2$. Square channels also use 192; that is an approximation, since the exact square-duct coefficient differs by a few per cent. $\kappa = 840$ for parallel plates, which is Aris's $W^2 U^2/(210 D)$ on the gap $W$ with $D_h = 2W$. Input: **Molecular diffusivity** $D_m$ [m² s⁻¹]. The deck checks refuse this model on a packed layer, and refuse Ruthven's packed-bed model on a structured layer. Switching a layer to a structured kind (`build::SetLayerKind`) installs it.

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
\lambda = 0.5\thinspace \frac{|u|}{\epsilon_b}\thinspace d_p\thinspace \rho_g\thinspace C_{p,g}
$$

### Yagi and Kunii (1957)

A stagnant-plus-flow correlation that accounts for the solid/gas conductivity ratio:

$$
\lambda = \lambda_{z}^{0} + 0.75\thinspace k_g\thinspace \mathrm{Pr}\thinspace \mathrm{Re}
$$

$$
\lambda_{z}^{0} = k_g \left(\frac{k_s}{k_g}\right)^{n}, \qquad
n = 0.28 - 0.757 \log_{10}(\epsilon_b) - 0.057 \log_{10}\negthinspace \left(\frac{k_s}{k_g}\right)
$$

| Symbol | Description | Unit |
|--------|-------------|------|
| $\lambda_z^0$ | Stagnant (no-flow) conductivity | W m⁻¹ K⁻¹ |
| $k_g, k_s$ | Gas and solid thermal conductivity | W m⁻¹ K⁻¹ |
| $\mathrm{Pr}, \mathrm{Re}$ | Prandtl and Reynolds numbers | - |

Ruthven's and Yagi and Kunii's correlations are for packed beds: the deck checks refuse them on a structured layer.

### Channel conduction (structured layers)

In a monolith or laminate the gas and the channel walls run side by side along the column, so axial heat flows through the two in parallel, each over its share of the cross-section (the parallel arrangement of the effective conductivity of a two-phase medium, which is exact for straight channels whose walls are continuous along the flow; see e.g. Kaviany 1995):

$$
\text{LTNE:}\quad \lambda_f = \phi\thinspace k_g, \qquad \lambda_s = (1-\phi)\thinspace k_s
$$

$$
\text{LTE:}\quad \lambda = \phi\thinspace k_g + (1-\phi)\thinspace k_s
$$

with $\phi$ the open frontal area and $k_s$ the layer's effective solid conductivity (the solid-volume-weighted conductivity of its solids). Under local thermal non-equilibrium $\lambda_f$ is the gas equation's conductivity and $\lambda_s$ conducts in the [solid energy balance](energy-balance.md#solid-conduction-in-structured-layers); under local thermal equilibrium the single temperature sees their sum. Thermal Taylor dispersion in the channels is left out, as it is small at the channel Péclet numbers of adsorption contactors. The model takes no input. The deck checks refuse it on a packed layer, and switching a layer to a structured kind (`build::SetLayerKind`) installs it.

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
h_w = \frac{k_g}{d_p}\left[2\epsilon_b + \frac{1-\epsilon_b}{\tfrac{1}{3}(k_g/k_s) + \theta_w}\right] + \frac{k_g}{d_p}\left(0.0835\thinspace \mathrm{Re}^{0.91}\right), \qquad \theta_w = 0.0024\left(\frac{D_c}{d_p}\right)^{1.58}
$$

then the effective radial conductivity,

$$
k_e = k_g\left[\epsilon_b + \frac{1-\epsilon_b}{\tfrac{2}{3}(k_g/k_s) + 0.22\thinspace \epsilon_b^2}\right] + \frac{|u|\thinspace \rho_g C_{p,g}\thinspace d_p}{\mathrm{Pe}_H}, \qquad \mathrm{Pe}_H = 8.65\left[1 + 19.4\left(\frac{d_p}{D_c}\right)^2\right]
$$

and finally the apparent internal coefficient via the wall Biot number $\mathrm{Bi} = h_w D_c / (2 k_e)$,

$$
\frac{1}{h_i} = \frac{1}{h_w} + \frac{D_c}{6 k_e}\cdot\frac{\mathrm{Bi} + 3}{\mathrm{Bi} + 4}
$$

where $D_c$ is the column (internal) diameter. The film term is valid for roughly $10 < \mathrm{Re} < 1200$.

Dixon's and Leva's correlations are for packed beds (they read the particle diameter): the deck checks refuse them on a structured layer, which uses a constant $h_i$. Switching a layer to a structured kind replaces either with a constant coefficient (default 0, an adiabatic contact, so set it); a constant coefficient keeps its value. At run start the solver warns about a structured layer whose $h_i$ is 0 while the column wall exchanges heat (an outer coefficient $h_o > 0$, a fixed-temperature wall end or a heating medium), naming `layer[ID].internal_htc`. In a structured layer the whole bed-to-wall flux passes through the solid ([Wall Energy Balance](wall-energy-balance.md)).

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
h_o = \frac{\mathrm{Nu}\thinspace k_{air}}{d_e}, \qquad \mathrm{Nu} = c\thinspace \mathrm{Re}_{air}^{\thinspace m}\thinspace \mathrm{Pr}_{air}^{\thinspace n}, \qquad \mathrm{Re}_{air} = \frac{d_e\thinspace v_{ext}\thinspace \rho_{air}}{\mu_{air}}
$$

Input: **External air flow velocity** $v_{ext}$. The air properties ($k_{air}$, $\rho_{air}$, $\mu_{air}$, $c_{p,air}$) and the correlation constants $c$, $m$, $n$ are built in; $d_e$ is the external column diameter.

---

## Particle Heat Transfer (LTNE only)

When the [energy balance](energy-balance.md) uses **local thermal non-equilibrium (LTNE)**, a gas–particle coefficient is needed. The **Constant** model converts an entered coefficient into a volumetric value using the specific surface area of the packing:

$$
h_p^{vol} = h_p\thinspace a_v, \qquad a_v = \frac{6(1-\epsilon_b)}{d_p}
$$

For a structured layer $a_v$ is the channel wall area per bed volume ($4a/p^2$ for a square monolith, $\pi a/p^2$ for a circular one, $2/(g+t)$ for a laminate).

Input: **Heat transfer coefficient** $h_p$ [W m⁻² K⁻¹]. Under local thermal equilibrium (LTE, the default) this model is not used.

The **Film** model computes the coefficient from correlations, $h_p^{vol} = h\thinspace a_v$:

- packed layer: $\mathrm{Nu} = h\thinspace d_p/k_g = 2 + 1.1\thinspace \mathrm{Pr}^{1/3}\mathrm{Re}^{0.6}$ (Wakao and Funazkri 1978), with Re on the superficial velocity and $d_p$;
- structured layer: $\mathrm{Nu} = h\thinspace D_h/k_g$, the channel's value below.

The gas–solid model is a setting of the bed, not of a layer. Switching a layer to a structured kind sets it to Film, and switching the last structured layer back to packed restores Constant.

## Film coefficients in structured layers

In a laminar channel the film Nusselt and Sherwood numbers take their fully developed values at constant wall temperature, $\mathrm{Nu}_T$ (Shah and London 1978), with $\mathrm{Sh} = \mathrm{Nu}$ by the heat–mass analogy:

| Channel | $\mathrm{Nu}_T = \mathrm{Sh}$ |
|---|---|
| Square | 2.976 |
| Circular | 3.657 |
| Parallel plates (laminate) | 7.541 |

With the layer's **entrance effects** switched on (`layer[ID].entrance_effects`), Hawthorn's (1974) length-averaged developing-flow form is used over the layer length $L$:

$$
\mathrm{Nu} = \mathrm{Nu}_T\left(1 + 0.095\thinspace \mathrm{Re}\thinspace \mathrm{Pr}\thinspace \frac{D_h}{L}\right)^{0.45}
$$

with Sc in place of Pr for Sh. The Reynolds number of a structured layer is the channel value $\mathrm{Re} = \rho\thinspace (|u|/\phi)\thinspace D_h/\mu$, and the solver warns at the end of a run when it passes 2000, where these laminar closures stop holding (the setting `settings.warn_channel_reynolds`). The film mass-transfer coefficient is $k_f = \mathrm{Sh}\thinspace D_m/D_h$ (Haynes rate model, [Kinetic Models](kinetic-models.md)), and the film heat-transfer coefficient is $h = \mathrm{Nu}\thinspace k_g/D_h$ (the Film gas–solid model above).

**Sources** (the equation and table numbers are still to be confirmed against the texts; the results are as stated above):

- R. K. Shah and A. L. London, *Laminar Flow Forced Convection in Ducts*, Academic Press, 1978: fully developed $\mathrm{Nu}_T$ for square and circular ducts and parallel plates.
- R. D. Hawthorn, AIChE Symp. Ser. 70(137) (1974) 428: developing-flow Nu in monolith channels.
- G. I. Taylor, Proc. R. Soc. Lond. A 219 (1953) 186–203: dispersion in a tube.
- R. Aris, Proc. R. Soc. Lond. A 235 (1956) 67–77: the molecular term and the parallel-plate result.
- N. Wakao and T. Funazkri, Chem. Eng. Sci. 33 (1978) 1375–1384: the packed-bed film correlation.
- E. Glueckauf, Trans. Faraday Soc. 51 (1955) 1540–1551: the linear driving force, [Kinetic Models](kinetic-models.md).
- M. Kaviany, *Principles of Heat Transfer in Porous Media*, 2nd ed., Springer, 1995: the parallel arrangement of the effective thermal conductivity, the channel conduction model above (the section is to be confirmed).

The verification tests `StructuredTransfer.*` check each closure against these results, and check a parallel-passage breakthrough against the moment analysis of the LDF model; `StructuredHeat.*` check the heat closures ([Verification](verification.md)).
