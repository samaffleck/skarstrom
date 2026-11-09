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
| $d_p$ | Particle diameter (a cylinder's extrudate diameter, used as an equivalent sphere; [Kinetic Models](kinetic-models.md#particle-shape)) | m |

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

For a packed layer of cylindrical particles ([Kinetic Models](kinetic-models.md#particle-shape)) $a_v = 4(1-\epsilon_b)/d_p$; the film correlations still use $d_p$. For a structured layer $a_v$ is the channel wall area per bed volume ($4a/p^2$ for a square monolith, $\pi a/p^2$ for a circular one, $2/(g+t)$ for a laminate).

Input: **Heat transfer coefficient** $h_p$ [W m⁻² K⁻¹]. Under local thermal equilibrium (LTE, the default) this model is not used.

The **Film** model computes the coefficient from correlations, $h_p^{vol} = h\thinspace a_v$:

- packed layer: $\mathrm{Nu} = h\thinspace d_p/k_g$ from the layer's [film correlation](#packed-bed-film-correlations) (default Wakao–Kaguei–Funazkri, $2 + 1.1\thinspace \mathrm{Pr}^{1/3}\mathrm{Re}^{0.6}$), with Re on the superficial velocity and $d_p$;
- structured layer: $\mathrm{Nu} = h\thinspace D_h/k_g$, the channel's value below.

The gas–solid model is a setting of the bed, not of a layer. Switching a layer to a structured kind sets it to Film, and switching the last structured layer back to packed restores Constant.

## Packed-bed film correlations

Each packed layer has one gas–particle **film correlation** (`layer[ID].film_correlation`; Python `layer.film_correlation = "gunn"`). By the heat–mass analogy the same correlation gives the particle Nusselt number with the Prandtl number and the Sherwood number with the Schmidt number:

- $\mathrm{Nu} = h\thinspace d_p/k_g$ feeds the **Film** gas–solid heat transfer model (above), $h_p^{vol} = h\thinspace a_v$;
- $\mathrm{Sh} = k_f\thinspace d_p/D_m$ feeds the external-film resistance of the mass-transfer coefficient ([Kinetic Models](kinetic-models.md#resistances); the `haynes` preset).

In both, $\mathrm{Re} = \rho\thinspace |u|\thinspace d_p/\mu$ on the superficial velocity $u$, $\mathrm{Pr} = c_p\thinspace\mu/k_g$ and $\mathrm{Sc} = \mu/(\rho\thinspace D_m)$ with the mole-fraction-weighted molecular diffusivity. A structured layer ignores the setting: its channels use the laminar Nu = Sh [below](#film-coefficients-in-structured-layers). Gas phase only.

| Choice | Correlation ($X$ = Pr for Nu, Sc for Sh) | Limit at Re → 0 |
|---|---|---|
| `wakao_funazkri` (default) | $2 + 1.1\thinspace X^{1/3}\thinspace \mathrm{Re}^{0.6}$ | 2 |
| `gunn` | $(7 - 10\epsilon_b + 5\epsilon_b^2)(1 + 0.7\thinspace \mathrm{Re}^{0.2} X^{1/3}) + (1.33 - 2.4\epsilon_b + 1.2\epsilon_b^2)\thinspace \mathrm{Re}^{0.7} X^{1/3}$ | $7 - 10\epsilon_b + 5\epsilon_b^2$ (2 at $\epsilon_b = 1$) |
| `ranz_marshall` | $2 + 0.6\thinspace \mathrm{Re}^{1/2} X^{1/3}$ (single sphere) | 2 |

**Wakao–Funazkri (mass) and Wakao–Kaguei–Funazkri (heat).** Fitted to packed-bed data re-reduced with the axial dispersion included, which matches the engine's dispersed plug-flow model. This is the Sherwood number the engine has always used for the Haynes coefficient, so decks keep their results bit for bit. (The cube root is evaluated as $X^{0.3333333}$, as before; it differs from $X^{1/3}$ by a relative $3.3\times10^{-8}\ln X$.) Validity: roughly $3 < \mathrm{Re} < 10\,000$ for mass and $15 < \mathrm{Re} < 8\,500$ for heat (*to be confirmed* against the papers).

**Gunn (1978).** Covers fixed and fluidised beds over the whole voidage range, and reduces to an isolated sphere's $\mathrm{Nu} = 2$ at $\epsilon_b = 1$, Re = 0. The void fraction is the layer's bed void fraction $\epsilon_b$. Re is on the superficial velocity and $d_p$, the basis of the engine's Re: Gunn's correlation is usually quoted with the superficial (empty-tube) slip velocity, which in CFD–DEM codes appears as $\epsilon\thinspace|u_f - u_p|$ (the velocity basis is *to be confirmed* against the paper). Validity: $0.35 \le \epsilon_b \le 1$, Re up to about $10^5$ (*to be confirmed*). At low Re in a packed bed ($\epsilon_b \approx 0.4$) its limit is $7 - 4 + 0.8 = 3.8$, above Wakao's 2.

**Ranz–Marshall (1952).** A correlation for a **single sphere** (evaporating drops), not for a packed bed: in a bed at moderate Re it under-predicts the coefficients. It is offered so that papers that use it can be reproduced, and the run log carries a warning for every packed layer that selects it (a warning, not an error). Re on the superficial velocity, as the papers that apply it to beds usually do.

Not included: older j-factor fits reduced without axial dispersion, which are not consistent with a dispersed plug-flow model.

**Sources.** Equation and table numbers are *to be confirmed* against the texts; the DOIs were checked against Crossref.

- N. Wakao and T. Funazkri, "Effect of fluid dispersion coefficients on particle-to-fluid mass transfer coefficients in packed beds", Chem. Eng. Sci. 33 (1978) 1375–1384, doi:10.1016/0009-2509(78)85120-3.
- N. Wakao, S. Kaguei and T. Funazkri, "Effect of fluid dispersion coefficients on particle-to-fluid heat transfer coefficients in packed beds", Chem. Eng. Sci. 34 (1979) 325–336, doi:10.1016/0009-2509(79)85064-2.
- N. Wakao and S. Kaguei, *Heat and Mass Transfer in Packed Beds*, Gordon and Breach, 1982 (both correlations).
- D. J. Gunn, "Transfer of heat or mass to particles in fixed and fluidised beds", Int. J. Heat Mass Transfer 21 (1978) 467–476, doi:10.1016/0017-9310(78)90080-7.
- W. E. Ranz and W. R. Marshall, "Evaporation from drops", Chem. Eng. Prog. 48 (1952) 141–146 and 173–180 (no DOI).

The verification tests `FilmCorrelation.*` check each correlation's analytic limits, hand-computed values, the heat–mass symmetry and $h \propto 1/d_p$, and that the Film heat model and the Haynes coefficient both follow the layer's choice ([Verification](verification.md)).

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
