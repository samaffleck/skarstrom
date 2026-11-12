## Kinetic Rate Models

<!-- doc-setup
import skarstrom
sim = skarstrom.examples.Perez2019_A.load()
ads = sim.adsorbents["13X"]
-->

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

The same symbols are used throughout this page: $k_i$ is the coefficient of the [rate model](#rate-models), whatever [resistances](#mass-transfer-coefficient) produce it; the dual-rate model below writes its two coefficients $k_{i,1}$ (fast) and $k_{i,2}$ (slow), and each comes from its own coefficient model.


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
| coefficient parameters (e.g. `k`) | Fast-site coefficient $k_{i,1}$ | as the [resistance](#resistances) | |
| `slow_<name>` (e.g. `slow_k`) | Slow-site coefficient $k_{i,2}$ | as the resistance | the fast site's values |

The slow site has its own coefficient model, with its own [switches](#setting-the-switches), capacity-factor mode and values, because published fits report both rate constants. With the lumped term alone, the registry paths are `adsorbent[A].kinetics[C].k` ($k_{i,1}$), `adsorbent[A].kinetics[C].slow_k` ($k_{i,2}$) and `adsorbent[A].kinetics[C].split` ($\varphi$). A preset (`set_rate_model`, `SetMassTransferKind`) applies to both sites. A new dual-rate model starts with the slow site a copy of the fast one, so it behaves like the LDF until you change `split` or the slow coefficient.

With $\varphi = 1$ the slow site holds nothing. The model is then exactly the first-order LDF with $k = k_{i,1}$, and the solver adds no unknown for it. The slow-site unknown exists only for components whose $\varphi$ is below 1. A warm start from a state taken with another $\varphi$, or with another rate model, gives a new slow site its share $(1 - \varphi)$ of the total loading. A slow site the deck no longer has is dropped, and the total loading keeps the inventory. $\varphi$ is one number per run; it may come from a variable but not from the local state.

```python
ads.set_rate_model("CO2", "dual_rate_ldf", "constant", k=0.05, slow_k=0.002, split=0.6)
```

---

## Mass Transfer Coefficient

The coefficient $k_i$ of the LDF models comes from one model: **mass-transfer resistances in series** (tag `series_resistance`). It has four main resistances, and two switches for the pore diffusivity of the macropore term. You switch each one on or off by itself, and the enabled terms add up:

$$
\frac{1}{k_i} = \underbrace{f_\text{film}\thinspace K\thinspace R_\text{film}}_\text{film}
 + \underbrace{K\thinspace R_\text{macro}}_\text{macropore (molecular and/or Knudsen)}
 + \underbrace{R_\text{micro}}_\text{micropore}
 + \underbrace{\frac{1}{k_{L,i}}}_\text{lumped}
$$

A term that is switched off adds nothing. $K$ is the capacity factor: it multiplies the gas-side terms (film and macropore) only. The [capacity-factor mode](#capacity-factor-modes) sets how the engine gets $K$. This is the series (Glueckauf) form of the LDF coefficient (Glueckauf, 1955; Ruthven, 1984). Stampi-Bombelli et al. (2024, *Ind. Eng. Chem. Res.* 63, 11637, doi:[10.1021/acs.iecr.4c01392](https://doi.org/10.1021/acs.iecr.4c01392)) use it for direct-air-capture breakthrough curves.

The coefficient kinds of earlier versions (`constant`, `arrhenius`, `macropore`, `macropore_knudsen`, `haynes`, `macropore_capacity`) are now [presets](#presets) of this model: a set of switches, a mode and a few fixed values.

!!! note "New in 0.9"
    The switch of each resistance and the capacity-factor mode are registry parameters: `adsorbent[A].kinetics[C].resistances.<tag>` and `adsorbent[A].kinetics[C].capacity_factor_mode`. The deck view also shows the [preset](#presets) the switches match and the [display equation](#display-equation). In 0.8 these were set only with the deck edits `SetMassTransferResistance` and `SetCapacityFactorMode`.

### Resistances

The geometry comes from the layer. Packed layers use spheres of radius $r_p = d_p/2$, or infinite cylinders of radius $r_f = d_p/2$ when the layer's [particle shape](#particle-shape) is `cylinder`. Structured layers ([Momentum Balance](momentum-balance.md#structured-layers-laminar-flow-in-channels)) use a slab of half-thickness $\ell$: half the monolith wall, half the laminate sheet, or the washcoat thickness when there is a washcoat. The factors 15 (sphere, on the radius), 8 (cylinder, on the radius) and 3 (slab, on the half-thickness) are Glueckauf's (1955) linear-driving-force approximation for each shape. For the slab terms, refer to Rezaei & Webley (2009, *Chem. Eng. Sci.* 64, Eq. 6 on p. 5184 and Table 1 on p. 5183).

| Switch (`resistances.<tag>`) | Term of $1/k_i$, sphere | Term of $1/k_i$, cylinder | Term of $1/k_i$, slab | Definitions |
| --- | --- | --- | --- | --- |
| `film` (external film) | $f_\text{film}\thinspace\dfrac{K\thinspace r_p}{3\thinspace k_f}$ | $f_\text{film}\thinspace\dfrac{K\thinspace r_f}{2\thinspace k_f}$ | $f_\text{film}\thinspace\dfrac{K\thinspace\ell}{k_f}$ | $k_f = \mathrm{Sh}\thinspace D_{m,i}/d_p$ (sphere, cylinder), $\mathrm{Sh}\thinspace D_{m,i}/D_h$ (slab) |
| `macropore` (macropore diffusion) | $\dfrac{K\thinspace\tau\thinspace r_p^2}{15\thinspace\varepsilon_p D}$ | $\dfrac{K\thinspace\tau\thinspace r_f^2}{8\thinspace\varepsilon_p D}$ | $\dfrac{K\thinspace\tau\thinspace\ell^2}{3\thinspace\varepsilon_p D}$ | $D = D_{m,i}$ (molecular only), $D_{K,i}$ (Knudsen only) or $D_{B,i}$ (both) |
| `molecular` (molecular diffusion; only with `macropore` on) | refines the macropore term | refines the macropore term | refines the macropore term | $D_{m,i}$: the local binary diffusivity, or the constant `molecular_diffusion_coefficient` |
| `knudsen` (Knudsen diffusion; only with `macropore` on) | refines the macropore term | refines the macropore term | refines the macropore term | $D_{K,i} = \tfrac{2}{3}\thinspace r_\text{pore}\sqrt{\dfrac{8RT}{\pi M_i}}$; with `molecular` on too, $\dfrac{1}{D_{B,i}} = \dfrac{1}{D_{K,i}} + \dfrac{1}{D_{m,i}}$ |
| `micropore` (micropore diffusion) | $\dfrac{r_c^2}{15\thinspace D_{c,i}}$ | the same | the same | $D_{c,i} = D_{c0}\thinspace e^{-E_{a,c}/RT}$ |
| `lumped` (lumped coefficient) | $\dfrac{1}{k_{L,i}}$ | the same | the same | $k_{L,i} = k\thinspace e^{-E_a/RT}$; $k_{L,i} = k$ when $E_a = 0$ |

!!! note "New in 0.12"
    The switch `molecular`. Before 0.12 the macropore term always held molecular diffusion; a deck from an earlier version loads with `molecular` on wherever `macropore` is on, so its results do not change. With `molecular` off and `knudsen` on, the pore diffusivity is $D_{K,i}$ alone.

The film coefficient $k_f$ uses the Sherwood number $\mathrm{Sh}$ of the layer's film correlation (default Wakao and Funazkri; [Transport Correlations](transport-correlations.md#packed-bed-film-correlations)). For a structured layer it uses the channel Sh and the hydraulic diameter $D_h$ ([Transport Correlations](transport-correlations.md#film-coefficients-in-structured-layers)). The micropore term is that of a spherical crystal in all geometries.

#### Particle shape

!!! note "New in 0.9"
    The layer leaf `layer[L].particle_shape` (`sphere`, `cylinder`), the deck edit `SetParticleShape` and the deck-view layer field `particle_shape`.

A packed layer holds spheres (`sphere`, the default and the shape of every deck from an earlier version) or extrudates and fibres (`cylinder`). The layer keeps one size field, `particle_diameter`. For a cylinder it is the extrudate (fibre) diameter $d_p = 2 r_f$. The engine takes the cylinder as infinitely long: the end faces are not in the diffusion geometry. This is a good approximation for extrudates that are several diameters long.

Only the diffusion geometry of the rate coefficient follows the shape:

| | Sphere (radius $r_p$) | Infinite cylinder (radius $r_f$) | Slab (half-thickness $\ell$) |
| --- | --- | --- | --- |
| LDF constant, $k = c\thinspace D/R^2$ | 15 | 8 | 3 |
| Volume-to-surface ratio $V/S$ (film term) | $r_p/3$ | $r_f/2$ | $\ell$ |

The constant comes from Glueckauf's (1955) assumption of a parabolic profile in the particle, $q(r) = a + b\thinspace r^2$. For a cylinder the mean is $\bar q = a + b R^2/2$ and the flux through the surface gives $d\bar q/dt = (2/R)\thinspace D\thinspace(\partial q/\partial r)_R = 4 D b$. With $q_s - \bar q = b R^2/2$, this is $d\bar q/dt = (8 D/R^2)(q_s - \bar q)$. The same steps for a sphere ($\bar q = a + 3bR^2/5$) give the factor 15, and for a slab give 3.

The pressure drop (Ergun, Kozeny–Carman, Blake–Kozeny), the axial dispersion, and the packed-bed film correlations (Sh and Nu) continue to use $d_p$ as an equivalent-sphere diameter, for the Reynolds and Sherwood length too. The specific surface area follows the shape: $a_v = 6(1-\varepsilon_b)/d_p$ for spheres and $a_v = 4(1-\varepsilon_b)/d_p$ for infinite cylinders of diameter $d_p$ (end faces not counted). The gas–solid heat transfer ($h\thinspace a_v$) and the packed hydraulic diameter $D_h = 4\varepsilon_b/a_v$ use it. A shape-specific pressure-drop correlation is not in the engine. A structured layer ignores the particle shape.

**Molecular diffusivity.** By default (`molecular_diffusion_coefficient` = 0), $D_{m,i}$ is the binary-diffusion value that the engine calculates from the local temperature and pressure (Chapman–Enskog). This value changes as $1/P$: at an evacuation pressure of 3–10 kPa it is 10–33 times its 1 bar value, and so is the macropore coefficient. A value above 0 [m² s⁻¹] is a constant $D_m$ for the macropore term, in each cell and at each pressure. Haghpanah et al. (2013) and Estupiñan Perez et al. (2019, Table 1) use $D_m = 1.6 \times 10^{-5}$ m² s⁻¹. The film term always uses the local value, because its Sherwood number comes from the Schmidt number of that same $D_m$.

The parameters are registry leaves with flat names, `adsorbent[A].kinetics[C].<name>`. The registry shows a leaf only while its resistance is on. The model keeps every value whatever its switch, so if you switch a term off and on again, its values come back.

| Leaf | Resistance | Symbol | Unit | Default |
| --- | --- | --- | --- | --- |
| `film_resistance_factor` | film | $f_\text{film}$ | – | 1 |
| `tortuosity` | macropore | $\tau$ | – | 3 |
| `macropore_void_fraction` | macropore | $\varepsilon_p$ | – | 0.35 |
| `molecular_diffusion_coefficient` | molecular | constant $D_m$ (0: local value) | m² s⁻¹ | 0 |
| `pore_radius` | knudsen | $r_\text{pore}$ | m | 10⁻⁹ |
| `crystal_radius` | micropore | $r_c$ | m | 10⁻⁶ |
| `micropore_diffusivity` | micropore | $D_{c0}$ | m² s⁻¹ | 10⁻¹³ |
| `micropore_activation_energy` | micropore | $E_{a,c}$ | J mol⁻¹ | 0 |
| `k` | lumped | $k$ (LDF2: kg mol⁻¹ s⁻¹) | s⁻¹ | 0.5 |
| `activation_energy` | lumped | $E_a$ | J mol⁻¹ | 0 |
| `concentration_floor` | modes slope, secant | $c_\text{floor}$ | mol m⁻³ | 10⁻³ |
| `capacity_factor` | mode fixed | $K_\text{fix}$ | – | 0 (choosing the fixed mode sets 0 to 1) |

At least one of film, macropore, micropore and lumped must stay on. The engine refuses a switch that turns the last one off. It also refuses molecular or Knudsen without macropore, and it refuses to switch off the last of molecular and Knudsen while macropore is on. Switching macropore on with neither on also switches molecular on. Switching macropore off also switches molecular and Knudsen off. When you switch a term on for the first time, a parameter of that term that is out of range (for example a zero from an old deck) gets the default.

Sources of the terms:

- **Film:** Haynes & Sarma (1973); Ruthven (1984).
- **Macropore:** Glueckauf (1955); Ruthven (1984). With the secant capacity factor, Haghpanah et al. (2013, *Ind. Eng. Chem. Res.* 52, p. 4253, Eqs. 40–41).
- **Knudsen:** Knudsen diffusion in the macropores, alone or in the Bosanquet combination with molecular diffusion; Ruthven (1984), Do (1998).
- **Micropore:** Kärger & Ruthven (1992); Farooq & Ruthven (*Chem. Eng. Sci.* 46, 1991, 2213).
- **Lumped:** a fitted LDF coefficient (Glueckauf, 1955), with an optional Arrhenius temperature dependence.
- **Slab forms:** Rezaei & Webley (2009, p. 5184, Eq. 6; Table 1, p. 5183).

The books are cited without page numbers.

### Capacity-factor modes

$K$ [-] is the ratio of the solid's capacity to the gas's capacity for the component. The film and macropore resistances are then the resistances of a solid-phase driving force. $q_i^*$ is in mol kg⁻¹, $c_i = y_i P/(RT)$ is in mol per m³ of gas, and $\rho_p$ is the particle density (`solid_density`). The mode is the leaf `capacity_factor_mode` (`none`, `slope`, `secant`, `fixed`). The registry shows this leaf only while film or macropore is on.

| Mode | $K$ | Notes |
| --- | --- | --- |
| `none` | $1$ | The gas-side terms as written. |
| `slope` | $\max\!\left(\rho_p\thinspace\dfrac{\partial q_i^*}{\partial c_i},\thinspace 1\right)$ | The local value from the Glueckauf derivation (Farooq & Ruthven, 1991). The engine uses a forward difference of the isotherm in $c_i$ (relative step 10⁻⁶). Where the slope is not positive, the engine uses the secant. |
| `secant` | $\max\!\left(\dfrac{\rho_p\thinspace q_i^*}{c_i},\thinspace 1\right)$ | The form of Haghpanah et al. (2013) and the papers that use their model (for example Estupiñan Perez et al. 2019, SI Eq. S4). |
| `fixed` | $K_\text{fix}$ | For example a published Henry constant. The engine does not floor this value. |

The slope and secant modes evaluate the isotherm at $\max(c_i, c_\text{floor})$, with each other component at its actual value (under the IAST mixing rule, the IAST mixture loading). As $c_i \to 0$, both $c_i$ and $q_i^*$ go to zero, and their ratio goes to the Henry limit. If the engine evaluated the state's own $c_i/q_i^*$ there, it would get a $0/0$ that drifts with the solver's error. The floor (default 10⁻³ mol m⁻³, approximately 2.5 Pa of partial pressure at room temperature) keeps the value finite and smooth. With the floor at 1, $K$ cannot make the coefficient larger than the $K = 1$ value. For a strongly adsorbed component (CO₂ on 13X at 0.15 bar, $\rho_p q^*/c \approx 500$), the coefficient is orders of magnitude below the value without the capacity factor.

### Display equation

!!! note "New in 0.9"
    The deck view's kinetics entry has `equation` (fast site) and `slow_equation` (slow site of a dual-rate model). Each is the LaTeX of the current coefficient model. The engine makes it only from the fragments of the kind descriptor (`MassTransferEquationLatex`), and it makes it again after each edit.

The equation is always in the series form, also for a single term. The definitions of the enabled terms and of the mode follow it, one definition on each line, aligned on the equals sign. With the film term on in a packed layer, the lines after $k_f$ are the Sherwood number of the layer's [film correlation](transport-correlations.md#packed-bed-film-correlations) and the Reynolds and Schmidt numbers it uses. For example, the `haynes` preset on packed spheres with the default Wakao–Funazkri correlation gives

$$
\begin{aligned}\frac{1}{k_i} &= f_{film}\,\frac{r_p}{3\,k_f} + \frac{\tau\,r_p^{2}}{15\,\varepsilon_p\,D_{m,i}} \\ k_f &= \frac{Sh\,D_{m,i}}{d_p} \\ Sh &= 2 + 1.1\,Sc^{1/3}\,Re^{0.6} \\ Re &= \frac{\rho_g\,u\,d_p}{\mu} \\ Sc &= \frac{\mu}{\rho_g\,D_{m,i}}\end{aligned}
$$

and `macropore_capacity` gives

$$
\begin{aligned}\frac{1}{k_i} &= \frac{K\,\tau\,r_p^{2}}{15\,\varepsilon_p\,D_{m,i}} \\ K &= \max\!\left(\frac{\rho_p\,q_i^{*}}{c_i},\,1\right)\end{aligned}
$$

If the adsorbent is in a packed layer, the deck view uses the fragments and the film correlation of that layer (sphere or cylinder; with several packed layers, the first that holds the adsorbent). If it is only in structured layers, it uses the slab fragments. If it is in layers of the two kinds, `equation` is the packed form and `equation_slab` is the slab form.

### Presets

<a id="constant"></a><a id="arrhenius"></a><a id="macropore-diffusion"></a><a id="macropore-diffusion-with-the-capacity-factor"></a><a id="macropore-knudsen-diffusion"></a><a id="haynes-1973"></a><a id="series-resistance-film-macropore-micropore"></a>

A preset is a name for a set of switches, a mode and some fixed values. The former coefficient kinds are presets. Each deck-view kinetics entry has the `preset` that its current switches match, or `custom`. The engine calculates the preset from the model each time, so it cannot be out of date: if you switch film on in a `macropore` model, it becomes `custom`, and if you switch film off again, it becomes `macropore` again.

| Preset (former kind) | Switches | Mode | Fixed values | Parameters it shows |
| --- | --- | --- | --- | --- |
| `constant` | lumped | – | $E_a = 0$ | `k` |
| `arrhenius` | lumped | – | ($E_a \ne 0$) | `k` ($k_0$; the former `pre_exp_factor`), `activation_energy` |
| `macropore` | macropore, molecular | none | | `tortuosity`, `macropore_void_fraction`, `molecular_diffusion_coefficient` |
| `macropore_knudsen` | macropore, molecular, knudsen | none | | the macropore leaves and `pore_radius` |
| `macropore_capacity` | macropore, molecular | secant | | the macropore leaves and `concentration_floor` |
| `haynes` | film, macropore, molecular | none | $f_\text{film} = 1$, constant $D_m$ = 0 | `film_resistance_factor`, the macropore leaves |
| `custom` | any other set | any | | the leaves of each enabled term and of the mode |

The single-term presets reduce to the forms of the former kinds:

- **constant:** $k_i = k$.
- **arrhenius:** $k_i = k_0\thinspace e^{-E_a/RT}$.
- **macropore:** $k_i = 60\thinspace D_\text{eff}/d_p^2$ with $D_\text{eff} = D\thinspace\varepsilon_p/\tau$ (cylinder: $8\thinspace D_\text{eff}/r_f^2$; slab: $3\thinspace D_\text{eff}/\ell^2$).
- **macropore_knudsen:** the same, with $D = D_{B,i}$.
- **macropore_capacity:** $k_i = \dfrac{15\thinspace\varepsilon_p D_{m,i}}{\tau\thinspace r_p^2}\thinspace\min\!\left(\dfrac{c_i}{\rho_p q_i^*},\thinspace 1\right)$.
- **haynes:** $k_i = \left(\dfrac{r_p}{3 k_f} + \dfrac{\tau\thinspace r_p^2}{15\thinspace\varepsilon_p D_{m,i}}\right)^{-1}$ (Haynes & Sarma, 1973).

The engine evaluates each preset with the operations of the former class, so the coefficient is the same double.

A preset tag is not a kind of its own. The builder, `set_rate_model` and `SetMassTransferKind` accept it as a shortcut: it sets the switches, the mode and the fixed values. The deck keeps only the switches and the values.

### Setting the switches

```python
ads.set_rate_model("CO2", "ldf", "haynes", tortuosity=3, macropore_void_fraction=0.35)   # a preset
kin = sim.adsorbents["13X"].kinetics["CO2"]
path = kin.resistances.path("micropore")   # adsorbent[A].kinetics[CarbonDioxide].resistances.micropore
sim[path] = True                           # switch on by path (or kin.resistances.micropore = True)
kin = sim.adsorbents["13X"].kinetics["CO2"]
kin.capacity_factor_mode = "slope"
kin.crystal_radius = 1e-6                  # listed now that micropore is on
print(path, kin.crystal_radius)
```

A custom set in one call:

```python
ads.set_rate_model("CO2", "ldf", "series_resistance", resistances=["film", "macropore", "micropore"],
                   capacity_factor_mode="slope", tortuosity=3, macropore_void_fraction=0.35,
                   crystal_radius=1e-6, micropore_diffusivity=2e-12, micropore_activation_energy=2e4)
ads.set_rate_resistance("CO2", "film", False)   # one switch; site="slow" for the slow site of a dual-rate model
```

In `resistances=`, `macropore` keeps molecular diffusion on (the form of earlier versions). For Knudsen diffusion alone, switch `molecular` off after the call:

```python
ads.set_rate_model("CO2", "ldf", "series_resistance", resistances=["macropore", "knudsen"], pore_radius=2e-9)
ads.set_rate_resistance("CO2", "molecular", False)   # D = D_K alone
```

!!! note "New in 0.9"
    `set_rate_model` accepts `resistances=` and `capacity_factor_mode=`, and `Adsorbent.set_rate_resistance` sets one switch. A model made from the tag `series_resistance` starts with the lumped term only ($k = 0.5$ s⁻¹), until you set its switches.

When you set a switch or the mode by path, the engine uses the deck edits `SetMassTransferResistance` and `SetCapacityFactorMode`, with their rules. A dual-rate model's slow site has its own switches, mode and parameters, with the prefix `slow_`: `kinetics[C].slow_resistances.<tag>`, `kinetics[C].slow_capacity_factor_mode` and `kinetics[C].slow_<name>`. The registry shows these leaves only while the model has a slow site.

The former path `kinetics[C].pre_exp_factor` (and `slow_pre_exp_factor`) of the `arrhenius` kind still resolves: it is the lumped `k`.

### Migration

Decks saved by earlier versions load unchanged, and their results do not change. A deck that stores one of the former coefficient kinds loads as the generic model with that preset's switches. A version-1 `series_resistance` model, in which a term was off when its parameter was zero, loads with the switches that those zeros meant. The deck is saved again in the generic form (cereal version 2).

---

## Which to use for direct air capture

Direct-air-capture sorbents (amine-functionalised resins, aluminas and silicas) at 400 ppm CO₂ sit on the steep, Henry-like part of their isotherm, where the capacity factor $K$ is large, and their breakthrough curves show a fast initial front followed by a slow tail. Two choices follow:

- **Rate model.** Use [dual-rate LDF](#dual-rate-linear-driving-force) when the measured breakthrough or uptake curve has the tail, with $\varphi$ the fraction of the capacity that loads fast and $k_{i,1} \gg k_{i,2}$. Stampi-Bombelli et al. (2024, [doi:10.1021/acs.iecr.4c01392](https://doi.org/10.1021/acs.iecr.4c01392)) fit breakthrough curves of amine-functionalised sorbents at DAC conditions with two rates. A single LDF is enough when the outlet reaches the feed in one front.
- **Coefficient.** Switch on the physical [resistances](#resistances) (film, macropore, micropore, with $K$ from the isotherm slope) when the rate constants should follow from particle and transport properties (film, macropore and micropore resistances, with $K$ from the isotherm) rather than being fitted, for instance to scale a laboratory column to another particle size or temperature. Use the lumped term alone (the `constant` [preset](#presets)) when you have fitted values of $k_{i,1}$ and $k_{i,2}$ directly. The slow site can use a different $D_{c,0}$ or $r_c$ from the fast one.

See also Stampi-Bombelli and Mazzotti (2024) on DAC process modelling. With $\varphi = 1$, the capacity factor set to 1 and the micropore term off, the two models reduce to the first-order LDF with the [`haynes` preset](#presets), bit for bit. The bundled deck `Illustrative_DAC_Breakthrough` (400 ppm CO₂ in N₂ on the Lewatit VP OC 1065 library entry) uses both together; its rate parameters are illustrative, not fitted.

---

Once loading and its rate are defined, see the [transport correlations](transport-correlations.md) for the dispersion and heat-transfer closures, and the [mass balance](mass-balance.md) for how the adsorption source enters the conservation equations.
