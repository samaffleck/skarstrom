## Isotherm Models

The **isotherm** describes the equilibrium adsorbed-phase loading of a component, $q_i^*$, as a function of the local gas state (temperature, pressure and composition). It is the thermodynamic ceiling that the kinetic [rate model](kinetic-models.md) drives the actual loading $q_i$ towards.

Each component in each adsorbent material is assigned its own isotherm. In a deck, every component of every adsorbent has an **isotherm type** and the parameters of that model. How the components' isotherms combine in a mixture is the adsorbent's [mixing rule](#mixture-loadings-the-mixing-rule-and-iast): each on its own (the default), or Ideal Adsorbed Solution Theory over their pure isotherms. To plot or check an isotherm outside a simulation, evaluate it with the engine's own code from Python: `skarstrom.evaluate_isotherm(kind, params, pressure, temperature)` (see [Isotherm Fitting](isotherm-fitting.md)).

> **Partial pressure vs. concentration.** Several models come in two forms. *Partial-pressure* models are written in terms of the partial pressure $P_i = y_i P$. *Concentration* models use the molar concentration $C_i = P_i / (R T)$, computed internally from the ideal gas law. The two forms are otherwise identical; choose whichever matches how your isotherm parameters were regressed.

> **Heat of adsorption.** Each component also has an isosteric heat of adsorption $\Delta H_i$, stored with the adsorbent (parameter `adsorbent[<id>].heat_of_adsorption[<component>]`). It couples the isotherm to the [energy balance](energy-balance.md); it is a separate input and is *not* part of the isotherm equation itself. The temperature-dependent parameters in the models below (the $q$, $B_2$, $\Delta H_0$ terms) describe how loading shifts with temperature but do not by themselves supply the heat released on adsorption.

The sections below give the exact formula implemented for each model together with its parameter names: they are the names in deck files, in the Python parameter paths (`adsorbent[<id>].isotherm[<component>].<name>`) and in `skarstrom.isotherm_parameters(kind)`. $R$ is the ideal gas constant and $T$ the temperature.

---

## Inert

A non-adsorbing component. Its loading is fixed at zero, so it acts purely as a carrier gas.

$$
q_i^* = 0
$$

---

## Henry

The linear (low-coverage) limit, with a temperature-dependent Henry constant.

$$
q_i^* = k \thinspace \exp\negthinspace \left(\frac{Q}{R T}\right) P_i
$$

| Symbol | Description | Unit |
|--------|-------------|------|
| $q_i^*$ | Equilibrium loading of component $i$ | mol kg⁻¹ |
| $k$ | Henry coefficient | mol kg⁻¹ Pa⁻¹ |
| $Q$ | Adsorption energy parameter | J mol⁻¹ |
| $P_i$ | Partial pressure of component $i$ | Pa |

---

## Single Site Langmuir

The classic monolayer model, available on a partial-pressure or concentration basis.

**Partial pressure form** (`Single Site Langmuir`):

$$
q_i^* = \frac{m_1 \thinspace b_1 \exp\negthinspace \left(\frac{q_1}{R T}\right) P_i}{1 + b_1 \exp\negthinspace \left(\frac{q_1}{R T}\right) P_i}
$$

**Concentration form** (`Single Site Langmuir (Concentration)`) — replace $P_i$ with $C_i = P_i/(RT)$:

$$
q_i^* = \frac{m_1 \thinspace b_1 \exp\negthinspace \left(\frac{q_1}{R T}\right) C_i}{1 + b_1 \exp\negthinspace \left(\frac{q_1}{R T}\right) C_i}
$$

| Symbol | Description | Unit (P-form / C-form) |
|--------|-------------|------------------------|
| $m_1$ | Saturation capacity | mol kg⁻¹ |
| $b_1$ | Affinity pre-factor | Pa⁻¹ / m³ mol⁻¹ |
| $q_1$ | Adsorption energy | J mol⁻¹ |

---

## Dual Site Langmuir

Two independent Langmuir sites, for heterogeneous surfaces. Available on a partial-pressure or concentration basis.

**Partial pressure form** (`Dual Site Langmuir`):

$$
q_i^* = \frac{m_1 \thinspace b_1 \exp\negthinspace \left(\frac{q_1}{R T}\right) P_i}{1 + b_1 \exp\negthinspace \left(\frac{q_1}{R T}\right) P_i} + \frac{m_2 \thinspace b_2 \exp\negthinspace \left(\frac{q_2}{R T}\right) P_i}{1 + b_2 \exp\negthinspace \left(\frac{q_2}{R T}\right) P_i}
$$

The concentration form (`Dual Site Langmuir (Concentration)`) substitutes $C_i$ for $P_i$ in both terms.

| Symbol | Description | Unit |
|--------|-------------|------|
| $m_1, m_2$ | Saturation capacity of each site | mol kg⁻¹ |
| $b_1, b_2$ | Affinity pre-factor of each site | Pa⁻¹ (P-form) |
| $q_1, q_2$ | Adsorption energy of each site | J mol⁻¹ |

---

## Extended (Competitive) Langmuir

The **Extended** variants make adsorption *competitive*: the loading of each component is reduced by the presence of every other adsorbing component through a shared denominator. Use these for multi-component mixtures where components compete for the same sites. They are available in single-site and dual-site forms, each on a partial-pressure or concentration basis (`Extended Single Site Langmuir`, `Extended Dual Site Langmuir`, and their `(Concentration)` counterparts).

Extended single-site (partial-pressure), with the sum running over every adsorbing component $j$:

$$
q_i^* = \frac{m_{1,i} \thinspace b_{1,i} \exp\negthinspace \left(\frac{q_{1,i}}{R T}\right) P_i}{1 + \sum_{j} b_{1,j} \exp\negthinspace \left(\frac{q_{1,j}}{R T}\right) P_j}
$$

The dual-site extended form adds a second competitive term with its own denominator $1 + \sum_j b_{2,j}\exp(q_{2,j}/RT)\thinspace P_j$. Parameters are entered per component exactly as for the ordinary Langmuir models; the coupling between components is applied automatically.

---

## Anti-Langmuir

For systems where loading accelerates with concentration (e.g. some cooperative-adsorption cases).

$$
q_i^* = \frac{a \thinspace P_i}{1 - b \thinspace P_i}
$$

| Symbol | Description | Unit |
|--------|-------------|------|
| $a$ | Linear coefficient | mol kg⁻¹ Pa⁻¹ |
| $b$ | Non-linearity coefficient | Pa⁻¹ |

---

## Freundlich

An empirical power-law model.

$$
q_i^* = \left(a \thinspace P_i\right)^{1/v}
$$

| Symbol | Description | Unit |
|--------|-------------|------|
| $a$ | Freundlich coefficient | mol kg⁻¹ Pa⁻¹ |
| $v$ | Heterogeneity exponent | - |

---

## Sips

A Langmuir–Freundlich hybrid with a heterogeneity exponent.

$$
q_i^* = q_{sat}\thinspace \frac{\left(b \thinspace P_i\right)^{1/v}}{1 + \left(b \thinspace P_i\right)^{1/v}}
$$

| Symbol | Description | Unit |
|--------|-------------|------|
| $q_{sat}$ | Saturation capacity | mol kg⁻¹ |
| $b$ | Affinity coefficient | Pa⁻¹ |
| $v$ | Heterogeneity exponent | - |

---

## Sips 2 (temperature-dependent)

A Sips form with explicit temperature dependence built into the capacity, affinity and exponent, referenced to a temperature $T_{ref}$.

$$
b_i = B_1 \exp\negthinspace \left(\frac{B_2}{R T}\right), \qquad
s_i = C_1 + C_2\negthinspace \left(\frac{T}{T_{ref}} - 1\right), \qquad
q_{max} = A_1 \exp\negthinspace \left[A_2\negthinspace \left(\frac{T}{T_{ref}} - 1\right)\right]
$$

$$
q_i^* = q_{max}\thinspace \frac{\left(b_i P_i\right)^{s_i}}{1 + \left(b_i P_i\right)^{s_i}}
$$

| Symbol | Description | Unit |
|--------|-------------|------|
| $A_1$ | Saturation capacity at $T_{ref}$ | mol kg⁻¹ |
| $A_2$ | Capacity temperature coefficient | - |
| $B_1$ | Affinity pre-factor | Pa⁻¹ |
| $B_2$ | Affinity energy | J mol⁻¹ |
| $C_1, C_2$ | Exponent and its temperature coefficient | - |
| $T_{ref}$ | Reference temperature | K |

An **Extended Sips 2** competitive multi-component form is also available (kind `extended_sips_2`); it shares the denominator $1 + \sum_j (b_j P_j)^{s_j}$ across components. Below a mole fraction of $10^{-5}$ each term is continued as described in [Competitive Sips kinds near zero mole fraction](#competitive-sips-kinds-near-zero-mole-fraction).

---

## Arrhenius Sips (Schell 2012)

The temperature-dependent Sips form of Schell et al. (*Adsorption* 18 (2012) 49, Eqs. 5–7), fitted there to CO₂, H₂ and N₂ on activated carbon: an Arrhenius saturation capacity and affinity, and an S-shaped (inverse-tangent) exponent. It is competitive, with the same shared denominator as Extended Sips 2 (Ruthven's empirical mixture Sips). Kind `arrhenius_sips`.

$$
n_{\infty,i} = a \exp\negthinspace \left(\frac{b}{R T}\right), \qquad
k_i = A \exp\negthinspace \left(\frac{B}{R T}\right), \qquad
c_i = c_{ref} + c_{amp}\thinspace \frac{2}{\pi} \arctan\negthinspace\left[\beta \left(T - T_c\right)\right]
$$

$$
q_i^* = n_{\infty,i}\thinspace \frac{\left(k_i P_i\right)^{c_i}}{1 + \sum_j \left(k_j P_j\right)^{c_j}}
$$

| Symbol | Parameter | Description | Unit |
|--------|-----------|-------------|------|
| $a$ | `n_inf_a` | Saturation capacity pre-factor | mol kg⁻¹ |
| $b$ | `n_inf_b` | Saturation capacity energy ($0$ for a constant capacity) | J mol⁻¹ |
| $A$ | `k_A` | Affinity pre-factor | Pa⁻¹ |
| $B$ | `k_B` | Affinity energy | J mol⁻¹ |
| $c_{ref}$ | `c_ref` | Exponent at $T = T_c$ | - |
| $c_{amp}$ | `c_amp` | Exponent amplitude ($0$ for a constant exponent) | - |
| $\beta$ | `c_beta` | Exponent steepness | K⁻¹ |
| $T_c$ | `c_Tref` | Exponent centre temperature | K |

Below a mole fraction of $10^{-5}$ each term is continued as described in the next section.

---

## Competitive Sips kinds near zero mole fraction

In Extended Sips 2 and Arrhenius Sips each component's term $z_i = (b_i P y_i)^{s_i}$ (and in [Extended Langmuir–Freundlich](#extended-langmuirfreundlich-temperature-dependent) the power $(P y_i / P_{ref})^{n_i}$ of $z_i = b_i (P y_i / P_{ref})^{n_i}$) has an infinite slope at $y_i = 0$ when $s_i < 1$. The integrator's Newton iteration cannot linearise that, and a component absent from part of the bed (ahead of its front, or after a purge) sits right on it. Below the mole fraction $y_0 = 10^{-5}$ the power is therefore replaced by the quadratic in $y_i$ that is zero at $y_i = 0$ and has the power's value and slope at $y_0$:

$$
z_i = \left(b_i P y_0\right)^{s_i}\thinspace t\left[(2 - s_i) - (1 - s_i)\thinspace t\right], \qquad t = \frac{y_i}{y_0}, \qquad 0 \le y_i < y_0
$$

- **Exact above $y_0$.** For $y_i \ge 10^{-5}$ (10 ppm) the published isotherm is used unchanged, so a trace component's loading is exact down to 10 ppm, at every pressure and temperature. Below $y_0$ the term has a finite Henry slope, $(2 - s_i)(b_i P y_0)^{s_i} / y_0$ at $y_i = 0$, where the power's is infinite; at $y_i = y_0/10$ it is 40 % of the power for the most singular exponent in the adsorbent library ($s = 0.425$, CO₂ on 13X).
- **Smooth and favourable.** The continuation matches the power's value and slope at $y_0$ and, for $s_i < 1$, is concave and increasing on $[0, y_0]$, so the isotherm has no inflection down to $y_i = 0$, where the loading is exactly 0. A slightly negative mole fraction, which the integrator can produce ahead of a front, gives $z_i = (2 - s_i)(b_i P y_0)^{s_i} \tanh t$, which continues the slope at 0 and stays bounded.
- $y_0$ is 1000 times the integrator's absolute tolerance on mole fractions ($10^{-8}$), so a Newton correction near $y_i = 0$ sees an almost linear term.
- For $s_i \ge 2$, where the quadratic would not be monotone (no published fit is near it), the power is used down to $y_i = 0$ and $z_i = 0$ below it.

The single-component evaluation of these kinds (`skarstrom.evaluate_isotherm`, whose `mole_fraction` defaults to 1, and the [isotherm fitter](isotherm-fitting.md)) uses the same function, so what is plotted or fitted is what the solver uses. Under the [IAST mixing rule](#mixture-loadings-the-mixing-rule-and-iast) the pure isotherms are evaluated in partial pressure and absent components take the infinite-dilution limit, so this continuation does not enter; a pure Sips term with $s_i < 1$ has its own continuation inside the IAST solve (see [Trace levels](#mixture-loadings-the-mixing-rule-and-iast) there).

> Before this was introduced (see the changelog), the two Sips kinds multiplied every term by $y_i^2 / (y_i^2 + 10^{-6})$, in the solver only. That halved a component's loading at $y_i = 10^{-3}$, cut it by 1 % at $10^{-2}$ and by 99 % at $10^{-4}$. Extended Langmuir–Freundlich used $y_i^2 / (y_i^2 + 10^{-10})$ in both paths: $10^{-4}$ at $y_i = 10^{-3}$, 1 % at $10^{-4}$, about 50 % at $10^{-5}$.

---

## Dual Site Sips

Two independent Sips sites.

$$
q_i^* = q_{sat,1}\thinspace \frac{\left(b_1 P_i\right)^{1/v_1}}{1 + \left(b_1 P_i\right)^{1/v_1}} + q_{sat,2}\thinspace \frac{\left(b_2 P_i\right)^{1/v_2}}{1 + \left(b_2 P_i\right)^{1/v_2}}
$$

| Symbol | Description | Unit |
|--------|-------------|------|
| $q_{sat,1}, q_{sat,2}$ | Saturation capacity of each site | mol kg⁻¹ |
| $b_1, b_2$ | Affinity of each site | Pa⁻¹ |
| $v_1, v_2$ | Heterogeneity exponent of each site | - |

---

## Redlich–Peterson

A three-parameter model bridging Henry and Langmuir behaviour.

$$
q_i^* = \frac{a \thinspace P_i}{1 + \left(b \thinspace P_i\right)^{v}}
$$

| Symbol | Description | Unit |
|--------|-------------|------|
| $a$ | Linear coefficient | mol kg⁻¹ Pa⁻¹ |
| $b$ | Affinity coefficient | Pa⁻¹ |
| $v$ | Exponent ($0 < v \le 1$) | - |

---

## Toth

A widely used heterogeneous model that reduces to Langmuir when $v = 1$.

$$
q_i^* = \frac{q_{sat}\thinspace b \thinspace P_i}{\left[1 + \left(b \thinspace P_i\right)^{v}\right]^{1/v}}
$$

| Symbol | Description | Unit |
|--------|-------------|------|
| $q_{sat}$ | Saturation capacity | mol kg⁻¹ |
| $b$ | Affinity coefficient | Pa⁻¹ |
| $v$ | Heterogeneity exponent | - |

---

## Toth 2 (temperature-dependent)

A Toth form in which the capacity, affinity and exponent all vary with temperature about a reference $T_0$.

$$
q_{s} = q_{sat}\exp\negthinspace \left[x\left(1 - \frac{T}{T_0}\right)\right], \qquad
b = b_0 \exp\negthinspace \left[\frac{\Delta H_0}{R T_0}\left(\frac{T_0}{T} - 1\right)\right], \qquad
w = w_0 + \alpha\left(1 - \frac{T_0}{T}\right)
$$

$$
q_i^* = \frac{q_{s}\thinspace b \thinspace P_i}{\left[1 + \left(b \thinspace P_i\right)^{w}\right]^{1/w}}
$$

| Symbol | Description | Unit |
|--------|-------------|------|
| $q_{sat}$ | Saturation capacity at $T_0$ | mol kg⁻¹ |
| $x$ | Capacity temperature coefficient | - |
| $T_0$ | Reference temperature | K |
| $b_0$ | Affinity at $T_0$ | Pa⁻¹ |
| $\Delta H_0$ | Isosteric heat parameter | J mol⁻¹ |
| $w_0$ | Exponent at $T_0$ | - |
| $\alpha$ | Exponent temperature coefficient | - |

---

## Extended Langmuir–Freundlich (temperature-dependent)

The competitive Langmuir–Freundlich form of Yang (*Gas Separation by Adsorption Processes*, 1987), with the capacity linear in $T$, an Arrhenius affinity and an exponent whose inverse is linear in $1/T$. It is the "extended Langmuir–Freundlich" isotherm of toPSAil, whose $k_1$–$k_6$ map one to one (with $P_{ref} = 10^5$ Pa for its bar-based $k_3$). Kind `extended_langmuir_freundlich`.

$$
q_{s,i} = k_1 + k_2\thinspace T, \qquad
b_i = k_3 \exp\negthinspace \left(\frac{k_4}{T}\right), \qquad
\frac{1}{n_i} = k_5 + \frac{k_6}{T}
$$

$$
q_i^* = q_{s,i}\thinspace \frac{b_i \left(P_i / P_{ref}\right)^{n_i}}{1 + \sum_j b_j \left(P_j / P_{ref}\right)^{n_j}}
$$

| Symbol | Description | Unit |
|--------|-------------|------|
| $k_1$ | Saturation capacity at 0 K (intercept) | mol kg⁻¹ |
| $k_2$ | Capacity temperature slope | mol kg⁻¹ K⁻¹ |
| $k_3$ | Affinity pre-factor (per $P_{ref}^{\,n}$) | - |
| $k_4$ | Affinity temperature | K |
| $k_5$ | Inverse exponent intercept | - |
| $k_6$ | Inverse exponent temperature | K |
| $P_{ref}$ | Pressure unit of the fit (default $10^5$ Pa) | Pa |

The affinity multiplies the power, $b P^n$, where Sips 2 raises the affinity with the pressure, $(b P)^n$; the two coincide only when $k_2 = k_6 = 0$ (then Extended Sips 2 with $A_1 = k_1$, $A_2 = 0$, $C_1 = 1/k_5$, $C_2 = 0$, $B_1 = k_3^{k_5} / P_{ref}$, $B_2 = R k_4 k_5$). Below a mole fraction of $10^{-5}$ the power $(P y_i / P_{ref})^{n_i}$ is continued exactly as in the competitive Sips kinds (with $b_i P$ there read as $P / P_{ref}$ and $s_i$ as $n_i$; see [Competitive Sips kinds near zero mole fraction](#competitive-sips-kinds-near-zero-mole-fraction)), so the isotherm is exact for $y_i \ge 10^{-5}$ at every temperature, whatever $n_i(T)$, and has a finite Henry slope below. Under the IAST mixing rule each component's term is a pure Sips isotherm with affinity $b_i^{1/n_i} / P_{ref}$ and exponent $n_i$.

---

## GAB and BET (water and other vapours)

The Guggenheim–Anderson–de Boer (GAB) isotherm describes multilayer adsorption of a condensable vapour, typically water, as a function of the **relative pressure** (relative humidity for water) $x_i = P_i / p_{sat,i}(T)$ with $P_i = y_i P$. Two kinds share the equation and differ only in how $C$ and $K$ depend on temperature.

$$
q_i^* = \frac{q_m\thinspace C\thinspace K\thinspace x_i}{(1 - K x_i)\thinspace(1 - K x_i + C K x_i)}
$$

$C$ measures how much more strongly the first layer binds than the bulk liquid, and $K$ how much the multilayer differs from it ($K = 1$: the multilayer is liquid). **BET** is the special case $K = 1$: use kind `gab` with $K_0 = 1$, $\Delta H_K = 0$. It has no kind of its own, so it shares the same code and the same safeguards below.

**Kind `gab`** (constants with an optional Arrhenius temperature dependence):

$$
q_m = q_{m0} \exp\negthinspace\left(\frac{\Delta H_m}{RT}\right), \qquad
C = C_0 \exp\negthinspace\left(\frac{\Delta H_C}{RT}\right), \qquad
K = K_0 \exp\negthinspace\left(\frac{\Delta H_K}{RT}\right)
$$

| Symbol | Name | Description | Unit |
|--------|------|-------------|------|
| $q_{m0}$ | `q_m0` | Monolayer capacity (its value when $\Delta H_m = 0$) | mol kg⁻¹ |
| $\Delta H_m$ | `dH_m` | Temperature coefficient of $q_m$ | J mol⁻¹ |
| $C_0$ | `C0` | GAB constant $C$ pre-factor | - |
| $\Delta H_C$ | `dH_C` | Temperature coefficient of $C$ | J mol⁻¹ |
| $K_0$ | `K0` | GAB constant $K$ pre-factor | - |
| $\Delta H_K$ | `dH_K` | Temperature coefficient of $K$ | J mol⁻¹ |

With $\Delta H_m = \Delta H_C = \Delta H_K = 0$, the loading depends on temperature only through $p_{sat}(T)$: the isotherms at different temperatures collapse onto one curve in relative humidity. That is the form of Stampi-Bombelli et al. (*Adsorption* 26 (2020) 1183, Eq. 6 and Table 2; their $c_m$, $c_G$, $K_{ads}$ are $q_{m0}$, $C_0$, $K_0$) and of Grossmann et al. (*Ind. Eng. Chem. Res.* 64 (2025), Eq. 12 and Table 5).

**Kind `gab_young`** (layer energies of Young et al., *Energy Environ. Sci.* 14 (2021) 5377, also used by Low et al., *J. Chem. Eng. Data* 68 (2023), Eqs. 17–22 and Table 9, and by pyDAC):

$$
C = \exp\negthinspace\left(\frac{E_1 - E_{10}}{RT}\right), \qquad
K = \exp\negthinspace\left(\frac{E_{2\text{–}9} - E_{10}}{RT}\right), \qquad
E_1 = C_E - e^{D T}, \quad E_{2\text{–}9} = F + G\thinspace T, \quad E_{10} = E_{10,A} + E_{10,B}\thinspace T
$$

| Symbol | Name | Description | Unit |
|--------|------|-------------|------|
| $q_m$ | `q_m` | Monolayer capacity | mol kg⁻¹ |
| $C_E$ | `C_E` | First-layer energy constant (the papers' $C$) | J mol⁻¹ |
| $D$ | `D` | First-layer energy exponent | K⁻¹ |
| $F$ | `F` | Energy of layers 2–9 at 0 K | J mol⁻¹ |
| $G$ | `G` | Temperature slope of the layers 2–9 energy | J mol⁻¹ K⁻¹ |
| $E_{10,A}$ | `E10_A` | Heat of condensation of the bulk liquid at 0 K (water: 57220) | J mol⁻¹ |
| $E_{10,B}$ | `E10_B` | Its temperature slope (water: −44.38) | J mol⁻¹ K⁻¹ |

The papers print $C_E$, $F$ and $G$ in kJ: enter them in J ($C_E = 46.83$ kJ mol⁻¹ is `46830`). The defaults are the Lewatit VP OC 1065 values of Young et al. as reprinted by Low et al. (SI Table S5).

**Saturation pressure.** $p_{sat,i}(T)$ comes from CoolProp for the component (its saturation ancillary equation: a closed-form fit of the full equation of state, about 0.1 µs per call). For water it agrees with the IAPWS-95 saturation pressure to within 0.02 % from 273 K to the critical point. Above the critical temperature $p_{sat}$ is held at the critical pressure. Below 0.75 times the triple-point temperature it is held at its value there; between that and the triple point, the curve over the (supercooled) liquid is extrapolated. Any CoolProp fluid can use a GAB kind, not only water.

**Safe limits.** The GAB has a pole at $K x = 1$, and a simulation can briefly produce $x > 1$ (supersaturation; condensation is not modelled) or $y_i < 0$ ahead of a front. So the loading is exact for $0 \le x \le x_c$ with $x_c = \min(1,\ 0.95/K)$, and continued so that it stays finite, smooth (continuous slope) and increasing everywhere:

- for $x < 0$ it follows the Henry line $q = q_m C K x$, the GAB's slope at $x = 0$;
- for $x > x_c$ it follows the tangent at $x_c$, $q = q(x_c) + q'(x_c)\thinspace(x - x_c)$.

With $K \le 0.95$, the usual case for water on amine sorbents, the GAB is therefore exact up to saturation. With $K$ near 1 (BET) it is exact up to $K x = 0.95$. GAB parameters are usually fitted to data up to $x \approx 0.9$, so loadings above that are extrapolations in any case.

**Mixtures.** A GAB component's loading depends only on its own relative pressure: the kinds are not competitive. Under the default mixing rule, water adsorbs independently of the other components, which keep their own isotherms. Under the IAST mixing rule a GAB kind takes part like any pure isotherm: its spreading pressure has the closed form $\psi = q_m \left[\ln(1 + (C - 1) K x) - \ln(1 - K x)\right]$ (and the matching integral of the tangent above $x_c$). Coupled CO₂–H₂O models, where water changes the CO₂ loading, are separate kinds: see [CO₂ with adsorbed water](#co2-with-adsorbed-water-toth-wadst-and-humid-toth).

**From Python.** `skarstrom.evaluate_isotherm`, `fit_isotherm` and `evaluate_iast` evaluate the GAB kinds for water unless `component=` (`components=` for IAST) names another vapour; see [Isotherm Fitting](isotherm-fitting.md). In a deck, each component uses its own saturation pressure.

**Heat of adsorption.** The temperature dependence implies an isosteric heat that varies with loading:

$$
q_{st} = R T^2 \left(\frac{\partial \ln P_i}{\partial T}\right)_{q} = \Delta H_{vap}(T) + R T^2 \left(\frac{\partial \ln x_i}{\partial T}\right)_{q},
\qquad \Delta H_{vap} = R T^2 \frac{d \ln p_{sat}}{dT}
$$

At low loading, $q \approx q_m C K x$, so for kind `gab` $q_{st} \to \Delta H_{vap} + \Delta H_m + \Delta H_C + \Delta H_K$ (water at 25 °C: $\Delta H_{vap} \approx 44.0$ kJ mol⁻¹). At high loading it tends towards $\Delta H_{vap} + \Delta H_m + \Delta H_K$, as the multilayer approaches liquid. The energy balance, as for every kind, takes the adsorbent's constant `heat_of_adsorption[<component>]` and does not derive it from the isotherm (under the default isosteric convention it is $q_{st}$; see [Model Conventions](model-conventions.md)). The C++ function `adsim::IsostericHeat(isotherm, component, T, P, y)` (`adsim/equations/isotherm_kinds.h`) evaluates the expression above for any kind at a given state, which helps choose that constant for the loadings of the process. For water, a value near the heat of condensation plus a few kJ mol⁻¹ is typical.

## CO2 with adsorbed water (Toth WADST and humid Toth)

On amine sorbents for direct air capture, adsorbed water raises the CO₂ capacity. Two kinds describe this for the CO₂ component. Its equilibrium loading depends on the **adsorbed water loading** $q_w$ on the same adsorbent in the same cell. That is the actual loading of the component named `Water`, which follows its own isotherm (for example `gab` or `gab_young`) and its own rate model. Water's isotherm does not depend on CO₂ in either model.

Both source papers use the water on the sorbent, not the equilibrium loading at the gas humidity. The difference matters when water uptake is slow, as it usually is in direct air capture: there the CO₂ capacity follows the water actually adsorbed, which lags the gas humidity. A deck without a `Water` component, or a cell with no water adsorbed, gives the dry isotherm. The deck checks warn when a deck uses one of these kinds without a `Water` component. A negative $q_w$ (a transient solver iterate) is taken as 0. At the start of a run the loadings are put in equilibrium in two passes, water first, so the CO₂ starts in equilibrium with the initial water loading.

**Kind `toth_wadst`**: the weighted-average dual-site Toth of Young et al., *Energy Environ. Sci.* 14 (2021) 5377, Eq. 14. A dry and a wet Toth site, weighted by the water loading:

$$
q_i^* = \left(1 - e^{-A/q_w}\right) q_{dry}(P_i, T) + e^{-A/q_w}\thinspace q_{wet}(P_i, T)
$$

with each site the temperature-dependent Toth of that paper (Eqs. 1–4):

$$
q_{site} = \frac{q_\infty\thinspace b\thinspace P_i}{\left[1 + (b P_i)^{\tau}\right]^{1/\tau}}, \qquad
q_\infty = q_{\infty,0} \exp\negthinspace\left[\chi\left(1 - \frac{T}{T_0}\right)\right], \qquad
b = b_0 \exp\negthinspace\left(\frac{\Delta H}{R T}\right), \qquad
\tau = \tau_0 + \alpha\left(1 - \frac{T_0}{T}\right)
$$

The affinity has no reference temperature, as in the paper, so $b_0$ is the affinity at infinite temperature. [Toth 2](#toth-2-temperature-dependent) writes the same site as $b_0' \exp[\Delta H/(R T_0)(T_0/T - 1)]$ with $b_0' = b_0 \exp[\Delta H/(R T_0)]$. At $q_w = 0$ the loading is exactly the dry site, and the weight stays negligible up to $q_w \approx A/30$. The parameter names are `T_0`, then `q_inf_dry`, `chi_dry`, `b0_dry`, `dH_dry`, `tau0_dry`, `alpha_dry`, the same six with `_wet`, and `A`.

| Symbol | Name | Description | Unit | Default (Lewatit VP OC 1065) |
|--------|------|-------------|------|------|
| $T_0$ | `T_0` | Reference temperature of $\chi$ and $\tau$ | K | 298.15 |
| $q_{\infty,0}$ | `q_inf_dry`, `q_inf_wet` | Site capacity at $T_0$ | mol kg⁻¹ | 4.86, 9.035 |
| $\chi$ | `chi_dry`, `chi_wet` | Capacity temperature coefficient | - | 0, 0 |
| $b_0$ | `b0_dry`, `b0_wet` | Affinity pre-exponential | Pa⁻¹ | 2.85×10⁻²¹, 1.230×10⁻¹⁸ |
| $\Delta H$ | `dH_dry`, `dH_wet` | Affinity energy ($-\Delta H_{ads}$, positive) | J mol⁻¹ | 117 798, 203 687 |
| $\tau_0$ | `tau0_dry`, `tau0_wet` | Exponent at $T_0$ | - | 0.209, 0.053 |
| $\alpha$ | `alpha_dry`, `alpha_wet` | Exponent temperature coefficient | - | 0.523, 0.053 |
| $A$ | `A` | Water loading scale of the wet-site weight | mol kg⁻¹ | 1.532 |

**Kind `toth_humid`**: the humid Toth of Stampi-Bombelli, Storione and Mazzotti, *Adsorption* 26 (2020) 1183–1197 (corrected 2022), Eqs. 11–14. This is [Toth 2](#toth-2-temperature-dependent), with the same parameter names, plus two water terms:

$$
q_s = \frac{q_{sat}\exp\negthinspace\left[x\left(1 - \frac{T}{T_0}\right)\right]}{1 - \gamma\thinspace q_w}, \qquad
b = b_0 \exp\negthinspace\left[\frac{\Delta H_0}{R T_0}\left(\frac{T_0}{T} - 1\right)\right]\left(1 + \beta\thinspace q_w\right), \qquad
w = w_0 + \alpha\left(1 - \frac{T_0}{T}\right)
$$

$$
q_i^* = \frac{q_s\thinspace b\thinspace P_i}{\left[1 + (b P_i)^{w}\right]^{1/w}}
$$

At $q_w = 0$ this is Toth 2 exactly. $\Delta H_0 > 0$ makes the affinity fall with temperature. The defaults are the paper's corrected APDES-NFC values: $q_{sat} = 2.38$ mol kg⁻¹, $x = 0$, $T_0 = 296$ K, $b_0 = 70.7$ kPa⁻¹ (0.0707 Pa⁻¹), $\Delta H_0 = 57\,047$ J mol⁻¹, $w_0 = 0.4148$, $\alpha = -1.606$, $\gamma = 0.0016$ kg mol⁻¹ and $\beta = 59.1$ kg mol⁻¹. Young et al. (2021) refit this form to Lewatit VP OC 1065 with $\gamma = -0.137$ and $\beta = 5.612$ kg mol⁻¹ on their dry Toth.

| Symbol | Name | Description | Unit |
|--------|------|-------------|------|
| $\gamma$ | `gamma` | Capacity enhancement by adsorbed water | kg mol⁻¹ |
| $\beta$ | `beta` | Affinity enhancement by adsorbed water | kg mol⁻¹ |

**Mixtures.** Water changes the CO₂ sites themselves; it does not compete through a shared denominator. Neither kind has a pure-component form, so the [IAST mixing rule](#mixture-loadings-the-mixing-rule-and-iast) rejects them with a deck-check error. Use the mixing rule `none`, the default. Give the other components (N₂, …) their own kinds as usual. These kinds read the water loading but do not change it.

**From Python and the fitter.** `skarstrom.evaluate_isotherm(..., water_loading=q_w)` evaluates either kind at a given water loading (default 0, the dry isotherm). `skarstrom.fit_isotherm(..., water_loading=...)` takes the water loading of each point and fits the way the papers do. Points at $q_w = 0$ fit the dry Toth block. Points at fixed $q_w > 0$ fit the water terms: $A$ and the wet site for WADST, $\gamma$ and $\beta$ for the humid Toth. A block that the data cannot identify is held. So with humid data only, pass the dry fit's values in `fixed=`; with both kinds of data, everything is fitted together. In C++ these are the `q_w` member of `isofit::Point` and the `q_water` argument of `adsim::EvalIsotherm`.

**Worked example: APDES-NFC at a humid air feed.**

!!! note "New in 0.9"
    The example and the [`Illustrative_Humid_DAC_Breakthrough`](examples.md) deck.

The library entry [APDES-NFC](adsorbent-library.md#apdes-nfc-stampi2020) carries a `gab` isotherm for water and a `toth_humid` isotherm for CO₂. At the feed of its source paper (400 ppm CO₂, $y_{H_2O} = 0.0115$, 20 °C, 1 bar, about 49 % relative humidity) the water loading is 2.82 mol/kg, and that water raises the CO₂ loading from 0.82 to 2.07 mol/kg:

```python
import skarstrom as sk

sim = sk.Simulation.new(title="APDES-NFC at the feed state")
for name in ("Nitrogen", "CarbonDioxide", "Water"):
    sim.add_component(name)
adsorbent_id = sim.add_library_adsorbent("apdes_nfc_stampi2020")
ads = sim.adsorbents.by_id(adsorbent_id)

def params(node):
    return {name: node.get(name) for name in node.fields}

T, P = 293.15, 1e5
q_w = sk.evaluate_isotherm("gab", params(ads.isotherm["Water"]), pressure=0.0115 * P, temperature=T)[0]
co2 = params(ads.isotherm["CarbonDioxide"])
dry = sk.evaluate_isotherm("toth_humid", co2, pressure=400e-6 * P, temperature=T)[0]
humid = sk.evaluate_isotherm("toth_humid", co2, pressure=400e-6 * P, temperature=T, water_loading=q_w)[0]
print(f"water {q_w:.2f}, CO2 dry {dry:.2f}, CO2 humid {humid:.2f} mol/kg")
# water 2.82, CO2 dry 0.82, CO2 humid 2.07 mol/kg
```

In a bed, the CO₂ isotherm reads the water actually adsorbed in each cell, so the humid capacity appears only behind the water front. The illustrative deck `Illustrative_Humid_DAC_Breakthrough` runs this feed through a column of APDES-NFC: the water front leaves the bed first, then the CO₂ front, with the humid loading behind it. See [Humid feeds](tsa-regeneration.md#6-humid-feeds) in the TSA tutorial.

---

## Mixture loadings: the mixing rule and IAST

Every adsorbent has a **mixing rule**, which says how its components' isotherms combine into mixture loadings. It is one switch per adsorbent: the registry parameter `adsorbent[<id>].mixing_rule`, `ads.mixing_rule` or `ads.set_mixing_rule(kind)` in Python.

- **`none`** (the default, and what every deck written before the switch existed loads as). Each component's isotherm is evaluated on its own at its partial pressure. Components compete only if they share an extended (competitive) kind, through that kind's own denominator.
- **`iast`**: **Ideal Adsorbed Solution Theory** (Myers and Prausnitz, *AIChE J.* 11 (1965) 121). The mixture loadings are predicted from the components' **pure-component** isotherms alone. Each component keeps its own kind and parameters, so the kinds can differ (Toth for one component and dual-site Langmuir for another, for example). No mixture data are needed.

Under IAST the adsorbed phase is treated as an ideal solution in which every component sits at the same reduced spreading pressure $\psi$. With $q_i^\circ$ the pure isotherm of component $i$ at the local temperature:

$$
\psi_i(P) = \int_0^{P} \frac{q_i^\circ(p)}{p}\thinspace dp, \qquad
P \thinspace y_i = x_i \thinspace P_i^\circ, \qquad
\psi_i(P_i^\circ) = \psi \ \text{for every } i, \qquad
\sum_i x_i = 1
$$

$$
\frac{1}{q_T} = \sum_i \frac{x_i}{q_i^\circ(P_i^\circ)}, \qquad q_i^* = x_i \thinspace q_T
$$

| Symbol | Description | Unit |
|--------|-------------|------|
| $\psi$ | Reduced spreading pressure of the mixture | mol kg⁻¹ |
| $P_i^\circ$ | Pressure at which pure $i$ reaches $\psi$ | Pa |
| $x_i$ | Mole fraction of $i$ in the adsorbed phase | - |
| $q_T$ | Total mixture loading | mol kg⁻¹ |

**How it is solved.** In every cell and at every residual evaluation, the engine finds $\psi$ by a bracketed Newton iteration on $\ln \sum_i x_i$, warm-started from the cell's previous value. The spreading-pressure integral has a closed form for Henry, the Langmuir and dual-site Langmuir kinds (both bases), Freundlich, Sips, dual-site Sips, Sips 2, Arrhenius Sips, anti-Langmuir and the GAB kinds. Toth, Toth 2 and Redlich–Peterson have no elementary integral, so they use a 48-point Gauss–Legendre rule in $\ln p$. The inverse $\psi_i \mapsto P_i^\circ$ is analytic for Henry, Langmuir, Freundlich, Sips and anti-Langmuir, and a Newton iteration otherwise. A kind without an IAST form is an error, never treated as inert; the CO₂-with-water kinds (`toth_wadst`, `toth_humid`) are rejected by the deck checks.

**Limits that hold exactly.** A single adsorbing component gives its pure isotherm at its partial pressure. Langmuir isotherms with **equal saturation capacities** give the extended Langmuir model. With unequal capacities, extended Langmuir is thermodynamically inconsistent and IAST is not.

**Absent components.** A component with zero or negative gas mole fraction, which the integrator produces ahead of a front, gets the infinite-dilution limit of IAST, $q_i^* = y_i P q_T / P_i^\circ(\psi)$. That is linear in $y_i$, so the loading is smooth through $y_i = 0$, as the extended kinds are. Inert components never adsorb.

**Trace levels.** Two cases near zero mole fraction are handled inside the IAST solve. A positive mole fraction at or below $10^{-150}$ (the decaying toe of a front reaches subnormal values, where the spreading pressure underflows) takes the infinite-dilution limit above. A Sips-family term (Sips, Sips 2, Arrhenius Sips, Extended Langmuir–Freundlich, and the extended Sips kinds in their pure limit) with exponent $n < 1$ has an infinite Henry slope, so the infinite-dilution loading $y_i P q_T / P_i^\circ(\psi)$ of an absent component would diverge as $\psi \to 0$, at the edge of a bed that starts in an inert gas. Below $p_0 = 10^{-4} P$ ($P$ the total pressure) such a term is replaced by $q = A t / (1 + B t)^k$, $t = p / p_0$, $0 < k \le 1$, whose three constants match the Sips value and slope at $p_0$ and the Sips spreading pressure $\psi(p_0) = (q_s / n) \ln(1 + z_0)$. The loading is monotone with a finite Henry slope, and $\psi$ above $p_0$ is the Sips closed form, so the IAST loadings are exact for every component with $y_i \ge 10^{-4}$. Only the IAST solve (and so `skarstrom.evaluate_iast`) uses this; `skarstrom.evaluate_isotherm` and the fitter show the published isotherm.

**Which kinds to use.** Give each component its pure-component kind. An extended kind under IAST is taken in its pure-component limit, so its own competition term is not applied. The solver log warns about such a combination at the start of a run.

**Limitations.** IAST assumes an ideal adsorbed phase. Strongly non-ideal mixtures need activity coefficients (Real Adsorbed Solution Theory, RAST), which the engine does not offer. Each component's isotherm also has to describe its pure loading well down to low pressure, because the spreading-pressure integral starts at zero. IAST costs more per residual evaluation than an extended kind: a few Newton iterations per cell, each inverting every component's spreading pressure. On the example decks with dual-site Langmuir isotherms, a residual evaluation takes about three times as long, while the number of integrator steps barely changes.

To check mixture loadings outside a simulation, call `skarstrom.evaluate_iast(kinds, params, pressure, temperature, mole_fraction)` with one pure kind and parameter set per component.

---

Once you have chosen an isotherm, pair it with a [kinetic rate model](kinetic-models.md) to control how quickly the loading approaches $q_i^*$.
