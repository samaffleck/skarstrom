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

An **Extended Sips 2** competitive multi-component form is also available; it shares the denominator $1 + \sum_j (b_j P_j)^{s_j}$ across components (with a small numerical regularisation applied near zero mole fraction for solver stability).

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

**How it is solved.** In every cell and at every residual evaluation, the engine finds $\psi$ by a bracketed Newton iteration on $\ln \sum_i x_i$, warm-started from the cell's previous value. The spreading-pressure integral has a closed form for Henry, the Langmuir and dual-site Langmuir kinds (both bases), Freundlich, Sips, dual-site Sips, Sips 2, Arrhenius Sips and anti-Langmuir. Toth, Toth 2 and Redlich–Peterson have no elementary integral, so they use a 48-point Gauss–Legendre rule in $\ln p$. The inverse $\psi_i \mapsto P_i^\circ$ is analytic for Henry, Langmuir, Freundlich, Sips and anti-Langmuir, and a Newton iteration otherwise.

**Limits that hold exactly.** A single adsorbing component gives its pure isotherm at its partial pressure. Langmuir isotherms with **equal saturation capacities** give the extended Langmuir model. With unequal capacities, extended Langmuir is thermodynamically inconsistent and IAST is not.

**Absent components.** A component with zero or negative gas mole fraction, which the integrator produces ahead of a front, gets the infinite-dilution limit of IAST, $q_i^* = y_i P q_T / P_i^\circ(\psi)$. That is linear in $y_i$, so the loading is smooth through $y_i = 0$, as the extended kinds are. Inert components never adsorb.

**Which kinds to use.** Give each component its pure-component kind. An extended kind under IAST is taken in its pure-component limit, so its own competition term is not applied. The solver log warns about such a combination at the start of a run.

**Limitations.** IAST assumes an ideal adsorbed phase. Strongly non-ideal mixtures need activity coefficients (Real Adsorbed Solution Theory, RAST), which the engine does not offer. Each component's isotherm also has to describe its pure loading well down to low pressure, because the spreading-pressure integral starts at zero. IAST costs more per residual evaluation than an extended kind: a few Newton iterations per cell, each inverting every component's spreading pressure. On the example decks with dual-site Langmuir isotherms, a residual evaluation takes about three times as long, while the number of integrator steps barely changes.

To check mixture loadings outside a simulation, call `skarstrom.evaluate_iast(kinds, params, pressure, temperature, mole_fraction)` with one pure kind and parameter set per component.

---

Once you have chosen an isotherm, pair it with a [kinetic rate model](kinetic-models.md) to control how quickly the loading approaches $q_i^*$.
