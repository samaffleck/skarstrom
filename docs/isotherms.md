## Isotherm Models

The **isotherm** describes the equilibrium adsorbed-phase loading of a component, $q_i^*$, as a function of the local gas state (temperature, pressure and composition). It is the thermodynamic ceiling that the kinetic [rate model](kinetic-models.md) drives the actual loading $q_i$ towards.

Each component in each adsorbent material is assigned its own isotherm. In a deck, every component of every adsorbent has an **isotherm type** and the parameters of that model. To plot or check an isotherm outside a simulation, evaluate it with the engine's own code from Python: `skarstrom.evaluate_isotherm(kind, params, pressure, temperature)` (see [Isotherm Fitting](isotherm-fitting.md)).

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

Once you have chosen an isotherm, pair it with a [kinetic rate model](kinetic-models.md) to control how quickly the loading approaches $q_i^*$.
