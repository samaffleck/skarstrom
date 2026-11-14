## Reactions

A **reaction** lives on an adsorbent (a solid), beside its isotherm and rate model. A catalyst is a solid with reactions, usually with an inert isotherm. Mix a catalyst and a sorbent in one layer by mass fraction to get a sorption-enhanced reactor. No other machinery is needed.

Each rate law gives the intrinsic rate $r$ per kilogram of the solid that carries the reaction, in mol of reaction extent per kg per second. Component $i$ is produced at $\nu_i\,\eta\, r$ per kg of that solid ($\nu_i < 0$ for a reactant, $> 0$ for a product; $\eta$ is the [effectiveness factor](#effectiveness-factor)). The rate laws other than `expression` are written in activities $a_i$: partial pressures $a_i = y_i P$ in **bar** (`partial_pressure`, the default) or concentrations $a_i = y_i P/(RT)$ in mol/m³ (`concentration`).

| Rate law | Rate |
|---|---|
| `power_law` | $r = k(T) \prod_i a_i^{n_i}$ |
| `reversible_power_law` | $r = k(T) \prod_i a_i^{n_i} (1 - Q / K(T))$, $Q = \prod_i a_i^{\nu_i}$ |
| `lhhw` | $r = k(T) \prod_i a_i^{n_i} (1 - Q/K(T)) / \mathrm{DEN}^m$, $\mathrm{DEN} = 1 + \sum_j K_j(T) \prod_i a_i^{e_{ji}}$ |
| `expression` | a user expression of `T`, `P`, `c` and `p_<C>`, `c_<C>`, `y_<C>` |

The temperature-dependent constants all have the Arrhenius form:

$$
k(T) = A\, e^{-E/RT}, \qquad K(T) = A_K\, e^{-E_K/RT}, \qquad K_j(T) = A_j\, e^{-E_j/RT}
$$

| Symbol | Parameter (`adsorbent[ID].reaction[ID].rate.*`) | Description | Unit |
|---|---|---|---|
| $A$, $E$ | `pre_exp_factor`, `activation_energy` | Rate constant | mol kg⁻¹ s⁻¹ (activity)$^{-\sum_i n_i}$; J mol⁻¹ |
| $n_i$ | `order[C]` | Reaction order; defaults to $-\nu_i$ for reactants and 0 for products | – |
| $A_K$, $E_K$ | `equilibrium_pre_exp_factor`, `equilibrium_energy` | Equilibrium constant (entered; there is no thermochemistry table) | (activity)$^{\sum_i \nu_i}$; J mol⁻¹ |
| $A_j$, $E_j$ | `adsorption[j].pre_exp_factor`, `adsorption[j].energy` | Adsorption constant of denominator term $j$ | (activity)$^{-\sum_i e_{ji}}$; J mol⁻¹ |
| $e_{ji}$ | `adsorption[j].exponent[C]` | Exponent of component $i$ in term $j$ | – |
| $m$ | `denominator_exponent` | Exponent of the LHHW denominator (default 1) | – |

With $K_j = A_j e^{-E_j/RT}$, an adsorption constant that falls with temperature has $E_j < 0$ (an adsorption enthalpy). An `lhhw` rate has the switch `rate.reversible` (default on). With it off, the factor $(1 - Q/K)$ is dropped and the rate is $k \prod_i a_i^{n_i} / \mathrm{DEN}^m$.

The reversible factor is evaluated as the difference of two monomials, $k\,[\prod_i a_i^{n_i} - K^{-1}\prod_i a_i^{n_i + \nu_i}]$, never through $Q$, so a zero product activity gives a finite rate. Each activity is clipped at 0, and an activity raised to a negative power is floored (see [Start-up of fast reactions](#start-up-of-fast-reactions)); the two forms agree wherever no activity is clipped or floored.

**The `expression` rate** gives $r$ in mol kg⁻¹ s⁻¹ directly and ignores the rate basis. Its variables are `T` [K], `P` [Pa], `c` = $P/(RT)$ [mol/m³] and, for each component, `p_<C>` = $y_C P$ in **Pa** (not bar), `c_<C>` [mol/m³] and `y_<C>` [–]. A partial-pressure law fitted in bar needs `p_<C>/1e5`.

**Temperature and composition.** Every rate is evaluated with the gas composition and pressure of the cell, at the catalyst temperature: the solid temperature under [local thermal non-equilibrium](energy-balance.md#local-thermal-non-equilibrium-ltne), the gas temperature otherwise.

### Effectiveness factor

The effectiveness factor $\eta$ multiplies the intrinsic rate. It is `constant` (the `effectiveness` value, default 1) or `thiele`. The `thiele` model computes $\eta$ in each cell from the generalised (Aris–Bischoff) Thiele modulus of a key reactant A, for an isothermal sphere without film resistance:

$$
\phi = \frac{R}{3}\,\frac{\rho_p\,|\nu_A|\,|r_s|}{\sqrt{2\,D_\text{eff}\displaystyle\int_{C_{A,\text{eq}}}^{C_{A,s}} \rho_p\,|\nu_A|\, r\, dC_A}}, \qquad
\eta = \frac{1}{\phi}\left[\frac{1}{\tanh 3\phi} - \frac{1}{3\phi}\right], \qquad
D_\text{eff} = \frac{\varepsilon_p}{\tau}\, D_{m,A}
$$

| Symbol | Description | Unit |
|---|---|---|
| $R$ | Particle radius, half the layer's `particle_diameter` | m |
| $\rho_p$ | Particle density of the solid (`solid_density`) | kg m⁻³ |
| $r_s$ | Intrinsic rate at the cell's gas state (the particle surface) | mol kg⁻¹ s⁻¹ |
| $C_{A,s}$ | Concentration of A in the cell's gas | mol m⁻³ |
| $\varepsilon_p$ | The adsorbent's `particle_void_fraction` | – |
| $\tau$ | The reaction's `tortuosity` (default 3) | – |
| $D_{m,A}$ | Molecular diffusivity of A in the cell's gas ([Equation of State](equation-of-state.md)) | m² s⁻¹ |

Inside the particle the other components follow the stoichiometric line $C_j = C_{j,s} + (\nu_j/\nu_A)(C_A - C_{A,s})$ (clipped at 0): every component diffuses with the same effective diffusivity, and temperature and pressure are those of the surface. $C_{A,\text{eq}}$ is the concentration of A on that line where the rate is zero (the equilibrium, found by regula falsi) or, for a rate that does not change sign, where a concentration on the line reaches zero. The integral uses 8-point Gauss–Legendre quadrature. The key reactant is `key_component` when it is a reactant of the reaction, else the first reactant. The particle is a sphere of radius $d_p/2$ also in a layer of cylinders or a structured layer. Where the modulus cannot be formed (no reactant, $r_s = 0$, a non-positive $R$, $\rho_p$ or $D_\text{eff}$), $\eta = 1$.

### Source terms and heat of reaction

The reactions of solid $k$ act at that solid's bulk density $\rho_{b,k}$ (its mass per bed volume, see [Mixed Adsorbent Layers](mixed-layers.md)). Component $i$ gains, per unit bed volume,

$$
S_i = \sum_k \rho_{b,k} \sum_r \nu_{r,i}\, \eta_{k,r}\, r_{k,r} \qquad [\text{mol m}^{-3}\,\text{s}^{-1}]
$$

$S_i$ enters the [component balance](mass-balance.md) of $i$, and the net $\sum_i S_i$ enters the overall mole balance, so a mole-changing reaction drives the velocity field as adsorption does. Each reaction also releases

$$
Q_\text{rxn} = \sum_k \rho_{b,k} \sum_r \left(-\Delta H_r(T)\right) \eta_{k,r}\, r_{k,r} \qquad [\text{W m}^{-3}]
$$

into the [energy balance](energy-balance.md) (the solid's under LTNE), with $T$ the catalyst temperature; it is 0 in an isothermal run. $\Delta H_r$ is per mol of extent, negative for an exothermic reaction. Two modes (`heat_mode`):

- `kirchhoff` (the default): $\Delta H_r(T) = \Delta H_\text{ref} + \sum_i \nu_i \int_{T_\text{ref}}^{T} C_{p,i}\, dT$, with the gas molar heat capacities the energy balance uses;
- `constant`: $\Delta H_r = \Delta H_\text{ref}$ at every temperature.

$\Delta H_\text{ref}$ is `heat_of_reaction` [J mol⁻¹] and $T_\text{ref}$ is `heat_reference_temperature` (default 298.15 K). Under the default mode, enter the heat of reaction at $T_\text{ref}$, not at the operating temperature.

### Building reactions in Python

!!! note "New in 0.9"
    `Adsorbent.add_reaction`, the `Reaction` setters and `sim.describe()["adsorbents"][i]["reactions"]`. `sk.to_script` writes them.

Create the structure with `ads.add_reaction(name, nu, kind=, basis=, effectiveness=, key_component=, adsorption_terms=)`. Then set every value through the returned `Reaction` object. The values are the registry paths `adsorbent[ID].reaction[ID].*`.

```python
import skarstrom as sk

sim = sk.Simulation.new(title="reformer")
for c in ("Methane", "Water", "CarbonMonoxide", "Hydrogen", "CarbonDioxide"):
    sim.add_component(c)
ni = sim.add_adsorbent("Ni catalyst", solid_density=2000, heat_capacity=900,
                       thermal_conductivity=0.5, particle_void_fraction=0.4)
sim.add_adsorbent("Sorbent", solid_density=1500, heat_capacity=900,
                  thermal_conductivity=0.5, particle_void_fraction=0.4)
sim.bed.set_shape("cylinder", diameter=0.05)
sim.wall.set_material("ss316l")
sim.wall.set_boundary("zero_flux")
sim.add_layer("Layer 1", composition={"Ni catalyst": 0.3, "Sorbent": 0.7},
              momentum="kozeny_carman", length=0.2, cells=10,
              bed_void_fraction=0.4, particle_diameter=3e-3)
sim.initial_conditions.pressure = 1e5
sim.initial_conditions.temperature = 800
sim.initial_conditions.y["Water"] = 1.0
feed = sim.add_step("Feed", time=2.0)
feed.bottom.set_streams("inlet", [sk.Stream("Feed inlet", flow=sk.velocity(0.05),
                                            y={"Methane": 0.2, "Water": 0.8}, temperature=800)])
feed.top.set_pressure("outlet", 1e5)

smr = ni.add_reaction("SMR", {"Methane": -1, "Water": -1, "CarbonMonoxide": 1, "Hydrogen": 3},
                      kind="lhhw", adsorption_terms=2, heat_of_reaction=206100.0)
smr.rate.pre_exp_factor = 1e3
smr.rate.activation_energy = 8e4
smr.rate.order["Hydrogen"] = 0.5
smr.rate.adsorption[2].exponent["Hydrogen"] = -1.0
wgs = ni.add_reaction("WGS", {"CarbonMonoxide": -1, "Water": -1, "CarbonDioxide": 1, "Hydrogen": 1},
                      kind="reversible_power_law", basis="concentration",
                      effectiveness="thiele", key_component="CarbonMonoxide")
sim.kpi.track_reaction("Methane", products=["Hydrogen"])
print(sk.to_script(sim)[:200])
```

The SMR heat of reaction, 206 100 J mol⁻¹, is the standard value at 298.15 K: with the default `kirchhoff` mode the engine corrects it to the bed temperature (800 K here) through the gas heat capacities, so it is not entered at 800 K.

The structure setters of a `Reaction` replace what was there. A new rate law keeps the parameters that it shares with the old one.

| Call | Effect |
|---|---|
| `ads.add_reaction(name, nu, kind="power_law", basis=None, effectiveness=None, key_component=None, adsorption_terms=0, **values)` | A new reaction. `kind` is one of `skarstrom.build_kinds()["reaction_rate"]`. `values` set fields such as `heat_of_reaction`. |
| `rx.set_kind(kind)` | A new rate law. |
| `rx.set_stoichiometry({component: nu})` | Replace the stoichiometry. |
| `rx.set_basis("partial_pressure" \| "concentration")` | The activity of the rate. |
| `rx.set_effectiveness("constant" \| "thiele", key_component=None)` | The effectiveness model. A Thiele model needs a reactant of the reaction as `key_component`. |
| `rx.add_adsorption_term()` / `rx.remove_adsorption_term(j)` | The denominator terms of an `lhhw` rate. They are numbered from 1. |
| `ads.remove_reaction(reaction)` | Remove a reaction. IDs are not reused. |
| `ads.reactions` | The reactions: `ads.reactions["SMR"]`, `ads.reactions.by_id(1)`. |

`skarstrom.to_script(sim)` writes the reactions, their rate-law values and their adsorption terms. Running the script gives the same deck. The raw operations (`sim._build("add_reaction", ...)`) are the C++ functions of `adsim/objects/deck_builder.h`. The reaction KPIs are in [Performance Metrics](performance-metrics.md).

### Start-up of fast reactions

!!! note "New in 0.9"
    The Jacobian treatment of reaction mole fractions below and the `REACTION_STIFF` run-health warning are not
    in 0.8. In 0.8 a fast reaction started from a bed without its products can fail with "too much work".

A catalyst often converts many times more than the flow brings: a steam-reforming catalyst at 450 °C can
react 10<sup>4</sup> to 10<sup>6</sup> times the methane fed, so the gas reaches equilibrium within
microseconds of meeting the catalyst. This is stiff, and the integrator (BDF) handles stiffness with long
steps only when its iteration matrix holds the reaction's coupling correctly.

Start-up is where that is hardest. A bed that holds no product yet (steam and argon before a methane feed)
has components whose mole fractions sit far below their absolute tolerance (10<sup>-8</sup>): the products,
and the reactants ahead of their front. The rate laws clip each activity at 0 and floor an activity raised
to a negative power (`kMinActivityBar`, 10<sup>-10</sup> bar; `kMinActivityConcentration`, 10<sup>-6</sup>
mol/m³), and several kinetics raise a product to a negative power (H<sub>2</sub> at −2.5 in the Xu-Froment
steam-reforming rate). The difference quotients of the Jacobian perturb each mole fraction by its error
weight, about 10<sup>-8</sup>, in the direction it is moving. For a mole fraction far below that:

- a falling mole fraction is perturbed below zero, where the rate is clipped, so even a linear term's slope
  is underestimated by the ratio of the mole fraction to the perturbation;
- a power below 1 has an unbounded slope at zero, and the secant across 10<sup>-8</sup> is far flatter than
  the tangent (25 times for y<sup>−2.5</sup> at y = 10<sup>-9</sup>).

The Newton iteration then fails to converge on the fast reaction and the integrator keeps its step at the
reaction's time scale: a run that should take a second takes minutes, or fails with "too much work".

The engine therefore treats the mole fraction of every component that a reaction in the cell involves (its
stoichiometry, an order or an exponent of its rate law) specially in the Jacobian: it is perturbed away from
zero whenever the usual perturbation would cross it, and, when the rate raises it to a power below 1, by
10<sup>-4</sup> of its own value (or of the activity floor's mole fraction for a negative power) when that is
smaller. Decks without reactions, and every other unknown, keep the usual perturbation bit for bit. The
`expression` kind cannot report its powers; its components are perturbed away from zero only.

On a Xu-Froment steam-reforming bed (three reversible LHHW rates, 3.94 % CH<sub>4</sub>, steam to carbon 4,
0.1 % H<sub>2</sub> in argon at 450 °C and 1.5 bar, a 4.2 cm bed fed for 1200 s,
`core/test/integration/test_reaction_startup.cpp`), started from steam and argon only:

| Start | Wall time | Steps | Newton convergence failures |
|---|---|---|---|
| Steam and argon, usual perturbation (0.8) | 52 s | 11 600 | 5 522 |
| Steam and argon, reaction treatment | 0.7 s | 509 | 20 |
| Gas at the equilibrium of the feed | 0.23 s | 123 | 0 |

The outlet methane conversion before breakthrough is the same from either start, and equals the gas-phase
equilibrium of the feed (0.529).

#### When the step still collapses

When a reaction still holds the step far below the flow's time scale, the run-health warning
[`REACTION_STIFF`](run-health.md#reaction_stiff) names the reaction and the component whose near-zero mole
fraction drives it. Remedies, from the least intrusive:

- **Start the bed nearer the reaction's state.** Give the initial gas (`initial.y[<Component>]`) a little of
  the product the warning names, or set it to the gas-phase equilibrium of the feed at the bed's temperature
  and pressure: the first breakthrough is then the reacted gas at once. Decks of fast reforming catalysts that
  fail from a product-free bed run in seconds from the equilibrium gas.
- **Put a trace of the product in the feed** (0.1 % H<sub>2</sub> in a reformer feed is common practice for
  this reason).
- **Check the rate law.** A negative order on a component that is absent at the inlet makes the rate grow
  without bound as the component goes to zero, up to the activity floor; make sure the orders, the
  denominator exponents and the units of the constants (bar, not Pa, on the partial-pressure basis) are those
  of the source.
