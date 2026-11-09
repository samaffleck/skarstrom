## Reactions

A **reaction** lives on an adsorbent (a solid), beside its isotherm and rate model. A catalyst is a solid with reactions, usually with an inert isotherm. Mix a catalyst and a sorbent in one layer by mass fraction to get a sorption-enhanced reactor. No other machinery is needed.

Each rate law gives the rate $r$ per kilogram of the solid that carries the reaction, in mol of reaction extent per kg per second. Component $i$ is produced at $\nu_i r$ ($\nu_i < 0$ for a reactant, $> 0$ for a product). The rate is written in activities: partial pressures in bar (`partial_pressure`, the default) or concentrations in mol/m³ (`concentration`).

| Rate law | Rate |
|---|---|
| `power_law` | $r = k(T) \prod_i a_i^{n_i}$ |
| `reversible_power_law` | $r = k(T) \prod_i a_i^{n_i} (1 - Q / K(T))$, $Q = \prod_i a_i^{\nu_i}$ |
| `lhhw` | $r = k(T) \prod_i a_i^{n_i} (1 - Q/K(T)) / \mathrm{DEN}^m$, $\mathrm{DEN} = 1 + \sum_j K_j(T) \prod_i a_i^{e_{ji}}$ |
| `expression` | a user expression of `T`, `P`, `c` and `p_<C>`, `c_<C>`, `y_<C>` |

$k(T) = A \exp(-E / RT)$. The orders $n_i$ default to $-\nu_i$ for reactants and 0 for products.

The effectiveness factor $\eta$ is `constant` (the `effectiveness` value) or `thiele` (from the Thiele modulus of a key reactant of the reaction).

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
