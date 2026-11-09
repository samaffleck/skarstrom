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
