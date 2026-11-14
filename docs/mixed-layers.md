# Mixed Adsorbent Layers

<!-- doc-setup
import skarstrom
sim = skarstrom.examples.Illustrative_Diluted_Bed.load()
-->

!!! note "New in 0.9"
    This page, the per-adsorbent series (`sim.series(..., adsorbent=)`, `sim.series_adsorbents()`, `skarstrom.series_manifest()`) and the bundled deck `Illustrative_Diluted_Bed` are new in 0.9. Mixed layers themselves are older.

A **layer** of the bed is packed with one *or several* solids, well mixed, in fixed **mass fractions**. This is how Skarstrom models a blend of two adsorbents, an adsorbent diluted with inert balls (a common way to damp the temperature rise of a lab breakthrough bed), or, later, a catalyst mixed with a sorbent. A layer of one solid is the special case with a single fraction of 1, and its results do not change.

## What a mixed layer is

The layer stores its **composition** as a map from adsorbent to mass fraction $w_k$, the mass of solid $k$ per mass of solid in the layer:

$$
\sum_k w_k = 1, \qquad 0 \le w_k \le 1 .
$$

A fraction may be a number or an expression over the deck's [global variables](global-variables.md) (`"x"` and `"1 - x"`), which is how a study varies one solid's share. The check is made on the evaluated values at run start, and the engine never rescales fractions for you. What `CheckLayerComposition` accepts and refuses:

| Layer composition | Result |
|---|---|
| fractions sum to 1 (to within $10^{-6}$), each in $[0, 1]$ | accepted (a fraction of 0 with the warning below) |
| a sum other than 1 | **error**: "mass fractions … sum to …; they must sum to 1" (the run does not start) |
| a negative or NaN fraction | **error**: "Mass fractions must be between 0 and 1" |
| a solid at mass fraction 0 | **warning**: it takes no part in the run; remove it from the composition |
| no solid, or an adsorbent that is not in the deck | **error** |
| a fraction expression that does not evaluate | **error**, with the evaluation message |
| every solid inert (an all-inert layer) | **accepted**: the layer is a spacer that holds gas and heat but adsorbs nothing |

An **inert solid** is an ordinary adsorbent whose components all have the `inert` isotherm (loading fixed at zero, see [Isotherms](isotherms.md#inert)). A new adsorbent starts that way for every component, so `sim.add_adsorbent("Glass beads", solid_density=2500, heat_capacity=840, thermal_conductivity=1.0, particle_void_fraction=0.0)` is an inert solid until you give it an isotherm. Its density, heat capacity and conductivity still enter the mixing rules below, so it dilutes the bed and adds thermal mass.

## Mixing rules

Each solid has its own solid density $\rho_k$, particle voidage $\varepsilon_{p,k}$, thermal conductivity $\lambda_k$ and specific heat capacity $C_{ps,k}$. The layer packs them at the **solid volume fractions**

$$
\phi_k = \frac{w_k / \rho_k}{\sum_j w_j / \rho_j}, \qquad \sum_k \phi_k = 1,
$$

and the effective properties of the solid are

$$
\rho_{mix} = \sum_k \phi_k \rho_k = \left(\sum_k \frac{w_k}{\rho_k}\right)^{-1}
\qquad
\varepsilon_p = \sum_k \phi_k\, \varepsilon_{p,k}
\qquad
\lambda_s = \sum_k \phi_k\, \lambda_k
\qquad
C_{ps} = \sum_k w_k\, C_{ps,k} .
$$

So the density is the **harmonic mean** weighted by mass fraction, the voidage and conductivity are weighted by **solid volume**, and the heat capacity, a per-mass property, by **mass fraction**. The conductivity is a plain volume average, an upper bound for a random packing of beads, so treat a strongly contrasting blend's axial heat conduction with care. Only the solids of the composition enter it: in a structured layer, a substrate that is not listed in the composition adds nothing to the solid conductivity.

The bed voidage $\varepsilon_b$, particle diameter and the pressure-drop, dispersion and heat-transfer models belong to the **layer**, not to a solid: a mixture shares one of each. The bulk density of the layer and of each solid in it (the mass of solid $k$ per bed volume) are

$$
\rho_b = (1-\varepsilon_b)\,\rho_{mix}, \qquad \rho_{b,k} = w_k\,\rho_b, \qquad \sum_k \rho_{b,k} = \rho_b .
$$

Each solid carries its **own loading** $q_{i,k}$ (mol per kg of solid $k$) in every cell, driven by its own isotherm and [rate model](kinetic-models.md), and the balances sum over the solids, each at its own bulk density. In the [mass balance](mass-balance.md) the adsorbed-phase term becomes $\sum_k \rho_{b,k}\,\partial q_{i,k}/\partial t$, and in the [energy balance](energy-balance.md) the loading heat capacity $\sum_k \rho_{b,k}\sum_i C_{pa,i}\,q_{i,k}$ and the heat of adsorption $\sum_k \rho_{b,k}\sum_i \Delta H^{\mathrm{eff}}_{i,k}\,\partial q_{i,k}/\partial t$ are summed the same way. $\Delta H^{\mathrm{eff}}_{i,k}$ is solid $k$'s stored heat of adsorption for component $i$, converted by the [heat of adsorption convention](energy-balance.md#heat-of-adsorption-convention) (the stored value itself for the default isosteric convention).

Each solid's isotherm and rate model see only that solid's loadings. A competitive isotherm, such as an extended Langmuir, divides by a denominator over the components on that one solid, so the solids do not compete with each other for adsorption sites: they interact only through the shared gas. A rate model that reads the solid's density, such as the capacity-factor macropore rate ([Matching Another Model](model-conventions.md#composition-used-by-the-rate-models)), uses that solid's own `solid_density` and its own equilibrium loading. An inert solid has $q_{i,k} = 0$, so it contributes thermal mass ($w_k C_{ps,k}$ in $C_{ps}$) and dilutes the bed (a smaller $\rho_{b,k}$ for the active solid), and nothing else.

!!! note "Volume fractions"
    Experimental papers on diluted beds often quote the *volume* share of the diluent. The engine stores *mass* fractions. With $\phi_k$ the solid-volume share of solid $k$ and $\rho_k$ its `solid_density`,

    $$
    w_k = \frac{\phi_k\,\rho_k}{\sum_j \phi_j\,\rho_j} .
    $$

    The mapping depends on what the paper's volume refers to (the solids only, or the packed bed including the voids between them) and on which density it uses for the particle, so check it against the paper before using it.

## Assumptions: one temperature, one gas

Within a cell all solids share **one solid temperature** (equal to the fluid temperature under local thermal equilibrium, or the single solid temperature of [LTNE](energy-balance.md)), and every solid takes up gas from the **same interparticle gas** (same pressure, composition and film conditions at the particle diameter of the layer). Solids do not exchange heat with each other as separate bodies: a hot adsorbent and a cold inert in the same cell are, by construction, at one temperature. The heat of a [reaction](reactions.md) on a catalyst in the blend is likewise released at that one solid temperature ([Energy Balance](energy-balance.md#heat-of-reaction)). This is accurate for a well-mixed bed of similar-sized particles, and tends to overstate the temperature damping of a diluent whose heat has to cross a poor contact.

## Mixture or thin alternating layers

Use a **mixture** when the solids really are intermingled at the particle scale, one particle size describes the packing, and you care about the average behaviour: a blended bed, or beads mixed through a bed. Use **separate layers** (one composition each, stacked along $z$) when the solids are segregated, when they differ in particle size, bed voidage or pressure-drop behaviour, or when the axial arrangement is itself the question (a guard layer, a layered desiccant and sorbent bed). Many thin alternating layers approximate a mixture only when each layer is a few cells and the cost of the extra interfaces is acceptable; they do not enforce a common temperature or gas within a pair, so for a truly well-mixed bed the mixture is both simpler and cheaper.

## Results for each solid

The solver log always holds the aggregate series, for example `CarbonDioxide Adsorbed Concentration`, a **mass-weighted sum over the solids**, $\sum_k w_k\, q_{i,k}$ per kg of solid in the layer; for a layer of one solid it is the solid's own loading. When a layer holds more than one solid the log also holds each solid's own series, in mol per kg of that solid: `<component> Adsorbed Concentration (<adsorbent>)`, and likewise for the equilibrium loading, the mass-transfer coefficient and, with a dual-rate model, the slow-site loading (`skarstrom.per_adsorbent_variables()`).

```python
result = sim.run(cycles=1)
sim.series_adsorbents()      # [{'id': 1, 'label': 'Activated Carbon', 'logged': True}, {'id': 2, 'label': 'Glass beads', ...}]
total = sim.series("adsorbed_concentration", "CarbonDioxide", z=0.1)                          # mass-weighted sum
carbon = sim.series("adsorbed_concentration", "CarbonDioxide", z=0.1, adsorbent="Activated Carbon")
beads = sim.series("adsorbed_concentration", "CarbonDioxide", z=0.1, adsorbent="Glass beads")  # all 0: inert
print(carbon.values[-1], 0.7 * carbon.values[-1], total.values[-1], beads.values[-1])
```

`adsorbent` takes a label from `sim.series_adsorbents()` (the name, with ` #<id>` added when two adsorbents share one), a deck ID or an adsorbent object such as `sim.adsorbents["Glass beads"]`; `sim.profile(..., adsorbent=)` works the same way. For a bed without a mixed layer nothing extra is logged and the per-adsorbent series is read off the aggregate. `skarstrom.series_manifest()` lists every bed series variable, its log name pattern, its unit and whether it is also logged per adsorbent.

```python
import skarstrom
per = [v["adsorbent_name"] for v in skarstrom.series_manifest()["variables"] if v["per_adsorbent"]]
print(per[0])
```

The CLI's `--csv out/` writes one file per series into `out/logs/`, with the per-adsorbent ones named after the log name, for example `CarbonDioxide_Adsorbed_Concentration__Activated_Carbon_.csv` and `CarbonDioxide_Adsorbed_Concentration__Glass_beads_.csv`.

## Setting a composition

!!! note "New in 0.9"
    `Layer.composition` and `Layer.set_composition()` are documented here for the first time, with the registry path below.

In Python a composition is given when the layer is created, or replaced on an existing layer (adsorbents by name, ID or object; a fraction is a number or an expression):

```python
beads = sim.adsorbents["Glass beads"]
layer = sim.layers[0]
print(dict(layer.composition))                          # adsorbent name -> mass fraction
layer.set_composition({"Activated Carbon": 0.5, "Glass beads": 0.5})
print(sim.get("layer[1].mass_fraction[2]"))             # the same value by its registry path
sim["layer[1].mass_fraction[2]"] = 0.3                  # one fraction: set the others to match before running
sim["layer[1].mass_fraction[1]"] = 0.7
```

When the layer is created, `sim.add_layer("Bed", composition={"13X": 0.7, "Beads": 0.3}, ...)` does the same. For a sweep, make the fractions expressions of one variable so that they sum to 1 for every value:

```python
sim.add_variable("w_carbon", 0.7, category="dimensionless")
layer.set_composition({"Activated Carbon": "w_carbon", "Glass beads": "1 - w_carbon"})
sim["variable[w_carbon]"] = 0.6
```

On the command line each fraction is the registry path `layer[<layer>].mass_fraction[<adsorbent>]`, repeated once per solid, and the fractions must sum to 1 when the run starts:

```bash
skarstrom examples run Illustrative_Diluted_Bed -q --json --set "layer[1].mass_fraction[1]=0.5" --set "layer[1].mass_fraction[2]=0.5"
```

## Example deck

`Illustrative_Diluted_Bed` is the [Ward 2022](examples.md) 100 ccm CO₂/He breakthrough on activated carbon with the adsorbent diluted by 30 % glass beads by mass in a single mixed layer, otherwise the same bed with the paper's values (it does not take the calibration of the `Ward2022_*` decks: ambient 291.4 K, q_s,2 6.04 mol/kg, the paper's dispersion and inside heat transfer coefficient). It is not a validated reproduction: the beads' density (2500 kg m⁻³), heat capacity (840 J kg⁻¹ K⁻¹) and conductivity (1.0 W m⁻¹ K⁻¹) are nominal soda-lime glass values. It runs in about a second. The carbon holds the same loading per kg as in the undiluted deck (the aggregate loading is 0.7 times as large), the beads hold none, and the breakthrough comes earlier because the bed holds 12 % less carbon per bed volume.
