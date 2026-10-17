# Tutorial: O₂ PSA (Arora 2021)

Reference: [Arora, A., Hasan, M.M.F. Flexible oxygen concentrators for medical applications. *Scientific Reports* 11, 14317 (2021).](https://doi.org/10.1038/s41598-021-93796-3)

**Goal.** Run a four-step Skarstrom cycle that makes oxygen from air on a LiLSX zeolite, learn how each step is wired, and read purity and recovery correctly when part of the product is sent back into the bed. The deck reproduces the PSA column of Table 4 in the paper (the case with at least 5 L/min and at least 85 % O₂).

**What you will see.** The O₂ purity of the production step comes out at 88.6 % against the paper's 84.7 %, but the net O₂ recovery is 4.2 % against the paper's 30.1 %. That recovery gap is an open disagreement, documented on the [Validation](validation.md#arora-hasan-2021) page. This tutorial shows where it comes from and what changes it.

---

## 1. Load the deck

```python
import numpy as np
import skarstrom

sim = skarstrom.examples.Arora2021_LiLSX_PSA.load()
print(sim.components, sim.steps.names)
# ['Oxygen', 'Nitrogen'] ['Production', 'Depressurization', 'Purge', 'Pressurization']
```

The companion deck *Arora2021_LiLSX_PVSA* is the vacuum version of the same cycle (desorption at 0.5 bar instead of 1 bar). `examples/arora_2021.py` runs both and compares them with the paper.

## 2. Walk through the setup

**Fluid.** Two components: *Oxygen* and *Nitrogen*; the feed is air at 21 % O₂.

**Adsorbents.** The deck carries three materials from Table 2 of the paper: *LiX Zeolite*, *LiLSX Zeolite* and *Zeolite 5A*. Only LiLSX is packed in the column; the other two are there so you can swap them in. *LiLSX Zeolite* has:

- a non-competitive **Dual Site Langmuir** isotherm for each component, in partial-pressure form. Nitrogen has the larger site energies, so it is the strongly adsorbed component and oxygen leaves as the light product;
- a **Linear Driving Force** rate model with constant coefficients, 104.07 s⁻¹ for O₂ and 42.84 s⁻¹ for N₂;
- heats of adsorption of 14.66 kJ/mol (O₂) and 25.45 kJ/mol (N₂).

The same material is also in the curated [adsorbent library](adsorbent-library.md#lilsx-arora2021) with its citations.

**Column.** A cylinder of 0.1 m diameter with one layer, *LiLSX Zeolite*: 0.127 m long, 30 cells, bed voidage 0.35 (from the paper's bulk density of 672.75 kg/m³ and particle density of 1035 kg/m³) and the particle diameter of 0.035 m as printed in the paper's SI. The momentum model is *Darcy (Kozeny–Carman)*, because this deck predates the Ergun option. Heat transfer to the wall uses constant coefficients of 70 W m⁻² K⁻¹ inside and out, with an axial thermal conductivity of 0.29 W m⁻¹ K⁻¹ (the paper's SI); the energy balance is non-isothermal with local thermal equilibrium.

**Cycle.** Four steps, 6.4 s per cycle, 50 cycles (`skarstrom params Arora2021_LiLSX_PSA.json` lists the flows and pressures with their parameter paths):

| Step | Time | Bottom (z = 0, feed end) | Top (z = L, product end) |
|------|------|--------------------------|--------------------------|
| Production | 3.0 s | **Flow** inlet: air (21 % O₂), 0.23 mol/s, 298 K | **Pressure** outlet at 3.5 bar: O₂ product |
| Depressurization | 1.0 s | **Pressure** outlet, linear ramp to 1 bar: N₂-rich vent | **Closed** |
| Purge | 1.4 s | **Pressure** outlet at 1 bar: waste | **Flow** inlet: 0.129 mol/s of 85 % O₂ product |
| Pressurization | 1.0 s | **Closed** | **Pressure** inlet, linear ramp to 3.5 bar with 85 % O₂ |

The purge flow is the paper's purge-to-feed velocity ratio (1.97) converted to a molar flow: 1.97 × 0.23 mol/s × (1 bar / 3.5 bar). The purge and pressurisation gas are typed in at 85 % O₂, the product specification; they are not connected to the production step, so the deck does not check that the product actually reaches that purity.

**Initial Conditions.** Air at 3.5 bar and 298 K.

**Cyclic steady state.** The method is *None*, so the run always does the full 50 cycles.

**KPI configuration.** Oxygen is tracked. The Feed is the air stream of Production. The Product boxes are ticked on three endpoints: the Production outlet, the Purge inlet and the Pressurization inlet. Product endpoints are summed as net outflow, so the two inlets *subtract* the O₂ sent back into the bed. That gives the net recovery the paper reports, but it also makes the total net product meaningless for purity (see below).

**Solver.** The deck caps the integrator step at 0.01 s (`settings.max_step`), kept from the validation test.

## 3. Run

```python
result = sim.run()
print(result.status, result.cycles_run, f"{result.wall_time_s:.0f} s")
# completed 50 6 s
```

The 50 cycles take about 6 s natively.

## 4. Interpret the results

The KPIs of the last logged cycle are cycle **49**: the solver numbers cycles from 0, so cycle 49 is the 50th. With the deck's KPI configuration:

```python
k = next(k for k in sim.kpis(cycle=-1) if k.component == "Oxygen")
print(k.purity, k.recovery, k.production_rate)      # nan 0.0424 0.00096
```

| KPI (cycle 49) | Value |
|---|---|
| Purity | `nan` (undefined) |
| Recovery | 4.24 % |
| Production rate | 0.96 × 10⁻³ mol/s |

**Why purity is undefined.** More gas goes back into the bed through the purge and pressurisation inlets than leaves through the production outlet, and that recycled gas is 15 % N₂ while the product is about 11 % N₂. The net nitrogen product is negative, the net total is not positive, and purity, the ratio of the two, is undefined. (In the PVSA deck the same configuration gives a "purity" above 100 %.) To read the purity of the oxygen you actually collect, integrate the production outlet over the last cycle from the stream logs:

```python
def moles(step, location, component=None):
    total = 0.0
    for e in sim.streams():
        if e["step"] == step and e["location"] == location:
            s = (sim.stream_series(step, e["name"], "molar_flow", component) if component
                 else sim.stream_series(step, e["name"], "molar_flow_rate"))
            total += np.sum(0.5 * (s.values[1:] + s.values[:-1]) * np.diff(s.time))
    return abs(total)

print(moles(0, "top", "Oxygen") / moles(0, "top"))  # 0.8856
```

That is a purity of **88.56 %** at cycle 49, against 84.72 % in the paper. Take recovery from the KPI configuration (net of the recycled oxygen, as the paper defines it) and purity from the collected product. The Validation page uses exactly these two numbers.

**Is it converged?** Not quite. The net recovery was 4.61 % at cycle 39 and 4.24 % at cycle 49, still falling by about 0.04 points per cycle:

```python
rec = [(k.cycle, k.recovery) for k in sim.kpis() if k.component == "Oxygen"]
```

Set a [cyclic steady state](cyclic-steady-state.md) method (`css.method`) and a larger number of cycles to run on to CSS.

**Product quality inside a step.** The O₂ mole fraction at the product end during *Production*:

```python
y = sim.stream_series(0, "Production Pressure Boundary Top", "mole_fraction", "Oxygen")
```

In the last cycle the product starts at 86 % O₂, rises to 92 % after about 2.4 s, then falls to 47 % at the end of the 3 s step as the nitrogen front breaks through. The production step is slightly too long for this feed rate.

**Where the recovery goes.** Each cycle feeds 0.23 mol/s × 3 s × 21 % ≈ 0.145 mol of O₂. The purge alone returns 0.129 mol/s × 1.4 s × 85 % ≈ 0.154 mol of O₂, more than the feed brings in, before counting pressurisation. With the purge flow fixed at the value derived from Table 4, almost all of the oxygen produced is spent regenerating the bed, which is why the net recovery is low. The paper's model reaches 30 %; we have not found the cause of the difference (the README of the test lists the input uncertainties, such as the production time t₁ = 3.0 s, which was not legible in the source table).

## 5. What to try next

- **Halve the purge.** Change the purge stream's molar flow from 0.1295 to 0.0647 mol/s: `sim["step[Purge].top[].molar_flow"] = 0.0647` (the stream has no name, hence `[]`; the ID path is `endpoint[5].molar_flow`). In our runs (50 cycles, same deck otherwise) the net recovery rose from 4.2 % to **34.1 %** and the production purity fell from 88.6 % to **73.0 %**: the classic purity–recovery trade-off. Sweep the purge flow in a [parametric study](parameterisation.md) to map the whole curve, or [optimise](parameterisation.md#optimisation) recovery subject to purity ≥ 0.85.
- **Shorten production** to stop before the nitrogen front breaks through, and watch the production-step purity.
- **Connect the recycle streams.** Replace the typed-in 85 % O₂ purge and pressurisation gas with [connections](connecting-beds.md) to the production outlet (the [light reflux](light-reflux-recycle.md) and [light product pressurisation](light-product-pressurization.md) pages show how), so the recycled gas has the composition the bed really produces.
- **Swap the adsorbent.** Change the layer's composition to *LiX Zeolite* or *Zeolite 5A* (a deck edit: `MixedAdsorbentLayer::composition_` in C++) and compare.
- **Try the vacuum cycle.** *Arora2021_LiLSX_PVSA* desorbs at 0.5 bar. The validation test reports a production purity of 94.06 % and a net recovery of 19.88 % for it, against 89.97 % and 65.51 % in the paper.
