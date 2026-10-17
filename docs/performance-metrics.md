# Performance Metrics

Skarstrom computes the standard PSA/TSA key performance indicators (KPIs) — purity, recovery, productivity, specific energy, and production rate — for a chosen product component over a chosen cycle. These are what you use to compare designs and to drive the [parametric study](parameterisation.md) and the [optimiser](parameterisation.md).

> KPIs are only meaningful at [cyclic steady state](cyclic-steady-state.md). Evaluate them on a converged cycle.

## Configuring KPIs

Before KPIs can be computed, the solver needs to know where fresh feed enters and where product is collected. This is the deck's **KPI configuration**, stored with the cycle:

1. Tick the component (for example *CarbonDioxide*) to **track** it. Only tracked components appear in the results.
2. The table lists every **endpoint** of every step: each stream of a Flow boundary, and each Pressure boundary, with the end of the column it sits on (Top or Bottom). Closed boundaries have no endpoints. Tick **Feed** and/or **Product** on each.

- **Feed** endpoints count the moles of the component that *enter* the bed through them.
- **Product** endpoints count the **net** moles that *leave* the bed through them: gas leaving counts positive, gas entering counts negative.

**Recycled product.** Because product endpoints are summed with their sign, ticking **Product** on an inlet that returns product to the bed, such as a purge or a [light-product pressurisation](light-product-pressurization.md) inlet, subtracts that recycled gas from the product. That gives a *net* recovery, $(\text{product} - \text{recycled})/\text{feed}$. Leaving the inlet unticked ignores it. If more gas is recycled than collected, the net total product can be zero or negative and **purity becomes undefined (shown as —) or exceeds 100 %**; read purity with only the product step ticked, as the [O₂ PSA tutorial](arora-2021-o2-psa.md) shows. The KPI helper in the workspace shows these definitions as equations.

## The metrics

For a product component, integrating the molar fluxes over the relevant steps of a cycle gives:

**Purity** — the fraction of the collected product that is the product component:

$$
\text{Purity} = \frac{n_{\text{product component, out}}}{n_{\text{total, out}}}
$$

**Recovery** — the net product component collected as a fraction of that component fed:

$$
\text{Recovery} = \frac{n_{\text{product, out}} - n_{\text{product, recycled}}}{n_{\text{product, fed}}}
$$

where the recycled term is the gas entering through endpoints ticked as Product. Recovery can be negative if more product is recycled than is collected, and above 100 % if the product endpoints include gas that was not counted as feed.

**Productivity (mass)** — net product mass produced per unit mass of adsorbent per unit time:

$$
\text{Productivity}_{\text{mass}} = \frac{n_{\text{product}}\, M}{t_{\text{cycle}}\; m_{\text{adsorbent}}} \quad \left[\frac{\text{kg}}{\text{kg}\cdot\text{s}}\right]
$$

**Productivity (volume)** — the same but per unit bed volume:

$$
\text{Productivity}_{\text{vol}} = \frac{n_{\text{product}}\, M}{t_{\text{cycle}}\; V_{\text{bed}}} \quad \left[\frac{\text{kg}}{\text{m}^3\cdot\text{s}}\right]
$$

**Energy Consumption** — the specific energy: total compression work over the cycle per unit mass of product:

$$
\text{Energy Consumption} = \frac{E_{\text{total}}}{n_{\text{product}}\, M} \quad \left[\frac{\text{J}}{\text{kg}}\right]
$$

where $E_{\text{total}}$ is the sum over all steps of the time-integrated compressor/pump power. Step [energy sources](energy-balance.md#energy-source-heating-cooling) (heaters and coolers, as in a TSA cycle) are **not** included; compute their duty per kilogram of product yourself, as in the [TSA tutorial](tsa-regeneration.md).

**Production Rate** — the net molar rate of product:

$$
\text{Production Rate} = \frac{n_{\text{product, out}} - n_{\text{product, recycled}}}{t_{\text{cycle}}} \quad \left[\frac{\text{mol}}{\text{s}}\right]
$$

where $M$ is the product molar mass, $t_{\text{cycle}}$ the cycle time, $m_{\text{adsorbent}}$ the total adsorbent mass, and $V_{\text{bed}}$ the total bed volume.

> **Energy Consumption requires compression power.** For the energy KPI to be non-zero, enable **Calculate pump power** (or use a pump device) on the [boundary conditions](boundary-conditions.md#pumps-compressors-and-vacuum-pumps) that do work — pressurisation, blowdown, evacuation — and set their efficiency and energy models.

## Reading the results

After a run, `sim.kpis()` returns one record per tracked component per logged cycle with **purity**, **recovery**, **productivity_mass**, **productivity_volume**, **energy_consumption** and **production_rate**; `sim.kpis(cycle=-1)` keeps the last cycle and `sim.kpi_table()` returns a pandas DataFrame. `skarstrom run deck.json` prints the KPIs of the last cycle, and `--csv DIR` writes every cycle to `kpis.csv`. Cycles are numbered from 0, so cycle 49 is the 50th, while the solver log counts from 1 (*CSS reached at cycle 50*). A KPI that is undefined for the run, for example a purity with no product drawn, is `nan` (`null` in JSON output).

To see how a KPI evolves as the process approaches cyclic steady state, plot it against the cycle number:

```python
import matplotlib.pyplot as plt
rows = [k for k in sim.kpis() if k.component == "CarbonDioxide"]
plt.plot([k.cycle + 1 for k in rows], [k.purity for k in rows])
```

In C++ the same numbers come from `adsim::Simulation::KPIs` (adsim/simulation.h) or the lower-level `KPIParams` (adsim/solver/kpis.h).
