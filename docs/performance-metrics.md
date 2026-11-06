# Performance Metrics

Skarstrom computes the standard PSA/TSA key performance indicators (KPIs) — purity, recovery, productivity, specific electrical energy, and production rate — for a chosen product component over a chosen cycle. These are what you use to compare designs and to drive the [parametric study](parameterisation.md) and the [optimiser](parameterisation.md).

> KPIs are only meaningful at [cyclic steady state](cyclic-steady-state.md). Evaluate them on a converged cycle.

## Configuring KPIs

Before KPIs can be computed, the solver needs to know where fresh feed enters and where product is collected. This is the deck's **KPI configuration**, stored with the cycle:

1. Tick the component (for example *CarbonDioxide*) to **track** it. Only tracked components appear in the results.
2. The table lists every **endpoint** of every step: each stream of a Flow boundary, and each Pressure boundary, with the end of the column it sits on (Top or Bottom). Closed boundaries have no endpoints. Tick **Feed** and/or **Product** on each.

- **Feed** endpoints count the moles of the component that *enter* the bed through them.
- **Product** endpoints count the **net** moles that *leave* the bed through them: gas leaving counts positive, gas entering counts negative.
- **Recycle** endpoints are inlets that return product to the bed, such as a purge or a [light-product pressurisation](light-product-pressurization.md) inlet fed from the product. They count the moles of the component that *enter* the bed through them.

**Recycled product.** Recycled gas is product the process uses itself, so it is subtracted from the product for recovery, productivity, production rate and specific energy: a *net* recovery, $(\text{product} - \text{recycled})/\text{feed}$. Purity is the purity of the gas collected at the product endpoints and does not depend on the recycle. In Python, pass the recycle endpoints to `sim.kpi.track(..., recycle=[...])`; in a deck they are the `Recycle_Endpoints` of the KPI configuration. The browser interface does not show the recycle role yet.

Do not tick a recycle inlet as **Product** as well: the signed sum would subtract it from the purity too, and when more gas is recycled than collected the "purity" becomes undefined (shown as —) or exceeds 100 %.

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

where the recycled term is the product component entering through the recycle endpoints. Recovery can be negative if more product is recycled than is collected, and above 100 % if the product endpoints include gas that was not counted as feed.

**Productivity (mass)** — net product mass produced per unit mass of adsorbent per unit time:

$$
\text{Productivity}_{\text{mass}} = \frac{n_{\text{product}}\thinspace M}{t_{\text{cycle}}\thickspace m_{\text{adsorbent}}} \quad \left[\frac{\text{kg}}{\text{kg}\cdot\text{s}}\right]
$$

**Productivity (volume)** — the same but per unit bed volume:

$$
\text{Productivity}_{\text{vol}} = \frac{n_{\text{product}}\thinspace M}{t_{\text{cycle}}\thickspace V_{\text{bed}}} \quad \left[\frac{\text{kg}}{\text{m}^3\cdot\text{s}}\right]
$$

**Electrical Energy** (`electrical_energy`, also `energy_consumption`) — the specific electrical energy: total compressor and vacuum-pump work over the cycle per unit mass of net product:

$$
\text{Electrical Energy} = \frac{E_{\text{total}}}{n_{\text{product}}\thinspace M} \quad \left[\frac{\text{J}}{\text{kg}}\right]
$$

where $E_{\text{total}}$ is the sum over all steps of the time-integrated compressor/pump power. Step [energy sources](energy-balance.md#energy-source-heating-cooling) (heaters and coolers, as in a TSA cycle) are **not** included: their heat is the thermal energy and cooling duty below.

`energy_consumption` is the original name of this KPI and keeps its meaning: it is the electrical energy, the same value as `electrical_energy`, in KPI records, sensitivity outputs (`Component.energy_consumption` or `Component.electrical_energy`), optimiser objectives and the cost estimate's electricity cost. Decks and scripts that use it are unchanged. It is undefined (`nan`) exactly when no net product is drawn.

**Thermal Energy** (`thermal_energy`) — the specific thermal energy: heat supplied to the bed over the cycle per unit mass of net product:

$$
\text{Thermal Energy} = \frac{\sum_{s} Q^{+}_{s}\thinspace t_s}{n_{\text{product}}\thinspace M} \quad \left[\frac{\text{J}}{\text{kg}}\right]
$$

where $Q^{+}_{s}$ is the total power of the step's heating [energy sources](energy-balance.md#energy-source-heating-cooling) (those with a positive power) and $t_s$ the step time. A heating medium's duty will be added to it once the engine has one.

**Cooling Duty** (`cooling_duty`) — the heat removed by the cooling energy sources (negative power) over the cycle per unit mass of net product, reported as a **positive** magnitude with the same formula using $|Q^{-}_{s}|$. It is for information: it is not subtracted from the thermal energy.

Both use the source powers the solver applied, so a disabled source counts nothing, and both are 0 in an isothermal run, where the sources do nothing. They are defined and undefined exactly as the electrical energy (`nan` when no net product is drawn). Each record also has `step_thermal_energy` and `step_cooling_duty`, the share of every step in the deck's step order (they sum to the cycle values). The electrical energy, thermal energy and cooling duty are separate figures; none includes another.

**Production Rate** — the net molar rate of product:

$$
\text{Production Rate} = \frac{n_{\text{product, out}} - n_{\text{product, recycled}}}{t_{\text{cycle}}} \quad \left[\frac{\text{mol}}{\text{s}}\right]
$$

where $n_{\text{product}} = n_{\text{product, out}} - n_{\text{product, recycled}}$ is the net product, $M$ is the product molar mass, $t_{\text{cycle}}$ the cycle time, $m_{\text{adsorbent}}$ the total adsorbent mass, and $V_{\text{bed}}$ the total bed volume.

> **Electrical Energy requires compression power.** For the electrical energy KPI to be non-zero, enable **Calculate pump power** (or use a pump device) on the [boundary conditions](boundary-conditions.md#pumps-compressors-and-vacuum-pumps) that do work — pressurisation, blowdown, evacuation — and set their efficiency and energy models.

### Energy: the three figures and their units

| KPI | What it counts | Field |
|---|---|---|
| Electrical energy | compressor and vacuum-pump work | `electrical_energy` (also `energy_consumption`) |
| Thermal energy | heat supplied by the heating energy sources | `thermal_energy` |
| Cooling duty | heat removed by the cooling energy sources (positive) | `cooling_duty` |

All three are in **J per kg of net product**. They are separate: add them only if your comparison means to (electricity and heat are not worth the same, so a plant study usually weighs them, for example by converting heat to primary energy or cost). To convert:

| To | From J/kg | Note |
|---|---|---|
| MJ/kg | divide by 10⁶ | |
| GJ/t (GJ per tonne) | divide by 10⁶ | 1 GJ/t = 10⁹ J / 10³ kg = 1 MJ/kg, so the number is the same as in MJ/kg |
| kWh/t | divide by 3600 | 1 kWh/t = 3.6 × 10⁶ J / 10³ kg = 3600 J/kg |

For example, the [TSA tutorial](tsa-regeneration.md) heater gives a thermal energy of 6.07 × 10⁶ J/kg of CO₂: 6.07 MJ/kg, 6.07 GJ/t, or 1685 kWh/t.

```python
k = sim.kpis(cycle=-1)[0]
print(f"electrical {k.electrical_energy / 1e6:.3f} GJ/t, thermal {k.thermal_energy / 1e6:.3f} GJ/t, "
      f"cooling {k.cooling_duty / 1e6:.3f} GJ/t")
```

## Reading the results

After a run, `sim.kpis()` returns one record per tracked component per cycle of the run (every cycle, not only the logged ones: the solver records the KPIs when each cycle completes) with **purity**, **recovery**, **productivity_mass**, **productivity_volume**, **energy_consumption**, **production_rate**, **electrical_energy** (the same value as energy_consumption), **thermal_energy**, **cooling_duty**, **step_thermal_energy** and **step_cooling_duty**; `sim.kpis(cycle=-1)` keeps the last cycle and `sim.kpi_table()` returns a pandas DataFrame. `skarstrom run deck.json` prints the KPIs of the last cycle, and `--csv DIR` writes every cycle to `kpis.csv` (its `energy [J/kg]` column is the electrical energy; `electrical energy [J/kg]`, `thermal energy [J/kg]` and `cooling duty [J/kg]` columns are appended). Cycles are numbered from 0, so cycle 49 is the 50th, while the solver log counts from 1 (*CSS reached at cycle 50*). A KPI that is undefined for the run, for example a purity with no product drawn, is `nan` (`null` in JSON output).

To see how a KPI evolves as the process approaches cyclic steady state, plot it against the cycle number:

```python
import matplotlib.pyplot as plt
rows = [k for k in sim.kpis() if k.component == "CarbonDioxide"]
plt.plot([k.cycle + 1 for k in rows], [k.purity for k in rows])
```

In C++ the same numbers come from `adsim::Simulation::KPIs` (adsim/simulation.h) or the lower-level `KPIParams` (adsim/solver/kpis.h).
