# Tutorial: CO₂ VPSA with light-product pressurisation (Perez 2019)

Reference: [Estupiñan Perez, L., Sarkar, P., Rajendran, A. Experimental validation of multi-objective optimization techniques for design of vacuum swing adsorption processes. *Separation and Purification Technology* 224, 553–563 (2019).](https://doi.org/10.1016/j.seppur.2019.05.039)

**Goal.** Run a five-step vacuum swing adsorption (VPSA) cycle that captures CO₂ from a 15 % CO₂ / 85 % N₂ flue gas on zeolite 13X, see how **light-product pressurisation** (LPP) is wired with a step connection, run it to **cyclic steady state**, and compare purity and recovery with the experiment. The deck is experiment D1 of Table 3 in the paper.

**What you will get.** CSS after 182 cycles, a CO₂ purity of 93.7 % and a recovery of 90.2 %, against 97.1 % and 86.1 % measured. See [Validation](validation.md#perez-2019) for all six Perez experiments.

---

## 1. Load the deck

```python
import skarstrom

sim = skarstrom.examples.Perez2019_D1.load()
print(sim.components, sim.steps.names)
# ['CarbonDioxide', 'Nitrogen'] ['Feed 1', 'Feed 2', 'Blowdown', 'Evacuation', 'Pressurisation']
```

*Perez2019_E* and *Perez2019_F* use the same LPP cycle with other step times and pressures; *A*, *B* and *C1* are the basic four-step cycle without LPP. `examples/perez_2019.py` runs all six and compares them with the paper.

## 2. Walk through the setup

**Fluid.** *CarbonDioxide* and *Nitrogen*.

**Adsorbents.** *13X*: the paper's competitive (extended) dual-site Langmuir isotherm in concentration form, a macropore-controlled linear driving force (tortuosity 3, particle voidage 0.35) and heats of adsorption of 30.85 kJ/mol (CO₂) and 15.93 kJ/mol (N₂). The same parameters are in the library as [Zeolite 13X (Perez 2019)](adsorbent-library.md#zeolite-13x-perez2019).

**Column.** The laboratory column: 0.30 m long, 0.032 m diameter, a 3.1 mm steel wall and an ambient temperature of 293 K. The layer *13X* has 30 cells, a bed voidage of 0.41, 1.5 mm particles and the *Darcy (Kozeny–Carman)* momentum model. Dispersion uses the Ruthven correlation; heat transfer uses constant coefficients of 8.6 W m⁻² K⁻¹ (bed to wall) and 10 W m⁻² K⁻¹ (wall to ambient) with an axial thermal conductivity of 0.0903 W m⁻¹ K⁻¹.

**Initial Conditions.** Pure N₂ at 1 bar and 298 K.

**Cycle.** Five steps, 360 s per cycle (`skarstrom params Perez2019_D1.json` lists the step times, flows and pressures with their parameter paths):

| Step | Time | Bottom (feed end) | Top (light end) |
|------|------|-------------------|-----------------|
| Feed 1 | 27.5 s | **Flow** inlet *Feed 1 inlet*: 15 % CO₂, 0.0328 m/s superficial, 298 K | **Flow** outlet with two streams: a **Valve** to 1 bar (C<sub>v</sub> = 10) and the *LPP Stream*, connected to Pressurisation |
| Feed 2 | 117 s | **Flow** inlet *Feed 2 inlet*, as Feed 1 | **Flow** outlet: the valve to 1 bar only |
| Blowdown | 21 s | **Closed** | **Pressure** outlet, exponential ramp to 7.2 kPa |
| Evacuation | 167 s | **Pressure** outlet, exponential ramp to 3 kPa: the CO₂ product | **Closed** |
| Pressurisation | 27.5 s | **Closed** | **Pressure** inlet, linear ramp to 101.3 kPa, connected to Feed 1 |

Blowdown removes N₂ from the light end at an intermediate vacuum, so the evacuation product is rich in CO₂. The CO₂ product is drawn from the feed end, where the bed is most saturated.

**Connections.** One connection, between the *LPP Stream* on the top of Feed 1 and the top pressure boundary of Pressurisation. This is the light-product pressurisation:

- the **Pressurisation** inlet takes its composition and temperature from the gas leaving the top of the bed during Feed 1, so the bed is re-pressurised with N₂-rich product, not with the placeholder composition typed into the step;
- the **LPP Stream** on Feed 1 withdraws, instant by instant, the molar flow that the Pressurisation step drew, so the gas used for pressurisation is taken out of the product.

The two steps represent two beds working at the same time, which is why they have the same duration, 27.5 s. The coupling lags by one cycle: each side uses what its partner did in the previous cycle, and on the first cycle the typed-in values (pure N₂, zero LPP flow) stand in. See [Connecting Beds](connecting-beds.md) and [Light Product Pressurisation](light-product-pressurization.md).

**Cyclic steady state.** *Mass Balance (Per-Component)* (`css.method` = 3) with a tolerance of 0.005: the run stops once, for every component, the moles entering the column over a cycle match the moles leaving to within 0.5 %. The minimum is 5 cycles (`css.min_cycles`) and the cycle limit is 250 (`cycle.number_of_cycles`).

**KPI configuration.** CO₂ is tracked. **Feed**: *Feed 1 inlet* and *Feed 2 inlet*. **Product**: *Evacuation Pressure Boundary Bottom*.

## 3. Run

```python
result = sim.run()
print(result.status, result.css_cycle + 1, f"{result.wall_time_s:.0f} s")
# reached_css 182 31 s
```

The run stops at cyclic steady state, at cycle 182. It took about 30 s natively on one core of a shared 4-core Linux machine. `skarstrom run Perez2019_D1.json` does the same from the command line.

## 4. Interpret the results

**KPIs.** `sim.kpis(cycle=-1)` gives the last logged cycle, 181 (the 182nd cycle, numbered from 0):

```python
k = next(k for k in sim.kpis(cycle=-1) if k.component == "CarbonDioxide")
print(k.purity, k.recovery, k.production_rate, k.productivity_mass)
```

| KPI | Simulated | Experiment D1 | Repeat D2 | Paper's model |
|---|---|---|---|---|
| CO₂ purity | 93.72 % | 97.1 % | 94.5 % | 96.1 % |
| CO₂ recovery | 90.22 % | 86.1 % | 83.4 % | 86.8 % |

Production rate is 5.79 × 10⁻⁵ mol/s and the mass productivity 1.57 × 10⁻⁵ kg CO₂ per kg adsorbent per second. The energy KPI is 0 because no pump has power calculation switched on (see below).

The model is 3.4 points low on purity and 4.1 points high on recovery against D1. The test README lists two known causes: the engine's macropore LDF coefficient omits the capacity factor the paper uses, which makes the kinetics faster, and the Kozeny–Carman constant is 180 against the paper's 150.

**Approach to CSS.** Plot the KPI against the cycle number (only the last `settings.max_cycles_to_keep` = 10 cycles are logged):

```python
import matplotlib.pyplot as plt
rows = [k for k in sim.kpis() if k.component == "CarbonDioxide"]
plt.plot([k.cycle + 1 for k in rows], [k.purity for k in rows], "o-")
print(sim.css_history()[-1])        # the CSS check of the last cycle
```

**Pressure history.** The axial pressure profile at the end of each step of the last cycle:

```python
for i, name in enumerate(sim.steps):
    p = sim.profile("pressure", cycle=-1, step=i)       # p.values: (time, z)
    print(f"{name:15s} {p.values[-1].min() / 1e3:7.2f} - {p.values[-1].max() / 1e3:7.2f} kPa")
```

In our run the bed sat at 1.000–1.001 bar during the feed steps, reached 8.1–8.2 kPa at the end of Blowdown (the exponential ramp aims at 7.2 kPa but has not got there after 21 s), 3.05 kPa at the end of Evacuation, and was back at 1.013 bar after Pressurisation.

**Fronts.** `sim.profile("mole_fraction", "CarbonDioxide", cycle=-1, step=1)` holds the CO₂ profile through Feed 2 (`values[i]` is the profile at `time[i]`): plot a few rows to see how far the CO₂ front advances towards the top of the bed, then do the same for Evacuation (`step=3`) to see the bed being emptied through the feed end.

## 5. What to try next

- **Accelerate CSS.** Switch on the Aitken (Graves-Morris) extrapolation: `sim["settings.acceleration_method"] = 1`. In our run CSS came at cycle 176 instead of 182 with the same purity and recovery to four significant figures. The 100-cycle minimum limits the gain; lower `css.min_cycles` to let it stop sooner. See [Numerics](numerics.md#cyclic-steady-state-detection-and-acceleration).
- **Add the vacuum pump energy.** Pump power is calculated on boundaries that have a pump model. The Python registry exposes the pump settings of boundaries that already have one (`endpoint[N].pump.*`); switching power calculation on for the Evacuation and Blowdown boundaries is a deck edit (in C++, or in the deck JSON). With a discharge pressure of 1 bar, the default constant efficiency of 0.7 and a polytropic constant of 1.4, the specific energy came out at about **0.51 MJ per kg of CO₂**; purity and recovery do not change, because power calculation is post-processing only.
- **Switch off LPP** by disconnecting the pair (`adsim::DisconnectPair` in `adsim/objects/deck_edit.h`), and compare with the four-step decks *Perez2019_A* to *C1*.
- **Map the trade-off.** Sweep the Feed 2 time (`step[Feed 2].time`) or the evacuation end pressure (`step[Evacuation].bottom.pressure_exponential_end`) in a [parametric study](parameterisation.md), or [optimise](parameterisation.md#optimisation) recovery subject to purity ≥ 0.95, as the paper does with its multi-objective optimiser.
