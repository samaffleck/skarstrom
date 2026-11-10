# Case Study: a CO₂ VSA Under Plant Constraints (Base Decks)

!!! note "New in 0.9"
    The decks `Perez2019_CaseStudy_Base` and `Perez2019_CaseStudy_Base_LPP` and the builder script `examples/perez_2019_case_study.py`.

This case study starts from the validated design of one CO₂/N₂ VSA paper, [Estupiñan Perez et al. (2019)](perez-2019-vpsa.md), adds the engineering constraints a real plant must meet, and re-optimises under them. This page describes its first step: the **base decks**. They run the paper's cycles at the paper's operating points, but the way a plant runs them:

- **Steps are driven by devices, not imposed profiles.** A blower delivers the feed, the light product leaves through a valve, the bed is pressurised through a valve that opens on a ramp, and one vacuum pump does the blowdown and the evacuation.
- **The cycle is an n-bed schedule.** The steps are timed as fractions of one cycle time, so the feed and the vacuum pump both run without a break.
- **Each constraint is one parameter.** Every device setting and every constraint value is a [global variable](global-variables.md), so a [study](sensitivity.md) or an [optimiser](parameterisation.md#optimisation) can change it by its registry path.

| Deck | Starts from | Cycle |
|---|---|---|
| `Perez2019_CaseStudy_Base` | `Perez2019_A` (four-step cycle) | Feed, Blowdown, Evacuation, Pressurisation from the feed header |
| `Perez2019_CaseStudy_Base_LPP` | `Perez2019_D1` (light-product pressurisation) | Feed 1, Feed 2, Blowdown, Evacuation, Pressurisation (LPP), Idle |

The column, the adsorbent, the heat transfer, the initial state and the solver settings are the validated decks' (`examples/perez_2019.py`). The momentum balance is Ergun (the validated decks use Kozeny–Carman), so the blower and the pump work against the full pressure drop.

```python
import skarstrom

sim = skarstrom.examples.Perez2019_CaseStudy_Base.load()
print(sim.steps.names, sim.cycle.total_time)
result = sim.run(cycles=3)
print(result.status, result.healthy)
```

`python examples/perez_2019_case_study.py` runs both decks to cyclic steady state and prints the comparison below.

---

## The schedule

The decks use fraction step-time mode: `cycle.total_time` and one fraction per step ([Cycles](cycles.md#step-timing)). The fractions are expressions of the variable `n_beds`, so the schedule follows the number of beds:

| Rule | Four-step deck | LPP deck |
|---|---|---|
| Continuous feed: the feed takes 1/n of the cycle | Feed: `1 / n_beds` | Feed 1: `feed1_share / n_beds`, Feed 2: `(1 - feed1_share) / n_beds` |
| Continuous vacuum pump: one bed is always in blowdown or evacuation | Blowdown: `blowdown_share / n_beds`, Evacuation: `(1 - blowdown_share) / n_beds` | the same |
| The LPP step is as long as the partner's first feed step | — | Pressurisation: `feed1_share / n_beds` |
| The rest of the cycle | Pressurisation: `1 - 2 / n_beds` | Idle (both ends closed): `1 - (2 + feed1_share) / n_beds` |

The cycle time is `n_beds * t_feed`. The feed time and velocity are the paper's, so with `n_beds` = 3 the cycle is 297 s (four-step) and 433.5 s (LPP). With fewer than three beds there is no time left for pressurisation.

**Where the single-bed model is an approximation.** Skarstrom solves one bed and represents its partners with [connections](connecting-beds.md):

- The LPP valve connects the top of Feed 1 with the top of Pressurisation. Each side reads the other bed's pressure from the previous cycle, so the two flows agree only at cyclic steady state.
- The feed header is held at a fixed pressure (`p_feed_header`). In a plant the blower also supplies the pressurisation gas while it feeds another bed. The feed flow does not respond to that extra draw here.
- The vacuum pump runs at all times, but the pump model has no inertia and no shared suction header. Blowdown and evacuation of different beds never overlap in this schedule, so this is exact for the schedule, not for a pump that serves two beds at once.
- The Idle step of the LPP deck is a closed bed. In a plant it is the bed waiting for its feed slot.
- The vent valves have no check valve. With a larger LPP valve the LPP draw exceeds what the feed bed produces and gas flows back in through the Feed 1 vent (the validated `Perez2019_D1` deck draws 6 % of its feed back this way). `cv_lpp` is chosen just below that point. Turning on the check valve (`endpoint[N].valve.check_valve`) in this single-bed model needs a feed that is a molar flow, not a velocity at the bed pressure; that is left for the optimisation step.

## The devices

| Step | Device | Model |
|---|---|---|
| Feed | Positive-displacement blower | A set feed velocity (`feed_velocity`), with the [compressor power](boundary-conditions.md#pumps-compressors-and-vacuum-pumps) calculated from the bed inlet pressure over the suction pressure `p_atm` (efficiency `blower_efficiency`) |
| Feed (top) | Vent valve to the atmosphere | A [valve stream](boundary-conditions.md#valves) to `p_atm`, C<sub>v</sub> = `cv_vent` |
| Blowdown (top), Evacuation (bottom) | One vacuum pump | A pump stream that discharges to `p_atm`. Its shaft power falls linearly with the suction pressure P, W = W<sub>max</sub> (1 − P / p<sub>atm</sub>), with W<sub>max</sub> = `vacuum_pump_power`. The efficiency curve is η = η<sub>max</sub> α p / (1 + α p), p in bar (`vp_eta_max`, `vp_alpha`) |
| Pressurisation (four-step) | Valve from the feed header | A valve stream from `p_feed_header`, C<sub>v</sub> = `cv_pressurisation`, that opens on a ramp (`valve_ramp_start`, `valve_ramp_time`; below) |
| Pressurisation (LPP) | Valve between the beds | A valve stream at each end, C<sub>v</sub> = `cv_lpp`, connected so that each reads the other bed's pressure. Both ends open on the same ramp |

The power falls to zero as the pressure ratio goes to 1. Near the atmosphere the pump then moves a finite volume flow, as a constant-speed pump does. With a constant power, the flow at the start of blowdown is unbounded: the pressure ratio is 1, so the power can move any flow. With a constant power the approach to fluidisation reached 19–31 at the start of blowdown. With the power curve it stays below 1 in blowdown and evacuation.

### The ramped valve opening

!!! note "New in 0.9"
    The [valve opening ramp](boundary-conditions.md#opening-ramp) and the variables `valve_ramp_start` and `valve_ramp_time` of both decks.

A pressurisation valve that opens fully at the start of the step sees the largest pressure difference it will ever see: the header is at 1 bar and the bed has just been evacuated to 4 kPa. Its flow is largest in that first instant, and it enters a bed full of thin gas, so the pressure gradient at the feed end is steep. With the four-step deck's valve fully open at once, the approach to fluidisation reaches 3.6 at the feed end for a moment. A smaller valve does not solve this: the approach falls below 1 only for C<sub>v</sub> below about 0.001, and that valve then cannot pressurise the bed in its step. An earlier search of 844 runs without the ramp found no design that meets the fluidisation constraint and still pressurises the bed.

A plant valve does not open in an instant. Actuated valves have a stroke time, and pressurisation valves are often given a slow opening to protect the bed. The decks now model this with the valve's [opening ramp](boundary-conditions.md#opening-ramp): the valve opens from `valve_ramp_start` of its C<sub>v</sub> to fully open over `valve_ramp_time` from the start of the step, along a smooth cubic. The first instant then sees a small valve, and the full valve acts only once the bed pressure, and with it the gas density at the feed end, has risen.

| Deck | Valve | `valve_ramp_start` | `valve_ramp_time` | Effect at the base point |
|---|---|---|---|---|
| Four-step | Feed header | 0.1 | 20 s | Approach to fluidisation in pressurisation 3.62 → 0.95. Pressure drop 14.1 → 3.7 kPa m⁻¹. The bed still comes within 2 % of 1 bar after 31.7 s of its 99 s step |
| LPP | Both ends of the LPP valve | 0.2 | 8 s | Pressure drop in the LPP step 16.2 → 5.7 kPa m⁻¹ (its approach to fluidisation was already 0.3: the LPP enters from the top). The bed ends the step at 99.84 kPa, 1.4 % below the validated deck |

The start opening sets the peak: at the base C<sub>v</sub>, an opening of 0.1 gives 0.95 and 0.2 gives 1.39, whatever the ramp time from 20 to 30 s. The ramp time sets how soon the bed fills. Both are decision variables of the optimisation.

## Matching the paper's operating points

The device settings were chosen so that the device-driven steps reach the pressures that the validated decks reach at cyclic steady state. These are not the nominal end values of the validated decks' imposed profiles, because those profiles do not reach their end values in the step time. The pump power and the blowdown share were solved together for the two vacuum pressures. The pressurisation valve was chosen to reach 1 bar within the paper's 30 s, and its opening ramp was then added (the pump power and the share were solved again with the ramp in place). The vent valve keeps the feed step within 0.7 kPa of 1 bar. The tolerance is 2 % on the step-end pressures.

| Setting | Four-step | LPP | Why |
|---|---|---|---|
| `vacuum_pump_power` (W<sub>max</sub>) | 10.2 W | 5.85 W | Blowdown and evacuation end pressures (solved with the share) |
| `blowdown_share` | 0.205 | 0.233 | as above |
| `vp_eta_max`, `vp_alpha` | 0.8, 19.55 | 0.8, 19.55 | The engine's default vacuum-pump curve |
| `cv_vent` | 0.5 mol s⁻¹ bar⁻¹ | 0.5 | Feed-step pressure within 0.7 kPa of 1 bar |
| `cv_pressurisation` | 0.006 mol s⁻¹ bar⁻¹ | — | 1 bar within the paper's 30 s |
| `cv_lpp` | — | 0.008 mol s⁻¹ bar⁻¹ | The largest valve that does not draw gas back in through the Feed 1 vent |
| `valve_ramp_start`, `valve_ramp_time` | 0.1, 20 s | 0.2, 8 s | Approach to fluidisation about 1 or less and a lower pressure drop, with the bed still pressurised in the step (above) |
| `blower_efficiency` | 0.7 | 0.7 | Typical. The blower uses less than 1 % of the energy |

At cyclic steady state (this engine version, 80 cells; specific energy = vacuum-pump and blower work per kg of CO₂, with the same pump curve added to the validated decks for comparison):

| Quantity | `Perez2019_A` | Base | `Perez2019_D1` | Base LPP |
|---|---|---|---|---|
| Cycle time | 205 s | 297 s | 360 s | 433.5 s |
| Blowdown: time / end pressure | 19 s / 7.21 kPa | 20.3 s / 7.19 kPa (−0.2 %) | 21 s / 8.16 kPa | 33.7 s / 8.08 kPa (−1.0 %) |
| Evacuation: time / end pressure | 57 s / 4.00 kPa | 78.7 s / 4.01 kPa (+0.4 %) | 167 s / 3.11 kPa | 110.8 s / 3.12 kPa (+0.0 %) |
| Pressurisation: time / end pressure | 30 s / 99.98 kPa | 99 s / 100.00 kPa (within 2 % after 31.7 s) | 27.5 s / 101.30 kPa | 27.5 s / 99.84 kPa (−1.4 %, within 2 % after 19.6 s) |
| CO₂ purity | 99.55 % | 99.57 % | 94.08 % | 94.24 % |
| CO₂ recovery | 60.1 % | 59.5 % | 89.5 % | 89.4 % |
| Productivity, 10⁻⁵ kg kg⁻¹ s⁻¹ | 3.68 | 2.53 | 1.56 | 1.29 |
| Specific energy, MJ kg⁻¹ | 0.752 | 0.743 | 0.821 | 0.810 |
| Cyclic steady state | cycle 48 | cycle 49, 43 s | cycle 92 | cycle 101, 59 s |

Purity, recovery and specific energy agree with the validated decks to within 0.6 points and 1.4 %. The productivity is lower because the n-bed schedule makes the cycle longer: the feed time is the paper's, but the evacuation and pressurisation get a full third of the cycle each. The LPP pressurisation ends just below the feed bed's pressure, at 99.8 kPa. The validated deck pressurises to 1 atm while its feed vents to 1 bar.

## The constraint parameters

These registry paths are the settings the case study varies:

| Path | Meaning | Base value |
|---|---|---|
| `variable[n_beds]` | Beds in the schedule | 3 |
| `variable[t_feed]` | Feed time, 1/n of the cycle | 99 s / 144.5 s |
| `variable[feed_velocity]` | Blower flow, as a superficial velocity | 0.0849 / 0.0328 m s⁻¹ |
| `variable[blowdown_share]` | Share of the pump's time spent in blowdown | above |
| `variable[vacuum_pump_power]` | Pump power at full vacuum, W | above |
| `variable[vp_eta_max]`, `variable[vp_alpha]` | Pump efficiency curve | 0.8, 19.55 |
| `variable[blower_efficiency]` | Blower efficiency | 0.7 |
| `variable[cv_vent]`, `variable[cv_pressurisation]`, `variable[cv_lpp]` | Valve coefficients, mol s⁻¹ bar⁻¹ | above |
| `variable[valve_ramp_start]`, `variable[valve_ramp_time]` | Pressurisation valve opening ramp: share of C<sub>v</sub> open at the start of the step, time to open fully | 0.1, 20 s / 0.2, 8 s |
| `variable[p_atm]`, `variable[p_feed_header]` | Atmosphere and feed-header pressure | 1 bar |
| `variable[feed1_share]` | Share of the feed time in Feed 1, the LPP supply (LPP deck) | 0.190 |
| `variable[p_evac_floor]` | Constraint: lowest allowed evacuation pressure | 10 kPa |
| `variable[purity_target]`, `variable[recovery_target]` | Constraint: CO₂ purity and recovery targets | 0.95, 0.90 |

The last three are constraint values. The run does not use them: the base decks reproduce the paper, whose evacuation goes below a 10 kPa floor. The optimisation step of the case study checks them.

```python
sim = skarstrom.examples.Perez2019_CaseStudy_Base_LPP.load()
for name in ("n_beds", "vacuum_pump_power", "cv_lpp", "p_evac_floor"):
    print(name, sim.get(f"variable[{name}]"))
sim["variable[n_beds]"] = 4          # a four-bed schedule: every fraction follows
print([round(s.fraction, 4) for s in sim.steps], sim.cycle.total_time)
```

The same from the command line:

```bash
skarstrom examples run Perez2019_CaseStudy_Base --json -q > base.json
```

---

## Re-optimising under constraints

!!! note "New in 0.9"
    The script `examples/perez_2019_case_study_optimise.py`.

`examples/perez_2019_case_study_optimise.py` searches the device settings and the schedule of the two base decks with NSGA-II ([pymoo](https://pymoo.org), `pip install pymoo`), in the way of `examples/optimise_nsga2.py` ([Optimisation](parameterisation.md#multi-objective-optimisation-nsga-ii)). Each candidate is one run to cyclic steady state. The candidates of a generation run in parallel processes.

<!-- no-run: needs a repository checkout and pymoo; each generation takes minutes -->
```bash
pip install pymoo
python examples/perez_2019_case_study_optimise.py --deck Base --base-only        # the base deck's own point
python examples/perez_2019_case_study_optimise.py                                 # energy vs productivity, both decks
python examples/perez_2019_case_study_optimise.py --objectives purity             # purity vs recovery, plant limits
python examples/perez_2019_case_study_optimise.py --objectives purity --unconstrained
python examples/perez_2019_case_study_optimise.py --pop 24 --gen 20 --workers 8 --seed 2 --out runs/
```

### Decision variables

The variables are registry paths of the decks' [global variables](global-variables.md). The bounds are in `VARIABLES` at the top of the script, with the reason for each. `n_beds` stays 3, so `variable[t_feed]` sets the cycle time (`n_beds * t_feed`).

| Variable | Four-step deck | LPP deck |
|---|---|---|
| `variable[t_feed]` | 40–200 s | 60–250 s |
| `variable[feed_velocity]` (blower flow, superficial) | 0.02–0.15 m s⁻¹ | 0.015–0.1 m s⁻¹ |
| `variable[vacuum_pump_power]` (W<sub>max</sub>) | 1–30 W | 0.5–20 W |
| `variable[blowdown_share]` | 0.08–0.6 | 0.08–0.6 |
| `variable[cv_pressurisation]` | 0.0003–0.02 mol s⁻¹ bar⁻¹ | — |
| `variable[cv_lpp]`, `variable[feed1_share]` | — | 0.002–0.016, 0.08–0.5 |
| `variable[valve_ramp_start]` (pressurisation valve) | 0.02–1 | 0.02–1 |
| `variable[valve_ramp_time]` | 1–60 s | 1–30 s |
| `variable[cv_vent]` | 0.05–2 mol s⁻¹ bar⁻¹ | 0.05–2 |

The base decks have no pressure set points: the pump power and the blowdown share set the blowdown and evacuation pressures, and the valves set the feed and pressurisation pressures. The pump power, the valve coefficients and the valve's start opening are searched in log10 (`LOG_SCALE`). Each spans one to two decades, and the plant constraints hold only near the small end, which a uniform search of a few generations seldom samples.

### Objectives

The script has two objective sets (`--objectives`). Each set is a separate run, because they answer different questions. With all four objectives in one run, the front is a four-dimensional surface that a small population cannot cover, and the purity and recovery targets would be objectives and constraints at the same time.

- `energy` (default): minimise the specific electrical energy, in kWh per tonne of CO₂ (the `electrical_energy` KPI: vacuum-pump and blower work per kg of CO₂, divided by 3600), and maximise the productivity, in kg of CO₂ per tonne of adsorbent per hour. All the constraints below apply. This is the design question.
- `purity`: maximise CO₂ purity and recovery, with the plant constraints but not the purity and recovery targets. This is the front that the paper optimised. With `--unconstrained`, only cyclic steady state and a healthy run are required: this is the paper's search with the devices of the case study.

### Constraints

A candidate is feasible when every margin is zero or more. Each margin is normalised by its bound. NSGA-II ranks infeasible candidates by their summed violation, so the search moves towards the feasible region before it has found it.

| Constraint | Value checked | Source |
|---|---|---|
| Purity ≥ `variable[purity_target]` (0.95) | CO₂ purity of the last cycle | KPI record (`energy` only) |
| Recovery ≥ `variable[recovery_target]` (0.90) | CO₂ recovery of the last cycle | KPI record (`energy` only) |
| No fluidisation | The largest approach to fluidisation, every packed face, every step, ≤ 1 | `sim.operating_limits()` |
| Bed pressure drop | The largest \|P<sub>top</sub> − P<sub>bottom</sub>\| / L over the logged times of every step ≤ `--max-dp-per-m` (10 kPa m⁻¹) | `sim.profile("pressure", ...)` |
| Evacuation floor | The lowest pressure in the column in blowdown and evacuation ≥ `variable[p_evac_floor]` (10 kPa) | `sim.operating_limits()` (`min_pressure`) |
| Product | CO₂ recovery ≥ 1 %. Without it a weak pump that evacuates nothing meets every plant limit, and NSGA-II, which ranks any feasible candidate above any infeasible one, fills its population with runs that make nothing | KPI record |
| Converged | The run reached cyclic steady state and `result.healthy` is true | `RunResult` |
| No outlet backflow (`--backflow-constraint`) | Gas in through an unconnected outlet ≤ 0.1 % of the cycle feed | `sim.operating_limits()` (`outlet_backflow` warning) |

The operating-limit values are those at the logged times. A peak between two log points is not seen. The backflow is always recorded. It is not enforced by default: the LPP base deck itself draws 12.4 % of its cycle feed back in through the Feed 1 vent at cyclic steady state, so every LPP candidate would be infeasible. That is a known issue of the connected LPP valves ([#241](https://github.com/samaffleck/skarstrom/issues/241)), not of the deck; the backflow becomes a default constraint when it is fixed. The vacuum floor is a constraint on the result, not a limit of the pump model.

The base points are infeasible. The table gives their values (`--base-only`):

| Constraint | `Perez2019_CaseStudy_Base` | `Perez2019_CaseStudy_Base_LPP` |
|---|---|---|
| Purity ≥ 95 % | 99.6 % | 94.2 % (fails) |
| Recovery ≥ 90 % | 59.5 % (fails) | 89.4 % (fails) |
| Approach to fluidisation ≤ 1 | 0.98, evacuation (0.95 in pressurisation; 3.62 without the ramp) | 0.30 |
| Pressure drop ≤ 10 kPa m⁻¹ | 3.7 kPa m⁻¹, pressurisation (14.1 without the ramp) | 5.7 kPa m⁻¹, LPP from the top (16.2 without the ramp) |
| Lowest vacuum pressure ≥ 10 kPa | 3.9 kPa (fails) | 3.1 kPa (fails) |
| Outlet backflow (recorded) | 0 | 12.4 % of the feed |

In the feed step the pressure drop is about 0.6 kPa m⁻¹. The pressure-drop bound acts on the pressurisation transients, where a valve opens onto an evacuated bed. With the opening ramps, both base points meet the fluidisation and pressure-drop constraints. They fail the vacuum floor and the targets.

### Warm starts, parallel runs and resuming

- **Warm starts.** Each worker process keeps a `skarstrom.optimisation.StateCache`. Each candidate starts from the cyclic steady state of the nearest candidate that the worker has run (`skarstrom.optimisation.run_warm`, distance in the variables scaled by their bounds). The runs below give the saving. `--no-warm-start` runs every candidate from the deck's initial state.
- **Seeds.** The first population is the base deck's own point and a uniform sample of the box from `--seed`. NSGA-II uses the same seed.
- **The results file.** Every evaluation is appended to `results.csv` in `--out` (default `examples/output/case_study/`): the deck, the study, the variables, the status, the CSS cycle, the cycles run, the run-health codes, whether the run was warm-started, the wall time, the four KPIs, each constraint's value, each margin and `feasible`.
- **Resuming.** Run the same command again. The candidates that are already in the file are read back and not run again. With the same seed, NSGA-II proposes the same candidates, so the search continues where it stopped. A larger `--gen` continues a finished search.

### The three fronts

One run of the six studies on a 4-core machine (2026-11-09, seed 1, warm starts on), with the opening ramps as decision variables. The `purity` studies had a population of 24 and ran 25 generations, the `energy` studies 24 and 20, the `--unconstrained` studies 16 and 15. That is 2640 runs in total, three times the 844 runs of an earlier search without the ramps, in about 4.8 hours. 2573 runs reached cyclic steady state. Each front below takes every converged evaluation of the deck, whatever its study, that satisfies the constraints named (`--summary`). The searches are still shorter than the paper's 4800 evaluations per front, so the fronts are lower bounds on what is possible.

| Front (purity / recovery, %) | Four-step | LPP |
|---|---|---|
| Paper (Table 3, the paper's model) | A 99.3 / 60.0, B 95.2 / 64.1, C1 83.4 / 68.2 | D1 96.1 / 86.8, E 90.4 / 91.1, F 79.9 / 94.6 |
| Unconstrained (CSS and run health only) | 99.6 / 80.2, 98.0 / 88.4, 97.1 / 90.3, 93.0 / 93.4, 84.8 / 96.0 | 98.8 / 88.4, 95.3 / 88.7, 94.7 / 93.9, 85.2 / 95.0 |
| Plant limits (fluidisation, pressure drop, 10 kPa floor) | 80.1 / 21.8, 72.9 / 26.0, 56.6 / 33.8, 45.1 / 39.2 | 89.2 / 13.9, 83.7 / 38.2, 60.2 / 59.8, 53.5 / 85.2 |
| Earlier search without the ramps, plant limits | 83.1 / 8.5, 57.4 / 13.8, 47.2 / 14.5 | 84.1 / 11.5, 80.5 / 47.8, 53.1 / 57.4 |

At the purity of each paper point, the best recovery on each front:

| Paper point | Paper recovery | Unconstrained | Plant limits |
|---|---|---|---|
| A (99.3 %) | 60.0 % | 83.7 % | not reached |
| B (95.2 %) | 64.1 % | 91.3 % | not reached |
| C1 (83.4 %) | 68.2 % | 96.0 % | not reached |
| D1 (96.1 %) | 86.8 % | 88.4 % | not reached |
| E (90.4 %) | 91.1 % | 94.4 % | not reached |
| F (79.9 %) | 94.6 % | 96.4 % | 39.8 % |

The cost of each practical constraint, alone and together, over the same evaluations (`constraint_cost.csv`):

| Constraints applied | Four-step: best recovery at 95 % purity | Four-step: best purity at 90 % recovery | LPP: best recovery at 95 % purity | LPP: best purity at 90 % recovery |
|---|---|---|---|---|
| None | 91.3 % | 97.1 % | 89.3 % | 94.7 % |
| Approach to fluidisation ≤ 1 | 84.2 % | 83.9 % | 89.3 % | 94.7 % |
| Pressure drop ≤ 10 kPa m⁻¹ | 89.7 % | 94.4 % | 89.3 % | 94.7 % |
| Evacuation ≥ 10 kPa | none reaches 95 % | none reaches 90 % | none reaches 95 % | 45.8 % |
| All three | none reaches 95 % | none reaches 90 % | none reaches 95 % | 45.8 % |

- **Unconstrained, the devices pass every paper point.** The three-bed schedule gives the evacuation a full third of the cycle. The four-step deck reaches 97.1 % purity at 90.3 % recovery, which meets both targets. It evacuates to about 1 kPa, and its pumps lift the bed (approach to fluidisation 4.3). The LPP deck comes within 0.3 points of the purity target (94.7 / 93.9) or 1.3 points of the recovery target (95.3 / 88.7), at 2–3 kPa.
- **With the plant limits, no design meets 95 % purity and 90 % recovery.** None of the 606 evaluations that meet every practical constraint reaches both targets. The `energy` studies therefore have no feasible point. Their output is the candidate closest to the targets: four-step 79.5 / 22.5 at 155 kWh t⁻¹, LPP 50.9 / 88.9 at 109 kWh t⁻¹.
- **The ramps open the fluidisation and pressure-drop constraints.** With the opening ramp, both constraints cost at most 7 points of recovery (the four-step deck's fluidisation, at 95 % purity), and nothing on the LPP deck. The plant-limit fronts moved up against the earlier search: four-step from 57.4 / 13.8 to 72.9 / 26.0 and 80.1 / 21.8, LPP to 53.5 / 85.2 from 53.1 / 57.4. The four-step designs on the plant-limit front open the feed-header valve from 2–7 % of its C<sub>v</sub>. The LPP designs mostly open their valve almost fully at once: the LPP enters from the top, where it does not lift the bed.
- **The 10 kPa floor is the constraint that decides.** Every point of the plant-limit fronts ends its evacuation within 0.9 kPa of the floor. At 10 kPa, 13X still holds much of its CO₂, so the working capacity collapses and the N₂ co-adsorbed in the feed step is not cleared. The floor alone removes every design at 95 % purity on both decks.

### What the laboratory column cannot reach

Under the practical constraints, the laboratory column of this case study (0.3 m of 13X, one vacuum pump, three beds on a 3-bed schedule) does not reach 95 % purity and 90 % recovery together. The best constrained designs are about 89 % purity at 14 % recovery and 54 % purity at 85 % recovery. The fluidisation and pressure-drop limits are met by slow-opening valves at little cost. The vacuum floor is not: the purity–recovery targets need an evacuation to 1–3 kPa, which is what the paper and the unconstrained search use.

These directions could change that. They are not results of this study:

- **A lower vacuum floor.** The floor is a constraint value (`variable[p_evac_floor]`), not a device limit. The unconstrained fronts show what a pump that holds 1–3 kPa at the column does. Whether a plant can hold that at scale, with pipe losses and a pump's real curve, is the question.
- **Pilot scale.** A longer column, a larger diameter and a real header change the ratio of the void and pipe volumes to the bed, the heat transfer through the wall and the pressure drop per metre. All three move the trade-off between the floor and the working capacity.
- **More beds or more steps.** A heavy reflux, a pressure equalisation or a second stage concentrates the CO₂ without a deeper vacuum. The paper's LPP cycle is one such step. Two-stage processes are the common answer for 95/90 at a moderate vacuum.

The summary tables are in [`examples/case_study_results/`](https://github.com/samaffleck/skarstrom/tree/main/examples/case_study_results): `fronts.csv` (the paper points and the two fronts, with the decision variables), `constraint_cost.csv`, `paper.csv`, `best.csv` and `counts.csv`. The full results file (2640 rows) is 1.4 MB, so it is not committed. The notebook [`examples/notebooks/case_study_vsa.ipynb`](https://github.com/samaffleck/skarstrom/blob/main/examples/notebooks/case_study_vsa.ipynb) plots the three fronts from `fronts.csv`.

The warm starts reduced the work by about half. A warm-started candidate ran 38 cycles to CSS on average, against 68 for the cold starts (23 s against 50 s of wall time).

To repeat the runs:

<!-- no-run: needs a repository checkout and pymoo; about five hours on 4 cores -->
```bash
cd examples
python perez_2019_case_study_optimise.py --objectives purity --pop 24 --gen 25 --out runs/
python perez_2019_case_study_optimise.py --objectives energy --pop 24 --gen 20 --out runs/
python perez_2019_case_study_optimise.py --objectives purity --unconstrained --pop 16 --gen 15 --out runs/
python perez_2019_case_study_optimise.py --out runs/ --summary case_study_results/
```
