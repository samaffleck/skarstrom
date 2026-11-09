# Case Study: a CO₂ VSA Under Plant Constraints (Base Decks)

!!! note "New in 0.9"
    The decks `Perez2019_CaseStudy_Base` and `Perez2019_CaseStudy_Base_LPP` and the builder script `examples/perez_2019_case_study.py`.

This case study starts from the validated design of one CO₂/N₂ VSA paper, [Estupiñan Perez et al. (2019)](perez-2019-vpsa.md), adds the engineering constraints a real plant must meet, and re-optimises under them. This page describes its first step: the **base decks**. They run the paper's cycles at the paper's operating points, but the way a plant runs them:

- **Steps are driven by devices, not imposed profiles.** A blower delivers the feed, the light product leaves through a valve, the bed is pressurised through a valve, and one vacuum pump does the blowdown and the evacuation.
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
| Pressurisation (four-step) | Valve from the feed header | A valve stream from `p_feed_header`, C<sub>v</sub> = `cv_pressurisation` |
| Pressurisation (LPP) | Valve between the beds | A valve stream at each end, C<sub>v</sub> = `cv_lpp`, connected so that each reads the other bed's pressure |

The power falls to zero as the pressure ratio goes to 1. Near the atmosphere the pump then moves a finite volume flow, as a constant-speed pump does. With a constant power, the flow at the start of blowdown is unbounded: the pressure ratio is 1, so the power can move any flow. With a constant power the approach to fluidisation reached 19–31 at the start of blowdown. With the power curve it stays at 1.05 or less in blowdown and evacuation. The four-step deck's pressurisation valve still reaches 3.6 for a moment at the feed end when it opens on the evacuated bed (the validated `Perez2019_A` reaches 1.1 with its linear ramp); a ramped or smaller valve is a constraint for the optimisation step.

## Matching the paper's operating points

The device settings were chosen so that the device-driven steps reach the pressures that the validated decks reach at cyclic steady state. These are not the nominal end values of the validated decks' imposed profiles, because those profiles do not reach their end values in the step time. The pump power and the blowdown share were solved together for the two vacuum pressures. The pressurisation valve was chosen to reach 1 bar within the paper's 30 s. The vent valve keeps the feed step within 0.7 kPa of 1 bar. The tolerance is 2 % on the step-end pressures.

| Setting | Four-step | LPP | Why |
|---|---|---|---|
| `vacuum_pump_power` (W<sub>max</sub>) | 10.4 W | 5.8 W | Blowdown and evacuation end pressures (solved with the share) |
| `blowdown_share` | 0.197 | 0.237 | as above |
| `vp_eta_max`, `vp_alpha` | 0.8, 19.55 | 0.8, 19.55 | The engine's default vacuum-pump curve |
| `cv_vent` | 0.5 mol s⁻¹ bar⁻¹ | 0.5 | Feed-step pressure within 0.7 kPa of 1 bar |
| `cv_pressurisation` | 0.006 mol s⁻¹ bar⁻¹ | — | 1 bar within the paper's 30 s |
| `cv_lpp` | — | 0.008 mol s⁻¹ bar⁻¹ | The largest valve that does not draw gas back in through the Feed 1 vent; 99 % of the end pressure after 16.8 s |
| `blower_efficiency` | 0.7 | 0.7 | Typical. The blower uses less than 1 % of the energy |

At cyclic steady state (this engine version, 80 cells; specific energy = vacuum-pump and blower work per kg of CO₂, with the same pump curve added to the validated decks for comparison):

| Quantity | `Perez2019_A` | Base | `Perez2019_D1` | Base LPP |
|---|---|---|---|---|
| Cycle time | 205 s | 297 s | 360 s | 433.5 s |
| Blowdown: time / end pressure | 19 s / 7.21 kPa | 19.5 s / 7.27 kPa (+0.8 %) | 21 s / 8.16 kPa | 34.2 s / 8.12 kPa (−0.5 %) |
| Evacuation: time / end pressure | 57 s / 4.00 kPa | 79.5 s / 3.99 kPa (−0.1 %) | 167 s / 3.11 kPa | 110.3 s / 3.14 kPa (+0.8 %) |
| Pressurisation: time / end pressure | 30 s / 99.98 kPa | 99 s / 100.00 kPa (99 % after 22.8 s) | 27.5 s / 101.30 kPa | 27.5 s / 100.01 kPa (−1.3 %, 99 % after 16.8 s) |
| CO₂ purity | 99.55 % | 99.59 % | 94.08 % | 94.06 % |
| CO₂ recovery | 60.1 % | 61.4 % | 89.5 % | 88.6 % |
| Productivity, 10⁻⁵ kg kg⁻¹ s⁻¹ | 3.68 | 2.61 | 1.56 | 1.28 |
| Specific energy, MJ kg⁻¹ | 0.752 | 0.735 | 0.821 | 0.807 |
| Cyclic steady state | cycle 48 | cycle 48, 41 s | cycle 92 | cycle 89, 55 s |

Purity, recovery and specific energy agree with the validated decks to within 1.3 points and 3 %. The productivity is lower because the n-bed schedule makes the cycle longer: the feed time is the paper's, but the evacuation and pressurisation get a full third of the cycle each. The LPP pressurisation ends just below the feed bed's pressure, at 100.0 kPa. The validated deck pressurises to 1 atm while its feed vents to 1 bar.

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
| `variable[cv_vent]` | 0.05–2 mol s⁻¹ bar⁻¹ | 0.05–2 |

The base decks have no pressure set points: the pump power and the blowdown share set the blowdown and evacuation pressures, and the valves set the feed and pressurisation pressures. The pump power and the valve coefficients are searched in log10 (`LOG_SCALE`). Each spans one to two decades, and the plant constraints hold only near the small end, which a uniform search of a few generations seldom samples.

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

The operating-limit values are those at the logged times. A peak between two log points is not seen. The backflow is always recorded. It is not enforced by default: the LPP base deck itself draws 14.4 % of its cycle feed back in through the Feed 1 vent at cyclic steady state, so every LPP candidate would be infeasible. That is a known issue of the connected LPP valves ([#241](https://github.com/samaffleck/skarstrom/issues/241)), not of the deck; the backflow becomes a default constraint when it is fixed. The vacuum floor is a constraint on the result, not a limit of the pump model.

The base points are infeasible. The table gives their values (`--base-only`):

| Constraint | `Perez2019_CaseStudy_Base` | `Perez2019_CaseStudy_Base_LPP` |
|---|---|---|
| Purity ≥ 95 % | 99.6 % | 94.1 % (fails) |
| Recovery ≥ 90 % | 61.4 % (fails) | 88.6 % (fails) |
| Approach to fluidisation ≤ 1 | 3.62, pressurisation valve opening (fails) | 0.30 |
| Pressure drop ≤ 10 kPa m⁻¹ | 14.1 kPa m⁻¹, pressurisation (fails) | 16.2 kPa m⁻¹, LPP from the top (fails) |
| Lowest vacuum pressure ≥ 10 kPa | 3.9 kPa (fails) | 3.1 kPa (fails) |
| Outlet backflow (recorded) | 0 | 14.4 % of the feed |

In the feed step the pressure drop is about 0.6 kPa m⁻¹. The pressure-drop bound acts on the pressurisation transients, where a valve opens onto an evacuated bed. The optimiser must use slower valves, a weaker pump or shorter evacuation, and keep purity and recovery.

### Warm starts, parallel runs and resuming

- **Warm starts.** Each worker process keeps a `skarstrom.optimisation.StateCache`. Each candidate starts from the cyclic steady state of the nearest candidate that the worker has run (`skarstrom.optimisation.run_warm`, distance in the variables scaled by their bounds). The runs below give the saving. `--no-warm-start` runs every candidate from the deck's initial state.
- **Seeds.** The first population is the base deck's own point and a uniform sample of the box from `--seed`. NSGA-II uses the same seed.
- **The results file.** Every evaluation is appended to `results.csv` in `--out` (default `examples/output/case_study/`): the deck, the study, the variables, the status, the CSS cycle, the cycles run, the run-health codes, whether the run was warm-started, the wall time, the four KPIs, each constraint's value, each margin and `feasible`.
- **Resuming.** Run the same command again. The candidates that are already in the file are read back and not run again. With the same seed, NSGA-II proposes the same candidates, so the search continues where it stopped. A larger `--gen` continues a finished search.

### The three fronts

A short search on one 4-core machine (2026-11-09, seed 1, warm starts on). The constrained studies had a population of 16: the `energy` studies ran 6 generations and the `purity` studies 12 (resumed from 6). The `--unconstrained` studies had a population of 12 and ran 4 generations. That is 844 runs in total, at about 13 s each. Each front below takes every converged evaluation of the deck, whatever its study, that satisfies the constraints named. These searches are far shorter than the paper's 4800 evaluations per front, so the fronts are lower bounds on what is possible.

| Front (purity / recovery) | Four-step | LPP |
|---|---|---|
| Paper (Table 3, the paper's model) | A 99.3 / 60.0, B 95.2 / 64.1, C1 83.4 / 68.2 | D1 96.1 / 86.8, E 90.4 / 91.1, F 79.9 / 94.6 |
| Unconstrained (CSS and run health only) | 100 / 78.5, 92.7 / 79.1, 87.6 / 85.1 | 99.2 / 72.9, 97.4 / 83.7, 94.9 / 87.8, 89.1 / 93.6, 86.9 / 95.1 |
| Plant limits (fluidisation, pressure drop, 10 kPa floor) | 83.1 / 8.5, 57.4 / 13.8, 47.2 / 14.5 | 84.1 / 11.5, 80.5 / 47.8, 53.1 / 57.4 |

- **Unconstrained, the devices reach and pass the paper's fronts.** At the paper's purity, the four-step deck recovers 78.5 % against A's 60.0 % and B's 64.1 %, and 85.1 % against C1's 68.2 %. The three-bed schedule gives the evacuation a full third of the cycle. The LPP deck passes F (95.1 % against 94.6 %), and it is 2.5 points below E and 3.1 points below D1 after 4 generations. Its unconstrained front evacuates to 2–4 kPa, the depth at which the paper's own search ended (its P<sub>low</sub> converged to its 3 kPa bound).
- **No evaluation meets 95 % purity and 90 % recovery.** This holds with and without the plant limits. The closest is the LPP deck's 94.9 / 87.8 at 229 kWh t⁻¹, which is unconstrained (and the LPP base point, 94.1 / 88.6 at 224 kWh t⁻¹). The four-step cycle stays below 86 % recovery at any purity. The `energy` studies therefore have no feasible point. Their output is the least-violating candidate.
- **The 10 kPa floor is the constraint that decides.** Every point of the plant-limit fronts ends its evacuation within 1 kPa of the floor (10.0–11.0 kPa). At that pressure 13X keeps much of its CO₂, so the working capacity collapses. The best constrained LPP point is 80.5 % purity at 47.8 % recovery, for 175 kWh t⁻¹. Fluidisation and pressure drop are met by slower valves (the four-step deck's pressurisation C<sub>v</sub> falls to about 0.0005), and they cost little.
- **The four-step deck's pressurisation valve.** The valve opens fully onto the evacuated bed, so the approach to fluidisation at the opening falls below 1 only for C<sub>v</sub> below about 0.001. That valve then cannot pressurise the bed in its step (95 kPa at C<sub>v</sub> = 0.001), and the feed must finish the pressurisation. A valve that opens on a ramp, a [log-linear pressure ramp](boundary-conditions.md#choosing-a-pressurisation-ramp) for example, is the device change that would remove this trade-off.

The warm starts reduced the work by about 40 %. A warm-started candidate ran 24 cycles to CSS on average, against 40 for the cold starts (12.7 s against 19.4 s of wall time). The results file of these runs is 416 kB, so this page gives the summary.

To repeat the runs:

<!-- no-run: needs a repository checkout and pymoo; about two hours on 4 cores -->
```bash
cd examples
python perez_2019_case_study_optimise.py --objectives energy --pop 16 --gen 6 --out runs/
python perez_2019_case_study_optimise.py --objectives purity --pop 16 --gen 12 --out runs/
python perez_2019_case_study_optimise.py --objectives purity --unconstrained --pop 12 --gen 4 --out runs/
```
