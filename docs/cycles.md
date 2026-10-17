# Cycles and Initial Conditions

A **cycle** is the sequence of operating steps a bed repeats — feed, blowdown, purge, pressurisation, and so on. In Skarstrom a cycle is a list of **steps**, each with its own duration and its own pair of [boundary conditions](boundary-conditions.md), executed in order and then repeated. A deck's cycle also holds the [connections between steps](connecting-beds.md) for multi-bed processes and the KPI configuration for the [performance metrics](performance-metrics.md).

> **Start from a template.** `data/cycles/` holds ready-made cycles: `breakthrough.json`, `4_step_vpsa.json`, `4_step_vpsa_lpp.json` (with light-product pressurisation), `skarstrom_cycle.json` and `6_step_vpsa_rinse.json`. They are written by `adsim/test/unit/test_cycle_designs.cpp`, which is also the clearest example of building a cycle in C++, and load with `adsim::LoadCycleFromJSON`. The [example decks](examples.md) are complete processes to adapt.

---

## Steps

Each step has a **name** and two boundary conditions:

- the **Top** boundary at $z = L$ (the far end of the column), and
- the **Bottom** boundary at $z = 0$ (the near end).

A step can also apply an [energy source](energy-balance.md) to a layer (heating or cooling). In C++ a cycle is edited with `Cycle::AddStep`, `Cycle::GetStep`, `DuplicateStep` (which copies a step with its settings and clears connections on the copy) and the boundary-condition factories; from Python the existing steps' values are parameters:

```python
sim.steps                          # ['Feed', 'Blowdown', 'Evacuation', 'Pressurisation']
sim["step[Feed].time"] = 60.0      # by name (alias) ...
sim["step[1].time"] = 60.0         # ... or by the step's stable ID
```

A few rules are enforced by the deck checks before a run:

- A cycle must contain at least one step, and step names must be unique.
- The two boundaries of a step cannot both be **inlets** — one must be an outlet (or a closed end).
- The Top and Bottom boundaries must be at different ends of the column.

## Number of cycles

A PSA/TSA process must run many cycles before it settles into a repeating pattern. The deck's **number of cycles** (`cycle.number_of_cycles`) is how many times the whole step sequence repeats; the total simulated time is the cycle time multiplied by the number of cycles. `sim.run(cycles=N)` and `skarstrom run --cycles N` override it for one run.

Rather than guessing the number of cycles, let the solver stop automatically once the bed stops changing from one cycle to the next — see [cyclic steady state](cyclic-steady-state.md). The number of cycles is then an upper limit.

## Step timing

The **step time mode** controls how step durations are given:

- **Absolute** — each step's duration in seconds (`step[...].time`). The cycle time is the sum of the step times.
- **Fraction** — a single **total cycle time** (`cycle.total_time`) and a fraction of it per step (`step[...].fraction`). The fractions must sum to 1. This makes it easy to rescale a whole cycle by changing one number.

Step times are expressions, so two steps can share a [global variable](global-variables.md) (`sim["step[Rinse].time"] = "t_feed"`).

---

## Initial conditions

The initial conditions set the state of the bed everywhere at $t = 0$, before the first step begins:

- **Pressure**
- **Bed temperature** (the initial gas and solid temperature)
- **Wall temperature**
- **Velocity**
- **Initial mole fraction** for each component — these must sum to 1.

Choose initial conditions that are physically consistent with the first step's boundary conditions; for example, if the first step feeds gas at high pressure into an evacuated bed, start from a low pressure. A poor initial guess can make the first solve harder to converge. (If the model reaches [cyclic steady state](cyclic-steady-state.md), the initial condition only affects how many cycles that takes, not the final answer.)

---

Next: define what happens at each end of the column in [Boundary Conditions](boundary-conditions.md).
