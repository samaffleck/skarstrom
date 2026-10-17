#!/usr/bin/env python3
"""Time to cyclic steady state for the bundled examples (the README's Performance table).

Runs each case three times and prints the median wall time. Cyclic steady
state is declared when the mass balance of every component over a cycle
closes within 0.5 % (css.method "mass_balance_component", at least 5 cycles),
with and without Aitken (Graves-Morris) acceleration. The per-component
balance is the strict test: on the light-product pressurisation cycles the
gas profiles and the total balance settle hundreds of cycles before the CO2
inventory does. Needs the skarstrom package
(pip install . from the repository root).

    python scripts/benchmark_css.py
"""

from __future__ import annotations

import statistics
import time

import skarstrom as sk

CYCLIC = ["Perez2019_E", "Perez2019_A", "Nguyen2023_Exp_VI"]
BREAKTHROUGH = ["Ward2022_Q100ccm", "Streb2021_A1"]
REPEATS = 3


def timed(sim: sk.Simulation):
    t0 = time.perf_counter()
    result = sim.run()
    return time.perf_counter() - t0, result


def main() -> None:
    print(f"skarstrom {sk.__version__}; median of {REPEATS} runs")
    for name in CYCLIC:
        for acceleration in ("none", "aitken"):
            times = []
            for _ in range(REPEATS):
                sim = sk.examples.get(name).load()
                sim.css.method = "mass_balance_component"
                sim.css.min_cycles = 5
                sim.css.mass_balance_tolerance = 0.005
                sim.settings.acceleration_method = acceleration
                dt, result = timed(sim)
                times.append(dt)
            cells = sum(layer.cells for layer in sim.layers)
            kpi = sim.kpis(cycle=-1)[0]
            print(f"{name:20s} {cells:4d} cells  acceleration={acceleration:6s}  {result.status:12s} "
                  f"CSS at cycle {result.css_cycle + 1:4d}  {statistics.median(times):6.2f} s  "
                  f"{kpi.component} purity {kpi.purity:.4f} recovery {kpi.recovery:.4f}")
    for name in BREAKTHROUGH:
        times = []
        for _ in range(REPEATS):
            dt, result = timed(sk.examples.get(name).load())
            times.append(dt)
        print(f"{name:20s} breakthrough  {result.status:12s} {statistics.median(times):6.2f} s")


if __name__ == "__main__":
    main()
