#!/usr/bin/env python3
"""Time to cyclic steady state for the bundled examples (the README's Performance table).

Runs each deck as bundled, three times, and prints the median wall time. The
decks declare cyclic steady state when the mass balance of every component
over a cycle closes within 0.5 % (css.method "mass_balance_component", at
least 5 cycles), without acceleration. Ward2024 is also run on a finer grid
(90 cells instead of 30) to show how the cost scales with the number of cells.
Needs the skarstrom package (pip install . from the repository root).

    python scripts/benchmark_css.py
"""

from __future__ import annotations

import statistics
import time

import skarstrom as sk

# (deck, cells per layer; None keeps the deck's grid)
CYCLIC = [("Perez2019_A", None), ("Nguyen2023_Exp_VI", None), ("Ferre2026_NaYF", None),
          ("Ward2024", None), ("Ward2024", 90)]
REPEATS = 3


def timed(sim: sk.Simulation):
    t0 = time.perf_counter()
    result = sim.run()
    return time.perf_counter() - t0, result


def main() -> None:
    print(f"skarstrom {sk.__version__}; median of {REPEATS} runs")
    for name, cells in CYCLIC:
        times = []
        for _ in range(REPEATS):
            sim = sk.examples.get(name).load()
            if cells is not None:
                for layer in sim.layers:
                    layer.cells = cells
            dt, result = timed(sim)
            times.append(dt)
        n_cells = sum(layer.cells for layer in sim.layers)
        print(f"{name:20s} {n_cells:4d} cells  {result.status:12s} "
              f"CSS at cycle {result.css_cycle + 1:4d}  {statistics.median(times):6.2f} s")


if __name__ == "__main__":
    main()
