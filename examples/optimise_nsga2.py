"""Multi-objective optimisation of a VSA cycle with NSGA-II (pymoo).

Searches the feed time and evacuation pressure of the basic four-step CO2/N2
VSA cycle of Perez et al. (2019, configuration A) for the trade-off between
CO2 purity and recovery at cyclic steady state. Each candidate is one full run
to CSS; candidates are evaluated in parallel processes, one simulation each.

Any optimiser that calls a Python function works the same way (SciPy, Optuna,
BoTorch, a grid): load the deck once, copy it per candidate, set the decision
variables, run, read the KPIs. Physics stays in the C++ engine; the search is
plain Python.

    pip install pymoo
    python examples/optimise_nsga2.py                       # 12 x 5 runs on all cores
    python examples/optimise_nsga2.py --pop 24 --gen 20 --workers 8
"""

from __future__ import annotations

import argparse
import os
import sys
import time
from concurrent.futures import ProcessPoolExecutor
from typing import List, Optional, Tuple

import numpy as np

import skarstrom as sk
from skarstrom import ureg as u  # pint units: 3.2 * u.cm, 40.98 * u.kJ / u.mol

# Decision variables: (path, lower bound, upper bound, unit of the bounds).
VARIABLES = [
    ("step[Feed].time", 60.0, 150.0, "s"),
    ("step[Evacuation].bottom.pressure_exponential_end", 2.0, 6.0, "kPa"),
]

_base: Optional[sk.Simulation] = None


def base_case() -> sk.Simulation:
    """The deck every candidate starts from, with a CSS test suited to a search."""
    sim = sk.examples.Perez2019_A.load()
    sim.css.method = "mass_balance_component"   # CSS when every component's balance closes...
    sim.css.mass_balance_tolerance = 0.005      # ...within 0.5 % over a cycle
    sim.css.min_cycles = 10
    sim.cycle.number_of_cycles = 300            # an upper bound; CSS stops the run
    sim.settings.acceleration_method = "aitken"
    return sim


def _init_worker(deck_json: str) -> None:
    global _base
    _base = sk.Simulation.from_json(deck_json)


def evaluate(x: np.ndarray) -> Tuple[float, float, float, str]:
    """Run one candidate to CSS: (CO2 purity, CO2 recovery, wall time [s], status); NaN if it failed."""
    sim = _base.copy()
    for (path, _, _, unit), value in zip(VARIABLES, x):
        sim.set(path, float(value), unit=unit)
    t0 = time.perf_counter()
    result = sim.run()
    if not result.ok:
        return float("nan"), float("nan"), time.perf_counter() - t0, result.status
    kpi = next(k for k in sim.kpis(cycle=-1) if k.component == "CarbonDioxide")
    return kpi.purity, kpi.recovery, time.perf_counter() - t0, result.status


def main(argv: Optional[List[str]] = None) -> int:
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--pop", type=int, default=12, help="population size (default 12)")
    ap.add_argument("--gen", type=int, default=5, help="generations (default 5)")
    ap.add_argument("--workers", type=int, default=os.cpu_count() or 1, help="parallel runs (default: all cores)")
    ap.add_argument("--seed", type=int, default=1)
    ap.add_argument("--cycles", type=int, default=None, help="cap the cycles per run (a quick look, not CSS)")
    ap.add_argument("--no-plot", action="store_true", help="do not save the Pareto plot")
    ap.add_argument("-q", "--quiet", action="store_true", help="only the summary")
    args = ap.parse_args(argv)

    try:
        from pymoo.algorithms.moo.nsga2 import NSGA2
        from pymoo.core.problem import Problem
        from pymoo.optimize import minimize
    except ImportError:
        print("this example needs pymoo: pip install pymoo", file=sys.stderr)
        return 1

    base = base_case()
    if args.cycles is not None:
        base.cycle.number_of_cycles = args.cycles
    deck = base.to_json()
    lower = np.array([v[1] for v in VARIABLES])
    upper = np.array([v[2] for v in VARIABLES])
    runs: List[Tuple[np.ndarray, float, float]] = []

    with ProcessPoolExecutor(max_workers=args.workers, initializer=_init_worker, initargs=(deck,)) as pool:

        class CycleDesign(Problem):
            def __init__(self) -> None:
                super().__init__(n_var=len(VARIABLES), n_obj=2, xl=lower, xu=upper)

            def _evaluate(self, X, out, *_, **__):
                results = list(pool.map(evaluate, X))
                F = []
                for x, (purity, recovery, seconds, status) in zip(X, results):
                    runs.append((x, purity, recovery))
                    if not args.quiet:
                        print(f"  feed {x[0]:6.1f} s  evacuation {x[1]:4.2f} kPa  ->  {status:11s} "
                              f"purity {purity:6.4f}  recovery {recovery:6.4f}  ({seconds:4.1f} s)")
                    # pymoo minimises; a failed run scores worst on both objectives.
                    F.append([1.0, 1.0] if np.isnan(purity) else [-purity, -recovery])
                out["F"] = np.array(F)

        t0 = time.perf_counter()
        res = minimize(CycleDesign(), NSGA2(pop_size=args.pop), ("n_gen", args.gen), seed=args.seed, verbose=False)
        elapsed = time.perf_counter() - t0

    finished = sum(1 for _, p, _ in runs if not np.isnan(p))
    print(f"\n{len(runs)} runs ({finished} completed) in {elapsed:.0f} s on {args.workers} worker(s). Pareto front:")
    print(f"  {'feed time':>10}  {'evacuation':>10}  {'purity':>7}  {'recovery':>8}")
    for x, f in sorted(zip(res.X, res.F), key=lambda r: r[1][0]):
        print(f"  {x[0]:8.1f} s  {x[1]:6.2f} kPa  {-f[0]:7.4f}  {-f[1]:8.4f}")

    if not args.no_plot:
        try:
            import matplotlib

            matplotlib.use("Agg")
            import matplotlib.pyplot as plt
        except ImportError:
            return 0
        done = [(p, r) for _, p, r in runs if not np.isnan(p)]
        fig, ax = plt.subplots(figsize=(5, 4))
        ax.scatter([r for _, r in done], [p for p, _ in done], s=12, color="0.7", label="all candidates")
        ax.scatter(-res.F[:, 1], -res.F[:, 0], s=24, color="C0", label="Pareto front")
        ax.set_xlabel("CO2 recovery [-]")
        ax.set_ylabel("CO2 purity [-]")
        ax.set_title("Perez 2019 A: NSGA-II")
        ax.legend()
        out_dir = os.environ.get("SKARSTROM_EXAMPLES_OUTPUT") or os.path.join(
            os.path.dirname(os.path.abspath(__file__)), "output")
        out = os.path.join(out_dir, "optimise_nsga2.png")
        os.makedirs(os.path.dirname(out), exist_ok=True)
        fig.tight_layout()
        fig.savefig(out, dpi=150)
        print(f"plot saved to {out}")
    return 0


if __name__ == "__main__":
    sys.exit(main())
