"""Calibrate a cyclic deck against per-cycle KPIs (purity and recovery against the cycle number).

A plant or pilot rig often reports one purity and one recovery per cycle rather
than curves in time. A dataset bound to a KPI observable
(``sk.Observable("recovery", "CarbonDioxide")``) holds such a history: x is the
cycle number and y the KPI of the deck's KPI configuration in that cycle, and
``skarstrom.calibrate`` matches them cycle by cycle.

This example is a twin experiment on ``Perez2019_A`` (CO2/N2 VSA on 13X), on a
coarse grid and a few cycles from the deck's initial state:

1. run the deck and record the CO2 purity and recovery of each cycle, with
   seeded noise of the stated uncertainty, as the "measured" histories;
2. move the parameters to be estimated away from the deck's values;
3. fit them back to the histories, and compare the estimates with the values
   that made the data.

Runs start from the deck's initial state, not from a cached cyclic steady state
(``warm_start=False``): the histories begin at cycle 1, so the transient is part
of what is measured. The CSS check is off (``css.method``) so every run makes all
the cycles.

    python examples/calibrate_cycle_kpis.py            # 8 cycles, 20 cells
    python examples/calibrate_cycle_kpis.py --quick    # 4 cycles, a few iterations
"""

from __future__ import annotations

import argparse

import numpy as np

import _common as c
import skarstrom as sk
from skarstrom import calibration as cal

DECK = "Perez2019_A"
COMPONENT = "CarbonDioxide"
CELLS = 20
SIGMA = {"purity": 0.002, "recovery": 0.005}   # stated uncertainty of each KPI [-] (fraction)

AFFINITY = "adsorbent[1].isotherm[CarbonDioxide].b1"   # CO2 site-1 affinity [m3/mol]
VOIDAGE = "layer[1].bed_void_fraction"                  # bed void fraction [-]
START = {AFFINITY: 1.5, VOIDAGE: 0.38}                  # the start: b1 x 1.5, voidage 0.38 (deck: 0.41)


def measured(cycles: int, seed: int = 1) -> dict:
    """The "measurements": the deck's own KPI history, with seeded noise of SIGMA."""
    sim = c.load(DECK, cycles)
    sim["layer[1].cells"] = CELLS
    sim["css.method"] = 0          # no CSS check: run every cycle
    result = sim.run()
    if not result.ok:
        raise RuntimeError(f"{DECK}: {result.message}")
    rows = [k for k in sim.kpis() if k.component == COMPONENT]
    rng = np.random.default_rng(seed)
    data = {"cycle": np.array([k.cycle + 1 for k in rows], dtype=float)}   # 1-based cycle numbers
    for kpi, sigma in SIGMA.items():
        data[kpi] = np.array([getattr(k, kpi) for k in rows]) + rng.normal(0.0, sigma, len(rows))
    data["truth"] = {AFFINITY: sim.get(AFFINITY), VOIDAGE: sim.get(VOIDAGE)}
    return data


def case(data: dict, cycles: int) -> cal.Case:
    """The deck with the two KPI histories bound as datasets."""
    sim = c.load(DECK)
    sim["layer[1].cells"] = CELLS
    for kpi, sigma in SIGMA.items():
        sim.experimental_data.add(
            f"{COMPONENT}_{kpi}", data["cycle"], data[kpi], x_label="cycle", y_label=f"{kpi} [-]",
            observable=sk.Observable(kpi, COMPONENT, cycle_offset=-1),   # x = 1-based cycle number
            x_unit="cycle", uncertainty=[sigma] * len(data["cycle"]),
            weight=float(len(data["cycle"])))   # weight = points: every point counts once, sigma^2 ~ 1
    return cal.Case(sim, [f"{COMPONENT}_{kpi}" for kpi in SIGMA], cycles=cycles,
                    settings={"css.method": 0}, name=DECK)


def main(argv=None) -> int:
    ap = argparse.ArgumentParser(description=__doc__.splitlines()[0])
    budget = ap.add_mutually_exclusive_group()
    budget.add_argument("--quick", action="store_true", help="4 cycles and at most 20 runs")
    budget.add_argument("--full", action="store_true", help="8 cycles, central differences (default)")
    ap.add_argument("--threads", type=int, help="simulations at once (default: the cores - 1)")
    ap.add_argument("--no-plot", action="store_true", help="do not save a plot")
    ap.add_argument("-q", "--quiet", action="store_true")
    args = ap.parse_args(argv)

    cycles = 4 if args.quick else 8
    options = (dict(differences="forward", max_simulations=20) if args.quick
               else dict(differences="central", max_simulations=120))

    print(f"Making the measured KPI histories: {DECK}, {cycles} cycles, {CELLS} cells")
    data = measured(cycles)
    truth = data["truth"]
    sim_case = case(data, cycles)
    params = [
        cal.Parameter(AFFINITY, truth[AFFINITY] / 4, truth[AFFINITY] * 4, scale="log",
                      start=truth[AFFINITY] * START[AFFINITY]),
        cal.Parameter(VOIDAGE, 0.30, 0.50, start=START[VOIDAGE]),
    ]

    def progress(info):
        if not args.quiet:
            print(f"  iteration {info['iteration']}: S = {info['objective']:.6g} "
                  f"({info['simulations']} simulations)", flush=True)

    result = sk.calibrate(sim_case, params, warm_start=False, threads=args.threads, progress=progress, **options)
    print()
    print(result.table())
    print()
    rows = []
    for p in result.parameters:
        t = truth[p.path]
        rows.append([p.path, c.fmt(t), c.fmt(p.start), c.fmt(p.value), c.fmt(p.std_error, 2),
                     f"{(p.value - t) / p.std_error:+.2f}" if p.std_error > 0 else "n/a"])
    print("The estimates against the values that made the data:")
    c.table(rows, ["parameter", "true", "start", "estimate", "std error", "(estimate - true) / std error"])
    print()
    if not result.ok:
        print(f"the fit failed: {result.message}")
        return 1
    check = result.calibrated()
    check.cycle.number_of_cycles = cycles
    check["css.method"] = 0
    print(f"calibrated deck: {check.run().status}")

    if not args.no_plot:
        plt = c.pyplot(args)
        if plt is not None:
            fig = result.plot()
            c.save(fig, "calibrate_cycle_kpis.png")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
