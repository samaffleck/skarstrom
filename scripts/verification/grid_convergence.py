#!/usr/bin/env python3
"""Grid convergence of cycle KPIs at cyclic steady state (observed order and GCI).

Runs each deck to cyclic steady state on a sequence of grids, each twice as
fine as the last (10, 20, 40, 80, 160 cells by default), with each advection
scheme, and reports for every tracked KPI (purity, recovery, mass
productivity, specific energy) the values on every grid, the observed order of
accuracy from the three finest grids (Richardson) and Roache's grid
convergence index (GCI) of the finest grid, with safety factor 1.25:

    p    = ln((f2 - f1) / (f3 - f2)) / ln 2           f1, f2, f3: coarse to fine
    GCI  = 1.25 |(f3 - f2) / f3| / (2^p - 1)          (relative, of f3)
    f_ex = f3 + (f3 - f2) / (2^p - 1)                 (Richardson extrapolation)

When the two differences change sign the convergence is oscillatory and p is
undefined; the table says so and gives the GCI with p = 1 (the scheme's
formal order at the bed ends) instead. Every run tightens the CSS state
tolerance (default 1e-6) so that what is left of the approach to CSS is small
next to the grid differences; the deck's own CSS method and minimum cycle
count are kept.

Needs the skarstrom package (pip install . from the repository root).

    python scripts/verification/grid_convergence.py                    # default decks, all schemes
    python scripts/verification/grid_convergence.py --decks Perez2019_A --schemes upwind \\
        --cells 10 20 40 80 --json grid.json --markdown grid.md
"""

from __future__ import annotations

import argparse
import json
import math
import sys
import time

import skarstrom as sk

DECKS = ["Perez2019_A", "Ward2024", "Illustrative_TSA_13X"]
SCHEMES = ["upwind", "weno", "van_leer"]
CELLS = [10, 20, 40, 80, 160]
KPIS = [("purity", 1.0, "purity"), ("recovery", 1.0, "recovery"),
        ("productivity_mass", 1.0, "productivity [kg/(kg s)]"),
        ("energy_consumption", 1e-3, "specific energy [kJ/kg]")]
SAFETY_FACTOR = 1.25


def run_case(deck: str, scheme: str, cells: int, state_tol: float, max_cycles: int) -> dict:
    sim = sk.examples.get(deck).load()
    for layer in sim.layers:
        layer.cells = cells
    sim.settings.advection_scheme = scheme
    sim.set("css.state_tolerance", state_tol)
    t0 = time.perf_counter()
    result = sim.run(cycles=max_cycles)
    wall = time.perf_counter() - t0
    row = {"deck": deck, "scheme": scheme, "cells": cells, "status": result.status,
           "css_cycle": result.css_cycle, "wall_s": wall, "kpis": {}}
    for k in sim.kpis(-1):
        row["kpis"][k.component] = {name: getattr(k, name) for name, _, _ in KPIS}
    return row


def richardson(f1: float, f2: float, f3: float, ratio: float = 2.0) -> dict:
    """Observed order, GCI of the finest grid and extrapolated value from three grids."""
    e21, e32 = f2 - f1, f3 - f2
    out = {"order": None, "gci": None, "extrapolated": None, "convergence": "monotone"}
    if any(not math.isfinite(v) for v in (f1, f2, f3)):
        out["convergence"] = "undefined"
        return out
    if e32 == 0.0:
        out.update(order=math.inf, gci=0.0, extrapolated=f3, convergence="converged")
        return out
    r = e21 / e32
    if r <= 0.0:
        # Oscillatory: the order is undefined. Report the GCI with p = 1.
        p = 1.0
        out["convergence"] = "oscillatory"
    else:
        p = math.log(r) / math.log(ratio)
        out["order"] = p
        if p <= 0.0:
            out["convergence"] = "diverging"
            return out
    denom = ratio**p - 1.0
    out["gci"] = SAFETY_FACTOR * abs(e32 / f3) / denom if f3 != 0.0 else None
    out["extrapolated"] = f3 + e32 / denom
    return out


def fmt(v, digits=5):
    if v is None:
        return "-"
    if isinstance(v, float) and not math.isfinite(v):
        return "n/a"
    return f"{v:.{digits}g}"


def tables(rows: list, cells: list) -> tuple[list, str]:
    """Per (deck, scheme, component, KPI): values on every grid and the Richardson estimates."""
    summary = []
    lines = []
    by_case = {}
    for row in rows:
        by_case.setdefault((row["deck"], row["scheme"]), {})[row["cells"]] = row
    for (deck, scheme), grid in by_case.items():
        comps = sorted({c for r in grid.values() for c in r["kpis"]})
        for comp in comps:
            for name, scale, label in KPIS:
                values = [grid[n]["kpis"].get(comp, {}).get(name, math.nan) if n in grid else math.nan
                          for n in cells]
                values = [v if v is not None else math.nan for v in values]
                if all(not math.isfinite(v) or v == 0.0 for v in values):
                    continue  # KPI not defined for this deck (e.g. no pump: no energy)
                est = richardson(*values[-3:]) if len(values) >= 3 else {}
                summary.append({"deck": deck, "scheme": scheme, "component": comp, "kpi": name,
                                "cells": cells, "values": values, **est})
                gci = est.get("gci")
                lines.append("| " + " | ".join(
                    [deck, scheme, f"{comp} {label}"] + [fmt(v * scale) for v in values]
                    + [fmt(est.get("order"), 3), est.get("convergence", "-"),
                       "-" if gci is None else f"{100 * gci:.2g} %",
                       fmt(None if est.get("extrapolated") is None else est["extrapolated"] * scale)]) + " |")
    header = ("| Deck | Scheme | KPI | " + " | ".join(f"N = {n}" for n in cells)
              + " | p (3 finest) | Convergence | GCI (finest) | Extrapolated |\n"
              + "|" + "---|" * (len(cells) + 7))
    return summary, header + "\n" + "\n".join(lines)


def main(argv=None) -> int:
    ap = argparse.ArgumentParser(description=__doc__.splitlines()[0])
    ap.add_argument("--decks", nargs="+", default=DECKS)
    ap.add_argument("--schemes", nargs="+", default=SCHEMES, choices=SCHEMES)
    ap.add_argument("--cells", nargs="+", type=int, default=CELLS,
                    help="grids, each twice as fine as the last")
    ap.add_argument("--state-tol", type=float, default=1e-6, help="css.state_tolerance of every run")
    ap.add_argument("--max-cycles", type=int, default=2000)
    ap.add_argument("--json", help="write the runs and the estimates here")
    ap.add_argument("--markdown", help="write the table here")
    args = ap.parse_args(argv)

    for a, b in zip(args.cells, args.cells[1:]):
        if b != 2 * a:
            ap.error("--cells must double from one grid to the next")

    rows = []
    for deck in args.decks:
        for scheme in args.schemes:
            for n in args.cells:
                row = run_case(deck, scheme, n, args.state_tol, args.max_cycles)
                rows.append(row)
                kp = "; ".join(f"{c}: " + ", ".join(f"{k}={fmt(v, 6)}" for k, v in d.items())
                               for c, d in row["kpis"].items())
                print(f"{deck:16s} {scheme:8s} N={n:4d} {row['status']:11s} CSS at {row['css_cycle']:4d} "
                      f"{row['wall_s']:7.1f} s  {kp}", flush=True)
                if row["status"] != "reached_css":
                    print(f"  warning: {deck} {scheme} N={n} did not reach CSS", file=sys.stderr)

    summary, md = tables(rows, args.cells)
    print()
    print(md)
    if args.json:
        with open(args.json, "w") as f:
            json.dump({"settings": vars(args), "runs": rows, "estimates": summary}, f, indent=1)
    if args.markdown:
        with open(args.markdown, "w") as f:
            f.write(md + "\n")
    return 0


if __name__ == "__main__":
    sys.exit(main())
