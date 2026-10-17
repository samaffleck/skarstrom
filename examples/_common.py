"""Shared helpers for the example scripts in this folder.

Not part of the ``skarstrom`` package: the scripts import it as a sibling
module. It finds the decks, the digitised experimental data (adsim/data/) and
the published reference values (adsim/test/papers/*/reference_values.csv) of a
repository checkout, and falls back to the decks bundled with the installed
package when the scripts are run outside the repository (without the
comparison data).
"""

from __future__ import annotations

import argparse
import csv
import math
import os
import time
from pathlib import Path
from typing import Dict, Iterable, List, Optional, Sequence, Tuple

import numpy as np

import skarstrom

HERE = Path(__file__).resolve().parent
REPO = HERE.parent
DECKS = REPO / "examples" / "decks"
DATA = REPO / "adsim" / "data"
PAPERS = REPO / "adsim" / "test" / "papers"
OUTPUT = Path(os.environ.get("SKARSTROM_EXAMPLES_OUTPUT", HERE / "output"))


# ---- decks -----------------------------------------------------------------------

def deck_path(name: str) -> str:
    """Path of an example deck: the repository copy if present, else the bundled one."""
    p = DECKS / name
    if p.is_file():
        return str(p)
    return skarstrom.example_path(name)


def load(name: str, cycles: Optional[int] = None) -> skarstrom.Simulation:
    sim = skarstrom.Simulation(deck_path(name))
    if cycles is not None:
        sim["cycle.number_of_cycles"] = int(cycles)
    return sim


def run(sim: skarstrom.Simulation, label: str, quiet: bool = False) -> skarstrom.RunResult:
    """Run a deck and print a one-line status (raises if the run fails)."""
    t0 = time.perf_counter()
    result = sim.run()
    dt = time.perf_counter() - t0
    css = f", CSS at cycle {result.css_cycle + 1}" if result.css_reached else ""
    print(f"[{label}] {result.status}: {result.cycles_run} cycle(s) in {dt:.1f} s{css}")
    if not result.ok:
        raise RuntimeError(f"{label}: run failed: {result.message}")
    return result


def css_status(sim: skarstrom.Simulation) -> str:
    """The deck's cyclic-steady-state check at the last logged cycle, as text."""
    hist = sim.css_history()
    if not hist or hist[-1].method == "None":
        return "CSS not checked (the deck's CSS method is None)"
    c = hist[-1]
    state = "converged" if c.converged else "not converged"
    return f"CSS {state} at cycle {c.cycle + 1} ({c.method}: error {c.error:.3g}, tolerance {c.tolerance:.3g})"


def kpi(sim: skarstrom.Simulation, component: str, cycle: int = -1):
    for k in sim.kpis(cycle=cycle):
        if k.component == component:
            return k
    raise KeyError(f"no KPIs for {component} (is it tracked in the deck's KPI configuration?)")


def print_kpis(sim: skarstrom.Simulation, indent: str = "  ") -> None:
    rows = sim.kpis(cycle=-1)
    if not rows:
        print(f"{indent}no KPIs (no component is tracked in the deck's KPI configuration)")
        return
    for k in rows:
        if not k.valid:
            print(f"{indent}{k.component}: no feed/product endpoints configured")
            continue
        print(f"{indent}{k.component}: purity {100 * k.purity:.2f} %, recovery {100 * k.recovery:.2f} %, "
              f"productivity {k.productivity_volume:.4g} kg/(m3 s), energy {k.energy_consumption:.4g} J/kg")


# ---- reference data ----------------------------------------------------------------

def have_reference_data() -> bool:
    return DATA.is_dir() and PAPERS.is_dir()


def load_csv(path: Path) -> Tuple[np.ndarray, np.ndarray]:
    """Two-column digitised curve (header row, then x,y)."""
    xs, ys = [], []
    with open(path, encoding="utf-8-sig", newline="") as f:
        r = csv.reader(f)
        next(r, None)
        for row in r:
            if len(row) >= 2 and row[0].strip():
                xs.append(float(row[0]))
                ys.append(float(row[1]))
    return np.asarray(xs), np.asarray(ys)


def reference(paper: str, case: str, quantity: str) -> Optional[float]:
    """A published value from adsim/test/papers/<paper>/reference_values.csv."""
    p = PAPERS / paper / "reference_values.csv"
    if not p.is_file():
        return None
    with open(p, encoding="utf-8", newline="") as f:
        for row in csv.DictReader(f):
            if row["case"] == case and row["quantity"] == quantity:
                return float(row["value"])
    return None


def nrmse(t_sim: np.ndarray, y_sim: np.ndarray, t_exp: np.ndarray, y_exp: np.ndarray,
          scale: float = 1.0) -> float:
    """RMSE of the simulation (interpolated at the data times) against the data, divided by scale."""
    y = np.interp(t_exp, t_sim, y_sim)
    return float(np.sqrt(np.mean((y - y_exp) ** 2)) / scale)


def crossing_time(t: np.ndarray, y: np.ndarray, level: float) -> float:
    """First time y reaches `level` (linear interpolation), or nan."""
    idx = np.nonzero(y >= level)[0]
    if idx.size == 0:
        return math.nan
    i = int(idx[0])
    if i == 0:
        return float(t[0])
    return float(t[i - 1] + (level - y[i - 1]) * (t[i] - t[i - 1]) / (y[i] - y[i - 1]))


def table(rows: Iterable[Sequence[str]], header: Sequence[str]) -> None:
    rows = [list(map(str, r)) for r in rows]
    widths = [len(h) for h in header]
    for r in rows:
        widths = [max(w, len(c)) for w, c in zip(widths, r)]
    print("  " + "  ".join(h.ljust(w) for h, w in zip(header, widths)))
    print("  " + "  ".join("-" * w for w in widths))
    for r in rows:
        print("  " + "  ".join(c.ljust(w) for c, w in zip(r, widths)))


def fmt(x: Optional[float], digits: int = 4) -> str:
    if x is None or (isinstance(x, float) and math.isnan(x)):
        return "-"
    return f"{x:.{digits}g}"


def endpoint_moles(sim: skarstrom.Simulation, step: int, location: str,
                   component: Optional[str] = None, cycle: int = -1) -> float:
    """Moles through every endpoint of one step and bed end over a cycle (absolute value).

    For decks without a KPI configuration: integrate the logged stream flow
    (``molar_flow`` of a component, or the total ``molar_flow_rate``).
    """
    total = 0.0
    for e in sim.streams():
        if e["step"] != step or e["location"] != location:
            continue
        s = (sim.stream_series(step, e["name"], "molar_flow", component, cycle=cycle) if component
             else sim.stream_series(step, e["name"], "molar_flow_rate", cycle=cycle))
        if len(s) > 1:
            total += float(np.sum(0.5 * (s.values[1:] + s.values[:-1]) * np.diff(s.time)))
    return abs(total)


def evacuation_product_kpis(sim: skarstrom.Simulation, component: str = "CarbonDioxide") -> Dict[str, float]:
    """Purity and recovery of a VSA heavy product drawn from the bottom in "Evacuation".

    The definition the paper tests use (Perez 2019, Nguyen 2023): the product
    is the evacuation outflow; the feed is every bottom inlet of the steps
    named "Feed..." plus a bottom-fed "Pressurisation". Unlike a KPI
    configuration that also counts a light-product-pressurisation inlet as
    (negative) product, this purity is always defined.
    """
    steps = sim.steps
    evac = steps.index("Evacuation")
    feeds = [i for i, n in enumerate(steps) if n.startswith("Feed")]
    if "Pressurisation" in steps:
        feeds.append(steps.index("Pressurisation"))
    fed = sum(endpoint_moles(sim, i, "bottom", component) for i in feeds)
    product = endpoint_moles(sim, evac, "bottom", component)
    return {f"{component}.purity": product / endpoint_moles(sim, evac, "bottom"),
            f"{component}.recovery": product / fed}


def compare_kpis(sim: skarstrom.Simulation, paper: str, case: str,
                 quantities: Sequence[Tuple[str, str, str]], values=None) -> List[Tuple[str, float, Optional[float]]]:
    """Print simulated purity/recovery [%] next to the published values.

    quantities: (label, "<Component>.<purity|recovery>", reference quantity name).
    values: optional callable sim -> {"<Component>.<field>": fraction} for decks
    without a KPI configuration; by default the deck's KPIs are used.
    Returns (label, simulated %, reference %) tuples.
    """
    out = []
    rows = []
    computed = values(sim) if values else None
    for label, what, ref_q in quantities:
        comp, field = what.split(".")
        value = 100.0 * (computed[what] if computed else getattr(kpi(sim, comp), field))
        ref = reference(paper, case, ref_q)
        out.append((label, value, ref))
        diff = "-" if ref is None else f"{value - ref:+.2f} pp"
        rows.append([label, f"{value:.2f} %", "-" if ref is None else f"{ref:.2f} %", diff])
    table(rows, ["quantity", "simulated", "published", "difference"])
    return out


# ---- command line and plotting -------------------------------------------------------

def parser(description: str, cases: Sequence[str], default: Optional[Sequence[str]] = None) -> argparse.ArgumentParser:
    ap = argparse.ArgumentParser(description=description)
    ap.add_argument("--case", action="append", choices=list(cases),
                    help=f"case(s) to run (repeatable; default: {', '.join(default or cases)})")
    ap.add_argument("--cycles", type=int, help="override the decks' number of cycles (quick look)")
    ap.add_argument("--no-plot", action="store_true", help="do not save a plot")
    ap.add_argument("-q", "--quiet", action="store_true")
    ap.set_defaults(_default_cases=list(default or cases))
    return ap


def selected(args: argparse.Namespace) -> List[str]:
    return args.case or args._default_cases


def pyplot(args: argparse.Namespace):
    """matplotlib.pyplot with a non-interactive backend, or None (no plot wanted or not installed)."""
    if getattr(args, "no_plot", False):
        return None
    try:
        import matplotlib
        matplotlib.use("Agg")
        import matplotlib.pyplot as plt
    except ImportError:
        print("(matplotlib is not installed: skipping the plot; pip install matplotlib)")
        return None
    return plt


def save(fig, name: str) -> Path:
    OUTPUT.mkdir(parents=True, exist_ok=True)
    p = OUTPUT / name
    fig.tight_layout()
    fig.savefig(p, dpi=120)
    print(f"plot saved to {p}")
    return p


def cyclic_family(argv, doc: str, paper: str, cases: Dict[str, Tuple[str, str]],
                  quantities: Sequence[Tuple[str, str, str]], plot_component: str, plot_name: str,
                  default: Optional[Sequence[str]] = None, values=None) -> int:
    """Run cyclic decks to CSS and compare their final KPIs with published values.

    cases: name -> (deck file, case name in reference_values.csv).
    quantities, values: see compare_kpis. The plot shows purity and recovery of
    plot_component from the deck's KPI configuration against the cycle number
    (the last logged cycles) for every case run; it is skipped for decks without
    a KPI configuration.
    """
    args = parser(doc.splitlines()[0], list(cases), default).parse_args(argv)
    plt = pyplot(args)
    fig, axes = (plt.subplots(1, 2, figsize=(10, 4)) if plt else (None, None))
    plotted = False
    for name in selected(args):
        deck, ref_case = cases[name]
        sim = load(deck, args.cycles)
        run(sim, f"{paper} {name}", args.quiet)
        print(f"  {css_status(sim)}")
        compare_kpis(sim, paper, ref_case, quantities, values)
        if axes is not None:
            rows = [k for k in sim.kpis() if k.component == plot_component and k.valid]
            cyc = [k.cycle + 1 for k in rows]
            if rows:
                plotted = True
                axes[0].plot(cyc, [100 * k.purity for k in rows], label=name)
                axes[1].plot(cyc, [100 * k.recovery for k in rows], label=name)
        print()
    if not have_reference_data():
        print("(run from a repository checkout to see the published values)")
    if fig is not None and plotted:
        axes[0].set(xlabel="cycle (last logged cycles)", ylabel=f"{plot_component} purity [%], deck KPI configuration")
        axes[1].set(xlabel="cycle", ylabel=f"{plot_component} recovery [%]")
        for a in axes:
            a.legend(fontsize=7)
        save(fig, plot_name)
    return 0
