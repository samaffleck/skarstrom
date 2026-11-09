"""Perez 2019 case study: re-optimise the device-driven VSA under plant constraints.

Searches the device settings and the schedule of the case-study base decks
(``Perez2019_CaseStudy_Base`` and ``_Base_LPP``, examples/perez_2019_case_study.py)
with NSGA-II (pymoo). Each candidate is one run to cyclic steady state (CSS);
candidates run in parallel processes, each warm-started from the CSS of the
nearest candidate that process has already run.

Two objective sets (``--objectives``):

* ``energy`` (default): minimise the specific electrical energy [kWh per tonne
  of CO2] and maximise the productivity, subject to every constraint below,
  including the purity and recovery targets. This is the design question.
* ``purity``: maximise CO2 purity and recovery, subject to the practical
  constraints only (not the targets). This front is the one the paper
  optimised, so it can be compared with the paper's operating points. Add
  ``--unconstrained`` to drop the practical constraints too (only CSS and a
  healthy run are required): the reproduction of the paper's search with the
  devices of the case study.

The two questions need different objectives, so they are two runs: with all
four objectives in one run the front is a 4-D surface that a small population
cannot cover, and the purity and recovery targets would be both objectives
and constraints.

Constraints (a candidate is feasible when every margin is >= 0; NSGA-II ranks
infeasible candidates by their summed normalised violation):

* CO2 purity >= ``variable[purity_target]`` and recovery >= ``variable[recovery_target]``
  (the KPI record of the last cycle; ``energy`` only);
* approach to fluidisation <= 1 on every packed face in every step
  (``sim.operating_limits()``, at the logged times);
* bed pressure drop per metre <= ``--max-dp-per-m`` (default 10 kPa/m): the
  largest |P_top - P_bottom| / L over the logged times of every step;
* lowest bed pressure in blowdown and evacuation >= ``variable[p_evac_floor]``
  (the operating-limit report's minimum over the column);
* with ``--backflow-constraint``: no gas drawn back in through an unconnected
  outlet (the operating-limit report's outlet-backflow warning, above 0.1 % of
  the cycle feed: the vents have no check valve, and gas that flows back
  carries the bed-end composition, not the outside gas). It is always
  recorded; it is not enforced by default because the LPP base deck itself
  draws gas back in through its feed-1 vent (see docs/case-study-vsa.md);
* the run reached CSS and its run-health checks passed (``result.healthy``).

Every evaluation is appended to ``results.csv`` (``--out``): the variables, the
objectives, each constraint's value and margin, the CSS cycle, the cycles run,
whether the run was warm-started, and the wall time. Running the same command
again resumes: the candidates already in the file are read back instead of
run, and NSGA-II (same ``--seed``) proposes the same candidates again, so the
search continues where it stopped.

    pip install pymoo
    python examples/perez_2019_case_study_optimise.py                         # both decks, energy
    python examples/perez_2019_case_study_optimise.py --deck Base --objectives purity
    python examples/perez_2019_case_study_optimise.py --deck Base --objectives purity --unconstrained
    python examples/perez_2019_case_study_optimise.py --pop 24 --gen 20 --workers 8 --seed 2
    python examples/perez_2019_case_study_optimise.py --deck Base --base-only  # the base deck's own point

Reference: L. Estupinan Perez, P. Sarkar, A. Rajendran, Separation and
Purification Technology 224 (2019) 553-563. doi:10.1016/j.seppur.2019.05.039
"""

from __future__ import annotations

import argparse
import csv
import math
import os
import sys
import time
from concurrent.futures import ProcessPoolExecutor
from pathlib import Path
from typing import Dict, List, Optional, Sequence, Tuple

import numpy as np

import _common as c
import skarstrom as sk
from skarstrom import optimisation as opt

DECKS = {"Base": "Perez2019_CaseStudy_Base", "Base_LPP": "Perez2019_CaseStudy_Base_LPP"}
CO2 = "CarbonDioxide"

# Decision variables per deck: (registry path, lower, upper, unit of the bounds).
# The cycle time is n_beds * t_feed (n_beds stays 3), so t_feed sets it. The base decks have
# no pressure set points left: the pump power and the blowdown share set the blowdown and
# evacuation pressures, the valves set the feed and pressurisation pressures.
VARIABLES: Dict[str, List[Tuple[str, float, float, str]]] = {
    "Base": [
        # Feed time, 1/3 of the cycle. The paper searched 20-200 s; below 40 s the pump's
        # third of the cycle is too short to evacuate a 0.3 m column. (Base 99 s.)
        ("variable[t_feed]", 40.0, 200.0, "s"),
        # Blower flow as a superficial velocity. The paper's 0.05-1.5 m/s interstitial range
        # is 0.02-0.6 superficial; above about 0.15 m/s the feed lifts the bed. (Base 0.085.)
        ("variable[feed_velocity]", 0.02, 0.15, "m/s"),
        # Vacuum pump power at full vacuum, W_max of W = W_max (1 - P/p_atm). Below 3 W the
        # bed is not evacuated in the time; above 30 W the evacuation passes any floor. (Base 10.4.)
        ("variable[vacuum_pump_power]", 3.0, 30.0, ""),
        # Share of the pump's time in blowdown: sets the intermediate pressure. (Base 0.197.)
        ("variable[blowdown_share]", 0.08, 0.6, ""),
        # Pressurisation valve from the feed header [mol/(s bar)]. It opens fully onto the evacuated
        # bed, so the first seconds set the approach to fluidisation: 0.001 still gives 1.15 (and
        # reaches 95 kPa), so the range reaches down to valves that leave the feed to finish the
        # pressurisation. Large: fast, lifts the bed at the feed end. (Base 0.006.)
        ("variable[cv_pressurisation]", 0.0003, 0.02, ""),
        # Vent valve [mol/(s bar)]: sets the feed pressure above 1 bar (a small valve gives a
        # pressure swing as well as a vacuum swing; the blower pays for it). (Base 0.5.)
        ("variable[cv_vent]", 0.05, 2.0, ""),
    ],
    "Base_LPP": [
        # Feed time (feed 1 + feed 2), 1/3 of the cycle. (Base 144.5 s.)
        ("variable[t_feed]", 60.0, 250.0, "s"),
        # Blower flow, superficial. (Base 0.0328 m/s.)
        ("variable[feed_velocity]", 0.015, 0.1, "m/s"),
        # Vacuum pump power at full vacuum [W]. (Base 5.8.)
        ("variable[vacuum_pump_power]", 2.0, 20.0, ""),
        # Share of the pump's time in blowdown. (Base 0.237.)
        ("variable[blowdown_share]", 0.08, 0.6, ""),
        # Share of the feed time in feed 1, which supplies the LPP; the LPP step is as long.
        # (Base 0.190.)
        ("variable[feed1_share]", 0.08, 0.5, ""),
        # LPP valve between the beds [mol/(s bar)]. Above about 0.008 the LPP draw exceeds what
        # the feed bed makes and gas flows back in through the vent (the backflow constraint).
        # (Base 0.008.)
        ("variable[cv_lpp]", 0.002, 0.016, ""),
        # Vent valve [mol/(s bar)]. (Base 0.5.)
        ("variable[cv_vent]", 0.05, 2.0, ""),
    ],
}
ALL_VARIABLES = sorted({path for v in VARIABLES.values() for path, *_ in v})

# The paper's operating points (Table 3, the paper's model), chosen from its purity-recovery
# fronts: the comparison for the "purity" fronts.
PAPER_POINTS = {
    "Base": {"A": (0.993, 0.600), "B": (0.952, 0.641), "C1": (0.834, 0.682)},
    "Base_LPP": {"D1": (0.961, 0.868), "E": (0.904, 0.911), "F": (0.799, 0.946)},
}

MAX_DP_PER_M = 10e3   # Pa/m, the default bound on the bed pressure drop per metre (--max-dp-per-m)
KWH_PER_TONNE = 3600.0   # J/kg in one kWh/t
FAILED = 10.0            # the normalised violation of each constraint of a failed run

CONSTRAINTS = ("purity", "recovery", "fluidisation", "pressure_drop", "evacuation_floor", "backflow", "converged")
PRACTICAL = ("fluidisation", "pressure_drop", "evacuation_floor")

COLUMNS = (["deck", "study", "seed", "index"] + ALL_VARIABLES
           + ["status", "css_reached", "css_cycle", "cycles_run", "healthy", "run_health", "warm", "start_cycles",
              "wall_time_s", "purity", "recovery", "energy_kwh_per_t", "productivity_kg_per_t_h",
              "max_fluidisation", "fluidisation_step", "max_dp_per_m", "dp_step", "min_evacuation_pressure",
              "backflow_fraction"]
           + [f"margin_{name}" for name in CONSTRAINTS] + ["feasible", "message"])

_base: Optional[sk.Simulation] = None
_cache = None
_warm = True


def base_case(deck: str, cycles: Optional[int], max_cycles: int) -> sk.Simulation:
    """The deck every candidate starts from: the base deck with a cycle cap (CSS stops a run earlier)."""
    sim = c.load(DECKS[deck])
    sim.cycle.number_of_cycles = cycles if cycles is not None else max_cycles
    return sim


def base_point(sim: sk.Simulation, deck: str) -> np.ndarray:
    """The base deck's own values of the decision variables, in the units of the bounds."""
    return np.array([sim.get(path, unit=unit) if unit else sim.get(path) for path, _, _, unit in VARIABLES[deck]])


def _init_worker(deck_json: str, lower: Sequence[float], upper: Sequence[float], warm: bool) -> None:
    global _base, _cache, _warm
    _base = sk.Simulation.from_json(deck_json)
    _cache = opt.StateCache(list(lower), list(upper))
    _warm = warm


def _finite(x: float) -> bool:
    return x is not None and isinstance(x, float) and math.isfinite(x)


def measure(sim: sk.Simulation, max_dp_per_m: float) -> Dict[str, object]:
    """The constraint values of the last logged cycle of a finished run."""
    out: Dict[str, object] = {}
    report = sim.operating_limits()
    fl = [(s.max_approach_to_fluidisation, s.step_name) for s in report.steps
          if math.isfinite(s.max_approach_to_fluidisation)]
    out["max_fluidisation"], out["fluidisation_step"] = max(fl) if fl else (float("nan"), "")
    vacuum = [s.min_pressure for s in report.steps if s.step_name in ("Blowdown", "Evacuation")]
    out["min_evacuation_pressure"] = min(vacuum) if vacuum else float("nan")
    # Gas in through an unconnected outlet, as a fraction of the cycle feed. The report names it only
    # above its threshold (1e-3, the bound used here): below it the fraction counts as 0.
    backflow = [w.value for w in report.warnings if w.limit == "outlet_backflow"]
    out["backflow_fraction"] = max(backflow) if backflow else 0.0
    dp = []
    for i, name in enumerate(sim.steps.names):
        p = sim.profile("pressure", cycle=-1, step=i)
        length = float(p.z[-1] - p.z[0])
        dp.append((float(np.abs(p.values[:, -1] - p.values[:, 0]).max()) / length, name))
    out["max_dp_per_m"], out["dp_step"] = max(dp)
    return out


def evaluate(x: Sequence[float], deck: str, targets: Dict[str, float], max_dp_per_m: float) -> Dict[str, object]:
    """Run one candidate to CSS and return its results-file row (without the bookkeeping columns)."""
    sim = _base.copy()
    for (path, _, _, unit), value in zip(VARIABLES[deck], x):
        if unit:
            sim.set(path, float(value), unit=unit)
        else:
            sim.set(path, float(value))
    row: Dict[str, object] = {path: float(v) for (path, *_), v in zip(VARIABLES[deck], x)}
    t0 = time.perf_counter()
    try:
        if _warm:
            result, info = opt.run_warm(sim, _cache, list(x))
            row.update(warm=bool(info["warm"]), start_cycles=info["start_cycles"])
        else:
            result = sim.run()
            row.update(warm=False, start_cycles=0)
    except Exception as exc:   # a deck the engine refuses: infeasible, not a crash of the search
        row.update(status="failed", message=str(exc)[:80], wall_time_s=time.perf_counter() - t0)
        return row
    row.update(status=result.status, message=result.message[:80].replace("\n", " "),
               css_reached=bool(result.css_reached), css_cycle=result.css_cycle, cycles_run=result.cycles_run,
               healthy=bool(result.healthy), run_health=" ".join(sorted({w.code for w in result.warnings})),
               wall_time_s=time.perf_counter() - t0)
    if result.status == "failed" or result.status == "cancelled":
        return row
    k = next((k for k in sim.kpis(cycle=-1) if k.component == CO2), None)
    if k is not None:
        row.update(purity=k.purity, recovery=k.recovery,
                   energy_kwh_per_t=k.electrical_energy / KWH_PER_TONNE,    # J/kg CO2 -> kWh/t CO2
                   productivity_kg_per_t_h=k.productivity_mass * 3.6e6)     # kg/(kg s) -> kg/(t h)
    row.update(measure(sim, max_dp_per_m))
    return row


def margins(row: Dict[str, object], targets: Dict[str, float], max_dp_per_m: float) -> Dict[str, float]:
    """Each constraint's margin, normalised by its bound (>= 0 is satisfied; NaN when not measured)."""
    def get(key):
        v = row.get(key)
        return float(v) if v not in (None, "") else float("nan")

    floor = targets["p_evac_floor"]
    return {
        "purity": (get("purity") - targets["purity_target"]) / targets["purity_target"],
        "recovery": (get("recovery") - targets["recovery_target"]) / targets["recovery_target"],
        "fluidisation": 1.0 - get("max_fluidisation"),
        "pressure_drop": (max_dp_per_m - get("max_dp_per_m")) / max_dp_per_m,
        "evacuation_floor": (get("min_evacuation_pressure") - floor) / floor,
        "backflow": (1e-3 - get("backflow_fraction")) / 1e-3,
        "converged": 0.0 if (str(row.get("css_reached")) == "True" and str(row.get("healthy")) == "True") else -1.0,
    }


def active_constraints(objectives: str, unconstrained: bool, backflow: bool = False) -> Tuple[str, ...]:
    names = ("purity", "recovery") if objectives == "energy" else ()
    if not unconstrained:
        names += PRACTICAL + (("backflow",) if backflow else ())
    return names + ("converged",)


def objectives_of(row: Dict[str, object], objectives: str) -> List[float]:
    """pymoo minimises: [energy, -productivity] or [-purity, -recovery]; the worst values when undefined."""
    def get(key, worst):
        v = row.get(key)
        v = float(v) if v not in (None, "") else float("nan")
        return v if math.isfinite(v) else worst

    if objectives == "energy":
        return [get("energy_kwh_per_t", 1e6), -get("productivity_kg_per_t_h", 0.0)]
    return [-get("purity", 0.0), -get("recovery", 0.0)]


def violations(m: Dict[str, float], names: Sequence[str]) -> List[float]:
    """pymoo inequality constraints G <= 0: the negated margins (FAILED when not measured)."""
    return [-m[n] if math.isfinite(m[n]) else FAILED for n in names]


# ---- the results file -------------------------------------------------------------------

def _key(deck: str, study: str, x: Sequence[float]) -> Tuple:
    return (deck, study) + tuple(round(float(v), 10) for v in x)


def read_results(path: Path) -> Dict[Tuple, Dict[str, str]]:
    if not path.exists():
        return {}
    done = {}
    with open(path, newline="") as f:
        for row in csv.DictReader(f):
            paths = [p for p, *_ in VARIABLES.get(row["deck"], [])]
            try:
                done[_key(row["deck"], row["study"], [float(row[p]) for p in paths])] = row
            except (KeyError, ValueError):
                continue
    return done


def append_result(path: Path, row: Dict[str, object]) -> None:
    new = not path.exists()
    path.parent.mkdir(parents=True, exist_ok=True)
    with open(path, "a", newline="") as f:
        w = csv.DictWriter(f, fieldnames=COLUMNS, extrasaction="ignore")
        if new:
            w.writeheader()
        w.writerow({k: _fmt(row.get(k, "")) for k in COLUMNS})


def _fmt(v: object) -> object:
    if isinstance(v, float):
        return "" if not math.isfinite(v) else repr(v)
    return v


# ---- the search -------------------------------------------------------------------------

class Search:
    """Evaluates batches of candidates of one deck and study: from the results file when they are
    there, otherwise in the process pool; appends every new evaluation to the file."""

    def __init__(self, deck: str, study: str, args: argparse.Namespace, targets: Dict[str, float], pool,
                 done: Dict[Tuple, Dict[str, str]]) -> None:
        self.deck, self.study, self.args, self.targets, self.pool, self.done = deck, study, args, targets, pool, done
        self.rows: List[Dict[str, object]] = []
        self.objectives = study.split("-")[0]
        self.constraints = active_constraints(self.objectives, study.endswith("unconstrained"),
                                              args.backflow_constraint)

    def run(self, X: np.ndarray) -> Tuple[np.ndarray, np.ndarray]:
        todo = [i for i, x in enumerate(X) if _key(self.deck, self.study, x) not in self.done]
        fresh = {}
        if todo:
            futures = [self.pool.submit(evaluate, X[i], self.deck, self.targets, self.args.max_dp_per_m)
                       for i in todo]
            for i, fut in zip(todo, futures):
                row = fut.result()
                row.update(deck=self.deck, study=self.study, seed=self.args.seed, index=len(self.rows) + len(fresh))
                m = margins(row, self.targets, self.args.max_dp_per_m)
                row.update({f"margin_{n}": v for n, v in m.items()})
                row["feasible"] = all(math.isfinite(m[n]) and m[n] >= 0 for n in self.constraints)
                append_result(self.args.results, row)
                self.done[_key(self.deck, self.study, X[i])] = {k: _fmt(v) for k, v in row.items()}
                fresh[i] = row
                if not self.args.quiet:
                    self._print(row)
        F, G = [], []
        for i, x in enumerate(X):
            row = fresh.get(i) or dict(self.done[_key(self.deck, self.study, x)])
            self.rows.append(row)
            m = margins(row, self.targets, self.args.max_dp_per_m)
            F.append(objectives_of(row, self.objectives))
            G.append(violations(m, self.constraints))
        return np.array(F), np.array(G)

    def _print(self, row: Dict[str, object]) -> None:
        xs = "  ".join(f"{p[9:-1]}={float(row[p]):.4g}" for p, *_ in VARIABLES[self.deck])
        print(f"  {xs}  -> {row.get('status', ''):11s} css@{row.get('css_cycle', '')} "
              f"pur {_pct(row.get('purity'))} rec {_pct(row.get('recovery'))} "
              f"E {_num(row.get('energy_kwh_per_t'), '.0f')} kWh/t  AF {_num(row.get('max_fluidisation'), '.2f')}  "
              f"Pmin {_num(row.get('min_evacuation_pressure'), '.0f', 1e-3)} kPa  "
              f"{'feasible' if row.get('feasible') else 'infeasible'}  ({float(row.get('wall_time_s', 0)):.0f} s"
              f"{', warm' if row.get('warm') else ''})", flush=True)


def _num(v, spec: str, scale: float = 1.0) -> str:
    try:
        v = float(v)
    except (TypeError, ValueError):
        return "-"
    return format(v * scale, spec) if math.isfinite(v) else "-"


def _pct(v) -> str:
    return _num(v, ".1f", 100.0) + "%"


def run_study(deck: str, study: str, args: argparse.Namespace, done) -> List[Dict[str, object]]:
    from pymoo.algorithms.moo.nsga2 import NSGA2
    from pymoo.core.problem import Problem
    from pymoo.optimize import minimize

    base = base_case(deck, args.cycles, args.max_cycles)
    targets = {name: base.get(f"variable[{name}]") for name in ("purity_target", "recovery_target", "p_evac_floor")}
    lower = np.array([v[1] for v in VARIABLES[deck]])
    upper = np.array([v[2] for v in VARIABLES[deck]])
    # The first population: the base deck's own point and a seeded uniform sample of the box.
    rng = np.random.default_rng(args.seed)
    X0 = np.vstack([np.clip(base_point(base, deck), lower, upper),
                    lower + rng.random((args.pop - 1, len(lower))) * (upper - lower)])

    with ProcessPoolExecutor(max_workers=args.workers, initializer=_init_worker,
                             initargs=(base.to_json(), lower, upper, not args.no_warm_start)) as pool:
        search = Search(deck, study, args, targets, pool, done)

        class Design(Problem):
            def __init__(self) -> None:
                super().__init__(n_var=len(lower), n_obj=2, n_ieq_constr=len(search.constraints), xl=lower, xu=upper)

            def _evaluate(self, X, out, *_, **__):
                out["F"], out["G"] = search.run(np.asarray(X))

        minimize(Design(), NSGA2(pop_size=args.pop, sampling=X0), ("n_gen", args.gen), seed=args.seed, verbose=False)
    return search.rows


def evaluate_base(deck: str, args: argparse.Namespace) -> Dict[str, object]:
    """The base deck's own point, through the same evaluation and results file (no pymoo needed)."""
    base = base_case(deck, args.cycles, args.max_cycles)
    targets = {name: base.get(f"variable[{name}]") for name in ("purity_target", "recovery_target", "p_evac_floor")}
    lower = [v[1] for v in VARIABLES[deck]]
    upper = [v[2] for v in VARIABLES[deck]]
    _init_worker(base.to_json(), lower, upper, False)
    row = evaluate(base_point(base, deck), deck, targets, args.max_dp_per_m)
    row.update(deck=deck, study="base", seed="", index=0)
    m = margins(row, targets, args.max_dp_per_m)
    row.update({f"margin_{n}": v for n, v in m.items()})
    row["feasible"] = all(math.isfinite(m[n]) and m[n] >= 0
                          for n in active_constraints("energy", False, args.backflow_constraint))
    append_result(args.results, row)
    return row


# ---- the summary ------------------------------------------------------------------------

def _f(row, key) -> float:
    try:
        return float(row.get(key))
    except (TypeError, ValueError):
        return float("nan")


def pareto(rows: List[Dict[str, object]], objectives: str) -> List[Dict[str, object]]:
    """The non-dominated feasible rows (unique), sorted by the first objective."""
    pts = [(objectives_of(r, objectives), r) for r in rows if str(r.get("feasible")) == "True"]
    front, seen = [], set()
    for f, r in pts:
        if any(all(g[j] <= f[j] for j in range(2)) and g != f for g, _ in pts):
            continue
        if tuple(f) in seen:
            continue
        seen.add(tuple(f))
        front.append((f, r))
    return [r for _, r in sorted(front, key=lambda t: t[0][0])]


def summarise(deck: str, study: str, rows: List[Dict[str, object]], quiet: bool) -> List[Dict[str, object]]:
    statuses = [str(r.get("status")) for r in rows]
    feasible = sum(1 for r in rows if str(r.get("feasible")) == "True")
    print(f"\n{deck} / {study}: {len(rows)} runs ({statuses.count('reached_css')} reached_css, "
          f"{statuses.count('completed')} completed, {statuses.count('failed')} failed), {feasible} feasible")
    objectives = study.split("-")[0]
    front = pareto(rows, objectives)
    feasible_front = bool(front)
    if not front:
        best = min(rows, key=lambda r: -sum(min(0.0, _f(r, f"margin_{n}")) if math.isfinite(_f(r, f"margin_{n}"))
                                            else -FAILED for n in active_constraints(objectives, "unconstrained" in study,
                                                                     "backflow" in study)), default=None)
        if best is not None:
            print("  no feasible candidate; the least-violating one:")
            front = [best]
    names = [p[9:-1] for p, *_ in VARIABLES[deck]]
    header = names + ["purity", "recovery", "kWh/t", "kg/(t h)", "AF", "dP kPa/m", "Pmin kPa"]
    c.table([[f"{_f(r, p):.4g}" for p, *_ in VARIABLES[deck]]
             + [_pct(r.get("purity")), _pct(r.get("recovery")), _num(r.get("energy_kwh_per_t"), ".0f"),
                _num(r.get("productivity_kg_per_t_h"), ".1f"), _num(r.get("max_fluidisation"), ".2f"),
                _num(r.get("max_dp_per_m"), ".1f", 1e-3), _num(r.get("min_evacuation_pressure"), ".1f", 1e-3)]
             for r in front], header)
    if objectives == "purity" and feasible_front:
        print("  against the paper's operating points (Table 3, model): best recovery at the same purity or higher")
        for name, (p, r) in PAPER_POINTS[deck].items():
            ok = [_f(x, "recovery") for x in front if _f(x, "purity") >= p]
            got = max(ok) if ok else float("nan")
            print(f"    {name:3s} purity {100 * p:.1f} %  recovery {100 * r:.1f} %  ->  "
                  f"{_pct(got) if math.isfinite(got) else 'not reached'}")
    return front if feasible_front else []


def plot(fronts: Dict[Tuple[str, str], Tuple[List, List]], args: argparse.Namespace) -> None:
    plt = c.pyplot(args)
    if plt is None:
        return
    fig, axes = plt.subplots(1, 2, figsize=(10, 4))
    for (deck, study), (rows, front) in fronts.items():
        if study.startswith("energy"):
            ax = axes[0]
            ax.scatter([_f(r, "productivity_kg_per_t_h") for r in front], [_f(r, "energy_kwh_per_t") for r in front],
                       label=deck, s=20)
        else:
            ax = axes[1]
            ax.scatter([_f(r, "recovery") for r in front], [_f(r, "purity") for r in front], s=20,
                       label=f"{deck} {study}")
    for deck in {d for d, _ in fronts}:
        pts = list(PAPER_POINTS[deck].values())
        axes[1].scatter([r for _, r in pts], [p for p, _ in pts], marker="x", color="k",
                        label="paper (Table 3)" if deck == "Base" else None)
    axes[0].set_xlabel("productivity [kg CO2 / (t h)]")
    axes[0].set_ylabel("specific energy [kWh / t CO2]")
    axes[0].set_title("feasible: purity, recovery, plant limits")
    axes[1].set_xlabel("CO2 recovery [-]")
    axes[1].set_ylabel("CO2 purity [-]")
    axes[1].set_title("purity-recovery fronts")
    for ax in axes:
        if ax.get_legend_handles_labels()[0]:
            ax.legend(fontsize=7)
    path = c.save(fig, "perez_2019_case_study_optimise.png")
    print(f"plot saved to {path}")


def main(argv: Optional[List[str]] = None) -> int:
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--deck", action="append", choices=list(DECKS), help="deck(s) (repeatable; default both)")
    ap.add_argument("--objectives", choices=("energy", "purity"), default="energy")
    ap.add_argument("--unconstrained", action="store_true",
                    help="purity objectives without the practical constraints (only CSS and run health)")
    ap.add_argument("--pop", type=int, default=12, help="population size (default 12)")
    ap.add_argument("--gen", type=int, default=5, help="generations (default 5)")
    ap.add_argument("--workers", type=int, default=os.cpu_count() or 1, help="parallel runs (default: all cores)")
    ap.add_argument("--seed", type=int, default=1)
    ap.add_argument("--max-cycles", type=int, default=200, help="cycle cap of a run; CSS stops it earlier (default 200)")
    ap.add_argument("--cycles", type=int, default=None, help="run exactly this many cycles (a quick look: no CSS)")
    ap.add_argument("--max-dp-per-m", type=float, default=MAX_DP_PER_M,
                    help=f"bound on the bed pressure drop per metre [Pa/m] (default {MAX_DP_PER_M:g})")
    ap.add_argument("--backflow-constraint", action="store_true",
                    help="make gas drawn back in through an unconnected outlet infeasible (always recorded)")
    ap.add_argument("--no-warm-start", action="store_true", help="run every candidate from the deck's initial state")
    ap.add_argument("--out", type=Path, default=None, help="directory of results.csv (default examples/output/case_study)")
    ap.add_argument("--base-only", action="store_true", help="evaluate the base deck's own point only (no pymoo)")
    ap.add_argument("--no-plot", action="store_true", help="do not save the plot")
    ap.add_argument("-q", "--quiet", action="store_true", help="only the summary")
    args = ap.parse_args(argv)
    if args.unconstrained and args.objectives != "purity":
        ap.error("--unconstrained applies to --objectives purity")
    args.results = (args.out or c.OUTPUT / "case_study") / "results.csv"
    decks = args.deck or list(DECKS)

    if args.base_only:
        for deck in decks:
            row = evaluate_base(deck, args)
            m = {n: _f(row, f"margin_{n}") for n in CONSTRAINTS}
            print(f"{deck}: {row.get('status')}  purity {_pct(row.get('purity'))}  recovery {_pct(row.get('recovery'))}  "
                  f"{_num(row.get('energy_kwh_per_t'), '.0f')} kWh/t  feasible {row['feasible']}")
            print("  margins: " + "  ".join(f"{n} {v:+.3f}" for n, v in m.items()))
        print(f"results appended to {args.results}")
        return 0

    try:
        import pymoo  # noqa: F401
    except ImportError:
        print("this example needs pymoo: pip install pymoo", file=sys.stderr)
        return 1
    study = args.objectives + ("-unconstrained" if args.unconstrained else "") \
        + ("-backflow" if args.backflow_constraint and not args.unconstrained else "")
    done = read_results(args.results)
    fronts = {}
    t0 = time.perf_counter()
    for deck in decks:
        rows = run_study(deck, study, args, done)
        fronts[(deck, study)] = (rows, summarise(deck, study, rows, args.quiet))
    print(f"\n{sum(len(r) for r, _ in fronts.values())} evaluations in {time.perf_counter() - t0:.0f} s on "
          f"{args.workers} worker(s); results in {args.results}")
    if not args.no_plot:
        plot(fronts, args)
    return 0


if __name__ == "__main__":
    sys.exit(main())
