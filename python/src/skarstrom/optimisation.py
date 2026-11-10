"""Single-objective optimisation of a deck, with warm-started candidates.

Searches registry paths within bounds for the lowest cost, where the cost is
built from the KPIs of each candidate's last cycle::

    import skarstrom
    from skarstrom import optimisation as opt

    sim = skarstrom.examples.Perez2019_A.load()
    result = opt.optimise(
        sim,
        {"step[Feed].time": (60, 150), "step[Evacuation].bottom.pressure_double_exponential_end": (2e3, 6e3)},
        objectives=[opt.maximise("CarbonDioxide", "recovery")],
        constraints=[opt.at_least("CarbonDioxide", "purity", 0.9, penalty=10)],
        evaluations=20,
    )
    result.x, result.cost, result.kpis      # the optimum, its cold cost and KPIs
    result.table()                         # every evaluation (pandas DataFrame)

The search is the engine's ``GlobalOptimizer`` (dlib's MaxLIPO with a trust
region; ``adsim/solver/optimisation.h``). Each candidate is run on a fresh copy
of the deck. With ``warm_start`` (default) it starts from the cyclic steady
state of the nearest candidate evaluated so far (Euclidean distance in the
variables scaled to [0, 1] by their bounds), which saves cycles; a state that
does not fit (a structural change) is skipped with a note. The optimum is then
re-run from the deck's initial conditions (``rerun_cold``, default) and that
cold value is the reported one; ``search_cost`` is the value the search saw.

The cost is minimised: ``+weight * KPI`` for :func:`minimise` terms,
``-weight * KPI`` for :func:`maximise` terms, plus ``penalty * violation`` per
constraint. A candidate whose run fails, or whose KPI is undefined, costs
1e10.

For a loop of your own (SciPy, Optuna, an ask/tell front end), use
:class:`StateCache` and :func:`run_warm` directly::

    cache = opt.StateCache(lower, upper)
    for x in candidates:
        cand = sim.copy(); ...set x...
        result, info = opt.run_warm(cand, cache, x)    # starts from the nearest cached state
"""

from __future__ import annotations

import math
from dataclasses import dataclass
from typing import TYPE_CHECKING, Any, Callable, Dict, List, Mapping, Optional, Sequence, Tuple

from . import _core, units
from ._core import StateCache

if TYPE_CHECKING:  # pragma: no cover
    from . import Simulation

__all__ = [
    "Constraint",
    "Objective",
    "OptimisationResult",
    "StateCache",
    "at_least",
    "at_most",
    "maximise",
    "minimise",
    "optimise",
    "run_warm",
]

KPIS = ("purity", "recovery", "productivity_mass", "productivity_volume", "energy_consumption",
        "electrical_energy", "production_rate", "thermal_energy", "cooling_duty")


@dataclass(frozen=True)
class Objective:
    """One objective term: a tracked component's KPI, minimised or maximised."""

    component: str
    kpi: str
    direction: str = "max"  # "min" or "max"
    weight: float = 1.0


@dataclass(frozen=True)
class Constraint:
    """One soft constraint: ``penalty * violation`` is added to the cost."""

    component: str
    kpi: str
    sense: str  # ">=" or "<="
    target: float
    penalty: float = 1.0


def minimise(component: str, kpi: str, weight: float = 1.0) -> Objective:
    return Objective(component, kpi, "min", float(weight))


def maximise(component: str, kpi: str, weight: float = 1.0) -> Objective:
    return Objective(component, kpi, "max", float(weight))


def at_least(component: str, kpi: str, target: float, penalty: float = 1.0) -> Constraint:
    return Constraint(component, kpi, ">=", float(target), float(penalty))


def at_most(component: str, kpi: str, target: float, penalty: float = 1.0) -> Constraint:
    return Constraint(component, kpi, "<=", float(target), float(penalty))


def run_warm(sim: "Simulation", cache: StateCache, x: Sequence[float], cycles: Optional[int] = None,
             progress: Optional[Callable] = None, log_level: str = "low",
             on_message: Optional[Callable] = None) -> Tuple[Any, Dict[str, Any]]:
    """Run ``sim`` from the state in ``cache`` nearest to ``x`` and cache its
    final state under ``x``. Returns ``(RunResult, info)``; ``info`` has
    ``warm``, ``distance``, ``start_cycles``, ``skipped``, ``retried_cold``,
    ``cached`` and ``note`` (why a warm start was skipped or repeated cold).
    ``x`` is in the units of the cache's bounds (or already scaled for a cache
    made without bounds)."""
    return _core.run_with_state_cache(sim._sim, cache, [float(v) for v in x], cycles, progress, log_level,
                                      on_message)


def _nan(x: Any) -> Any:
    return None if isinstance(x, float) and not math.isfinite(x) else x


class OptimisationResult:
    """The outcome of :func:`optimise`."""

    def __init__(self, paths: List[str], raw: Dict[str, Any]) -> None:
        self.paths = paths
        self.raw = raw
        #: The optimum: path -> value (SI units).
        self.x: Dict[str, float] = dict(zip(paths, raw["x"]))
        #: The optimum's cost: the cold re-run's, or the search's without it.
        self.cost: float = raw["cost"]
        #: The optimum's cost as the (warm-started) search evaluated it.
        self.search_cost: float = raw["search_cost"]
        self.best_index: int = raw["best_index"]
        #: The cold re-run of the optimum (a dict like each evaluation), or None.
        self.cold: Optional[Dict[str, Any]] = raw["cold"]
        self.evaluations: List[Dict[str, Any]] = raw["evaluations"]
        self.cancelled: bool = raw["cancelled"]

    @property
    def kpis(self) -> Dict[str, Dict[str, float]]:
        """The optimum's KPIs per referenced component (from the cold re-run when there is one)."""
        if self.cold is not None:
            return self.cold["kpis"]
        if 0 <= self.best_index < len(self.evaluations):
            return self.evaluations[self.best_index]["kpis"]
        return {}

    def apply(self, sim: "Simulation") -> None:
        """Set the optimum's values in ``sim``."""
        for path, value in self.x.items():
            sim.set(path, value)

    def table(self):
        """Every evaluation, then the cold re-run (index -1): values, cost, KPIs, warm start, cycles,
        wall time (pandas DataFrame)."""
        import pandas as pd

        rows = []
        for e in self.evaluations + ([self.cold] if self.cold is not None else []):
            row: Dict[str, Any] = {"index": e["index"]}
            row.update(dict(zip(self.paths, e["x"])))
            row["cost"] = e["cost"]
            for comp, values in e["kpis"].items():
                for k, v in values.items():
                    row[f"{comp}.{k}"] = _nan(v)
            row.update(status=e["status"], accepted=e["accepted"], warm=e["warm"], distance=_nan(e["distance"]),
                       cycles_run=e["cycles_run"], css_reached=e["css_reached"], wall_time_s=e["wall_time_s"],
                       note=e["note"], message=e["message"])
            rows.append(row)
        return pd.DataFrame(rows)

    def __repr__(self) -> str:
        return (f"<OptimisationResult cost={self.cost:.6g} (search {self.search_cost:.6g}), "
                f"{len(self.evaluations)} evaluations, x={self.x}>")


def _bound(sim: "Simulation", path: str, value: Any) -> float:
    if units.is_quantity(value):
        info = sim._index().find(path)
        if info is not None and info.unit_checked:
            return units._quantity_to(value, info.unit, path)
        return float(value.to_base_units().magnitude)
    return float(value)


def optimise(
    sim: "Simulation",
    variables: Mapping[str, Tuple[Any, Any]],
    objectives: Sequence[Objective],
    constraints: Sequence[Constraint] = (),
    evaluations: int = 35,
    init_samples: int = 5,
    warm_start: bool = True,
    rerun_cold: bool = True,
    require_css: bool = False,
    cycles: Optional[int] = None,
    settings: Optional[Mapping[str, Any]] = None,
    cache_capacity: int = 64,
    progress: Optional[Callable[[Dict[str, Any]], Any]] = None,
) -> OptimisationResult:
    """Minimise the cost of ``objectives`` and ``constraints`` over ``variables``.

    ``variables`` maps registry paths to ``(lower, upper)`` bounds (SI numbers
    or pint Quantities). ``evaluations`` is the total number of runs of the
    search (``init_samples`` of them spread over the box first), plus one
    cold re-run of the optimum. ``warm_start`` starts each candidate from the
    nearest cached cyclic steady state; ``rerun_cold`` re-runs the optimum
    from the deck's initial conditions and reports that value. ``require_css``
    penalises candidates that stop at the cycle limit; ``cycles`` overrides
    the deck's cycle limit; ``settings`` (path -> value) are applied to the
    deck before any run. ``progress(evaluation)`` is called after every run
    with the evaluation's dict. ``sim`` is not changed (see
    :meth:`OptimisationResult.apply`). Ctrl-C stops the search.
    """
    if not variables:
        raise ValueError("an optimisation needs at least one variable")
    if not objectives:
        raise ValueError("an optimisation needs at least one objective")
    for term in list(objectives) + list(constraints):
        if term.kpi not in KPIS:
            raise ValueError(f"unknown KPI {term.kpi!r} ({', '.join(KPIS)})")
    paths = list(variables)
    var_list = []
    for path in paths:
        lo, hi = variables[path]
        var_list.append((path, _bound(sim, path, lo), _bound(sim, path, hi)))
    setting_text = {}
    for path, value in (settings or {}).items():
        if units.is_quantity(value):
            setting_text[path] = f"{float(value.magnitude)!r} [{value.units:~}]"
        elif isinstance(value, bool):
            setting_text[path] = "1" if value else "0"
        elif isinstance(value, float):
            setting_text[path] = repr(value)
        else:
            setting_text[path] = str(value)
    options = {
        "evaluations": int(evaluations),
        "init_samples": int(init_samples),
        "warm_start": bool(warm_start),
        "rerun_cold": bool(rerun_cold),
        "require_css": bool(require_css),
        "cycles": cycles,
        "cache_capacity": int(cache_capacity),
        "settings": setting_text,
    }
    raw = _core.optimise(
        sim.to_json(), var_list,
        [{"component": o.component, "kpi": o.kpi, "direction": o.direction, "weight": o.weight} for o in objectives],
        [{"component": c.component, "kpi": c.kpi, "sense": c.sense, "target": c.target, "penalty": c.penalty}
         for c in constraints],
        options, progress)
    return OptimisationResult(paths, raw)
