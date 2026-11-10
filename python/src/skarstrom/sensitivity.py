"""Sensitivity analysis and KPI error bars.

Which uncertain inputs move a deck's KPIs, and how uncertain do they make
them? A study has inputs (registry paths, each with a distribution) and
outputs (``"Component.kpi"``, the KPIs at the last cycle), and answers in two
steps::

    import skarstrom
    from skarstrom import sensitivity as sa

    sim = skarstrom.examples.Perez2019_A.load()
    study = sa.Study(sim, [
        sa.relative("step[Feed].time", 0.10),                        # +-10 %, uniform
        sa.normal("adsorbent[1].isotherm[CarbonDioxide].b1", 3.934e-7, 0.2e-7, lower=0),
        sa.uniform("layer[1].bed_void_fraction", 0.39, 0.43),
    ], workers=4)
    screen = study.screen()               # Morris: which inputs matter (r (k + 1) runs)
    screen.table()                        # mu*, mu, sigma per output and input
    spread = study.propagate(inputs=screen.important())   # Latin hypercube Monte Carlo
    spread.table()                        # mean, sd, p05 / p50 / p95 (the error bar), converged
    spread.shares()                       # each input's share of the spread

A study's set-up (inputs, outputs, run policy, seed and sizes, not its
results) can be saved with the deck and read back::

    study.save_to(sim); sim.save("deck.json")
    study = sa.Study.from_deck(skarstrom.Simulation("deck.json"), workers=4)

The engine plans the runs and analyses them (``adsim/solver/sensitivity.h``);
this module runs them, in parallel with ``workers > 1`` (separate processes).
Every run must reach cyclic steady state (``require_css``); runs that fail or
do not settle are counted and shown, never silently dropped. Ctrl-C stops a
study. Full guide: https://samaffleck.github.io/skarstrom/sensitivity/
"""

from __future__ import annotations

import json
import math
import os
from concurrent.futures import ProcessPoolExecutor
from dataclasses import asdict, dataclass, replace
from typing import TYPE_CHECKING, Any, Callable, Dict, List, Mapping, Optional, Sequence, Union

from . import _core, units

if TYPE_CHECKING:  # pragma: no cover
    from . import Simulation

__all__ = [
    "Input",
    "Propagation",
    "Screening",
    "Study",
    "kpi_names",
    "step_output_names",
    "normal",
    "relative",
    "triangular",
    "uniform",
]

Number = Union[float, int, Any]  # a number in SI units or a pint Quantity
Progress = Callable[[int, int], Any]


@dataclass(frozen=True)
class Input:
    """One uncertain input: a parameter path and its distribution.

    Build one with :func:`uniform`, :func:`relative`, :func:`triangular` or
    :func:`normal`. Numbers are in the parameter's SI unit; a pint Quantity is
    converted when the study is created.
    """

    path: str
    distribution: str = "uniform"
    lower: Optional[Number] = None
    upper: Optional[Number] = None
    mode: Optional[Number] = None
    mean: Optional[Number] = None
    sd: Optional[Number] = None
    relative: Optional[float] = None  # +- fraction of the deck's value (resolved by Study)

    def as_dict(self) -> Dict[str, Any]:
        """The input as the engine takes it (needs plain numbers)."""
        d: Dict[str, Any] = {"path": self.path, "distribution": self.distribution}
        for key in ("lower", "upper", "mode", "mean", "sd"):
            value = getattr(self, key)
            d[key] = None if value is None else float(value)
        return d

    def quantile(self, u: Union[float, Sequence[float]]) -> Any:
        """The value(s) at probability ``u`` in (0, 1)."""
        import numpy as np

        return _core.sensitivity_quantile(self.as_dict(), list(np.atleast_1d(np.asarray(u, dtype=float))))


def uniform(path: str, lower: Number, upper: Number) -> Input:
    """Equally likely anywhere in [lower, upper]."""
    return Input(path, "uniform", lower=lower, upper=upper)


def relative(path: str, fraction: float) -> Input:
    """Uniform over the deck's value +- ``fraction`` (0.1 for +-10 %)."""
    return Input(path, "uniform", relative=float(fraction))


def triangular(path: str, lower: Number, mode: Number, upper: Number) -> Input:
    """Between lower and upper, most likely at ``mode``."""
    return Input(path, "triangular", lower=lower, mode=mode, upper=upper)


def normal(path: str, mean: Number, sd: Number, lower: Optional[Number] = None,
           upper: Optional[Number] = None) -> Input:
    """A best value and its standard deviation (an isotherm fit's standard error,
    say), optionally truncated to [lower, upper] (e.g. ``lower=0`` for a positive
    property)."""
    return Input(path, "normal", lower=lower, upper=upper, mean=mean, sd=sd)


def kpi_names() -> List[str]:
    """The KPIs an output can name (``"CarbonDioxide.purity"``...)."""
    return list(_core.sensitivity_kpi_names())


def step_output_names() -> List[str]:
    """The step outputs an output can name after a step reference
    (``"step[Evacuation].end_pressure_bottom"``; the step by deck ID or name):
    the pressure the step ends at, at each bed end and over the bed
    (:meth:`Simulation.step_end_pressures`)."""
    return list(_core.step_output_names())


# ---- running ------------------------------------------------------------------------

# Warm-start caches of this process, one per (deck, input paths, policy): the
# cyclic steady states of the runs it made, keyed by their probability
# coordinates. Worker processes keep their own. Only the most recent few are
# kept, so a session running many studies does not hold every state.
_CACHES: Dict[Any, Any] = {}
_MAX_CACHES = 4


def _cache_for(deck: str, inputs: List[Dict[str, Any]], policy: Dict[str, Any]) -> Any:
    key = (hash(deck), tuple(i["path"] for i in inputs), json.dumps(policy, sort_keys=True, default=str))
    cache = _CACHES.pop(key, None)
    if cache is None:
        cache = _core.StateCache()
    _CACHES[key] = cache  # most recent last
    while len(_CACHES) > _MAX_CACHES:
        _CACHES.pop(next(iter(_CACHES)))
    return cache


def _evaluate(deck: str, inputs: List[Dict[str, Any]], values: List[float], outputs: List[str],
              policy: Dict[str, Any], key: Optional[List[float]] = None) -> Dict[str, Any]:
    """One sample (top level, so worker processes can run it). With ``key``
    (the sample's probability coordinates) the run warm-starts from the
    nearest state this process has cached for the same study."""
    if key:
        cache = _cache_for(deck, inputs, policy)
        return _core.sensitivity_evaluate(deck, inputs, values, outputs, policy, cache, list(key))
    return _core.sensitivity_evaluate(deck, inputs, values, outputs, policy)


def _nan_to_none(x: Any) -> Any:
    if isinstance(x, float) and not math.isfinite(x):
        return None
    if isinstance(x, dict):
        return {k: _nan_to_none(v) for k, v in x.items()}
    if isinstance(x, (list, tuple)):
        return [_nan_to_none(v) for v in x]
    return x


class Study:
    """A sensitivity study of one deck.

    ``inputs`` are :class:`Input` objects; ``outputs`` are ``"Component.kpi"``
    names (default: purity, recovery, mass productivity and specific energy of
    every component the deck tracks; see :func:`kpi_names`) or step outputs
    ``"step[<ID or name>].end_pressure_bottom"`` (see :func:`step_output_names`).

    Every run starts from the deck as it is now. ``settings`` are applied to
    every run (paths and values as for :meth:`Simulation.set`), typically a
    tighter CSS tolerance: a run that is not fully settled wobbles, and the
    wobble would read as sensitivity. ``require_css`` (default) only counts
    runs that reach cyclic steady state; ``cycles`` overrides the deck's
    cycle limit. ``workers`` runs that many processes. ``seed`` makes the
    study repeatable.

    ``warm_start`` (default) starts each run from the cyclic steady state of
    the nearest run made so far (in the inputs' probability coordinates; the
    nominal run counts, at the median), which saves cycles. Warm and cold
    runs stop at the same cyclic steady state to within the CSS tolerance's
    bound, not bit for bit; with ``workers > 1`` each process keeps its own
    states, so which run starts from which depends on scheduling. Set
    ``warm_start=False`` for runs that are exactly repeatable. Each sample
    records ``warm``, ``cycles_run`` and ``note`` (why a warm start was
    skipped). ``trajectories``, ``levels`` and ``threshold`` are the
    screening's defaults, ``runs``, ``batch``, ``max_runs`` and ``rel_tol``
    the propagation's, and ``bootstrap`` both analyses' (see :meth:`screen`
    and :meth:`propagate`).

    :meth:`from_deck` builds the study saved with a deck and :meth:`save_to`
    saves this one's set-up into a deck (not its results).
    """

    def __init__(
        self,
        sim: "Simulation",
        inputs: Sequence[Input],
        outputs: Optional[Sequence[str]] = None,
        settings: Optional[Mapping[str, Any]] = None,
        require_css: bool = True,
        cycles: Optional[int] = None,
        seed: int = 1,
        workers: int = 1,
        trajectories: int = 10,
        levels: int = 4,
        threshold: float = 0.1,
        runs: int = 100,
        batch: int = 50,
        max_runs: int = 400,
        rel_tol: float = 0.1,
        bootstrap: int = 1000,
        warm_start: bool = True,
    ) -> None:
        if not inputs:
            raise ValueError("a study needs at least one input")
        self.inputs: List[Input] = [self._resolve(sim, i) for i in inputs]
        self.settings: Dict[str, str] = {k: self._setting_text(sim, k, v) for k, v in (settings or {}).items()}
        self.require_css = bool(require_css)
        self.cycles = None if cycles is None else int(cycles)
        self.seed = int(seed)
        self.workers = max(1, int(workers))
        self.trajectories = int(trajectories)
        self.levels = int(levels)
        self.threshold = float(threshold)
        self.runs = int(runs)
        self.batch = int(batch)
        self.max_runs = int(max_runs)
        self.rel_tol = float(rel_tol)
        self.bootstrap = int(bootstrap)
        self.warm_start = bool(warm_start)
        self.deck = sim.to_json()
        self.title = sim.get_text("info.title")
        self._named_outputs: List[str] = list(outputs) if outputs else []
        self.outputs: List[str] = self._named_outputs or list(_core.sensitivity_default_outputs(self.deck))
        _core.sensitivity_check(self.deck, self._input_dicts(), self.outputs, self.policy)
        self._nominal: Optional[Dict[str, Any]] = None

    # ---- saved with the deck ----------------------------------------------------------
    @classmethod
    def from_deck(cls, sim: "Simulation", workers: int = 1, **overrides: Any) -> "Study":
        """The study saved with ``sim``'s deck. ``overrides`` replace saved
        values (any constructor argument, e.g. ``seed=2`` or
        ``require_css=False``). Raises ValueError when the deck has no study."""
        saved = sim._sim.sensitivity_study()
        if not saved["inputs"]:
            raise ValueError("the deck has no sensitivity study (save one with Study.save_to)")
        policy = saved["policy"]
        kwargs: Dict[str, Any] = {
            "inputs": [Input(**i) for i in saved["inputs"]],
            "outputs": saved["outputs"] or None,
            "settings": policy["settings"],
            "require_css": policy["require_css"],
            "cycles": policy["cycles"],
            **{k: saved[k] for k in ("seed", "trajectories", "levels", "threshold", "runs", "batch", "max_runs",
                                      "rel_tol", "bootstrap")},
        }
        kwargs.update(overrides)
        return cls(sim, workers=workers, **kwargs)

    def setup(self) -> Dict[str, Any]:
        """The study's set-up as saved with a deck (see :meth:`save_to`)."""
        return {
            "inputs": [i.as_dict() for i in self.inputs],
            "outputs": list(self._named_outputs),
            "policy": self.policy,
            "seed": self.seed,
            "trajectories": self.trajectories,
            "levels": self.levels,
            "threshold": self.threshold,
            "runs": self.runs,
            "batch": self.batch,
            "max_runs": self.max_runs,
            "rel_tol": self.rel_tol,
            "bootstrap": self.bootstrap,
        }

    def save_to(self, sim: "Simulation") -> None:
        """Save this study's set-up into ``sim``'s deck (then ``sim.save(path)``
        writes it to a file). Outputs left to the default stay the default, so
        they follow the deck's tracked components. Results are not saved."""
        sim._sim.set_sensitivity_study(self.setup())

    # ---- set-up --------------------------------------------------------------------
    @staticmethod
    def _resolve(sim: "Simulation", inp: Input) -> Input:
        info = sim._index().find(inp.path)
        if info is None:
            raise KeyError(f"input {inp.path!r}: no such parameter in the deck (see sim.parameters())")
        if inp.relative is not None:
            nominal = float(sim.get(inp.path))
            f = inp.relative
            if not (f > 0) or nominal == 0 or not math.isfinite(nominal):
                raise ValueError(f"input {inp.path!r}: a relative range needs a fraction > 0 and a non-zero "
                                 f"value (the deck has {nominal}); give bounds instead")
            a, b = nominal * (1 - f), nominal * (1 + f)
            return Input(inp.path, "uniform", lower=min(a, b), upper=max(a, b))

        def si(value: Any, key: str) -> Optional[float]:
            if value is None:
                return None
            if units.is_quantity(value):
                if key == "sd" and str(value.units) in ("degree_Celsius", "degree_Fahrenheit"):
                    raise ValueError(f"input {inp.path!r}: give a temperature spread as a difference "
                                     "(ureg.delta_degC) or a plain number in K")
                if not info.unit_checked:
                    return float(value.to_base_units().magnitude)
                return units._quantity_to(value, info.unit, f"{inp.path} ({key})")
            return float(value)

        return replace(inp, **{k: si(getattr(inp, k), k) for k in ("lower", "upper", "mode", "mean", "sd")})

    @staticmethod
    def _setting_text(sim: "Simulation", path: str, value: Any) -> str:
        if sim._index().find(path) is None:
            raise KeyError(f"setting {path!r}: no such parameter in the deck")
        if units.is_quantity(value):
            return f"{float(value.magnitude)!r} [{value.units:~}]"
        if isinstance(value, bool):
            return "1" if value else "0"
        if isinstance(value, float):
            return repr(value)
        return str(value)

    @property
    def policy(self) -> Dict[str, Any]:
        return {"require_css": self.require_css, "cycles": self.cycles, "settings": dict(self.settings)}

    def _input_dicts(self, inputs: Optional[Sequence[Input]] = None) -> List[Dict[str, Any]]:
        return [i.as_dict() for i in (inputs if inputs is not None else self.inputs)]

    def _select(self, inputs: Optional[Sequence[Union[str, Input]]]) -> List[Input]:
        if inputs is None:
            return list(self.inputs)
        by_path = {i.path: i for i in self.inputs}
        out = []
        for i in inputs:
            path = i.path if isinstance(i, Input) else i
            if path not in by_path:
                raise KeyError(f"{path!r} is not an input of this study")
            out.append(by_path[path])
        if not out:
            raise ValueError("no inputs selected")
        return out

    def describe(self) -> Dict[str, Any]:
        """The study's set-up (inputs, outputs, run policy, seed)."""
        return {
            "title": self.title,
            "inputs": [_nan_to_none(asdict(i)) for i in self.inputs],
            "outputs": list(self.outputs),
            "policy": self.policy,
            "seed": self.seed,
            "screening": {"trajectories": self.trajectories, "levels": self.levels, "threshold": self.threshold},
            "propagation": {"runs": self.runs, "batch": self.batch, "max_runs": self.max_runs,
                            "rel_tol": self.rel_tol},
            "bootstrap": self.bootstrap,
            "warm_start": self.warm_start,
        }

    # ---- running ---------------------------------------------------------------------
    def _run(self, inputs: Sequence[Input], rows: List[List[float]], progress: Optional[Progress],
             offset: int = 0, total: Optional[int] = None,
             keys: Optional[List[List[float]]] = None) -> List[Dict[str, Any]]:
        args = (self.deck, self._input_dicts(inputs))
        total = total if total is not None else len(rows)
        if not self.warm_start or keys is None:
            keys = [None] * len(rows)  # type: ignore[list-item]
        out: List[Dict[str, Any]] = []
        if self.workers == 1 or len(rows) == 1:
            for i, (values, key) in enumerate(zip(rows, keys)):
                out.append(_evaluate(*args, values, self.outputs, self.policy, key))
                if progress:
                    progress(offset + i + 1, total)
            return out
        with ProcessPoolExecutor(max_workers=min(self.workers, len(rows))) as pool:
            futures = [pool.submit(_evaluate, *args, values, self.outputs, self.policy, key)
                       for values, key in zip(rows, keys)]
            try:
                for i, f in enumerate(futures):
                    out.append(f.result())
                    if progress:
                        progress(offset + i + 1, total)
            except BaseException:
                for f in futures:
                    f.cancel()
                raise
        return out

    def nominal(self) -> Dict[str, Any]:
        """The unperturbed run (with the study's settings): status, message,
        accepted, outputs (dict by name) and wall_time_s. Run once and cached."""
        if self._nominal is None:
            key = [0.5] * len(self.inputs) if self.warm_start else None
            r = _evaluate(self.deck, self._input_dicts(), [], self.outputs, self.policy, key)
            r["outputs"] = dict(zip(self.outputs, r["outputs"]))
            self._nominal = r
        return self._nominal

    def screen(self, trajectories: Optional[int] = None, levels: Optional[int] = None,
               inputs: Optional[Sequence[Union[str, Input]]] = None, progress: Optional[Progress] = None,
               bootstrap: Optional[int] = None) -> "Screening":
        """Morris screening: ``trajectories * (inputs + 1)`` runs.

        Sizes left out are the study's (``Study(trajectories=10, levels=4)``
        by default). ``progress(done, total)`` is called after every run.
        """
        trajectories = self.trajectories if trajectories is None else trajectories
        levels = self.levels if levels is None else levels
        bootstrap = self.bootstrap if bootstrap is None else bootstrap
        chosen = self._select(inputs)
        plan = _core.sensitivity_morris_plan(self._input_dicts(chosen), trajectories, levels, self.seed)
        samples = self._run(chosen, plan["values"], progress, keys=plan["u"])
        y = [s["outputs"] for s in samples]
        result = _core.sensitivity_analyse_morris(self._input_dicts(chosen), plan, self.outputs, y, bootstrap)
        return Screening(self, chosen, plan, samples, result)

    def propagate(self, runs: Optional[int] = None, batch: Optional[int] = None, max_runs: Optional[int] = None,
                  rel_tol: Optional[float] = None, inputs: Optional[Sequence[Union[str, Input]]] = None,
                  progress: Optional[Progress] = None, bootstrap: Optional[int] = None) -> "Propagation":
        """Monte Carlo propagation with Latin hypercube sampling.

        Starts with ``runs`` runs and adds batches of ``batch`` until every
        output's percentiles are known to within ``rel_tol`` of its 5-95 %
        band, or ``max_runs`` is reached (then ``converged`` is False).
        Inputs left out (``inputs``, e.g. :meth:`Screening.important`) stay at
        the deck's value. ``progress(done, total)`` counts runs against
        ``max_runs``. Sizes left out are the study's (``runs=100, batch=50,
        max_runs=400, rel_tol=0.1`` by default).
        """
        runs = self.runs if runs is None else runs
        batch = self.batch if batch is None else batch
        max_runs = self.max_runs if max_runs is None else max_runs
        rel_tol = self.rel_tol if rel_tol is None else rel_tol
        bootstrap = self.bootstrap if bootstrap is None else bootstrap
        if runs < 2 or batch < 1 or max_runs < runs:
            raise ValueError("need runs >= 2, batch >= 1 and max_runs >= runs")
        chosen = self._select(inputs)
        dicts = self._input_dicts(chosen)
        values: List[List[float]] = []
        samples: List[Dict[str, Any]] = []
        batches = 0
        n = runs
        while True:
            plan = _core.sensitivity_monte_carlo_plan(dicts, n, self.seed, batches)
            samples += self._run(chosen, plan["values"], progress, offset=len(values), total=max_runs,
                                 keys=plan["u"])
            values += plan["values"]
            batches += 1
            result = _core.sensitivity_analyse_monte_carlo(
                dicts, values, self.outputs, [s["outputs"] for s in samples], self.seed, bootstrap, rel_tol)
            if result["converged"] or len(values) >= max_runs:
                break
            n = min(batch, max_runs - len(values))
        return Propagation(self, chosen, values, samples, result, batches, rel_tol)


class _StudyResult:
    def __init__(self, study: Study, inputs: List[Input], samples: List[Dict[str, Any]],
                 values: List[List[float]], result: Dict[str, Any]) -> None:
        self.study = study
        self.inputs = inputs
        self.samples = samples
        self.values = values
        self.result = result

    @property
    def runs(self) -> int:
        return len(self.samples)

    @property
    def failed(self) -> int:
        """Runs that did not count (failed, or did not reach CSS)."""
        return sum(1 for s in self.samples if not s["accepted"])

    @property
    def warnings(self) -> List[str]:
        out = []
        if self.failed:
            reasons: Dict[str, int] = {}
            for s in self.samples:
                if not s["accepted"]:
                    key = s["message"] or s["status"]
                    reasons[key] = reasons.get(key, 0) + 1
            detail = "; ".join(f"{n} x {m}" for m, n in sorted(reasons.items(), key=lambda x: -x[1]))
            out.append(f"{self.failed} of {self.runs} runs did not count ({detail}). Failures often cluster in one "
                       "corner of the ranges, so the results describe the rest of it.")
        return out

    def samples_table(self):
        """Every run: input values, outputs, status and message (pandas DataFrame)."""
        import pandas as pd

        rows = []
        for values, s in zip(self.values, self.samples):
            row = {i.path: v for i, v in zip(self.inputs, values)}
            row.update(dict(zip(self.study.outputs, s["outputs"])))
            row.update(status=s["status"], accepted=s["accepted"], message=s["message"], wall_time_s=s["wall_time_s"],
                       cycles_run=s.get("cycles_run"), warm=s.get("warm"), note=s.get("note"))
            rows.append(row)
        return pd.DataFrame(rows)

    def to_dict(self) -> Dict[str, Any]:
        raise NotImplementedError

    def save(self, path: Union[str, os.PathLike]) -> None:
        """Write the study, every run and the analysis as JSON (NaN as null)."""
        with open(path, "w", encoding="utf-8") as f:
            json.dump(_nan_to_none(self.to_dict()), f, indent=1)


class Screening(_StudyResult):
    """Morris screening results. ``mu*`` ranks the inputs (in the output's units:
    the typical change when the input moves across its range); a ``sigma``
    comparable to ``mu*`` means the input's effect depends on the others or is
    curved."""

    def __init__(self, study: Study, inputs: List[Input], plan: Dict[str, Any], samples: List[Dict[str, Any]],
                 result: Dict[str, Any]) -> None:
        super().__init__(study, inputs, samples, plan["values"], result)
        self.plan = plan

    def table(self):
        """mu*, its 95 % interval, mu and sigma per output and input (pandas DataFrame)."""
        import pandas as pd

        rows = []
        for o in self.result["outputs"]:
            for e in o["effects"]:
                rows.append({
                    "output": o["output"], "input": e["input"], "mu_star": e["mu_star"],
                    "mu_star_low": e["mu_star_interval"]["low"], "mu_star_high": e["mu_star_interval"]["high"],
                    "mu": e["mu"], "sigma": e["sigma"], "effects": e["effects"],
                })
        return pd.DataFrame(rows)

    def ranking(self, output: Optional[str] = None) -> List[tuple]:
        """(input, mu*) pairs, most influential first, for one output (default: the first)."""
        name = output or self.study.outputs[0]
        for o in self.result["outputs"]:
            if o["output"] == name:
                effects = [(e["input"], e["mu_star"]) for e in o["effects"]]
                return sorted(effects, key=lambda x: -x[1] if math.isfinite(x[1]) else math.inf)
        raise KeyError(f"{name!r} is not an output of this study")

    def important(self, threshold: Optional[float] = None) -> List[str]:
        """Inputs whose mu* reaches ``threshold`` (default: the study's, 0.1
        unless set) times the largest mu* of at least one output: the ones to
        keep for :meth:`Study.propagate`."""
        threshold = self.study.threshold if threshold is None else threshold
        keep = set()
        for o in self.result["outputs"]:
            finite = [e["mu_star"] for e in o["effects"] if math.isfinite(e["mu_star"])]
            if not finite or max(finite) <= 0:
                continue
            top = max(finite)
            keep.update(e["input"] for e in o["effects"]
                        if math.isfinite(e["mu_star"]) and e["mu_star"] >= threshold * top)
        return [i.path for i in self.inputs if i.path in keep]

    def to_dict(self) -> Dict[str, Any]:
        return {"method": "morris", "study": self.study.describe(), "inputs": [i.path for i in self.inputs],
                "plan": self.plan, "samples": self.samples, "result": self.result, "warnings": self.warnings}

    def __repr__(self) -> str:
        return f"<Screening {len(self.inputs)} inputs, {self.runs} runs, {self.failed} failed>"


class Propagation(_StudyResult):
    """Monte Carlo results: per output the mean, sd and the 5th / 50th / 95th
    percentiles (the 90 % band is the error bar), and each input's share of
    the spread."""

    def __init__(self, study: Study, inputs: List[Input], values: List[List[float]], samples: List[Dict[str, Any]],
                 result: Dict[str, Any], batches: int, rel_tol: float) -> None:
        super().__init__(study, inputs, samples, values, result)
        self.batches = batches
        self.rel_tol = rel_tol

    @property
    def converged(self) -> bool:
        return bool(self.result["converged"])

    def table(self, nominal: bool = False):
        """One row per output (pandas DataFrame). ``nominal=True`` adds the
        unperturbed run's value (one extra run)."""
        import pandas as pd

        base = self.study.nominal()["outputs"] if nominal else {}
        rows = []
        for o in self.result["outputs"]:
            row = {"output": o["output"]}
            if nominal:
                row["nominal"] = base.get(o["output"], math.nan)
            row.update({k: o[k] for k in ("mean", "sd", "p05", "p50", "p95")})
            for k in ("p05", "p95"):
                row[f"{k}_low"] = o[f"{k}_interval"]["low"]
                row[f"{k}_high"] = o[f"{k}_interval"]["high"]
            row.update({k: o[k] for k in ("min", "max", "valid", "missing", "converged", "r2")})
            rows.append(row)
        return pd.DataFrame(rows)

    def shares(self):
        """Each input's rank correlation, standardised rank regression
        coefficient and share of each output's spread (pandas DataFrame). The
        shares add up to the output's ``r2``; a low ``r2`` means interactions."""
        import pandas as pd

        rows = []
        for o in self.result["outputs"]:
            for s in o["shares"]:
                rows.append({"output": o["output"], "input": s["input"], "share": s["share"],
                             "srrc": s["srrc"], "spearman": s["spearman"], "r2": o["r2"]})
        return pd.DataFrame(rows)

    def error_bars(self) -> Dict[str, tuple]:
        """``{output: (p05, p50, p95)}``."""
        return {o["output"]: (o["p05"], o["p50"], o["p95"]) for o in self.result["outputs"]}

    def to_dict(self) -> Dict[str, Any]:
        return {"method": "monte_carlo", "study": self.study.describe(), "inputs": [i.path for i in self.inputs],
                "batches": self.batches, "rel_tol": self.rel_tol, "values": self.values, "samples": self.samples,
                "result": self.result, "warnings": self.warnings}

    def __repr__(self) -> str:
        state = "converged" if self.converged else "NOT converged"
        return f"<Propagation {len(self.inputs)} inputs, {self.runs} runs ({state}), {self.failed} failed>"
