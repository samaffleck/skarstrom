"""Calibration: fit deck parameters to measured data stored in the decks.

The engine's calibration object (``adsim/solver/calibration.h``) from Python.
One deck and the datasets in it::

    import skarstrom

    sim = skarstrom.examples.Illustrative_Diluted_Bed.load()
    result = sim.calibrate(["adsorbent[1].isotherm[CarbonDioxide].m2"], max_simulations=20)
    print(result.table())          # estimates, standard errors, fit per dataset
    result.apply(sim)              # write the estimates into the deck

Several experiments (decks) with shared and per-case parameters::

    from skarstrom import calibration as cal

    result = skarstrom.calibrate(
        [cal.Case(sim_50, ["outlet_CO2"]), cal.Case(sim_100, ["outlet_CO2"])],
        [cal.Parameter("adsorbent[1].isotherm[CarbonDioxide].m2", 3, 8),                # one value
         cal.Parameter("layer[1].film_coefficient", 1, 500, scale="log", shared=False)])  # one per case

Each dataset must be bound to the quantity it measures
(``sim.experimental_data[name].observable``, :class:`skarstrom.Observable`).
The objective is a weighted least-squares sum over every measured point,
scaled by each point's uncertainty (else by the dataset's y range) and by the
dataset's weight over its number of points; see ``docs/parameterisation.md``.

The runs happen in the engine's threads (``threads``); the GIL is released
while they run. Ctrl-C or :meth:`CancelToken.cancel` stops the fit: Ctrl-C
raises ``KeyboardInterrupt`` afterwards, ``cancel()`` returns the best point
so far (``result.cancelled``). The decks passed are never changed by the fit
(:meth:`CalibrationResult.apply` writes the estimates).
"""

from __future__ import annotations

import json
import math
from dataclasses import dataclass, field
from typing import TYPE_CHECKING, Any, Callable, Dict, List, Optional, Sequence, Tuple, Union

import numpy as np

from . import _core, units

if TYPE_CHECKING:  # pragma: no cover
    from . import Simulation

__all__ = [
    "CalibrationResult",
    "CancelToken",
    "Case",
    "DatasetFit",
    "Estimate",
    "Parameter",
    "calibrate",
]

#: A flag that stops a running calibration from any thread (or from a
#: ``progress`` callback): ``token = CancelToken(); calibrate(..., cancel=token)``,
#: then ``token.cancel()``. The fit returns its best point so far.
CancelToken = _core.CalibrationCancel

METHODS = ("least_squares", "global")
DIFFERENCES = ("forward", "central")
SCALES = ("linear", "log")

DatasetSpec = Union[str, int, Tuple[Union[str, int], Optional[float]]]


@dataclass
class Parameter:
    """One parameter to estimate.

    ``path`` is a tunable registry path or alias (``sim.parameters(tunable_only=True)``).
    ``lower``, ``upper`` and ``start`` are in the parameter's SI unit, or
    carry their own: a pint Quantity, a text such as ``"2 min"`` or
    ``"2 [min]"``, or a number with ``unit="min"``. Bounds left as ``None``
    take the registry's suggested bounds. ``start=None`` starts from the
    deck's value (moved inside the bounds when outside, with a note).

    ``scale="log"`` searches the logarithm (for positive values that span
    decades: rate coefficients, affinities); ``shared=False`` estimates one
    value per case. For a per-case parameter, ``case_bounds`` (a ``(lower,
    upper)`` per case, ``None`` entries keep the bounds above) and
    ``case_start`` (a start per case) are optional.
    """

    path: str
    lower: Any = None
    upper: Any = None
    scale: str = "linear"
    shared: bool = True
    start: Any = None
    unit: Optional[str] = None
    case_bounds: Optional[Sequence[Optional[Tuple[Any, Any]]]] = None
    case_start: Optional[Sequence[Any]] = None

    def __post_init__(self) -> None:
        if self.scale not in SCALES:
            raise ValueError(f"scale is 'linear' or 'log'; got {self.scale!r}")


class Case:
    """One experiment: a deck and the datasets in it to reproduce.

    ``datasets``: names (or 0-based indices) of the deck's experimental data,
    each optionally as ``(name, weight)``; ``None`` fits every bound dataset.
    ``weight`` is the weight of the listed datasets that give none (default:
    each dataset's own ``weight``). ``cycles`` caps the cycles of every run
    (default: the deck's), ``require_css`` counts a run only when it reaches
    cyclic steady state (off by default: calibration data are often
    breakthrough runs), and ``settings`` (path -> value) are applied to every
    run before the parameters (e.g. a tighter ``css.state_tolerance``).

    The deck is copied when the fit starts; later changes to ``sim`` do not
    reach a fit that is running.
    """

    def __init__(self, sim: "Simulation", datasets: Optional[Sequence[DatasetSpec]] = None,
                 weight: Optional[float] = None, name: Optional[str] = None, cycles: Optional[int] = None,
                 require_css: bool = False, settings: Optional[Dict[str, Any]] = None) -> None:
        self.sim = sim
        self.datasets = list(datasets) if datasets is not None else None
        self.weight = None if weight is None else float(weight)
        self.name = name
        self.cycles = cycles
        self.require_css = bool(require_css)
        self.settings = dict(settings or {})

    @classmethod
    def of(cls, item: Any) -> "Case":
        """A case from a :class:`Case`, a Simulation, or a ``(Simulation, datasets[, weight])`` tuple."""
        if isinstance(item, Case):
            return item
        if isinstance(item, tuple):
            if not 1 <= len(item) <= 3:
                raise TypeError("a case tuple is (Simulation, [dataset names], weight)")
            return cls(*item)
        if hasattr(item, "_sim") and hasattr(item, "to_json"):
            return cls(item)
        raise TypeError(f"a case is a Simulation, a (Simulation, datasets, weight) tuple or a Case; got {item!r}")

    def _dataset_refs(self) -> Optional[List[Tuple[Union[str, int], Optional[float]]]]:
        names: Optional[List[Any]] = self.datasets
        if names is None:
            if self.weight is None:
                return None  # every bound dataset, with its own weight
            names = [s.name for s in self.sim.experimental_data if s.observable is not None]
        refs = []
        for d in names:
            if isinstance(d, tuple):
                key, w = d
                refs.append((key, None if w is None else float(w)))
            else:
                refs.append((d, self.weight))
        return refs

    def _to_engine(self) -> Dict[str, Any]:
        d: Dict[str, Any] = {"deck_json": self.sim.to_json(), "require_css": self.require_css,
                             "cycles": self.cycles, "settings": {k: _setting_text(v) for k, v in self.settings.items()}}
        if self.name is not None:
            d["name"] = self.name
        refs = self._dataset_refs()
        if refs is not None:
            d["datasets"] = refs
        return d

    def __repr__(self) -> str:
        what = "every bound dataset" if self.datasets is None else ", ".join(map(str, self.datasets))
        label = f"{self.name!r}: " if self.name else ""
        return f"<Case {label}{what}>"


def _setting_text(value: Any) -> str:
    if units.is_quantity(value):
        return f"{float(value.magnitude)!r} [{value.units:~}]"
    if isinstance(value, bool):
        return "1" if value else "0"
    if isinstance(value, float):
        return repr(value)
    return str(value)


def _to_si(sim: "Simulation", path: str, value: Any, unit: Optional[str], what: str) -> Optional[float]:
    """A bound or start in the parameter's SI unit: a number (in ``unit`` when
    given), a pint Quantity, or a text "2 min" / "2 [min]"."""
    if value is None:
        return None
    info = sim._index().find(path)
    if info is None:
        raise KeyError(f"unknown parameter: {path}; list them with sim.parameters(tunable_only=True)")
    si = info.unit if info.unit_checked else None
    if units.is_quantity(value):
        if si is None:
            return float(value.to_base_units().magnitude)
        return units._quantity_to(value, si, f"{path} ({what})")
    if isinstance(value, str):
        split = _core.split_quantity(value.strip())
        if split is None:
            try:
                number, text_unit = float(value), None
            except ValueError:
                raise ValueError(f"{path}: the {what} {value!r} is not a number or a quantity such as '2 min'") from None
        else:
            number, text_unit = split
            text_unit = text_unit or None
        if text_unit is not None and unit is not None:
            raise ValueError(f"{path}: the {what} {value!r} has a unit; do not also give unit={unit!r}")
        value, unit = number, text_unit or unit
    if isinstance(value, bool) or not isinstance(value, (int, float, np.integer, np.floating)):
        raise TypeError(f"{path}: the {what} must be a number, a quantity or a text such as '2 min'; got {value!r}")
    if unit is None or si is None:
        return float(value)
    return units.convert(float(value), unit, si)


def _parameter(item: Any) -> Parameter:
    if isinstance(item, Parameter):
        return item
    if isinstance(item, str):
        return Parameter(item)
    if isinstance(item, tuple) and 1 <= len(item) <= 4:
        return Parameter(*item)
    raise TypeError(f"a parameter is a path, a (path, lower, upper[, scale]) tuple or a Parameter; got {item!r}")


def _parameter_to_engine(p: Parameter, sim: "Simulation") -> Dict[str, Any]:
    lower = _to_si(sim, p.path, p.lower, p.unit, "lower bound")
    upper = _to_si(sim, p.path, p.upper, p.unit, "upper bound")
    if lower is None or upper is None:
        info = sim._index().find(p.path)
        bounds = info.bounds if info is not None else None
        if bounds is None:
            raise ValueError(f"{p.path} has no suggested bounds; give them: Parameter({p.path!r}, lower, upper)")
        lower = bounds[0] if lower is None else lower
        upper = bounds[1] if upper is None else upper
    d: Dict[str, Any] = {"path": p.path, "lower": lower, "upper": upper, "scale": p.scale, "shared": bool(p.shared),
                         "initial": _to_si(sim, p.path, p.start, p.unit, "start")}
    if p.case_bounds is not None:
        d["case_bounds"] = [None if b is None else (_to_si(sim, p.path, b[0], p.unit, "lower bound"),
                                                    _to_si(sim, p.path, b[1], p.unit, "upper bound"))
                            for b in p.case_bounds]
    if p.case_start is not None:
        d["case_initial"] = [_to_si(sim, p.path, v, p.unit, "start") for v in p.case_start]
    return d


# ---- the result -----------------------------------------------------------------

def _none_to_nan(x: Any) -> float:
    return math.nan if x is None else float(x)


def _nan_to_none(x: Any) -> Any:
    if isinstance(x, (float, np.floating)):
        return None if not math.isfinite(float(x)) else float(x)
    if isinstance(x, np.ndarray):
        return [_nan_to_none(v) for v in x.tolist()]
    if isinstance(x, dict):
        return {k: _nan_to_none(v) for k, v in x.items()}
    if isinstance(x, (list, tuple)):
        return [_nan_to_none(v) for v in x]
    return x


def _array(x: Any) -> np.ndarray:
    return np.asarray([math.nan if v is None else v for v in x], dtype=float)


@dataclass
class Estimate:
    """One estimated value: a shared parameter, or a per-case parameter in one case.

    ``value``, ``std_error``, ``start``, ``lower`` and ``upper`` are in the
    parameter's SI ``unit``. ``std_error`` is NaN when the data cannot fix the
    value (see ``CalibrationResult.notes``); at a bound it is only a local
    statement (``at_bound``).
    """

    name: str            # the path, with " [<case name>]" for a per-case value
    path: str
    case: Optional[int]  # None: shared
    value: float
    std_error: float
    start: float
    lower: float
    upper: float
    at_lower: bool
    at_upper: bool
    scale: str
    unit: str = ""

    @property
    def at_bound(self) -> bool:
        return self.at_lower or self.at_upper

    @property
    def relative_error(self) -> float:
        """``std_error / |value|`` (NaN at zero)."""
        return self.std_error / abs(self.value) if self.value else math.nan


@dataclass
class DatasetFit:
    """How the calibrated model reproduces one dataset (in the dataset's units).

    ``x``, ``measured``, ``simulated`` and ``residuals`` (simulated - measured)
    are at the points in the objective; ``weighted_residuals`` are the
    objective's entries. ``series_x`` and ``series_y`` are the whole simulated
    series at the estimate.
    """

    case: int
    case_name: str
    name: str
    weight: float
    points: int
    skipped: int
    rmse: float
    nrmse: float
    max_abs_error: float
    objective: float
    x: np.ndarray = field(repr=False)
    measured: np.ndarray = field(repr=False)
    simulated: np.ndarray = field(repr=False)
    residuals: np.ndarray = field(repr=False)
    weighted_residuals: np.ndarray = field(repr=False)
    series_x: np.ndarray = field(repr=False)
    series_y: np.ndarray = field(repr=False)


_ARRAYS = ("x", "measured", "simulated", "residuals", "weighted_residuals", "series_x", "series_y")
_FLOATS = ("weight", "rmse", "nrmse", "max_abs_error", "objective")


class CalibrationResult:
    """The outcome of :func:`calibrate`.

    ``parameters`` (one :class:`Estimate` per estimated value), ``correlation``
    (NumPy matrix in the same order; NaN when unavailable), ``datasets`` (one
    :class:`DatasetFit` each), ``objective`` (S at the estimate), ``sigma2``
    (S / (N - p)), ``objective_trace``, ``simulations``, ``failed_simulations``,
    ``iterations``, ``wall_time_s``, ``method``, ``stop_reason``, ``message``,
    ``notes``, ``ok`` and ``cancelled``.
    """

    def __init__(self, data: Dict[str, Any], cases: Optional[List[Case]] = None) -> None:
        self._data = data
        self._cases = cases
        self.parameters: List[Estimate] = [
            Estimate(name=p["label"], path=p["path"], case=None if p["case_index"] < 0 else p["case_index"],
                     value=_none_to_nan(p["estimate"]), std_error=_none_to_nan(p["standard_error"]),
                     start=_none_to_nan(p["initial"]), lower=_none_to_nan(p["lower"]),
                     upper=_none_to_nan(p["upper"]), at_lower=bool(p["at_lower"]), at_upper=bool(p["at_upper"]),
                     scale=p["scale"], unit=p.get("unit", ""))
            for p in data["parameters"]]
        n = len(self.parameters)
        corr = data["correlation"]
        self.correlation: np.ndarray = (np.asarray([[_none_to_nan(v) for v in row] for row in corr], dtype=float)
                                        if len(corr) else np.full((n, n), math.nan)).reshape(n, n)
        self.datasets: List[DatasetFit] = []
        for d in data["datasets"]:
            kw = {k: _array(d[k]) for k in _ARRAYS}
            kw.update({k: _none_to_nan(d[k]) for k in _FLOATS})
            self.datasets.append(DatasetFit(case=int(d["case_index"]), case_name=d["case_name"], name=d["name"],
                                            points=int(d["points"]), skipped=int(d["skipped"]), **kw))
        self.trace: List[Dict[str, Any]] = [
            {"solver": t["solver"], "iteration": t["iteration"], "objective": _none_to_nan(t["objective"]),
             "values": _array(t["values"]), "simulations": t["simulations"]} for t in data["trace"]]
        self.objective: float = _none_to_nan(data["objective"])
        self.sigma2: float = _none_to_nan(data["sigma2"])
        self.residual_count: int = int(data["residuals"])
        self.simulations: int = int(data["simulations"])
        self.failed_simulations: int = int(data["failed_simulations"])
        self.iterations: int = int(data["iterations"])
        self.wall_time_s: float = float(data["wall_time_s"])
        self.method: str = data["solver"]
        self.refined: bool = bool(data["refined"])
        self.stop_reason: str = data["stop_reason"]
        self.message: str = data["message"]
        self.notes: List[str] = list(data["notes"])
        self.cancelled: bool = bool(data["cancelled"])
        self.ok: bool = bool(data["ok"])
        self.calibrated_decks: List[str] = list(data["calibrated_decks"])
        self.case_names: List[str] = list(data.get("case_names", []))

    # -- reading -----------------------------------------------------------------------
    @property
    def objective_trace(self) -> np.ndarray:
        """The objective S after each iteration (or global candidate)."""
        return np.asarray([t["objective"] for t in self.trace], dtype=float)

    @property
    def values(self) -> Dict[str, float]:
        """The estimates by name (SI units)."""
        return {p.name: p.value for p in self.parameters}

    def __getitem__(self, name: str) -> Estimate:
        """An estimate by name (``path`` or ``path [case name]``), or by path for a shared parameter."""
        for p in self.parameters:
            if p.name == name:
                return p
        matches = [p for p in self.parameters if p.path == name]
        if len(matches) == 1:
            return matches[0]
        if matches:
            raise KeyError(f"{name} is per case; ask for one of {[p.name for p in matches]}")
        raise KeyError(f"no estimate {name!r}; have {[p.name for p in self.parameters]}")

    def dataset(self, name: str, case: Optional[int] = None) -> DatasetFit:
        """The fit of one dataset by name (``case``: its case index when several cases have it)."""
        matches = [d for d in self.datasets if d.name == name and (case is None or d.case == case)]
        if len(matches) == 1:
            return matches[0]
        if not matches:
            raise KeyError(f"no dataset {name!r}; have {[d.name for d in self.datasets]}")
        raise KeyError(f"{len(matches)} cases have dataset {name!r}; give case=")

    # -- using the estimates ----------------------------------------------------------
    def _case_values(self, case: int) -> Dict[str, float]:
        out: Dict[str, float] = {}
        for p in self.parameters:
            if p.case is None or p.case == case:
                out[p.path] = p.value
        return out

    def apply(self, target: Union[int, "Simulation", None] = None, case: Optional[int] = None) -> "Simulation":
        """Write the estimates into a deck, as ordinary parameter edits.

        ``apply()`` or ``apply(i)``: case ``i``'s own Simulation (the one
        passed to :func:`calibrate`; default 0). ``apply(sim)``: any deck,
        with the shared values and the per-case values of ``case`` (default:
        the case whose Simulation ``sim`` is, else 0). Returns the deck.
        """
        if target is None or isinstance(target, (int, np.integer)):
            index = int(target or 0)
            if self._cases is None:
                raise ValueError("this result has no decks (it was read from JSON): pass the Simulation, "
                                 "apply(sim, case=i), or use calibrated()")
            if not 0 <= index < len(self._cases):
                raise IndexError(f"the fit has {len(self._cases)} cases")
            sim = self._cases[index].sim
            case = index if case is None else case
        else:
            sim = target
            if case is None:
                case = next((i for i, c in enumerate(self._cases or []) if c.sim is sim), 0)
        if not self.ok:
            raise ValueError(f"the fit has no estimate to apply: {self.message}")
        for path, value in self._case_values(case).items():
            sim.set(path, value)
        return sim

    def calibrated(self, case: Optional[int] = None) -> Union["Simulation", List["Simulation"]]:
        """New Simulations with the estimates applied: case ``case``'s, or, by
        default, the only case's (a list when the fit had several cases). The
        case's ``settings`` are not applied."""
        from . import Simulation

        if not self.calibrated_decks:
            raise ValueError(f"the fit has no calibrated decks: {self.message}")
        if case is not None:
            return Simulation.from_json(self.calibrated_decks[case])
        sims = [Simulation.from_json(text) for text in self.calibrated_decks]
        return sims[0] if len(sims) == 1 else sims

    # -- text ------------------------------------------------------------------------
    def table(self) -> str:
        """The estimates, the fit of each dataset and how the fit stopped, as text."""
        from .cli import _fmt

        def column(rows: List[List[str]], header: List[str]) -> List[str]:
            widths = [max(len(h), *(len(r[i]) for r in rows)) if rows else len(h) for i, h in enumerate(header)]
            lines = ["  ".join(h.ljust(w) for h, w in zip(header, widths)).rstrip(),
                     "  ".join("-" * w for w in widths).rstrip()]
            lines += ["  ".join(c.ljust(w) for c, w in zip(r, widths)).rstrip() for r in rows]
            return lines

        rows = []
        for p in self.parameters:
            rel = f"{100 * p.relative_error:.3g} %" if math.isfinite(p.relative_error) else "n/a"
            flag = "at lower bound" if p.at_lower else "at upper bound" if p.at_upper else ""
            rows.append([p.name, _fmt(p.value), _fmt(p.std_error), rel, p.unit, _fmt(p.start),
                         f"{_fmt(p.lower)} .. {_fmt(p.upper)}", p.scale, flag])
        lines = column(rows, ["parameter", "estimate", "std error", "rel.", "unit", "start", "bounds", "scale", ""])
        lines.append("")
        rows = [[d.name, d.case_name, str(d.points) + (f" (+{d.skipped} skipped)" if d.skipped else ""),
                 _fmt(d.rmse), _fmt(d.nrmse), _fmt(d.max_abs_error), _fmt(d.weight)] for d in self.datasets]
        lines += column(rows, ["dataset", "case", "points", "rmse", "nrmse", "max |error|", "weight"])
        lines.append("")
        n = len(self.parameters)
        if n > 1 and np.isfinite(self.correlation).all():
            off = np.abs(self.correlation - np.eye(n))
            i, j = np.unravel_index(int(np.argmax(off)), off.shape)
            lines.append(f"largest correlation: {self.correlation[i, j]:+.3f} "
                         f"({self.parameters[i].name}, {self.parameters[j].name})")
        failed = f" ({self.failed_simulations} failed)" if self.failed_simulations else ""
        lines.append(f"objective S = {_fmt(self.objective)}, sigma^2 = {_fmt(self.sigma2)}, "
                     f"{self.residual_count} residuals")
        lines.append(f"{self.method}{' + least_squares' if self.refined else ''}: {self.stop_reason} after "
                     f"{self.iterations} iterations, {self.simulations} simulations{failed}, "
                     f"{self.wall_time_s:.1f} s")
        lines.append(self.message)
        lines += [f"note: {n}" for n in self.notes]
        return "\n".join(lines)

    def __str__(self) -> str:
        return self.table()

    def __repr__(self) -> str:
        values = ", ".join(f"{p.name}={p.value:.6g}" for p in self.parameters)
        return f"<CalibrationResult {self.stop_reason}, S={self.objective:.6g}, {values}>"

    # -- serialisation -------------------------------------------------------------
    def to_dict(self) -> Dict[str, Any]:
        """Every field as JSON-ready values (NaN as ``None``, arrays as lists),
        including the calibrated decks. :meth:`from_dict` reads it back."""
        return _nan_to_none({
            "parameters": [{"label": p.name, "path": p.path, "case_index": -1 if p.case is None else p.case,
                            "estimate": p.value, "standard_error": p.std_error, "initial": p.start,
                            "lower": p.lower, "upper": p.upper, "at_lower": p.at_lower, "at_upper": p.at_upper,
                            "scale": p.scale, "unit": p.unit} for p in self.parameters],
            "correlation": self.correlation,
            "datasets": [{"case_index": d.case, "case_name": d.case_name, "name": d.name, "points": d.points,
                          "skipped": d.skipped, **{k: getattr(d, k) for k in _FLOATS},
                          **{k: getattr(d, k) for k in _ARRAYS}} for d in self.datasets],
            "trace": [dict(t) for t in self.trace],
            "objective": self.objective, "sigma2": self.sigma2, "residuals": self.residual_count,
            "simulations": self.simulations, "failed_simulations": self.failed_simulations,
            "iterations": self.iterations, "wall_time_s": self.wall_time_s, "solver": self.method,
            "refined": self.refined, "stop_reason": self.stop_reason, "message": self.message,
            "notes": self.notes, "cancelled": self.cancelled, "ok": self.ok,
            "case_names": self.case_names, "case_settings": self._data.get("case_settings", []),
            "calibrated_decks": self.calibrated_decks,
        })

    def to_json(self, path: Optional[str] = None, indent: Optional[int] = 1) -> str:
        """:meth:`to_dict` as JSON text, also written to ``path`` when given."""
        text = json.dumps(self.to_dict(), indent=indent)
        if path is not None:
            with open(path, "w", encoding="utf-8") as f:
                f.write(text)
        return text

    @classmethod
    def from_dict(cls, data: Dict[str, Any]) -> "CalibrationResult":
        return cls(data)

    @classmethod
    def from_json(cls, text: str) -> "CalibrationResult":
        """A result from :meth:`to_json` text (or a path to such a file)."""
        if not text.lstrip().startswith("{"):
            with open(text, encoding="utf-8") as f:
                text = f.read()
        return cls(json.loads(text))

    # -- plot ------------------------------------------------------------------------
    def _start_series(self, d: DatasetFit) -> Optional[Tuple[np.ndarray, np.ndarray]]:
        """The dataset's simulated counterpart at the starting values (one run of its case)."""
        from . import Simulation

        if not self.calibrated_decks:
            return None
        sim = Simulation.from_json(self.calibrated_decks[d.case])
        settings = self._data.get("case_settings") or []
        if d.case < len(settings):
            for path, text in (settings[d.case].get("settings") or {}).items():
                sim.set(path, text)
            if settings[d.case].get("cycles"):
                sim.cycle.number_of_cycles = settings[d.case]["cycles"]
        for p in self.parameters:
            if p.case is None or p.case == d.case:
                sim.set(p.path, p.start)
        if not sim.run().ok:
            return None
        x, y = sim.experimental_data[d.name].simulated(sim)
        return np.asarray(x), np.asarray(y)

    def plot(self, dataset: Union[str, int, None] = None, start: bool = True, axes: Any = None):
        """Plot the measured points, the simulated curve at the estimate and
        (``start=True``: one more run per dataset's case) at the start, with
        the residuals below. ``dataset`` picks one (name or index; default:
        every dataset, one column each). Needs matplotlib. Returns the figure."""
        import matplotlib.pyplot as plt

        if dataset is None:
            chosen = self.datasets
        elif isinstance(dataset, int):
            chosen = [self.datasets[dataset]]
        else:
            chosen = [d for d in self.datasets if d.name == dataset]
            if not chosen:
                raise KeyError(f"no dataset {dataset!r}; have {[d.name for d in self.datasets]}")
        if axes is None:
            fig, axes = plt.subplots(2, len(chosen), squeeze=False, sharex="col",
                                     gridspec_kw={"height_ratios": [3, 1]}, figsize=(5 * len(chosen), 5))
        else:
            axes = np.asarray(axes).reshape(2, len(chosen))
            fig = axes[0, 0].figure
        for k, d in enumerate(chosen):
            top, bottom = axes[0, k], axes[1, k]
            top.plot(d.x, d.measured, "o", ms=3, label="measured")
            if start:
                s = self._start_series(d)
                if s is not None:
                    top.plot(s[0], s[1], "--", label="start")
            top.plot(d.series_x, d.series_y, "-", label="calibrated")
            top.set_title(d.name if len(self.case_names) <= 1 else f"{d.name} [{d.case_name}]", fontsize="small")
            top.legend(fontsize="small")
            bottom.axhline(0.0, color="0.6", lw=0.8)
            bottom.plot(d.x, d.residuals, "o", ms=3)
            bottom.set_ylabel("residual")
        fig.tight_layout()
        return fig


# ---- the fit ----------------------------------------------------------------------

def calibrate(
    cases: Union[Any, Sequence[Any]],
    parameters: Sequence[Union[str, Parameter, tuple]],
    method: str = "least_squares",
    refine: bool = False,
    threads: Optional[int] = None,
    max_simulations: int = 500,
    max_iterations: int = 50,
    ftol: float = 1e-6,
    xtol: float = 1e-6,
    differences: str = "forward",
    fd_step: float = 1e-3,
    global_evaluations: int = 30,
    warm_start: bool = True,
    standard_errors: bool = True,
    progress: Optional[Callable[[Dict[str, Any]], Any]] = None,
    cancel: Optional[CancelToken] = None,
) -> CalibrationResult:
    """Fit deck parameters to the measured datasets stored in the decks.

    ``cases``: one experiment or a list of them; each a Simulation (every
    bound dataset), a ``(Simulation, [dataset names], weight)`` tuple or a
    :class:`Case`. ``parameters``: registry paths (bounds from the registry's
    suggested bounds), ``(path, lower, upper[, scale])`` tuples or
    :class:`Parameter` objects (shared or per case, linear or log scale).

    ``method="least_squares"`` (bounded Levenberg-Marquardt, the default) is
    fast and gives standard errors; ``"global"`` (MaxLIPO over the bounds,
    ``global_evaluations`` candidates) is for a poor start or several minima,
    and ``refine=True`` polishes its best point with least squares.
    ``differences="central"`` doubles the cost of a Jacobian but is needed
    for strongly correlated parameters. ``max_simulations`` caps the runs
    (one case run is one simulation); ``ftol`` / ``xtol`` stop least squares
    on a small relative change of the objective / of the normalised values.
    ``threads`` runs that many simulations at once (default: the cores - 1).
    ``warm_start`` starts a cyclic deck's runs from the nearest cached
    steady state.

    ``progress(info)`` is called after each iteration with ``solver``,
    ``iteration``, ``objective``, ``values`` (per estimate, SI) and
    ``simulations``; returning ``False`` stops the fit. ``cancel`` is a
    :class:`CancelToken`. Raises ``ValueError`` for an invalid problem (an
    unknown path, a dataset that is missing or not bound to an observable,
    bad bounds), before any run.
    """
    if method not in METHODS:
        raise ValueError(f"method is 'least_squares' or 'global'; got {method!r}")
    if differences not in DIFFERENCES:
        raise ValueError(f"differences is 'forward' or 'central'; got {differences!r}")
    if isinstance(cases, (Case, tuple)) or hasattr(cases, "to_json"):
        cases = [cases]
    case_list = [Case.of(c) for c in cases]
    if not case_list:
        raise ValueError("a calibration needs at least one case")
    if isinstance(parameters, (str, Parameter)):
        parameters = [parameters]
    params = [_parameter(p) for p in parameters]
    if not params:
        raise ValueError("a calibration needs at least one parameter")
    first = case_list[0].sim
    engine_cases = [c._to_engine() for c in case_list]
    engine_params = [_parameter_to_engine(p, first) for p in params]
    options = {"method": method, "refine": bool(refine), "threads": 0 if threads is None else int(threads),
               "max_simulations": int(max_simulations), "max_iterations": int(max_iterations),
               "ftol": float(ftol), "xtol": float(xtol), "differences": differences, "fd_step": float(fd_step),
               "global_evaluations": int(global_evaluations), "warm_start": bool(warm_start),
               "standard_errors": bool(standard_errors)}
    raw = _core.calibrate(engine_cases, engine_params, options, progress, cancel)
    unit_of = {}
    for p in params:
        info = first._index().find(p.path)
        unit_of[p.path] = info.unit if info is not None and info.unit not in ("-",) else ""
    for p in raw["parameters"]:
        p["unit"] = unit_of.get(p["path"], "")
    names = {d["case_index"]: d["case_name"] for d in raw["datasets"]}
    raw["case_names"] = [names.get(i, f"case {i + 1}") for i in range(len(case_list))]
    raw["case_settings"] = [{"settings": c["settings"], "cycles": c["cycles"]} for c in engine_cases]
    return CalibrationResult(raw, case_list)
