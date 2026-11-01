"""Skarstrom: pressure and temperature swing adsorption cycle simulation.

Load a deck (a JSON simulation file; :mod:`skarstrom.examples` holds the
bundled ones), change it, run it to cyclic steady state and read the results::

    import skarstrom
    from skarstrom import ureg as u                # pint's unit registry

    sim = skarstrom.examples.Perez2019_A.load()   # or skarstrom.Simulation("deck.json")
    sim.cycle.number_of_cycles = 20
    sim.steps["Feed"].time = 2 * u.min             # a step by name; a pint quantity...
    sim.steps["Feed"].set("time", 2, "min")        # ...or a number and its unit...
    sim.steps["Feed"].bottom.temperature = 298.15  # ...or a bare number, in SI (K, s, Pa, m, mol/s)
    sim.layers[0].length = 1.2 * u.m
    result = sim.run()
    print(result.status, result.css_cycle)
    for k in sim.kpis(cycle=-1):
        print(k.component, k.purity, k.recovery)
    y = sim.series("mole_fraction", "CarbonDioxide", location="top")
    y.time, y.values                               # NumPy arrays

Every value also has a stable string path (``sim["step[Feed].time"] = 120``)
for sweeps, the optimiser and the CLI; ``sim.parameters()`` lists them. A bare
number is always in SI units; see :mod:`skarstrom.units`. Full guide:
https://samaffleck.github.io/skarstrom/python/
"""

from __future__ import annotations

import functools
import json
import math
import os
from dataclasses import dataclass
from typing import Any, Callable, Dict, List, Mapping, Optional, Sequence, Union

import numpy as np

from . import _core, examples, sensitivity, units
from .units import Quantity, ureg
from ._core import (
    AxialProfile,
    CSSRecord,
    KPIRecord,
    OperatingLimitsReport,
    OperatingWarning,
    RunResult,
    StepOperatingSummary,
    TimeSeries,
    WallDesignResult,
)
from ._example import Example
from .build import (
    Flow,
    KPIConfig,
    PolynomialPower,
    Profile,
    Pump,
    PumpStream,
    Stream,
    ValveStream,
    build_kinds,
    constant,
    exponential,
    linear,
    log_linear,
    double_exponential,
    mass_flow,
    molar_flow,
    polynomial_power,
    to_script,
    velocity,
    volumetric_flow,
)
from .model import (
    CSS,
    Adsorbent,
    Bed,
    BufferTank,
    Collection,
    Cycle,
    EnergyBalance,
    EnergySource,
    Endpoint,
    Fluid,
    Info,
    InitialConditions,
    Isotherm,
    Kinetics,
    Layer,
    NodeMap,
    Settings,
    Step,
    StepEnd,
    Valve,
    ValueMap,
    VoidRegion,
    Wall,
)
from . import build
from . import model as _model

__version__: str = _core.__version__

__all__ = [
    "Adsorbent",
    "AxialProfile",
    "Bed",
    "BufferTank",
    "CSS",
    "CSSRecord",
    "OperatingLimitsReport",
    "OperatingWarning",
    "StepOperatingSummary",
    "WallDesignResult",
    "Collection",
    "Cycle",
    "EnergyBalance",
    "EnergySource",
    "ExperimentalData",
    "ExperimentalSeries",
    "Endpoint",
    "Example",
    "Flow",
    "Fluid",
    "Info",
    "InitialConditions",
    "Isotherm",
    "KPIConfig",
    "KPIRecord",
    "Kinetics",
    "Layer",
    "NodeMap",
    "Parameter",
    "PolynomialPower",
    "Profile",
    "Pump",
    "PumpStream",
    "Quantity",
    "RunResult",
    "Settings",
    "Simulation",
    "Step",
    "StepEnd",
    "Stream",
    "TimeSeries",
    "Valve",
    "ValueMap",
    "ValveStream",
    "VoidRegion",
    "Wall",
    "__version__",
    "adsorbent_library",
    "adsorbent_library_dir",
    "adsorbent_path",
    "build",
    "build_kinds",
    "bulk_variables",
    "component_variables",
    "constant",
    "evaluate_iast",
    "evaluate_isotherm",
    "examples",
    "exponential",
    "fit_isotherm",
    "isotherm_kinds",
    "isotherm_parameters",
    "linear",
    "log_linear",
    "double_exponential",
    "mass_flow",
    "molar_flow",
    "polynomial_power",
    "stream_component_variables",
    "stream_variables",
    "to_script",
    "units",
    "ureg",
    "velocity",
    "volumetric_flow",
]

PathLike = Union[str, "os.PathLike[str]"]


@dataclass(frozen=True)
class Parameter:
    """One settable deck value, as listed by :meth:`Simulation.parameters`.

    ``path`` is the stable, ID-based address (it survives renames); ``alias``
    is the same address written with display names. Either can be passed to
    :meth:`Simulation.get` / :meth:`Simulation.set`.
    """

    path: str
    alias: str
    category: str
    subcategory: str
    name: str
    unit: str  # SI unit ("-" dimensionless); bare numbers are in this unit
    kind: str  # "real", "expression", "integer", "boolean" or "string"
    tunable: bool
    lower: float
    upper: float
    value: Union[float, bool, str]
    text: str
    unit_checked: bool = True  # False: a custom-unit user variable (units are not checked)
    choices: Sequence[str] = ()  # names of a model/mode switch's values (assign a name or its number)

    @property
    def bounds(self) -> Optional[tuple]:
        """Suggested (lower, upper) search bounds, or None."""
        if math.isnan(self.lower) or math.isnan(self.upper):
            return None
        return (self.lower, self.upper)


def _data_dir() -> str:
    return os.path.join(os.path.dirname(__file__), "data")


def adsorbent_library_dir() -> str:
    """Directory of the curated adsorbent library shipped with the package.

    It holds ``library.json`` (the catalogue: citation, DOI and per-parameter
    provenance of every entry) and one adsorbent file per entry
    (``<key>.json``, a single serialised adsorbent).
    """
    return os.path.join(_data_dir(), "adsorbents")


def adsorbent_library() -> Dict[str, Any]:
    """The curated adsorbent library catalogue (``library.json``) as a dict.

    ``adsorbent_library()["entries"]`` lists every adsorbent with its ``key``,
    ``label``, components, sources and properties; ``["sources"]`` maps source
    ids to citations and DOIs.
    """
    with open(os.path.join(adsorbent_library_dir(), "library.json"), encoding="utf-8") as f:
        return json.load(f)


def adsorbent_path(key: str) -> str:
    """Absolute path of one library adsorbent file, by its library key."""
    for e in adsorbent_library()["entries"]:
        if e["key"] == key:
            return os.path.join(adsorbent_library_dir(), e["file"])
    raise KeyError(f"no adsorbent {key!r} in the library")


_NUMERIC_KINDS = ("real", "expression")


@functools.lru_cache(maxsize=256)
def _dimensional(unit: str) -> bool:
    """Whether a parameter unit has a physical dimension (so a bare number is ambiguous)."""
    if unit in ("", "-"):
        return False
    try:
        return not units.parse(unit).dimensionless
    except ValueError:
        return True


def _is_number(value: Any) -> bool:
    return isinstance(value, (int, float, np.integer, np.floating)) and not isinstance(value, bool)


@dataclass(frozen=True)
class ExperimentalSeries:
    """One set of measured (x, y) points stored in the deck."""

    name: str
    x: np.ndarray
    y: np.ndarray
    x_label: str = ""
    y_label: str = ""

    def compare(self, x: Any, y: Any, normalise: Optional[float] = None) -> Dict[str, float]:
        """Compare a simulated curve (``x``, ``y``, e.g. a :class:`TimeSeries`'s ``time`` and
        ``values``) with these points.

        The simulation is interpolated linearly at the measured x values that it
        covers. Returns ``rmse``, ``nrmse`` (the RMSE over ``normalise``, by
        default the range of the measured y), ``max_abs_error`` and ``points``
        (how many measured points were compared).
        """
        xs = np.asarray(x, dtype=float)
        ys = np.asarray(y, dtype=float)
        inside = (self.x >= xs.min()) & (self.x <= xs.max())
        if not inside.any():
            raise ValueError(f"the simulated x range [{xs.min()}, {xs.max()}] does not cover '{self.name}'")
        err = np.interp(self.x[inside], xs, ys) - self.y[inside]
        rmse = float(np.sqrt(np.mean(err**2)))
        scale = float(normalise) if normalise is not None else float(np.ptp(self.y))
        return {"rmse": rmse, "nrmse": rmse / scale if scale > 0 else math.nan,
                "max_abs_error": float(np.max(np.abs(err))), "points": int(inside.sum())}


class ExperimentalData:
    """Measured data stored in the deck (``sim.experimental_data``): breakthrough
    curves, temperature histories, KPIs per cycle... It is saved with the deck,
    so a case and the data it is checked against travel together.

    Index by name or position; ``add``, ``load_csv``, ``remove`` and ``clear``
    change it.
    """

    def __init__(self, sim: "Simulation") -> None:
        self._sim = sim

    def _all(self) -> List[ExperimentalSeries]:
        return [ExperimentalSeries(d["name"], np.asarray(d["x"]), np.asarray(d["y"]), d["x_label"], d["y_label"])
                for d in self._sim._sim.experimental_data()]

    def add(self, name: str, x: Sequence[float], y: Sequence[float], x_label: str = "",
            y_label: str = "") -> ExperimentalSeries:
        """Store a set of points; ``x_label`` and ``y_label`` describe the columns (e.g. ``"time [s]"``)."""
        self._sim._sim.add_experimental_data(name, [float(v) for v in x], [float(v) for v in y], x_label, y_label)
        return self[-1]

    def load_csv(self, path: PathLike, name: Optional[str] = None, x_column: Union[int, str] = 0,
                 y_column: Union[int, str] = 1, delimiter: Optional[str] = None) -> ExperimentalSeries:
        """Read two columns of a delimited text file (comma, semicolon, tab or
        whitespace; a header row is detected and gives the labels).

        Columns are chosen by position or by header name. ``name`` defaults to
        the file name without its extension.
        """
        import csv

        text = open(os.fspath(path), newline="").read()
        if delimiter is None:
            try:
                delimiter = csv.Sniffer().sniff(text.splitlines()[0], delimiters=",;\t ").delimiter
            except (csv.Error, IndexError):
                delimiter = ","
        rows = [r for r in csv.reader(text.splitlines(), delimiter=delimiter, skipinitialspace=True) if r]

        def number(cell: str) -> Optional[float]:
            try:
                return float(cell)
            except ValueError:
                return None

        header = rows[0] if rows and any(number(c) is None for c in rows[0]) else None
        body = rows[1:] if header else rows

        def index(column: Union[int, str]) -> int:
            if isinstance(column, int):
                return column
            if header is None or column not in header:
                raise KeyError(f"no column {column!r} in {os.fspath(path)}")
            return header.index(column)

        ix, iy = index(x_column), index(y_column)
        pts = [(number(r[ix]), number(r[iy])) for r in body if len(r) > max(ix, iy)]
        pts = [(a, b) for a, b in pts if a is not None and b is not None]
        if not pts:
            raise ValueError(f"no numeric rows in {os.fspath(path)}")
        label = (lambda i: header[i] if header else "")
        return self.add(name or os.path.splitext(os.path.basename(os.fspath(path)))[0],
                        [a for a, _ in pts], [b for _, b in pts], label(ix), label(iy))

    def remove(self, key: Union[int, str]) -> None:
        """Delete one set (by name or position)."""
        self._sim._sim.remove_experimental_data(self._index(key))

    def clear(self) -> None:
        for _ in range(len(self)):
            self._sim._sim.remove_experimental_data(0)

    def names(self) -> List[str]:
        return [s.name for s in self._all()]

    def _index(self, key: Union[int, str]) -> int:
        items = self._all()
        if isinstance(key, int):
            if not -len(items) <= key < len(items):
                raise IndexError(f"experimental data has {len(items)} sets")
            return key % len(items)
        for i, s in enumerate(items):
            if s.name == key:
                return i
        raise KeyError(f"no experimental data named {key!r}; have {[s.name for s in items]}")

    def __getitem__(self, key: Union[int, str]) -> ExperimentalSeries:
        return self._all()[self._index(key)]

    def __iter__(self):
        return iter(self._all())

    def __len__(self) -> int:
        return len(self._sim._sim.experimental_data())

    def __contains__(self, name: object) -> bool:
        return name in self.names()

    def __repr__(self) -> str:
        rows = [f"  {s.name!r}: {len(s.x)} points ({s.x_label or 'x'} vs {s.y_label or 'y'})" for s in self._all()]
        return "<experimental data" + (":\n" + "\n".join(rows) if rows else ": none") + ">"


class Simulation:
    """One deck and its solver.

    ``Simulation(deck)`` accepts a path to a deck file, an :class:`Example`
    (``skarstrom.examples.Perez2019_A``), the deck's JSON text, or nothing (an
    empty deck, mainly useful with :meth:`load`).

    The deck's parts are objects: :attr:`steps`, :attr:`layers`,
    :attr:`adsorbents`, :attr:`variables`, :attr:`cycle`, :attr:`bed`,
    :attr:`wall`, :attr:`settings`, :attr:`css`, :attr:`energy_balance` and
    :attr:`info` (see :mod:`skarstrom.model`). Every value also has a string
    path for :meth:`get` / :meth:`set` and ``sim[path]``.

    With ``require_units=True`` a bare number for a dimensional parameter is
    rejected (``TypeError``): every time, temperature, pressure... must carry
    its unit (``2 * ureg.minute``, ``set(path, 2, unit="min")``).
    """

    def __init__(self, deck: Optional[Union[PathLike, str, Example]] = None, require_units: bool = False) -> None:
        self._sim = _core.Simulation()
        self._index_cache: Optional[_model._Index] = None
        self._epoch = 0  # bumped by every structural change; deck objects re-read themselves
        self.require_units = require_units
        if deck is not None:
            text = os.fspath(deck)
            if text.lstrip().startswith("{"):
                self._sim.load_json(text)
            else:
                self._sim.load(text)

    # ---- deck --------------------------------------------------------------
    @classmethod
    def from_file(cls, path: PathLike, require_units: bool = False) -> "Simulation":
        sim = cls(require_units=require_units)
        sim._sim.load(os.fspath(path))
        return sim

    @classmethod
    def from_json(cls, text: str, require_units: bool = False) -> "Simulation":
        sim = cls(require_units=require_units)
        sim._sim.load_json(text)
        return sim

    def load(self, path: PathLike) -> None:
        self._sim.load(os.fspath(path))
        self._invalidate()

    def save(self, path: PathLike) -> None:
        """Write the deck (with any parameter changes) as a deck file."""
        self._sim.save(os.fspath(path))

    def to_json(self) -> str:
        return self._sim.to_json()

    def copy(self) -> "Simulation":
        """An independent copy of the deck (results are not copied)."""
        return Simulation.from_json(self.to_json(), require_units=self.require_units)

    def check(self) -> tuple:
        """(ok, message): evaluates every expression and runs the deck checks."""
        return self._sim.check()

    # ---- the deck as objects ---------------------------------------------------
    def _index(self) -> "_model._Index":
        if self._index_cache is None:
            self._index_cache = _model._Index(self.parameters(), self._sim.steps, self._sim.streams())
        return self._index_cache

    def _invalidate(self) -> None:
        self._index_cache = None
        self._epoch += 1

    @property
    def steps(self) -> Collection[Step]:
        """The cycle's steps in order: ``sim.steps["Feed"]``, ``sim.steps[0]``."""
        return _model._steps(self)

    @property
    def layers(self) -> Collection[Layer]:
        """The bed's adsorbent layers, feed end first."""
        return _model._layers(self)

    @property
    def adsorbents(self) -> Collection[Adsorbent]:
        """The adsorbents the layers use: ``sim.adsorbents["13X"].isotherm["CO2"].b1``."""
        return _model._adsorbents(self)

    @property
    def variables(self) -> ValueMap:
        """User variables by name (the names expressions refer to)."""
        return _model._variables(self)

    @property
    def cycle(self) -> Cycle:
        """Cycle-level values: number of cycles, total cycle time."""
        return _model._group(self, Cycle, "cycle")

    @property
    def bed(self) -> Bed:
        """Column geometry and ambient temperature."""
        return _model._group(self, Bed, "bed")

    @property
    def wall(self) -> Wall:
        """Column wall properties."""
        return _model._group(self, Wall, "wall")

    @property
    def settings(self) -> Settings:
        """Integrator and logging settings."""
        return _model._group(self, Settings, "settings")

    @property
    def css(self) -> CSS:
        """Cyclic-steady-state detection settings."""
        return _model._group(self, CSS, "css")

    @property
    def energy_balance(self) -> EnergyBalance:
        """Energy-balance switches (isothermal or not)."""
        return _model._group(self, EnergyBalance, "energy_balance")

    @property
    def initial_conditions(self) -> InitialConditions:
        """The state of the bed at the start of the first cycle."""
        return _model._group(self, InitialConditions, "initial")

    @property
    def fluid(self) -> Fluid:
        """Gas and adsorbed-phase properties (adsorbed-phase heat capacity, heat of adsorption convention)."""
        return _model._group(self, Fluid, "fluid")

    @property
    def experimental_data(self) -> ExperimentalData:
        """Measured data stored in the deck (see :class:`ExperimentalData`)."""
        return ExperimentalData(self)

    @property
    def kpi(self) -> KPIConfig:
        """KPI configuration: ``sim.kpi.track("CO2", feed=[...], product=[...])``."""
        return KPIConfig(self)

    @property
    def info(self) -> Info:
        """Title, description and author of the deck."""
        return _model._group(self, Info, "info")

    # ---- parameters by path ------------------------------------------------------
    def parameters(self, tunable_only: bool = False, search: Optional[str] = None) -> List[Parameter]:
        """Every settable parameter of this deck.

        ``tunable_only`` keeps the numeric model values the optimiser may vary
        (each with suggested bounds). ``search`` keeps entries whose path,
        alias or name contains the text (case-insensitive).
        """
        out = [Parameter(**d) for d in self._sim.parameters(tunable_only)]
        if search:
            s = search.lower()
            out = [p for p in out if s in p.path.lower() or s in p.alias.lower() or s in p.name.lower()]
        return out

    def get(self, path: str, unit: Optional[str] = None) -> Union[float, int, str]:
        """A parameter's value: float in its SI unit (or converted to ``unit``),
        int for integers, bool for switches, str for text.

        Raises KeyError for a path the deck does not have.
        """
        try:
            value = self._sim.get(path)
        except IndexError as e:
            raise KeyError(str(e)) from None
        if unit is None:
            return value
        info = self._index().find(path)
        if info is None or info.kind not in _NUMERIC_KINDS or not info.unit_checked:
            raise TypeError(f"{path} has no unit to convert from")
        return units.convert(float(value), info.unit, unit)

    def get_text(self, path: str) -> str:
        """A parameter's text: the expression for expression parameters."""
        try:
            return self._sim.get_text(path)
        except IndexError as e:
            raise KeyError(str(e)) from None

    def set(self, path: str, value: Union[float, int, str, Any], unit: Optional[str] = None) -> None:
        """Set a parameter.

        * a number is in the parameter's SI unit (``unit="min"`` says otherwise);
        * a pint Quantity is converted with pint;
        * a string is an expression of user variables (``"2*feed_time"``), a
          deck-syntax quantity (``"2 [min]"``, also ``"2 min"``), or the text
          of a text parameter.

        A unit of the wrong dimension raises ValueError, as does an invalid
        expression (the old value is kept).
        """
        info = self._index().find(path)
        try:
            self._set(path, info, value, unit)
        except IndexError as e:
            raise KeyError(str(e)) from None
        if info is None or info.kind not in _NUMERIC_KINDS:
            self._invalidate()  # a model switch or a renamed object can change the deck's structure

    def _set(self, path: str, info: Optional[Parameter], value: Any, unit: Optional[str]) -> None:
        if isinstance(value, bool):
            value = int(value)
        if unit is not None:
            if not _is_number(value):
                raise TypeError(f"unit= gives the unit of a plain number; got {type(value).__name__} for {path}")
            self._sim.set_text(path, f"{float(value)!r} [{unit}]")
        elif units.is_quantity(value):
            if info is None:
                raise KeyError(f"unknown parameter: {path}")
            if info.kind not in _NUMERIC_KINDS:
                raise TypeError(f"{path} is {info.kind}, not a dimensional value")
            if info.unit_checked:
                number = units._quantity_to(value, info.unit, path)
            else:
                number = float(value.to_base_units().magnitude)
            self._sim.set_value(path, number)
        elif isinstance(value, str):
            if self.require_units and info is not None and self._needs_unit(info):
                try:
                    float(value)
                except ValueError:
                    pass
                else:
                    raise TypeError(self._unit_required(info, value))
            self._sim.set_text(path, value)
        elif _is_number(value):
            if self.require_units and info is not None and self._needs_unit(info):
                raise TypeError(self._unit_required(info, value))
            self._sim.set_value(path, float(value))
        else:
            raise TypeError(f"cannot set {path} to a {type(value).__name__}")

    @staticmethod
    def _needs_unit(info: Parameter) -> bool:
        return info.kind in _NUMERIC_KINDS and info.unit_checked and _dimensional(info.unit)

    @staticmethod
    def _unit_required(info: Parameter, value: Any) -> str:
        return (f"{info.path} is a {units.dimension_name(info.unit)} [{info.unit}]; with require_units=True "
                f"give the unit, e.g. '{value} {info.unit}' or set(..., unit='{info.unit}')")

    def update(self, values: Mapping[str, Union[float, int, str, Any]]) -> None:
        """Set several parameters: ``sim.update({"step[1].time": 50 * ureg.s, ...})``."""
        for k, v in values.items():
            self.set(k, v)

    def __getitem__(self, path: str) -> Union[float, int, str]:
        return self.get(path)

    def __setitem__(self, path: str, value: Union[float, int, str, Any]) -> None:
        self.set(path, value)

    # ---- run -------------------------------------------------------------------
    def run(
        self,
        cycles: Optional[int] = None,
        progress: Optional[Callable[[Dict[str, Any]], Optional[bool]]] = None,
        log_level: str = "low",
        on_message: Optional[Callable[[str, str], None]] = None,
    ) -> RunResult:
        """Run the deck (blocking) and return the run summary.

        ``cycles`` overrides the deck's number of cycles. ``progress`` is
        called after every step with a dict (cycle, step, total_cycles,
        total_steps, fraction); return False to cancel. ``on_message`` gets
        the solver's messages up to ``log_level`` ("low", "medium", "high").
        Ctrl-C cancels the run and raises KeyboardInterrupt.
        """
        return self._sim.run(cycles, progress, log_level, on_message)

    @property
    def last_result(self) -> RunResult:
        return self._sim.last_result

    @property
    def messages(self) -> List[tuple]:
        """Solver messages of the last run as (level, text) pairs."""
        return self._sim.messages

    # ---- structural edits ---------------------------------------------------------
    # ---- building from scratch (skarstrom.build) --------------------------------------
    @classmethod
    def new(cls, title: str = "", description: str = "", author: str = "",
            require_units: bool = False) -> "Simulation":
        """An empty deck to build in code: add components, adsorbents, layers and steps."""
        sim = cls(require_units=require_units)
        for key, value in (("title", title), ("description", description), ("author", author)):
            if value:
                sim.set(f"info.{key}", value)
        return sim

    def _build(self, op: str, *args: Any) -> Any:
        try:
            result = self._sim.build(op, list(args))
        finally:
            self._invalidate()
        return result

    def _component(self, name: str) -> str:
        return _model._resolve_key(name, self.components, "components")

    def _adsorbent_id(self, adsorbent: Any) -> int:
        if isinstance(adsorbent, Adsorbent):
            return adsorbent.id
        if isinstance(adsorbent, int):
            return adsorbent
        return self.adsorbents[adsorbent].id

    def _layer_id(self, layer: Any) -> int:
        if isinstance(layer, Layer):
            return layer.id
        if isinstance(layer, int):
            return layer
        return self.layers[layer].id

    def _endpoint(self, endpoint_id: int) -> Endpoint:
        return _model._endpoint(self, endpoint_id)

    def describe(self) -> Dict[str, Any]:
        """The deck's structure: the model kind of every object, how each step end is set up,
        connections and KPI endpoints (endpoint IDs). Values are in :meth:`parameters`."""
        return self._sim.describe()

    def add_adsorbent(self, name: str, **values: Any) -> Adsorbent:
        """A new adsorbent (default physics for every component: inert isotherm, LDF with a
        constant coefficient); ``values`` set its fields (solid_density, heat_capacity,
        thermal_conductivity, particle_void_fraction)."""
        adsorbent_id = self._build("add_adsorbent", name)
        ads = self.adsorbents.by_id(adsorbent_id)
        _model._apply(ads, values)
        return ads

    def add_layer(self, name: str, composition: Optional[Mapping[Any, float]] = None, momentum: str = "ergun",
                  **values: Any) -> Layer:
        """A new layer at the product end of the bed.

        ``composition`` maps adsorbents (names, IDs or Adsorbent objects) to
        mass fractions; ``momentum`` is ``"ergun"`` or ``"kozeny_carman"``;
        ``values`` set its fields (length, cells, bed_void_fraction,
        particle_diameter...). The transport models start as constant
        dispersion and constant heat transfer; change them with
        ``layer.set_dispersion(...)`` and the like.
        """
        layer_id = self._build("add_layer", name)
        if composition:
            pairs = [[self._adsorbent_id(a), float(f)] for a, f in composition.items()]
            self._build("set_layer_composition", layer_id, pairs)
        self._build("set_momentum", layer_id, momentum)
        layer = self.layers.by_id(layer_id)
        _model._apply(layer, values)
        return layer

    def add_step(self, name: str, time: Any = None, fraction: Any = None) -> Step:
        """A new step at the end of the cycle, both ends closed until set
        (``step.bottom.set_streams(...)``, ``step.top.set_pressure(...)``)."""
        step_id = self._build("add_step", name)
        step = self.steps.by_id(step_id)
        if time is not None:
            step.time = time
        if fraction is not None:
            step.fraction = fraction
        return step

    def connect(self, a: Any, b: Any) -> None:
        """Connect two endpoints so one supplies the other (light-product pressurisation, a purge
        drawn from another step's product). Each is an Endpoint, a single-endpoint step end, or an ID."""
        ids = build._endpoint_ids(a) + build._endpoint_ids(b)
        if len(ids) != 2:
            raise ValueError("connect() takes exactly two endpoints (choose one stream of a multi-stream end)")
        self._build("connect", ids[0], ids[1])

    def add_variable(self, name: str, value: Any = None, category: str = "dimensionless",
                     formula: Optional[str] = None, description: str = "") -> None:
        """A user variable that other values' expressions can use (``sim.steps["Feed"].time = "2 * t_feed"``).

        ``category`` is its unit category (``skarstrom.build_kinds()["variable_category"]``:
        ``"time"``, ``"pressure"``...), which fixes the SI unit of ``value`` and lets a
        quantity (``30 * ureg.s``) be converted. Give a ``value`` or a ``formula`` of other
        variables. A variable of the same name is replaced.
        """
        self._build("add_variable", name, category, formula or "", description)
        if value is not None:
            self.set(f"variable[{name}]", value)

    def add_component(self, name: str) -> None:
        """Add a fluid component (a CoolProp name such as ``"Argon"``).

        The component gets a zero mole fraction in every inlet, boundary and
        initial condition, and default physics on every adsorbent (inert
        isotherm, LDF rate model, zero heat of adsorption) for you to set.
        """
        self._build("add_component", name)

    def remove_component(self, name: str) -> None:
        """Remove a component from the fluid, every adsorbent and every mole-fraction map."""
        self._sim.remove_component(name)
        self._invalidate()

    def add_library_adsorbent(self, key: str) -> int:
        """Add an adsorbent from the bundled library by its key; returns its adsorbent ID.

        Components of the deck that the entry does not cover get default
        physics. A layer with no adsorbent yet gets it automatically; otherwise
        assign it to a layer in C++ (its parameters appear in
        :meth:`parameters` once a layer uses it). See :func:`adsorbent_library`.
        """
        adsorbent_id = int(self._sim.add_library_adsorbent(adsorbent_library_dir(), key))
        self._invalidate()
        return adsorbent_id

    # ---- structure -------------------------------------------------------------
    @property
    def components(self) -> List[str]:
        """Fluid components (CoolProp names), in deck order."""
        return self._sim.components

    def streams(self) -> List[Dict[str, Any]]:
        """Boundary endpoints: dicts with step (index), step_name, name, endpoint_id and location."""
        return self._sim.streams()

    # ---- results ---------------------------------------------------------------
    @property
    def has_results(self) -> bool:
        return self._sim.has_results

    @property
    def cycles(self) -> List[int]:
        """Logged cycle indices (0-based)."""
        return self._sim.cycles

    @property
    def last_cycle(self) -> int:
        return self._sim.last_cycle

    @property
    def positions(self) -> np.ndarray:
        """Cell-centre positions along the bed [m]."""
        return self._sim.positions

    def kpis(self, cycle: Optional[int] = None) -> List[KPIRecord]:
        """KPIs of every tracked component: all logged cycles, or one (-1 = last)."""
        return self._sim.kpis(cycle)

    def kpi_table(self, cycle: Optional[int] = None):
        """:meth:`kpis` as a pandas DataFrame (needs pandas)."""
        import pandas as pd

        return pd.DataFrame([k.as_dict() for k in self.kpis(cycle)])

    def css_history(self) -> List[CSSRecord]:
        """The deck's cyclic-steady-state check at every logged cycle."""
        return self._sim.css_history()

    def operating_limits(self, cycle: int = -1) -> OperatingLimitsReport:
        """Pressure, wall-temperature and fluidisation extremes of each step of
        one cycle (-1 = last), with ``warnings`` where the cycle passes the
        wall's design pressure or temperature, goes below atmospheric without
        vacuum service, or the approach to fluidisation passes 1. The design
        checks run only when the wall has a design (a positive
        ``wall.design_pressure`` or ``wall.vacuum_service``)."""
        return self._sim.operating_limits(cycle)

    def wall_design(self) -> WallDesignResult:
        """The wall's minimum thickness [m] for its pressure design
        (``wall.design_code``, ``design_pressure``, ``design_temperature``,
        ``weld_joint_efficiency``, ``corrosion_allowance``,
        ``vacuum_service``): ``required_thickness``, or ``valid`` False with
        a ``message``. Set ``sim.wall.thickness`` to it to size the wall."""
        return self._sim.wall_design()

    def series(
        self,
        variable: str,
        component: Optional[str] = None,
        location: str = "top",
        z: Optional[float] = None,
        cycle: int = -1,
    ) -> TimeSeries:
        """A log variable over one cycle at one point of the bed.

        ``location`` is "top", "bottom" or "position"; passing ``z`` [m]
        implies "position". Time is measured from the start of the cycle.
        """
        if z is not None:
            location = "position"
        return self._sim.series(variable, component or "", location, 0.0 if z is None else float(z), cycle)

    def profile(self, variable: str, component: Optional[str] = None, cycle: int = -1, step: int = 0) -> AxialProfile:
        """Axial profiles over one step: ``values`` has shape (len(time), len(z))."""
        return self._sim.profile(variable, component or "", cycle, step)

    def stream_series(
        self, step: int, stream: str, variable: str, component: Optional[str] = None, cycle: int = -1
    ) -> TimeSeries:
        """A boundary-stream log variable over one step (time from the step start)."""
        return self._sim.stream_series(step, stream, variable, component or "", cycle)

    def endpoint_series(
        self, endpoint_id: int, variable: str, component: Optional[str] = None, cycle: int = -1
    ) -> TimeSeries:
        """:meth:`stream_series` for the endpoint with this stable ID (see :meth:`streams`)."""
        return self._sim.endpoint_series(endpoint_id, variable, component or "", cycle)

    def save_logs(self, directory: PathLike) -> None:
        """Write the full equation log as CSV files."""
        self._sim.save_logs(os.fspath(directory))

    def __repr__(self) -> str:
        return f"<skarstrom.Simulation {len(self.steps)} steps, components {self.components}>"


def bulk_variables() -> List[str]:
    """Names accepted by :meth:`Simulation.series` / :meth:`Simulation.profile` without a component."""
    return _core.bulk_variables()


def component_variables() -> List[str]:
    """Names accepted by :meth:`Simulation.series` / :meth:`Simulation.profile` with a component."""
    return _core.component_variables()


def stream_variables() -> List[str]:
    return _core.stream_variables()


def stream_component_variables() -> List[str]:
    return _core.stream_component_variables()


# ---- isotherms -------------------------------------------------------------------

def isotherm_kinds() -> List[str]:
    """Isotherm kinds that can be evaluated and fitted."""
    return _core.isotherm_kinds()


def isotherm_parameters(kind: str) -> List[str]:
    """Parameter names of an isotherm kind (SI units, as stored in decks)."""
    return _core.isotherm_parameters(kind)


def _as_list(x: Union[float, Sequence[float], np.ndarray], n: Optional[int] = None) -> List[float]:
    a = np.atleast_1d(np.asarray(x, dtype=float))
    if n is not None and a.size == 1 and n != 1:
        a = np.full(n, float(a[0]))
    return a.tolist()


def evaluate_isotherm(
    kind: str,
    params: Mapping[str, float],
    pressure: Union[float, Sequence[float], np.ndarray],
    temperature: Union[float, Sequence[float], np.ndarray],
    mole_fraction: Union[float, Sequence[float], np.ndarray] = 1.0,
) -> np.ndarray:
    """Equilibrium loading q [mol/kg] from the engine's isotherm code.

    ``pressure`` [Pa] is the total pressure; the component's partial pressure
    is ``mole_fraction * pressure``. Parameters left out keep the engine
    default.
    """
    p = _as_list(pressure)
    return _core.evaluate_isotherm(kind, dict(params), p, _as_list(temperature, len(p)), _as_list(mole_fraction, len(p)))


def evaluate_iast(
    kinds: Sequence[str],
    params: Sequence[Mapping[str, float]],
    pressure: float,
    temperature: float,
    mole_fraction: Sequence[float],
) -> np.ndarray:
    """Mixture loadings q_i [mol/kg] by Ideal Adsorbed Solution Theory.

    One pure isotherm per component: ``kinds[i]`` (see ``isotherm_kinds()``)
    with its parameters ``params[i]`` (those left out keep the engine default),
    at total ``pressure`` [Pa], ``temperature`` [K] and gas mole fractions
    ``mole_fraction``. An inert kind or a zero mole fraction gives zero loading.
    The engine's ``adsim/equations/iast.h``; see the isotherms page of the docs.
    """
    return _core.evaluate_iast(list(kinds), [dict(p) for p in params], float(pressure), float(temperature),
                               _as_list(mole_fraction, len(kinds)))


def fit_isotherm(
    kind: str,
    pressure: Union[Sequence[float], np.ndarray],
    loading: Union[Sequence[float], np.ndarray],
    temperature: Union[float, Sequence[float], np.ndarray],
    mole_fraction: Union[float, Sequence[float], np.ndarray] = 1.0,
    fixed: Optional[Mapping[str, float]] = None,
    initial: Optional[Mapping[str, float]] = None,
    bounds: Optional[Mapping[str, tuple]] = None,
    max_global_evals: int = 0,
) -> Dict[str, Any]:
    """Fit an isotherm kind to (pressure [Pa], loading [mol/kg], temperature [K]) data.

    Uses the engine's fitter (adsim/solver/isotherm_fit.h: a global MaxLIPO
    search, then bounded Levenberg-Marquardt). ``fixed`` holds parameters at given values; ``initial`` and
    ``bounds`` override the data-derived starting guesses and search bounds
    (for a pre-exponential factor paired with an energy, b in b*exp(Q/RT),
    the bounds apply to its value at the reference temperature t_ref).
    Returns a dict with params, std_errors, r2, rmse, sse, n, dof, t_ref,
    q_fit (model at each point) and any warning.
    """
    p = _as_list(pressure)
    q = _as_list(loading)
    t = _as_list(temperature, len(p))
    y = _as_list(mole_fraction, len(p))
    specs = None
    if fixed or initial or bounds:
        specs, _t_ref = _core.default_fit_specs(kind, p, q, t)
        names = {s["name"] for s in specs}
        for group in (fixed or {}), (initial or {}), (bounds or {}):
            for k in group:
                if k not in names:
                    raise ValueError(f"unknown parameter {k!r} for {kind}; expected one of {sorted(names)}")
        for s in specs:
            n = s["name"]
            if bounds and n in bounds:
                s["lower"], s["upper"] = (float(b) for b in bounds[n])
            if initial and n in initial:
                s["value"] = float(initial[n])
            if fixed and n in fixed:
                s["value"] = float(fixed[n])  # held: bounds are not used
                s["fixed"] = True
    return _core.fit_isotherm(kind, p, q, t, y, specs, 0.0, max_global_evals)
