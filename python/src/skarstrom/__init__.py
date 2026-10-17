"""Skarstrom: pressure and temperature swing adsorption cycle simulation.

Load a deck (a JSON simulation file; see :func:`examples` for the bundled
ones), change parameters, run it to cyclic steady state and read the results::

    import skarstrom

    sim = skarstrom.Simulation("deck.json")
    sim["cycle.number_of_cycles"] = 20
    sim["step[Feed].time"] = 120.0          # by name, or by ID path: "step[1].time"
    result = sim.run()
    print(result.status, result.css_cycle)
    for k in sim.kpis(cycle=-1):
        print(k.component, k.purity, k.recovery)
    y = sim.series("mole_fraction", "CarbonDioxide", location="top")
    y.time, y.values                         # NumPy arrays

Parameters are addressed by stable paths from the engine's parameter
registry; ``sim.parameters()`` lists them. See
https://samaffleck.github.io/skarstrom/python/ for the full guide.
"""

from __future__ import annotations

import json
import math
import os
import shutil
from dataclasses import dataclass
from typing import Any, Callable, Dict, Iterable, List, Mapping, Optional, Sequence, Union

import numpy as np

from . import _core
from ._core import AxialProfile, CSSRecord, KPIRecord, RunResult, TimeSeries

__version__: str = _core.__version__

__all__ = [
    "AxialProfile",
    "CSSRecord",
    "KPIRecord",
    "Parameter",
    "RunResult",
    "Simulation",
    "TimeSeries",
    "__version__",
    "adsorbent_library",
    "adsorbent_library_dir",
    "adsorbent_path",
    "bulk_variables",
    "copy_example",
    "component_variables",
    "evaluate_isotherm",
    "example_path",
    "examples",
    "fit_isotherm",
    "isotherm_kinds",
    "isotherm_parameters",
    "stream_component_variables",
    "stream_variables",
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
    unit: str
    kind: str  # "real", "expression", "integer", "boolean" or "string"
    tunable: bool
    lower: float
    upper: float
    value: Union[float, bool, str]
    text: str

    @property
    def bounds(self) -> Optional[tuple]:
        """Suggested (lower, upper) search bounds, or None."""
        if math.isnan(self.lower) or math.isnan(self.upper):
            return None
        return (self.lower, self.upper)


def _data_dir() -> str:
    return os.path.join(os.path.dirname(__file__), "data")


def _examples_dir() -> str:
    return os.path.join(_data_dir(), "decks")


def examples() -> List[str]:
    """File names of the example decks shipped with the package (sorted).

    The decks reproduce published PSA/VSA/TSA studies (see the examples page of
    the documentation); every one is an ordinary deck file that
    :class:`Simulation` loads.
    """
    d = _examples_dir()
    if not os.path.isdir(d):
        return []
    return sorted(f for f in os.listdir(d) if f.endswith(".json"))


def example_path(name: str) -> str:
    """Absolute path of a shipped example deck (see :func:`examples`).

    ``name`` may omit the ``.json`` suffix.
    """
    names = examples()
    if name not in names and name + ".json" in names:
        name += ".json"
    p = os.path.join(_examples_dir(), name)
    if not os.path.isfile(p):
        raise FileNotFoundError(f"no example deck {name!r}; available: {', '.join(names)}")
    return p


def copy_example(name: str, destination: PathLike = ".") -> str:
    """Copy a shipped example deck to ``destination`` (a directory or file path).

    Returns the path of the copy, which you can edit and run.
    """
    src = example_path(name)
    dst = os.fspath(destination)
    if os.path.isdir(dst):
        dst = os.path.join(dst, os.path.basename(src))
    shutil.copyfile(src, dst)
    return dst


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


class Simulation:
    """One deck and its solver.

    ``Simulation(deck)`` accepts a path to a deck file, the deck's JSON text,
    or nothing (an empty deck, mainly useful with :meth:`load`).
    """

    def __init__(self, deck: Optional[Union[PathLike, str]] = None) -> None:
        self._sim = _core.Simulation()
        if deck is not None:
            text = os.fspath(deck)
            if text.lstrip().startswith("{"):
                self._sim.load_json(text)
            else:
                self._sim.load(text)

    # ---- deck --------------------------------------------------------------
    @classmethod
    def from_file(cls, path: PathLike) -> "Simulation":
        sim = cls()
        sim._sim.load(os.fspath(path))
        return sim

    @classmethod
    def from_json(cls, text: str) -> "Simulation":
        sim = cls()
        sim._sim.load_json(text)
        return sim

    def load(self, path: PathLike) -> None:
        self._sim.load(os.fspath(path))

    def save(self, path: PathLike) -> None:
        """Write the deck (with any parameter changes) as a deck file."""
        self._sim.save(os.fspath(path))

    def to_json(self) -> str:
        return self._sim.to_json()

    def copy(self) -> "Simulation":
        """An independent copy of the deck (results are not copied)."""
        return Simulation.from_json(self.to_json())

    def check(self) -> tuple:
        """(ok, message): evaluates every expression and runs the deck checks."""
        return self._sim.check()

    # ---- parameters ----------------------------------------------------------
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

    def get(self, path: str) -> Union[float, int, str]:
        """A parameter's value (float, int for integers, str for text).

        Raises KeyError for a path the deck does not have.
        """
        try:
            return self._sim.get(path)
        except IndexError as e:
            raise KeyError(str(e)) from None

    def get_text(self, path: str) -> str:
        """A parameter's text: the expression for expression parameters."""
        try:
            return self._sim.get_text(path)
        except IndexError as e:
            raise KeyError(str(e)) from None

    def set(self, path: str, value: Union[float, int, str]) -> None:
        """Set a parameter.

        A number sets the value. A string sets an expression (``"2*feed"``,
        ``"5 [bar]"``) for expression parameters, or the text for text ones.
        """
        if isinstance(value, bool):
            value = int(value)
        try:
            if isinstance(value, (int, float, np.integer, np.floating)):
                self._sim.set_value(path, float(value))
            elif isinstance(value, str):
                self._sim.set_text(path, value)
            else:
                raise TypeError(f"cannot set {path} to a {type(value).__name__}")
        except IndexError as e:
            raise KeyError(str(e)) from None

    def update(self, values: Mapping[str, Union[float, int, str]]) -> None:
        """Set several parameters: ``sim.update({"step[1].time": 50, ...})``."""
        for k, v in values.items():
            self.set(k, v)

    def __getitem__(self, path: str) -> Union[float, int, str]:
        return self.get(path)

    def __setitem__(self, path: str, value: Union[float, int, str]) -> None:
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
    def add_component(self, name: str) -> None:
        """Add a fluid component (a CoolProp name such as ``"Argon"``).

        The component gets a zero mole fraction in every inlet, boundary and
        initial condition, and default physics on every adsorbent (inert
        isotherm, LDF rate model, zero heat of adsorption) for you to set.
        """
        self._sim.add_component(name)

    def remove_component(self, name: str) -> None:
        """Remove a component from the fluid, every adsorbent and every mole-fraction map."""
        self._sim.remove_component(name)

    def add_library_adsorbent(self, key: str) -> int:
        """Add an adsorbent from the bundled library by its key; returns its adsorbent ID.

        Components of the deck that the entry does not cover get default
        physics. A layer with no adsorbent yet gets it automatically; otherwise
        assign it to a layer in C++ (its parameters appear in
        :meth:`parameters` once a layer uses it). See :func:`adsorbent_library`.
        """
        return int(self._sim.add_library_adsorbent(adsorbent_library_dir(), key))

    # ---- structure -------------------------------------------------------------
    @property
    def components(self) -> List[str]:
        return self._sim.components

    @property
    def steps(self) -> List[str]:
        return self._sim.steps

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
