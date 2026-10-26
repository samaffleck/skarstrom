"""Objects for the parts of a deck: cycle steps and their streams, layers,
adsorbents, user variables and the settings groups.

They are views onto the engine's parameter registry, the same one behind
string paths (``sim["step[Feed].time"]``), the CLI and the optimiser. Reading
an attribute reads the deck; assigning one sets it, with the same checks.
Nothing is copied, so a view stays valid while the deck changes (and the IDs
it holds survive renames)::

    from skarstrom import ureg as u                   # pint's unit registry

    sim.steps["Feed"].time = 2 * u.min                # a step, by name
    sim.steps[0].bottom.set("temperature", 25, "degC")  # its feed-end stream
    sim.steps["Feed"].bottom.y["CO2"] = 0.15          # a mole fraction
    sim.layers[0].length = 1.2                        # metres (SI)
    sim.adsorbents["13X"].isotherm["CO2"].b1          # an isotherm parameter
    sim.variables["feed_flow"] = 2 * u.L / u.min      # a user variable
    sim.settings.reltol = 1e-6

Rules:

* Attributes are the fields every object of that kind has (``time``,
  ``length``); brackets select things you named (steps, streams, layers,
  adsorbents) or components. An integer in brackets is a position (0-based);
  ``.by_id(n)`` selects by the ID stored in the deck.
* A bare number is in the field's SI unit (see :mod:`skarstrom.units`); give
  a pint Quantity (``2 * u.min``) or ``set(name, 2, "min")`` to use another
  unit. Text is an expression (``"2 * t_feed"``) or a deck-syntax quantity.
* A field the deck does not have (``fraction`` when steps are timed, ``velocity``
  on a pressure boundary) raises ``AttributeError`` naming what it does have;
  so does a misspelling, so a typo never creates a new attribute silently.
* Component names are CoolProp names (``"CarbonDioxide"``); the usual formulas
  (``"CO2"``, ``"N2"``, ``"CH4"``...) are accepted too.
"""

from __future__ import annotations

import difflib
from typing import (
    TYPE_CHECKING,
    Any,
    Callable,
    Dict,
    Generic,
    Iterator,
    List,
    Optional,
    Tuple,
    TypeVar,
    Union,
    overload,
)

if TYPE_CHECKING:
    from . import Parameter, Simulation

__all__ = [
    "Adsorbent",
    "Bed",
    "CSS",
    "Collection",
    "Cycle",
    "EnergyBalance",
    "EnergySource",
    "Endpoint",
    "Fluid",
    "Info",
    "InitialConditions",
    "Isotherm",
    "Kinetics",
    "Layer",
    "NodeMap",
    "Pump",
    "Settings",
    "Step",
    "StepEnd",
    "BufferTank",
    "Valve",
    "ValueMap",
    "VoidRegion",
    "Wall",
]

#: What a field accepts: a number (SI), a pint Quantity, or a string (an
#: expression such as "2*feed_time", or a deck-syntax quantity such as "2 [min]").
Value = Union[float, int, str, Any]

# Formula-style names for the CoolProp components decks use most.
_COMPONENT_ALIASES = {
    "CO2": "CarbonDioxide",
    "N2": "Nitrogen",
    "O2": "Oxygen",
    "H2": "Hydrogen",
    "CH4": "Methane",
    "He": "Helium",
    "Ar": "Argon",
    "H2O": "Water",
    "CO": "CarbonMonoxide",
    "C2H6": "Ethane",
    "C2H4": "Ethylene",
    "C3H8": "n-Propane",
    "C3H6": "Propylene",
    "H2S": "HydrogenSulfide",
    "SO2": "SulfurDioxide",
    "NH3": "Ammonia",
    "Kr": "Krypton",
    "Xe": "Xenon",
    "Ne": "Neon",
}


# ---------------------------------------------------------------------------
# The registry, as a tree
# ---------------------------------------------------------------------------


def _split_path(path: str) -> List[Tuple[str, Optional[str]]]:
    """``"adsorbent[1].isotherm[CO2].b1"`` -> [("adsorbent", "1"), ("isotherm", "CO2"), ("b1", None)]."""
    parts: List[str] = []
    depth = start = 0
    for i, c in enumerate(path):
        if c == "[":
            depth += 1
        elif c == "]":
            depth = max(0, depth - 1)
        elif c == "." and depth == 0:
            parts.append(path[start:i])
            start = i + 1
    parts.append(path[start:])
    out: List[Tuple[str, Optional[str]]] = []
    for p in parts:
        b = p.find("[")
        if b > 0 and p.endswith("]"):
            out.append((p[:b], p[b + 1 : -1]))
        else:
            out.append((p, None))
    return out


class _Tree:
    """The registry paths under one object, grouped by shape."""

    __slots__ = ("leaves", "maps", "groups", "keyed")

    def __init__(self) -> None:
        self.leaves: Dict[str, str] = {}  # time -> "step[1].time"
        self.maps: Dict[str, Dict[str, str]] = {}  # y -> {Nitrogen: "endpoint[1].y[Nitrogen]"}
        self.groups: Dict[str, _Tree] = {}  # pump -> subtree
        self.keyed: Dict[str, Dict[str, _Tree]] = {}  # isotherm -> {CarbonDioxide: subtree}

    def insert(self, segments: List[Tuple[str, Optional[str]]], path: str) -> None:
        node = self
        for i, (name, key) in enumerate(segments):
            last = i == len(segments) - 1
            if last and key is None:
                node.leaves[name] = path
            elif last:
                node.maps.setdefault(name, {})[key] = path
            elif key is None:
                node = node.groups.setdefault(name, _Tree())
            else:
                node = node.keyed.setdefault(name, {}).setdefault(key, _Tree())

    def names(self) -> List[str]:
        return [*self.leaves, *self.maps, *self.groups, *self.keyed]

    def paths(self) -> Iterator[str]:
        yield from self.leaves.values()
        for m in self.maps.values():
            yield from m.values()
        for g in self.groups.values():
            yield from g.paths()
        for k in self.keyed.values():
            for t in k.values():
                yield from t.paths()


_EMPTY = _Tree()


class _Index:
    """One snapshot of the registry: parameters by path, the tree, object names."""

    def __init__(self, params: List["Parameter"], step_names: List[str], streams: List[Dict[str, Any]]) -> None:
        self.by_path: Dict[str, "Parameter"] = {p.path: p for p in params}
        self.by_alias: Dict[str, List["Parameter"]] = {}
        for p in params:
            if p.alias != p.path:
                self.by_alias.setdefault(p.alias, []).append(p)
        self.root = _Tree()
        self.names: Dict[Tuple[str, str], str] = {}
        for p in params:
            if " #" in p.path:  # a duplicate the registry disambiguated; string paths only
                continue
            segments = _split_path(p.path)
            self.root.insert(segments, p.path)
            kind, key = segments[0]
            if key is not None and kind in ("step", "layer", "adsorbent"):
                rest = p.path[len(kind) + len(key) + 2 :]  # after "kind[key]"
                if p.alias.startswith(kind + "[") and p.alias.endswith(rest) and p.alias != p.path:
                    self.names.setdefault((kind, key), p.alias[len(kind) + 1 : len(p.alias) - len(rest) - 1])
        self.step_names = step_names
        self.streams = streams

    def find(self, path: str) -> Optional["Parameter"]:
        p = self.by_path.get(path)
        if p is not None:
            return p
        matches = self.by_alias.get(path, [])
        return matches[0] if len(matches) == 1 else None


# ---------------------------------------------------------------------------
# Fields
# ---------------------------------------------------------------------------

T = TypeVar("T")


class _Field(Generic[T]):
    """A declared field: documents it (with its SI unit) for editors and help().

    Reading returns the deck's value; assigning sets it (see :data:`Value`).
    """

    def __init__(self) -> None:
        self.name = ""

    def __set_name__(self, owner: type, name: str) -> None:
        self.name = name

    @overload
    def __get__(self, obj: None, owner: Any = None) -> "_Field[T]": ...

    @overload
    def __get__(self, obj: "_Node", owner: Any = None) -> T: ...

    def __get__(self, obj: Optional["_Node"], owner: Any = None) -> Any:
        if obj is None:
            return self
        return obj._read(self.name)

    def __set__(self, obj: "_Node", value: Value) -> None:
        obj._write(self.name, value)


class Real(_Field[float]):
    """A number (or an expression, which reads back as its value)."""


class Integer(_Field[int]):
    """A whole number: a count or a model switch."""


class Boolean(_Field[bool]):
    """A switch."""


class Text(_Field[str]):
    """Free text."""


class Choice(_Field[str]):
    """One of a fixed set of names (a model or mode switch). Reads give the
    name; assign the name (or its number). ``parameters()`` lists the choices."""


def _format(value: Any, unit: str) -> str:
    if isinstance(value, bool) or isinstance(value, str):
        return repr(value)
    if isinstance(value, float):
        text = f"{value:.6g}"
    else:
        text = str(value)
    return f"{text} {unit}" if unit and unit != "-" else text


def _suggest(name: str, options: List[str]) -> str:
    close = difflib.get_close_matches(name, options, n=3)
    if close:
        return f"; did you mean {', '.join(close)}?"
    if options:
        return f"; it has: {', '.join(options)}"
    return ""


# ---------------------------------------------------------------------------
# Nodes
# ---------------------------------------------------------------------------


class _Node:
    """Base of every deck object: attribute access to its registry fields."""

    __slots__ = ("_sim", "_t", "_label", "_epoch")

    def __init__(self, sim: "Simulation", tree: Optional[_Tree], label: str) -> None:
        object.__setattr__(self, "_sim", sim)
        object.__setattr__(self, "_t", tree)
        object.__setattr__(self, "_label", label)
        object.__setattr__(self, "_epoch", sim._epoch)

    # ---- lookup ---------------------------------------------------------------
    @property
    def _tree(self) -> Optional[_Tree]:
        """The object's fields, refreshed when a structural edit changed the deck
        (a step, layer or adsorbent held in a variable stays usable)."""
        if self._epoch != self._sim._epoch:
            object.__setattr__(self, "_epoch", self._sim._epoch)
            try:
                fresh = self._refetch()
            except KeyError:  # the object was removed: keep what it had
                fresh = None
            if fresh is not None:
                for cls in type(self).__mro__:
                    for slot in getattr(cls, "__slots__", ()):
                        if slot not in ("_sim", "_epoch"):
                            object.__setattr__(self, slot, getattr(fresh, slot))
        return self._t

    def _refetch(self) -> Optional["_Node"]:
        """This object, read again from the current deck (None: keep the old fields)."""
        if type(self) in _TOP_LEVEL:
            return _group(self._sim, type(self), self._label)
        return None

    def _tree_or_raise(self) -> _Tree:
        tree = self._tree
        if tree is None:
            raise AttributeError(self._no_tree())
        return tree

    def _no_tree(self) -> str:
        return f"{self._label} has no parameters"

    def _missing(self, name: str) -> str:
        tree = self._tree if self._tree is not None else _EMPTY
        declared = isinstance(getattr(type(self), name, None), _Field)
        where = " in this deck" if declared else ""
        return f"{self._label} has no {name!r}{where}" + _suggest(name, tree.names())

    def _path_of(self, name: str) -> str:
        path = self._tree_or_raise().leaves.get(name)
        if path is None:
            raise AttributeError(self._missing(name))
        return path

    def _read(self, name: str) -> Any:
        return self._value(self._path_of(name))

    def _value(self, path: str) -> Any:
        """A registry value; the name for a choice (an enum-valued integer)."""
        value = self._sim.get(path)
        info = self._sim._index().by_path.get(path)
        if info is not None and info.choices and isinstance(value, int) and 0 <= value < len(info.choices):
            return info.choices[value]
        return value

    def _write(self, name: str, value: Value) -> None:
        self._sim.set(self._path_of(name), value)

    def __getattr__(self, name: str) -> Any:
        if name.startswith("_"):
            raise AttributeError(name)
        tree = self._tree_or_raise()
        if name in tree.leaves:
            return self._value(tree.leaves[name])
        if name in tree.maps:
            return ValueMap(self._sim, tree.maps[name], f"{self._label}.{name}")
        if name in tree.groups:
            return Group(self._sim, tree.groups[name], f"{self._label}.{name}")
        if name in tree.keyed:
            return NodeMap(self._sim, tree.keyed[name], f"{self._label}.{name}", Group)
        raise AttributeError(self._missing(name))

    def __setattr__(self, name: str, value: Any) -> None:
        if name.startswith("_"):
            object.__setattr__(self, name, value)
            return
        attr = getattr(type(self), name, None)
        if isinstance(attr, _Field):
            attr.__set__(self, value)
            return
        if isinstance(attr, property):
            raise AttributeError(f"{self._label}.{name} cannot be assigned")
        tree = self._tree_or_raise()
        if name in tree.leaves:
            self._sim.set(tree.leaves[name], value)
            return
        if name in tree.maps or name in tree.groups or name in tree.keyed:
            raise AttributeError(f"{self._label}.{name} holds several values; assign its entries, e.g. "
                                 f"{self._label}.{name}[...]")
        raise AttributeError(self._missing(name))

    def __dir__(self) -> List[str]:
        cls = type(self)
        own = [n for n in dir(cls) if not n.startswith("_") and not isinstance(getattr(cls, n, None), _Field)]
        tree = self._tree if self._tree is not None else _EMPTY
        return sorted(set(own) | set(tree.names()))

    # ---- public helpers ---------------------------------------------------------
    @property
    def fields(self) -> List[str]:
        """Names of the values this object has in this deck."""
        tree = self._tree if self._tree is not None else _EMPTY
        return list(tree.leaves)

    def get(self, name: str, unit: Optional[str] = None) -> Any:
        """A field's value, in its SI unit or converted to ``unit`` (``get("time", unit="min")``)."""
        if unit is None:
            return self._read(name)
        return self._sim.get(self._path_of(name), unit=unit)

    def set(self, name: str, value: Value, unit: Optional[str] = None) -> None:
        """Set a field; ``unit`` gives the unit of a plain number (``set("time", 2, unit="min")``)."""
        self._sim.set(self._path_of(name), value, unit=unit)

    def path(self, name: str) -> str:
        """The stable string path of a field (for ``sim[...]``, sweeps and the optimiser)."""
        return self._path_of(name)

    def text(self, name: str) -> str:
        """A field's text: its expression (``"2*feed_time"``, ``"25 [degC]"``) or formatted value."""
        return self._sim.get_text(self._path_of(name))

    def parameters(self) -> List["Parameter"]:
        """Every parameter under this object, with units, bounds and paths."""
        tree = self._tree if self._tree is not None else _EMPTY
        mine = set(tree.paths())
        return [p for p in self._sim.parameters() if p.path in mine]

    def _header(self) -> str:
        return f"<{self._label}>"

    def __repr__(self) -> str:
        lines = [self._header()]
        tree = self._tree
        if tree is None:
            lines.append(f"  ({self._no_tree()})")
            return "\n".join(lines)
        index = self._sim._index()
        rows: List[Tuple[str, str, str]] = []

        def row(label: str, path: str) -> None:
            p = index.by_path.get(path)
            unit = p.unit if p is not None else ""
            try:
                value = _format(self._sim.get(path), unit)
            except Exception as e:  # noqa: BLE001 - a repr must not raise
                value = f"<{e}>"
            rows.append((label, value, path))

        for name, path in tree.leaves.items():
            row(name, path)
        for name, entries in tree.maps.items():
            for key, path in entries.items():
                row(f"{name}[{key!r}]", path)
        for name in tree.groups:
            rows.append((name, f"(.{name})", ""))
        for name, entries in tree.keyed.items():
            rows.append((f"{name}[...]", ", ".join(repr(k) for k in entries), ""))
        width = max((len(r[0]) for r in rows), default=0)
        vwidth = max((len(r[1]) for r in rows if r[2]), default=0)
        for label, value, path in rows:
            lines.append(f"  {label:<{width}}  {value:<{vwidth}}  {path}".rstrip())
        return "\n".join(lines)


class Group(_Node):
    """A group of values with no fixed fields (fields depend on the model)."""


def _alias_key(alias: str, head: str) -> str:
    """The display name in ``...<head>[<name>]...`` of an alias (the last such bracket)."""
    i = alias.rfind(head + "[")
    if i < 0:
        return alias
    j, depth = i + len(head) + 1, 1
    k = j
    while k < len(alias) and depth:
        depth += {"[": 1, "]": -1}.get(alias[k], 0)
        k += 1
    return alias[j : k - 1]


def _named_paths(sim: "Simulation", paths: Dict[str, str], head: str) -> Dict[str, str]:
    """``paths`` (keyed by ID) re-keyed by the display name in each path's alias."""
    index = sim._index()
    out: Dict[str, str] = {}
    for key, path in paths.items():
        p = index.by_path.get(path)
        out[_alias_key(p.alias, head) if p is not None else key] = path
    return out


def _named_trees(sim: "Simulation", trees: Dict[str, "_Tree"], head: str) -> Dict[str, "_Tree"]:
    index = sim._index()
    out: Dict[str, _Tree] = {}
    for key, tree in trees.items():
        name = key
        for path in tree.paths():
            p = index.by_path.get(path)
            if p is not None:
                name = _alias_key(p.alias, head)
            break
        out[name] = tree
    return out


def _apply(node: "_Node", values: Dict[str, Any]) -> None:
    """Set several fields of a node (``None`` values are skipped)."""
    for name, value in values.items():
        if value is not None:
            node.set(name, value)


# ---------------------------------------------------------------------------
# Mappings
# ---------------------------------------------------------------------------


def _resolve_key(key: str, keys: List[str], what: str) -> str:
    if key in keys:
        return key
    alias = _COMPONENT_ALIASES.get(key)
    if alias in keys:
        return alias  # type: ignore[return-value]
    folded = [k for k in keys if k.lower() == str(key).lower()]
    if len(folded) == 1:
        return folded[0]
    raise KeyError(f"{what} has no {key!r}" + _suggest(str(key), keys))


class ValueMap:
    """Values keyed by name: mole fractions ``y[...]``, heats of adsorption,
    user variables. ``map[key]`` reads, ``map[key] = value`` sets."""

    __slots__ = ("_sim", "_paths", "_label")

    def __init__(self, sim: "Simulation", paths: Dict[str, str], label: str) -> None:
        self._sim = sim
        self._paths = paths
        self._label = label

    def _path(self, key: str) -> str:
        return self._paths[_resolve_key(key, list(self._paths), self._label)]

    def __getitem__(self, key: str) -> Any:
        return self._sim.get(self._path(key))

    def __setitem__(self, key: str, value: Value) -> None:
        self._sim.set(self._path(key), value)

    def get(self, key: str, unit: Optional[str] = None) -> Any:
        """The value of ``key``, in SI or converted to ``unit``."""
        return self._sim.get(self._path(key), unit=unit)

    def set(self, key: str, value: Value, unit: Optional[str] = None) -> None:
        """Set ``key``; ``unit`` gives the unit of a plain number."""
        self._sim.set(self._path(key), value, unit=unit)

    def path(self, key: str) -> str:
        """The stable string path of one entry."""
        return self._path(key)

    def text(self, key: str) -> str:
        return self._sim.get_text(self._path(key))

    def update(self, values: Dict[str, Value]) -> None:
        for k, v in values.items():
            self[k] = v

    def keys(self) -> List[str]:
        return list(self._paths)

    def values(self) -> List[Any]:
        return [self[k] for k in self._paths]

    def items(self) -> List[Tuple[str, Any]]:
        return [(k, self[k]) for k in self._paths]

    def to_dict(self) -> Dict[str, Any]:
        return dict(self.items())

    def __iter__(self) -> Iterator[str]:
        return iter(list(self._paths))

    def __len__(self) -> int:
        return len(self._paths)

    def __contains__(self, key: object) -> bool:
        try:
            _resolve_key(str(key), list(self._paths), self._label)
        except KeyError:
            return False
        return True

    def _ipython_key_completions_(self) -> List[str]:
        return list(self._paths)

    def __repr__(self) -> str:
        index = self._sim._index()
        lines = [f"<{self._label}>"]
        rows = []
        for key, path in self._paths.items():
            p = index.by_path.get(path)
            rows.append((f"[{key!r}]", _format(self._sim.get(path), p.unit if p is not None else ""), path))
        width = max((len(r[0]) for r in rows), default=0)
        vwidth = max((len(r[1]) for r in rows), default=0)
        for key, value, path in rows:
            lines.append(f"  {key:<{width}}  {value:<{vwidth}}  {path}")
        return "\n".join(lines)


N = TypeVar("N", bound=_Node)


class NodeMap(Generic[N]):
    """Objects keyed by component: ``adsorbent.isotherm["CO2"]``, ``adsorbent.kinetics["N2"]``."""

    __slots__ = ("_sim", "_trees", "_label", "_make")

    def __init__(self, sim: "Simulation", trees: Dict[str, _Tree], label: str, make: Callable[..., N]) -> None:
        self._sim = sim
        self._trees = trees
        self._label = label
        self._make = make

    def __getitem__(self, key: str) -> N:
        k = _resolve_key(key, list(self._trees), self._label)
        return self._make(self._sim, self._trees[k], f"{self._label}[{k!r}]")

    def keys(self) -> List[str]:
        return list(self._trees)

    def __iter__(self) -> Iterator[str]:
        return iter(list(self._trees))

    def __len__(self) -> int:
        return len(self._trees)

    def __contains__(self, key: object) -> bool:
        try:
            _resolve_key(str(key), list(self._trees), self._label)
        except KeyError:
            return False
        return True

    def _ipython_key_completions_(self) -> List[str]:
        return list(self._trees)

    def __repr__(self) -> str:
        return "\n\n".join(repr(self[k]) for k in self._trees) or f"<{self._label}: empty>"


# ---------------------------------------------------------------------------
# Collections of named objects
# ---------------------------------------------------------------------------


class Collection(Generic[N]):
    """Named deck objects (steps, layers, adsorbents, a step end's streams).

    ``c["Feed"]`` selects by name, ``c[0]`` by position (negative counts from
    the end), ``c.by_id(3)`` by the ID stored in the deck. Iterating yields the
    objects; ``c.names`` lists the names.
    """

    __slots__ = ("_items", "_label")

    def __init__(self, items: List[N], label: str) -> None:
        self._items = items
        self._label = label

    def __getitem__(self, key: Union[int, str]) -> N:
        if isinstance(key, bool):
            raise TypeError(f"{self._label} are selected by name (str) or position (int)")
        if isinstance(key, int):
            try:
                return self._items[key]
            except IndexError:
                raise IndexError(f"{self._label} has {len(self._items)} entries; no position {key}") from None
        if isinstance(key, str):
            matches = [x for x in self._items if x.name == key]
            if not matches:
                matches = [x for x in self._items if x.name.lower() == key.lower()]
            if len(matches) == 1:
                return matches[0]
            if len(matches) > 1:
                ids = ", ".join(str(x.id) for x in matches)
                raise KeyError(f"{self._label}: {len(matches)} entries are named {key!r} (IDs {ids}); "
                               f"select one by position or with .by_id()")
            raise KeyError(f"{self._label} has no {key!r}" + _suggest(key, self.names))
        raise TypeError(f"{self._label} are selected by name (str) or position (int), not {type(key).__name__}")

    def by_id(self, id: int) -> N:
        """The entry with this deck ID (the number in its string path, e.g. ``step[3]``)."""
        for x in self._items:
            if x.id == int(id):
                return x
        raise KeyError(f"{self._label} has no ID {id}; IDs: {', '.join(str(x.id) for x in self._items)}")

    @property
    def names(self) -> List[str]:
        return [x.name for x in self._items]

    def __iter__(self) -> Iterator[N]:
        return iter(list(self._items))

    def __len__(self) -> int:
        return len(self._items)

    def __contains__(self, key: object) -> bool:
        if isinstance(key, str):
            return any(x.name == key for x in self._items)
        return any(x is key for x in self._items)

    def _ipython_key_completions_(self) -> List[str]:
        return [x.name for x in self._items if x.name]

    def __repr__(self) -> str:
        lines = [f"<{self._label}: {len(self._items)}>"]
        for i, x in enumerate(self._items):
            lines.append(f"  [{i}] {x.name!r:<36} id {x.id}")
        return "\n".join(lines)


# ---------------------------------------------------------------------------
# Typed objects. Each declared field documents its SI unit in brackets; the
# test suite checks the brackets against the engine registry.
# ---------------------------------------------------------------------------


class Info(_Node):
    """The deck's descriptive information."""

    title = Text()
    """Title of the deck."""
    description = Text()
    """What the deck models and where it comes from."""
    author = Text()
    """Who made the deck."""


class Settings(_Node):
    """Integrator (SUNDIALS IDA) and logging settings."""

    init_step = Real()
    """Initial integrator step [s]."""
    min_step = Real()
    """Minimum integrator step [s]."""
    max_step = Real()
    """Maximum integrator step [s]."""
    reltol = Real()
    """Relative tolerance [-]."""
    nonlin_conv_coef = Real()
    """Nonlinear convergence coefficient [-]."""
    max_order = Integer()
    """Maximum BDF order."""
    max_nonlin_iters = Integer()
    """Maximum nonlinear iterations per step."""
    max_conv_fails = Integer()
    """Maximum convergence failures per step."""
    max_err_test_fails = Integer()
    """Maximum error-test failures per step."""
    max_num_steps = Integer()
    """Maximum internal steps per call."""
    max_iters_per_step = Integer()
    """Maximum iterations per step."""
    advection_scheme = Choice()
    """Advection discretisation: "upwind", "weno" or "van_leer"."""
    max_cycles_to_keep = Integer()
    """Cycles kept in the log."""
    max_log_points_per_step = Integer()
    """Maximum logged points per step."""
    warn_pressure_step_change = Boolean()
    """Warn at the start of a run about a constant pressure boundary far from the bed pressure."""
    warn_wall_thickness = Boolean()
    """Warn at the start of a run when the wall is thinner than its pressure design requires."""
    warn_design_pressure = Boolean()
    """Warn when the last cycle's pressure passes the wall's design pressure."""
    warn_vacuum = Boolean()
    """Warn when the last cycle goes below atmospheric and the wall is not in vacuum service."""
    warn_design_temperature = Boolean()
    """Warn when the last cycle's wall temperature passes the design temperature."""
    warn_fluidisation = Boolean()
    """Warn when the last cycle's approach to fluidisation passes 1."""
    acceleration_method = Choice()
    """Cyclic-steady-state acceleration: "none" or "aitken" (Graves-Morris extrapolation)."""


class CSS(_Node):
    """Cyclic-steady-state detection."""

    method = Choice()
    """CSS criterion: "none", "molar_profile", "mass_balance_total" or "mass_balance_component"."""
    min_cycles = Integer()
    """Minimum number of cycles before CSS can be declared."""
    molar_absolute_tolerance = Real()
    """Absolute tolerance on the cycle-to-cycle molar change [-]."""
    molar_relative_tolerance = Real()
    """Relative tolerance on the cycle-to-cycle molar change [-]."""
    mass_balance_tolerance = Real()
    """Mass-balance closure tolerance [-]."""


class EnergyBalance(_Node):
    """Energy-balance switches."""

    energy_balance_mode = Choice()
    """"non_isothermal" (solve the energy balances) or "isothermal"."""
    isothermal_use_initial_temperature = Boolean()
    """Isothermal runs take their temperature from the initial conditions."""
    isothermal_temperature = Real()
    """Temperature of an isothermal run [K]."""


class Cycle(_Node):
    """Cycle-level values."""

    number_of_cycles = Integer()
    """Number of cycles to run."""
    step_time_mode = Choice()
    """"time" (each step has a duration in seconds) or "fraction" (steps are fractions of total_time)."""
    total_time = Real()
    """Total cycle time when steps are given as fractions [s]."""


class Pump(_Node):
    """The pump, compressor or vacuum pump on a stream or boundary."""

    efficiency = Real()
    """Isentropic efficiency [-] (constant efficiency model)."""
    polytropic_index = Real()
    """Polytropic index of the compression [-]."""
    suction_pressure = Real()
    """Compressor suction pressure [Pa]."""
    discharge_pressure = Real()
    """Vacuum-pump discharge pressure [Pa]."""
    eta_max = Real()
    """Peak efficiency of the vacuum-pump efficiency curve [-]."""
    alpha = Real()
    """Decay coefficient of the vacuum-pump efficiency curve [-]."""
    poly_a = Real()
    """Efficiency polynomial coefficient a [-] (eta = a + b P + c P^2 + d P^3 + e P^4 + f P^5)."""
    poly_b = Real()
    """Efficiency polynomial coefficient b, of P [1/Pa]."""
    poly_c = Real()
    """Efficiency polynomial coefficient c, of P^2 [1/Pa^2]."""
    poly_d = Real()
    """Efficiency polynomial coefficient d, of P^3 [1/Pa^3]."""
    poly_e = Real()
    """Efficiency polynomial coefficient e, of P^4 [1/Pa^4]."""
    poly_f = Real()
    """Efficiency polynomial coefficient f, of P^5 [1/Pa^5]."""


class BufferTank(_Node):
    """The buffer tank an inlet draws from (:meth:`Endpoint.set_buffer_tank`)."""

    __slots__ = ("_source_of",)

    percentage = Real()
    """Share of the source endpoint's moles (from its step's last run) the inlet draws [%]."""

    def _refetch(self) -> Optional[_Node]:
        return _endpoint(self._sim, self._source_of).buffer_tank

    @property
    def source(self) -> "Endpoint":
        """The endpoint whose product fills the tank."""
        for s in self._sim.describe()["steps"]:
            for loc in ("bottom", "top"):
                end = s[loc]
                for e in [end] + end.get("streams", []):
                    if e.get("endpoint") == self._source_of and e.get("buffer_tank"):
                        return self._sim._endpoint(e["buffer_tank"])
        raise AttributeError(f"{self._label} has no source")


class Valve(_Node):
    """The valve of a valve stream."""

    cv = Real()
    """Flow coefficient Cv (US gpm per psi^0.5) [-]."""
    xt = Real()
    """Pressure-drop ratio factor xT [-]."""


class Endpoint(_Node):
    """A flow stream or a pressure boundary at one end of a step."""

    __slots__ = ("_id", "_name", "_location", "_step")

    temperature = Real()
    """Inlet temperature [K]."""
    temperature_linear_start = Real()
    """Inlet temperature at the start of a linear ramp [K]."""
    temperature_linear_end = Real()
    """Inlet temperature at the end of a linear ramp [K]."""
    temperature_exponential_start = Real()
    """Inlet temperature at the start of an exponential ramp [K]."""
    temperature_exponential_end = Real()
    """Inlet temperature at the end of an exponential ramp [K]."""
    temperature_linear_time = Real()
    """Duration of the inlet temperature linear ramp [s]."""
    temperature_exponential_time = Real()
    """Duration of the inlet temperature exponential ramp [s]."""
    velocity_linear_start = Real()
    """Superficial velocity at the start of a linear ramp [m/s]."""
    velocity_linear_end = Real()
    """Superficial velocity at the end of a linear ramp [m/s]."""
    velocity_linear_time = Real()
    """Duration of the velocity linear ramp [s]."""
    velocity_exponential_start = Real()
    """Superficial velocity at the start of an exponential ramp [m/s]."""
    velocity_exponential_end = Real()
    """Superficial velocity at the end of an exponential ramp [m/s]."""
    velocity_exponential_time = Real()
    """Duration of the velocity exponential ramp [s]."""
    molar_flow_linear_start = Real()
    """Molar flow at the start of a linear ramp [mol/s]."""
    molar_flow_linear_end = Real()
    """Molar flow at the end of a linear ramp [mol/s]."""
    molar_flow_linear_time = Real()
    """Duration of the molar flow linear ramp [s]."""
    molar_flow_exponential_start = Real()
    """Molar flow at the start of an exponential ramp [mol/s]."""
    molar_flow_exponential_end = Real()
    """Molar flow at the end of an exponential ramp [mol/s]."""
    molar_flow_exponential_time = Real()
    """Duration of the molar flow exponential ramp [s]."""
    volumetric_flow_linear_start = Real()
    """Volumetric flow at the start of a linear ramp [m3/s]."""
    volumetric_flow_linear_end = Real()
    """Volumetric flow at the end of a linear ramp [m3/s]."""
    volumetric_flow_linear_time = Real()
    """Duration of the volumetric flow linear ramp [s]."""
    volumetric_flow_exponential_start = Real()
    """Volumetric flow at the start of an exponential ramp [m3/s]."""
    volumetric_flow_exponential_end = Real()
    """Volumetric flow at the end of an exponential ramp [m3/s]."""
    volumetric_flow_exponential_time = Real()
    """Duration of the volumetric flow exponential ramp [s]."""
    mass_flow_linear_start = Real()
    """Mass flow at the start of a linear ramp [kg/s]."""
    mass_flow_linear_end = Real()
    """Mass flow at the end of a linear ramp [kg/s]."""
    mass_flow_linear_time = Real()
    """Duration of the mass flow linear ramp [s]."""
    mass_flow_exponential_start = Real()
    """Mass flow at the start of an exponential ramp [kg/s]."""
    mass_flow_exponential_end = Real()
    """Mass flow at the end of an exponential ramp [kg/s]."""
    mass_flow_exponential_time = Real()
    """Duration of the mass flow exponential ramp [s]."""
    velocity = Real()
    """Inlet superficial velocity [m/s]."""
    molar_flow = Real()
    """Inlet molar flow [mol/s]."""
    volumetric_flow = Real()
    """Inlet volumetric flow [m3/s]."""
    mass_flow = Real()
    """Inlet mass flow [kg/s]."""
    pressure = Real()
    """Constant boundary pressure [Pa]."""
    pressure_linear_end = Real()
    """Boundary pressure at the end of a linear ramp [Pa]."""
    pressure_exponential_end = Real()
    """Boundary pressure at the end of an exponential ramp [Pa]."""
    power = Real()
    """Shaft power of a pump-driven stream (constant power model) [W]."""
    power_poly_a = Real()
    """Pump power polynomial a, in a + b P + ... + f P^5 with P the boundary pressure [W]."""
    power_poly_b = Real()
    """Pump power polynomial b [W/Pa]."""
    power_poly_c = Real()
    """Pump power polynomial c [W/Pa^2]."""
    power_poly_d = Real()
    """Pump power polynomial d [W/Pa^3]."""
    power_poly_e = Real()
    """Pump power polynomial e [W/Pa^4]."""
    power_poly_f = Real()
    """Pump power polynomial f [W/Pa^5]."""

    def __init__(self, sim: "Simulation", tree: Optional[_Tree], label: str, id: int, name: str,
                 location: str, step: str) -> None:
        super().__init__(sim, tree, label)
        object.__setattr__(self, "_id", id)
        object.__setattr__(self, "_name", name)
        object.__setattr__(self, "_location", location)
        object.__setattr__(self, "_step", step)

    @property
    def id(self) -> int:
        """The endpoint ID stored in the deck (``endpoint[<id>]`` in string paths)."""
        return self._id

    def _refetch(self) -> Optional[_Node]:
        return _endpoint(self._sim, self._id)

    @property
    def name(self) -> str:
        return self._name

    @property
    def location(self) -> str:
        """``"bottom"`` (the feed end, z = 0) or ``"top"``."""
        return self._location

    @property
    def step(self) -> str:
        """Name of the step this endpoint belongs to."""
        return self._step

    @property
    def y(self) -> ValueMap:
        """Inlet mole fractions by component [-]."""
        return self._map("y")

    @property
    def pump(self) -> Pump:
        """The pump on this stream or boundary."""
        tree = self._tree_or_raise().groups.get("pump")
        if tree is None:
            raise AttributeError(self._missing("pump"))
        return Pump(self._sim, tree, f"{self._label}.pump")

    @property
    def valve(self) -> Valve:
        """The valve of a valve stream."""
        tree = self._tree_or_raise().groups.get("valve")
        if tree is None:
            raise AttributeError(self._missing("valve"))
        return Valve(self._sim, tree, f"{self._label}.valve")

    @property
    def buffer_tank(self) -> Optional[BufferTank]:
        """The buffer tank this inlet draws from, or ``None``."""
        tree = self._tree_or_raise().groups.get("buffer_tank")
        if tree is None:
            return None
        tank = BufferTank(self._sim, tree, f"{self._label}.buffer_tank")
        object.__setattr__(tank, "_source_of", self._id)
        return tank

    def set_buffer_tank(self, source: Optional["Endpoint"], percentage: Any = 100) -> "Endpoint":
        """Feed this inlet from a buffer tank holding what ``source`` (another
        endpoint, usually a product outlet) produced in its step's last run.

        The inlet draws ``percentage`` % of those moles, spread evenly over its
        own step, at their average temperature and composition: a well-mixed
        surge tank between two beds. Until the source step has run once, the
        inlet uses its own temperature, composition and flow. ``None`` removes
        the tank. An endpoint with a connection cannot have a buffer tank.
        """
        self._sim._build("set_buffer_tank", self._id, 0 if source is None else int(getattr(source, "id", source)))
        if source is not None:
            self._sim.set(f"endpoint[{self._id}].buffer_tank.percentage", percentage)
        return self._sim._endpoint(self._id)

    def set_pump(self, pump: Any) -> "Endpoint":
        """Fit a pump (a :class:`skarstrom.Pump`) to this stream or boundary; ``None`` removes it."""
        from . import build

        build._apply_pump(self._sim, self._id, pump)
        return self._sim._endpoint(self._id)

    def _map(self, name: str) -> ValueMap:
        paths = self._tree_or_raise().maps.get(name)
        if paths is None:
            raise AttributeError(self._missing(name))
        return ValueMap(self._sim, paths, f"{self._label}.{name}")

    def _header(self) -> str:
        return f"<{self._label} (endpoint id {self._id})>"


class StepEnd(Endpoint):
    """One end of a step (``step.bottom`` or ``step.top``): its streams, or its pressure boundary.

    ``end["Feed gas"]`` / ``end[0]`` selects one endpoint. When the end has
    exactly one, its fields are available directly: ``step.bottom.temperature``.
    """

    __slots__ = ("_eps", "_step_id")

    def __init__(self, sim: "Simulation", endpoints: List[Endpoint], label: str, step: str, location: str,
                 step_id: int = -1) -> None:
        only = endpoints[0] if len(endpoints) == 1 else None
        super().__init__(sim, only._tree if only else None, label, only.id if only else -1,
                         only.name if only else "", location, step)
        object.__setattr__(self, "_eps", endpoints)
        object.__setattr__(self, "_step_id", step_id)

    def _refetch(self) -> Optional[_Node]:
        return _steps(self._sim).by_id(self._step_id)._end(self._location)

    @property
    def _endpoints(self) -> List[Endpoint]:
        self._tree  # re-read after a structural edit
        return self._eps

    # ---- building ----------------------------------------------------------------
    def set_closed(self) -> None:
        """Close this end of the step (no flow)."""
        self._sim._build("set_closed", self._step_id, self._location)

    def set_pressure(self, direction: str, pressure: Any, temperature: Any = None, y: Optional[Dict[str, Any]] = None,
                     pump: Any = None) -> Endpoint:
        """Make this end a pressure boundary.

        ``direction`` is ``"inlet"`` or ``"outlet"``; ``pressure`` a value or a
        profile (:func:`skarstrom.constant`, :func:`skarstrom.linear`,
        :func:`skarstrom.exponential`). ``temperature`` and ``y`` are the gas
        that enters (for an outlet, when the flow reverses); ``pump`` a
        :class:`skarstrom.Pump`. Returns the new boundary.
        """
        from . import build

        return build._set_pressure(self._sim, self._step_id, self._location, direction, pressure, temperature, y, pump)

    def set_streams(self, direction: str, streams: List[Any]) -> List[Endpoint]:
        """Make this end a flow boundary with these streams (:class:`skarstrom.Stream`,
        :class:`skarstrom.ValveStream`). Returns the new streams, in order."""
        from . import build

        return build._set_streams(self._sim, self._step_id, self._location, direction, streams)

    def _no_tree(self) -> str:
        if not self._endpoints:
            return f"{self._label} is closed (no stream or pressure boundary)"
        names = ", ".join(repr(e.name) for e in self._endpoints)
        return (f"{self._label} has {len(self._endpoints)} streams ({names}); "
                f"choose one: {self._label}[{self._endpoints[0].name!r}] or {self._label}[0]")

    @property
    def endpoints(self) -> Collection[Endpoint]:
        """Every stream or boundary at this end."""
        return Collection(list(self._endpoints), self._label)

    def __getitem__(self, key: Union[int, str]) -> Endpoint:
        return self.endpoints[key]

    def __len__(self) -> int:
        return len(self._endpoints)

    def __iter__(self) -> Iterator[Endpoint]:
        return iter(list(self._endpoints))

    def _ipython_key_completions_(self) -> List[str]:
        return [e.name for e in self._endpoints if e.name]

    def __repr__(self) -> str:
        if len(self._endpoints) == 1:
            return repr(self._endpoints[0])
        return repr(self.endpoints) if self._endpoints else f"<{self._no_tree()}>"


class Step(_Node):
    """One step of the cycle."""

    __slots__ = ("_id", "_name", "_index")

    time = Real()
    """Duration of the step [s] (when steps are timed)."""
    fraction = Real()
    """Share of the total cycle time [-] (when steps are fractions of cycle.total_time)."""

    def __init__(self, sim: "Simulation", tree: Optional[_Tree], label: str, id: int, name: str, index: int) -> None:
        super().__init__(sim, tree, label)
        object.__setattr__(self, "_id", id)
        object.__setattr__(self, "_name", name)
        object.__setattr__(self, "_index", index)

    @property
    def id(self) -> int:
        """The step ID stored in the deck (``step[<id>]`` in string paths)."""
        return self._id

    def _refetch(self) -> Optional[_Node]:
        return _steps(self._sim).by_id(self._id)

    @property
    def name(self) -> str:
        return self._name

    @property
    def index(self) -> int:
        """Position of the step in the cycle (0-based), as used by results (``sim.profile(step=...)``)."""
        return self._index

    def _end(self, location: str) -> StepEnd:
        index = self._sim._index()
        root = index.root.keyed.get("endpoint", {})
        endpoints = []
        for s in index.streams:
            if s["step"] != self._index or s["location"] != location:
                continue
            eid = int(s["endpoint_id"])
            label = f"steps[{self._name!r}].{location}[{s['name']!r}]"
            endpoints.append(Endpoint(self._sim, root.get(str(eid), _Tree()), label, eid, s["name"], location,
                                      self._name))
        return StepEnd(self._sim, endpoints, f"steps[{self._name!r}].{location}", self._name, location, self._id)

    @property
    def bottom(self) -> StepEnd:
        """The feed end (z = 0): its stream(s) or pressure boundary."""
        return self._end("bottom")

    @property
    def top(self) -> StepEnd:
        """The product end (z = L): its stream(s) or pressure boundary."""
        return self._end("top")

    @property
    def energy_sources(self) -> "NodeMap[EnergySource]":
        """Heaters and coolers acting on a layer during this step, by layer name."""
        trees = self._tree_or_raise().keyed.get("energy", {})
        return NodeMap(self._sim, _named_trees(self._sim, trees, "energy"), f"{self._label}.energy_sources",
                       EnergySource)

    def add_energy_source(self, layer: Union["Layer", str, int], power: Value = None, enabled: bool = True) -> None:
        """Heat (positive ``power``) or cool (negative) a layer during this step [W]."""
        layer_id = self._sim._layer_id(layer)
        self._sim._build("set_energy_source", self._id, layer_id, bool(enabled))
        if power is not None:
            self._sim.set(f"step[{self._id}].energy[{layer_id}].power", power)

    def _header(self) -> str:
        return f"<{self._label} (position {self._index}, id {self._id})>"


class EnergySource(_Node):
    """A heater or cooler acting on one layer during one step."""

    enabled = Boolean()
    """Whether the source is on."""
    power = Real()
    """Total power into the layer [W] (negative = cooling)."""


class Bed(_Node):
    """Column geometry and surroundings."""

    diameter = Real()
    """Inner diameter of a cylindrical column [m]."""
    width = Real()
    """Width of a square or rectangular column [m]."""
    length = Real()
    """Depth of a rectangular column [m] (bed length is per layer: layers[i].length)."""
    ambient_temperature = Real()
    """Ambient temperature outside the wall [K]."""
    thermal_equilibrium_mode = Choice()
    """"lte" (gas and solid at one temperature) or "ltne" (separate gas and solid temperatures)."""
    particle_htc = Real()
    """Gas-to-particle heat transfer coefficient (non-equilibrium model) [W/m2-K]."""

    @property
    def bottom_void(self) -> "VoidRegion":
        """The empty volume below the bed."""
        return VoidRegion(self._sim, self._tree_or_raise().groups.get("bottom_void", _Tree()), "bed.bottom_void")

    @property
    def top_void(self) -> "VoidRegion":
        """The empty volume above the bed."""
        return VoidRegion(self._sim, self._tree_or_raise().groups.get("top_void", _Tree()), "bed.top_void")

    def set_shape(self, kind: str, **values: Value) -> "Bed":
        """``"cylinder"`` (diameter), ``"square"`` (width) or ``"rectangle"`` (width, length)."""
        self._sim._build("set_bed_shape", kind)
        bed = self._sim.bed
        _apply(bed, values)
        return bed

    def set_thermal_equilibrium(self, kind: str) -> None:
        """``"lte"`` (gas and solid at one temperature) or ``"ltne"``."""
        self._sim._build("set_thermal_equilibrium", kind)

    def set_particle_heat_transfer(self, kind: str = "constant", **values: Value) -> None:
        """The gas-to-particle heat transfer model (``"constant"``: particle_htc)."""
        self._sim._build("set_particle_heat_transfer", kind)
        _apply(self._sim.bed, values)


class VoidRegion(_Node):
    """An empty volume at one end of the bed."""

    enabled = Boolean()
    """Whether the void volume is modelled."""
    volume = Real()
    """Void volume [m3]."""
    resistance_coefficient = Real()
    """Flow resistance through the void [Pa.s/m]."""


class Wall(_Node):
    """The column wall."""

    density = Real()
    """Wall density [kg/m3]."""
    heat_capacity = Real()
    """Wall specific heat capacity [J/(kg-K)]."""
    thermal_conductivity = Real()
    """Wall thermal conductivity [W/(m-K)]."""
    thickness = Real()
    """Wall thickness [m]."""
    fixed_temperature = Real()
    """Wall temperature when it is held fixed [K]."""
    design_code = Choice()
    """Pressure design code: ``"asme_viii_1"`` (UG-27) or ``"asme_b31_3"`` (304.1.2); see
    :meth:`Simulation.wall_design`."""
    design_pressure = Real()
    """Design (maximum internal) pressure, gauge [Pa]; 0 = no pressure design."""
    design_temperature = Real()
    """Design temperature, for the allowable stress [K]."""
    weld_joint_efficiency = Real()
    """Weld joint efficiency E [-]: 1.0 seamless or fully radiographed, 0.85 spot, 0.70 none."""
    corrosion_allowance = Real()
    """Corrosion allowance added to the required thickness [m]."""
    vacuum_service = Boolean()
    """Also design the wall for full vacuum (external pressure)."""

    def set_material(self, kind: str, **values: Value) -> "Wall":
        """A wall material (``"ss316l"``, ``"carbon_steel_a106"``...; ``"custom"``): sets its
        density, heat capacity and conductivity, which ``values`` may then override."""
        self._sim._build("set_wall_material", kind)
        wall = self._sim.wall
        _apply(wall, values)
        return wall

    def set_boundary(self, kind: str, **values: Value) -> "Wall":
        """The outside of the wall: ``"zero_flux"`` (insulated by the ambient model) or ``"fixed"``
        (``fixed_temperature``)."""
        self._sim._build("set_wall_boundary", kind)
        wall = self._sim.wall
        _apply(wall, values)
        return wall


class InitialConditions(_Node):
    """The state of the bed at the start of the first cycle."""

    pressure = Real()
    """Initial pressure [Pa]."""
    temperature = Real()
    """Initial gas and solid temperature [K]."""
    wall_temperature = Real()
    """Initial wall temperature [K]."""
    velocity = Real()
    """Initial superficial velocity [m/s]."""

    @property
    def y(self) -> ValueMap:
        """Initial gas mole fractions by component [-]."""
        return ValueMap(self._sim, self._tree_or_raise().maps.get("y", {}), f"{self._label}.y")


class Fluid(_Node):
    """Properties of the gas and the adsorbed phase."""

    adsorbed_heat_capacity_mode = Choice()
    """"constant" (adsorbed_heat_capacity per component) or "gas_phase" (equal to the gas phase)."""
    heat_of_adsorption_convention = Choice()
    """Heats of adsorption are "isosteric" heats or follow the "internal_energy" convention."""

    @property
    def adsorbed_heat_capacity(self) -> ValueMap:
        """Adsorbed-phase molar heat capacity by component [J/(mol-K)]."""
        return ValueMap(self._sim, self._tree_or_raise().maps.get("adsorbed_heat_capacity", {}),
                        f"{self._label}.adsorbed_heat_capacity")


class Layer(_Node):
    """An adsorbent layer of the bed (feed end first)."""

    __slots__ = ("_id", "_name")

    length = Real()
    """Layer length [m]."""
    cells = Integer()
    """Number of finite-volume cells."""
    bed_void_fraction = Real()
    """Interparticle (bed) void fraction [-]."""
    particle_diameter = Real()
    """Particle diameter [m]."""
    momentum_model = Choice()
    """Pressure-drop model: "kozeny_carman" (Darcy) or "ergun"."""
    dispersion_coefficient = Real()
    """Constant axial dispersion coefficient [m2/s]."""
    molecular_diffusion_coefficient = Real()
    """Molecular diffusivity used by the dispersion correlation [m2/s]."""
    heat_dispersion_coefficient = Real()
    """Constant effective axial thermal conductivity lambda of the energy balance [W/(m-K)]."""
    internal_htc = Real()
    """Constant bed-to-wall heat transfer coefficient [W/m2-K]."""
    leva_eta1 = Real()
    """Leva correlation Re coefficient [-]."""
    leva_eta2 = Real()
    """Leva correlation Re exponent [-]."""
    leva_exp_coeff = Real()
    """Leva correlation wall-decay coefficient [-]."""
    external_htc = Real()
    """Constant wall-to-ambient heat transfer coefficient [W/m2-K]."""
    holman_c = Real()
    """Holman correlation coefficient [-]."""
    holman_m = Real()
    """Holman correlation Re exponent [-]."""
    holman_n = Real()
    """Holman correlation Pr exponent [-]."""
    external_velocity = Real()
    """Air velocity outside the column for the Holman correlation [m/s]."""

    def __init__(self, sim: "Simulation", tree: Optional[_Tree], label: str, id: int, name: str) -> None:
        super().__init__(sim, tree, label)
        object.__setattr__(self, "_id", id)
        object.__setattr__(self, "_name", name)

    @property
    def id(self) -> int:
        """The layer ID stored in the deck (``layer[<id>]`` in string paths)."""
        return self._id

    @property
    def composition(self) -> ValueMap:
        """Mass fraction of each adsorbent in the layer, by adsorbent name [-]."""
        paths = self._tree_or_raise().maps.get("mass_fraction", {})
        return ValueMap(self._sim, _named_paths(self._sim, paths, "mass_fraction"), f"{self._label}.composition")

    # ---- building ----------------------------------------------------------------
    def _fresh(self) -> "Layer":
        return self._sim.layers.by_id(self._id)

    _refetch = _fresh

    def set_composition(self, composition: Dict[Any, float]) -> "Layer":
        """The layer's adsorbents (names, IDs or Adsorbent objects) and their mass fractions."""
        pairs = [[self._sim._adsorbent_id(a), float(f)] for a, f in composition.items()]
        self._sim._build("set_layer_composition", self._id, pairs)
        return self._fresh()

    def set_momentum(self, kind: str) -> "Layer":
        """The pressure-drop model: ``"ergun"`` or ``"kozeny_carman"``."""
        self._sim._build("set_momentum", self._id, kind)
        return self._fresh()

    def _set_model(self, op: str, kind: str, values: Dict[str, Value]) -> "Layer":
        self._sim._build(op, self._id, kind)
        layer = self._fresh()
        _apply(layer, values)
        return layer

    def set_dispersion(self, kind: str, **values: Value) -> "Layer":
        """Axial dispersion: ``"constant"`` (dispersion_coefficient) or ``"ruthven"``
        (molecular_diffusion_coefficient)."""
        return self._set_model("set_dispersion", kind, values)

    def set_heat_dispersion(self, kind: str, **values: Value) -> "Layer":
        """Axial heat dispersion: ``"constant"`` (heat_dispersion_coefficient), ``"ruthven"`` or ``"yagi_kunii"``."""
        return self._set_model("set_heat_dispersion", kind, values)

    def set_internal_heat_transfer(self, kind: str, **values: Value) -> "Layer":
        """Bed-to-wall heat transfer: ``"constant"`` (internal_htc), ``"dixon"`` or ``"leva"``
        (leva_eta1, leva_eta2, leva_exp_coeff)."""
        return self._set_model("set_internal_heat_transfer", kind, values)

    def set_external_heat_transfer(self, kind: str, **values: Value) -> "Layer":
        """Wall-to-ambient heat transfer: ``"constant"`` (external_htc) or ``"holman"``
        (holman_c, holman_m, holman_n, external_velocity)."""
        return self._set_model("set_external_heat_transfer", kind, values)

    @property
    def name(self) -> str:
        return self._name

    def _header(self) -> str:
        return f"<{self._label} (id {self._id})>"


class Isotherm(_Node):
    """Isotherm parameters of one component. The fields depend on the isotherm
    kind (``b1``, ``q1``, ``m1`` for dual-site Langmuir; ``q_sat``, ``b``, ``v``
    for Sips...); ``fields`` lists them and the repr shows them with units."""


class Kinetics(_Node):
    """Mass-transfer (rate) model parameters of one component."""

    k = Real()
    """Mass transfer coefficient [1/s] ([kg/(mol-s)] for the LDF2 model)."""
    pre_exp_factor = Real()
    """Arrhenius pre-exponential factor of the mass transfer coefficient [1/s] ([kg/(mol-s)] for LDF2)."""
    activation_energy = Real()
    """Arrhenius activation energy of the mass transfer coefficient [J/mol]."""
    tortuosity = Real()
    """Particle tortuosity [-]."""
    macropore_void_fraction = Real()
    """Macropore void fraction [-]."""
    pore_radius = Real()
    """Mean pore radius [m]."""


class Adsorbent(_Node):
    """An adsorbent used by the bed's layers."""

    __slots__ = ("_id", "_name")

    solid_density = Real()
    """Particle (pellet) density [kg/m3]."""
    heat_capacity = Real()
    """Specific heat capacity [J/(kg-K)]."""
    thermal_conductivity = Real()
    """Thermal conductivity [W/(m-K)]."""
    particle_void_fraction = Real()
    """Intraparticle void fraction [-]."""

    def __init__(self, sim: "Simulation", tree: Optional[_Tree], label: str, id: int, name: str) -> None:
        super().__init__(sim, tree, label)
        object.__setattr__(self, "_id", id)
        object.__setattr__(self, "_name", name)

    @property
    def id(self) -> int:
        """The adsorbent ID stored in the deck (``adsorbent[<id>]`` in string paths)."""
        return self._id

    @property
    def name(self) -> str:
        return self._name

    @property
    def heat_of_adsorption(self) -> ValueMap:
        """Heat of adsorption by component [J/mol]."""
        paths = self._tree_or_raise().maps.get("heat_of_adsorption")
        if paths is None:
            raise AttributeError(self._missing("heat_of_adsorption"))
        return ValueMap(self._sim, paths, f"{self._label}.heat_of_adsorption")

    @property
    def isotherm(self) -> NodeMap[Isotherm]:
        """Isotherm parameters by component: ``isotherm["CO2"].b1``."""
        return NodeMap(self._sim, self._tree_or_raise().keyed.get("isotherm", {}), f"{self._label}.isotherm",
                       Isotherm)

    @property
    def kinetics(self) -> NodeMap[Kinetics]:
        """Rate-model parameters by component: ``kinetics["CO2"].k``."""
        return NodeMap(self._sim, self._tree_or_raise().keyed.get("kinetics", {}), f"{self._label}.kinetics",
                       Kinetics)

    # ---- building ----------------------------------------------------------------
    def _fresh(self) -> "Adsorbent":
        return self._sim.adsorbents.by_id(self._id)

    _refetch = _fresh

    def set_isotherm(self, kind: str, params: Optional[Dict[str, Dict[str, Value]]] = None,
                     components: Optional[List[str]] = None) -> "Adsorbent":
        """One isotherm of ``kind`` (see ``skarstrom.isotherm_kinds()``) for several components.

        ``params`` maps each component to its parameter values
        (``{"CO2": {"m1": 3.1, "b1": 2e-7 * u.m**3 / u.mol, ...}}``); ``components``
        defaults to the keys of ``params``. A competitive kind couples the
        components given together.
        """
        params = params or {}
        names = [self._sim._component(c) for c in (components if components is not None else params)]
        self._sim._build("set_isotherm", self._id, kind, names)
        ads = self._fresh()
        for component, values in params.items():
            _apply(ads.isotherm[self._sim._component(component)], values)
        return ads

    def set_rate_model(self, component: str, model: str = "ldf", coefficient: str = "constant",
                       **values: Value) -> "Adsorbent":
        """The mass-transfer model of one component: ``model`` ``"ldf"`` or ``"ldf2"``, with a
        ``coefficient`` model ``"constant"`` (k), ``"arrhenius"`` (pre_exp_factor, activation_energy),
        ``"macropore"`` / ``"haynes"`` (tortuosity, macropore_void_fraction) or
        ``"macropore_knudsen"`` (+ pore_radius)."""
        name = self._sim._component(component)
        self._sim._build("set_rate_model", self._id, name, model, coefficient)
        ads = self._fresh()
        _apply(ads.kinetics[name], values)
        return ads

    def _header(self) -> str:
        return f"<{self._label} (id {self._id})>"


# ---------------------------------------------------------------------------
# Builders used by Simulation
# ---------------------------------------------------------------------------

# The deck-wide groups (sim.bed, sim.cycle...), read again by name when stale.
_TOP_LEVEL = (Info, Settings, CSS, EnergyBalance, Cycle, Bed, Wall, InitialConditions, Fluid)


def _group(sim: "Simulation", cls: Callable[..., N], name: str) -> N:
    tree = sim._index().root.groups.get(name, _Tree())
    return cls(sim, tree, name)


def _steps(sim: "Simulation") -> Collection[Step]:
    index = sim._index()
    trees = index.root.keyed.get("step", {})
    items = []
    for i, (key, tree) in enumerate(trees.items()):
        name = index.names.get(("step", key))
        if name is None:
            name = index.step_names[i] if i < len(index.step_names) else key
        items.append(Step(sim, tree, f"steps[{name!r}]", int(key), name, i))
    return Collection(items, "steps")


def _layers(sim: "Simulation") -> Collection[Layer]:
    index = sim._index()
    items = []
    for key, tree in index.root.keyed.get("layer", {}).items():
        name = index.names.get(("layer", key), key)
        items.append(Layer(sim, tree, f"layers[{name!r}]", int(key), name))
    return Collection(items, "layers")


def _adsorbents(sim: "Simulation") -> Collection[Adsorbent]:
    index = sim._index()
    items = []
    for key, tree in index.root.keyed.get("adsorbent", {}).items():
        name = index.names.get(("adsorbent", key), key)
        items.append(Adsorbent(sim, tree, f"adsorbents[{name!r}]", int(key), name))
    return Collection(items, "adsorbents")


def _endpoint(sim: "Simulation", endpoint_id: int) -> Endpoint:
    """The endpoint (stream or pressure boundary) with this deck ID."""
    for step in _steps(sim):
        for end in (step.bottom, step.top):
            for e in end:
                if e.id == endpoint_id:
                    return e
    raise KeyError(f"no stream or pressure boundary with endpoint ID {endpoint_id}")


def _variables(sim: "Simulation") -> ValueMap:
    return ValueMap(sim, sim._index().root.maps.get("variable", {}), "variables")
