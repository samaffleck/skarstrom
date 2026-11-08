"""Display data the engine defines: the deck view's schema and the manifests.

Front ends read these instead of keeping their own copies (see
``adsim/objects/deck_view.h`` and ``adsim/objects/manifests.h``):

    import skarstrom as sk

    sim.deck_view()               # structure + every editable value (registry paths)
    sk.deck_view_schema()         # JSON Schema of the view
    sk.kpi_manifest()             # KPI key, label, unit, family, description
    sk.unit_manifest()            # unit categories and suggested units (parsed by the engine)
    sk.check_unit("kPa", "pressure")
    sk.warning_manifest()         # advisory warnings: switch, default, message template
    sk.component_manifest()       # formula, molar mass, CAS, heat-capacity source
    sk.example_manifest()         # example decks: title, family, hidden, summary
    sk.series_manifest()          # bed series variables: key, log name, unit, per-adsorbent
"""

from __future__ import annotations

import functools
import json
import os
from typing import Any, Dict, List, Optional

from . import _core

__all__ = [
    "check_unit",
    "component_aliases",
    "component_manifest",
    "deck_view_schema",
    "example_manifest",
    "kpi_manifest",
    "series_manifest",
    "unit_manifest",
    "warning_manifest",
]


def deck_view_schema() -> Dict[str, Any]:
    """The JSON Schema (draft 4) of :meth:`skarstrom.Simulation.deck_view`."""
    return json.loads(_core.deck_view_schema())


def kpi_manifest() -> List[Dict[str, Any]]:
    """Every KPI the engine reports, in report order: ``key`` (the :class:`KPIRecord` field),
    ``label``, ``short_label``, ``unit`` (SI), ``family``, ``description``, ``display``
    (shown in summary tables) and ``per_step`` (a list per step)."""
    return json.loads(_core.kpi_manifest())


def series_manifest() -> Dict[str, Any]:
    """The bed series variables :meth:`skarstrom.Simulation.series` and ``profile`` take:
    ``variables`` (``key``, ``label``, log ``name`` pattern such as ``"{component} Adsorbed
    Concentration"``, ``unit``, ``kind`` bulk/component, ``porous_only``, ``per_adsorbent`` and
    ``adsorbent_name``, the per-adsorbent log name pattern ``"{component} Adsorbed
    Concentration ({adsorbent})"``) and ``per_adsorbent_logged_when``."""
    return json.loads(_core.series_manifest())


def unit_manifest() -> Dict[str, Any]:
    """Unit categories (``key``, ``name``, ``si_unit``, suggested ``units`` with ``factor`` and
    ``offset`` from the engine's unit parser), the parser's symbol ``overrides`` and ``syntax``
    notes. Any unit the parser accepts can be used; check one with :func:`check_unit`."""
    return json.loads(_core.unit_manifest())


def check_unit(unit: str, category: Optional[str] = None) -> Dict[str, Any]:
    """Check unit text with the engine's parser: ``{ok, category, dimension, factor, offset,
    message}``. With ``category`` (e.g. ``"pressure"``), a unit of another category fails."""
    return dict(_core.check_unit(unit, category or ""))


def warning_manifest() -> List[Dict[str, Any]]:
    """The advisory warnings: ``key``, ``label``, ``phase`` (``start`` or ``end`` of a run),
    ``setting`` (the registry path of its switch, or None), ``default_on`` and ``message``."""
    return json.loads(_core.warning_manifest())


def component_manifest() -> Dict[str, Any]:
    """The components a deck can use: ``{"cp_reference": ..., "components": [...]}``. Each component has
    ``name`` (CoolProp), ``formula``, ``molar_mass`` [kg/mol], ``cas``, ``properties_source``,
    ``aliases``, ``ideal_gas_heat_capacity``, ``heat_capacity_held`` and the polynomial's provenance
    ``cp_species`` and ``cp_source`` (empty when the heat capacity is held); ``cp_reference`` is the
    table the polynomials come from."""
    return json.loads(_core.component_manifest())


@functools.lru_cache(maxsize=1)
def _component_aliases() -> Dict[str, str]:
    table: Dict[str, str] = {}
    comps = component_manifest()["components"]
    for pass_keys in (("formula",), ("aliases",)):  # formulas win over CoolProp aliases
        for c in comps:
            names = [c[pass_keys[0]]] if pass_keys[0] == "formula" else c[pass_keys[0]]
            for n in names:
                if n and n != c["name"]:
                    table.setdefault(n, c["name"])
    return table


def component_aliases() -> Dict[str, str]:
    """Formula / alias -> CoolProp name (``"CO2"`` -> ``"CarbonDioxide"``), read from the
    component manifest once per process. Decks accept these wherever a component is named."""
    return dict(_component_aliases())


def _data_dir() -> str:
    return os.path.join(os.path.dirname(__file__), "data")


def _example_manifest_path() -> Optional[str]:
    for path in (
        os.path.join(_data_dir(), "example_decks.json"),
        # A source checkout (the package next to examples/).
        os.path.join(os.path.dirname(__file__), "..", "..", "..", "examples", "example_decks.json"),
    ):
        if os.path.isfile(path):
            return os.path.normpath(path)
    return None


def example_manifest() -> List[Dict[str, Any]]:
    """The example decks: ``file``, ``name``, ``title``, ``family`` (display label, e.g.
    ``"Perez (2019)"``), ``family_key`` (``"Perez2019"``), ``hidden``, ``summary`` and ``description``. Raises FileNotFoundError when the package has no manifest."""
    path = _example_manifest_path()
    decks = os.path.join(_data_dir(), "decks")
    if path is None or not os.path.isdir(decks):
        raise FileNotFoundError("no example deck manifest is installed with this package")
    return json.loads(_core.example_deck_manifest(path, decks))
