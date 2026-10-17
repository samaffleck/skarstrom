"""Compare two decks field by field (for the builder tests).

A deck is serialised by cereal; two decks that describe the same simulation
can still differ in bookkeeping: object IDs, connection pair IDs, the order of
unordered maps, polymorphic type numbering, embedded results, experimental
data, cost data and GUI display settings. :func:`normalise` removes those, and
the settings the engine ignores (an unused pump's values, an inactive buffer
tank, the start and duration of a pressure ramp, which the engine takes from
the bed state and the step time). What remains must be equal.
"""

from __future__ import annotations

import json
import math
import re
from typing import Any, Dict, List, Tuple

DROP = {
    "cereal_class_version", "DisplayUnitIndex", "Cycle_Endpoint_Summary", "Experimental_Data", "Cost_Estimator",
    "Next_ID", "Next_Layer_ID", "ID", "Endpoint_ID", "Pair_ID", "Step_ID", "Version",
    "Design_Pressure", "Design_Temperature", "Weld_Joint_Efficiency", "Vacuum_Service", "Proportional_Gain",
}
# Values the engine overrides at run time for pressure profiles (bed pressure, step time).
PRESSURE_PROFILE_KEYS = {"Boundary_Pressure", "Boundary_Value"}
RUNTIME_PROFILE_FIELDS = {"Start_Value", "Ramp_Time", "Step_Time"}
# Per-component maps where a missing component means 0 (the engine fills it in).
ZERO_DEFAULT_MAPS = {"Boundary_Mole_Fraction", "Adsorbed_Heat_Capacity"}
ZERO_DEFAULT_PATHS = (".y[", "adsorbed_heat_capacity[")


_QUANTITY = re.compile(r"^\s*[-+]?(\d+\.?\d*|\.\d+)([eE][-+]?\d+)?\s*\[[^\]]+\]\s*$")


def _root(text: str) -> Dict[str, Any]:
    d = json.loads(text)
    return d.get("value0", d.get("simulation", d))


class _Normaliser:
    def __init__(self) -> None:
        self.poly: Dict[int, str] = {}
        self.adsorbent_names: Dict[Any, str] = {}
        self.layer_names: Dict[Any, str] = {}

    def value(self, o: Any, key: str = "") -> Any:
        if isinstance(o, dict):
            if "polymorphic_id" in o:
                pid = o["polymorphic_id"]
                if "polymorphic_name" in o:
                    self.poly[pid & 0x7FFFFFFF] = o["polymorphic_name"]
                name = self.poly.get(pid & 0x7FFFFFFF, f"?{pid}")
                inner = o.get("ptr_wrapper", {}).get("data", {})
                out = self.value(inner, key)
                if isinstance(out, dict):
                    out = {"type": name, **out}
                return out
            if "ptr_wrapper" in o:
                return self.value(o["ptr_wrapper"].get("data"), key)
            if "Value" in o and "Expression" in o:
                text = o["Expression"]
                try:
                    float(text)
                    return o["Value"]
                except (TypeError, ValueError):
                    if _QUANTITY.match(text):  # "30.85 [kJ/mol]": a number written with a unit
                        return o["Value"]
                    return ("expr", text)
            if set(o) == {"index", "data"}:  # std::variant
                return {"variant": o["index"], **(self.value(o["data"]) or {})}
            out = {}
            for k, v in o.items():
                if k in DROP:
                    continue
                if key in PRESSURE_PROFILE_KEYS and k in RUNTIME_PROFILE_FIELDS:
                    continue
                out[k] = self.value(v, k)
            if out.get("Use_Pump") is False:  # an unused pump: its settings do nothing
                return {"Use_Pump": False}
            if out.get("Active") is False and "Percentage" in out:  # an inactive buffer tank
                return {"Active": False}
            return out
        if isinstance(o, list):
            if key in ZERO_DEFAULT_MAPS:
                pairs = ((str(x["key"]), self.value(x["value"], key)) for x in o)
                return {k: v for k, v in pairs if v != 0}
            if o and all(isinstance(x, dict) and set(x) == {"key", "value"} for x in o):
                return {str(self.key(x["key"], key)): self.value(x["value"], key) for x in o}
            return [self.value(x, key) for x in o]
        return o

    def key(self, k: Any, map_name: str) -> Any:
        if map_name == "Composition":
            return self.adsorbent_names.get(k, k)
        if map_name == "Energy_Sources":
            return self.layer_names.get(k, k)
        return k


def normalise(text: str, ignore_info: bool = False) -> Dict[str, Any]:
    """A deck's JSON with the bookkeeping removed (see the module docstring).

    ``ignore_info`` also drops the deck's title, description and author.
    """
    root = dict(_root(text))
    if ignore_info:
        root.pop("System_Information", None)
    n = _Normaliser()
    for entry in root.get("Adsorbent_Directory", {}).get("Adsorbents", []):
        n.adsorbent_names[entry["key"]] = entry["value"]["Name"]
    for layer in root.get("Bed", {}).get("Layers", {}).get("Layers", []):
        n.layer_names[layer.get("ID")] = layer.get("Name")
    used = set()
    for layer in root.get("Bed", {}).get("Layers", {}).get("Layers", []):
        used |= {str(x["key"]) for x in layer.get("Composition", []) if x["value"] > 0}
    out = n.value(root)
    ads = out.get("Adsorbent_Directory", {}).get("Adsorbents", {})
    if isinstance(ads, dict):  # key adsorbents by name, not ID
        for key, a in ads.items():
            if key not in used:  # no layer holds it, so the engine never reads its adsorption models
                a.pop("Adsorption_Source", None)
        out["Adsorbent_Directory"]["Adsorbents"] = {v["Name"]: v for v in ads.values()}
    return out


def _close(a: Any, b: Any, rel: float) -> bool:
    if isinstance(a, bool) or isinstance(b, bool):
        return a == b
    if isinstance(a, (int, float)) and isinstance(b, (int, float)):
        if math.isnan(a) and math.isnan(b):
            return True
        return a == b or abs(a - b) <= rel * max(abs(a), abs(b))
    return a == b


def diff(a: Any, b: Any, path: str = "", rel: float = 1e-12) -> List[Tuple[str, Any, Any]]:
    """Every (path, a, b) where two normalised decks differ."""
    out: List[Tuple[str, Any, Any]] = []
    if isinstance(a, dict) and isinstance(b, dict):
        for k in sorted(set(a) | set(b), key=str):
            if k not in a or k not in b:
                out.append((f"{path}/{k}", a.get(k, "<missing>"), b.get(k, "<missing>")))
            else:
                out += diff(a[k], b[k], f"{path}/{k}", rel)
    elif isinstance(a, list) and isinstance(b, list):
        if len(a) != len(b):
            out.append((f"{path}[len]", len(a), len(b)))
        for i, (x, y) in enumerate(zip(a, b)):
            out += diff(x, y, f"{path}[{i}]", rel)
    elif not _close(a, b, rel):
        out.append((path, a, b))
    return out


def endpoint_alias(sim: Any, endpoint_id: int) -> str:
    """A name for an endpoint that survives renumbering: step, end, stream position."""
    for s in sim.describe()["steps"]:
        for loc in ("bottom", "top"):
            end = s[loc]
            if end.get("endpoint") == endpoint_id:
                return f"{s['name']}/{loc}"
            for i, st in enumerate(end.get("streams", [])):
                if st["endpoint"] == endpoint_id:
                    return f"{s['name']}/{loc}/{i}"
    return f"?{endpoint_id}"


def connections(sim: Any) -> set:
    return {frozenset((endpoint_alias(sim, a), endpoint_alias(sim, b))) for a, b in sim.describe()["connections"]}


def parameter_values(sim: Any) -> Dict[str, Any]:
    """Every registry value keyed by alias (display names, so independent of IDs).

    Zero entries of per-component maps where a missing component means 0
    (backflow mole fractions, adsorbed-phase heat capacities) are left out.
    """
    d = sim.describe()
    used = {i for layer in d["layers"] for i, frac in layer["composition"] if frac > 0}
    unused = tuple(f"adsorbent[{a['name']}].{m}" for a in d["adsorbents"] if a["id"] not in used
                   for m in ("isotherm", "kinetics", "heat_of_adsorption"))
    out = {}
    for p in sim.parameters():
        value = p.text if p.kind == "string" else p.value
        if value == 0 and any(s in p.alias for s in ZERO_DEFAULT_PATHS) and not p.alias.startswith("initial."):
            continue
        if unused and p.alias.startswith(unused):
            continue
        out[p.alias] = value
    return out


def parameter_diff(a: Any, b: Any, rel: float = 1e-12) -> List[Tuple[str, Any, Any]]:
    """Every (alias, a, b) where two simulations' registry values differ."""
    pa, pb = parameter_values(a), parameter_values(b)
    return [(k, pa.get(k, "<missing>"), pb.get(k, "<missing>")) for k in sorted(set(pa) | set(pb), key=str)
            if k not in pa or k not in pb or not _close(pa[k], pb[k], rel)]
