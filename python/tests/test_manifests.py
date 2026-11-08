"""The deck view, its schema and the display manifests (adsim/objects/deck_view.h, manifests.h)."""

import math
import os

import pytest

import skarstrom as sk
from skarstrom import cli


def _values(node):
    """Every registry value in a view (objects with path, kind and unit)."""
    if isinstance(node, dict):
        if {"path", "kind", "unit"} <= set(node):
            yield node
            return
        for v in node.values():
            yield from _values(v)
    elif isinstance(node, list):
        for v in node:
            yield from _values(v)


def test_deck_view_lists_every_registry_value_once():
    sim = sk.examples.Perez2019_A.load()
    view = sim.deck_view()
    assert view["view_version"] >= 1
    paths = [v["path"] for v in _values(view)]
    assert len(paths) == len(set(paths))
    assert set(paths) == {p.path for p in sim.parameters()}
    assert view["other"]["values"] == []
    assert view["steps"] and view["layers"] and view["adsorbents"]
    assert {c["name"] for c in view["components"]} == set(sim._sim.components)


def test_deck_view_set_back_by_path_leaves_the_deck_unchanged():
    sim = sk.examples.Perez2019_A.load()
    for v in _values(sim.deck_view()):
        if v["kind"] in ("expression", "string"):
            sim.set(v["path"], v["text"])
        elif v["kind"] == "boolean":
            sim.set(v["path"], 1.0 if v["value"] else 0.0)
        elif isinstance(v["value"], (int, float)):
            sim.set(v["path"], v["value"])
    after = sim.deck_view()
    fresh = sk.examples.Perez2019_A.load().deck_view()
    assert [(v["path"], v["text"]) for v in _values(after)] == [(v["path"], v["text"]) for v in _values(fresh)]


def test_deck_view_schema_names_the_kinds():
    schema = sk.deck_view_schema()
    assert schema["title"] == "DeckView"
    iso = schema["definitions"]["AdsorbentComponent"]["properties"]["isotherm"]["enum"]
    assert set(sk.build_kinds()["isotherm"]) <= set(iso)


def test_kpi_manifest_matches_kpi_records():
    kpis = sk.kpi_manifest()
    scalar = [k["key"] for k in kpis if not k["per_step"] and not k["diagnostic"] and not k["reaction"]]
    assert scalar == list(sk.sensitivity.kpi_names())
    # Reported KPIs first, then the diagnostics.
    flags = [k["diagnostic"] for k in kpis]
    assert flags == sorted(flags)
    record_fields = {f for f in dir(sk.KPIRecord) if not f.startswith("_")}
    for k in kpis:
        assert k["key"] in record_fields, k["key"]
        assert sk.check_unit(k["unit"])["ok"], k["unit"]
    # The CLI's columns come from the manifest; kpis.csv keeps its headers.
    assert [f for f, _ in cli.KPI_CSV_FIELDS] == scalar
    assert ("purity", "purity [mol/mol]") in cli.KPI_FIELDS
    assert "energy_consumption" not in [f for f, _ in cli.KPI_FIELDS]


def test_check_unit():
    ok = sk.check_unit("kPa", "pressure")
    assert ok["ok"] and ok["category"] == "pressure" and math.isclose(ok["factor"], 1000.0)
    assert not sk.check_unit("bar", "time")["ok"]
    bad = sk.check_unit("blargs")
    assert not bad["ok"] and bad["message"]
    units = sk.unit_manifest()
    for cat in units["categories"]:
        for u in cat["units"]:
            assert sk.check_unit(u["symbol"], cat["key"])["ok"], u["symbol"]


def test_warning_and_component_manifests():
    warnings = {w["key"]: w for w in sk.warning_manifest()}
    assert warnings["channel_reynolds"]["setting"] == "settings.warn_channel_reynolds"
    sim = sk.Simulation.new(title="t")
    for w in warnings.values():
        if w["setting"]:
            assert sim.get(w["setting"]) in (0.0, 1.0)
    manifest = sk.component_manifest()
    assert "NASA" in manifest["cp_reference"]
    comps = {c["name"]: c for c in manifest["components"]}
    for c in comps.values():
        assert bool(c["cp_species"]) == c["ideal_gas_heat_capacity"], c["name"]
    assert comps["CarbonDioxide"]["formula"] == "CO2"
    assert comps["Nitrogen"]["molar_mass"] == pytest.approx(0.028013, rel=1e-3)


def test_example_manifest_covers_the_bundled_decks():
    rows = sk.example_manifest()
    assert {r["name"] for r in rows} == {e.name for e in sk.examples.ALL}
    perez = next(r for r in rows if r["name"] == "Perez2019_A")
    assert perez["family"] == "Perez (2019)" and perez["family_key"] == "Perez2019"
    assert perez["summary"] and perez["hidden"] is False
    assert next(r for r in rows if r["name"] == "Illustrative_TSA_13X")["hidden"] is True


def test_kpi_manifest_covers_every_emitted_kpi_field():
    sim = sk.examples.Perez2019_A.load()
    sim.run(cycles=1)
    record = sim.kpis(cycle=-1)
    record = record[0] if isinstance(record, (list, tuple)) else record
    emitted = set(record.as_dict()) - {"cycle", "component", "valid"}
    assert emitted <= {k["key"] for k in sk.kpi_manifest()}


def test_manifests_are_not_empty():
    assert sk.kpi_manifest()
    assert sk.warning_manifest()
    assert sk.component_manifest()["components"]
    assert sk.example_manifest()
    units = sk.unit_manifest()
    assert units["categories"]
    # Every category lists suggested units except "custom" (any unit text the parser accepts).
    assert all(c["units"] for c in units["categories"] if c["key"] != "custom")


def test_every_warning_has_a_registry_switch():
    sim = sk.Simulation.new(title="t")
    registry = {p.path for p in sim.parameters()}
    switches = {p for p in registry if p.startswith("settings.warn_")}
    seen = set()
    for w in sk.warning_manifest():
        assert w["phase"] in ("start", "end"), w["key"]
        assert w["message"] and w["label"], w["key"]
        if w["setting"] is None:  # a warning without a switch is always on
            assert w["default_on"], w["key"]
            continue
        assert w["setting"] == f"settings.warn_{w['key']}", w["key"]
        assert w["setting"] in registry, w["key"]
        assert sim.get(w["setting"]) == (1.0 if w["default_on"] else 0.0), w["key"]
        seen.add(w["setting"])
    # No switch exists without a manifest entry.
    assert switches == seen


def test_every_example_file_in_the_manifest_exists():
    rows = sk.example_manifest()
    assert len({r["name"] for r in rows}) == len(rows)
    for r in rows:
        assert r["file"] == r["name"] + ".json"
        assert sk.examples.get(r["name"]).path.endswith(r["file"])
        assert os.path.isfile(sk.examples.get(r["name"]).path), r["file"]
        assert r["title"] and r["family"] and r["family_key"], r["name"]


def test_component_aliases_come_from_the_component_manifest():
    comps = {c["name"]: c for c in sk.component_manifest()["components"]}
    aliases = sk.component_aliases()
    assert aliases["CO2"] == "CarbonDioxide" and aliases["N2"] == "Nitrogen"
    for formula in ("CO2", "N2", "O2", "H2", "CH4", "He", "Ar", "H2O", "CO", "C2H6", "C2H4",
                    "C3H8", "C3H6", "H2S", "SO2", "NH3", "Kr", "Xe", "Ne"):
        assert aliases[formula] in comps, formula
    for alias, name in aliases.items():
        assert name in comps, alias
    # The typed objects resolve a formula through the same table.
    sim = sk.examples.Perez2019_A.load()
    assert sim.steps[0].bottom.y["CO2"] == sim.steps[0].bottom.y["CarbonDioxide"]


def test_python_keeps_no_copy_of_the_manifests():
    # The CLI's KPI columns and the sensitivity outputs are the manifest's, in its order.
    scalar = [k["key"] for k in sk.kpi_manifest() if not k["per_step"] and not k["diagnostic"] and not k["reaction"]]
    assert [f for f, _ in cli.KPI_CSV_FIELDS] == scalar
    from skarstrom import model
    assert not hasattr(model, "_COMPONENT_ALIASES")
