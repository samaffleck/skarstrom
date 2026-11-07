"""The deck view, its schema and the display manifests (adsim/objects/deck_view.h, manifests.h)."""

import math

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
    scalar = [k["key"] for k in kpis if not k["per_step"] and not k["diagnostic"]]
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
