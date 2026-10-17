"""The object layer (skarstrom.model), units (skarstrom.units) and skarstrom.examples."""

import ast
import json
import math
import os
import pathlib
import re

import pytest

import skarstrom
from skarstrom import model

EX = skarstrom.examples


@pytest.fixture
def perez():
    return EX.Perez2019_A.load()


# ---- examples -------------------------------------------------------------------


def test_examples_module_lists_every_deck_as_an_attribute():
    decks = sorted(p[:-5] for p in os.listdir(os.path.dirname(EX.Perez2019_A.path)) if p.endswith(".json"))
    assert EX.names() == decks and len(decks) >= 30
    for name in decks:
        assert name.isidentifier()
        ex = getattr(EX, name)
        assert isinstance(ex, skarstrom.Example) and ex.name == name and os.path.isfile(ex.path)
    assert [e.name for e in EX.ALL] == decks


def test_example_objects(tmp_path):
    ex = EX.Perez2019_A
    assert "Perez" in ex.description and ex.title
    assert os.fspath(ex) == ex.path
    assert EX.get("Perez2019_A") is ex and EX.get("Perez2019_A.json") is ex
    with pytest.raises(KeyError, match="did you mean Perez2019_A"):
        EX.get("Perez2019A")
    assert skarstrom.Simulation(ex).steps.names == ex.load().steps.names
    copy = ex.copy_to(tmp_path)
    assert os.path.basename(copy) == "Perez2019_A.json" and os.path.isfile(copy)


# ---- collections ----------------------------------------------------------------


def test_steps_by_name_position_and_id(perez):
    steps = perez.steps
    assert steps.names == ["Feed", "Blowdown", "Evacuation", "Pressurisation"]
    assert len(steps) == 4 and "Feed" in steps
    assert steps["Feed"].id == steps[0].id == steps.by_id(steps[0].id).id
    assert steps[-1].name == "Pressurisation" and steps[2].index == 2
    assert [s.name for s in steps] == steps.names
    assert steps["feed"].name == "Feed"  # case-insensitive fallback
    with pytest.raises(KeyError, match="did you mean Feed"):
        steps["Fed"]
    with pytest.raises(IndexError):
        steps[9]
    with pytest.raises(KeyError):
        steps.by_id(999)
    assert steps._ipython_key_completions_() == steps.names
    assert "Feed" in repr(steps)


def test_step_fields_read_and_write_the_deck(perez):
    feed = perez.steps["Feed"]
    feed.time = 120
    assert feed.time == 120.0 == perez["step[Feed].time"] == perez.get(feed.path("time"))
    assert feed.path("time") == "step[1].time"
    assert feed.fields == ["time"]
    assert "time" in dir(feed) and "fraction" not in dir(feed)


def test_typos_and_absent_fields_raise_instead_of_creating_attributes(perez):
    feed = perez.steps["Feed"]
    with pytest.raises(AttributeError, match="did you mean time"):
        feed.tiem = 3
    with pytest.raises(AttributeError, match="in this deck"):
        feed.fraction  # the cycle uses step times
    assert not hasattr(feed, "fraction")
    with pytest.raises(AttributeError, match="cannot be assigned"):
        feed.name = "x"
    with pytest.raises(AttributeError):
        perez.layers[0].lenght


# ---- units ------------------------------------------------------------------------


def test_bare_numbers_are_si_and_quantities_convert(perez):
    feed = perez.steps["Feed"]
    feed.time = "2 min"
    assert feed.time == 120.0 and feed.text("time") == "2 [min]"
    feed.time = "3 [min]"
    assert feed.time == 180.0
    feed.set("time", 1.5, unit="min")
    assert feed.time == 90.0
    assert feed.get("time", unit="min") == pytest.approx(1.5)
    feed.bottom.temperature = "25 degC"
    assert feed.bottom.temperature == pytest.approx(298.15)
    assert feed.bottom.get("temperature", unit="degC") == pytest.approx(25.0)
    perez.adsorbents[0].heat_capacity = "0.92 kJ/(kg*K)"
    assert perez.adsorbents[0].heat_capacity == pytest.approx(920.0)
    perez.set("endpoint[2].pressure", "1.2 bar")
    assert perez.get("endpoint[2].pressure", unit="bar") == pytest.approx(1.2)


def test_wrong_dimension_is_rejected_and_keeps_the_old_value(perez):
    feed = perez.steps["Feed"]
    feed.time = 100
    with pytest.raises(ValueError, match="time"):
        feed.time = "2 bar"
    with pytest.raises(ValueError):
        feed.set("time", 2, unit="K")
    assert feed.time == 100.0
    with pytest.raises(ValueError):
        feed.get("time", unit="Pa")


def test_require_units(perez):
    strict = skarstrom.Simulation(EX.Perez2019_A, require_units=True)
    feed = strict.steps["Feed"]
    with pytest.raises(TypeError, match="require_units"):
        feed.time = 120
    with pytest.raises(TypeError, match="require_units"):
        feed.time = "120"
    feed.time = "120 s"
    feed.set("time", 2, unit="min")
    assert feed.time == 120.0
    strict.settings.reltol = 1e-6  # dimensionless: a bare number is unambiguous
    strict.cycle.number_of_cycles = 2
    assert strict.copy().require_units


def test_pint_quantities():
    pint = pytest.importorskip("pint")
    ureg = pint.UnitRegistry()
    sim = EX.Perez2019_A.load()
    feed = sim.steps["Feed"]
    feed.time = 2 * ureg.minute
    assert feed.time == pytest.approx(120.0)
    feed.bottom.temperature = ureg.Quantity(25, ureg.degC)
    assert feed.bottom.temperature == pytest.approx(298.15)
    sim.layers[0].length = 120 * ureg.cm
    assert sim.layers[0].length == pytest.approx(1.2)
    sim.adsorbents[0].isotherm["CO2"].b1 = 1e-3 * ureg("m**3/mol")
    with pytest.raises(ValueError, match="time"):
        feed.time = 2 * ureg.bar
    strict = skarstrom.Simulation(EX.Perez2019_A, require_units=True)
    strict.steps[0].time = 30 * ureg.s


def test_units_module():
    u = skarstrom.units
    assert u.convert(25, "degC", "K") == pytest.approx(298.15)
    assert u.convert(1, "1/kPa", "1/bar") == pytest.approx(100.0)
    assert u.dimension_name("W/m2-K") == "heat transfer coefficient"
    assert u.parse("kJ/mol").factor == 1e3
    with pytest.raises(ValueError):
        u.convert(1, "bar", "s")
    with pytest.raises(ValueError):
        u.parse("furlong")


# ---- step ends, streams, components -------------------------------------------------


def test_step_ends_and_streams():
    sim = EX.FiveStepProductPurge.load()
    feed = sim.steps["Feed1"]
    assert feed.bottom.name == "Feed 1 inlet"
    feed.bottom.temperature = 300  # one stream: fields are direct
    assert feed.bottom["Feed 1 inlet"].temperature == 300.0
    top = feed.top
    assert len(top) == 2 and top.endpoints.names == ["Ambient Stream", "Feed 1 Press Stream"]
    with pytest.raises(AttributeError, match="choose one"):
        top.temperature
    assert top[0].name == "Ambient Stream" and top["Feed 1 Press Stream"].location == "top"
    assert all(e.step == "Feed1" for e in top)


def test_mole_fractions_accept_component_formulas(perez):
    y = perez.steps["Feed"].bottom.y
    assert set(y.keys()) == {"Nitrogen", "CarbonDioxide"}
    y["CO2"] = 0.2
    y["N2"] = "80 %"
    assert y["CarbonDioxide"] == pytest.approx(0.2) and y["Nitrogen"] == pytest.approx(0.8)
    assert "CO2" in y and "Argon" not in y
    with pytest.raises(KeyError):
        y["Argon"]
    assert y.to_dict() == {"Nitrogen": pytest.approx(0.8), "CarbonDioxide": pytest.approx(0.2)}


def test_layers_adsorbents_isotherms_and_kinetics(perez):
    layer = perez.layers[0]
    assert layer.name == perez.layers.names[0] and layer.length > 0
    layer.length = "120 cm"
    assert layer.length == pytest.approx(1.2)
    ads = perez.adsorbents[0]
    assert set(ads.isotherm.keys()) == {"CarbonDioxide", "Nitrogen"}
    iso = ads.isotherm["CO2"]
    assert "b1" in iso.fields and "b1" in dir(iso)
    iso.b1 = 2e-3
    assert iso.b1 == 2e-3 == perez[iso.path("b1")]
    with pytest.raises(AttributeError, match="did you mean"):
        iso.bb1 = 1
    ads.heat_of_adsorption["CO2"] = "35 kJ/mol"
    assert ads.heat_of_adsorption["CarbonDioxide"] == pytest.approx(35e3)
    assert ads.kinetics["N2"].tortuosity > 0
    assert len(ads.parameters()) > 10


def test_groups(perez):
    perez.cycle.number_of_cycles = 3
    assert perez.cycle.number_of_cycles == 3 and isinstance(perez.cycle.number_of_cycles, int)
    perez.settings.reltol = 1e-6
    assert perez.settings.reltol == 1e-6
    perez.energy_balance.isothermal_use_initial_temperature = True
    assert perez.energy_balance.isothermal_use_initial_temperature is True
    perez.info.title = "mine"
    assert perez.info.title == "mine"
    assert perez.bed.diameter > 0 and perez.wall.thickness > 0
    for node in (perez.cycle, perez.settings, perez.css, perez.energy_balance, perez.info, perez.bed,
                 perez.wall, perez.steps[0], perez.steps[0].bottom, perez.layers[0], perez.adsorbents[0],
                 perez.adsorbents[0].isotherm, perez.steps[0].bottom.y):
        assert repr(node)


def test_user_variables(tmp_path):
    deck = json.loads(open(EX.Perez2019_A.path, encoding="utf-8").read())
    root = deck.get("value0", deck)
    root["Variable_Registry"]["Variables"] = [
        {"cereal_class_version": 3, "Name": "t_feed", "Value": 100.0, "ValueExpr": "", "Description": "",
         "UnitCategory": 1},
        {"Name": "k_any", "Value": 2.0, "ValueExpr": "", "Description": "", "UnitCategory": 8},
    ]
    sim = skarstrom.Simulation(json.dumps(deck))
    v = sim.variables
    assert v.keys() == ["t_feed", "k_any"]
    v["t_feed"] = "2 min"
    assert v["t_feed"] == 120.0
    with pytest.raises(ValueError):
        v["t_feed"] = "2 bar"
    v["k_any"] = 5  # a custom unit: not checked
    sim.steps["Feed"].time = "t_feed"
    assert sim.steps["Feed"].time == 120.0


# ---- the object layer covers the whole registry ----------------------------------------------


def _walk(node, seen):
    """Every registry path reachable from a node, through attributes only."""
    if isinstance(node, model.ValueMap):
        for k in node.keys():
            seen.add(node.path(k))
        return
    if isinstance(node, model.NodeMap):
        for k in node.keys():
            _walk(node[k], seen)
        return
    if isinstance(node, model.StepEnd):
        for e in node:
            _walk(e, seen)
        return
    for name in node.fields:
        seen.add(node.path(name))
        getattr(node, name)  # every field reads
    tree = node._tree
    for name in [*tree.maps, *tree.groups, *tree.keyed]:
        _walk(getattr(node, name), seen)
    if isinstance(node, model.Step):
        _walk(node.bottom, seen)
        _walk(node.top, seen)


@pytest.mark.parametrize("example", EX.ALL, ids=lambda e: e.name)
def test_every_registry_parameter_is_reachable_as_an_attribute(example):
    sim = example.load()
    seen = set()
    for group in (sim.cycle, sim.settings, sim.css, sim.energy_balance, sim.info, sim.bed, sim.wall):
        _walk(group, seen)
    for coll in (sim.steps, sim.layers, sim.adsorbents):
        for item in coll:
            _walk(item, seen)
    _walk(sim.variables, seen)
    paths = {p.path for p in sim.parameters() if " #" not in p.path}
    assert paths - seen == set(), "parameters the object layer cannot reach"
    assert seen <= paths


def _declared_field_docs():
    """(class name, field) -> docstring, from the attribute docstrings in model.py."""
    tree = ast.parse(pathlib.Path(model.__file__).read_text(encoding="utf-8"))
    out = {}
    for cls in [n for n in tree.body if isinstance(n, ast.ClassDef)]:
        body = cls.body
        for a, b in zip(body, body[1:]):
            if (isinstance(a, ast.Assign) and isinstance(a.value, ast.Call)
                    and isinstance(a.value.func, ast.Name) and a.value.func.id in ("Real", "Integer", "Boolean", "Text")):
                doc = b.value.value if isinstance(b, ast.Expr) and isinstance(b.value, ast.Constant) else ""
                out[(cls.name, a.targets[0].id)] = doc
    return out


def test_declared_fields_document_the_registry_unit():
    docs = _declared_field_docs()
    assert docs, "no declared fields found"
    owners = {"Step": lambda s: list(s.steps), "Layer": lambda s: list(s.layers),
              "Adsorbent": lambda s: list(s.adsorbents), "Bed": lambda s: [s.bed], "Wall": lambda s: [s.wall],
              "Cycle": lambda s: [s.cycle], "Settings": lambda s: [s.settings], "CSS": lambda s: [s.css],
              "EnergyBalance": lambda s: [s.energy_balance], "Info": lambda s: [s.info],
              "Endpoint": lambda s: [e for st in s.steps for end in (st.bottom, st.top) for e in end],
              "Pump": lambda s: [e.pump for st in s.steps for end in (st.bottom, st.top) for e in end
                                 if "pump" in e._tree.groups],
              "Kinetics": lambda s: [a.kinetics[k] for a in s.adsorbents for k in a.kinetics.keys()]}
    checked = set()
    for example in EX.ALL:
        sim = example.load()
        index = sim._index()
        for cls_name, getter in owners.items():
            for node in getter(sim):
                for field in node.fields:
                    doc = docs.get((cls_name, field))
                    assert doc is not None, f"{cls_name}.{field} is in the registry but not declared in model.py"
                    unit = index.by_path[node.path(field)].unit
                    if unit not in ("", "-") and index.by_path[node.path(field)].kind in ("real", "expression"):
                        assert f"[{unit}]" in doc, f"{cls_name}.{field}: registry unit {unit!r}, doc {doc!r}"
                    checked.add((cls_name, field))
    assert len(checked) > 40


def test_object_helpers_do_not_shadow_registry_fields():
    reserved = {"get", "set", "path", "text", "fields", "parameters", "name", "id", "index", "location", "step",
                "endpoints", "bottom", "top", "keys", "values", "items", "update"}
    for example in EX.ALL:
        for p in example.load().parameters():
            leaf = re.sub(r"\[.*?\]", "", p.path.split(".")[-1])
            assert leaf not in reserved, f"{p.path} would be hidden by a helper"


def test_index_refreshes_after_structural_edits():
    sim = EX.Perez2019_A.load()
    assert "Argon" not in sim.steps["Feed"].bottom.y
    sim.add_component("Argon")
    assert "Argon" in sim.steps["Feed"].bottom.y
    sim.remove_component("Argon")
    assert "Argon" not in sim.steps["Feed"].bottom.y


def test_nan_free_reads():
    sim = EX.Perez2019_A.load()
    for p in sim.parameters():
        if p.kind in ("real", "expression"):
            assert isinstance(sim.get(p.path), float) and not math.isnan(sim.get(p.path))
