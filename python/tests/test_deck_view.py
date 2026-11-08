"""Reference-client checks of the deck view for every bundled deck: it validates against the
engine's schema, and each value in it can be set back through the registry by its path."""

import json

import pytest

import skarstrom as sk

DECKS = [pytest.param(e, id=e.name) for e in sk.examples.ALL]


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


def _deck(sim):
    d = json.loads(sim.to_json())["simulation"]
    d.pop("System_Information", None)  # the save counter moves on every save
    return d


@pytest.mark.parametrize("example", DECKS)
def test_view_validates_against_the_schema(example):
    # jsonschema validates the view against the engine's schema; without it only the top level is checked.
    view = example.load().deck_view()
    schema = sk.deck_view_schema()
    try:
        import jsonschema
    except ImportError:
        assert set(schema["required"]) <= set(view)
        return
    jsonschema.validate(view, schema)


@pytest.mark.parametrize("example", DECKS)
def test_view_values_are_the_registry(example):
    sim = example.load()
    view = sim.deck_view()
    paths = [v["path"] for v in _values(view)]
    assert len(paths) == len(set(paths))
    assert set(paths) == {p.path for p in sim.parameters()}


@pytest.mark.parametrize("example", DECKS)
def test_every_view_value_sets_back_and_the_deck_is_unchanged(example):
    sim = example.load()
    before = _deck(sim)
    view = sim.deck_view()
    for v in _values(view):
        if v["kind"] in ("expression", "string"):
            sim.set(v["path"], v["text"])
        elif v["kind"] == "boolean":
            sim.set(v["path"], 1.0 if v["value"] else 0.0)
        else:
            assert isinstance(v["value"], (int, float)), v["path"]
            sim.set(v["path"], v["value"])
    assert _deck(sim) == before
    # The saved deck loads to the same deck. Compared by path: the engine's maps are unordered,
    # so a reloaded deck can list its component-keyed values in another order.
    again = sk.Simulation(sim.to_json())
    assert {v["path"]: (v["text"], v["unit"]) for v in _values(again.deck_view())} == {
        v["path"]: (v["text"], v["unit"]) for v in _values(view)}
    ok, msg = sim.check()
    assert ok, msg
