"""sim.edit: the engine's structural cascades (adsim/objects/deck_edit.h), and skarstrom.model_kinds()."""

import json

import pytest

import skarstrom

EX = skarstrom.examples


@pytest.fixture
def perez():
    return EX.Perez2019_A.load()


def test_model_kinds_lists_every_family_with_parameters():
    kinds = skarstrom.model_kinds()
    for family in ("rate_coefficient", "dispersion", "heat_dispersion", "internal_heat_transfer",
                   "external_heat_transfer", "particle_heat_transfer", "profile", "flow", "pump",
                   "pump_efficiency", "pump_power", "stream", "boundary", "isotherm"):
        assert kinds[family], family
    coeff = {k["tag"]: k for k in kinds["rate_coefficient"]}
    assert [p["name"] for p in coeff["macropore"]["parameters"]] == [
        "tortuosity", "macropore_void_fraction", "molecular_diffusion_coefficient"]
    iso = {k["tag"]: k for k in kinds["isotherm"]}
    assert iso["dual_site_langmuir_p"]["paired"] == "extended_dual_site_langmuir_p"
    assert iso["toth_wadst"]["water_coupled"] and not iso["toth_wadst"]["iast"]
    # The builder's tag lists are the descriptors' tags.
    build = skarstrom.build_kinds()
    for family in ("rate_coefficient", "dispersion", "flow", "pump", "pump_efficiency", "pump_power", "profile"):
        # Legacy kinds (the former rate coefficients, now presets) are not offered.
        assert build[family] == [k["tag"] for k in kinds[family] if not k.get("legacy")], family


def test_layers_and_steps(perez):
    sim = perez
    n_layers, n_steps = len(sim.layers), len(sim.steps)
    layer = sim.edit.add_layer()
    assert len(sim.layers) == n_layers + 1
    assert sim.layers.by_id(layer).name == f"Layer {n_layers + 1}"
    sim.edit.move_layer(layer, -1)
    with pytest.raises(ValueError):
        sim.edit.move_layer(layer, 5)
    sim.edit.remove_layer(layer)
    assert len(sim.layers) == n_layers

    step = sim.edit.add_step("Rinse")
    assert sim.steps[-1].name == "Rinse" and sim.steps[-1].id == step
    copy = sim.edit.duplicate_step(step)
    assert sim.steps[-1].id == copy
    sim.edit.move_step(copy, -1)
    assert sim.steps[-2].id == copy
    sim.edit.remove_step(sim.steps["Rinse"])
    sim.edit.remove_step(copy)
    assert len(sim.steps) == n_steps
    ok, msg = sim.check()
    assert ok, msg


def test_streams_and_kind_switches(perez):
    sim = perez
    step = sim.edit.add_step()
    valve = sim.edit.add_stream(step, "bottom", "valve")
    assert sim[f"endpoint[{valve}].valve.cv"] == pytest.approx(100.0)
    with pytest.raises(ValueError, match="flow boundary"):
        sim.edit.add_stream(step, "top", "flow")
    with pytest.raises(ValueError, match="unknown stream kind"):
        sim.edit.add_stream(step, "bottom", "pipe")

    # A stream keeps its endpoint ID through a device switch.
    sim.edit.set_stream_kind(step, "bottom", 1, "pump")
    assert f"endpoint[{valve}].power" in [p.path for p in sim.parameters()]
    sim.edit.set_pump_power_kind(valve, "polynomial")
    assert sim.get(f"endpoint[{valve}].power_poly_a") == 0.0
    sim.edit.remove_stream(step, "bottom", 1)

    # Pressure profile: the constant's value becomes the ramp's end.
    top = sim.steps.by_id(step).top.id
    sim[f"endpoint[{top}].pressure"] = 2.5e5
    sim.edit.set_profile_kind(top, "pressure", "linear")
    assert sim[f"endpoint[{top}].pressure_linear_end"] == pytest.approx(2.5e5)

    # Isotherm switch to the paired kind keeps every parameter.
    ads = sim.adsorbents[0]
    comp = sim.components[0]
    sim.edit.set_isotherm_kind(ads, comp, "dual_site_langmuir_p")
    sim[f"adsorbent[{ads.id}].isotherm[{comp}].m1"] = 4.25
    sim.edit.set_isotherm_kind(ads, comp, "extended_dual_site_langmuir_p")
    assert sim[f"adsorbent[{ads.id}].isotherm[{comp}].m1"] == pytest.approx(4.25)

    sim.edit.set_mass_transfer_kind(ads, comp, "ldf", "macropore")
    sim[f"adsorbent[{ads.id}].kinetics[{comp}].tortuosity"] = 4.5
    sim.edit.set_mass_transfer_kind(ads, comp, "ldf", "haynes")
    assert sim[f"adsorbent[{ads.id}].kinetics[{comp}].tortuosity"] == pytest.approx(4.5)

    layer = sim.layers[0]
    sim.edit.set_dispersion_kind(layer, "ruthven")
    sim.edit.set_boundary_kind(step, "top", "closed")
    sim.edit.ensure_ids()
    assert sim.edit.repair_connections() == 0
    sim.edit.fill_gaps()
    ok, msg = sim.check()
    assert ok, msg


def _deck(sim):
    d = json.loads(sim.to_json())["simulation"]
    d.pop("System_Information")  # the save counter moves on every save
    return d


def test_refused_edits_change_nothing(perez):
    sim = perez
    before = _deck(sim)
    with pytest.raises(ValueError, match="no step"):
        sim.edit.remove_step(9999)
    with pytest.raises(ValueError):
        sim.edit.set_isotherm_kind(sim.adsorbents[0], sim.components[0], "no_such_kind")
    with pytest.raises(ValueError, match="already the first"):
        sim.edit.move_step(sim.steps[0], -1)
    assert _deck(sim) == before


# ---- every cascade leaves a deck that checks and round-trips ---------------------------------
# The kinds are read from model_kinds(), never listed here, so a new kind is covered at once.


def _paths(sim):
    return {p.path: p.text for p in sim.parameters()}


def _assert_sound(sim, model_mismatch=None):
    """The deck passes check() and a save/load round trip keeps every registry value and the
    check's verdict. ``model_mismatch`` names a text the check may report instead: a kind that
    does not suit the bed (Taylor-Aris dispersion in a packed bed) is the user's modelling
    choice, which the check reports, not a broken deck."""
    ok, msg = sim.check()
    assert ok or (model_mismatch and model_mismatch in msg), msg
    again = skarstrom.Simulation(sim.to_json())
    assert _paths(again) == _paths(sim)
    assert again.check() == (ok, msg)


def _tags(family):
    return [k["tag"] for k in skarstrom.model_kinds()[family]]


@pytest.mark.parametrize("family, setter", [
    ("dispersion", "set_dispersion_kind"),
    ("heat_dispersion", "set_heat_dispersion_kind"),
    ("internal_heat_transfer", "set_internal_heat_transfer_kind"),
    ("external_heat_transfer", "set_external_heat_transfer_kind"),
])
def test_layer_kind_switches_round_trip(perez, family, setter):
    layer = perez.layers[0]
    for tag in _tags(family):
        getattr(perez.edit, setter)(layer, tag)
        _assert_sound(perez, model_mismatch="packed bed")


def test_particle_heat_transfer_kind_switches_round_trip(perez):
    for tag in _tags("particle_heat_transfer"):
        perez.edit.set_particle_heat_transfer_kind(tag)
        _assert_sound(perez)


def test_every_isotherm_kind_round_trips(perez):
    ads, comp = perez.adsorbents[0], perez.components[0]
    for tag in _tags("isotherm"):
        perez.edit.set_isotherm_kind(ads, comp, tag)
        perez.edit.fill_gaps()
        _assert_sound(perez)


def test_every_mass_transfer_kind_round_trips(perez):
    ads, comp = perez.adsorbents[0], perez.components[0]
    for model in _tags("rate_model"):
        for coeff in _tags("rate_coefficient"):
            perez.edit.set_mass_transfer_kind(ads, comp, model, coeff)
            perez.edit.fill_gaps()
            _assert_sound(perez)


def test_boundary_stream_flow_pump_and_profile_switches_round_trip(perez):
    step = perez.edit.add_step("Edit")
    for tag in _tags("boundary"):
        perez.edit.set_boundary_kind(step, "top", tag)
        _assert_sound(perez)
    # A flow end holds streams: every device kind, and every flow quantity of a set-flow stream.
    perez.edit.set_boundary_kind(step, "bottom", "flow")
    for tag in _tags("stream"):
        endpoint = perez.edit.add_stream(step, "bottom", tag)
        _assert_sound(perez)
        if tag == "flow":
            for flow in _tags("flow"):
                perez.edit.set_flow_kind(endpoint, flow)
                _assert_sound(perez)
        if tag == "pump":
            for pump in _tags("pump_power"):
                perez.edit.set_pump_power_kind(endpoint, pump)
                _assert_sound(perez)
    for index in reversed(range(3)):
        perez.edit.remove_stream(step, "bottom", index)
        _assert_sound(perez)
    # A pressure end: every profile kind and every pump kind.
    perez.edit.set_boundary_kind(step, "top", "pressure")
    top = perez.steps.by_id(step).top.id
    for tag in _tags("profile"):
        perez.edit.set_profile_kind(top, "pressure", tag)
        _assert_sound(perez)
    for tag in _tags("pump"):
        perez.edit.set_pump_kind(top, tag)
        _assert_sound(perez)


def test_add_remove_connect_cascades_round_trip(perez):
    before_layers, before_steps = len(perez.layers), len(perez.steps)
    layer = perez.edit.add_layer("Second")
    _assert_sound(perez)
    a = perez.edit.add_step("A")
    b = perez.edit.duplicate_step(a)
    _assert_sound(perez)
    perez.edit.move_step(b, -1)
    perez.edit.move_layer(layer, -1)
    _assert_sound(perez)
    # Connect two flow ends, then disconnect by removing a stream (its link goes with it).
    perez.edit.set_boundary_kind(a, "top", "flow")
    perez.edit.set_boundary_kind(b, "top", "flow")
    s1 = perez.edit.add_stream(a, "top", "flow")
    s2 = perez.edit.add_stream(b, "top", "flow")
    perez.edit.connect(s1, s2)
    _assert_sound(perez)
    perez.edit.remove_stream(a, "top", 0)
    assert perez.edit.repair_connections() == 0
    _assert_sound(perez)
    perez.edit.remove_step(a)
    perez.edit.remove_step(b)
    perez.edit.remove_layer(layer)
    assert (len(perez.layers), len(perez.steps)) == (before_layers, before_steps)
    perez.edit.sync_kpi_config()
    _assert_sound(perez)


@pytest.mark.parametrize("example", [pytest.param(e, id=e.name) for e in EX.ALL])
def test_housekeeping_cascades_on_every_example(example):
    sim = example.load()
    before = _paths(sim)
    sim.edit.ensure_ids()
    assert sim.edit.repair_connections() == 0
    assert _paths(sim) == before
    # fill_gaps only adds entries a deck lacks (older decks), never changes one, and is idempotent.
    sim.edit.fill_gaps()
    filled = _paths(sim)
    assert {k: v for k, v in filled.items() if k in before} == before
    sim.edit.fill_gaps()
    assert _paths(sim) == filled
    ok, msg = sim.check()
    assert ok, msg
