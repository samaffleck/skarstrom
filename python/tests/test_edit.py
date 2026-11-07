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
        assert build[family] == [k["tag"] for k in kinds[family]], family


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
