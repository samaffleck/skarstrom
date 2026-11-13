"""Every older rate-coefficient tag is still accepted, as a preset of the generic resistance model.

``set_rate_model`` and the ``set_mass_transfer_kind`` edit take each tag; the
reported coefficient is the tag given, the deck round-trips through JSON, and
the old parameter names keep working.
"""

import pytest

import skarstrom as sk

TAGS = ["constant", "macropore", "macropore_knudsen", "haynes", "arrhenius", "macropore_capacity",
        "series_resistance"]


def _sim():
    sim = sk.Simulation.new(title="presets")
    sim.add_component("CarbonDioxide")
    sim.add_component("Nitrogen")
    ads = sim.add_adsorbent("13X", solid_density=1138, heat_capacity=900, thermal_conductivity=0.63,
                            particle_void_fraction=0.35)
    return sim, ads


def test_one_kind_and_the_old_tags_as_presets():
    assert sk.build_kinds()["rate_coefficient"] == ["series_resistance"]
    assert sorted(sk.build_kinds()["rate_coefficient_preset"]) == sorted(TAGS[:-1])
    kinds = sk.model_kinds()["rate_coefficient"]
    assert [k["tag"] for k in kinds] == TAGS
    assert [k["legacy"] for k in kinds] == [True] * 6 + [False]
    generic = kinds[-1]
    assert [r["tag"] for r in generic["resistances"]] == ["film", "macropore", "molecular", "knudsen", "micropore", "lumped"]
    assert [m["tag"] for m in generic["capacity_modes"]] == ["none", "slope", "secant", "fixed"]
    assert {p["tag"] for p in generic["presets"]} == set(TAGS[:-1])


@pytest.mark.parametrize("rate", ["ldf", "ldf2", "dual_rate_ldf"])
@pytest.mark.parametrize("tag", TAGS)
def test_set_rate_model_accepts_each_tag(tag, rate):
    sim, ads = _sim()
    ads.set_rate_model("CarbonDioxide", rate, tag)
    path = f"adsorbent[{ads._id}].kinetics[CarbonDioxide]."
    names = {p.path[len(path):] for p in sim.parameters() if p.path.startswith(path)}
    # The preset's parameters are listed (pre_exp_factor is the lumped k); a
    # new series_resistance model has the lumped term only.
    expected = (["k", "resistances.lumped"] if tag == "series_resistance" else
                [{"pre_exp_factor": "k"}.get(p["name"], p["name"])
                 for p in sk.model_kinds()["rate_coefficient"][TAGS.index(tag)]["parameters"]])
    for name in expected:
        assert name in names, (tag, name, sorted(names))
    back = sk.Simulation(sim.to_json())
    assert _values(back) == _values(sim)
    assert _values(sk.Simulation(back.to_json())) == _values(back)


def _values(sim):
    return {p.path: sim.get_text(p.path) for p in sim.parameters() if not p.path.startswith("info.")}


def test_old_parameter_names_keep_working():
    sim, ads = _sim()
    ads.set_rate_model("CarbonDioxide", "ldf", "constant", k=0.3)
    ads.set_rate_model("Nitrogen", "ldf", "macropore", tortuosity=2.5, macropore_void_fraction=0.3)
    a = f"adsorbent[{ads._id}]."
    assert sim.get(a + "kinetics[CarbonDioxide].k") == pytest.approx(0.3)
    assert sim.get(a + "kinetics[Nitrogen].tortuosity") == pytest.approx(2.5)
    ads.set_rate_model("CarbonDioxide", "ldf", "arrhenius", pre_exp_factor=50.0, activation_energy=2e4)
    assert sim.get(a + "kinetics[CarbonDioxide].pre_exp_factor") == pytest.approx(50.0)
    assert sim.get(a + "kinetics[CarbonDioxide].k") == pytest.approx(50.0)
    # series_resistance: the switches are explicit; a parameter of a
    # switched-off term is not there.
    with pytest.raises(AttributeError):
        ads.set_rate_model("CarbonDioxide", "ldf", "series_resistance", crystal_radius=1e-6)
    ads.set_rate_model("CarbonDioxide", "ldf", "series_resistance", resistances=["micropore"],
                       crystal_radius=1e-6, micropore_diffusivity=2e-12)
    assert sim.get(a + "kinetics[CarbonDioxide].crystal_radius") == pytest.approx(1e-6)
    assert sim.get(a + "kinetics[CarbonDioxide].resistances.micropore") is True
    assert sim.get(a + "kinetics[CarbonDioxide].resistances.lumped") is False
    assert '"Micropore": true' in sim.to_json()


def test_switches_by_helper_and_path():
    sim, ads = _sim()
    ads.set_rate_model("CarbonDioxide", "dual_rate_ldf", "haynes", tortuosity=3.0)
    kin = sim.adsorbents["13X"].kinetics["CarbonDioxide"]
    assert kin.resistances.film is True and kin.resistances.lumped is False
    ads.set_rate_resistance("CarbonDioxide", "micropore", True, site="slow")
    kin = sim.adsorbents["13X"].kinetics["CarbonDioxide"]
    assert kin.slow_resistances.micropore is True and kin.resistances.micropore is False
    assert "slow_crystal_radius" in kin.fields and "crystal_radius" not in kin.fields
    ads.set_rate_resistance("CarbonDioxide", "film", False)
    with pytest.raises(ValueError):  # the last term of the fast site
        ads.set_rate_resistance("CarbonDioxide", "macropore", False)
    kin = sim.adsorbents["13X"].kinetics["CarbonDioxide"]
    kin.capacity_factor_mode = "slope"
    assert sim.adsorbents["13X"].kinetics["CarbonDioxide"].capacity_factor_mode == "slope"
    comp = sim.deck_view()["adsorbents"][0]["components"][0]
    assert comp["rate_coefficient"] == "series_resistance"
    assert comp["preset"] == "custom" and comp["slow_preset"] == "custom"
    assert comp["equation"].startswith("\\begin{aligned}\\frac{1}{k_i} &= ")


def test_to_script_reproduces_custom_switches():
    sim, ads = _sim()
    ads.set_rate_model("CarbonDioxide", "dual_rate_ldf", "series_resistance",
                       resistances=["film", "macropore", "micropore"], capacity_factor_mode="slope",
                       tortuosity=2.0, crystal_radius=2e-6, split=0.6)
    ads.set_rate_resistance("CarbonDioxide", "film", False, site="slow")
    sim.edit.set_capacity_factor_mode(ads.id, "CarbonDioxide", "fixed", site="slow")
    sim.adsorbents["13X"].kinetics["CarbonDioxide"].slow_capacity_factor = 25.0
    ns = {}
    exec(compile(sk.to_script(sim), "<script>", "exec"), ns)
    assert _values(ns["sim"]) == _values(sim)
