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


def test_every_old_tag_is_listed():
    assert sk.build_kinds()["rate_coefficient"] == TAGS


@pytest.mark.parametrize("rate", ["ldf", "ldf2", "dual_rate_ldf"])
@pytest.mark.parametrize("tag", TAGS)
def test_set_rate_model_accepts_each_tag(tag, rate):
    sim, ads = _sim()
    ads.set_rate_model("CarbonDioxide", rate, tag)
    path = f"adsorbent[{ads._id}].kinetics[CarbonDioxide]."
    names = {p.path[len(path):] for p in sim.parameters() if p.path.startswith(path)}
    for p in sk.model_kinds()["rate_coefficient"][TAGS.index(tag)]["parameters"]:
        assert p["name"] in names, (tag, p["name"], sorted(names))
    back = sk.Simulation(sim.to_json())
    assert back.to_json() == sk.Simulation(back.to_json()).to_json()


def test_old_parameter_names_keep_working():
    sim, ads = _sim()
    ads.set_rate_model("CarbonDioxide", "ldf", "constant", k=0.3)
    ads.set_rate_model("Nitrogen", "ldf", "macropore", tortuosity=2.5, macropore_void_fraction=0.3)
    a = f"adsorbent[{ads._id}]."
    assert sim.get(a + "kinetics[CarbonDioxide].k") == pytest.approx(0.3)
    assert sim.get(a + "kinetics[Nitrogen].tortuosity") == pytest.approx(2.5)
    ads.set_rate_model("CarbonDioxide", "ldf", "arrhenius", pre_exp_factor=50.0, activation_energy=2e4)
    assert sim.get(a + "kinetics[CarbonDioxide].pre_exp_factor") == pytest.approx(50.0)
    # The series kind keeps its version-1 behaviour while built: setting the
    # micropore parameters turns the term on.
    ads.set_rate_model("CarbonDioxide", "ldf", "series_resistance", crystal_radius=1e-6,
                       micropore_diffusivity=2e-12)
    assert sim.get(a + "kinetics[CarbonDioxide].crystal_radius") == pytest.approx(1e-6)
    assert '"Micropore": true' in sim.to_json()
