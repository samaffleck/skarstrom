"""What a stored dataset measures (its observable binding): set and read it,
save it with the deck, extract the simulated counterpart and compare.

Runs are cut to 20 cells (Ward 2022, one breakthrough) and 3 cycles (Perez
2019 A): they check the plumbing, not the decks' results.
"""

import numpy as np
import pytest

import skarstrom as sk
from skarstrom import ureg as u


@pytest.fixture(scope="module")
def ward():
    sim = sk.examples.Ward2022_Q100ccm.load()
    sim["layer[1].cells"] = 20
    assert sim.run().ok
    return sim


def test_observable_infers_its_kind():
    assert sk.Observable("temperature").kind == "bulk"
    assert sk.Observable("mole_fraction", "CarbonDioxide").kind == "component"
    assert sk.Observable("purity", "CarbonDioxide", cycle_offset=-1).kind == "kpi"
    o = sk.Observable("temperature", location=7.15 * u.cm)
    assert o.location == "position" and o.z == pytest.approx(0.0715)
    assert sk.Observable.from_dict(o.to_dict()) == o
    with pytest.raises(ValueError):
        sk.Observable("temperatures")
    with pytest.raises(ValueError):
        sk.Observable("temperature", location="middle")
    with pytest.raises(AttributeError):
        o.z = 1.0


def test_bundled_ward_deck_is_bound():
    data = sk.examples.Ward2022_Q100ccm.load().experimental_data
    by_name = {s.name.split("_Qin")[0]: s.observable for s in data}
    assert by_name["outlet_CO2_mole_fraction[-]"] == sk.Observable("mole_fraction", "CarbonDioxide")
    assert by_name["outlet_volumetric_flowrate"] == sk.Observable("volumetric_flow_rate")
    assert by_name["temperature[K]_z=0.0715[m]"] == sk.Observable("temperature", location=0.0715)


def test_binding_is_checked_written_back_and_saved():
    sim = sk.examples.Ward2022_Q100ccm.load()
    data = sim.experimental_data
    s = data.add("T", [0, 1, 2], [20.0, 21.0, 22.0], observable=sk.Observable("temperature", location="bottom"),
                 x_unit="min", y_unit="degC", uncertainty=[0.5, 0.5, 0.5], weight=2)
    assert s.x_unit == "min" and s.y_unit == "degC" and s.si_unit == "K" and s.weight == 2.0
    assert s.uncertainty.tolist() == [0.5, 0.5, 0.5]
    with pytest.raises(ValueError, match="no component"):
        data.add("bad", [0, 1], [0, 1], observable=sk.Observable("mole_fraction", "Nitrogen"))
    with pytest.raises(ValueError, match="y unit"):
        data.add("bad", [0, 1], [0, 1], observable=sk.Observable("temperature"), y_unit="bar")
    with pytest.raises(ValueError, match="uncertainty"):
        data.add("bad", [0, 1], [0, 1], observable=sk.Observable("temperature"), uncertainty=[1.0])
    assert "bad" not in data

    s.weight = 3.0                                                     # written back
    with pytest.raises(ValueError, match="y unit"):                   # y_unit is still degC
        s.observable = sk.Observable("mole_fraction", "CarbonDioxide")
    s.bind(sk.Observable("mole_fraction", "CO2", location="top"), y_unit="%")  # formula -> deck name
    assert data["T"].weight == 3.0
    assert data["T"].observable == sk.Observable("mole_fraction", "CarbonDioxide")
    with pytest.raises(ValueError):
        s.y_unit = "degC"                                              # a mole fraction is not a temperature
    assert s.y_unit == "%"
    s.uncertainty = None

    copy = sk.Simulation.from_json(sim.to_json())
    t = copy.experimental_data["T"]
    assert t.observable == sk.Observable("mole_fraction", "CarbonDioxide")
    assert (t.x_unit, t.y_unit, t.uncertainty, t.weight) == ("min", "%", None, 3.0)
    t.observable = None
    assert copy.experimental_data["T"].observable is None


def test_deck_view_and_edit_ops_carry_the_binding():
    sim = sk.examples.Ward2022_Q100ccm.load()
    n = len(sim.experimental_data)
    index = sim._sim.edit("add_experimental_data", ["P", [0.0, 1.0], [1.0, 1.1], "t", "P", {"kind": "bulk",
                                                    "variable": "pressure", "y_unit": "bar"}])
    assert index == n
    view = sim.deck_view()
    entry = view["experimental_data"][n]
    assert entry["name"] == "P" and entry["y"] == [1.0, 1.1]
    assert entry["observable"]["kind"] == "bulk" and entry["observable"]["y_unit"] == "bar"
    assert entry["observable"]["si_unit"] == "Pa" and entry["observable"]["x_axis"] == "time"
    sim._sim.edit("set_experimental_data_binding", [n, {"kind": "none"}])
    assert sim.deck_view()["experimental_data"][n]["observable"]["kind"] == "none"
    with pytest.raises(ValueError):
        sim._sim.edit("set_experimental_data_binding", [n, {"kind": "bulk", "variable": "nope"}])
    sim._sim.edit("remove_experimental_data", [n])
    assert len(sim.experimental_data) == n


def test_simulated_and_compare_with_a_run(ward):
    outlet = next(s for s in ward.experimental_data if s.name.startswith("outlet_CO2"))
    t, y = outlet.simulated(ward)
    ref = ward.series("mole_fraction", "CarbonDioxide", location="top")
    np.testing.assert_array_equal(t, ref.time)
    np.testing.assert_array_equal(y, ref.values)

    fit = outlet.compare(ward)
    assert set(fit) == {"rmse", "nrmse", "max_abs_error", "points", "weighted_sse"}
    assert 0 < fit["points"] <= len(outlet.x) and fit["nrmse"] < 0.3
    # The engine and the (x, y) form agree on the same curve.
    legacy = outlet.compare(ref.time, ref.values)
    assert legacy["points"] == fit["points"]
    assert legacy["rmse"] == pytest.approx(fit["rmse"], rel=1e-12)

    every = ward.experimental_data.compare(ward)
    assert len(every) == 3 and all(r["points"] > 0 for r in every.values())

    # Units: the same thermocouple data in minutes and degC.
    tc = next(s for s in ward.experimental_data if s.name.startswith("temperature"))
    c = ward.experimental_data.add("tc (degC, min)", tc.x / 60.0, tc.y - 273.15, observable=tc.observable,
                                   x_unit="min", y_unit="degC")
    assert c.compare(ward)["rmse"] == pytest.approx(tc.compare(ward)["rmse"], rel=1e-9)
    t_c, y_c = c.simulated(ward)
    t_k, y_k = tc.simulated(ward)
    np.testing.assert_allclose(t_c, t_k / 60.0)
    np.testing.assert_allclose(y_c, y_k - 273.15)

    unbound = ward.experimental_data.add("free", [0, 1], [0, 1])
    with pytest.raises(ValueError, match="not bound"):
        unbound.compare(ward)


def test_per_cycle_kpi_dataset():
    sim = sk.examples.Perez2019_A.load()
    sim["layer[1].cells"] = 20
    assert sim.run(cycles=3).ok
    purity = [k.purity for k in sim.kpis() if k.component == "CarbonDioxide"]
    s = sim.experimental_data.add("purity", [1, 2, 3], [100 * p for p in purity],
                                  observable=sk.Observable("purity", "CarbonDioxide", cycle_offset=-1), y_unit="%")
    assert s.x_unit == "cycle"
    x, y = s.simulated(sim)
    assert x.tolist() == [1, 2, 3]
    np.testing.assert_allclose(y, [100 * p for p in purity])
    fit = s.compare(sim)
    assert fit["points"] == 3 and fit["rmse"] == pytest.approx(0, abs=1e-9)
