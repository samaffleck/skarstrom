"""Tests for the skarstrom Python package (run against the installed wheel)."""

import json
import math
import os

import numpy as np
import pytest

import skarstrom

PEREZ = skarstrom.examples.Perez2019_A      # four-step CO2 VSA, CO2 tracked: ~0.2 s per cycle
WARD = skarstrom.examples.Ward2022_Q100ccm  # one-step breakthrough: < 1 s


@pytest.fixture
def perez():
    sim = skarstrom.Simulation(PEREZ)
    sim["cycle.number_of_cycles"] = 3
    return sim


@pytest.fixture(scope="module")
def perez_run():
    sim = skarstrom.Simulation(PEREZ)
    sim["cycle.number_of_cycles"] = 3
    # Upwind keeps the outlet mole fractions in [0, 1] to round-off, which
    # test_run_results checks (to 1e-9); the deck's WENO undershoots by about
    # the absolute tolerance (-2e-8 at the top end; docs/numerics.md).
    sim.settings.advection_scheme = "upwind"
    result = sim.run()
    assert result.ok, result.message
    return sim


def test_version_and_examples(tmp_path):
    assert skarstrom.__version__.count(".") == 2
    names = skarstrom.examples.names()
    assert PEREZ.name in names and WARD.name in names
    assert len(names) >= 30, "every deck in examples/decks/ ships with the package"
    with pytest.raises(KeyError):
        skarstrom.examples.get("nope")
    copy = WARD.copy_to(tmp_path)
    assert os.path.basename(copy) == "Ward2022_Q100ccm.json" and os.path.isfile(copy)


def test_adsorbent_library_ships_with_the_package():
    lib = skarstrom.adsorbent_library()
    assert lib["schema_version"] == 3
    keys = [e["key"] for e in lib["entries"]]
    assert "zeolite_13x_perez2019" in keys and len(keys) == len(set(keys))
    for e in lib["entries"]:
        assert os.path.isfile(skarstrom.adsorbent_path(e["key"])), e["key"]
        if e["status"] == "sourced":
            assert e["sources"], e["key"]
    with pytest.raises(KeyError):
        skarstrom.adsorbent_path("no_such_adsorbent")
    # A material id gives its recommended entry.
    for mid, mat in lib["materials"].items():
        rec = next(e for e in lib["entries"] if e["key"] == mat["recommended"])
        assert rec["material_id"] == mid
        assert skarstrom.adsorbent_path(mid) == skarstrom.adsorbent_path(rec["key"])


def test_structural_edits():
    sim = skarstrom.Simulation(PEREZ)
    sim.add_component("Argon")
    assert "Argon" in sim.components
    assert any("y[Argon]" in p.path for p in sim.parameters())
    with pytest.raises(ValueError):
        sim.add_component("Unobtainium")
    sim.remove_component("Argon")
    assert "Argon" not in sim.components
    def n_adsorbents():
        return len(json.loads(sim.to_json())["simulation"]["Adsorbent_Directory"]["Adsorbents"])
    before = n_adsorbents()
    ads = sim.add_library_adsorbent("zeolite_13x_perez2019")
    assert ads > 0 and n_adsorbents() == before + 1
    assert sim.add_library_adsorbent("zeolite-13x") > ads  # a material id: its recommended entry
    assert n_adsorbents() == before + 2
    with pytest.raises(ValueError):
        sim.add_library_adsorbent("no_such_key")
    sim["cycle.number_of_cycles"] = 1
    assert sim.run().ok


def test_licences_ship_with_the_package():
    d = os.path.join(os.path.dirname(skarstrom.__file__), "licenses")
    for f in ("LICENSE", "NOTICE", "THIRD_PARTY_NOTICES.md"):
        assert os.path.isfile(os.path.join(d, f)), f


def test_load_from_path_and_json(tmp_path):
    a = skarstrom.Simulation(WARD)
    b = skarstrom.Simulation(a.to_json())
    c = skarstrom.Simulation.from_file(WARD)
    assert a.components == b.components == c.components == ["CarbonDioxide", "Helium"]
    assert a.steps.names == b.steps.names == ["Feed"]
    a["settings.reltol"] = 2e-6
    out = tmp_path / "deck.json"
    a.save(out)
    assert skarstrom.Simulation(out)["settings.reltol"] == pytest.approx(2e-6)
    assert a.copy()["settings.reltol"] == pytest.approx(2e-6)
    with pytest.raises(RuntimeError):
        skarstrom.Simulation(str(tmp_path / "missing.json"))


def test_parameters_list_get_set(perez):
    params = perez.parameters()
    paths = [p.path for p in params]
    assert len(paths) == len(set(paths))
    assert "settings.reltol" in paths and "info.title" in paths and "step[1].time" in paths

    tunable = perez.parameters(tunable_only=True)
    assert tunable and all(p.tunable and p.bounds is not None for p in tunable)
    assert all(p.kind in ("real", "expression") for p in tunable)

    feed = perez.steps[0].name
    assert perez["step[1].time"] == perez[f"step[{feed}].time"]
    perez["step[1].time"] = 50
    assert perez["step[1].time"] == 50.0
    perez["step[1].time"] = "1 [min]"
    assert perez["step[1].time"] == pytest.approx(60.0)
    assert perez.get_text("step[1].time") == "1 [min]"
    with pytest.raises(ValueError):
        perez["step[1].time"] = "1 +* 2"
    assert perez["step[1].time"] == pytest.approx(60.0)

    assert isinstance(perez["cycle.number_of_cycles"], int)
    flag = "energy_balance.isothermal_use_initial_temperature"
    perez[flag] = False
    assert perez[flag] is False
    perez[flag] = "true"
    assert perez[flag] is True
    with pytest.raises(ValueError):
        perez[flag] = "maybe"
    perez["info.title"] = "from python"
    assert perez["info.title"] == "from python"
    perez.update({"css.min_cycles": 2, "settings.max_order": 4})
    assert perez["css.min_cycles"] == 2 and perez["settings.max_order"] == 4

    with pytest.raises(KeyError):
        perez["no.such.parameter"]
    with pytest.raises(KeyError):
        perez["no.such.parameter"] = 1.0
    with pytest.raises(TypeError):
        perez["step[1].time"] = [1, 2]

    found = perez.parameters(search="TIME")
    assert found and all("time" in (p.path + p.alias + p.name).lower() for p in found)


def test_run_results(perez_run):
    sim = perez_run
    r = sim.last_result
    assert r.status in ("completed", "reached_css")
    assert r.cycles_run == 3 and sim.cycles == [0, 1, 2] and sim.last_cycle == 2
    assert r.stats["steps"] > 0 and r.wall_time_s > 0

    kpis = sim.kpis()
    assert len(kpis) == 3
    last = sim.kpis(cycle=-1)
    assert len(last) == 1 and last[0].cycle == 2 and last[0].component == "CarbonDioxide"
    assert last[0].valid and 0 < last[0].purity <= 1 and 0 < last[0].recovery <= 1

    css = sim.css_history()
    assert [c.cycle for c in css] == [0, 1, 2]
    # Mass Balance (Per-Component): the balance plus the end-of-cycle state test,
    # on three consecutive cycles; three cycles are far from CSS.
    assert sim.css.state_tolerance == pytest.approx(2e-5) and sim.css.consecutive_cycles == 3
    assert not any(c.converged for c in css)
    assert all(c.required_cycles == 3 and c.state_tolerance == pytest.approx(2e-5) for c in css)
    assert math.isnan(css[0].state_change) and css[2].state_change > 2e-5
    # The remaining-distance diagnostic needs at least three decaying state
    # changes; three cycles give two.
    assert not any(c.distance_estimated for c in css)
    assert all(c.distance_note == "history too short" for c in css)
    assert all(math.isnan(c.distance) and math.isnan(c.distance_ratio) and math.isnan(c.convergence_rate)
               for c in css)

    y = sim.series("mole_fraction", "CarbonDioxide", location="top")
    assert isinstance(y.time, np.ndarray) and y.time.shape == y.values.shape
    assert y.time[0] == pytest.approx(0.0) and np.all(np.diff(y.time) >= 0)
    assert np.all((y.values >= -1e-9) & (y.values <= 1 + 1e-9))
    assert y.unit == "mol/mol"
    p = sim.series("pressure", z=0.1, cycle=0)
    assert p.unit == "Pa" and len(p) > 0

    prof = sim.profile("temperature", step=1)
    assert prof.values.shape == (prof.time.size, prof.z.size)
    assert prof.z.size == sim.positions.size and not prof.faces
    np.testing.assert_array_equal(prof.z, sim.positions)
    # Cell centres sit halfway between the faces (no ghost cells); faces span 0..L.
    zf = sim.face_positions
    assert zf.size == sim.positions.size + 1 and zf[0] == 0.0
    np.testing.assert_allclose(sim.positions, 0.5 * (zf[:-1] + zf[1:]))
    u = sim.profile("velocity", step=1)
    assert u.faces and u.values.shape == (u.time.size, zf.size)
    np.testing.assert_array_equal(u.z, zf)

    streams = sim.streams()
    assert streams
    assert len({s["endpoint_id"] for s in streams}) == len(streams)
    s = streams[0]
    f = sim.stream_series(s["step"], s["name"], "molar_flow_rate")
    assert f.time.shape == f.values.shape
    g = sim.endpoint_series(s["endpoint_id"], "molar_flow_rate")
    np.testing.assert_array_equal(f.values, g.values)
    assert "pressure" in skarstrom.stream_variables()

    # The bed ends are the endpoint logs: for a step with one endpoint at an
    # end, series(location=...) over that step equals stream_series.
    t_all = sim.series("pressure", location="bottom").time
    for s in streams:
        if sum(1 for o in streams if o["step"] == s["step"] and o["location"] == s["location"]) != 1:
            continue
        n = sim.stream_series(s["step"], s["name"], "pressure")
        rows = [sim.profile("pressure", step=k).time.size for k in range(len(sim.steps))]
        start = sum(rows[: s["step"]])
        for var, comp in (("pressure", None), ("temperature", None), ("mole_fraction", "CarbonDioxide")):
            bed = sim.series(var, comp, location=s["location"]).values[start : start + rows[s["step"]]]
            log = sim.stream_series(s["step"], s["name"], var, comp).values
            np.testing.assert_array_equal(bed, log, err_msg=f"{s['name']} {var}")
        assert n.unit == "Pa" and t_all.size == sum(rows)

    with pytest.raises(ValueError):
        sim.series("not_a_variable")
    with pytest.raises(ValueError):
        sim.series("mole_fraction")  # needs a component
    with pytest.raises(IndexError):
        sim.series("pressure", cycle=99)


def test_kpi_table(perez_run):
    pd = pytest.importorskip("pandas")
    df = perez_run.kpi_table()
    assert list(df["cycle"]) == [0, 1, 2]
    assert {"purity", "recovery", "energy_consumption", "electrical_energy"} <= set(df.columns)
    # The electrical energy is energy_consumption under its explicit name.
    assert df["electrical_energy"].equals(df["energy_consumption"])


def test_progress_cancel_and_messages(perez):
    calls = []

    def progress(p):
        calls.append(p)
        return len(calls) < 3

    messages = []
    r = perez.run(progress=progress, on_message=lambda lvl, msg: messages.append((lvl, msg)))
    assert r.status == "cancelled"
    assert len(calls) == 3
    assert calls[0]["total_steps"] == 4 and 0 < calls[0]["fraction"] <= 1
    assert messages and perez.messages


def test_progress_exception_propagates(perez):
    def boom(p):
        raise ZeroDivisionError("from callback")

    with pytest.raises(ZeroDivisionError):
        perez.run(progress=boom)


def test_cycles_override(perez):
    r = perez.run(cycles=1)
    assert r.ok and r.cycles_run == 1


def test_isotherm_evaluate_and_fit():
    kinds = skarstrom.isotherm_kinds()
    kind = "single_site_langmuir_p"
    assert kind in kinds
    names = skarstrom.isotherm_parameters(kind)
    assert {"m1", "b1", "q1"} <= set(names)

    truth = {"m1": 4.0, "b1": 2e-9, "q1": 25000.0}
    P = np.linspace(1e3, 5e5, 25)
    T = np.repeat([288.15, 308.15, 328.15], 25)
    PP = np.tile(P, 3)
    q = skarstrom.evaluate_isotherm(kind, truth, PP, T)
    # closed form: q = m1 b P / (1 + b P), b = b1 exp(q1 / (R T))
    R = 8.31446261815324
    b = truth["b1"] * np.exp(truth["q1"] / (R * T))
    np.testing.assert_allclose(q, truth["m1"] * b * PP / (1 + b * PP), rtol=1e-6)

    fit = skarstrom.fit_isotherm(kind, PP, q, T)
    assert fit["r2"] > 0.9999
    for k, v in truth.items():
        assert fit["params"][k] == pytest.approx(v, rel=1e-3)
    held = skarstrom.fit_isotherm(kind, PP, q, T, fixed={"q1": 25000.0})
    assert "q1" in held["fixed"] and held["params"]["q1"] == 25000.0
    with pytest.raises(ValueError):
        skarstrom.fit_isotherm(kind, PP, q, T, fixed={"bogus": 1.0})
    with pytest.raises(ValueError):
        skarstrom.evaluate_isotherm("not-a-kind", {}, [1e5], 298.15)


def _gab(qm, c, k, x):
    return qm * c * k * x / ((1 - k * x) * (1 - k * x + c * k * x))


def test_gab_isotherm_component_evaluate_fit_and_iast():
    """The GAB kinds work in relative pressure x = p / p_sat(T) of a named vapour (default water)."""
    kinds = skarstrom.isotherm_kinds()
    assert {"gab", "gab_young"} <= set(kinds)
    assert skarstrom.isotherm_parameters("gab") == ["q_m0", "dH_m", "C0", "dH_C", "K0", "dH_K"]

    truth = {"q_m0": 3.51, "dH_m": 0.0, "C0": 2.73, "dH_C": 0.0, "K0": 0.867, "dH_K": 0.0}
    ps = 3169.89  # water at 298.15 K [Pa] (CoolProp saturation ancillary)
    x = np.linspace(0.05, 0.85, 12)
    q = skarstrom.evaluate_isotherm("gab", truth, x * ps, 298.15)
    np.testing.assert_allclose(q, _gab(3.51, 2.73, 0.867, x), rtol=1e-4)
    np.testing.assert_allclose(skarstrom.evaluate_isotherm("gab", truth, x * ps, 298.15, component="Water"), q, rtol=0)

    # Another vapour: the same relative pressures of ethanol give the same loadings;
    # read against water's saturation pressure they do not.
    pe = 7870.0  # ethanol at 298.15 K, approximately [Pa]
    q_eth = skarstrom.evaluate_isotherm("gab", truth, x * pe, 298.15, component="Ethanol")
    np.testing.assert_allclose(q_eth, q, rtol=2e-2)
    assert not np.allclose(skarstrom.evaluate_isotherm("gab", truth, x * pe, 298.15), q_eth, rtol=0.1)

    fit = skarstrom.fit_isotherm("gab", x * ps, q, 298.15)
    assert fit["r2"] > 0.9999
    for k in ("q_m0", "C0", "K0"):
        assert fit["params"][k] == pytest.approx(truth[k], rel=1e-3)
    assert {"dH_m", "dH_C", "dH_K"} <= set(fit["fixed"])  # one temperature
    fit_eth = skarstrom.fit_isotherm("gab", x * pe, q_eth, 298.15, component="Ethanol")
    assert fit_eth["params"]["K0"] == pytest.approx(0.867, rel=1e-3)

    # IAST: water (GAB) with an inert component gets its pure loading.
    yw = 0.5 * ps / 1e5
    qi = skarstrom.evaluate_iast(["gab", "inert"], [truth, {}], pressure=1e5, temperature=298.15,
                                 mole_fraction=[yw, 1 - yw])
    assert qi[0] == pytest.approx(_gab(3.51, 2.73, 0.867, 0.5), rel=1e-4)
    assert qi[1] == 0.0
    qe = skarstrom.evaluate_iast(["gab", "inert"], [truth, {}], pressure=1e5, temperature=298.15,
                                 mole_fraction=[0.5 * pe / 1e5, 1 - 0.5 * pe / 1e5], components=["Ethanol", None])
    assert qe[0] == pytest.approx(qi[0], rel=2e-2)
    with pytest.raises(ValueError):
        skarstrom.evaluate_iast(["gab", "inert"], [truth, {}], 1e5, 298.15, [0.01, 0.99], components=["Water"])
