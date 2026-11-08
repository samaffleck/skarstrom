"""Reaction KPIs (``sim.kpi.track_reaction``, the ``conversion`` ... fields of KPIRecord) and the
reaction-rate series (``series(..., reaction=)``, ``profile(..., reaction=)``), through the
Python API and the CLI."""

import csv
import json
import math

import numpy as np
import pytest

import skarstrom as sk
from skarstrom.cli import KPI_CSV_FIELDS, main

N2, BUT, IBUT = "Nitrogen", "Butane", "Isobutane"


def build_reactor(reaction=True, cells=20, step_time=5.0):
    """A packed catalyst bed fed with 10 % butane; first-order isomerisation butane -> isobutane."""
    sim = sk.Simulation.new(title="isomerisation")
    for c in (N2, BUT, IBUT):
        sim.add_component(c)
    ads = sim.add_adsorbent("Catalyst", solid_density=1000, heat_capacity=900, thermal_conductivity=0.5,
                            particle_void_fraction=0.4)
    sim.bed.set_shape("cylinder", diameter=0.05)
    sim.wall.set_material("ss316l")
    sim.wall.set_boundary("zero_flux")
    sim.add_layer("Bed", composition={"Catalyst": 1.0}, momentum="kozeny_carman", length=0.2, cells=cells,
                  bed_void_fraction=0.4, particle_diameter=5e-3)
    sim.initial_conditions.pressure = 1e5
    sim.initial_conditions.temperature = 400
    sim.initial_conditions.y[N2] = 1.0
    feed = sim.add_step("Feed", time=step_time)
    feed.bottom.set_streams("inlet", [sk.Stream("Feed inlet", flow=sk.velocity(0.1), y={N2: 0.9, BUT: 0.1},
                                                temperature=400)])
    feed.top.set_pressure("outlet", 1e5)
    if reaction:
        sim._build("add_reaction", ads.id, "R1", {BUT: -1.0, IBUT: 1.0}, "power_law")
        sim._build("set_reaction_basis", ads.id, 1, "concentration")
        # Da = rho_b k L / u = 1
        sim["adsorbent[Catalyst].reaction[R1].rate.pre_exp_factor"] = 0.1 / (600 * 0.2)
    sim["energy_balance.energy_balance_mode"] = "isothermal"
    sim["css.method"] = 0
    sim.kpi.track(BUT, feed=[feed.bottom], product=[feed.top])
    return sim


@pytest.fixture(scope="module")
def reactor():
    sim = build_reactor()
    sim.kpi.track_reaction(BUT, products=[IBUT])
    assert sim.run(cycles=2).ok
    return sim


def test_reaction_kpi_fields(reactor):
    rows = {k.component: k for k in reactor.kpis(cycle=1)}
    assert set(rows) == {BUT, IBUT}  # isobutane: a product, not tracked
    a, b = rows[BUT], rows[IBUT]
    assert a.reaction_kpis and b.reaction_kpis
    assert 0.55 < a.conversion < 0.65  # 1 - exp(-1) = 0.632, upwind on 20 cells 0.623
    assert math.isnan(a.selectivity) and math.isnan(b.conversion)
    assert b.selectivity == pytest.approx(1.0, rel=1e-9)
    assert b.reaction_yield == pytest.approx(a.conversion, rel=1e-9)
    assert b.reaction_generation == pytest.approx(-a.reaction_generation, rel=1e-9)
    assert b.space_time_yield > 0 and b.space_time_yield_unit == "kg/(kg s)"
    assert a.peak_temperature == pytest.approx(400.0, abs=1e-6)  # isothermal
    assert not b.valid and math.isnan(b.purity)
    d = b.as_dict()
    assert {"conversion", "selectivity", "reaction_yield", "space_time_yield", "peak_temperature"} <= set(d)
    # Every emitted KPI field is in the manifest (the flag and the unit are not KPIs).
    keys = {k["key"] for k in sk.kpi_manifest()}
    assert set(d) - {"cycle", "component", "valid", "reaction_kpis", "space_time_yield_unit"} <= keys
    reaction = [k["key"] for k in sk.kpi_manifest() if k["reaction"]]
    assert {"conversion", "selectivity", "reaction_yield", "space_time_yield"} <= set(reaction)


def test_kpi_table_has_the_reaction_columns(reactor):
    pytest.importorskip("pandas")
    table = reactor.kpi_table(cycle=1)
    assert "conversion" in table.columns and "selectivity" in table.columns


def test_registry_and_basis(reactor):
    info = next(p for p in reactor.parameters(search="kpi.reaction.key_reactant"))
    assert info.choices[reactor.get("kpi.reaction.key_reactant")] == BUT
    assert reactor.get("kpi.reaction.product[Isobutane]")
    assert not reactor.get("kpi.reaction.product[Nitrogen]")
    reactor["kpi.reaction.basis"] = "bed_volume"
    try:
        b = next(k for k in reactor.kpis(cycle=1) if k.component == IBUT)
        assert b.space_time_yield_unit == "kg/(m3 s)"
        volume = math.pi * 0.05 ** 2 / 4 * 0.2
        mass = 600 * volume
        b_mass = b.space_time_yield * volume / mass  # the same generation per kg of catalyst
        reactor["kpi.reaction.basis"] = "catalyst_mass"
        b2 = next(k for k in reactor.kpis(cycle=1) if k.component == IBUT)
        assert b2.space_time_yield == pytest.approx(b_mass, rel=1e-9)
    finally:
        reactor["kpi.reaction.basis"] = "catalyst_mass"
    with pytest.raises(ValueError):
        reactor["kpi.reaction.key_reactant"] = "Unobtainium"
    with pytest.raises(ValueError):
        reactor.kpi.track_reaction(BUT, products=[IBUT], basis="per_hour")
    assert "track_reaction('Butane', products=['Isobutane'])" in sk.to_script(reactor)


def test_reaction_rate_series_and_profile(reactor):
    (r,) = reactor.series_reactions()
    assert r["name"] == "R1" and r["adsorbent"] == "Catalyst" and r["variable_name"] == "R1 Rate (Catalyst)"
    assert sk.reaction_variables() == ["reaction_rate"]
    p = reactor.profile("reaction_rate", reaction="R1")
    assert p.reaction == "R1" and p.unit == "mol/(kg.s)" and p.values.shape[1] == 20
    by_id = reactor.profile("reaction_rate", reaction=r["id"], adsorbent="Catalyst")
    np.testing.assert_array_equal(by_id.values, p.values)
    top = reactor.series("reaction_rate", reaction="R1")
    bottom = reactor.series("reaction_rate", reaction="R1", location="bottom")
    assert top.reaction == "R1" and "R1" in repr(top)
    assert bottom.values[-1] > top.values[-1] > 0
    assert top.values[-1] == p.values[-1, -1]  # the last packed cell
    with pytest.raises(ValueError):
        reactor.series("reaction_rate")
    with pytest.raises(ValueError):
        reactor.series("pressure", reaction="R1")
    with pytest.raises(ValueError, match="unknown reaction"):
        reactor.profile("reaction_rate", reaction="R9")


def test_cli_reaction_columns(capsys, tmp_path, reactor):
    deck = tmp_path / "reactor.json"
    reactor.save(deck)
    out = tmp_path / "out"
    code = main(["run", str(deck), "--cycles", "1", "--csv", str(out), "--json", "-q"])
    stdout, err = capsys.readouterr()
    assert code == 0, err
    payload = json.loads(stdout)
    a = next(k for k in payload["kpis"] if k["component"] == BUT)
    assert a["reaction_kpis"] and 0 < a["conversion"] < 1 and a["selectivity"] is None
    with open(out / "kpis.csv") as f:
        header = next(csv.reader(f))
    assert "conversion [mol/mol]" in header and "space-time yield [kg/(kg s)]" in header
    assert (out / "logs" / "R1_Rate__Catalyst_.csv").exists()


def test_no_reactions_keep_their_kpis_and_csv(capsys, tmp_path):
    # A deck without reactions: configuring a key reactant adds no record, field or column.
    sim = build_reactor(reaction=False, cells=10, step_time=2.0)
    files = []
    for key in (None, BUT):
        if key:
            sim.kpi.track_reaction(key, products=[IBUT])
        deck = tmp_path / f"plain_{key}.json"
        sim.save(deck)
        out = tmp_path / f"out_{key}"
        code = main(["run", str(deck), "--cycles", "1", "--csv", str(out), "--json", "-q"])
        stdout, err = capsys.readouterr()
        assert code == 0, err
        payload = json.loads(stdout)
        assert [k["component"] for k in payload["kpis"]] == [BUT]
        assert "conversion" not in payload["kpis"][0]
        with open(out / "kpis.csv") as f:
            rows = list(csv.reader(f))
        assert rows[0] == ["cycle", "component", "valid"] + [h for _, h in KPI_CSV_FIELDS]
        files.append(rows)
    assert files[0] == files[1]
    assert sim.run(cycles=1).ok
    (k,) = sim.kpis(cycle=-1)
    assert not k.reaction_kpis and math.isnan(k.conversion)
    assert sim.series_reactions() == []
