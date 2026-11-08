"""Compositions of mixed layers through the Python API, the registry and the CLI."""

import json
import subprocess
import sys

import pytest

import skarstrom as sk


@pytest.fixture()
def sim():
    return sk.examples.Illustrative_Diluted_Bed.load()


def test_composition_reads_back_by_name(sim):
    comp = sim.layers[0].composition
    assert comp["Activated Carbon"] == pytest.approx(0.7)
    assert comp["Glass beads"] == pytest.approx(0.3)


def test_set_composition_and_registry_path(sim):
    sim.layers[0].set_composition({"Activated Carbon": 0.5, "Glass beads": 0.5})
    assert sim.get("layer[1].mass_fraction[2]") == pytest.approx(0.5)
    sim["layer[1].mass_fraction[2]"] = 0.25
    sim["layer[1].mass_fraction[1]"] = 0.75
    assert sim.layers[0].composition["Glass beads"] == pytest.approx(0.25)


def test_fractions_not_summing_to_one_are_refused_at_run(sim):
    sim["layer[1].mass_fraction[1]"] = 0.5
    result = sim.run(cycles=1)
    assert result.status == "failed" and "sum to" in result.message


def test_inert_solid_holds_nothing_and_aggregate_is_mass_weighted(sim):
    assert sim.run(cycles=1).ok
    assert [a["label"] for a in sim.series_adsorbents()] == ["Activated Carbon", "Glass beads"]
    carbon = sim.series("adsorbed_concentration", "CarbonDioxide", z=0.1, adsorbent="Activated Carbon").values[-1]
    beads = sim.series("adsorbed_concentration", "CarbonDioxide", z=0.1, adsorbent="Glass beads").values[-1]
    total = sim.series("adsorbed_concentration", "CarbonDioxide", z=0.1).values[-1]
    assert carbon > 0 and beads == pytest.approx(0, abs=1e-12)
    assert total == pytest.approx(0.7 * carbon, rel=1e-6)


def test_cli_set_composition_by_path(tmp_path):
    deck = sk.examples.Illustrative_Diluted_Bed.path
    base = [sys.executable, "-m", "skarstrom", "run", str(deck), "-q", "--json"]
    ok = subprocess.run(base + ["--set", "layer[1].mass_fraction[1]=0.5", "--set", "layer[1].mass_fraction[2]=0.5"],
                        capture_output=True, text=True)
    assert ok.returncode == 0, ok.stderr
    assert json.loads(ok.stdout)["status"] == "completed"
    bad = subprocess.run(base + ["--set", "layer[1].mass_fraction[1]=0.5"], capture_output=True, text=True)
    assert bad.returncode == 1
    assert json.loads(bad.stdout)["status"] == "failed"
