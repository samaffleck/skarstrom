"""Per-adsorbent series of mixed layers: ``series(..., adsorbent=)``, ``profile(..., adsorbent=)``,
``series_adsorbents()``, the series manifest and the CLI's CSV logs."""

import numpy as np
import pytest

import skarstrom as sk
from skarstrom.cli import main
from test_build import build_vsa

CO2 = "CarbonDioxide"


@pytest.fixture(scope="module")
def mixed():
    """The from-scratch VSA deck: one layer of 13X (w = 0.9) with inert beads (0.1)."""
    sim, _ = build_vsa()
    assert sim.run(cycles=1).ok
    return sim


def test_series_and_profile_per_adsorbent(mixed):
    ads = {a["label"]: a for a in mixed.series_adsorbents()}
    assert set(ads) == {"13X", "Beads"}
    assert all(a["logged"] for a in ads.values())

    agg = mixed.profile("adsorbed_concentration", CO2)
    q13x = mixed.profile("adsorbed_concentration", CO2, adsorbent="13X")
    beads = mixed.profile("adsorbed_concentration", CO2, adsorbent=ads["Beads"]["id"])
    assert q13x.adsorbent == "13X" and beads.adsorbent == "Beads"
    assert q13x.values.shape == agg.values.shape
    assert np.max(np.abs(q13x.values)) > 0 and np.all(beads.values == 0)
    np.testing.assert_allclose(agg.values, 0.9 * q13x.values, rtol=1e-12, atol=1e-15)

    # By adsorbent object too; series at a position.
    s = mixed.series("equilibrium_adsorbed_concentration", CO2, z=0.1, adsorbent=mixed.adsorbents["13X"])
    s_agg = mixed.series("equilibrium_adsorbed_concentration", CO2, z=0.1)
    assert s.adsorbent == "13X" and "(13X)" in repr(s)
    np.testing.assert_allclose(s_agg.values, 0.9 * s.values, rtol=1e-12, atol=1e-15)

    with pytest.raises(ValueError):
        mixed.series("pressure", adsorbent="13X")
    with pytest.raises(ValueError):
        mixed.series("mole_fraction", CO2, adsorbent="13X")
    with pytest.raises(ValueError, match="unknown adsorbent"):
        mixed.series("adsorbed_concentration", CO2, adsorbent="Alumina")


def test_single_solid_reads_the_aggregate():
    sim = sk.examples.Ward2022_Q100ccm.load()
    assert sim.run(cycles=1).ok
    (a,) = sim.series_adsorbents()
    assert not a["logged"]
    agg = sim.profile("adsorbed_concentration", sim.components[0])
    one = sim.profile("adsorbed_concentration", sim.components[0], adsorbent=a["label"])
    np.testing.assert_array_equal(one.values, agg.values)


def test_series_manifest_lists_the_per_adsorbent_variables():
    m = sk.series_manifest()
    keys = [v["key"] for v in m["variables"]]
    assert keys == sk.bulk_variables() + sk.component_variables()
    per = [v["key"] for v in m["variables"] if v["per_adsorbent"]]
    assert per == sk.per_adsorbent_variables()
    assert set(per) == {"adsorbed_concentration", "equilibrium_adsorbed_concentration",
                        "mass_transfer_coefficient", "slow_adsorbed_concentration"}
    v = next(v for v in m["variables"] if v["key"] == "adsorbed_concentration")
    assert v["adsorbent_name"] == "{component} Adsorbed Concentration ({adsorbent})"


def test_cli_csv_logs_have_the_per_adsorbent_series(capsys, tmp_path):
    sim, _ = build_vsa()
    deck = tmp_path / "mixed.json"
    sim.save(deck)
    out = tmp_path / "out"
    code = main(["run", str(deck), "--cycles", "1", "--csv", str(out), "-q"])
    _, err = capsys.readouterr()
    assert code == 0, err
    logs = out / "logs"
    assert (logs / f"{CO2}_Adsorbed_Concentration.csv").exists()
    for name in (f"{CO2}_Adsorbed_Concentration__13X_.csv", f"{CO2}_Adsorbed_Concentration__Beads_.csv",
                 f"{CO2}_Equilibrium_Adsorbed_Concentration__13X_.csv", f"{CO2}_Mass_Transfer_Coefficient__13X_.csv"):
        assert (logs / name).exists(), name
    with open(logs / f"{CO2}_Adsorbed_Concentration__13X_.csv") as f:
        header = f.readline().strip().split(",")
    assert header[:4] == ["Cycle Index", "Step Index", "time", "cell_0"]
