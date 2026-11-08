"""The APDES-NFC library entry carries the water and humid CO2 parameters of
Stampi-Bombelli et al. (2020, corrected 2022).

The reference numbers are derived with the engine at 20 degC, p_CO2 = 40 Pa and
a relative humidity of 0.494 (p_sat(20 degC) = 2339.2 Pa); the paper does not
print them.
"""
import json
from pathlib import Path

import pytest

import skarstrom as sk

T = 293.15
P_WATER = 0.494 * 2339.2
DRY = dict(q_sat=2.38, x=0, T_0=296, b0=0.0707, dH0=57047, w0=0.4148, alpha=-1.606)
HUMID = dict(DRY, gamma=0.0016, beta=59.1)
LIB_DIR = Path(__file__).resolve().parents[2] / "data" / "adsorbents"
GAB = dict(q_m0=36.48, dH_m=0, C0=0.1489, dH_C=0, K0=0.5751, dH_K=0)


def test_water_dry_and_humid_loadings():
    q_w = float(sk.evaluate_isotherm("gab", GAB, pressure=P_WATER, temperature=T)[0])
    assert q_w == pytest.approx(2.84, abs=0.01)
    dry = float(sk.evaluate_isotherm("toth_humid", HUMID, pressure=40, temperature=T)[0])
    humid = float(sk.evaluate_isotherm("toth_humid", HUMID, pressure=40, temperature=T, water_loading=q_w)[0])
    toth2 = float(sk.evaluate_isotherm("toth_2", DRY, pressure=40, temperature=T)[0])
    assert dry == pytest.approx(0.82, abs=0.005)
    assert dry == pytest.approx(toth2, rel=1e-12)  # q_w = 0 is the dry Toth
    assert humid == pytest.approx(2.07, abs=0.005)


def test_library_entry_has_water_and_humid_kinds(monkeypatch):
    # The repository library, not the one bundled with an installed wheel.
    monkeypatch.setattr(sk, "adsorbent_library_dir", lambda: str(LIB_DIR))
    catalogue = json.loads((LIB_DIR / "library.merged.json").read_text())
    entry = next(e for e in catalogue["entries"] if e["key"] == "apdes_nfc_stampi2020")
    assert entry["components"] == ["CarbonDioxide", "Water"]
    assert entry["heat_of_adsorption"]["Water"]["value"] == 49000.0
    sim = sk.Simulation.new(title="humid")
    sim.add_component("CarbonDioxide")
    sim.add_component("Water")
    adsorbent_id = sim.add_library_adsorbent("apdes_nfc_stampi2020")
    ads = sim.adsorbents.by_id(adsorbent_id)
    co2, water = ads.isotherm["CarbonDioxide"], ads.isotherm["Water"]
    assert co2.gamma == pytest.approx(0.0016) and co2.beta == pytest.approx(59.1)
    assert water.q_m0 == pytest.approx(36.48) and water.C0 == pytest.approx(0.1489)
    assert water.K0 == pytest.approx(0.5751)
    assert ads.kinetics["Water"].k == pytest.approx(2e-3)
