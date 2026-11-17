"""The SUPERSATURATION run-health code reaches Python unchanged.

The bundled humid DAC breakthrough runs at about 0.49 relative saturation
(1.15 % water at 293.15 K and 1 atm); raising its feed to 3 % water (1.3 times
p_sat at 293.15 K) takes it past saturation at the wall, which the model does not handle.
"""
import math

import skarstrom as sk


def _humid_dac(y_water):
    sim = sk.examples.Illustrative_Humid_DAC_Breakthrough.load()
    sim["layer[1].cells"] = 10
    sim["step[1].time"] = 300.0
    sim["endpoint[1].y[Water]"] = y_water
    sim["endpoint[1].y[Nitrogen]"] = 1.0 - y_water - 0.0004
    return sim


def test_humid_feed_below_saturation_is_silent():
    sim = _humid_dac(0.0115)
    result = sim.run()
    assert result.healthy and result.warnings == []
    cycle = sim.run_health().cycles[-1]
    assert 0.3 < cycle.supersaturation_gas < 1.0
    assert 0.3 < cycle.supersaturation_wall < 1.0


def test_supersaturated_feed_gives_the_code():
    sim = _humid_dac(0.03)
    assert sim["settings.warn_supersaturation"] is True
    result = sim.run()
    assert result.status in ("completed", "reached_css")  # a warning, not a failure
    assert not result.healthy
    (w,) = [w for w in result.warnings if w.code == "SUPERSATURATION"]
    assert w.component == "Water" and w.threshold == 1.0
    # The heat of water sorption warms the gas at the front, so the cooler wall
    # passes saturation first (the cold-wall measure).
    assert 1.0 < w.value < 1.4
    assert w.as_dict()["code"] == "SUPERSATURATION"
    assert "supersaturated" in w.message

    off = _humid_dac(0.03)
    off["settings.warn_supersaturation"] = False
    result = off.run()
    assert all(w.code != "SUPERSATURATION" for w in result.warnings)
    assert math.isnan(off.run_health().cycles[-1].supersaturation_gas)
