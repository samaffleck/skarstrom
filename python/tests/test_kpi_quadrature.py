"""KPI quadratures: the reported KPIs are integrated by IDAS alongside the
solution; the trapezoid KPIs are kept beside them (trapezoid_*), with the
cycle mole-balance error of each."""
import math

import pytest

import skarstrom as sk

FIELDS = ["purity", "recovery", "productivity_mass", "productivity_volume", "production_rate"]


def test_trapezoid_kpis_sit_beside_the_reported_kpis():
    sim = sk.examples.Perez2019_A.load()
    for layer in sim.layers:
        layer.cells = 20
    result = sim.run(cycles=2)
    assert result.ok, result.message
    for k in sim.kpis():
        d = k.as_dict()
        for name in FIELDS:
            assert d["trapezoid_" + name] == pytest.approx(d[name], rel=1e-4), name
        assert "trapezoid_energy_consumption" in d
        assert "quadrature_purity" not in d
        assert abs(k.mass_balance_error) < 1e-6  # component flows on the face basis
        assert math.isfinite(k.trapezoid_mass_balance_error)
