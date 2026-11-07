"""KPI quadratures: the KPIs integrated by IDAS alongside the solution, reported
beside the (trapezoid) KPIs, and the cycle mole-balance error."""
import math

import pytest

import skarstrom as sk

FIELDS = ["purity", "recovery", "productivity_mass", "productivity_volume", "production_rate"]


def test_quadrature_kpis_sit_beside_the_reported_kpis():
    sim = sk.examples.Perez2019_A.load()
    for layer in sim.layers:
        layer.cells = 20
    result = sim.run(cycles=2)
    assert result.ok, result.message
    for k in sim.kpis():
        d = k.as_dict()
        for name in FIELDS:
            assert d["quadrature_" + name] == pytest.approx(d[name], rel=1e-4), name
        assert "quadrature_energy_consumption" in d
        assert abs(k.quadrature_mass_balance_error) < 1e-3
        assert math.isfinite(k.mass_balance_error)
