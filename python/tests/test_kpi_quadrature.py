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


@pytest.mark.parametrize("name, expect_inflow", [("Perez2019_A", False), ("Perez2019_B", False),
                                                 ("Perez2019_C1", False), ("Perez2019_D1", True)])
def test_uncounted_inflow_diagnostic(name, expect_inflow):
    """Gross inflow through unconnected endpoints that are not feed or recycle: negligible
    (about 1e-6 of the feed, a start-of-step transient at the Feed outlet) for the Perez 2019 A,
    B and C1 decks, and about 0.4 % of the feed on D1, from its vent valves."""
    sim = getattr(sk.examples, name).load()
    for layer in sim.layers:
        layer.cells = 20
    result = sim.run(cycles=2)
    assert result.ok, result.message
    k = sim.kpis(cycle=-1)[0]
    d = k.as_dict()
    assert "uncounted_inflow" in d and "uncounted_inflow_fraction" in d
    if expect_inflow:
        assert k.uncounted_inflow > 0
        assert k.uncounted_inflow_fraction > 1e-3  # 0.42 % at CSS on the full grid; larger in two cycles on 20 cells
    else:
        assert 0 <= k.uncounted_inflow_fraction < 1e-5
