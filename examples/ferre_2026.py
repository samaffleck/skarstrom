"""Ferre et al. (2026): biogas upgrading (CO2/CH4) by 4-step PVSA on zeolite NaY.

Runs the NaY (Feng isotherm data) deck to cyclic steady state and compares
CO2 and CH4 purity and recovery with the paper's detailed-model result
(Table 5, recovery-recovery case).

    python examples/ferre_2026.py

Reference: A. Ferre, C. Wedler, A. Ward, J. Burger, R. Pini, C. Petit,
Separation and Purification Technology 406 (2026) 138821.
doi:10.1016/j.seppur.2026.138821
"""

from __future__ import annotations

import sys

import _common as c

CASES = {"NaY_F": ("Ferre2026_NaYF.json", "NaY_F")}
QUANTITIES = [
    ("CO2 purity", "CarbonDioxide.purity", "CO2_purity"),
    ("CO2 recovery", "CarbonDioxide.recovery", "CO2_recovery"),
    ("CH4 purity", "Methane.purity", "CH4_purity"),
    ("CH4 recovery", "Methane.recovery", "CH4_recovery"),
]


def kpis(sim):
    """The deck has no KPI configuration; integrate the stream logs of the last cycle.

    Steps: 0 Feed, 1 Blowdown, 2 Evacuation, 3 Pressurisation. CO2 is the heavy
    product (evacuation, bottom), CH4 the light product (feed raffinate, top);
    the feed is the feed and pressurisation inlets (bottom), as in the test.
    """
    m = c.endpoint_moles
    out = {}
    for comp, (step, loc) in (("CarbonDioxide", (2, "bottom")), ("Methane", (0, "top"))):
        fed = m(sim, 0, "bottom", comp) + m(sim, 3, "bottom", comp)
        out[f"{comp}.purity"] = m(sim, step, loc, comp) / m(sim, step, loc)
        out[f"{comp}.recovery"] = m(sim, step, loc, comp) / fed
    return out


def main(argv=None) -> int:
    return c.cyclic_family(argv, __doc__, "ferre_2026", CASES, QUANTITIES, "Methane", "ferre_2026.png", values=kpis)


if __name__ == "__main__":
    sys.exit(main())
