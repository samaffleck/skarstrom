"""Arora & Hasan (2021): medical oxygen from air on LiLSX, 4-step PSA and PVSA.

Runs both decks (50 cycles each) and compares the O2 purity of the collected
product and the net O2 recovery with the paper's optimised model results
(Table 4). The deck's KPI configuration subtracts the O2 product recycled to
purge and pressurisation, so the recovery is net, as in the paper.

    python examples/arora_2021.py
    python examples/arora_2021.py --case LiLSX_PVSA

Reference: A. Arora, M.M.F. Hasan, Flexible oxygen concentrators for medical
applications, Scientific Reports 11 (2021) 14317. doi:10.1038/s41598-021-93796-3
"""

from __future__ import annotations

import sys

import _common as c

CASES = {
    "LiLSX_PSA": ("Arora2021_LiLSX_PSA.json", "LiLSX_PSA"),
    "LiLSX_PVSA": ("Arora2021_LiLSX_PVSA.json", "LiLSX_PVSA"),
}
QUANTITIES = [("O2 purity", "Oxygen.purity", "O2_purity"), ("O2 recovery", "Oxygen.recovery", "O2_recovery")]


def kpis(sim):
    """Purity of the collected product, recovery net of the recycled oxygen.

    The deck's KPI configuration counts the purge and pressurisation inlets as
    (negative) product, which gives the paper's net recovery but makes the
    purity of that net product meaningless (undefined, or above 100 % for the
    PVSA deck). As in the validation test, the purity is that of the
    production outlet (step 0, top), integrated over the last cycle.
    """
    o2 = c.endpoint_moles(sim, 0, "top", "Oxygen")
    return {"Oxygen.purity": o2 / c.endpoint_moles(sim, 0, "top"),
            "Oxygen.recovery": c.kpi(sim, "Oxygen").recovery}


def main(argv=None) -> int:
    return c.cyclic_family(argv, __doc__, "arora_hasan_2021", CASES, QUANTITIES, "Oxygen", "arora_2021.png",
                           values=kpis)


if __name__ == "__main__":
    sys.exit(main())
