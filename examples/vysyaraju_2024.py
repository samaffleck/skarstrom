"""Vysyaraju et al. (2024): helium recovery from He/N2 by 5-step PSA on binderless 5A.

Runs the experiment-1 deck to cyclic steady state and compares the He purity
and recovery with the measured values and the paper's model (Table 3).

    python examples/vysyaraju_2024.py

Reference: R. Vysyaraju, J.A. Sawada, L. Estupinan Perez, A. Rajendran,
Chemical Engineering Journal 499 (2024) 155929. doi:10.1016/j.cej.2024.155929
"""

from __future__ import annotations

import sys

import _common as c

CASES = {"PSA_Exp1": ("Vysyaraju2024_PSA_exp1.json", "PSA_Exp1")}
QUANTITIES = [
    ("He purity vs experiment", "Helium.purity", "He_purity_exp"),
    ("He recovery vs experiment", "Helium.recovery", "He_recovery_exp"),
    ("He purity vs paper model", "Helium.purity", "He_purity_model"),
    ("He recovery vs paper model", "Helium.recovery", "He_recovery_model"),
]


def kpis(sim):
    """The deck has no KPI configuration; integrate the stream logs of the last cycle.

    He purity is the He fraction of the adsorption-step top product (step 1);
    He recovery is 1 - (He lost in the blowdown, step 3 bottom) / (He fed in
    steps 0 and 1, bottom), the conservation form the validation test uses.
    """
    m = c.endpoint_moles
    fed = m(sim, 0, "bottom", "Helium") + m(sim, 1, "bottom", "Helium")
    return {
        "Helium.purity": m(sim, 1, "top", "Helium") / m(sim, 1, "top"),
        "Helium.recovery": 1.0 - m(sim, 3, "bottom", "Helium") / fed,
    }


def main(argv=None) -> int:
    return c.cyclic_family(argv, __doc__, "vysyarajua_2024", CASES, QUANTITIES, "Helium", "vysyaraju_2024.png",
                           values=kpis)


if __name__ == "__main__":
    sys.exit(main())
