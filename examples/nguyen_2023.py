"""Nguyen et al. (2023): CO2/N2 VPSA on structured CALF-20, experiments I-X.

Experiments I-V use the basic 4-step cycle, VI-X the 4-step cycle with
light-product pressurisation. Each deck runs to cyclic steady state; the CO2
purity and recovery of the evacuation product (as the validation test defines
them) are compared with the measured values and with the
paper's detailed model (Table 2).

    python examples/nguyen_2023.py                 # all ten (several minutes)
    python examples/nguyen_2023.py --case I --case VI

Reference: T.T.T. Nguyen, G.K.H. Shimizu, A. Rajendran, Chemical Engineering
Journal 452 (2023) 139550. doi:10.1016/j.cej.2022.139550
"""

from __future__ import annotations

import sys

import _common as c

NUMERALS = ["I", "II", "III", "IV", "V", "VI", "VII", "VIII", "IX", "X"]
CASES = {n: (f"Nguyen2023_Exp_{n}.json", f"Experiment{n}") for n in NUMERALS}
QUANTITIES = [
    ("CO2 purity vs experiment", "CarbonDioxide.purity", "CO2_purity_exp"),
    ("CO2 recovery vs experiment", "CarbonDioxide.recovery", "CO2_recovery_exp"),
    ("CO2 purity vs paper model", "CarbonDioxide.purity", "CO2_purity_model"),
    ("CO2 recovery vs paper model", "CarbonDioxide.recovery", "CO2_recovery_model"),
]


def main(argv=None) -> int:
    return c.cyclic_family(argv, __doc__, "nguyen_2023", CASES, QUANTITIES, "CarbonDioxide", "nguyen_2023.png",
                           values=c.evacuation_product_kpis)


if __name__ == "__main__":
    sys.exit(main())
