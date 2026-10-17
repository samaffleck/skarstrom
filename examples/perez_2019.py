"""Estupinan Perez et al. (2019): CO2/N2 VSA on zeolite 13X, configurations A-F.

A, B and C1 use the basic 4-step cycle; D1, E and F add light-product
pressurisation. Each deck runs to cyclic steady state; the CO2 purity and
recovery of the evacuation product (as the validation test defines them) are
compared with the measured values and the paper's detailed model
(Table 3).

    python examples/perez_2019.py                  # all six
    python examples/perez_2019.py --case A

Reference: L. Estupinan Perez, P. Sarkar, A. Rajendran, Separation and
Purification Technology 224 (2019) 553-563. doi:10.1016/j.seppur.2019.05.039
"""

from __future__ import annotations

import sys

import _common as c

CASES = {n: (f"Perez2019_{n}.json", f"Config{n}") for n in ["A", "B", "C1", "D1", "E", "F"]}
QUANTITIES = [
    ("CO2 purity vs experiment", "CarbonDioxide.purity", "CO2_purity_exp"),
    ("CO2 recovery vs experiment", "CarbonDioxide.recovery", "CO2_recovery_exp"),
    ("CO2 purity vs paper model", "CarbonDioxide.purity", "CO2_purity_model"),
    ("CO2 recovery vs paper model", "CarbonDioxide.recovery", "CO2_recovery_model"),
]


def main(argv=None) -> int:
    return c.cyclic_family(argv, __doc__, "perez_2019", CASES, QUANTITIES, "CarbonDioxide", "perez_2019.png",
                           values=c.evacuation_product_kpis)


if __name__ == "__main__":
    sys.exit(main())
