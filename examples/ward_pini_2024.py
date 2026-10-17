"""Ward & Pini (2024): CO2/N2 4-step VSA on Mg-MOF-74.

Runs the deck to cyclic steady state and compares CO2 purity and recovery
with the paper's optimised detailed-model result (Table 4, maximum separation
effectiveness, single adsorbent).

    python examples/ward_pini_2024.py

Reference: A. Ward, R. Pini, Industrial & Engineering Chemistry Research 63
(2024) 13787-13800. doi:10.1021/acs.iecr.4c01345
"""

from __future__ import annotations

import sys

import _common as c

CASES = {"FourStepPSA": ("Ward2024.json", "FourStepPSA")}
QUANTITIES = [
    ("CO2 purity", "CarbonDioxide.purity", "CO2_purity"),
    ("CO2 recovery", "CarbonDioxide.recovery", "CO2_recovery"),
]


def main(argv=None) -> int:
    return c.cyclic_family(argv, __doc__, "ward_pini_2024", CASES, QUANTITIES, "CarbonDioxide", "ward_pini_2024.png")


if __name__ == "__main__":
    sys.exit(main())
