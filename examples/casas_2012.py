"""Casas et al. (2012): CO2/H2 breakthrough on activated carbon AP3-60 at 15 bar, 25 C.

Runs the three feed compositions (50/50 of Fig. 4, 25/75 and 75/25 of Fig. 6)
and compares the H2 and CO2 breakthrough times (outlet mole fraction above
0.002, the paper's Sect. 4.4 definition) with the paper's simulated values
(Table 5). The isotherms are the Sips fits of Schell et al. (2012).

    python examples/casas_2012.py
    python examples/casas_2012.py --case Fig4_50_50

Reference: N. Casas, J. Schell, R. Pini, M. Mazzotti, Adsorption 18 (2012)
143-161. doi:10.1007/s10450-012-9389-z
"""

from __future__ import annotations

import sys

import _common as c

CASES = {
    "Fig4_50_50": ("Casas2012_Fig4_50_50_15bar_25C.json", "Fig4_Reference_50_50_15bar_25C"),
    "Fig6_25_75": ("Casas2012_Fig6_25_75_15bar_25C.json", "Fig6_25_75_15bar_25C"),
    "Fig6_75_25": ("Casas2012_Fig6_75_25_15bar_25C.json", "Fig6_75_25_15bar_25C"),
}
THRESHOLD = 0.002


def main(argv=None) -> int:
    args = c.parser(__doc__.splitlines()[0], list(CASES)).parse_args(argv)
    plt = c.pyplot(args)
    fig, ax = (plt.subplots(figsize=(6, 4)) if plt else (None, None))
    rows = []
    for name in c.selected(args):
        deck, ref_case = CASES[name]
        sim = c.load(deck, args.cycles)
        c.run(sim, f"Casas 2012 {name}", args.quiet)
        for comp, label in (("Hydrogen", "H2"), ("CarbonDioxide", "CO2")):
            s = sim.series("mole_fraction", comp, location="top", cycle=0)
            t = c.crossing_time(s.time, s.values, THRESHOLD)
            ref = c.reference("casas_2012", ref_case, f"{label}_breakthrough_time")
            err = "-" if ref is None else f"{100 * (t - ref) / ref:+.1f} %"
            rows.append([name, label, c.fmt(t), c.fmt(ref), err])
            if ax is not None:
                ax.plot(s.time, s.values, label=f"{name} {label}")
    print()
    c.table(rows, ["case", "component", "t_br simulated [s]", "t_br paper model [s]", "difference"])
    if fig is not None:
        ax.set(xlabel="time [s]", ylabel="outlet mole fraction [-]")
        ax.legend(fontsize=7)
        c.save(fig, "casas_2012.png")
    return 0


if __name__ == "__main__":
    sys.exit(main())
