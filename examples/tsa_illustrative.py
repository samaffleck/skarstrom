"""Illustrative TSA: CO2 capture on zeolite 13X with indirect heating and cooling.

Not a reproduction of a paper: a four-step temperature swing cycle (feed,
heat, purge, cool) driven by heater and cooler energy sources, built for the
TSA tutorial in the documentation. The script runs it to cyclic steady state,
prints the CSS check and the KPIs, and plots the bed temperature at the top
of the column over the last cycle.

    python examples/tsa_illustrative.py
"""

from __future__ import annotations

import sys

import _common as c


def main(argv=None) -> int:
    args = c.parser(__doc__.splitlines()[0], ["TSA"]).parse_args(argv)
    sim = c.load("Illustrative_TSA_13X.json", args.cycles)
    c.run(sim, "Illustrative TSA", args.quiet)
    print(f"  {c.css_status(sim)}")
    c.print_kpis(sim)
    plt = c.pyplot(args)
    if plt:
        fig, ax = plt.subplots(figsize=(6, 4))
        for loc in ("bottom", "top"):
            T = sim.series("temperature", location=loc, cycle=-1)
            ax.plot(T.time, T.values, label=f"bed temperature, {loc}")
        ax.set(xlabel="time in the last cycle [s]", ylabel="temperature [K]")
        ax.legend(fontsize=7)
        c.save(fig, "tsa_illustrative.png")
    return 0


if __name__ == "__main__":
    sys.exit(main())
