"""Six-step CO2/N2 VSA with product purge (template deck, no literature source).

Feed, blowdown, evacuation, pressurisation and product purge on a 13X-type
isotherm. The deck comes from an integration test and has no KPI
configuration, so the script reports the run statistics, the CSS check and
the outlet CO2 mole fraction over the last cycle. Use it as a starting point
for your own cycles, not as data.

    python examples/five_step_product_purge.py
"""

from __future__ import annotations

import sys

import _common as c


def main(argv=None) -> int:
    args = c.parser(__doc__.splitlines()[0], ["FiveStepProductPurge"]).parse_args(argv)
    sim = c.load("FiveStepProductPurge.json", args.cycles)
    result = c.run(sim, "FiveStepProductPurge", args.quiet)
    print(f"  steps: {', '.join(sim.steps)}")
    print(f"  {c.css_status(sim)}")
    print(f"  integrator: {result.stats['steps']} steps, {result.stats['residual_evaluations']} residual evaluations")
    c.print_kpis(sim)
    plt = c.pyplot(args)
    if plt:
        fig, ax = plt.subplots(figsize=(6, 4))
        for loc in ("bottom", "top"):
            y = sim.series("mole_fraction", "CarbonDioxide", location=loc, cycle=-1)
            ax.plot(y.time, y.values, label=f"y_CO2 at the {loc}")
        ax.set(xlabel="time in the last cycle [s]", ylabel="CO2 mole fraction [-]")
        ax.legend(fontsize=7)
        c.save(fig, "five_step_product_purge.png")
    return 0


if __name__ == "__main__":
    sys.exit(main())
