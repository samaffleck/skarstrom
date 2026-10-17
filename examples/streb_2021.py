"""Streb & Mazzotti (2021): CO2/CH4 breakthrough on zeolite 13X with thermocouples.

Experiments A1 (1 bar) and A3 (10 bar, the reference experiment) feed a
50:50 CO2/CH4 mixture into a helium-filled column. The outlet mole fractions
and the five centre-line thermocouples (z = 10, 35, 60, 85, 110 cm) are
compared with the digitised experiments (adsim/data/streb_2021/).

    python examples/streb_2021.py
    python examples/streb_2021.py --case A3

Reference: A. Streb, M. Mazzotti, Adsorption 27 (2021) 541-558.
doi:10.1007/s10450-021-00306-y
"""

from __future__ import annotations

import math
import sys

import numpy as np

import _common as c

CASES = {"A1": "streb_2021_A1.json", "A3": "streb_2021_A3.json"}
Z_CM = [10, 35, 60, 85, 110]


def main(argv=None) -> int:
    args = c.parser(__doc__.splitlines()[0], list(CASES)).parse_args(argv)
    plt = c.pyplot(args)
    rows = []
    for name in c.selected(args):
        sim = c.load(CASES[name], args.cycles)
        c.run(sim, f"Streb 2021 {name}", args.quiet)
        fig, axes = (plt.subplots(1, 2, figsize=(10, 4)) if plt else (None, None))
        d = c.DATA / "streb_2021" / name
        row = [name]
        for comp, label in (("CarbonDioxide", "CO2"), ("Methane", "CH4")):
            y = sim.series("mole_fraction", comp, location="top", cycle=0)
            t50 = c.crossing_time(y.time, y.values, 0.25)
            err = "-"
            f = d / f"{name}_mole_fraction_{label}.csv"
            if f.is_file():
                t_e, y_e = c.load_csv(f)
                err = c.fmt(c.nrmse(y.time, y.values, t_e, y_e, 0.5))
                if axes is not None:
                    axes[0].plot(t_e, y_e, "o", ms=3, label=f"{label} experiment")
            if axes is not None:
                axes[0].plot(y.time, y.values, "-", label=f"{label} Skarstrom")
            row += [c.fmt(t50), err]
        sq, n = 0.0, 0
        for z in Z_CM:
            T = sim.series("temperature", z=z / 100.0, cycle=0)
            f = d / f"{name}_temperature_z={z}[cm].csv"
            if f.is_file():
                t_e, T_e = c.load_csv(f)
                sq += float(np.sum((np.interp(t_e, T.time, T.values) - T_e) ** 2))
                n += len(t_e)
                if axes is not None:
                    axes[1].plot(t_e, T_e, "o", ms=2)
            if axes is not None:
                axes[1].plot(T.time, T.values, "-", label=f"z = {z} cm")
        row.append(c.fmt(math.sqrt(sq / n)) if n else "-")
        rows.append(row)
        if fig is not None:
            axes[0].set(xlabel="time [s]", ylabel="outlet mole fraction [-]")
            axes[1].set(xlabel="time [s]", ylabel="temperature [K] (lines: Skarstrom, markers: experiment)")
            for a in axes:
                a.legend(fontsize=7)
            c.save(fig, f"streb_2021_{name}.png")
    print()
    c.table(rows, ["case", "t50 CO2 [s]", "NRMSE y_CO2", "t50 CH4 [s]", "NRMSE y_CH4", "RMSE T, 5 TCs [K]"])
    return 0


if __name__ == "__main__":
    sys.exit(main())
