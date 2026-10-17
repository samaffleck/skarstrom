"""Ward & Pini (2022): CO2 breakthrough on activated carbon at 50, 100 and 200 ccm.

Runs the three single-step breakthrough decks (pure CO2 into a bed initially
full of air) and compares the outlet CO2 mole fraction and the thermocouple
at z = 0.0715 m with the digitised experiments (adsim/data/ward_2022/). The
breakthrough time (outlet y_CO2 = 0.5) is set against the 95 % credibility
interval of the paper's own model (Table 10).

    python examples/ward_2022.py                 # all three flow rates
    python examples/ward_2022.py --case 100      # the fastest single case

Reference: A. Ward, R. Pini, Adsorption 28 (2022) 161-183.
doi:10.1007/s10450-022-00361-z
"""

from __future__ import annotations

import sys

import _common as c

CASES = {"50": "Ward_2022_Q=50[ccm].json", "100": "Ward_2022_Q=100[ccm].json", "200": "Ward_2022_Q=200[ccm].json"}
Z_TC = 0.0715  # thermocouple position [m]


def main(argv=None) -> int:
    args = c.parser(__doc__.splitlines()[0], list(CASES)).parse_args(argv)
    plt = c.pyplot(args)
    fig, axes = (plt.subplots(1, 2, figsize=(10, 4)) if plt else (None, None))
    rows = []
    for q in c.selected(args):
        sim = c.load(CASES[q], args.cycles)
        c.run(sim, f"Ward 2022, {q} ccm", args.quiet)
        y = sim.series("mole_fraction", "CarbonDioxide", location="top", cycle=0)
        T = sim.series("temperature", z=Z_TC, cycle=0)
        t_br = c.crossing_time(y.time, y.values, 0.5)
        lo = c.reference("ward_2022", f"Q{q}ccm", "breakthrough_time_CI_low")
        hi = c.reference("ward_2022", f"Q{q}ccm", "breakthrough_time_CI_high")
        row = [f"{q} ccm", c.fmt(t_br), "-" if lo is None else f"[{lo:g}, {hi:g}]", "-", "-"]
        data_dir = c.DATA / "ward_2022" / f"Qin={q}"
        if data_dir.is_dir():
            t_e, y_e = c.load_csv(data_dir / f"outlet_CO2_mole_fraction[-]_Qin={q}[ccm].csv")
            tT_e, T_e = c.load_csv(data_dir / f"temperature[K]_z=0.0715[m]_Qin={q}[ccm].csv")
            row[3] = c.fmt(c.nrmse(y.time, y.values, t_e, y_e))  # feed is pure CO2: already normalised
            row[4] = c.fmt(c.nrmse(T.time, T.values, tT_e, T_e))
            if axes is not None:
                axes[0].plot(t_e, y_e, "o", ms=3, label=f"{q} ccm, experiment")
                axes[1].plot(tT_e, T_e, "o", ms=3, label=f"{q} ccm, experiment")
        if axes is not None:
            axes[0].plot(y.time, y.values, "-", label=f"{q} ccm, Skarstrom")
            axes[1].plot(T.time, T.values, "-", label=f"{q} ccm, Skarstrom")
        rows.append(row)
    print()
    c.table(rows, ["case", "t_br (y=0.5) [s]", "paper model 95% CI [s]", "NRMSE y_CO2", "RMSE T [K]"])
    if not c.have_reference_data():
        print("(run from a repository checkout to compare with the digitised data)")
    if fig is not None:
        axes[0].set(xlabel="time [s]", ylabel="outlet CO2 mole fraction [-]")
        axes[1].set(xlabel="time [s]", ylabel=f"temperature at z = {Z_TC} m [K]")
        for a in axes:
            a.legend(fontsize=7)
        c.save(fig, "ward_2022.png")
    return 0


if __name__ == "__main__":
    sys.exit(main())
