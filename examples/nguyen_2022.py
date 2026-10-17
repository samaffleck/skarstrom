"""Nguyen et al. (2022): CO2 breakthrough on CALF-20 (Figs. 4 and 5).

Fig. 4 is a 70 % CO2 feed and Fig. 5 a 15 % CO2 feed (balance helium) into a
bed initially full of helium. The outlet CO2 mole fraction is compared with the
digitised experiments (adsim/data/nguyen_2022/): NRMSE normalised by the feed
mole fraction, and the time at which the outlet reaches half the feed value.

    python examples/nguyen_2022.py
    python examples/nguyen_2022.py --case Fig4

Reference: T.T.T. Nguyen, J.-B. Lin, G.K.H. Shimizu, A. Rajendran, Chemical
Engineering Journal 442 (2022) 136263. doi:10.1016/j.cej.2022.136263
"""

from __future__ import annotations

import sys

import _common as c

CASES = {"Fig4": ("Nguyen2022_Fig4.json", 70), "Fig5": ("Nguyen2022_Fig5.json", 15)}


def main(argv=None) -> int:
    args = c.parser(__doc__.splitlines()[0], list(CASES)).parse_args(argv)
    plt = c.pyplot(args)
    fig, ax = (plt.subplots(figsize=(6, 4)) if plt else (None, None))
    rows = []
    for name in c.selected(args):
        deck, feed_pct = CASES[name]
        y_feed = feed_pct / 100.0
        sim = c.load(deck, args.cycles)
        c.run(sim, f"Nguyen 2022 {name}", args.quiet)
        y = sim.series("mole_fraction", "CarbonDioxide", location="top", cycle=0)
        t50 = c.crossing_time(y.time, y.values, 0.5 * y_feed)
        row = [name, f"{feed_pct} %", c.fmt(t50), "-", "-"]
        data = c.DATA / "nguyen_2022" / name / f"Y_CO2_{feed_pct}%.csv"
        if data.is_file():
            t_e, y_e = c.load_csv(data)
            row[3] = c.fmt(c.crossing_time(t_e, y_e, 0.5 * y_feed))
            row[4] = c.fmt(c.nrmse(y.time, y.values, t_e, y_e, y_feed))
            if ax is not None:
                ax.plot(t_e, y_e, "o", ms=3, label=f"{name} experiment")
        if ax is not None:
            ax.plot(y.time, y.values, "-", label=f"{name} Skarstrom")
        rows.append(row)
    print()
    c.table(rows, ["case", "feed CO2", "t50 simulated [s]", "t50 experiment [s]", "NRMSE y_CO2"])
    if fig is not None:
        ax.set(xlabel="time [s]", ylabel="outlet CO2 mole fraction [-]")
        ax.legend(fontsize=7)
        c.save(fig, "nguyen_2022.png")
    return 0


if __name__ == "__main__":
    sys.exit(main())
