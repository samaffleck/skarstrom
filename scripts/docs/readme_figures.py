"""Draw the README figures from a CLI sweep of the Perez2019_A example.

The sweep is the loop in the README's quick start: one ``skarstrom examples
run Perez2019_A --json -q`` per feed time and evacuation pressure, each
summary saved as ``run_<P>kPa_<T>s.json``. This script reads those files and
writes a light and a dark PNG for the README:

    python scripts/docs/readme_figures.py SWEEP_DIR            # -> docs/assets/readme-*.png
    python scripts/docs/readme_figures.py SWEEP_DIR --out DIR

Needs matplotlib. Not run in CI: the figures are committed images.
"""

from __future__ import annotations

import argparse
import json
import re
from pathlib import Path

import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt  # noqa: E402

ROOT = Path(__file__).resolve().parents[2]
PURITY_TARGET = 0.99

THEMES = {
    "light": dict(surface="#fcfcfb", text="#0b0b0b", muted="#52514e", grid="#e4e3df",
                  series=["#2a78d6", "#eb6834", "#1baf7a"]),
    "dark": dict(surface="#1a1a19", text="#f0efec", muted="#c3c2b7", grid="#34332f",
                 series=["#3987e5", "#d95926", "#199e70"]),
}


def load(sweep: Path):
    runs = {}
    for f in sweep.glob("run_*kPa_*s.json"):
        p, t = map(float, re.match(r"run_([\d.]+)kPa_([\d.]+)s\.json", f.name).groups())
        summary = json.loads(f.read_text())
        if summary["status"] not in ("completed", "reached_css"):
            continue
        kpi = next(k for k in summary["kpis"] if k["cycle"] == summary["cycles_run"] - 1
                   and k["component"] == "CarbonDioxide")
        runs.setdefault(p, []).append((t, kpi["purity"], kpi["recovery"]))
    return {p: sorted(v) for p, v in sorted(runs.items())}


def best(runs):
    ok = [(r, p, t, y) for p, v in runs.items() for t, y, r in v if y >= PURITY_TARGET]
    return max(ok)


def draw(runs, theme: str, out: Path):
    c = THEMES[theme]
    plt.rcParams.update({"font.family": "DejaVu Sans", "font.size": 11, "text.color": c["text"],
                         "axes.labelcolor": c["muted"], "xtick.color": c["muted"],
                         "ytick.color": c["muted"], "axes.edgecolor": c["grid"]})
    fig, axes = plt.subplots(1, 2, figsize=(10, 4.2), facecolor=c["surface"], sharex=True)
    r_best, p_best, t_best, y_best = best(runs)
    for ax, key, label in ((axes[0], 2, "CO₂ recovery [%]"), (axes[1], 1, "CO₂ purity [%]")):
        ax.set_facecolor(c["surface"])
        ax.grid(True, axis="y", color=c["grid"], linewidth=0.8)
        ax.set_axisbelow(True)
        for s in ("top", "right", "left"):
            ax.spines[s].set_visible(False)
        ax.tick_params(length=0)
        for colour, (p, v) in zip(c["series"], runs.items()):
            t = [x[0] for x in v]
            y = [100 * x[key] for x in v]
            ax.plot(t, y, color=colour, linewidth=2, marker="o", markersize=7,
                    markeredgecolor=c["surface"], markeredgewidth=2, label=f"{p:g} kPa vacuum")
            if key == 2:
                ax.annotate(f"{p:g} kPa", (t[-1], y[-1]), xytext=(8, 0), textcoords="offset points",
                            va="center", color=c["text"], fontsize=10)
        y_mark = 100 * (r_best if key == 2 else y_best)
        ax.plot([t_best], [y_mark], marker="o", markersize=15, markerfacecolor="none",
                markeredgecolor=c["text"], markeredgewidth=1.6)
        ax.set_xlabel("feed step time [s]")
        ax.set_title(label, loc="left", color=c["text"], fontsize=11)
    axes[0].set_xlim(right=max(t for v in runs.values() for t, _, _ in v) + 22)
    axes[0].annotate(f"best with ≥ {100 * PURITY_TARGET:g} % purity:\n{p_best:g} kPa, {t_best:g} s feed "
                     f"→ {100 * r_best:.0f} % recovery",
                     (t_best, 100 * r_best), xytext=(t_best - 25, 100 * r_best - 10), textcoords="data",
                     ha="center", va="center", color=c["text"], fontsize=9.5,
                     arrowprops=dict(arrowstyle="-", color=c["muted"], linewidth=0.8))
    axes[1].axhline(100 * PURITY_TARGET, color=c["muted"], linewidth=1, linestyle=(0, (4, 3)))
    axes[1].annotate(f"{100 * PURITY_TARGET:g} % target", (axes[1].get_xlim()[1], 100 * PURITY_TARGET),
                     xytext=(-4, -12), textcoords="offset points", ha="right", color=c["muted"], fontsize=9)
    leg = axes[1].legend(loc="lower right", frameon=False, fontsize=9.5)
    for text in leg.get_texts():
        text.set_color(c["text"])
    fig.suptitle("Perez 2019 VSA (15 % CO₂ in N₂, zeolite 13X): 18 runs to cyclic steady state",
                 x=0.012, ha="left", color=c["text"], fontsize=12.5, fontweight="bold")
    fig.tight_layout(rect=(0, 0, 1, 0.97))
    fig.savefig(out, dpi=160, facecolor=c["surface"])
    plt.close(fig)


def main() -> None:
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("sweep", type=Path, help="directory with the run_<P>kPa_<T>s.json summaries")
    ap.add_argument("--out", type=Path, default=ROOT / "docs" / "assets")
    args = ap.parse_args()
    runs = load(args.sweep)
    for theme in THEMES:
        path = args.out / f"readme-sweep-{theme}.png"
        draw(runs, theme, path)
        print(path)


if __name__ == "__main__":
    main()
