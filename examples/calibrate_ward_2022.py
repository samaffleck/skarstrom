"""Calibrate the Ward 2022 breakthrough decks against the measured curves with skarstrom.calibrate.

The shipped Ward 2022 decks (``Ward2022_Q50ccm``, ``_Q100ccm``, ``_Q200ccm``) were
calibrated by hand: a weighted least-squares fit of the three breakthrough runs
together to the digitised outlet CO2 mole fraction and the thermocouple at
z = 0.0715 m (core/test/papers/ward_2022/README.md, "Calibration"). This script
does the same fit with the calibration API, so it is the reproducible method of
that calibration and a template for a fit of your own:

- three cases, one per flow rate (50, 100, 200 cm3/min), each with its two bound
  datasets (the decks carry them, ``sim.experimental_data``);
- the residual weights of the hand fit: an uncertainty of 0.03 on each outlet
  mole fraction and 0.5 K on each thermocouple reading (the assumed digitising
  error), every point counted once (dataset weight = its number of points);
- shared parameters: q_s,2 of the CO2 isotherm and the factor on the paper's
  dispersion correlation (Eq. 4), the factor as a user variable that the three
  decks' dispersion coefficients use; the ambient and wall-end temperature, also
  a shared user variable, is held at 293 K as in the hand calibration
  (``--fit-ambient`` fits it too);
- one inside heat transfer coefficient h_in per case;
- the start at the paper's values (q_s,2 6.04 mol/kg, factor 4, h_in 77.5 / 57.9 /
  45.2 W/(m2 K)), so the fit has to find the calibration, not confirm it.

    python examples/calibrate_ward_2022.py            # the full fit (60 cells, central differences)
    python examples/calibrate_ward_2022.py --quick    # coarse grid and budget: a look in seconds

The bundled decks are never changed. The calibrated decks are written to a new
temporary folder (or ``--out DIR``).

Reference: A. Ward, R. Pini, Adsorption 28 (2022) 161-183.
doi:10.1007/s10450-022-00361-z
"""

from __future__ import annotations

import argparse
import tempfile
from pathlib import Path

import _common as c
import skarstrom as sk
from skarstrom import calibration as cal

FLOWS = ("50", "100", "200")
DECKS = {q: f"Ward2022_Q{q}ccm" for q in FLOWS}
OUTLET = "outlet_CO2_mole_fraction[-]_Qin={q}[ccm]"
THERMOCOUPLE = "temperature[K]_z=0.0715[m]_Qin={q}[ccm]"

SIGMA_Y = 0.03   # assumed digitising error of the outlet mole fraction [-]
SIGMA_T = 0.5    # assumed digitising error of the thermocouple [K]

# The paper's values (the start of the fit) and the hand calibration (README, "Calibration").
PAPER = {"q_s2": 6.04, "factor": 4.0, "h_in": {"50": 77.5, "100": 57.9, "200": 45.2}}
HAND = {"q_s2": (5.89, 0.01), "factor": (3.75, 0.12),
        "h_in": {"50": (97.3, 6.6), "100": (81.5, 2.4), "200": (69.0, 1.1)}}
DECK_FACTOR = 3.75   # the factor the shipped decks' dispersion coefficients carry

Q_S2 = "adsorbent[1].isotherm[CarbonDioxide].m2"
FACTOR = "variable[dispersion_factor]"
AMBIENT = "variable[ambient_temperature]"
H_IN = "layer[1].internal_htc"


def case(q: str, cells: int | None = None) -> cal.Case:
    """One flow rate: the shipped deck, made ready for the fit (in memory only)."""
    sim = c.load(DECKS[q])
    # The deck's dispersion coefficient is 3.75 x Eq. 4 at this flow. Write it as
    # dispersion_factor x Eq. 4, so one shared variable is the factor of all three decks.
    eq4 = sim.get("layer[1].dispersion_coefficient") / DECK_FACTOR
    sim.add_variable("dispersion_factor", DECK_FACTOR)
    sim["layer[1].dispersion_coefficient"] = f"dispersion_factor * {eq4!r}"
    # The ambient and the wall ends (held at a fixed temperature, Table 4) take one value.
    sim.add_variable("ambient_temperature", sim.get("bed.ambient_temperature"), category="temperature")
    sim["bed.ambient_temperature"] = "ambient_temperature"
    sim["wall.fixed_temperature"] = "ambient_temperature"
    if cells is not None:
        sim["layer[1].cells"] = cells
    # The weights of the hand fit: r = (simulated - measured) / sigma for every point.
    # The objective divides a dataset's weight by its number of points, so a weight
    # equal to that number counts every point once, as the hand fit did.
    datasets = []
    for name, sigma in ((OUTLET.format(q=q), SIGMA_Y), (THERMOCOUPLE.format(q=q), SIGMA_T)):
        ds = sim.experimental_data[name]
        ds.uncertainty = [sigma] * len(ds.x)
        datasets.append((name, float(len(ds.x))))
    return cal.Case(sim, datasets, name=f"{q} ccm")


def parameters(fit_ambient: bool) -> list:
    params = [
        cal.Parameter(Q_S2, 5.0, 7.0, start=PAPER["q_s2"]),                       # mol/kg, shared
        cal.Parameter(FACTOR, 1.0, 8.0, start=PAPER["factor"]),                   # -, shared
        cal.Parameter(H_IN, 10.0, 500.0, shared=False,                            # W/(m2 K), per case
                      case_start=[PAPER["h_in"][q] for q in FLOWS]),
    ]
    if fit_ambient:
        params.append(cal.Parameter(AMBIENT, 289.0, 296.0, start=293.0))         # K, shared
    return params


def comparison(result: cal.CalibrationResult) -> list:
    """Rows: value, hand calibration (± 1 sigma), API (± 1 sigma)."""
    def api(name):
        p = result[name]
        return p.value, p.std_error

    rows = [("q_s,2 [mol/kg]", HAND["q_s2"], api(Q_S2), "{:.2f}", "{:.2f}"),
            ("factor on Eq. 4 [-]", HAND["factor"], api(FACTOR), "{:.2f}", "{:.2f}")]
    for i, q in enumerate(FLOWS):
        rows.append((f"h_in, {q} ccm [W/(m2 K)]", HAND["h_in"][q], api(f"{H_IN} [{q} ccm]"), "{:.1f}", "{:.1f}"))
    if any(p.path == AMBIENT for p in result.parameters):
        rows.append(("ambient [K]", (293.0, float("nan")), api(AMBIENT), "{:.1f}", "{:.2f}"))
    return [(label, f"{fv.format(h[0])} ± {fe.format(h[1])}", f"{fv.format(a[0])} ± {fe.format(a[1])}")
            for label, h, a, fv, fe in rows]


def main(argv=None) -> int:
    ap = argparse.ArgumentParser(description=__doc__.splitlines()[0])
    budget = ap.add_mutually_exclusive_group()
    budget.add_argument("--quick", action="store_true",
                        help="20 cells, forward differences and at most 40 runs: shows the method in seconds")
    budget.add_argument("--full", action="store_true", help="the decks' 60 cells, central differences (default)")
    ap.add_argument("--fit-ambient", action="store_true",
                    help="also fit the ambient and wall-end temperature (the hand calibration set it to 293 K)")
    ap.add_argument("--threads", type=int, help="simulations at once (default: the cores - 1)")
    ap.add_argument("--out", type=Path, help="folder for the calibrated decks (default: a new temporary folder)")
    ap.add_argument("--no-plot", action="store_true", help="do not save a plot")
    ap.add_argument("-q", "--quiet", action="store_true")
    args = ap.parse_args(argv)

    cells = 20 if args.quick else None
    options = (dict(differences="forward", max_simulations=40, max_iterations=4) if args.quick
               else dict(differences="central", max_simulations=600))
    cases = [case(q, cells) for q in FLOWS]

    def progress(info):
        if not args.quiet:
            print(f"  iteration {info['iteration']}: S = {info['objective']:.6g} "
                  f"({info['simulations']} simulations)", flush=True)

    print(f"Calibrating Ward 2022 at 50, 100 and 200 ccm ({'quick' if args.quick else 'full'} budget)")
    result = sk.calibrate(cases, parameters(args.fit_ambient), threads=args.threads, progress=progress, **options)
    print()
    print(result.table())
    print()

    if args.quick:
        print("(a quick fit: coarse grid and few iterations; run without --quick to reproduce the hand calibration)")
    print("Hand calibration (README) and the API:")
    c.table(comparison(result), ["value", "hand calibration", "skarstrom.calibrate"])
    print()

    if not result.ok:
        print(f"the fit failed: {result.message}")
        return 1

    out = args.out or Path(tempfile.mkdtemp(prefix="ward2022_calibrated_"))
    out.mkdir(parents=True, exist_ok=True)
    rows = []
    for i, q in enumerate(FLOWS):
        sim = result.calibrated(i)
        path = out / f"{DECKS[q]}_calibrated.json"
        sim.save(str(path))
        status = sim.run().status   # a check that the calibrated deck runs on its own
        fits = [d for d in result.datasets if d.case == i]
        rows.append([f"{q} ccm", status] + [c.fmt(d.rmse, 3) for d in fits] + [str(path)])
    c.table(rows, ["case", "calibrated deck", "y_CO2 RMSE [-]", "T RMSE [K]", "written to"])

    if not args.no_plot:
        plt = c.pyplot(args)
        if plt is not None:
            fig = result.plot(start=False)
            c.save(fig, "calibrate_ward_2022.png")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
