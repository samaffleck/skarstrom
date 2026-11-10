"""Illustrative DAC breakthrough: 400 ppm CO2 in N2 on Lewatit VP OC 1065.

Not a reproduction of a paper: an isothermal breakthrough of air-like feed (400 ppm
CO2 in N2) through a small packed bed of the library adsorbent "Lewatit VP OC 1065",
built to show the two direct-air-capture kinetics together:

* the dual-rate LDF (``dual_rate_ldf``): a fast site holding a fraction of the
  equilibrium loading and a slow site holding the rest, so the outlet shows an early
  front and a slow tail (cf. Stampi-Bombelli et al. 2024, doi:10.1021/acs.iecr.4c01392);
* the series-resistance coefficient (``series_resistance``): external film, macropore and
  Arrhenius micropore resistances in series, the first two scaled by the capacity factor
  taken from the isotherm slope. The slow site has its own micropore diffusivity.

The rate parameters are round numbers chosen to show the behaviour, not fitted values.
The script builds the deck from scratch, runs it and prints the times at which the outlet
reaches 5 %, 50 % and 95 % of the feed.

    python examples/illustrative_dac_breakthrough.py

``build()`` returns the configured simulation without running it
(python/tests/test_examples_build.py checks it against examples/decks/).
"""

from __future__ import annotations

import sys

import numpy as np

import _common as c
import skarstrom as sk
from skarstrom import ureg as u

FEED_CO2 = 400e-6          # mole fraction


def build(case: str = "DAC") -> sk.Simulation:
    """The illustrative DAC breakthrough, ready to run (there is one case, "DAC")."""
    sim = sk.Simulation.new(
        title="Illustrative DAC breakthrough - 400 ppm CO2 on Lewatit VP OC 1065",
        description="ILLUSTRATIVE, not a validated reproduction. Isothermal breakthrough of 400 ppm CO2 in N2 on the "
                    "library adsorbent Lewatit VP OC 1065, with the dual-rate LDF (fast and slow site) and the "
                    "series-resistance coefficient (film, macropore and micropore resistances with the capacity "
                    "factor). Rate parameters are round numbers chosen to show an early front and a slow tail.",
        author="Skarstrom examples",
    )
    sim.add_component("Nitrogen")
    sim.add_component("CarbonDioxide")

    # ---- adsorbent: the library entry, with the DAC kinetics for CO2 -----------------------
    ads_id = sim.add_library_adsorbent("lewatit_vp_oc_1065_low2023")
    ads = sim.adsorbents.by_id(ads_id)
    ads.set_rate_model("Nitrogen", "ldf", "constant", k=0.1 / u.s)
    # fast site: 70 % of the capacity; the slow site has a twenty times smaller micropore diffusivity
    ads.set_rate_model("CarbonDioxide", "dual_rate_ldf", "series_resistance", split=0.7,
                       resistances=["film", "macropore", "micropore"], capacity_factor_mode="slope",
                       tortuosity=3.0, macropore_void_fraction=0.338,
                       crystal_radius=1e-5 * u.m, micropore_diffusivity=2e-10 * u.m**2 / u.s,
                       micropore_activation_energy=20 * u.kJ / u.mol,
                       slow_tortuosity=3.0, slow_macropore_void_fraction=0.338,
                       slow_crystal_radius=1e-5 * u.m, slow_micropore_diffusivity=1e-12 * u.m**2 / u.s,
                       slow_micropore_activation_energy=20 * u.kJ / u.mol)

    # ---- column ---------------------------------------------------------------------------
    sim.bed.set_shape("cylinder", diameter=2 * u.cm)
    # The wall does not matter in an isothermal run.
    sim.wall.set_material("ss316l", density=7990 * u.kg / u.m**3, heat_capacity=500 * u.J / (u.kg * u.K),
                          thermal_conductivity=16.3 * u.W / (u.m * u.K))
    sim.wall.thickness = 3.1 * u.mm
    sim.wall.set_boundary("zero_flux")
    layer = sim.add_layer("Lewatit", composition={ads.name: 1.0}, momentum="ergun", length=5 * u.cm, cells=20,
                          bed_void_fraction=0.4, particle_diameter=0.6 * u.mm)
    layer.set_internal_heat_transfer("constant", internal_htc=0 * u.W / (u.m**2 * u.K))
    layer.set_external_heat_transfer("constant", external_htc=0 * u.W / (u.m**2 * u.K))
    sim.energy_balance.energy_balance_mode = "isothermal"

    sim.initial_conditions.pressure = 1 * u.bar
    sim.initial_conditions.temperature = 298.15 * u.K
    sim.initial_conditions.wall_temperature = 298.15 * u.K
    sim.initial_conditions.velocity = 0 * u.m / u.s
    sim.initial_conditions.y["Nitrogen"] = 1.0
    sim.initial_conditions.y["CarbonDioxide"] = 0.0

    sim.cycle.number_of_cycles = 1
    feed = sim.add_step("Feed", time=6e4 * u.s)
    feed.bottom.set_streams("inlet", [
        sk.Stream("Air", flow=sk.velocity(0.1 * u.m / u.s), y={"Nitrogen": 1 - FEED_CO2, "CarbonDioxide": FEED_CO2},
                  temperature=298.15 * u.K)])
    feed.top.set_pressure("outlet", 1 * u.bar)
    sim.kpi.track("CarbonDioxide", feed=[feed.bottom], product=[feed.top])
    return sim


DECKS = {"DAC": "Illustrative_DAC_Breakthrough"}       # the shipped deck


def main(argv=None) -> int:
    args = c.parser(__doc__.splitlines()[0], ["DAC"]).parse_args(argv)
    sim = build()
    c.run(sim, "Illustrative DAC breakthrough", args.quiet)
    y = sim.series("mole_fraction", "CarbonDioxide", location="top")
    frac = np.asarray(y.values) / FEED_CO2
    for level in (0.05, 0.5, 0.95):
        i = int(np.argmax(frac >= level)) if (frac >= level).any() else None
        print(f"  outlet reaches {level:4.0%} of the feed at " + (f"{y.time[i]:.0f} s" if i is not None else "-"))
    plt = c.pyplot(args)
    if plt:
        fig, ax = plt.subplots(figsize=(6, 4))
        ax.plot(np.asarray(y.time) / 3600, frac)
        ax.set(xlabel="time [h]", ylabel="outlet CO2 / feed CO2")
        c.save(fig, "illustrative_dac_breakthrough.png")
    return 0


if __name__ == "__main__":
    sys.exit(main())
