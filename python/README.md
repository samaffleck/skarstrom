# skarstrom

Python interface to [Skarstrom](https://github.com/samaffleck/skarstrom), an
open-source simulator for pressure, vacuum and temperature swing adsorption
(PSA / VSA / TSA) cycles. The package wraps the C++ engine: load a deck (a
JSON simulation file), change parameters, run it to cyclic steady state and
read KPIs, time series and axial profiles as NumPy arrays.

```python
import skarstrom

sim = skarstrom.examples.Perez2019_A.load()          # or skarstrom.Simulation("deck.json")
sim.cycle.number_of_cycles = 10
sim.steps["Feed"].time = "2 min"                     # bare numbers are SI; strings carry units
result = sim.run()
print(result.status, result.message)
for k in sim.kpis(cycle=-1):
    print(k.component, k.purity, k.recovery)
y = sim.series("mole_fraction", "CarbonDioxide", location="top")   # y.time, y.values
```

Command line:

```bash
skarstrom examples                                     # the bundled decks
skarstrom examples run "Ward2022_Q100ccm" --json   # run one, JSON summary
skarstrom run deck.json --cycles 20 --set "step[Feed].time=120" --csv out/
skarstrom params deck.json                             # settable parameters
```

The package bundles every example deck of the repository (`skarstrom.examples`)
and the curated adsorbent library (`skarstrom.adsorbent_library()`), and it
exposes the engine's isotherm evaluation and fitter
(`skarstrom.evaluate_isotherm`, `skarstrom.fit_isotherm`).

- Documentation: <https://samaffleck.github.io/skarstrom/>
- Source and issues: <https://github.com/samaffleck/skarstrom>

Skarstrom is licensed under the Apache License 2.0 (see `LICENSE` and
`NOTICE`, installed in `skarstrom/licenses/`). The wheels also contain the
third-party components listed in `skarstrom/licenses/THIRD_PARTY_NOTICES.md`.
