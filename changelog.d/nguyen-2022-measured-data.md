### Fixed — Nguyen 2022: compare with the measured data, not the paper's model
- The files in `core/data/nguyen_2022/Fig4/` and `Fig5/` were the paper's model curves, digitised; the paper test and the `Nguyen2022_Fig4` and `Nguyen2022_Fig5` decks used them as the experiment. They are now `core/data/nguyen_2022/paper_model/` (the test's "vs paper model" check). The measured markers are re-digitised from the paper's figures (`core/data/nguyen_2022/digitise_figures.py`), including the outlet N2 of Fig. 5.
- The decks follow the paper and the author's thesis (Table A2.1): Fig. 5 feeds CO2 in nitrogen (the decks fed CO2 in helium), each run has its own temperature and pressure, and the decks use the constant gas heat capacities of Table 2, the capacity-factor LDF of Eq. 5 with D_m = 1.496e-5 m2/s, Eq. 6 dispersion with the same D_m, the Blake–Kozeny pressure drop and no gas holdup in the pores.
- The measured thermocouple 1.82 cm from the outlet is now compared too.
- Nothing is fitted. Against the measurement, the CO2 fronts are within 3–11 % (the paper's own model: 5–15 %). The Fig. 5 thermocouple RMSE is 1.0–2.7 K. The Fig. 4 temperature peaks are 6–8 K too high; this is an open item. `core/test/papers/nguyen_2022/README.md` has the provenance, the attribution and the open items.

### Added — Nguyen 2022: all six breakthrough runs as example decks
- `Nguyen2022_Fig4_50pct`, `Nguyen2022_Fig4_100pct`, `Nguyen2022_Fig5_5pct` and `Nguyen2022_Fig5_30pct` join `Nguyen2022_Fig4` (70 % CO2) and `Nguyen2022_Fig5` (15 % CO2). `examples/nguyen_2022.py` builds all six.
- The test cases are renamed `Nguyen2022Test.Fig4_100pct` … `Fig5_5pct`.

### Fixed — registry coverage test across map entries
- `RegistryCoverage` intersects its three passes per map entry (one component's value at a time), not per generic path. A deck whose adsorbent mixed a capacity-factor component with a lumped one was reported as uncovered, although each value is covered in one of the passes.
