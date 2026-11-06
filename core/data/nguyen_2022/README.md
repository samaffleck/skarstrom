# Experimental data: Nguyen et al. (2022)

Breakthrough curves from Nguyen, Lin, Shimizu, Rajendran, *Separation of CO₂
and N₂ on a hydrophobic metal organic framework CALF-20*, Chemical
Engineering Journal 442:136263 (2022),
doi:[10.1016/j.cej.2022.136263](https://doi.org/10.1016/j.cej.2022.136263),
on a column packed with CALF-20:

- `Fig4/`: binary CO₂/He breakthrough at feed CO₂ fractions of 50, 70 and
  100 % (file suffix = feed fraction).
- `Fig5/`: breakthrough at feed CO₂ fractions of 5, 15 and 30 %, described as
  CO₂/N₂/He experiments (the decks feed CO₂ in He; see the paper README).

Each feed fraction has an outlet CO₂ mole fraction curve (`Y_CO2_*`) and a
temperature curve (`T_CO2_*`).

**These are digitised figure points, not measured data files.** The paper
text has never been in the repository, so whether the paper's SI or a data
repository holds the measurements has not been checked. If they exist, they
replace these files (protocol step 4 in `core/test/papers/README.md`).

## Source figures

| Item | Record |
|---|---|
| Figure | Fig. 4 and Fig. 5 (folder names and the labels in `nguyen_2022.cpp`; not re-checked against the paper) |
| Panel per file | not recorded |
| Figure axis ranges | not recorded |
| Feed flow per curve | not recorded; the tests assume 50 cm³/min for every Fig. 4 curve and 100 cm³/min for every Fig. 5 curve (the flows of the two bundled cases) |
| Thermocouple position of `T_CO2_*` | not recorded (not stated in the material available to the repository) |
| Digitising tool and operator | not recorded (the files predate the repository's history, which starts with the open-source release of 2026-10-17) |

## Files

Columns: *Points* counts data rows; *Time range* and *Data range* are the
minimum and maximum of the points in the file (the figure's own axis ranges
were not recorded when the curves were digitised). Every file is two columns,
time in s and the quantity named in the header, comma-separated.

| File | Quantity (header) | Points | Time range [s] | Data range | Note |
|---|---|---|---|---|---|
| `Fig4/T_CO2_100%.csv` | `Time [s], T [K]` | 32 | 0–7069 | 295–327.3 | UTF-8 BOM before the header |
| `Fig4/T_CO2_50%.csv` | `Time [s], T [K]` | 30 | 0–7069 | 295.4–315 | UTF-8 BOM before the header |
| `Fig4/T_CO2_70%.csv` | `Time [s], T [K]` | 35 | 0–7069 | 295.5–320.6 | UTF-8 BOM before the header |
| `Fig4/Y_CO2_100%.csv` | `Time [s], Y_CO2 [-]` | 15 | 0–7069 | 0–1 | UTF-8 BOM before the header |
| `Fig4/Y_CO2_50%.csv` | `Time [s], Y_CO2 [-]` | 24 | 0–7069 | 0–0.4984 | UTF-8 BOM before the header |
| `Fig4/Y_CO2_70%.csv` | `Time [s], Y_CO2 [-]` | 22 | 0–7069 | 0–0.7091 | UTF-8 BOM before the header |
| `Fig5/T_CO2_15%.csv` | `Time [s], T [K]` | 37 | 0–1.178e+04 | 293.8–309.9 |  |
| `Fig5/T_CO2_30%.csv` | `Time [s], T [K]` | 49 | 0–1.178e+04 | 295–322.9 |  |
| `Fig5/T_CO2_5%.csv` | `Time [s], T [K]` | 23 | 0–1.178e+04 | 296.7–300.1 |  |
| `Fig5/Y_CO2_15%.csv` | `Time [s], Y_CO2 [-]` | 31 | 0–1.178e+04 | 0–0.15 |  |
| `Fig5/Y_CO2_30%.csv` | `Time [s], Y_CO2 [-]` | 32 | 0–1.178e+04 | 0–0.3 |  |
| `Fig5/Y_CO2_5%.csv` | `Time [s], Y_CO2 [-]` | 16 | 0–1.176e+04 | 0–0.04649 |  |

## Use

- `core/test/papers/nguyen_2022/nguyen_2022.cpp`: outlet CO₂ curve NRMSE
  (normalised by the feed fraction) and front time at half the feed fraction,
  for all six `Y_CO2_*` curves. The temperature curves are not scored because
  the thermocouple position is unknown.
- The initial temperature of the `T_CO2_*` curves (294.8–296.7 K) is the
  basis of the decks' assumed 295 K ambient temperature (`nguyen_2022.h`).
- Embedded as experimental data in the bundled `Nguyen2022_Fig4` and
  `Nguyen2022_Fig5` decks (`examples/decks/`), all files of the folder.

## Known caveats

- The Fig. 4 files start with a UTF-8 byte-order mark; the loaders strip it.
- `Fig5/Y_CO2_5%.csv` ends at y = 0.0465, below the 0.05 feed fraction: the
  digitised record stops before the plateau is reached.
- `Fig5/T_CO2_15%.csv` falls to 293.8 K near t = 11 100 s, 1 K below its starting value;
  whether that is in the figure or a digitising error is not recorded.
