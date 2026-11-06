# Experimental data: Streb & Mazzotti (2021)

Breakthrough curves from Streb & Mazzotti, *Adsorption for efficient low
carbon hydrogen production: part 1 — adsorption equilibrium and breakthrough
studies for H₂/CO₂/CH₄ on zeolite 13X*, Adsorption 27:541–558 (2021),
doi:[10.1007/s10450-021-00306-y](https://doi.org/10.1007/s10450-021-00306-y).
One folder per experiment of the paper's Table 2: the binary CO₂:CH₄ = 50:50
experiments A1–A7 and the ternary H₂:CO₂:CH₄ = 75:20:5 experiments B1 and B2
(the column starts under helium). Each folder holds the outlet mole fraction
of every fed component and the five centre-line thermocouples at z = 10, 35,
60, 85 and 110 cm (Sect. 2.2).

| Experiment | T [°C] | P [bar] | Flow [cm³/s at P, T] | Feed |
|---|---|---|---|---|
| A1 | 25 | 1 | 20 | CO₂:CH₄ 50:50 |
| A2 | 25 | 10 | 10 | CO₂:CH₄ 50:50 |
| A3 | 25 | 10 | 20 | CO₂:CH₄ 50:50 (reference experiment) |
| A4 | 25 | 10 | 40 | CO₂:CH₄ 50:50 |
| A5 | 25 | 25 | 20 | CO₂:CH₄ 50:50 |
| A6 | 45 | 10 | 20 | CO₂:CH₄ 50:50 |
| A7 | 65 | 10 | 20 | CO₂:CH₄ 50:50 |
| B1 | 25 | 1 | 20 | H₂:CO₂:CH₄ 75:20:5 |
| B2 | 25 | 10 | 10 | H₂:CO₂:CH₄ 75:20:5 |

(Conditions from Table 2 as transcribed in `core/test/papers/streb_2021/README.md`.)

**These are digitised figure points, not measured data files.** No SI table
or data repository with the measurements is recorded. If one exists, it
replaces these files (protocol step 4 in `core/test/papers/README.md`).

## Source figures

| Item | Record |
|---|---|
| Figure | A3: Fig. 6; A1, A2, A4–A7: Fig. 7; B1, B2: Fig. 11 (the labels the test source `streb_2021.cpp` attaches to each comparison; not re-checked against the paper) |
| Panel per experiment | not recorded |
| Figure axis ranges | not recorded |
| Digitising tool and operator | not recorded (the files predate the repository's history, which starts with the open-source release of 2026-10-17) |

## Files

Columns: *Points* counts data rows; *Time range* and *Data range* are the
minimum and maximum of the points in the file (the figure's own axis ranges
were not recorded when the curves were digitised). Every file is two columns,
time in s and the quantity named in the header, comma-separated.

| File | Quantity (header) | Points | Time range [s] | Data range | Note |
|---|---|---|---|---|---|
| `A1/A1_mole_fraction_CH4.csv` | `Time [s], Y_CH4 [-]` | 17 | 0–4042 | 0–1.002 |  |
| `A1/A1_mole_fraction_CO2.csv` | `Time [s], Y_CO2 [-]` | 13 | 0–4041 | 0–0.439 |  |
| `A1/A1_temperature_z=10[cm].csv` | `Time [s], T [K]` | 39 | 0–4041 | 297.3–365.5 | 1 point(s) out of time order (the metrics helpers sort by time) |
| `A1/A1_temperature_z=110[cm].csv` | `Time [s], T [K]` | 32 | 0–4039 | 297.6–359.9 |  |
| `A1/A1_temperature_z=35[cm].csv` | `Time [s], T [K]` | 35 | 0–4047 | 297.4–360.7 |  |
| `A1/A1_temperature_z=60[cm].csv` | `Time [s], T [K]` | 34 | 0–4046 | 297.5–360.3 |  |
| `A1/A1_temperature_z=85[cm].csv` | `Time [s], T [K]` | 39 | 0–4043 | 297.6–359.4 |  |
| `A2/A2_mole_fraction_CH4.csv` | `Time [s], Y_CH4 [-]` | 27 | 0–1451 | 0–1.003 |  |
| `A2/A2_mole_fraction_CO2.csv` | `Time [s], Y_CO2 [-]` | 21 | 0–1451 | 0–0.5016 |  |
| `A2/A2_temperature_z=10[cm].csv` | `Time [s], T [K]` | 41 | 0–1449 | 298–398.9 |  |
| `A2/A2_temperature_z=110[cm].csv` | `Time [s], T [K]` | 39 | 0–1449 | 297.8–379.6 |  |
| `A2/A2_temperature_z=35[cm].csv` | `Time [s], T [K]` | 36 | 0–1452 | 298–390.2 |  |
| `A2/A2_temperature_z=60[cm].csv` | `Time [s], T [K]` | 39 | 0–1451 | 297.8–385 |  |
| `A2/A2_temperature_z=85[cm].csv` | `Time [s], T [K]` | 38 | 0–1449 | 297.7–381.8 |  |
| `A3/A3_mole_fraction_CH4.csv` | `Time [s], Y_CH4 [-]` | 33 | 0–1181 | 0–1 |  |
| `A3/A3_mole_fraction_CO2.csv` | `Time [s], Y_CO2 [-]` | 26 | 0–1179 | 0–0.4993 |  |
| `A3/A3_temperature_z=10[cm].csv` | `Time [s], T [K]` | 39 | 0–1180 | 297.4–397.7 |  |
| `A3/A3_temperature_z=110[cm].csv` | `Time [s], T [K]` | 37 | 0–1180 | 298–383 |  |
| `A3/A3_temperature_z=35[cm].csv` | `Time [s], T [K]` | 42 | 0–1181 | 297.8–393.4 |  |
| `A3/A3_temperature_z=60[cm].csv` | `Time [s], T [K]` | 42 | 0–1180 | 297.6–389.2 |  |
| `A3/A3_temperature_z=85[cm].csv` | `Time [s], T [K]` | 36 | 0–1180 | 297.8–385.6 |  |
| `A4/A4_mole_fraction_CH4.csv` | `Time [s], Y_CH4 [-]` | 31 | 0–920.1 | 0–1.002 |  |
| `A4/A4_mole_fraction_CO2.csv` | `Time [s], Y_CO2 [-]` | 19 | 0–920.1 | 0–0.5025 |  |
| `A4/A4_temperature_z=10[cm].csv` | `Time [s], T [K]` | 30 | 0–920.1 | 297.9–391 |  |
| `A4/A4_temperature_z=110[cm].csv` | `Time [s], T [K]` | 41 | 0–920.2 | 298.5–381.6 |  |
| `A4/A4_temperature_z=35[cm].csv` | `Time [s], T [K]` | 39 | 0–920.2 | 297.9–390.3 |  |
| `A4/A4_temperature_z=60[cm].csv` | `Time [s], T [K]` | 45 | 0–920.6 | 297.9–387.8 |  |
| `A4/A4_temperature_z=85[cm].csv` | `Time [s], T [K]` | 41 | 0–921.6 | 298.2–384.7 |  |
| `A5/A5_mole_fraction_CH4.csv` | `Time [s], Y_CH4 [-]` | 34 | 0–694.7 | 0–1.001 |  |
| `A5/A5_mole_fraction_CO2.csv` | `Time [s], Y_CO2` | 19 | 0–694.6 | 0–0.5111 |  |
| `A5/A5_temperature_z=10[cm].csv` | `Time [s], T [K]` | 27 | 0–694.3 | 297.3–395.7 |  |
| `A5/A5_temperature_z=110[cm].csv` | `Time [s], T [K]` | 42 | 0–692.6 | 297.8–386.4 |  |
| `A5/A5_temperature_z=35[cm].csv` | `Time [s], T [K]` | 32 | 0–694.5 | 297.3–396.8 |  |
| `A5/A5_temperature_z=60[cm].csv` | `Time [s], T [K]` | 35 | 0–694.4 | 297.3–393.6 |  |
| `A5/A5_temperature_z=85[cm].csv` | `Time [s], T [K]` | 32 | 0–694.5 | 297.6–389.7 |  |
| `A6/A6_mole_fraction_CH4.csv` | `Time [s], Y_CH4 [-]` | 26 | 0–1005 | 0–1.001 |  |
| `A6/A6_mole_fraction_CO2.csv` | `Time [s], Y_CO2 [-]` | 19 | 0–1005 | 0–0.4986 |  |
| `A6/A6_temperature_z=10[cm].csv` | `Time [s], T [K]` | 32 | 0–1005 | 309–410.7 |  |
| `A6/A6_temperature_z=110[cm].csv` | `Time [s], T [K]` | 29 | 0–1006 | 317.7–397.4 |  |
| `A6/A6_temperature_z=35[cm].csv` | `Time [s], T [K]` | 38 | 0–1006 | 316.8–406.5 |  |
| `A6/A6_temperature_z=60[cm].csv` | `Time [s], T [K]` | 34 | 0–1005 | 318.1–402.5 |  |
| `A6/A6_temperature_z=85[cm].csv` | `Time [s], T [K]` | 33 | 0–1005 | 318–399.8 |  |
| `A7/A7_mole_fraction_CH4.csv` | `Time [s], Y_CH4` | 28 | 0–1061 | 0–1.001 |  |
| `A7/A7_mole_fraction_CO2.csv` | `Time [s], Y_CO2` | 17 | 0–1060 | 0–0.4962 |  |
| `A7/A7_temperature_z=10[cm].csv` | `Time [s], T [K]` | 29 | 0–1059 | 319.9–423 |  |
| `A7/A7_temperature_z=110[cm].csv` | `Time [s], T [K]` | 26 | 0–1059 | 337.5–412.6 |  |
| `A7/A7_temperature_z=35[cm].csv` | `Time [s], T [K]` | 35 | 0–1059 | 335–419.8 |  |
| `A7/A7_temperature_z=60[cm].csv` | `Time [s], T [K]` | 31 | 0–1060 | 337.6–416.4 |  |
| `A7/A7_temperature_z=85[cm].csv` | `Time [s], T [K]` | 33 | 0–1059 | 337.5–413.9 |  |
| `B1/B1_mole_fraction_CH4.csv` | `Time [s], Y_CH4 [-]` | 13 | 0–1.034e+04 | 0–0.0633 |  |
| `B1/B1_mole_fraction_CO2.csv` | `Time [s], Y_CO2 [-]` | 14 | 0–1.035e+04 | 0–0.2023 |  |
| `B1/B1_mole_fraction_H2.csv` | `Time [s], Y_H2 [-]` | 25 | 0–1.034e+04 | 0–1.001 | 1 point(s) out of time order (the metrics helpers sort by time) |
| `B1/B1_temperature_z=10[cm].csv` | `Time [s], T [K]` | 20 | 0–1.034e+04 | 298.1–341.1 |  |
| `B1/B1_temperature_z=110[cm].csv` | `Time [s], T [K]` | 18 | 0–1.035e+04 | 297.9–340.1 |  |
| `B1/B1_temperature_z=35[cm].csv` | `Time [s], T [K]` | 18 | 0–1.035e+04 | 298.1–340.3 |  |
| `B1/B1_temperature_z=60[cm].csv` | `Time [s], T [K]` | 18 | 0–1.034e+04 | 298.1–340.1 |  |
| `B1/B1_temperature_z=85[cm].csv` | `Time [s], T [K]` | 18 | 0–1.035e+04 | 298.2–340.4 |  |
| `B2/B2_mole_fraction_CH4.csv` | `Time [s], Y_CH4 [-]` | 17 | 0–3310 | 0–0.08164 |  |
| `B2/B2_mole_fraction_CO2.csv` | `Time [s], Y_CO2 [-]` | 10 | 0–3308 | 0–0.1992 |  |
| `B2/B2_mole_fraction_H2.csv` | `Time [s], Y_H2 [-]` | 24 | 0–3309 | 0–1.001 |  |
| `B2/B2_temperature_z=10[cm].csv` | `Time [s], T [K]` | 30 | 0–3308 | 297.4–385.4 |  |
| `B2/B2_temperature_z=110[cm].csv` | `Time [s], T [K]` | 30 | 0–3309 | 297.9–382.5 |  |
| `B2/B2_temperature_z=35[cm].csv` | `Time [s], T [K]` | 34 | 0–3308 | 297.7–382.9 |  |
| `B2/B2_temperature_z=60[cm].csv` | `Time [s], T [K]` | 35 | 0–3308 | 297.7–382.7 |  |
| `B2/B2_temperature_z=85[cm].csv` | `Time [s], T [K]` | 31 | 0–3307 | 298–382.4 |  |

## Use

- `core/test/papers/streb_2021/streb_2021.cpp`: outlet curve NRMSE (normalised
  by the feed mole fraction) and front time at half the feed fraction for CO₂
  and CH₄ (and the H₂ curve for B1, B2); thermocouple RMSE pooled over the
  five positions.
- Embedded as experimental data in the bundled `Streb2021_A1` and
  `Streb2021_A3` decks (`examples/decks/`).

## Known caveats

- **Sampling line.** The measured compositions are taken at the mass
  spectrometer after about 10.5 m of sampling line (Table 1), which the paper
  models and the decks do not; the curves therefore include a transport delay
  of order seconds to tens of seconds.
- **Thermocouples on the centre line.** The paper notes that the centre line
  cools more slowly than its 1D model predicts.
- **Out-of-order points.** Two files have a point earlier in time than the
  one before it: `A1/A1_temperature_z=10[cm].csv` (rows 15–16, near 260 s on
  the steep rise) and `B1/B1_mole_fraction_H2.csv` (rows 3–4, near 58 s).
  The validation helpers sort the points by time before use.
- Header spelling differs between files (`Y_CO2` without `[-]` in A5 and A7,
  `Y_CH4` in A7); the loaders do not read the header.
- Values slightly above the physical limit (y up to 1.003) are digitising
  error.
