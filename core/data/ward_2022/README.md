# Experimental data: Ward & Pini (2022)

Curves for the three single-component CO₂ breakthrough experiments used by
Ward & Pini, *Integrated uncertainty quantification and sensitivity analysis
of single-component dynamic column breakthrough experiments*, Adsorption
28:161–183 (2022), doi:[10.1007/s10450-022-00361-z](https://doi.org/10.1007/s10450-022-00361-z).
The measurements are those of the paper's ref. [5] (Pini, Joss, Hejazi 2021):
pure CO₂ fed at 293 K and 1 bar into a column of activated carbon initially
under helium, outlet CO₂ mole fraction, outlet volumetric flow and one
thermocouple at z = L/4 = 0.0715 m.

**These are digitised figure points, not measured data files.** The original
measurements of ref. [5] were not available to the repository. If they become
available (an SI table or a data repository), they replace these files
(protocol step 4 in `core/test/papers/README.md`).

## Source figure

| Item | Record |
|---|---|
| Figure | Fig. 2 according to the first version of this README (2026-10-31); the test source (`ward_2022.cpp`) labels the comparison "Fig. 3 / data of ref. [5]". Which is right has not been re-checked against the paper. |
| Panel per file | not recorded |
| Figure axis ranges | not recorded |
| Digitising tool and operator | not recorded (the files predate the repository's history, which starts with the open-source release of 2026-10-17) |

One folder per feed flow (`Qin=50`, `Qin=100`, `Qin=200` cm³/min).

## Files

Columns: *Points* counts data rows; *Time range* and *Data range* are the
minimum and maximum of the points in the file (the figure's own axis ranges
were not recorded when the curves were digitised). Every file is two columns,
time in s and the quantity named in the header, comma-separated.

| File | Quantity (header) | Points | Time range [s] | Data range | Note |
|---|---|---|---|---|---|
| `Qin=100/outlet_CO2_mole_fraction[-]_Qin=100[ccm].csv` | `Time [s], y_CO2 [-]` | 29 | 0–3000 | 0–0.9961 |  |
| `Qin=100/outlet_volumetric_flowrate_Qin=100[ccm].csv` | `Time [s], Outlet volumetric flowrate [ccm]` | 26 | 0–3002 | 3.025e-08–1.686e-06 | values are in m³/s, not the cm³/min of the header: the plateau 1.686e-06 m³/s = 101.2 cm³/min |
| `Qin=100/temperature[K]_z=0.0715[m]_Qin=100[ccm].csv` | `Time [s], Temperature [K]` | 37 | 0–2996 | 291.4–306.1 |  |
| `Qin=200/outlet_CO2_mole_fraction[-]_Qin=200[ccm].csv` | `Time [s], y_CO2 [-]` | 28 | 0–2992 | 0–0.9974 |  |
| `Qin=200/outlet_volumetric_flowrate_Qin=200[ccm].csv` | `Time [s], Outlet volumetric flowrate [ccm]` | 30 | 0–2998 | 8.2e-08–3.333e-06 | values are in m³/s, not the cm³/min of the header: the plateau 3.333e-06 m³/s = 200 cm³/min |
| `Qin=200/temperature[K]_z=0.0715[m]_Qin=200[ccm].csv` | `Time [s], Temperature [K]` | 42 | 0–3008 | 291.7–313.3 |  |
| `Qin=50/outlet_CO2_mole_fraction[-]_Qin=50[ccm].csv` | `Time [s], y_CO2 [-]` | 27 | 0–2994 | 0–0.9961 |  |
| `Qin=50/outlet_volumetric_flowrate_Qin=50[ccm].csv` | `Time [s], Outlet volumetric Flowrate [ccm]` | 23 | 0–2998 | 2.163e-08–8.323e-07 | values are in m³/s, not the cm³/min of the header: the plateau 8.323e-07 m³/s = 49.94 cm³/min |
| `Qin=50/temperature[K]_z=0.0715[m]_Qin=50[ccm].csv` | `Time [s], Temperature [K]` | 31 | 0–3000 | 291.7–299.4 |  |

## Use

- `core/test/papers/ward_2022/ward_2022.cpp`: outlet CO₂ curve NRMSE, front
  time at y_out = 0.5 (Table 9's t_br) and thermocouple RMSE. The outlet-flow
  series is not compared.
- Embedded as experimental data in the bundled `Ward2022_Q50ccm`,
  `Ward2022_Q100ccm` and `Ward2022_Q200ccm` decks (`examples/decks/`), all
  three series per folder.

## Known caveats

- **Outlet-flow unit.** The header says `[ccm]` but the values are m³/s
  (plateaus of 8.32 × 10⁻⁷, 1.69 × 10⁻⁶ and 3.33 × 10⁻⁶ m³/s are 49.9, 101
  and 200 cm³/min). The experimental-data loader does not parse units from
  headers, so the decks show the mislabelled header with SI values.
- **Front timing.** A digitised point on the steep part of a front carries a
  time uncertainty of a few seconds, below the differences the validation
  reports.
- The first outlet-flow points (about 10⁻⁸ m³/s, before breakthrough) are
  non-zero digitising offsets, not a measured flow.
