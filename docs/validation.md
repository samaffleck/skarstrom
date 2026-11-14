<!-- GENERATED — do not edit. Built by scripts/docs/gen_validation_page.py from
     core/test/papers/validation_summary.json, which the paper tests in core/test/papers/ write.
     Regenerate: run the paper tests (ctest), then
     python scripts/docs/gen_validation_page.py --merge <build>/validation -->

# Validation

Skarstrom is checked against published adsorption studies by automated tests. Every comparison on this page is produced by a test in [`core/test/papers/`](https://github.com/samaffleck/skarstrom/tree/main/core/test/papers) that runs the bundled input deck, computes an explicit error metric against numbers transcribed from the paper (`reference_values.csv` in each paper folder) or digitised from its figures (`core/data/`), and writes the result to a machine-readable summary. This page is generated from that summary — it is never edited by hand.

**How to read the tables.**

- **Against** says what the reference is: an *experiment* (measured data — this is validation), the *paper model* (the authors' own simulation — a code-to-code comparison), or an isotherm fit.
- **NRMSE** is the root-mean-square difference between the simulated and the digitised outlet mole fraction, divided by the feed mole fraction of that component. **RMSE** on temperature is in kelvin. Purity and recovery errors are absolute, in percentage points (pp); productivity and breakthrough-time errors are relative.
- **Tolerance** is a regression bound set from the measured error of the current model, rounded up — not an acceptance criterion chosen in advance. A case "passes" when it reproduces its recorded error; a large tolerance means a large, documented disagreement, not a good fit. Rows marked *info* are recorded but not asserted (e.g. approximate values read from the paper text).
- Most input decks are *predictions* with the paper's published parameters. Two are calibrated, and their findings say so: Ward 2022 (fitted to the curves it is compared with) and Streb 2021 (calibrated on experiment A3, with a sampling line and per-experiment feed flows). Known causes of disagreement are listed per paper below and in each folder's README.

## Summary

68 cases from 13 publications; 255 asserted comparisons (0 outside their recorded tolerance) and 130 informational ones. Summary generated on 2026-11-14.

| Publication | Cases | Against | Comparisons | Median error | Worst error |
|---|---|---|---|---|---|
| [Arora & Hasan (2021)](#arora-hasan-2021) | 2 | paper model | 4 | 5.35 pp (abs. error) | 46.6 pp (abs. error) |
| [Casas et al. (2012) / Schell et al. (2012)](#casas-2012) | 4 | paper model, isotherm fit | 8 | 3.35 % (rel. error)<br/>0.595 % (max rel. error) | 16.8 % (rel. error)<br/>6.99 % (max rel. error) |
| [Ferré et al. (2026)](#ferre-2026) | 4 | paper model | 20 | 0.787 pp (abs. error)<br/>6.74 % (rel. error) | 11.5 pp (abs. error)<br/>59.4 % (rel. error) |
| [Haghpanah et al. (2013) / Glover et al. (2026)](#haghpanah-2013) | 7 | reference | 23 | 4.10e-4 (abs. error)<br/>0.61 pp (abs. error)<br/>2.77 % (rel. error) | 4.10e-4 (abs. error)<br/>1.59 pp (abs. error)<br/>4.11 % (rel. error) |
| [Nguyen et al. (2022)](#nguyen-2022) | 6 | experiment, paper model | 33 | 11.8 % (NRMSE)<br/>4.01 % (rel. error)<br/>1.72 K (RMSE) | 29.4 % (NRMSE)<br/>10.9 % (rel. error)<br/>7.97 K (RMSE) |
| [Nguyen et al. (2023)](#nguyen-2023) | 10 | experiment, paper model | 60 | 2.05 pp (abs. error)<br/>8.42 % (rel. error) | 8.39 pp (abs. error)<br/>15.9 % (rel. error) |
| [Nguyen et al. (2025)](#nguyen-2024) | 2 | experiment | 2 | 0.931 pp (abs. error)<br/>2.01 % (rel. error) | 0.931 pp (abs. error)<br/>2.01 % (rel. error) |
| [Estupiñan Perez et al. (2019)](#perez-2019) | 6 | experiment, paper model | 24 | 2.02 pp (abs. error) | 5.81 pp (abs. error) |
| [Streb & Mazzotti (2021)](#streb-2021) | 18 | experiment | 47 | 10.5 % (NRMSE)<br/>0.691 % (rel. error)<br/>7.37 K (RMSE) | 31.2 % (NRMSE)<br/>4.7 % (rel. error)<br/>11.4 K (RMSE) |
| [Vysyaraju et al. (2024)](#vysyarajua-2024) | 3 | experiment, isotherm (derived), paper model | 7 | 1.89 % (rel. error)<br/>0.251 K (abs. error)<br/>1.93 pp (abs. error) | 3.92 % (rel. error)<br/>0.251 K (abs. error)<br/>33.5 pp (abs. error) |
| [Ward & Pini (2022)](#ward-2022) | 3 | experiment | 9 | 2.76 % (NRMSE)<br/>0.0209 % (rel. error)<br/>0.675 K (RMSE) | 3.89 % (NRMSE)<br/>0.404 % (rel. error)<br/>1.4 K (RMSE) |
| [Ward & Pini (2024)](#ward-pini-2024) | 1 | paper model | 3 | 2.34 pp (abs. error)<br/>1.86 % (rel. error) | 2.83 pp (abs. error)<br/>1.86 % (rel. error) |
| [Halabi (2011)](#halabi-2011) | 2 | paper model | 15 | 0.61 pp (abs. error) | 1.95 pp (abs. error) |

## Arora & Hasan (2021) {#arora-hasan-2021}

A. Arora, M.M.F. Hasan, Flexible oxygen concentrators for medical applications, Scientific Reports 11 (2021) 14317. [doi:10.1038/s41598-021-93796-3](https://doi.org/10.1038/s41598-021-93796-3)

*System:* O2/N2 (air) on LiLSX, 4-step PSA and PVSA. *Compared against:* paper's optimised model results (Table 4).

**Findings.**

- O2 purity is within about 4 percentage points of Table 4, but recovery is far below the paper (4 % vs 30 % for PSA, 20 % vs 66 % for PVSA). This is an open disagreement, not a tuned result.
- The paper's recovery is a model value from an optimiser; the t1/t4 rows of the PSA column of Table 4 were lost in the text extraction (t1 = 3.0 s is an estimate). The isotherm form is the non-competitive dual-site Langmuir; the paper does not state whether its DSL is competitive, so it was left as is.
- The PSA deck now vents the depressurisation step counter-currently (feed end), as the paper states for every Table 4 case; it previously vented co-currently.

| Case | Quantity | Against | Simulated | Reference | Error | Tolerance | Source |
|---|---|---|---|---|---|---|---|
| LiLSX_PSA | O2 purity at cycle 50 | paper model | 90.07 % | 84.72 % | 5.35 pp abs. error | 6.5 pp | Table 4, PSA cycle with LiLSX adsorbent, column PC_O2 &gt;= 5 L/min, P_O2 &gt;= 85% |
| LiLSX_PSA | O2 recovery at cycle 50 | paper model | 7.328 % | 30.14 % | 22.8 pp abs. error | 30 pp | Table 4, PSA cycle with LiLSX adsorbent, column PC_O2 &gt;= 5 L/min, P_O2 &gt;= 85% |
| LiLSX_PVSA | O2 purity at cycle 50 | paper model | 93.84 % | 89.97 % | 3.87 pp abs. error | 4.7 pp | Table 4, PVSA cycle with LiLSX adsorbent, column PC_O2 &gt;= 5 L/min, P_O2 &gt;= 90% |
| LiLSX_PVSA | O2 recovery at cycle 50 | paper model | 18.94 % | 65.51 % | 46.6 pp abs. error | 53 pp | Table 4, PVSA cycle with LiLSX adsorbent, column PC_O2 &gt;= 5 L/min, P_O2 &gt;= 90% |

??? note "Notes on these comparisons"

    - paper's optimised model result (simulation, not experiment)
    - paper's optimised model result; recovery net of the O2 product recycled to purge and pressurisation

Input decks, assumptions and known discrepancies: [`core/test/papers/arora_hasan_2021/README.md`](https://github.com/samaffleck/skarstrom/blob/main/core/test/papers/arora_hasan_2021/README.md).

## Casas et al. (2012) / Schell et al. (2012) {#casas-2012}

N. Casas, J. Schell, R. Pini, M. Mazzotti, Fixed bed adsorption of CO2/H2 mixtures on activated carbon: experiments and modeling, Adsorption 18 (2012) 143-161; J. Schell, N. Casas, R. Pini, M. Mazzotti, Pure and binary adsorption of CO2, H2, and N2 on activated carbon, Adsorption 18 (2012) 49-65. [doi:10.1007/s10450-012-9389-z](https://doi.org/10.1007/s10450-012-9389-z); [doi:10.1007/s10450-011-9382-y](https://doi.org/10.1007/s10450-011-9382-y)

*System:* CO2/H2 breakthrough on activated carbon AP3-60. *Compared against:* paper's simulated breakthrough times (Table 5); isotherm fits (Schell Table 3).

**Findings.**

- CO2 breakthrough times agree with the paper's simulated Table 5 values to 0.5-1.3 % for the 50/50 and 25/75 feeds and 13 % for 75/25; H2 breakthrough is 15-19 % early for the two CO2-lean feeds.
- Table 5 is the paper's simulation (feed velocity fitted per experiment, piping to the mass spectrometer modelled); adsim uses the 10 cm3/s set-point and no piping. The experimental curves (Figs. 4, 6) have not been digitised.
- The competitive Arrhenius-Sips isotherm reproduces the paper's per-temperature Sips fits (Schell et al., Table 3) to 7 % (CO2) and 0.6 % (H2).

| Case | Quantity | Against | Simulated | Reference | Error | Tolerance | Source |
|---|---|---|---|---|---|---|---|
| Fig4_Reference_50_50_15bar_25C | H2 breakthrough time (y &gt; 0.002) | paper model | 101.5 s | 105 s | 3.35 % rel. error | 4.1 % | Casas et al. 2012, Table 5 (50/50, 25 C, 15 bar) |
| Fig4_Reference_50_50_15bar_25C | CO2 breakthrough time (y &gt; 0.002) | paper model | 449.4 s | 440 s | 2.15 % rel. error | 2.4 % | Casas et al. 2012, Table 5 (50/50, 25 C, 15 bar) |
| Fig4_Reference_50_50_15bar_25C | H2 breakthrough time vs experiment (approx.) | experiment | 101.5 s | 110 s | 7.75 % rel. error | info | Casas et al. 2012, Sect. 4.1 text (Fig. 4) |
| Fig6_25_75_15bar_25C | H2 breakthrough time (y &gt; 0.002) | paper model | 75 s | 80 s | 6.25 % rel. error | 7.6 % | Casas et al. 2012, Table 5 (25/75, 25 C, 15 bar) |
| Fig6_25_75_15bar_25C | CO2 breakthrough time (y &gt; 0.002) | paper model | 627.7 s | 610 s | 2.91 % rel. error | 3.6 % | Casas et al. 2012, Table 5 (25/75, 25 C, 15 bar) |
| Fig6_75_25_15bar_25C | H2 breakthrough time (y &gt; 0.002) | paper model | 151.9 s | 130 s | 16.8 % rel. error | 21 % | Casas et al. 2012, Table 5 (75/25, 25 C, 15 bar) |
| Fig6_75_25_15bar_25C | CO2 breakthrough time (y &gt; 0.002) | paper model | 360.3 s | 320 s | 12.6 % rel. error | 15 % | Casas et al. 2012, Table 5 (75/25, 25 C, 15 bar) |
| Isotherm_Table3 | CO2 loading vs Table 3 Sips (25-140 C, 0.5-40 bar), max deviation | isotherm fit | — | — | 6.99 % max rel. error | 8 % | Schell et al. 2012, Table 3 |
| Isotherm_Table3 | H2 loading vs Table 3 Sips (25-65 C, 0.5-40 bar), max deviation | isotherm fit | — | — | 0.595 % max rel. error | 2 % | Schell et al. 2012, Table 3 |

??? note "Notes on these comparisons"

    - paper's SIMULATED breakthrough time; Sect. 4.4 defines breakthrough as the time the outlet mole fraction first exceeds 0.002. The paper's model also simulates the piping between column and MS (Sect. 3.2); adsim does not
    - paper's simulated breakthrough time (y &gt; 0.002), see note above
    - approximate experimental value: 'the temperature peak at 110 cm is measured at about 110 s, which corresponds to the H2 breakthrough time observed in Fig. 4'; text value is approximate ('about 110 s')
    - ArrheniusSipsIsotherm with the casas_2012.h constants vs the paper's per-temperature Sips fits

Input decks, assumptions and known discrepancies: [`core/test/papers/casas_2012/README.md`](https://github.com/samaffleck/skarstrom/blob/main/core/test/papers/casas_2012/README.md).

## Ferré et al. (2026) {#ferre-2026}

A. Ferré, C. Wedler, A. Ward, J. Burger, R. Pini, C. Petit, Co-production of high-purity carbon dioxide and methane via a single-bed pressure-vacuum swing adsorption process, Separation and Purification Technology 406 (2026) 138821. [doi:10.1016/j.seppur.2026.138821](https://doi.org/10.1016/j.seppur.2026.138821)

*System:* CO2/CH4 (biogas) 4-step PVSA, four adsorbents. *Compared against:* paper's optimised detailed-model results (Table 5).

**Findings.**

- Product purities agree with Table 5 within 0.1-1.8 percentage points and CO2 recovery within 0.6-3 pp; CH4 recovery is 5-7 pp low for the zeolites and 11 pp high for MIL-160(Al), whose CH4 productivity is 59 % high.
- The decks now use the extended (competitive) DSL isotherm the paper states (Sect. 3.1); they previously used the non-competitive form. adsim uses a constant LDF coefficient and a single heat of adsorption per component instead of the paper's loading-dependent forms.
- Specific energy is not compared: the decks drive the pressure steps with pressure boundary conditions, not pump models.

| Case | Quantity | Against | Simulated | Reference | Error | Tolerance | Source |
|---|---|---|---|---|---|---|---|
| MIL160Al | CO2 purity (heavy product) | paper model | 97.19 % | 95 % | 2.19 pp abs. error | 2.2 pp | Table 5, MIL-160(Al), recovery-recovery case, first row |
| MIL160Al | CO2 recovery | paper model | 91.77 % | 94.1 % | 2.33 pp abs. error | 3.5 pp | Table 5, MIL-160(Al), recovery-recovery case, first row |
| MIL160Al | CH4 purity (light product) | paper model | 96.81 % | 97.6 % | 0.787 pp abs. error | 1.4 pp | Table 5, MIL-160(Al), recovery-recovery case, first row |
| MIL160Al | CH4 recovery | paper model | 47.33 % | 35.8 % | 11.5 pp abs. error | 13 pp | Table 5, MIL-160(Al), recovery-recovery case, first row |
| MIL160Al | CH4 productivity | paper model | 1.435 mol/(m3 s) | 0.9 mol/(m3 s) | 59.4 % rel. error | 68 % | Table 5, MIL-160(Al), recovery-recovery case, first row |
| NaY_F | CO2 purity (heavy product) | paper model | 94.77 % | 95.1 % | 0.335 pp abs. error | 1.1 pp | Table 5, NaY (F), recovery-recovery (Re_CO2-Re_CH4) case, first row |
| NaY_F | CO2 recovery | paper model | 97.61 % | 97.8 % | 0.191 pp abs. error | 0.25 pp | Table 5, NaY (F), recovery-recovery case, first row |
| NaY_F | CH4 purity (light product) | paper model | 98.76 % | 98.8 % | 0.0412 pp abs. error | 0.097 pp | Table 5, NaY (F), recovery-recovery case, first row |
| NaY_F | CH4 recovery | paper model | 62.47 % | 66 % | 3.53 pp abs. error | 5.6 pp | Table 5, NaY (F), recovery-recovery case, first row |
| NaY_F | CH4 productivity | paper model | 3.104 mol/(m3 s) | 3.2 mol/(m3 s) | 3.01 % rel. error | 4.6 % | Table 5, NaY (F), recovery-recovery case, first row |
| NaY_W | CO2 purity (heavy product) | paper model | 93.79 % | 95 % | 1.21 pp abs. error | 1.8 pp | Table 5, NaY (W), recovery-recovery case, first row |
| NaY_W | CO2 recovery | paper model | 96.48 % | 95.9 % | 0.578 pp abs. error | 0.89 pp | Table 5, NaY (W), recovery-recovery case, first row |
| NaY_W | CH4 purity (light product) | paper model | 98.15 % | 97.6 % | 0.551 pp abs. error | 0.65 pp | Table 5, NaY (W), recovery-recovery case, first row |
| NaY_W | CH4 recovery | paper model | 59.95 % | 65.4 % | 5.45 pp abs. error | 7.6 pp | Table 5, NaY (W), recovery-recovery case, first row |
| NaY_W | CH4 productivity | paper model | 2.891 mol/(m3 s) | 3.1 mol/(m3 s) | 6.74 % rel. error | 8.5 % | Table 5, NaY (W), recovery-recovery case, first row |
| Zeolite13X | CO2 purity (heavy product) | paper model | 93.66 % | 95 % | 1.34 pp abs. error | 1.5 pp | Table 5, 13X, recovery-recovery case, first row |
| Zeolite13X | CO2 recovery | paper model | 95.64 % | 95.9 % | 0.262 pp abs. error | 0.35 pp | Table 5, 13X, recovery-recovery case, first row |
| Zeolite13X | CH4 purity (light product) | paper model | 97.99 % | 98 % | 0.0139 pp abs. error | 0.076 pp | Table 5, 13X, recovery-recovery case, first row |
| Zeolite13X | CH4 recovery | paper model | 43.63 % | 48.6 % | 4.97 pp abs. error | 5.7 pp | Table 5, 13X, recovery-recovery case, first row |
| Zeolite13X | CH4 productivity | paper model | 1.188 mol/(m3 s) | 1.3 mol/(m3 s) | 8.6 % rel. error | 11 % | Table 5, 13X, recovery-recovery case, first row |

??? note "Notes on these comparisons"

    - paper's detailed-model result
    - CH4 per bed volume per cycle time
    - paper's optimised detailed-model result (simulation)
    - CH4 produced per bed volume per total cycle time (Sect. 4.2)

Input decks, assumptions and known discrepancies: [`core/test/papers/ferre_2026/README.md`](https://github.com/samaffleck/skarstrom/blob/main/core/test/papers/ferre_2026/README.md).

## Haghpanah et al. (2013) / Glover et al. (2026) {#haghpanah-2013}

R. Haghpanah, A. Majumder, R. Nilam, A. Rajendran, S. Farooq, I.A. Karimi, M. Amanullah, Multiobjective optimization of a four-step adsorption process for postcombustion CO2 capture via finite volume simulation, Industrial & Engineering Chemistry Research 52 (2013) 4249-4265; T. Glover, T. Papathanasiou, R. Pini, Accelerating simulation and optimisation of cyclic adsorption processes with differentiable programming, arXiv:2606.21523 (2026). [doi:10.1021/ie302658y](https://doi.org/10.1021/ie302658y)

*System:* CO2/N2 4-step PVSA on zeolite 13X (pressurisation, adsorption, forward blowdown, reverse evacuation). *Compared against:* paper's own model results (Sect. 6.1, Table 6); Glover et al.'s re-implementation of the same model (Table 4) and its CSS multiplier.

**Findings.**

- On the five minimum-energy designs of Table 6, CO2 purity and recovery agree with Glover et al.'s verification runs within 1.6 percentage points and the productivity within 3.4 %. The specific energy is 0.7-6 % above Haghpanah's Table 6 for evacuation pressures of 0.05-0.1 bar but 12 % above at 0.02-0.03 bar, where nearly all of it is vacuum work.
- The energy gap and the Sect. 6.1 purity (79.7 % against 88.74 %; recovery matches within 1.6 pp) are open. Grid, heat-of-adsorption value and convention, wall heat transfer and dispersion closure move them by under 0.3 pp or 1 %; the unit basis of q* in the paper's loading-dependent LDF coefficient (Eq. 40) moves the energy into agreement only at the cost of the Sect. 6.1 recovery (README).
- The dominant multiplier of the cycle map at Glover et al.'s benchmark design, fitted from the decay of the end-of-cycle state change on 10 cells, is 0.9513 against Glover's 0.9509.

| Case | Quantity | Against | Simulated | Reference | Error | Tolerance | Source |
|---|---|---|---|---|---|---|---|
| Fig5Benchmark | Dominant CSS multiplier (fitted decay) | reference | 0.9513 | 0.9509 | 4.10e-4 abs. error | 0.0015 | Glover et al. 2026, p. 27 |
| Fig5Benchmark | Dominant CSS multiplier (engine rate, median) | reference | 0.9513 | 0.9509 | 4.45e-4 abs. error | info | Glover et al. 2026, p. 27 |
| Fig5Benchmark | Cycles to a 1e-5 state change, N2 start | reference | 99 cycles | 115 cycles | 13.9 % rel. error | info | Glover et al. 2026, pp. 23-24 |
| Fig5Benchmark | Cycles to a 1e-5 state change, feed start | reference | 72 cycles | 82 cycles | 12.2 % rel. error | info | Glover et al. 2026, pp. 23-24 |
| PL002 | Energy vs Haghpanah Table 6 | reference | 167.5 kWh/t | 149 kWh/t | 12.4 % rel. error | info | Haghpanah et al. 2013, Table 6, P_L = 0.02 bar |
| PL002 | CO2 purity vs Glover | reference | 90.86 % | 90.6 % | 0.256 pp abs. error | 0.35 pp | Glover et al. 2026, Table 4, this work |
| PL002 | CO2 recovery vs Glover | reference | 87.86 % | 89.4 % | 1.54 pp abs. error | 1.8 pp | Glover et al. 2026, Table 4, this work |
| PL002 | Productivity vs Glover | reference | 0.2888 mol/(m3 s) | 0.297 mol/(m3 s) | 2.77 % rel. error | 3.3 % | Glover et al. 2026, Table 4, this work |
| PL002 | Productivity vs Glover Ref. | reference | 0.2888 mol/(m3 s) | 0.312 mol/(m3 s) | 7.44 % rel. error | info | Glover et al. 2026, Table 4, Ref. |
| PL002 | Energy vs Glover | reference | 167.5 kWh/t | 143.6 kWh/t | 16.6 % rel. error | info | Glover et al. 2026, Table 4, this work |
| PL002 | Productivity vs Haghpanah | reference | 0.4584 mol/(m3 adsorbent s) | 0.49 mol/(m3 adsorbent s) | 6.45 % rel. error | info | Haghpanah et al. 2013, Conclusions, p. 4262 |
| PL003 | Energy vs Haghpanah Table 6 | reference | 225.3 kWh/t | 201.5 kWh/t | 11.8 % rel. error | info | Haghpanah et al. 2013, Table 6, P_L = 0.03 bar |
| PL003 | CO2 purity vs Glover | reference | 90.94 % | 90.6 % | 0.342 pp abs. error | 0.45 pp | Glover et al. 2026, Table 4, this work |
| PL003 | CO2 recovery vs Glover | reference | 89.17 % | 90.7 % | 1.53 pp abs. error | 1.8 pp | Glover et al. 2026, Table 4, this work |
| PL003 | Productivity vs Glover | reference | 0.491 mol/(m3 s) | 0.508 mol/(m3 s) | 3.36 % rel. error | 4 % | Glover et al. 2026, Table 4, this work |
| PL003 | Productivity vs Glover Ref. | reference | 0.491 mol/(m3 s) | 0.512 mol/(m3 s) | 4.11 % rel. error | 4.8 % | Glover et al. 2026, Table 4, Ref. |
| PL003 | Energy vs Glover | reference | 225.3 kWh/t | 194.5 kWh/t | 15.8 % rel. error | info | Glover et al. 2026, Table 4, this work |
| PL005 | Energy vs Haghpanah Table 6 | reference | 317.8 kWh/t | 298.9 kWh/t | 6.33 % rel. error | info | Haghpanah et al. 2013, Table 6, P_L = 0.05 bar |
| PL005 | CO2 purity vs Glover | reference | 89.69 % | 90.3 % | 0.61 pp abs. error | 0.75 pp | Glover et al. 2026, Table 4, this work |
| PL005 | CO2 recovery vs Glover | reference | 90.73 % | 90.9 % | 0.171 pp abs. error | 0.25 pp | Glover et al. 2026, Table 4, this work |
| PL005 | Productivity vs Glover | reference | 0.3305 mol/(m3 s) | 0.333 mol/(m3 s) | 0.746 % rel. error | 1.2 % | Glover et al. 2026, Table 4, this work |
| PL005 | Productivity vs Glover Ref. | reference | 0.3305 mol/(m3 s) | 0.348 mol/(m3 s) | 5.02 % rel. error | info | Glover et al. 2026, Table 4, Ref. |
| PL005 | Energy vs Glover | reference | 317.8 kWh/t | 289.3 kWh/t | 9.87 % rel. error | info | Glover et al. 2026, Table 4, this work |
| PL007 | Energy vs Haghpanah Table 6 | reference | 402.6 kWh/t | 388.7 kWh/t | 3.57 % rel. error | 4.2 % | Haghpanah et al. 2013, Table 6, P_L = 0.07 bar |
| PL007 | CO2 purity vs Glover | reference | 90.57 % | 91 % | 0.431 pp abs. error | 0.55 pp | Glover et al. 2026, Table 4, this work |
| PL007 | CO2 recovery vs Glover | reference | 90.78 % | 91.4 % | 0.615 pp abs. error | 0.75 pp | Glover et al. 2026, Table 4, this work |
| PL007 | Productivity vs Glover | reference | 0.691 mol/(m3 s) | 0.704 mol/(m3 s) | 1.84 % rel. error | 2.3 % | Glover et al. 2026, Table 4, this work |
| PL007 | Productivity vs Glover Ref. | reference | 0.691 mol/(m3 s) | 0.688 mol/(m3 s) | 0.442 % rel. error | 1 % | Glover et al. 2026, Table 4, Ref. |
| PL007 | Energy vs Glover | reference | 402.6 kWh/t | 373.2 kWh/t | 7.87 % rel. error | info | Glover et al. 2026, Table 4, this work |
| PL010 | Energy vs Haghpanah Table 6 | reference | 503.2 kWh/t | 499.6 kWh/t | 0.716 % rel. error | 1.2 % | Haghpanah et al. 2013, Table 6, P_L = 0.1 bar |
| PL010 | CO2 purity vs Glover | reference | 90.65 % | 91.8 % | 1.15 pp abs. error | 1.35 pp | Glover et al. 2026, Table 4, this work |
| PL010 | CO2 recovery vs Glover | reference | 91.57 % | 91.1 % | 0.47 pp abs. error | 0.6 pp | Glover et al. 2026, Table 4, this work |
| PL010 | Productivity vs Glover | reference | 0.7213 mol/(m3 s) | 0.721 mol/(m3 s) | 0.0354 % rel. error | 0.5 % | Glover et al. 2026, Table 4, this work |
| PL010 | Productivity vs Glover Ref. | reference | 0.7213 mol/(m3 s) | 0.694 mol/(m3 s) | 3.93 % rel. error | 4.6 % | Glover et al. 2026, Table 4, Ref. |
| PL010 | Energy vs Glover | reference | 503.2 kWh/t | 486.8 kWh/t | 3.37 % rel. error | 4 % | Glover et al. 2026, Table 4, this work |
| Sec61 | CO2 purity | reference | 79.68 % | 88.74 % | 9.06 pp abs. error | info | Haghpanah et al. 2013, Sect. 6.1, p. 4257 |
| Sec61 | CO2 recovery | reference | 37.42 % | 35.83 % | 1.59 pp abs. error | 1.9 pp | Haghpanah et al. 2013, Sect. 6.1, p. 4257 |

??? note "Notes on these comparisons"

    - dominant Floquet multiplier of the cycle map at the Figure 5 design, Glover's model on 10 cells; least-squares slope of ln(state change) over the geometric window, N2-saturated start
    - dominant Floquet multiplier of the cycle map at the Figure 5 design, Glover's model on 10 cells; median of the CSS check's fitted convergence rate over the window; spread 0.001225
    - successive substitution to a 1e-5 relative state change, bed saturated with N2; engine state-change norm, not Glover's max norm
    - successive substitution to a 1e-5 relative state change, bed saturated with feed gas; engine state-change norm, not Glover's max norm
    - minimum-energy point of the fixed-P_L Pareto front; paper model on 30 cells
    - Glover's model (constant LDF 26.25 1/s, 10 cells) on Haghpanah's design; printed as the fraction 0.906
    - Glover's model on Haghpanah's design; printed as the fraction 0.894
    - per column volume; per column volume
    - per column volume; not printed in Haghpanah et al. 2013; per column volume; Glover's reference value, not printed by Haghpanah
    - Glover's model on Haghpanah's design
    - per volume of adsorbent; at 149 kWh/t, P_L = 0.02 bar; per adsorbent volume (1 - eps) A L
    - Glover's model on Haghpanah's design; printed as the fraction 0.906
    - Glover's model on Haghpanah's design; printed as the fraction 0.907
    - Glover's model on Haghpanah's design; printed as the fraction 0.903
    - Glover's model on Haghpanah's design; printed as the fraction 0.909
    - Glover's model on Haghpanah's design; printed as the fraction 0.910
    - Glover's model on Haghpanah's design; printed as the fraction 0.914
    - Glover's model on Haghpanah's design; printed as the fraction 0.918
    - Glover's model on Haghpanah's design; printed as the fraction 0.911
    - paper model, WENO on 1000 cells, 500 cycles; design of Table 4

Input decks, assumptions and known discrepancies: [`core/test/papers/haghpanah_2013/README.md`](https://github.com/samaffleck/skarstrom/blob/main/core/test/papers/haghpanah_2013/README.md).

## Nguyen et al. (2022) {#nguyen-2022}

T.T.T. Nguyen, J.-B. Lin, G.K.H. Shimizu, A. Rajendran, Separation of CO2 and N2 on a hydrophobic metal organic framework CALF-20, Chemical Engineering Journal 442 (2022) 136263. [doi:10.1016/j.cej.2022.136263](https://doi.org/10.1016/j.cej.2022.136263)

*System:* CO2/He and CO2/N2 breakthrough on CALF-20 (six runs of Figs. 4 and 5). *Compared against:* measured markers re-digitised from Figs. 4 and 5 (outlet y_CO2, y_N2, thermocouple); the paper's model curves (code to code).

**Findings.**

- The curves used until 2026-11-12 as the experiment were the paper's model lines, digitised; the measured markers were re-digitised from the 300 dpi figures, with each run's flow, temperature and pressure from the author's thesis (Table A2.1). Fig. 5 is a CO2/N2 feed (the decks fed CO2 in helium before).
- With the paper's parameters and model forms (nothing fitted), the CO2 fronts are +6.5, +3.2, +10.9 % (Fig. 4: 100, 70, 50 % CO2) and -7.2, +5.6, -9.7 % (Fig. 5: 30, 15, 5 %) from the measurements; the paper's own model is +7.6, +5.0, +15.4, -9.6, +10.0, -5.8 %. Code to code the fronts agree with the paper's model within 0.1-4.1 %.
- Fig. 5 thermocouple RMSE 1.0-2.7 K and outlet N2 RMSE 0.04-0.28; the Fig. 4 temperature peaks are 6-8 K above the measurement and 5-6 K above the paper's model (open item). Calibrating the thermal values on Fig. 4, as the paper did, makes the Fig. 5 temperatures worse, so no calibration is shipped. The particle diameter (1.5 mm, pellets 1-3 mm) and the wall-end condition are assumptions.

| Case | Quantity | Against | Simulated | Reference | Error | Tolerance | Source |
|---|---|---|---|---|---|---|---|
| Fig4_100pct | outlet y_CO2 | experiment | curve | digitised data | 29.4 % NRMSE | 33 % | Fig. 4 (measured markers, re-digitised) |
| Fig4_100pct | CO2 front time (y_CO2 = y_feed/2) | experiment | 1236 s | 1160 s | 6.52 % rel. error | 7.5 % | Fig. 4 (measured markers, re-digitised) |
| Fig4_100pct | thermocouple T (1.82 cm from the outlet) | experiment | curve | digitised data | 7.97 K RMSE | 9 K | Fig. 4 (measured markers, re-digitised) |
| Fig4_100pct | outlet y_CO2 vs paper model | paper model | curve | digitised data | 15.3 % NRMSE | 17 % | Fig. 4 (paper's model curve, digitised) |
| Fig4_100pct | CO2 front time (y_CO2 = y_feed/2) vs paper model | paper model | 1236 s | 1258 s | 1.79 % rel. error | 2 % | Fig. 4 (paper's model curve, digitised) |
| Fig4_50pct | outlet y_CO2 | experiment | curve | digitised data | 18.8 % NRMSE | 21 % | Fig. 4 (measured markers, re-digitised) |
| Fig4_50pct | CO2 front time (y_CO2 = y_feed/2) | experiment | 2346 s | 2116 s | 10.9 % rel. error | 12 % | Fig. 4 (measured markers, re-digitised) |
| Fig4_50pct | thermocouple T (1.82 cm from the outlet) | experiment | curve | digitised data | 6.28 K RMSE | 7 K | Fig. 4 (measured markers, re-digitised) |
| Fig4_50pct | outlet y_CO2 vs paper model | paper model | curve | digitised data | 10.7 % NRMSE | 12 % | Fig. 4 (paper's model curve, digitised) |
| Fig4_50pct | CO2 front time (y_CO2 = y_feed/2) vs paper model | paper model | 2346 s | 2441 s | 3.89 % rel. error | 4.5 % | Fig. 4 (paper's model curve, digitised) |
| Fig4_70pct | outlet y_CO2 | experiment | curve | digitised data | 11.8 % NRMSE | 13 % | Fig. 4 (measured markers, re-digitised) |
| Fig4_70pct | CO2 front time (y_CO2 = y_feed/2) | experiment | 1718 s | 1664 s | 3.25 % rel. error | 4 % | Fig. 4 (measured markers, re-digitised) |
| Fig4_70pct | thermocouple T (1.82 cm from the outlet) | experiment | curve | digitised data | 5.23 K RMSE | 6 K | Fig. 4 (measured markers, re-digitised) |
| Fig4_70pct | outlet y_CO2 vs paper model | paper model | curve | digitised data | 8.88 % NRMSE | 10 % | Fig. 4 (paper's model curve, digitised) |
| Fig4_70pct | CO2 front time (y_CO2 = y_feed/2) vs paper model | paper model | 1718 s | 1747 s | 1.65 % rel. error | 2 % | Fig. 4 (paper's model curve, digitised) |
| Fig5_15pct | outlet y_CO2 | experiment | curve | digitised data | 13.4 % NRMSE | 15 % | Fig. 5 (measured markers, re-digitised) |
| Fig5_15pct | CO2 front time (y_CO2 = y_feed/2) | experiment | 3229 s | 3059 s | 5.56 % rel. error | 6.5 % | Fig. 5 (measured markers, re-digitised) |
| Fig5_15pct | thermocouple T (1.82 cm from the outlet) | experiment | curve | digitised data | 1.03 K RMSE | 1.25 K | Fig. 5 (measured markers, re-digitised) |
| Fig5_15pct | outlet y_N2 | experiment | curve | digitised data | 0.139 RMSE | 0.16 | Fig. 5 (measured markers, re-digitised) |
| Fig5_15pct | outlet y_CO2 vs paper model | paper model | curve | digitised data | 8.49 % NRMSE | 10 % | Fig. 5 (paper's model curve, digitised) |
| Fig5_15pct | CO2 front time (y_CO2 = y_feed/2) vs paper model | paper model | 3229 s | 3364 s | 4.01 % rel. error | 4.5 % | Fig. 5 (paper's model curve, digitised) |
| Fig5_30pct | outlet y_CO2 | experiment | curve | digitised data | 15.2 % NRMSE | 17 % | Fig. 5 (measured markers, re-digitised) |
| Fig5_30pct | CO2 front time (y_CO2 = y_feed/2) | experiment | 1706 s | 1839 s | 7.25 % rel. error | 8 % | Fig. 5 (measured markers, re-digitised) |
| Fig5_30pct | thermocouple T (1.82 cm from the outlet) | experiment | curve | digitised data | 2.66 K RMSE | 3 K | Fig. 5 (measured markers, re-digitised) |
| Fig5_30pct | outlet y_N2 | experiment | curve | digitised data | 0.282 RMSE | 0.32 | Fig. 5 (measured markers, re-digitised) |
| Fig5_30pct | outlet y_CO2 vs paper model | paper model | curve | digitised data | 9.94 % NRMSE | 11 % | Fig. 5 (paper's model curve, digitised) |
| Fig5_30pct | CO2 front time (y_CO2 = y_feed/2) vs paper model | paper model | 1706 s | 1708 s | 0.146 % rel. error | 0.5 % | Fig. 5 (paper's model curve, digitised) |
| Fig5_5pct | outlet y_CO2 | experiment | curve | digitised data | 9.49 % NRMSE | 11 % | Fig. 5 (measured markers, re-digitised) |
| Fig5_5pct | CO2 front time (y_CO2 = y_feed/2) | experiment | 7091 s | 7849 s | 9.66 % rel. error | 11 % | Fig. 5 (measured markers, re-digitised) |
| Fig5_5pct | thermocouple T (1.82 cm from the outlet) | experiment | curve | digitised data | 1.72 K RMSE | 2 K | Fig. 5 (measured markers, re-digitised) |
| Fig5_5pct | outlet y_N2 | experiment | curve | digitised data | 0.0355 RMSE | 0.04 | Fig. 5 (measured markers, re-digitised) |
| Fig5_5pct | outlet y_CO2 vs paper model | paper model | curve | digitised data | 13 % NRMSE | 15 % | Fig. 5 (paper's model curve, digitised) |
| Fig5_5pct | CO2 front time (y_CO2 = y_feed/2) vs paper model | paper model | 7091 s | 7397 s | 4.14 % rel. error | 5 % | Fig. 5 (paper's model curve, digitised) |

??? note "Notes on these comparisons"

    - RMSE over 13 digitised points; normalised by 1 -
    - first time the curve reaches 0.5 (reference time interpolated from the digitised points)
    - RMSE over 25 digitised points
    - RMSE over 15 digitised points; normalised by 1 -
    - RMSE over 64 digitised points; normalised by 0.5 -
    - first time the curve reaches 0.25 (reference time interpolated from the digitised points)
    - RMSE over 31 digitised points
    - RMSE over 24 digitised points; normalised by 0.5 -
    - RMSE over 70 digitised points; normalised by 0.7 -
    - first time the curve reaches 0.35 (reference time interpolated from the digitised points)
    - RMSE over 34 digitised points
    - RMSE over 22 digitised points; normalised by 0.7 -
    - RMSE over 45 digitised points; normalised by 0.15 -
    - first time the curve reaches 0.075 (reference time interpolated from the digitised points)
    - RMSE over 38 digitised points
    - RMSE over 28 digitised points
    - RMSE over 31 digitised points; normalised by 0.15 -
    - RMSE over 26 digitised points; normalised by 0.3 -
    - first time the curve reaches 0.15 (reference time interpolated from the digitised points)
    - RMSE over 29 digitised points
    - RMSE over 14 digitised points
    - RMSE over 32 digitised points; normalised by 0.3 -
    - RMSE over 81 digitised points; normalised by 0.05 -
    - first time the curve reaches 0.025 (reference time interpolated from the digitised points)
    - RMSE over 69 digitised points
    - RMSE over 48 digitised points
    - RMSE over 16 digitised points; normalised by 0.05 -

Input decks, assumptions and known discrepancies: [`core/test/papers/nguyen_2022/README.md`](https://github.com/samaffleck/skarstrom/blob/main/core/test/papers/nguyen_2022/README.md).

## Nguyen et al. (2023) {#nguyen-2023}

T.T.T. Nguyen, G.K.H. Shimizu, A. Rajendran, CO2/N2 separation by vacuum swing adsorption using a metal-organic framework, CALF-20: Multi-objective optimization and experimental validation, Chemical Engineering Journal 452 (2023) 139550. [doi:10.1016/j.cej.2022.139550](https://doi.org/10.1016/j.cej.2022.139550)

*System:* CO2/N2 VSA on CALF-20 (4-step and LPP). *Compared against:* experimental and model purity/recovery/productivity (Table 2).

**Findings.**

- Four-step cycles (I-V): CO2 purity within 1.7-4.9 pp and recovery within 1.2-2.7 pp of the experiments; productivity within 5-9 %.
- LPP cycles (VI-X): purity within 0.2-2.7 pp of the experiments, but recovery is 7-10 pp high and exceeds 100 % for IX and X, at the cyclic steady state of the per-component mass balance (VIII-X need about 300-500 cycles to reach it). The single-bed emulation couples the light-product pressurisation to the feed step by composition, not by moles, so CO2 enters through the LPP stream that is not drawn from the raffinate; this is a limitation of the deck/connection model, recorded rather than hidden.
- The decks now use the competitive DSL form (as the paper's ref. [10]); the paper's column length (0.267 m) is inconsistent with its adsorbent mass and particle density, and the deck keeps the mass (L = 0.304 m).

| Case | Quantity | Against | Simulated | Reference | Error | Tolerance | Source |
|---|---|---|---|---|---|---|---|
| ExperimentI | CO2 purity vs experiment | experiment | 95.71 % | 97.7 % | 1.99 pp abs. error | 2.4 pp | Table 2, Exp. I, Experimental Pu |
| ExperimentI | CO2 recovery vs experiment | experiment | 69.4 % | 66.9 % | 2.5 pp abs. error | 2.8 pp | Table 2, Exp. I, Experimental Re |
| ExperimentI | CO2 productivity vs experiment | experiment | 0.608 mol CO2/(m3 s) | 0.65 mol CO2/(m3 s) | 6.46 % rel. error | 7.6 % | Table 2, Exp. I, Experimental Pr |
| ExperimentI | CO2 purity vs paper model | paper model | 95.71 % | 96.7 % | 0.99 pp abs. error | 1.4 pp | Table 2, Exp. I, Model Pu |
| ExperimentI | CO2 recovery vs paper model | paper model | 69.4 % | 70.1 % | 0.695 pp abs. error | 1.1 pp | Table 2, Exp. I, Model Re |
| ExperimentI | CO2 productivity vs paper model | paper model | 0.608 mol CO2/(m3 s) | 0.63 mol CO2/(m3 s) | 3.49 % rel. error | 4.6 % | Table 2, Exp. I, Model Pr |
| ExperimentII | CO2 purity vs experiment | experiment | 89.36 % | 91.2 % | 1.84 pp abs. error | 2.2 pp | Table 2, Exp. II, Experimental Pu |
| ExperimentII | CO2 recovery vs experiment | experiment | 73 % | 71.7 % | 1.3 pp abs. error | 1.6 pp | Table 2, Exp. II, Experimental Re |
| ExperimentII | CO2 productivity vs experiment | experiment | 0.4715 mol CO2/(m3 s) | 0.52 mol CO2/(m3 s) | 9.33 % rel. error | 11 % | Table 2, Exp. II, Experimental Pr |
| ExperimentII | CO2 purity vs paper model | paper model | 89.36 % | 90.3 % | 0.94 pp abs. error | 1.3 pp | Table 2, Exp. II, Model Pu |
| ExperimentII | CO2 recovery vs paper model | paper model | 73 % | 72.3 % | 0.704 pp abs. error | 0.94 pp | Table 2, Exp. II, Model Re |
| ExperimentII | CO2 productivity vs paper model | paper model | 0.4715 mol CO2/(m3 s) | 0.54 mol CO2/(m3 s) | 12.7 % rel. error | 15 % | Table 2, Exp. II, Model Pr |
| ExperimentIII | CO2 purity vs experiment | experiment | 83.92 % | 85.6 % | 1.68 pp abs. error | 2.1 pp | Table 2, Exp. III, Experimental Pu |
| ExperimentIII | CO2 recovery vs experiment | experiment | 74.36 % | 71.8 % | 2.56 pp abs. error | 2.9 pp | Table 2, Exp. III, Experimental Re |
| ExperimentIII | CO2 productivity vs experiment | experiment | 0.4396 mol CO2/(m3 s) | 0.48 mol CO2/(m3 s) | 8.42 % rel. error | 9.8 % | Table 2, Exp. III, Experimental Pr |
| ExperimentIII | CO2 purity vs paper model | paper model | 83.92 % | 85.7 % | 1.78 pp abs. error | 2.2 pp | Table 2, Exp. III, Model Pu |
| ExperimentIII | CO2 recovery vs paper model | paper model | 74.36 % | 73.4 % | 0.956 pp abs. error | 1.2 pp | Table 2, Exp. III, Model Re |
| ExperimentIII | CO2 productivity vs paper model | paper model | 0.4396 mol CO2/(m3 s) | 0.51 mol CO2/(m3 s) | 13.8 % rel. error | 16 % | Table 2, Exp. III, Model Pr |
| ExperimentIV | CO2 purity vs experiment | experiment | 79.15 % | 82.3 % | 3.15 pp abs. error | 3.7 pp | Table 2, Exp. IV, Experimental Pu |
| ExperimentIV | CO2 recovery vs experiment | experiment | 75.22 % | 72.5 % | 2.72 pp abs. error | 3.1 pp | Table 2, Exp. IV, Experimental Re |
| ExperimentIV | CO2 productivity vs experiment | experiment | 0.428 mol CO2/(m3 s) | 0.46 mol CO2/(m3 s) | 6.96 % rel. error | 8.1 % | Table 2, Exp. IV, Experimental Pr |
| ExperimentIV | CO2 purity vs paper model | paper model | 79.15 % | 81.3 % | 2.15 pp abs. error | 2.6 pp | Table 2, Exp. IV, Model Pu |
| ExperimentIV | CO2 recovery vs paper model | paper model | 75.22 % | 74.3 % | 0.921 pp abs. error | 1.2 pp | Table 2, Exp. IV, Model Re |
| ExperimentIV | CO2 productivity vs paper model | paper model | 0.428 mol CO2/(m3 s) | 0.5 mol CO2/(m3 s) | 14.4 % rel. error | 17 % | Table 2, Exp. IV, Model Pr |
| ExperimentIX | CO2 purity vs experiment | experiment | 84.4 % | 85.67 % | 1.27 pp abs. error | 1.6 pp | Table 2, Exp. IX, Experimental Pu |
| ExperimentIX | CO2 recovery vs experiment | experiment | 98.89 % | 90.5 % | 8.39 pp abs. error | 12 pp | Table 2, Exp. IX, Experimental Re |
| ExperimentIX | CO2 productivity vs experiment | experiment | 0.1431 mol CO2/(m3 s) | 0.13 mol CO2/(m3 s) | 10 % rel. error | 14 % | Table 2, Exp. IX, Experimental Pr |
| ExperimentIX | CO2 purity vs paper model | paper model | 84.4 % | 88.3 % | 3.9 pp abs. error | 4.2 pp | Table 2, Exp. IX, Model Pu |
| ExperimentIX | CO2 recovery vs paper model | paper model | 98.89 % | 91.6 % | 7.29 pp abs. error | 11 pp | Table 2, Exp. IX, Model Re |
| ExperimentIX | CO2 productivity vs paper model | paper model | 0.1431 mol CO2/(m3 s) | 0.14 mol CO2/(m3 s) | 2.18 % rel. error | 4.7 % | Table 2, Exp. IX, Model Pr |
| ExperimentV | CO2 purity vs experiment | experiment | 73.55 % | 78.4 % | 4.85 pp abs. error | 5.7 pp | Table 2, Exp. V, Experimental Pu |
| ExperimentV | CO2 recovery vs experiment | experiment | 76.83 % | 74.2 % | 2.63 pp abs. error | 3 pp | Table 2, Exp. V, Experimental Re |
| ExperimentV | CO2 productivity vs experiment | experiment | 0.3805 mol CO2/(m3 s) | 0.4 mol CO2/(m3 s) | 4.88 % rel. error | 6 % | Table 2, Exp. V, Experimental Pr |
| ExperimentV | CO2 purity vs paper model | paper model | 73.55 % | 75.6 % | 2.05 pp abs. error | 2.4 pp | Table 2, Exp. V, Model Pu |
| ExperimentV | CO2 recovery vs paper model | paper model | 76.83 % | 75.4 % | 1.43 pp abs. error | 1.7 pp | Table 2, Exp. V, Model Re |
| ExperimentV | CO2 productivity vs paper model | paper model | 0.3805 mol CO2/(m3 s) | 0.44 mol CO2/(m3 s) | 13.5 % rel. error | 16 % | Table 2, Exp. V, Model Pr |
| ExperimentVI | CO2 purity vs experiment | experiment | 96.21 % | 97.2 % | 0.99 pp abs. error | 1.3 pp | Table 2, Exp. VI, Experimental Pu |
| ExperimentVI | CO2 recovery vs experiment | experiment | 89.83 % | 83.1 % | 6.73 pp abs. error | 7.8 pp | Table 2, Exp. VI, Experimental Re |
| ExperimentVI | CO2 productivity vs experiment | experiment | 0.2801 mol CO2/(m3 s) | 0.3 mol CO2/(m3 s) | 6.62 % rel. error | 7.6 % | Table 2, Exp. VI, Experimental Pr |
| ExperimentVI | CO2 purity vs paper model | paper model | 96.21 % | 96.7 % | 0.49 pp abs. error | 0.74 pp | Table 2, Exp. VI, Model Pu |
| ExperimentVI | CO2 recovery vs paper model | paper model | 89.83 % | 84.8 % | 5.03 pp abs. error | 5.9 pp | Table 2, Exp. VI, Model Re |
| ExperimentVI | CO2 productivity vs paper model | paper model | 0.2801 mol CO2/(m3 s) | 0.28 mol CO2/(m3 s) | 0.0501 % rel. error | 0.26 % | Table 2, Exp. VI, Model Pr |
| ExperimentVII | CO2 purity vs experiment | experiment | 93.46 % | 94.3 % | 0.843 pp abs. error | 1.1 pp | Table 2, Exp. VII.a, Experimental Pu |
| ExperimentVII | CO2 recovery vs experiment | experiment | 95.24 % | 88.4 % | 6.84 pp abs. error | 9 pp | Table 2, Exp. VII.a, Experimental Re |
| ExperimentVII | CO2 productivity vs experiment | experiment | 0.158 mol CO2/(m3 s) | 0.14 mol CO2/(m3 s) | 12.8 % rel. error | 16 % | Table 2, Exp. VII.a, Experimental Pr |
| ExperimentVII | CO2 purity vs paper model | paper model | 93.46 % | 95.5 % | 2.04 pp abs. error | 2.3 pp | Table 2, Exp. VII, Model Pu |
| ExperimentVII | CO2 recovery vs paper model | paper model | 95.24 % | 88.1 % | 7.14 pp abs. error | 9.3 pp | Table 2, Exp. VII, Model Re |
| ExperimentVII | CO2 productivity vs paper model | paper model | 0.158 mol CO2/(m3 s) | 0.14 mol CO2/(m3 s) | 12.8 % rel. error | 16 % | Table 2, Exp. VII, Model Pr |
| ExperimentVII | CO2 purity vs repeat experiment VII.b | experiment | 93.46 % | 93.8 % | 0.343 pp abs. error | info | Table 2, Exp. VII.b, Experimental Pu |
| ExperimentVII | CO2 recovery vs repeat experiment VII.b | experiment | 95.24 % | 87.5 % | 7.74 pp abs. error | info | Table 2, Exp. VII.b, Experimental Re |
| ExperimentVIII | CO2 purity vs experiment | experiment | 89.9 % | 92.7 % | 2.8 pp abs. error | 3.1 pp | Table 2, Exp. VIII, Experimental Pu |
| ExperimentVIII | CO2 recovery vs experiment | experiment | 96.14 % | 87.8 % | 8.34 pp abs. error | 11 pp | Table 2, Exp. VIII, Experimental Re |
| ExperimentVIII | CO2 productivity vs experiment | experiment | 0.1551 mol CO2/(m3 s) | 0.17 mol CO2/(m3 s) | 8.74 % rel. error | 11 % | Table 2, Exp. VIII, Experimental Pr |
| ExperimentVIII | CO2 purity vs paper model | paper model | 89.9 % | 93.6 % | 3.7 pp abs. error | 4.1 pp | Table 2, Exp. VIII, Model Pu |
| ExperimentVIII | CO2 recovery vs paper model | paper model | 96.14 % | 89.9 % | 6.24 pp abs. error | 8.3 pp | Table 2, Exp. VIII, Model Re |
| ExperimentVIII | CO2 productivity vs paper model | paper model | 0.1551 mol CO2/(m3 s) | 0.16 mol CO2/(m3 s) | 3.04 % rel. error | 4 % | Table 2, Exp. VIII, Model Pr |
| ExperimentX | CO2 purity vs experiment | experiment | 80.56 % | 80.6 % | 0.0355 pp abs. error | 0.13 pp | Table 2, Exp. X, Experimental Pu |
| ExperimentX | CO2 recovery vs experiment | experiment | 99.43 % | 92.4 % | 7.03 pp abs. error | 9.7 pp | Table 2, Exp. X, Experimental Re |
| ExperimentX | CO2 productivity vs experiment | experiment | 0.1507 mol CO2/(m3 s) | 0.15 mol CO2/(m3 s) | 0.458 % rel. error | 2.9 % | Table 2, Exp. X, Experimental Pr |
| ExperimentX | CO2 purity vs paper model | paper model | 80.56 % | 80.6 % | 0.0355 pp abs. error | 0.13 pp | Table 2, Exp. X, Model Pu |
| ExperimentX | CO2 recovery vs paper model | paper model | 99.43 % | 93.6 % | 5.83 pp abs. error | 8.3 pp | Table 2, Exp. X, Model Re |
| ExperimentX | CO2 productivity vs paper model | paper model | 0.1507 mol CO2/(m3 s) | 0.13 mol CO2/(m3 s) | 15.9 % rel. error | 21 % | Table 2, Exp. X, Model Pr |

??? note "Notes on these comparisons"

    - measured at CSS
    - per m3 of adsorbent per second (Sect. 2.3)
    - paper's detailed model
    - measured at CSS (LPP cycle)
    - repeat VII.b: 93.8 %
    - repeat VII.b: 87.5 %
    - repeat VII.b: 0.13
    - repeat experiment; repeat run of VII

Input decks, assumptions and known discrepancies: [`core/test/papers/nguyen_2023/README.md`](https://github.com/samaffleck/skarstrom/blob/main/core/test/papers/nguyen_2023/README.md).

## Nguyen et al. (2025) {#nguyen-2024}

T.T.T. Nguyen, G.K.H. Shimizu, A. Rajendran, Experimental demonstration of humid post-combustion CO2 capture by vacuum swing adsorption using CALF-20, Adsorption 31 (2025) 45. [doi:10.1007/s10450-025-00600-z](https://doi.org/10.1007/s10450-025-00600-z)

*System:* CO2/N2 VSA with light-product pressurisation on CALF-20 (dry cases of a humid-feed study). *Compared against:* experimental purity/recovery/productivity of the dry cases (Table 1).

**Findings.**

- Dry case B0 (the operating point of Nguyen 2023 Exp. VII): CO2 purity within 0.9 pp and productivity within 2 % of Table 1, inside the printed uncertainty; the productivity agreement is partly compensation, because the engine's recovery is 11.6 pp high. Dry case A0: purity 3.4 pp low (recorded).
- Recovery is about 11 pp high in both cases, the light-product over-recovery of the Nguyen 2023 LPP decks; none of the printed parameters closes it (README attribution: the SI's fitted blowdown and evacuation decay rates lower it by up to 5 pp and raise the purity, the measured 0.267 m column raises it by about 2 pp). A0's printed productivity is inconsistent with its own cycle time and is recorded, not asserted.
- The decks are the Nguyen2023_Exp_VII deck at 0.96 bar and 22 C. The main text was not read here (the SI was); the humid cases are not reproduced: the paper prints no water isotherm, and the engine has no competitive CO2/H2O kind for a physisorbent (#272).

| Case | Quantity | Against | Simulated | Reference | Error | Tolerance | Source |
|---|---|---|---|---|---|---|---|
| CaseA0 | CO2 purity vs experiment | experiment | 94.89 % | 98.3 % | 3.41 pp abs. error | info | Table 1, p. 6, Case 1, 0 % RH (A0), purity |
| CaseA0 | CO2 recovery vs experiment | experiment | 81.98 % | 71.2 % | 10.8 pp abs. error | info | Table 1, p. 6, Case 1, 0 % RH (A0), recovery |
| CaseA0 | CO2 productivity vs experiment | experiment | 0.4671 mol CO2/(m3 s) | 0.577 mol CO2/(m3 s) | 19.1 % rel. error | info | Table 1, p. 6, Case 1, 0 % RH (A0), productivity |
| CaseB0 | CO2 purity vs experiment | experiment | 93.33 % | 92.4 % | 0.931 pp abs. error | 1 pp | Table 1, p. 6, Case 2, 0 % RH (B0), purity |
| CaseB0 | CO2 productivity vs experiment | experiment | 0.15 mol CO2/(m3 s) | 0.147 mol CO2/(m3 s) | 2.01 % rel. error | 2.5 % | Table 1, p. 6, Case 2, 0 % RH (B0), productivity |
| CaseB0 | CO2 recovery vs experiment | experiment | 94.87 % | 83.3 % | 11.6 pp abs. error | info | Table 1, p. 6, Case 2, 0 % RH (B0), recovery |

??? note "Notes on these comparisons"

    - measured at CSS, dry basis (Eq. 2); printed uncertainty +/- 0.88; recorded: 3.4 pp below, outside the printed uncertainty (README)
    - measured at CSS (Eq. 3); printed uncertainty +/- 1.46; recorded: LPP over-recovery (README)
    - per m3 of adsorbent per second (Eq. 4); printed uncertainty +/- 0.014; recorded: the printed value is inconsistent with the cycle time (README)
    - measured at CSS, dry basis (Eq. 2); printed uncertainty +/- 1.31
    - per m3 of adsorbent per second (Eq. 4); printed uncertainty +/- 0.005
    - measured at CSS (Eq. 3); printed uncertainty +/- 1.97; recorded: LPP over-recovery (README)

Input decks, assumptions and known discrepancies: [`core/test/papers/nguyen_2024/README.md`](https://github.com/samaffleck/skarstrom/blob/main/core/test/papers/nguyen_2024/README.md).

## Estupiñan Perez et al. (2019) {#perez-2019}

L. Estupiñan Perez, P. Sarkar, A. Rajendran, Experimental validation of multi-objective optimization techniques for design of vacuum swing adsorption processes, Separation and Purification Technology 224 (2019) 553-563. [doi:10.1016/j.seppur.2019.05.039](https://doi.org/10.1016/j.seppur.2019.05.039)

*System:* CO2/N2 VSA on zeolite 13X (4-step and LPP). *Compared against:* experimental and model purity/recovery (Table 3).

**Findings.**

- The decks use the paper's model forms (WP10 audit of the paper and SI, 2026-10-31): the two-term exponential evacuation pressure of Sect. 4.3 (0.019 and 0.101 1/s; the split of the terms is not printed and an equal split is used, chosen on A), wall ends at ambient (SI Table S2), the Blake-Kozeny pressure drop (Ergun) and the capacity-factor LDF of SI Eq. S4. Four-step cycles: CO2 purity within 0.2-0.6 pp and recovery within 0.7-1.1 pp of experiments A and B (recovery was 8-10 pp high); C1 is 4.5 pp below its experiment and 12 pp below the repeat C2, which differ from each other by 7 pp. LPP cycles: purity 1.5-3.3 pp low and recovery 1.9-3.3 pp high.
- Attribution (README): the four-step recovery excess was the evacuation pressure history. The paper's profile keeps the column above P_low for most of the 57 s step where the deck's single exponential reached it; the capacity-factor kinetics, the wall ends and the pressure-drop constant are each worth under 0.5 pp. The LPP cycles' 167 s evacuation reaches P_low whatever the profile; their remaining recovery excess is the composition-coupled LPP connection.
- Open: the split of the two evacuation terms (Fig. 6c would pin it); the experiment ended steps on pressure triggers whose durations are not printed; the one-second hold between steps and about 5 % dead volume are not modelled.

| Case | Quantity | Against | Simulated | Reference | Error | Tolerance | Source |
|---|---|---|---|---|---|---|---|
| ConfigA | CO2 purity vs experiment | experiment | 99.55 % | 99.7 % | 0.149 pp abs. error | 0.3 pp | Table 3, 4-step, Experiment A, Experimental Purity |
| ConfigA | CO2 recovery vs experiment | experiment | 60.14 % | 58.7 % | 1.44 pp abs. error | 1.5 pp | Table 3, 4-step, Experiment A, Experimental Recovery |
| ConfigA | CO2 purity vs paper model | paper model | 99.55 % | 99.3 % | 0.251 pp abs. error | 0.3 pp | Table 3, 4-step, Experiment A, Model Purity |
| ConfigA | CO2 recovery vs paper model | paper model | 60.14 % | 60 % | 0.144 pp abs. error | 0.43 pp | Table 3, 4-step, Experiment A, Model Recovery |
| ConfigB | CO2 purity vs experiment | experiment | 96.96 % | 96.1 % | 0.861 pp abs. error | 0.9 pp | Table 3, 4-step, Experiment B, Experimental Purity |
| ConfigB | CO2 recovery vs experiment | experiment | 62.03 % | 61.8 % | 0.234 pp abs. error | 0.3 pp | Table 3, 4-step, Experiment B, Experimental Recovery |
| ConfigB | CO2 purity vs paper model | paper model | 96.96 % | 95.2 % | 1.76 pp abs. error | 1.8 pp | Table 3, 4-step, Experiment B, Model Purity |
| ConfigB | CO2 recovery vs paper model | paper model | 62.03 % | 64.1 % | 2.07 pp abs. error | 4.3 pp | Table 3, 4-step, Experiment B, Model Recovery |
| ConfigC1 | CO2 purity vs experiment | experiment | 88.43 % | 88.1 % | 0.325 pp abs. error | 0.4 pp | Table 3, 4-step, Experiment C1, Experimental Purity |
| ConfigC1 | CO2 recovery vs experiment | experiment | 62.39 % | 64.6 % | 2.21 pp abs. error | 2.9 pp | Table 3, 4-step, Experiment C1, Experimental Recovery |
| ConfigC1 | CO2 purity vs paper model | paper model | 88.43 % | 83.4 % | 5.03 pp abs. error | 5.1 pp | Table 3, 4-step, Experiment C1/C2, Model Purity |
| ConfigC1 | CO2 recovery vs paper model | paper model | 62.39 % | 68.2 % | 5.81 pp abs. error | 9.7 pp | Table 3, 4-step, Experiment C1/C2, Model Recovery |
| ConfigC1 | CO2 purity vs repeat experiment C2 | experiment | 88.43 % | 86.9 % | 1.53 pp abs. error | info | Table 3, 4-step, Experiment C2, Experimental Purity |
| ConfigC1 | CO2 recovery vs repeat experiment C2 | experiment | 62.39 % | 71.9 % | 9.51 pp abs. error | info | Table 3, 4-step, Experiment C2, Experimental Recovery |
| ConfigD1 | CO2 purity vs experiment | experiment | 94.08 % | 97.1 % | 3.02 pp abs. error | 3.9 pp | Table 3, LPP, Experiment D1, Experimental Purity |
| ConfigD1 | CO2 recovery vs experiment | experiment | 89.48 % | 86.1 % | 3.38 pp abs. error | 4.5 pp | Table 3, LPP, Experiment D1, Experimental Recovery |
| ConfigD1 | CO2 purity vs paper model | paper model | 94.08 % | 96.1 % | 2.02 pp abs. error | 2.7 pp | Table 3, LPP, Experiment D1/D2, Model Purity |
| ConfigD1 | CO2 recovery vs paper model | paper model | 89.48 % | 86.8 % | 2.68 pp abs. error | 3.7 pp | Table 3, LPP, Experiment D1/D2, Model Recovery |
| ConfigD1 | CO2 purity vs repeat experiment D2 | experiment | 94.08 % | 94.5 % | 0.416 pp abs. error | info | Table 3, LPP, Experiment D2, Experimental Purity |
| ConfigD1 | CO2 recovery vs repeat experiment D2 | experiment | 89.48 % | 83.4 % | 6.08 pp abs. error | info | Table 3, LPP, Experiment D2, Experimental Recovery |
| ConfigE | CO2 purity vs experiment | experiment | 87.25 % | 88.5 % | 1.25 pp abs. error | 2.1 pp | Table 3, LPP, Experiment E, Experimental Purity |
| ConfigE | CO2 recovery vs experiment | experiment | 93.91 % | 91.7 % | 2.21 pp abs. error | 3.2 pp | Table 3, LPP, Experiment E, Experimental Recovery |
| ConfigE | CO2 purity vs paper model | paper model | 87.25 % | 90.4 % | 3.15 pp abs. error | 4.3 pp | Table 3, LPP, Experiment E, Model Purity |
| ConfigE | CO2 recovery vs paper model | paper model | 93.91 % | 91.1 % | 2.81 pp abs. error | 4 pp | Table 3, LPP, Experiment E, Model Recovery |
| ConfigF | CO2 purity vs experiment | experiment | 76.78 % | 77.7 % | 0.916 pp abs. error | 2.6 pp | Table 3, LPP, Experiment F, Experimental Purity |
| ConfigF | CO2 recovery vs experiment | experiment | 97.2 % | 95.6 % | 1.6 pp abs. error | 2.6 pp | Table 3, LPP, Experiment F, Experimental Recovery |
| ConfigF | CO2 purity vs paper model | paper model | 76.78 % | 79.9 % | 3.12 pp abs. error | 4.4 pp | Table 3, LPP, Experiment F, Model Purity |
| ConfigF | CO2 recovery vs paper model | paper model | 97.2 % | 94.6 % | 2.6 pp abs. error | 3.8 pp | Table 3, LPP, Experiment F, Model Recovery |

??? note "Notes on these comparisons"

    - measured at CSS
    - paper's detailed model
    - repeat C2: 86.9 %
    - repeat C2: 71.9 %
    - repeat of C1; repeat run
    - repeat D2: 94.5 %
    - repeat D2: 83.4 %
    - repeat of D1; repeat run

Input decks, assumptions and known discrepancies: [`core/test/papers/perez_2019/README.md`](https://github.com/samaffleck/skarstrom/blob/main/core/test/papers/perez_2019/README.md).

## Streb & Mazzotti (2021) {#streb-2021}

A. Streb, M. Mazzotti, Adsorption for efficient low carbon hydrogen production: part 1 - adsorption equilibrium and breakthrough studies for H2/CO2/CH4 on zeolite 13X, Adsorption 27 (2021) 541-558. [doi:10.1007/s10450-021-00306-y](https://doi.org/10.1007/s10450-021-00306-y)

*System:* CO2/CH4 and H2/CO2/CH4 breakthrough on zeolite 13X. *Compared against:* digitised experimental outlet compositions and thermocouples (Figs. 6, 7, 11); decks calibrated on A3.

**Findings.**

- The decks are calibrated on A3 alone: CO2 and CH4 LDF coefficients 0.0592 and 0.137 1/s (Table 4: 0.06, 0.2), wall-to-thermofluid coefficient 49.7 W/(m2 K) (paper 220) and a 0.142 L plug-flow sampling line to the MS (an inert, nearly empty layer above the bed). Each other case is a prediction except its feed flow, which is fitted to the speed of its own CO2 temperature front (as the paper fits the velocity per experiment). CO2 fronts are within 1.4 % of the experiments and CH4 fronts within 2.4 % (B1 4.7 % early); outlet NRMSE 1.3-16 % (CO2) and 4.7-31 % (CH4); thermocouple RMSE 5.5-11.4 K. Before: 0.3-11 %, 0.2-18 %, 2.4-45 %, 13-84 % and 8.7-28.5 K.
- Six flow factors exceed 2 % of the set point (A1 +5.1, A2 -6.1, A4 -9.4, A6 -3.5, B1 +5.7, B2 -11.3 %). The paper's bound is 2 % of the MFC full scale, which is not printed: with a full scale of at least 31 NL/min every case but A4 is inside it. Held to 2 % of the set point, A1, A2, A4, B1 and B2 keep front errors of 3-10 % (README).
- The engine's void region is one well-mixed cell, which smears the CH4 front, so the line is a multi-cell inert layer. The decks keep extended Sips competition; IAST, recorded on the same calibrated decks, puts the CO2 fronts 0.3-11 % late.

| Case | Quantity | Against | Simulated | Reference | Error | Tolerance | Source |
|---|---|---|---|---|---|---|---|
| ExperimentA1 | outlet y_CO2 | experiment | curve | digitised data | 16.5 % NRMSE | 19 % | Fig. 7 (digitised) |
| ExperimentA1 | outlet y_CH4 | experiment | curve | digitised data | 20.3 % NRMSE | 23 % | Fig. 7 (digitised) |
| ExperimentA1 | CO2 front time (y_CO2 = y_feed/2) | experiment | 3686 s | 3663 s | 0.639 % rel. error | 1 % | Fig. 7 (digitised) |
| ExperimentA1 | CH4 front time (y_CH4 = y_feed/2) | experiment | 464.5 s | 458.6 s | 1.3 % rel. error | 1.5 % | Fig. 7 (digitised) |
| ExperimentA1 | thermocouple T (5 positions, pooled) | experiment | curve | digitised data | 6.91 K RMSE | 8 K | Fig. 7 (digitised) |
| ExperimentA1_IAST | outlet y_CO2 vs experiment (IAST) | experiment | curve | digitised data | 32 % NRMSE | info | Fig. 7 (digitised) |
| ExperimentA1_IAST | outlet y_CH4 vs experiment (IAST) | experiment | curve | digitised data | 34.3 % NRMSE | info | Fig. 7 (digitised) |
| ExperimentA1_IAST | CO2 front time (y_CO2 = y_feed/2) vs experiment (IAST) | experiment | 3773 s | 3663 s | 3.01 % rel. error | info | Fig. 7 (digitised) |
| ExperimentA1_IAST | CH4 front time (y_CH4 = y_feed/2) vs experiment (IAST) | experiment | 457.3 s | 458.6 s | 0.277 % rel. error | info | Fig. 7 (digitised) |
| ExperimentA1_IAST | thermocouple T (5 positions, pooled) vs experiment (IAST) | experiment | curve | digitised data | 13.6 K RMSE | info | Fig. 7 (digitised) |
| ExperimentA2 | outlet y_CO2 | experiment | curve | digitised data | 10.5 % NRMSE | 12 % | Fig. 7 (digitised) |
| ExperimentA2 | outlet y_CH4 | experiment | curve | digitised data | 26.4 % NRMSE | 30 % | Fig. 7 (digitised) |
| ExperimentA2 | CO2 front time (y_CO2 = y_feed/2) | experiment | 973.2 s | 965.7 s | 0.783 % rel. error | 1 % | Fig. 7 (digitised) |
| ExperimentA2 | CH4 front time (y_CH4 = y_feed/2) | experiment | 371.4 s | 365 s | 1.76 % rel. error | 2 % | Fig. 7 (digitised) |
| ExperimentA2 | thermocouple T (5 positions, pooled) | experiment | curve | digitised data | 8.51 K RMSE | 9.5 K | Fig. 7 (digitised) |
| ExperimentA2_IAST | outlet y_CO2 vs experiment (IAST) | experiment | curve | digitised data | 44.7 % NRMSE | info | Fig. 7 (digitised) |
| ExperimentA2_IAST | outlet y_CH4 vs experiment (IAST) | experiment | curve | digitised data | 44.8 % NRMSE | info | Fig. 7 (digitised) |
| ExperimentA2_IAST | CO2 front time (y_CO2 = y_feed/2) vs experiment (IAST) | experiment | 1072 s | 965.7 s | 11 % rel. error | info | Fig. 7 (digitised) |
| ExperimentA2_IAST | CH4 front time (y_CH4 = y_feed/2) vs experiment (IAST) | experiment | 364 s | 365 s | 0.277 % rel. error | info | Fig. 7 (digitised) |
| ExperimentA2_IAST | thermocouple T (5 positions, pooled) vs experiment (IAST) | experiment | curve | digitised data | 22 K RMSE | info | Fig. 7 (digitised) |
| ExperimentA3 | outlet y_CO2 | experiment | curve | digitised data | 8.26 % NRMSE | 10 % | Fig. 6 (digitised) |
| ExperimentA3 | outlet y_CH4 | experiment | curve | digitised data | 17.7 % NRMSE | 20 % | Fig. 6 (digitised) |
| ExperimentA3 | CO2 front time (y_CO2 = y_feed/2) | experiment | 420.8 s | 414.9 s | 1.43 % rel. error | 2 % | Fig. 6 (digitised) |
| ExperimentA3 | CH4 front time (y_CH4 = y_feed/2) | experiment | 160.3 s | 156.5 s | 2.39 % rel. error | 3 % | Fig. 6 (digitised) |
| ExperimentA3 | thermocouple T (5 positions, pooled) | experiment | curve | digitised data | 5.45 K RMSE | 6 K | Fig. 6 (digitised) |
| ExperimentA3 | CH4 breakthrough time (y &gt; 0.01) | experiment (text, approx.) | 158.9 s | 140 s | 18.9 s abs. error | info | Sect. 4.2.1 text (Fig. 6, reference experiment A3) |
| ExperimentA3 | CO2 breakthrough time (y &gt; 0.01) | experiment (text, approx.) | 407.8 s | 400 s | 7.79 s abs. error | info | Sect. 4.2.1 text (Fig. 6, reference experiment A3) |
| ExperimentA3_IAST | outlet y_CO2 vs experiment (IAST) | experiment | curve | digitised data | 27.8 % NRMSE | info | Fig. 6 (digitised) |
| ExperimentA3_IAST | outlet y_CH4 vs experiment (IAST) | experiment | curve | digitised data | 31.3 % NRMSE | info | Fig. 6 (digitised) |
| ExperimentA3_IAST | CO2 front time (y_CO2 = y_feed/2) vs experiment (IAST) | experiment | 456.5 s | 414.9 s | 10 % rel. error | info | Fig. 6 (digitised) |
| ExperimentA3_IAST | CH4 front time (y_CH4 = y_feed/2) vs experiment (IAST) | experiment | 156.8 s | 156.5 s | 0.194 % rel. error | info | Fig. 6 (digitised) |
| ExperimentA3_IAST | thermocouple T (5 positions, pooled) vs experiment (IAST) | experiment | curve | digitised data | 13.9 K RMSE | info | Fig. 6 (digitised) |
| ExperimentA3_IAST | CH4 breakthrough time (y &gt; 0.01) | experiment (text, approx.) | 155.5 s | 140 s | 15.5 s abs. error | info | Sect. 4.2.1 text (Fig. 6, reference experiment A3) |
| ExperimentA3_IAST | CO2 breakthrough time (y &gt; 0.01) | experiment (text, approx.) | 446.2 s | 400 s | 46.2 s abs. error | info | Sect. 4.2.1 text (Fig. 6, reference experiment A3) |
| ExperimentA4 | outlet y_CO2 | experiment | curve | digitised data | 2.34 % NRMSE | 3 % | Fig. 7 (digitised) |
| ExperimentA4 | outlet y_CH4 | experiment | curve | digitised data | 28.9 % NRMSE | 32 % | Fig. 7 (digitised) |
| ExperimentA4 | CO2 front time (y_CO2 = y_feed/2) | experiment | 225.4 s | 226.2 s | 0.339 % rel. error | 0.5 % | Fig. 7 (digitised) |
| ExperimentA4 | CH4 front time (y_CH4 = y_feed/2) | experiment | 86.25 s | 85.66 s | 0.691 % rel. error | 1 % | Fig. 7 (digitised) |
| ExperimentA4 | thermocouple T (5 positions, pooled) | experiment | curve | digitised data | 7.07 K RMSE | 8 K | Fig. 7 (digitised) |
| ExperimentA4_IAST | outlet y_CO2 vs experiment (IAST) | experiment | curve | digitised data | 23.6 % NRMSE | info | Fig. 7 (digitised) |
| ExperimentA4_IAST | outlet y_CH4 vs experiment (IAST) | experiment | curve | digitised data | 51.5 % NRMSE | info | Fig. 7 (digitised) |
| ExperimentA4_IAST | CO2 front time (y_CO2 = y_feed/2) vs experiment (IAST) | experiment | 242.7 s | 226.2 s | 7.3 % rel. error | info | Fig. 7 (digitised) |
| ExperimentA4_IAST | CH4 front time (y_CH4 = y_feed/2) vs experiment (IAST) | experiment | 84.13 s | 85.66 s | 1.79 % rel. error | info | Fig. 7 (digitised) |
| ExperimentA4_IAST | thermocouple T (5 positions, pooled) vs experiment (IAST) | experiment | curve | digitised data | 9 K RMSE | info | Fig. 7 (digitised) |
| ExperimentA5 | outlet y_CO2 | experiment | curve | digitised data | 4.99 % NRMSE | 6 % | Fig. 7 (digitised) |
| ExperimentA5 | outlet y_CH4 | experiment | curve | digitised data | 31.2 % NRMSE | 35 % | Fig. 7 (digitised) |
| ExperimentA5 | CO2 front time (y_CO2 = y_feed/2) | experiment | 186.8 s | 189.5 s | 1.4 % rel. error | 2 % | Fig. 7 (digitised) |
| ExperimentA5 | CH4 front time (y_CH4 = y_feed/2) | experiment | 94.02 s | 94.98 s | 1.01 % rel. error | 1.5 % | Fig. 7 (digitised) |
| ExperimentA5 | thermocouple T (5 positions, pooled) | experiment | curve | digitised data | 11.4 K RMSE | 13 K | Fig. 7 (digitised) |
| ExperimentA5_IAST | outlet y_CO2 vs experiment (IAST) | experiment | curve | digitised data | 27.3 % NRMSE | info | Fig. 7 (digitised) |
| ExperimentA5_IAST | outlet y_CH4 vs experiment (IAST) | experiment | curve | digitised data | 49.6 % NRMSE | info | Fig. 7 (digitised) |
| ExperimentA5_IAST | CO2 front time (y_CO2 = y_feed/2) vs experiment (IAST) | experiment | 210.1 s | 189.5 s | 10.9 % rel. error | info | Fig. 7 (digitised) |
| ExperimentA5_IAST | CH4 front time (y_CH4 = y_feed/2) vs experiment (IAST) | experiment | 92.26 s | 94.98 s | 2.87 % rel. error | info | Fig. 7 (digitised) |
| ExperimentA5_IAST | thermocouple T (5 positions, pooled) vs experiment (IAST) | experiment | curve | digitised data | 14.3 K RMSE | info | Fig. 7 (digitised) |
| ExperimentA6 | outlet y_CO2 | experiment | curve | digitised data | 1.32 % NRMSE | 2 % | Fig. 7 (digitised) |
| ExperimentA6 | outlet y_CH4 | experiment | curve | digitised data | 11.8 % NRMSE | 13 % | Fig. 7 (digitised) |
| ExperimentA6 | CO2 front time (y_CO2 = y_feed/2) | experiment | 447 s | 447.9 s | 0.186 % rel. error | 0.5 % | Fig. 7 (digitised) |
| ExperimentA6 | CH4 front time (y_CH4 = y_feed/2) | experiment | 160.8 s | 160.6 s | 0.172 % rel. error | 0.5 % | Fig. 7 (digitised) |
| ExperimentA6 | thermocouple T (5 positions, pooled) | experiment | curve | digitised data | 7.46 K RMSE | 8.5 K | Fig. 7 (digitised) |
| ExperimentA6_IAST | outlet y_CO2 vs experiment (IAST) | experiment | curve | digitised data | 30 % NRMSE | info | Fig. 7 (digitised) |
| ExperimentA6_IAST | outlet y_CH4 vs experiment (IAST) | experiment | curve | digitised data | 45.9 % NRMSE | info | Fig. 7 (digitised) |
| ExperimentA6_IAST | CO2 front time (y_CO2 = y_feed/2) vs experiment (IAST) | experiment | 477.8 s | 447.9 s | 6.69 % rel. error | info | Fig. 7 (digitised) |
| ExperimentA6_IAST | CH4 front time (y_CH4 = y_feed/2) vs experiment (IAST) | experiment | 157.3 s | 160.6 s | 2.01 % rel. error | info | Fig. 7 (digitised) |
| ExperimentA6_IAST | thermocouple T (5 positions, pooled) vs experiment (IAST) | experiment | curve | digitised data | 13 K RMSE | info | Fig. 7 (digitised) |
| ExperimentA7 | outlet y_CO2 | experiment | curve | digitised data | 1.38 % NRMSE | 2 % | Fig. 7 (digitised) |
| ExperimentA7 | outlet y_CH4 | experiment | curve | digitised data | 25.9 % NRMSE | 29 % | Fig. 7 (digitised) |
| ExperimentA7 | CO2 front time (y_CO2 = y_feed/2) | experiment | 424.7 s | 426.5 s | 0.406 % rel. error | 0.5 % | Fig. 7 (digitised) |
| ExperimentA7 | CH4 front time (y_CH4 = y_feed/2) | experiment | 145.8 s | 146.2 s | 0.286 % rel. error | 0.5 % | Fig. 7 (digitised) |
| ExperimentA7 | thermocouple T (5 positions, pooled) | experiment | curve | digitised data | 7.16 K RMSE | 8 K | Fig. 7 (digitised) |
| ExperimentA7_IAST | outlet y_CO2 vs experiment (IAST) | experiment | curve | digitised data | 20.6 % NRMSE | info | Fig. 7 (digitised) |
| ExperimentA7_IAST | outlet y_CH4 vs experiment (IAST) | experiment | curve | digitised data | 51 % NRMSE | info | Fig. 7 (digitised) |
| ExperimentA7_IAST | CO2 front time (y_CO2 = y_feed/2) vs experiment (IAST) | experiment | 448.5 s | 426.5 s | 5.18 % rel. error | info | Fig. 7 (digitised) |
| ExperimentA7_IAST | CH4 front time (y_CH4 = y_feed/2) vs experiment (IAST) | experiment | 142.6 s | 146.2 s | 2.48 % rel. error | info | Fig. 7 (digitised) |
| ExperimentA7_IAST | thermocouple T (5 positions, pooled) vs experiment (IAST) | experiment | curve | digitised data | 9.86 K RMSE | info | Fig. 7 (digitised) |
| ExperimentB1 | outlet y_CO2 | experiment | curve | digitised data | 4.51 % NRMSE | 5 % | Fig. 11 (digitised) |
| ExperimentB1 | outlet y_CH4 | experiment | curve | digitised data | 6.51 % NRMSE | 8 % | Fig. 11 (digitised) |
| ExperimentB1 | CO2 front time (y_CO2 = y_feed/2) | experiment | 8466 s | 8474 s | 0.0969 % rel. error | 0.5 % | Fig. 11 (digitised) |
| ExperimentB1 | CH4 front time (y_CH4 = y_feed/2) | experiment | 321.8 s | 337.7 s | 4.7 % rel. error | 5.5 % | Fig. 11 (digitised) |
| ExperimentB1 | outlet y_H2 | experiment | curve | digitised data | 24.9 % NRMSE | 28 % | Fig. 11 (digitised) |
| ExperimentB1 | thermocouple T (5 positions, pooled) | experiment | curve | digitised data | 7.37 K RMSE | 8.5 K | Fig. 11 (digitised) |
| ExperimentB1_IAST | outlet y_CO2 vs experiment (IAST) | experiment | curve | digitised data | 8.93 % NRMSE | info | Fig. 11 (digitised) |
| ExperimentB1_IAST | outlet y_CH4 vs experiment (IAST) | experiment | curve | digitised data | 6.9 % NRMSE | info | Fig. 11 (digitised) |
| ExperimentB1_IAST | CO2 front time (y_CO2 = y_feed/2) vs experiment (IAST) | experiment | 8499 s | 8474 s | 0.291 % rel. error | info | Fig. 11 (digitised) |
| ExperimentB1_IAST | CH4 front time (y_CH4 = y_feed/2) vs experiment (IAST) | experiment | 321.8 s | 337.7 s | 4.7 % rel. error | info | Fig. 11 (digitised) |
| ExperimentB1_IAST | outlet y_H2 vs experiment (IAST) | experiment | curve | digitised data | 25 % NRMSE | info | Fig. 11 (digitised) |
| ExperimentB1_IAST | thermocouple T (5 positions, pooled) vs experiment (IAST) | experiment | curve | digitised data | 6.59 K RMSE | info | Fig. 11 (digitised) |
| ExperimentB2 | outlet y_CO2 | experiment | curve | digitised data | 4.96 % NRMSE | 6 % | Fig. 11 (digitised) |
| ExperimentB2 | outlet y_CH4 | experiment | curve | digitised data | 4.68 % NRMSE | 6 % | Fig. 11 (digitised) |
| ExperimentB2 | CO2 front time (y_CO2 = y_feed/2) | experiment | 2693 s | 2694 s | 0.0639 % rel. error | 0.5 % | Fig. 11 (digitised) |
| ExperimentB2 | CH4 front time (y_CH4 = y_feed/2) | experiment | 725.2 s | 713.3 s | 1.67 % rel. error | 2 % | Fig. 11 (digitised) |
| ExperimentB2 | outlet y_H2 | experiment | curve | digitised data | 21.7 % NRMSE | 24 % | Fig. 11 (digitised) |
| ExperimentB2 | thermocouple T (5 positions, pooled) | experiment | curve | digitised data | 9.13 K RMSE | 10.5 K | Fig. 11 (digitised) |
| ExperimentB2_IAST | outlet y_CO2 vs experiment (IAST) | experiment | curve | digitised data | 20.9 % NRMSE | info | Fig. 11 (digitised) |
| ExperimentB2_IAST | outlet y_CH4 vs experiment (IAST) | experiment | curve | digitised data | 9.43 % NRMSE | info | Fig. 11 (digitised) |
| ExperimentB2_IAST | CO2 front time (y_CO2 = y_feed/2) vs experiment (IAST) | experiment | 2755 s | 2694 s | 2.26 % rel. error | info | Fig. 11 (digitised) |
| ExperimentB2_IAST | CH4 front time (y_CH4 = y_feed/2) vs experiment (IAST) | experiment | 715.6 s | 713.3 s | 0.33 % rel. error | info | Fig. 11 (digitised) |
| ExperimentB2_IAST | outlet y_H2 vs experiment (IAST) | experiment | curve | digitised data | 22.2 % NRMSE | info | Fig. 11 (digitised) |
| ExperimentB2_IAST | thermocouple T (5 positions, pooled) vs experiment (IAST) | experiment | curve | digitised data | 16.2 K RMSE | info | Fig. 11 (digitised) |

??? note "Notes on these comparisons"

    - RMSE over 13 digitised points; normalised by 0.5 -
    - RMSE over 17 digitised points; normalised by 0.5 -
    - first time the curve reaches 0.25 (reference time interpolated from the digitised points)
    - RMSE pooled over 179 digitised points of the centre-line thermocouples at z = 10, 35, 60, 85, 110 cm
    - RMSE over 21 digitised points; normalised by 0.5 -
    - RMSE over 27 digitised points; normalised by 0.5 -
    - RMSE pooled over 192 digitised points of the centre-line thermocouples at z = 10, 35, 60, 85, 110 cm
    - RMSE over 26 digitised points; normalised by 0.5 -
    - RMSE over 33 digitised points; normalised by 0.5 -
    - RMSE pooled over 196 digitised points of the centre-line thermocouples at z = 10, 35, 60, 85, 110 cm
    - approximate experimental value: 'After 140 s ... CH4 breaks through'; simulated: first time outlet y_CH4 &gt; 0.01
    - approximate experimental value: 'After 400 s ... CO2 breakthrough occurs'; simulated: first time outlet y_CO2 &gt; 0.01
    - RMSE over 19 digitised points; normalised by 0.5 -
    - RMSE over 31 digitised points; normalised by 0.5 -
    - RMSE over 34 digitised points; normalised by 0.5 -
    - RMSE pooled over 168 digitised points of the centre-line thermocouples at z = 10, 35, 60, 85, 110 cm
    - RMSE pooled over 166 digitised points of the centre-line thermocouples at z = 10, 35, 60, 85, 110 cm
    - RMSE over 28 digitised points; normalised by 0.5 -
    - RMSE pooled over 154 digitised points of the centre-line thermocouples at z = 10, 35, 60, 85, 110 cm
    - RMSE over 14 digitised points; normalised by 0.2 -
    - RMSE over 13 digitised points; normalised by 0.05 -
    - first time the curve reaches 0.1 (reference time interpolated from the digitised points)
    - first time the curve reaches 0.025 (reference time interpolated from the digitised points)
    - RMSE over 24 digitised points; normalised by 0.75 -
    - RMSE pooled over 92 digitised points of the centre-line thermocouples at z = 10, 35, 60, 85, 110 cm
    - RMSE over 10 digitised points; normalised by 0.2 -
    - RMSE over 17 digitised points; normalised by 0.05 -
    - RMSE pooled over 160 digitised points of the centre-line thermocouples at z = 10, 35, 60, 85, 110 cm

Input decks, assumptions and known discrepancies: [`core/test/papers/streb_2021/README.md`](https://github.com/samaffleck/skarstrom/blob/main/core/test/papers/streb_2021/README.md).

## Vysyaraju et al. (2024) {#vysyarajua-2024}

R. Vysyaraju, J.A. Sawada, L. Estupiñan Perez, A. Rajendran, Simulation and experimental demonstration of helium purification from He/N2 mixtures by pressure swing adsorption, Chemical Engineering Journal 499 (2024) 155929. [doi:10.1016/j.cej.2024.155929](https://doi.org/10.1016/j.cej.2024.155929)

*System:* He/N2 on zeolite 5A: breakthrough and 5-step PSA. *Compared against:* experimental values stated in the text and Table 3.

**Findings.**

- The N2 breakthrough loadings match the reported value (ambient) and the DSL value (11.2 bar) to 2-4 %, and the peak bed temperature to 0.3 K.
- He purity is 0.54 pp below the paper model (1.4 pp above the experiment), and He recovery is 63 % against 31 %: the deck is a single bed with self-equalisation, the experiment and the paper's model are a 3-bed cycle.

| Case | Quantity | Against | Simulated | Reference | Error | Tolerance | Source |
|---|---|---|---|---|---|---|---|
| DCB_Ambient_N2 | N2 equilibrium loading | experiment | 0.5381 mol/kg | 0.56 mol/kg | 3.92 % rel. error | 4.9 % | Sect. 4.2 text (Fig. 4a) |
| DCB_Ambient_N2 | peak bed temperature | experiment | 304.3 K | 304 K | 0.251 K abs. error | 0.56 K | Sect. 4.2 text (Fig. 4) |
| DCB_HighPressure_N2 | N2 equilibrium loading | isotherm (derived) | 2.126 mol/kg | 2.167 mol/kg | 1.89 % rel. error | 2.9 % | Table 2 DSL parameters evaluated at 11.2 bar, 298 K (Sect. 4.2, Fig. 5b) |
| PSA_Exp1 | He purity vs experiment | experiment | 99.93 % | 98 % | 1.93 pp abs. error | 2.4 pp | Table 3, Exp. 1, Exp Pu_He |
| PSA_Exp1 | He purity vs paper model | paper model | 99.93 % | 99.9 % | 0.0318 pp abs. error | 0.084 pp | Table 3, Exp. 1, Model Pu_He |
| PSA_Exp1 | He recovery vs experiment | experiment | 64.08 % | 30.6 % | 33.5 pp abs. error | 36 pp | Table 3, Exp. 1, Exp Re_He |
| PSA_Exp1 | He recovery vs paper model | paper model | 64.08 % | 31.3 % | 32.8 pp abs. error | 35 pp | Table 3, Exp. 1, Model Re_He |

??? note "Notes on these comparisons"

    - N2 equilibrium loading from the mass balance of the 100 % N2 adsorption DCB at 0.95 bar, 295 K
    - 'The temperature inside the bed reached a maximum of 304 K from 295 K during adsorption'
    - DERIVED, not printed in the paper: the paper states the measured high-pressure DCB loading 'was within 5% of the loading obtained from the DSL model'; this is that DSL loading computed from Table 2
    - measured at CSS (3-bed rig)
    - paper's detailed model
    - measured at CSS (3-bed rig); single-bed model with self-equalisation vs the paper's 3-bed rig; recovery by He conservation
    - paper's detailed model; single-bed model with self-equalisation vs the paper's 3-bed model

Input decks, assumptions and known discrepancies: [`core/test/papers/vysyarajua_2024/README.md`](https://github.com/samaffleck/skarstrom/blob/main/core/test/papers/vysyarajua_2024/README.md).

## Ward & Pini (2022) {#ward-2022}

A. Ward, R. Pini, Integrated uncertainty quantification and sensitivity analysis of single-component dynamic column breakthrough experiments, Adsorption 28 (2022) 161-183. [doi:10.1007/s10450-022-00361-z](https://doi.org/10.1007/s10450-022-00361-z)

*System:* CO2 breakthrough on activated carbon at 50/100/200 ccm. *Compared against:* digitised experimental outlet composition and thermocouple.

**Findings.**

- The decks use the paper's model forms (WP10 audit, 2026-10-31: helium carrier, the paper's dispersion correlation Eq. 4 with D_m = 0.605e-5 m2/s, -dU_1 + RT as the heat of adsorption, the wall ends held at ambient (Table 4), the Ergun pressure drop) with four values calibrated against the digitised curves of all three flows (2026-11): the ambient and wall-end temperature 293 K (Table 1: 291.4 K; the thermocouple settles at 293 K after the front), q_s,2 5.89 mol/kg (Table 2: 6.04, inside the paper's 95 % interval), the factor 3.75 on Eq. 4 (the paper fits 4) and h_in 97.3 / 81.5 / 69.0 W/(m2 K) (Table 6: 77.5 / 57.9 / 45.2). The CO2 front is within 0.4 % of the measurement, the outlet NRMSE 2.7-3.9 % and the thermocouple RMSE 0.18-1.40 K (before the calibration: 1.7-4.5 % late, 14-35 % and 1.15-2.56 K).
- Attribution (one change at a time, README): the 293 K ambient alone brings the front to within 0.7-1.6 %, q_s,2 closes the rest, the larger h_in takes out the extra heat and the dispersion factor sets the front width. The calibration is fitted to the same digitised curves it is compared with, so it is not an independent validation. Of the audit's changes, the fitted dispersion factor halves the temperature rise and the wall end condition moves the front from early to late; the carrier gas, the RT term and the pressure-drop constant are each worth under 1 % in the front time.
- Every deck value is traced to the paper in the README's provenance table. The deck keeps the particle voidage at 0 in the gas holdup (the paper's balances count the bed voidage only) and uses 60 cells to the paper's 30.

| Case | Quantity | Against | Simulated | Reference | Error | Tolerance | Source |
|---|---|---|---|---|---|---|---|
| Q100ccm | outlet y_CO2 | experiment | curve | digitised data | 2.72 % NRMSE | 3.2 % | Fig. 3 / data of ref. [5] (digitised) |
| Q100ccm | CO2 front time (y_CO2 = 0.5) | experiment | 1259 s | 1259 s | 0.0193 % rel. error | 0.2 % | Fig. 3 / data of ref. [5] (digitised) |
| Q100ccm | thermocouple T at z = 0.0715 m | experiment | curve | digitised data | 0.675 K RMSE | 0.78 K | Fig. 3 / data of ref. [5] (digitised) |
| Q100ccm | breakthrough time (y_CO2,out = 0.5) | paper model (95% CI) | 1259 s | 1250 s | 9.21 s abs. error | info | Table 10, adsorption, f_in = 100 ccm, t_br 95% CI |
| Q100ccm | temperature rise at z = 0.0715 m | paper model (95% CI) | 13.09 K | 15.05 K | 1.96 K abs. error | info | Table 10, adsorption, f_in = 100 ccm, Delta T_ads 95% CI |
| Q200ccm | outlet y_CO2 | experiment | curve | digitised data | 3.89 % NRMSE | 4.5 % | Fig. 3 / data of ref. [5] (digitised) |
| Q200ccm | CO2 front time (y_CO2 = 0.5) | experiment | 593.5 s | 591.1 s | 0.404 % rel. error | 0.47 % | Fig. 3 / data of ref. [5] (digitised) |
| Q200ccm | thermocouple T at z = 0.0715 m | experiment | curve | digitised data | 1.4 K RMSE | 1.6 K | Fig. 3 / data of ref. [5] (digitised) |
| Q200ccm | breakthrough time (y_CO2,out = 0.5) | paper model (95% CI) | 593.5 s | 575 s | 18.5 s abs. error | info | Table 10, adsorption, f_in = 200 ccm, t_br 95% CI |
| Q200ccm | temperature rise at z = 0.0715 m | paper model (95% CI) | 20.87 K | 23.6 K | 2.73 K abs. error | info | Table 10, adsorption, f_in = 200 ccm, Delta T_ads 95% CI |
| Q50ccm | outlet y_CO2 | experiment | curve | digitised data | 2.76 % NRMSE | 3.3 % | Fig. 3 / data of ref. [5] (digitised) |
| Q50ccm | CO2 front time (y_CO2 = 0.5) | experiment | 2574 s | 2574 s | 0.0209 % rel. error | 0.2 % | Fig. 3 / data of ref. [5] (digitised) |
| Q50ccm | thermocouple T at z = 0.0715 m | experiment | curve | digitised data | 0.181 K RMSE | 0.22 K | Fig. 3 / data of ref. [5] (digitised) |
| Q50ccm | breakthrough time (y_CO2,out = 0.5) | paper model (95% CI) | 2574 s | 2540 s | 33.7 s abs. error | info | Table 10, adsorption, f_in = 50 ccm, t_br 95% CI |
| Q50ccm | temperature rise at z = 0.0715 m | paper model (95% CI) | 6.666 K | 8.445 K | 1.78 K abs. error | info | Table 10, adsorption, f_in = 50 ccm, Delta T_ads 95% CI |

??? note "Notes on these comparisons"

    - RMSE over 29 digitised points; normalised by 1 -
    - first time the curve reaches 0.5 (reference time interpolated from the digitised points)
    - RMSE over 37 digitised points
    - 95% credibility interval of the paper's model prediction; reference = CI midpoint; simulated value inside the paper's 95% CI [1048, 1451]; t_br = the time y_CO2,out reaches 0.5 (Table 9)
    - 95% credibility interval of the paper's model prediction; reference = CI midpoint; simulated value inside the paper's 95% CI [12.9, 17.2]; simulated = max T at the thermocouple minus T_F (Table 9)
    - RMSE over 28 digitised points; normalised by 1 -
    - RMSE over 39 digitised points
    - 95% credibility interval of the paper's model prediction; reference = CI midpoint; simulated value inside the paper's 95% CI [482, 668]; t_br = the time y_CO2,out reaches 0.5 (Table 9)
    - 95% credibility interval of the paper's model prediction; reference = CI midpoint; simulated value inside the paper's 95% CI [20, 27.2]; simulated = max T at the thermocouple minus T_F (Table 9)
    - RMSE over 27 digitised points; normalised by 1 -
    - RMSE over 31 digitised points
    - 95% credibility interval of the paper's MODEL prediction (not a measurement); t_br = time y_out reaches 0.5 (Table 9); reference = CI midpoint; simulated value inside the paper's 95% CI [2132, 2949]; t_br = the time y_CO2,out reaches 0.5 (Table 9)
    - 95% credibility interval of the paper's model prediction; reference = CI midpoint; simulated value outside the paper's 95% CI [7.36, 9.53]; simulated = max T at the thermocouple minus T_F (Table 9)

Input decks, assumptions and known discrepancies: [`core/test/papers/ward_2022/README.md`](https://github.com/samaffleck/skarstrom/blob/main/core/test/papers/ward_2022/README.md).

## Ward & Pini (2024) {#ward-pini-2024}

A. Ward, R. Pini, Design and Performance Evaluation of Multisorbent Vacuum-Swing Adsorption Processes for Postcombustion Carbon Capture, Industrial & Engineering Chemistry Research 63 (2024) 13787-13800. [doi:10.1021/acs.iecr.4c01345](https://doi.org/10.1021/acs.iecr.4c01345)

*System:* CO2/N2 4-step VSA on Mg-MOF-74. *Compared against:* paper's optimised detailed-model results (Table 4).

**Findings.**

- CO2 purity within 2.3 pp, recovery within 1.8 pp and productivity within 0.7 % of Table 4.
- The deck now uses the extended (competitive) DSL the paper states (Sect. 2.1). Energy is not compared (no pump model in the deck).

| Case | Quantity | Against | Simulated | Reference | Error | Tolerance | Source |
|---|---|---|---|---|---|---|---|
| FourStepPSA | CO2 purity | paper model | 96.16 % | 98.5 % | 2.34 pp abs. error | 2.7 pp | Table 4, Maximum separation effectiveness, Mg-MOF-74 (single-adsorbent) |
| FourStepPSA | CO2 recovery | paper model | 91.83 % | 89 % | 2.83 pp abs. error | 3.4 pp | Table 4, Maximum separation effectiveness, Mg-MOF-74 (single-adsorbent) |
| FourStepPSA | CO2 productivity | paper model | 0.8434 mol/(m3 s) | 0.828 mol/(m3 s) | 1.86 % rel. error | 2.7 % | Table 4, Maximum separation effectiveness, Mg-MOF-74 (single-adsorbent) |

??? note "Notes on these comparisons"

    - paper's optimised detailed-model result (simulation)
    - paper's detailed-model result
    - CO2 from the evacuation step per bed volume (V_bed = L A) per total cycle time (Sect. 4.2)

Input decks, assumptions and known discrepancies: [`core/test/papers/ward_pini_2024/README.md`](https://github.com/samaffleck/skarstrom/blob/main/core/test/papers/ward_pini_2024/README.md).

## Halabi (2011) {#halabi-2011}

M.H. Halabi, Sorption enhanced catalytic reforming of methane for pure hydrogen production: experimental and modeling, PhD thesis, Eindhoven University of Technology, 2011 (chapters 8 and 9). [doi:10.6100/IR709035](https://doi.org/10.6100/IR709035)

*System:* Sorption-enhanced steam methane reforming: Rh/Ce0.6Zr0.4O2 catalyst mixed with K2CO3-promoted hydrotalcite, 450 C, 1.5-6 bar. *Compared against:* thesis's computed gas-phase equilibria (Fig. 7 legend, Table 6); measured conversion and H2 (Table 6, recorded only).

**Findings.**

- With an inert sorbent the chapter-8 rate model (LHHW, the two site balances expanded into one denominator) reaches the thesis's printed equilibria: Fig. 7 within 0.46 pp, the Table 6 equilibrium column within 1.95 pp. The model equals a direct solution of the two equilibria with the thesis's K values; the Table 6 column lies 0.3-2.0 pp below it.
- The sorption-enhanced runs are recorded, not asserted. The model's conversion has no plateau before breakthrough, and CO2 breaks through after 1.6-3.4 min against the thesis's 4-15 min. A 60 s reading seemed to agree at 1.5 bar, but that is an artefact of the read time.
- The open input is the sorbent capacity below 0.05 bar CO2: the thesis's Freundlich fit (0.13 mol/kg at 0.024 bar), its Table 3 and its Henry constant disagree by up to 2x (skarstrom#245).

| Case | Quantity | Against | Simulated | Reference | Error | Tolerance | Source |
|---|---|---|---|---|---|---|---|
| GasPhaseEquilibrium | CH4 equilibrium conversion vs paper model, Fig. 7 feed, 1.5 bar | paper model | 57.15 % | 57.1 % | 0.0528 pp abs. error | 1 pp | Ch. 9, Fig. 7 legend (p. 205): 'Equilibrium conversion' |
| GasPhaseEquilibrium | CH4 equilibrium conversion vs paper model, Fig. 7 feed, 4.5 bar | paper model | 40.14 % | 40.1 % | 0.0372 pp abs. error | 1 pp | Ch. 9, Fig. 7 legend (p. 205): 'Equilibrium conversion' |
| GasPhaseEquilibrium | CH4 equilibrium conversion vs paper model, Fig. 7 feed, 6 bar | paper model | 36.34 % | 36.8 % | 0.455 pp abs. error | 1 pp | Ch. 9, Fig. 7 legend (p. 205): 'Equilibrium conversion' |
| GasPhaseEquilibrium | CH4 equilibrium conversion vs paper model, Table 6, S/C 3, 1.5 bar | paper model | 47.25 % | 45.3 % | 1.95 pp abs. error | 2.5 pp | Ch. 9, Table 6 (p. 204), CH4 equil. conversion |
| GasPhaseEquilibrium | CH4 equilibrium conversion vs paper model, Table 6, S/C 3, 4.5 bar | paper model | 32.71 % | 32.1 % | 0.61 pp abs. error | 2.5 pp | Ch. 9, Table 6 (p. 204), CH4 equil. conversion |
| GasPhaseEquilibrium | CH4 equilibrium conversion vs paper model, Table 6, S/C 3, 6 bar | paper model | 29.55 % | 29.2 % | 0.349 pp abs. error | 2.5 pp | Ch. 9, Table 6 (p. 204), CH4 equil. conversion |
| GasPhaseEquilibrium | CH4 equilibrium conversion vs paper model, Table 6, S/C 4, 1.5 bar | paper model | 52.86 % | 51.2 % | 1.66 pp abs. error | 2.5 pp | Ch. 9, Table 6 (p. 204), CH4 equil. conversion |
| GasPhaseEquilibrium | CH4 equilibrium conversion vs paper model, Table 6, S/C 4, 4.5 bar | paper model | 36.76 % | 36.1 % | 0.656 pp abs. error | 2.5 pp | Ch. 9, Table 6 (p. 204), CH4 equil. conversion |
| GasPhaseEquilibrium | CH4 equilibrium conversion vs paper model, Table 6, S/C 4, 6 bar | paper model | 33.22 % | 31.8 % | 1.42 pp abs. error | 2.5 pp | Ch. 9, Table 6 (p. 204), CH4 equil. conversion |
| GasPhaseEquilibrium | CH4 equilibrium conversion vs paper model, Table 6, S/C 5, 1.5 bar | paper model | 57.41 % | 55.6 % | 1.81 pp abs. error | 2.5 pp | Ch. 9, Table 6 (p. 204), CH4 equil. conversion |
| GasPhaseEquilibrium | CH4 equilibrium conversion vs paper model, Table 6, S/C 5, 4.5 bar | paper model | 40.13 % | 39.6 % | 0.532 pp abs. error | 2.5 pp | Ch. 9, Table 6 (p. 204), CH4 equil. conversion |
| GasPhaseEquilibrium | CH4 equilibrium conversion vs paper model, Table 6, S/C 5, 6 bar | paper model | 36.3 % | 35.9 % | 0.402 pp abs. error | 2.5 pp | Ch. 9, Table 6 (p. 204), CH4 equil. conversion |
| GasPhaseEquilibrium | CH4 equilibrium conversion vs paper model, Table 6, S/C 6, 1.5 bar | paper model | 61.21 % | 59.4 % | 1.81 pp abs. error | 2.5 pp | Ch. 9, Table 6 (p. 204), CH4 equil. conversion |
| GasPhaseEquilibrium | CH4 equilibrium conversion vs paper model, Table 6, S/C 6, 4.5 bar | paper model | 43.04 % | 42.4 % | 0.645 pp abs. error | 2.5 pp | Ch. 9, Table 6 (p. 204), CH4 equil. conversion |
| GasPhaseEquilibrium | CH4 equilibrium conversion vs paper model, Table 6, S/C 6, 6 bar | paper model | 38.97 % | 38.5 % | 0.471 pp abs. error | 2.5 pp | Ch. 9, Table 6 (p. 204), CH4 equil. conversion |
| SorptionEnhanced | CH4 conversion at 60 s vs experiment, S/C 3, 1.5 bar | experiment | 81.25 % | 82 % | 0.753 pp abs. error | info | Ch. 9, Table 6 (p. 204), CH4 conversion |
| SorptionEnhanced | max CH4 conversion (t &gt;= 5 s) vs experiment, S/C 3, 1.5 bar | experiment | 94.15 % | 82 % | 12.2 pp abs. error | info | Ch. 9, Table 6 (p. 204), CH4 conversion |
| SorptionEnhanced | H2 dry Ar-free at 60 s vs experiment, S/C 3, 1.5 bar | experiment | 94.35 % | 94.7 % | 0.349 pp abs. error | info | Ch. 9, Table 6 (p. 204), H2 concentration |
| SorptionEnhanced | CO2 breakthrough time vs experiment, S/C 3, 1.5 bar | experiment | 1.628 min | 4 min | 2.37 min abs. error | info | Ch. 9, p. 203: "average breakthrough times were determined in the range of 4-15 min" |
| SorptionEnhanced | CH4 conversion at 60 s vs experiment, S/C 4, 1.5 bar | experiment | 90.52 % | 92.7 % | 2.18 pp abs. error | info | Ch. 9, Table 6 (p. 204), CH4 conversion |
| SorptionEnhanced | max CH4 conversion (t &gt;= 5 s) vs experiment, S/C 4, 1.5 bar | experiment | 99.15 % | 92.7 % | 6.45 pp abs. error | info | Ch. 9, Table 6 (p. 204), CH4 conversion |
| SorptionEnhanced | H2 dry Ar-free at 60 s vs experiment, S/C 4, 1.5 bar | experiment | 97.3 % | 98 % | 0.7 pp abs. error | info | Ch. 9, Table 6 (p. 204), H2 concentration |
| SorptionEnhanced | CO2 breakthrough time vs experiment, S/C 4, 1.5 bar | experiment | 1.677 min | 4 min | 2.32 min abs. error | info | Ch. 9, p. 203: "average breakthrough times were determined in the range of 4-15 min" |
| SorptionEnhanced | CH4 conversion at 60 s vs experiment, S/C 5, 1.5 bar | experiment | 96.59 % | 98 % | 1.41 pp abs. error | info | Ch. 9, Table 6 (p. 204), CH4 conversion |
| SorptionEnhanced | max CH4 conversion (t &gt;= 5 s) vs experiment, S/C 5, 1.5 bar | experiment | 99.91 % | 98 % | 1.91 pp abs. error | info | Ch. 9, Table 6 (p. 204), CH4 conversion |
| SorptionEnhanced | H2 dry Ar-free at 60 s vs experiment, S/C 5, 1.5 bar | experiment | 99.05 % | 99.5 % | 0.449 pp abs. error | info | Ch. 9, Table 6 (p. 204), H2 concentration |
| SorptionEnhanced | CO2 breakthrough time vs experiment, S/C 5, 1.5 bar | experiment | 1.792 min | 4 min | 2.21 min abs. error | info | Ch. 9, p. 203: "average breakthrough times were determined in the range of 4-15 min" |
| SorptionEnhanced | CH4 conversion at 60 s vs experiment, S/C 6, 1.5 bar | experiment | 99.28 % | 99.6 % | 0.321 pp abs. error | info | Ch. 9, Table 6 (p. 204), CH4 conversion |
| SorptionEnhanced | max CH4 conversion (t &gt;= 5 s) vs experiment, S/C 6, 1.5 bar | experiment | 99.97 % | 99.6 % | 0.366 pp abs. error | info | Ch. 9, Table 6 (p. 204), CH4 conversion |
| SorptionEnhanced | H2 dry Ar-free at 60 s vs experiment, S/C 6, 1.5 bar | experiment | 99.8 % | 99.8 % | 0.00317 pp abs. error | info | Ch. 9, Table 6 (p. 204), H2 concentration |
| SorptionEnhanced | CO2 breakthrough time vs experiment, S/C 6, 1.5 bar | experiment | 1.923 min | 4 min | 2.08 min abs. error | info | Ch. 9, p. 203: "average breakthrough times were determined in the range of 4-15 min" |
| SorptionEnhanced | CH4 conversion at 60 s vs experiment, S/C 3, 4.5 bar | experiment | 72.1 % | 74.9 % | 2.8 pp abs. error | info | Ch. 9, Table 6 (p. 204), CH4 conversion |
| SorptionEnhanced | max CH4 conversion (t &gt;= 5 s) vs experiment, S/C 3, 4.5 bar | experiment | 83.67 % | 74.9 % | 8.77 pp abs. error | info | Ch. 9, Table 6 (p. 204), CH4 conversion |
| SorptionEnhanced | H2 dry Ar-free at 60 s vs experiment, S/C 3, 4.5 bar | experiment | 91.18 % | 91.9 % | 0.724 pp abs. error | info | Ch. 9, Table 6 (p. 204), H2 concentration |
| SorptionEnhanced | CO2 breakthrough time vs experiment, S/C 3, 4.5 bar | experiment | 2.922 min | 4 min | 1.08 min abs. error | info | Ch. 9, p. 203: "average breakthrough times were determined in the range of 4-15 min" |
| SorptionEnhanced | CH4 conversion at 60 s vs experiment, S/C 4, 4.5 bar | experiment | 80.27 % | 85.2 % | 4.93 pp abs. error | info | Ch. 9, Table 6 (p. 204), CH4 conversion |
| SorptionEnhanced | max CH4 conversion (t &gt;= 5 s) vs experiment, S/C 4, 4.5 bar | experiment | 91.04 % | 85.2 % | 5.84 pp abs. error | info | Ch. 9, Table 6 (p. 204), CH4 conversion |
| SorptionEnhanced | H2 dry Ar-free at 60 s vs experiment, S/C 4, 4.5 bar | experiment | 94.19 % | 95.6 % | 1.41 pp abs. error | info | Ch. 9, Table 6 (p. 204), H2 concentration |
| SorptionEnhanced | CO2 breakthrough time vs experiment, S/C 4, 4.5 bar | experiment | 2.854 min | 4 min | 1.15 min abs. error | info | Ch. 9, p. 203: "average breakthrough times were determined in the range of 4-15 min" |
| SorptionEnhanced | CH4 conversion at 60 s vs experiment, S/C 5, 4.5 bar | experiment | 86.39 % | 93.4 % | 7.01 pp abs. error | info | Ch. 9, Table 6 (p. 204), CH4 conversion |
| SorptionEnhanced | max CH4 conversion (t &gt;= 5 s) vs experiment, S/C 5, 4.5 bar | experiment | 95.28 % | 93.4 % | 1.88 pp abs. error | info | Ch. 9, Table 6 (p. 204), CH4 conversion |
| SorptionEnhanced | H2 dry Ar-free at 60 s vs experiment, S/C 5, 4.5 bar | experiment | 96.18 % | 98.2 % | 2.02 pp abs. error | info | Ch. 9, Table 6 (p. 204), H2 concentration |
| SorptionEnhanced | CO2 breakthrough time vs experiment, S/C 5, 4.5 bar | experiment | 2.864 min | 4 min | 1.14 min abs. error | info | Ch. 9, p. 203: "average breakthrough times were determined in the range of 4-15 min" |
| SorptionEnhanced | CH4 conversion at 60 s vs experiment, S/C 6, 4.5 bar | experiment | 90.99 % | 97.5 % | 6.51 pp abs. error | info | Ch. 9, Table 6 (p. 204), CH4 conversion |
| SorptionEnhanced | max CH4 conversion (t &gt;= 5 s) vs experiment, S/C 6, 4.5 bar | experiment | 97.6 % | 97.5 % | 0.0954 pp abs. error | info | Ch. 9, Table 6 (p. 204), CH4 conversion |
| SorptionEnhanced | H2 dry Ar-free at 60 s vs experiment, S/C 6, 4.5 bar | experiment | 97.55 % | 99.2 % | 1.65 pp abs. error | info | Ch. 9, Table 6 (p. 204), H2 concentration |
| SorptionEnhanced | CO2 breakthrough time vs experiment, S/C 6, 4.5 bar | experiment | 2.91 min | 4 min | 1.09 min abs. error | info | Ch. 9, p. 203: "average breakthrough times were determined in the range of 4-15 min" |
| SorptionEnhanced | CH4 conversion at 60 s vs experiment, S/C 3, 6 bar | experiment | 69.12 % | 72 % | 2.88 pp abs. error | info | Ch. 9, Table 6 (p. 204), CH4 conversion |
| SorptionEnhanced | max CH4 conversion (t &gt;= 5 s) vs experiment, S/C 3, 6 bar | experiment | 79.82 % | 72 % | 7.82 pp abs. error | info | Ch. 9, Table 6 (p. 204), CH4 conversion |
| SorptionEnhanced | H2 dry Ar-free at 60 s vs experiment, S/C 3, 6 bar | experiment | 89.98 % | 91.3 % | 1.32 pp abs. error | info | Ch. 9, Table 6 (p. 204), H2 concentration |
| SorptionEnhanced | CO2 breakthrough time vs experiment, S/C 3, 6 bar | experiment | 3.379 min | 4 min | 0.621 min abs. error | info | Ch. 9, p. 203: "average breakthrough times were determined in the range of 4-15 min" |
| SorptionEnhanced | CH4 conversion at 60 s vs experiment, S/C 4, 6 bar | experiment | 76.91 % | 83.4 % | 6.49 pp abs. error | info | Ch. 9, Table 6 (p. 204), CH4 conversion |
| SorptionEnhanced | max CH4 conversion (t &gt;= 5 s) vs experiment, S/C 4, 6 bar | experiment | 87.28 % | 83.4 % | 3.88 pp abs. error | info | Ch. 9, Table 6 (p. 204), CH4 conversion |
| SorptionEnhanced | H2 dry Ar-free at 60 s vs experiment, S/C 4, 6 bar | experiment | 93.02 % | 95.3 % | 2.28 pp abs. error | info | Ch. 9, Table 6 (p. 204), H2 concentration |
| SorptionEnhanced | CO2 breakthrough time vs experiment, S/C 4, 6 bar | experiment | 3.283 min | 4 min | 0.717 min abs. error | info | Ch. 9, p. 203: "average breakthrough times were determined in the range of 4-15 min" |
| SorptionEnhanced | CH4 conversion at 60 s vs experiment, S/C 5, 6 bar | experiment | 82.8 % | 91.4 % | 8.6 pp abs. error | info | Ch. 9, Table 6 (p. 204), CH4 conversion |
| SorptionEnhanced | max CH4 conversion (t &gt;= 5 s) vs experiment, S/C 5, 6 bar | experiment | 92 % | 91.4 % | 0.598 pp abs. error | info | Ch. 9, Table 6 (p. 204), CH4 conversion |
| SorptionEnhanced | H2 dry Ar-free at 60 s vs experiment, S/C 5, 6 bar | experiment | 95.05 % | 97.5 % | 2.45 pp abs. error | info | Ch. 9, Table 6 (p. 204), H2 concentration |
| SorptionEnhanced | CO2 breakthrough time vs experiment, S/C 5, 6 bar | experiment | 3.282 min | 4 min | 0.718 min abs. error | info | Ch. 9, p. 203: "average breakthrough times were determined in the range of 4-15 min" |
| SorptionEnhanced | CH4 conversion at 60 s vs experiment, S/C 6, 6 bar | experiment | 87.36 % | 95.2 % | 7.84 pp abs. error | info | Ch. 9, Table 6 (p. 204), CH4 conversion |
| SorptionEnhanced | max CH4 conversion (t &gt;= 5 s) vs experiment, S/C 6, 6 bar | experiment | 95.01 % | 95.2 % | 0.193 pp abs. error | info | Ch. 9, Table 6 (p. 204), CH4 conversion |
| SorptionEnhanced | H2 dry Ar-free at 60 s vs experiment, S/C 6, 6 bar | experiment | 96.49 % | 98.8 % | 2.31 pp abs. error | info | Ch. 9, Table 6 (p. 204), H2 concentration |
| SorptionEnhanced | CO2 breakthrough time vs experiment, S/C 6, 6 bar | experiment | 3.29 min | 4 min | 0.71 min abs. error | info | Ch. 9, p. 203: "average breakthrough times were determined in the range of 4-15 min" |

??? note "Notes on these comparisons"

    - thesis's computed equilibrium, feed 3 vol% CH4, 0.1 vol% H2, S/C 4, 75 nml/min, 450 C (caption and p. 205)
    - thesis's computed equilibrium, S/C 3, 1.5 bar, 450 C, feed 3.94 vol% CH4, 0.1 vol% H2 (p. 203)
    - thesis's computed equilibrium, S/C 3, 4.5 bar, 450 C, feed 3.94 vol% CH4, 0.1 vol% H2 (p. 203)
    - thesis's computed equilibrium, S/C 3, 6 bar, 450 C, feed 3.94 vol% CH4, 0.1 vol% H2 (p. 203)
    - thesis's computed equilibrium, S/C 4, 1.5 bar, 450 C, feed 3.94 vol% CH4, 0.1 vol% H2 (p. 203)
    - thesis's computed equilibrium, S/C 4, 4.5 bar, 450 C, feed 3.94 vol% CH4, 0.1 vol% H2 (p. 203)
    - thesis's computed equilibrium, S/C 4, 6 bar, 450 C, feed 3.94 vol% CH4, 0.1 vol% H2 (p. 203)
    - thesis's computed equilibrium, S/C 5, 1.5 bar, 450 C, feed 3.94 vol% CH4, 0.1 vol% H2 (p. 203)
    - thesis's computed equilibrium, S/C 5, 4.5 bar, 450 C, feed 3.94 vol% CH4, 0.1 vol% H2 (p. 203)
    - thesis's computed equilibrium, S/C 5, 6 bar, 450 C, feed 3.94 vol% CH4, 0.1 vol% H2 (p. 203)
    - thesis's computed equilibrium, S/C 6, 1.5 bar, 450 C, feed 3.94 vol% CH4, 0.1 vol% H2 (p. 203)
    - thesis's computed equilibrium, S/C 6, 4.5 bar, 450 C, feed 3.94 vol% CH4, 0.1 vol% H2 (p. 203)
    - thesis's computed equilibrium, S/C 6, 6 bar, 450 C, feed 3.94 vol% CH4, 0.1 vol% H2 (p. 203)
    - measured, sorption-enhanced regime; the thesis does not say when it was read (recorded against the model at 60 s and at its peak); recorded only: no plateau before breakthrough, early breakthrough (README)
    - measured; taken as dry, Ar-free (the thesis reports effluent composition on a dry basis, Fig. 3, p. 201); recorded only: no plateau before breakthrough, early breakthrough (README)
    - low end of the range the thesis gives for all Table 6 runs; no per-row value is printed; recorded only: the thesis gives a range of 4-15 min for all rows; the reference is its low end

Input decks, assumptions and known discrepancies: [`core/test/papers/halabi_2011/README.md`](https://github.com/samaffleck/skarstrom/blob/main/core/test/papers/halabi_2011/README.md).

## Reproducing this page

<!-- no-run: builds the engine and runs the paper tests -->
```bash
cmake --preset linux-release
cmake --build --preset linux-release
ctest --preset linux-release -R "(Arora2021|Casas2012|Ferre2026|Nguyen2022|Nguyen2023|Perez2019|Streb2021|Vysyaraju2024|Ward2022|WardPini2024)" -j2
python scripts/docs/gen_validation_page.py --merge build/linux-release/validation
```

Each test writes `<build>/validation/<paper>__<case>.json`; the script merges them into `core/test/papers/validation_summary.json` and regenerates this page.
