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
- Every input deck is a *prediction* with the paper's published parameters: adsim does not refit transport coefficients or feed velocities per experiment, and it does not model the sampling lines between column and detector. Known causes of disagreement are listed per paper below and in each folder's README.

## Summary

48 cases from 10 publications; 194 asserted comparisons (0 outside their recorded tolerance) and 15 informational ones. Summary generated on 2026-11-07.

| Publication | Cases | Against | Comparisons | Median error | Worst error |
|---|---|---|---|---|---|
| [Arora & Hasan (2021)](#arora-hasan-2021) | 2 | paper model | 4 | 5.35 pp (abs. error) | 46.6 pp (abs. error) |
| [Casas et al. (2012) / Schell et al. (2012)](#casas-2012) | 4 | paper model, isotherm fit | 8 | 3.35 % (rel. error)<br/>0.595 % (max rel. error) | 16.8 % (rel. error)<br/>6.99 % (max rel. error) |
| [Ferré et al. (2026)](#ferre-2026) | 4 | paper model | 20 | 0.787 pp (abs. error)<br/>6.74 % (rel. error) | 11.5 pp (abs. error)<br/>59.4 % (rel. error) |
| [Nguyen et al. (2022)](#nguyen-2022) | 6 | experiment | 12 | 9.71 % (NRMSE)<br/>0.751 % (rel. error) | 13.7 % (NRMSE)<br/>4.81 % (rel. error) |
| [Nguyen et al. (2023)](#nguyen-2023) | 10 | experiment, paper model | 60 | 2.05 pp (abs. error)<br/>8.42 % (rel. error) | 8.39 pp (abs. error)<br/>15.9 % (rel. error) |
| [Estupiñan Perez et al. (2019)](#perez-2019) | 6 | experiment, paper model | 24 | 2.31 pp (abs. error) | 5.97 pp (abs. error) |
| [Streb & Mazzotti (2021)](#streb-2021) | 9 | experiment | 47 | 37.7 % (NRMSE)<br/>6.48 % (rel. error)<br/>16.4 K (RMSE) | 84.1 % (NRMSE)<br/>18 % (rel. error)<br/>28.5 K (RMSE) |
| [Vysyaraju et al. (2024)](#vysyarajua-2024) | 3 | experiment, isotherm (derived), paper model | 7 | 1.89 % (rel. error)<br/>0.251 K (abs. error)<br/>1.93 pp (abs. error) | 3.92 % (rel. error)<br/>0.251 K (abs. error)<br/>33.5 pp (abs. error) |
| [Ward & Pini (2022)](#ward-2022) | 3 | experiment | 9 | 26.5 % (NRMSE)<br/>3.56 % (rel. error)<br/>1.24 K (RMSE) | 35.2 % (NRMSE)<br/>4.48 % (rel. error)<br/>2.56 K (RMSE) |
| [Ward & Pini (2024)](#ward-pini-2024) | 1 | paper model | 3 | 2.39 pp (abs. error)<br/>2.12 % (rel. error) | 2.76 pp (abs. error)<br/>2.12 % (rel. error) |

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
| LiLSX_PVSA | O2 recovery at cycle 50 | paper model | 18.95 % | 65.51 % | 46.6 pp abs. error | 53 pp | Table 4, PVSA cycle with LiLSX adsorbent, column PC_O2 &gt;= 5 L/min, P_O2 &gt;= 90% |

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
| NaY_F | CO2 recovery | paper model | 97.61 % | 97.8 % | 0.19 pp abs. error | 0.25 pp | Table 5, NaY (F), recovery-recovery case, first row |
| NaY_F | CH4 purity (light product) | paper model | 98.76 % | 98.8 % | 0.0412 pp abs. error | 0.097 pp | Table 5, NaY (F), recovery-recovery case, first row |
| NaY_F | CH4 recovery | paper model | 62.47 % | 66 % | 3.53 pp abs. error | 5.6 pp | Table 5, NaY (F), recovery-recovery case, first row |
| NaY_F | CH4 productivity | paper model | 3.104 mol/(m3 s) | 3.2 mol/(m3 s) | 3.02 % rel. error | 4.6 % | Table 5, NaY (F), recovery-recovery case, first row |
| NaY_W | CO2 purity (heavy product) | paper model | 93.79 % | 95 % | 1.21 pp abs. error | 1.8 pp | Table 5, NaY (W), recovery-recovery case, first row |
| NaY_W | CO2 recovery | paper model | 96.48 % | 95.9 % | 0.578 pp abs. error | 0.89 pp | Table 5, NaY (W), recovery-recovery case, first row |
| NaY_W | CH4 purity (light product) | paper model | 98.15 % | 97.6 % | 0.551 pp abs. error | 0.65 pp | Table 5, NaY (W), recovery-recovery case, first row |
| NaY_W | CH4 recovery | paper model | 59.95 % | 65.4 % | 5.45 pp abs. error | 7.6 pp | Table 5, NaY (W), recovery-recovery case, first row |
| NaY_W | CH4 productivity | paper model | 2.891 mol/(m3 s) | 3.1 mol/(m3 s) | 6.74 % rel. error | 8.5 % | Table 5, NaY (W), recovery-recovery case, first row |
| Zeolite13X | CO2 purity (heavy product) | paper model | 93.66 % | 95 % | 1.34 pp abs. error | 1.5 pp | Table 5, 13X, recovery-recovery case, first row |
| Zeolite13X | CO2 recovery | paper model | 95.64 % | 95.9 % | 0.261 pp abs. error | 0.35 pp | Table 5, 13X, recovery-recovery case, first row |
| Zeolite13X | CH4 purity (light product) | paper model | 97.99 % | 98 % | 0.0139 pp abs. error | 0.076 pp | Table 5, 13X, recovery-recovery case, first row |
| Zeolite13X | CH4 recovery | paper model | 43.63 % | 48.6 % | 4.97 pp abs. error | 5.7 pp | Table 5, 13X, recovery-recovery case, first row |
| Zeolite13X | CH4 productivity | paper model | 1.188 mol/(m3 s) | 1.3 mol/(m3 s) | 8.6 % rel. error | 11 % | Table 5, 13X, recovery-recovery case, first row |

??? note "Notes on these comparisons"

    - paper's detailed-model result
    - CH4 per bed volume per cycle time
    - paper's optimised detailed-model result (simulation)
    - CH4 produced per bed volume per total cycle time (Sect. 4.2)

Input decks, assumptions and known discrepancies: [`core/test/papers/ferre_2026/README.md`](https://github.com/samaffleck/skarstrom/blob/main/core/test/papers/ferre_2026/README.md).

## Nguyen et al. (2022) {#nguyen-2022}

T.T.T. Nguyen, J.-B. Lin, G.K.H. Shimizu, A. Rajendran, Separation of CO2 and N2 on a hydrophobic metal organic framework CALF-20, Chemical Engineering Journal 442 (2022) 136263. [doi:10.1016/j.cej.2022.136263](https://doi.org/10.1016/j.cej.2022.136263)

*System:* CO2 breakthrough on CALF-20. *Compared against:* digitised experimental breakthrough curves (Figs. 4, 5).

**Findings.**

- Outlet CO2 curves agree with the digitised Figs. 4 and 5 to an NRMSE of 8-12 % and the CO2 front position (half the feed fraction) to 0.7-5 %.
- The paper text is not in the repo. The four additional feed fractions assume the same flow as the bundled case of the same figure (50 ccm for Fig. 4, 100 ccm for Fig. 5); Fig. 5 is described as a CO2/N2/He experiment but the deck feeds CO2 in He. Particle diameter and ambient temperature are assumptions (see the README).

| Case | Quantity | Against | Simulated | Reference | Error | Tolerance | Source |
|---|---|---|---|---|---|---|---|
| Fig4_BinaryBreakthrough | outlet y_CO2 | experiment | curve | digitised data | 8.9 % NRMSE | 10 % | Fig. 4 (digitised) |
| Fig4_BinaryBreakthrough | CO2 front time (y_CO2 = y_feed/2) | experiment | 1734 s | 1747 s | 0.751 % rel. error | 2.7 % | Fig. 4 (digitised) |
| Fig4_Binary_100pct | outlet y_CO2 | experiment | curve | digitised data | 10.4 % NRMSE | 12 % | Fig. 4 (digitised) |
| Fig4_Binary_100pct | CO2 front time (y_CO2 = y_feed/2) | experiment | 1251 s | 1258 s | 0.58 % rel. error | 1.7 % | Fig. 4 (digitised) |
| Fig4_Binary_50pct | outlet y_CO2 | experiment | curve | digitised data | 9.71 % NRMSE | 14 % | Fig. 4 (digitised) |
| Fig4_Binary_50pct | CO2 front time (y_CO2 = y_feed/2) | experiment | 2364 s | 2441 s | 3.15 % rel. error | 6 % | Fig. 4 (digitised) |
| Fig5_TernaryBreakthrough | outlet y_CO2 | experiment | curve | digitised data | 8.49 % NRMSE | 9.7 % | Fig. 5 (digitised) |
| Fig5_TernaryBreakthrough | CO2 front time (y_CO2 = y_feed/2) | experiment | 3202 s | 3364 s | 4.81 % rel. error | 5.9 % | Fig. 5 (digitised) |
| Fig5_Ternary_30pct | outlet y_CO2 | experiment | curve | digitised data | 10.1 % NRMSE | 13 % | Fig. 5 (digitised) |
| Fig5_Ternary_30pct | CO2 front time (y_CO2 = y_feed/2) | experiment | 1728 s | 1708 s | 1.15 % rel. error | 3.9 % | Fig. 5 (digitised) |
| Fig5_Ternary_5pct | outlet y_CO2 | experiment | curve | digitised data | 13.7 % NRMSE | 15 % | Fig. 5 (digitised) |
| Fig5_Ternary_5pct | CO2 front time (y_CO2 = y_feed/2) | experiment | 7401 s | 7397 s | 0.0572 % rel. error | 0.26 % | Fig. 5 (digitised) |

??? note "Notes on these comparisons"

    - RMSE over 21 digitised points; normalised by 0.7 -
    - first time the curve reaches 0.35 (reference time interpolated from the digitised points)
    - RMSE over 15 digitised points; normalised by 1 -
    - first time the curve reaches 0.5 (reference time interpolated from the digitised points)
    - RMSE over 24 digitised points; normalised by 0.5 -
    - first time the curve reaches 0.25 (reference time interpolated from the digitised points)
    - RMSE over 31 digitised points; normalised by 0.15 -
    - first time the curve reaches 0.075 (reference time interpolated from the digitised points)
    - RMSE over 32 digitised points; normalised by 0.3 -
    - first time the curve reaches 0.15 (reference time interpolated from the digitised points)
    - RMSE over 16 digitised points; normalised by 0.05 -
    - first time the curve reaches 0.025 (reference time interpolated from the digitised points)

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
| ExperimentIII | CO2 recovery vs paper model | paper model | 74.36 % | 73.4 % | 0.955 pp abs. error | 1.2 pp | Table 2, Exp. III, Model Re |
| ExperimentIII | CO2 productivity vs paper model | paper model | 0.4396 mol CO2/(m3 s) | 0.51 mol CO2/(m3 s) | 13.8 % rel. error | 16 % | Table 2, Exp. III, Model Pr |
| ExperimentIV | CO2 purity vs experiment | experiment | 79.15 % | 82.3 % | 3.15 pp abs. error | 3.7 pp | Table 2, Exp. IV, Experimental Pu |
| ExperimentIV | CO2 recovery vs experiment | experiment | 75.22 % | 72.5 % | 2.72 pp abs. error | 3.1 pp | Table 2, Exp. IV, Experimental Re |
| ExperimentIV | CO2 productivity vs experiment | experiment | 0.428 mol CO2/(m3 s) | 0.46 mol CO2/(m3 s) | 6.97 % rel. error | 8.1 % | Table 2, Exp. IV, Experimental Pr |
| ExperimentIV | CO2 purity vs paper model | paper model | 79.15 % | 81.3 % | 2.15 pp abs. error | 2.6 pp | Table 2, Exp. IV, Model Pu |
| ExperimentIV | CO2 recovery vs paper model | paper model | 75.22 % | 74.3 % | 0.921 pp abs. error | 1.2 pp | Table 2, Exp. IV, Model Re |
| ExperimentIV | CO2 productivity vs paper model | paper model | 0.428 mol CO2/(m3 s) | 0.5 mol CO2/(m3 s) | 14.4 % rel. error | 17 % | Table 2, Exp. IV, Model Pr |
| ExperimentIX | CO2 purity vs experiment | experiment | 84.4 % | 85.67 % | 1.27 pp abs. error | 1.6 pp | Table 2, Exp. IX, Experimental Pu |
| ExperimentIX | CO2 recovery vs experiment | experiment | 98.89 % | 90.5 % | 8.39 pp abs. error | 12 pp | Table 2, Exp. IX, Experimental Re |
| ExperimentIX | CO2 productivity vs experiment | experiment | 0.1431 mol CO2/(m3 s) | 0.13 mol CO2/(m3 s) | 10 % rel. error | 14 % | Table 2, Exp. IX, Experimental Pr |
| ExperimentIX | CO2 purity vs paper model | paper model | 84.4 % | 88.3 % | 3.9 pp abs. error | 4.2 pp | Table 2, Exp. IX, Model Pu |
| ExperimentIX | CO2 recovery vs paper model | paper model | 98.89 % | 91.6 % | 7.29 pp abs. error | 11 pp | Table 2, Exp. IX, Model Re |
| ExperimentIX | CO2 productivity vs paper model | paper model | 0.1431 mol CO2/(m3 s) | 0.14 mol CO2/(m3 s) | 2.19 % rel. error | 4.7 % | Table 2, Exp. IX, Model Pr |
| ExperimentV | CO2 purity vs experiment | experiment | 73.55 % | 78.4 % | 4.85 pp abs. error | 5.7 pp | Table 2, Exp. V, Experimental Pu |
| ExperimentV | CO2 recovery vs experiment | experiment | 76.83 % | 74.2 % | 2.63 pp abs. error | 3 pp | Table 2, Exp. V, Experimental Re |
| ExperimentV | CO2 productivity vs experiment | experiment | 0.3805 mol CO2/(m3 s) | 0.4 mol CO2/(m3 s) | 4.89 % rel. error | 6 % | Table 2, Exp. V, Experimental Pr |
| ExperimentV | CO2 purity vs paper model | paper model | 73.55 % | 75.6 % | 2.05 pp abs. error | 2.4 pp | Table 2, Exp. V, Model Pu |
| ExperimentV | CO2 recovery vs paper model | paper model | 76.83 % | 75.4 % | 1.43 pp abs. error | 1.7 pp | Table 2, Exp. V, Model Re |
| ExperimentV | CO2 productivity vs paper model | paper model | 0.3805 mol CO2/(m3 s) | 0.44 mol CO2/(m3 s) | 13.5 % rel. error | 16 % | Table 2, Exp. V, Model Pr |
| ExperimentVI | CO2 purity vs experiment | experiment | 96.21 % | 97.2 % | 0.99 pp abs. error | 1.3 pp | Table 2, Exp. VI, Experimental Pu |
| ExperimentVI | CO2 recovery vs experiment | experiment | 89.83 % | 83.1 % | 6.73 pp abs. error | 7.8 pp | Table 2, Exp. VI, Experimental Re |
| ExperimentVI | CO2 productivity vs experiment | experiment | 0.2801 mol CO2/(m3 s) | 0.3 mol CO2/(m3 s) | 6.62 % rel. error | 7.6 % | Table 2, Exp. VI, Experimental Pr |
| ExperimentVI | CO2 purity vs paper model | paper model | 96.21 % | 96.7 % | 0.49 pp abs. error | 0.74 pp | Table 2, Exp. VI, Model Pu |
| ExperimentVI | CO2 recovery vs paper model | paper model | 89.83 % | 84.8 % | 5.03 pp abs. error | 5.9 pp | Table 2, Exp. VI, Model Re |
| ExperimentVI | CO2 productivity vs paper model | paper model | 0.2801 mol CO2/(m3 s) | 0.28 mol CO2/(m3 s) | 0.0509 % rel. error | 0.26 % | Table 2, Exp. VI, Model Pr |
| ExperimentVII | CO2 purity vs experiment | experiment | 93.46 % | 94.3 % | 0.843 pp abs. error | 1.1 pp | Table 2, Exp. VII.a, Experimental Pu |
| ExperimentVII | CO2 recovery vs experiment | experiment | 95.25 % | 88.4 % | 6.85 pp abs. error | 9 pp | Table 2, Exp. VII.a, Experimental Re |
| ExperimentVII | CO2 productivity vs experiment | experiment | 0.158 mol CO2/(m3 s) | 0.14 mol CO2/(m3 s) | 12.8 % rel. error | 16 % | Table 2, Exp. VII.a, Experimental Pr |
| ExperimentVII | CO2 purity vs paper model | paper model | 93.46 % | 95.5 % | 2.04 pp abs. error | 2.3 pp | Table 2, Exp. VII, Model Pu |
| ExperimentVII | CO2 recovery vs paper model | paper model | 95.25 % | 88.1 % | 7.15 pp abs. error | 9.3 pp | Table 2, Exp. VII, Model Re |
| ExperimentVII | CO2 productivity vs paper model | paper model | 0.158 mol CO2/(m3 s) | 0.14 mol CO2/(m3 s) | 12.8 % rel. error | 16 % | Table 2, Exp. VII, Model Pr |
| ExperimentVII | CO2 purity vs repeat experiment VII.b | experiment | 93.46 % | 93.8 % | 0.343 pp abs. error | info | Table 2, Exp. VII.b, Experimental Pu |
| ExperimentVII | CO2 recovery vs repeat experiment VII.b | experiment | 95.25 % | 87.5 % | 7.75 pp abs. error | info | Table 2, Exp. VII.b, Experimental Re |
| ExperimentVIII | CO2 purity vs experiment | experiment | 89.9 % | 92.7 % | 2.8 pp abs. error | 3.1 pp | Table 2, Exp. VIII, Experimental Pu |
| ExperimentVIII | CO2 recovery vs experiment | experiment | 96.14 % | 87.8 % | 8.34 pp abs. error | 11 pp | Table 2, Exp. VIII, Experimental Re |
| ExperimentVIII | CO2 productivity vs experiment | experiment | 0.1551 mol CO2/(m3 s) | 0.17 mol CO2/(m3 s) | 8.74 % rel. error | 11 % | Table 2, Exp. VIII, Experimental Pr |
| ExperimentVIII | CO2 purity vs paper model | paper model | 89.9 % | 93.6 % | 3.7 pp abs. error | 4.1 pp | Table 2, Exp. VIII, Model Pu |
| ExperimentVIII | CO2 recovery vs paper model | paper model | 96.14 % | 89.9 % | 6.24 pp abs. error | 8.3 pp | Table 2, Exp. VIII, Model Re |
| ExperimentVIII | CO2 productivity vs paper model | paper model | 0.1551 mol CO2/(m3 s) | 0.16 mol CO2/(m3 s) | 3.03 % rel. error | 4 % | Table 2, Exp. VIII, Model Pr |
| ExperimentX | CO2 purity vs experiment | experiment | 80.56 % | 80.6 % | 0.035 pp abs. error | 0.13 pp | Table 2, Exp. X, Experimental Pu |
| ExperimentX | CO2 recovery vs experiment | experiment | 99.43 % | 92.4 % | 7.03 pp abs. error | 9.7 pp | Table 2, Exp. X, Experimental Re |
| ExperimentX | CO2 productivity vs experiment | experiment | 0.1507 mol CO2/(m3 s) | 0.15 mol CO2/(m3 s) | 0.463 % rel. error | 2.9 % | Table 2, Exp. X, Experimental Pr |
| ExperimentX | CO2 purity vs paper model | paper model | 80.56 % | 80.6 % | 0.035 pp abs. error | 0.13 pp | Table 2, Exp. X, Model Pu |
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

## Estupiñan Perez et al. (2019) {#perez-2019}

L. Estupiñan Perez, P. Sarkar, A. Rajendran, Experimental validation of multi-objective optimization techniques for design of vacuum swing adsorption processes, Separation and Purification Technology 224 (2019) 553-563. [doi:10.1016/j.seppur.2019.05.039](https://doi.org/10.1016/j.seppur.2019.05.039)

*System:* CO2/N2 VSA on zeolite 13X (4-step and LPP). *Compared against:* experimental and model purity/recovery (Table 3).

**Findings.**

- The decks use the paper's model forms (WP10 audit of the paper and SI, 2026-10-31): the two-term exponential evacuation pressure of Sect. 4.3 (0.019 and 0.101 1/s; the split of the terms is not printed and an equal split is used, chosen on A), wall ends at ambient (SI Table S2), the Blake-Kozeny pressure drop (Ergun) and the capacity-factor LDF of SI Eq. S4. Four-step cycles: CO2 purity within 0.2-0.6 pp and recovery within 0.7-1.1 pp of experiments A and B (recovery was 8-10 pp high); C1 is 4.5 pp below its experiment and 12 pp below the repeat C2, which differ from each other by 7 pp. LPP cycles: purity 1.5-3.3 pp low and recovery 1.9-3.3 pp high.
- Attribution (README): the four-step recovery excess was the evacuation pressure history. The paper's profile keeps the column above P_low for most of the 57 s step where the deck's single exponential reached it; the capacity-factor kinetics, the wall ends and the pressure-drop constant are each worth under 0.5 pp. The LPP cycles' 167 s evacuation reaches P_low whatever the profile; their remaining recovery excess is the composition-coupled LPP connection.
- Open: the split of the two evacuation terms (Fig. 6c would pin it); the experiment ended steps on pressure triggers whose durations are not printed; the one-second hold between steps and about 5 % dead volume are not modelled.

| Case | Quantity | Against | Simulated | Reference | Error | Tolerance | Source |
|---|---|---|---|---|---|---|---|
| ConfigA | CO2 purity vs experiment | experiment | 99.51 % | 99.7 % | 0.192 pp abs. error | 0.3 pp | Table 3, 4-step, Experiment A, Experimental Purity |
| ConfigA | CO2 recovery vs experiment | experiment | 59.64 % | 58.7 % | 0.944 pp abs. error | 1.2 pp | Table 3, 4-step, Experiment A, Experimental Recovery |
| ConfigA | CO2 purity vs paper model | paper model | 99.51 % | 99.3 % | 0.208 pp abs. error | 0.3 pp | Table 3, 4-step, Experiment A, Model Purity |
| ConfigA | CO2 recovery vs paper model | paper model | 59.64 % | 60 % | 0.356 pp abs. error | 0.43 pp | Table 3, 4-step, Experiment A, Model Recovery |
| ConfigB | CO2 purity vs experiment | experiment | 96.79 % | 96.1 % | 0.693 pp abs. error | 0.7 pp | Table 3, 4-step, Experiment B, Experimental Purity |
| ConfigB | CO2 recovery vs experiment | experiment | 61.66 % | 61.8 % | 0.144 pp abs. error | 0.2 pp | Table 3, 4-step, Experiment B, Experimental Recovery |
| ConfigB | CO2 purity vs paper model | paper model | 96.79 % | 95.2 % | 1.59 pp abs. error | 1.8 pp | Table 3, 4-step, Experiment B, Model Purity |
| ConfigB | CO2 recovery vs paper model | paper model | 61.66 % | 64.1 % | 2.44 pp abs. error | 4.3 pp | Table 3, 4-step, Experiment B, Model Recovery |
| ConfigC1 | CO2 purity vs experiment | experiment | 87.95 % | 88.1 % | 0.147 pp abs. error | 0.2 pp | Table 3, 4-step, Experiment C1, Experimental Purity |
| ConfigC1 | CO2 recovery vs experiment | experiment | 62.23 % | 64.6 % | 2.37 pp abs. error | 2.9 pp | Table 3, 4-step, Experiment C1, Experimental Recovery |
| ConfigC1 | CO2 purity vs paper model | paper model | 87.95 % | 83.4 % | 4.55 pp abs. error | 4.8 pp | Table 3, 4-step, Experiment C1/C2, Model Purity |
| ConfigC1 | CO2 recovery vs paper model | paper model | 62.23 % | 68.2 % | 5.97 pp abs. error | 9.7 pp | Table 3, 4-step, Experiment C1/C2, Model Recovery |
| ConfigC1 | CO2 purity vs repeat experiment C2 | experiment | 87.95 % | 86.9 % | 1.05 pp abs. error | info | Table 3, 4-step, Experiment C2, Experimental Purity |
| ConfigC1 | CO2 recovery vs repeat experiment C2 | experiment | 62.23 % | 71.9 % | 9.67 pp abs. error | info | Table 3, 4-step, Experiment C2, Experimental Recovery |
| ConfigD1 | CO2 purity vs experiment | experiment | 93.79 % | 97.1 % | 3.31 pp abs. error | 3.9 pp | Table 3, LPP, Experiment D1, Experimental Purity |
| ConfigD1 | CO2 recovery vs experiment | experiment | 89.83 % | 86.1 % | 3.73 pp abs. error | 4.5 pp | Table 3, LPP, Experiment D1, Experimental Recovery |
| ConfigD1 | CO2 purity vs paper model | paper model | 93.79 % | 96.1 % | 2.31 pp abs. error | 2.7 pp | Table 3, LPP, Experiment D1/D2, Model Purity |
| ConfigD1 | CO2 recovery vs paper model | paper model | 89.83 % | 86.8 % | 3.03 pp abs. error | 3.7 pp | Table 3, LPP, Experiment D1/D2, Model Recovery |
| ConfigD1 | CO2 purity vs repeat experiment D2 | experiment | 93.79 % | 94.5 % | 0.706 pp abs. error | info | Table 3, LPP, Experiment D2, Experimental Purity |
| ConfigD1 | CO2 recovery vs repeat experiment D2 | experiment | 89.83 % | 83.4 % | 6.43 pp abs. error | info | Table 3, LPP, Experiment D2, Experimental Recovery |
| ConfigE | CO2 purity vs experiment | experiment | 86.65 % | 88.5 % | 1.85 pp abs. error | 2.1 pp | Table 3, LPP, Experiment E, Experimental Purity |
| ConfigE | CO2 recovery vs experiment | experiment | 94.36 % | 91.7 % | 2.66 pp abs. error | 3.2 pp | Table 3, LPP, Experiment E, Experimental Recovery |
| ConfigE | CO2 purity vs paper model | paper model | 86.65 % | 90.4 % | 3.75 pp abs. error | 4.3 pp | Table 3, LPP, Experiment E, Model Purity |
| ConfigE | CO2 recovery vs paper model | paper model | 94.36 % | 91.1 % | 3.26 pp abs. error | 4 pp | Table 3, LPP, Experiment E, Model Recovery |
| ConfigF | CO2 purity vs experiment | experiment | 75.65 % | 77.7 % | 2.05 pp abs. error | 2.6 pp | Table 3, LPP, Experiment F, Experimental Purity |
| ConfigF | CO2 recovery vs experiment | experiment | 97.69 % | 95.6 % | 2.09 pp abs. error | 2.6 pp | Table 3, LPP, Experiment F, Experimental Recovery |
| ConfigF | CO2 purity vs paper model | paper model | 75.65 % | 79.9 % | 4.25 pp abs. error | 4.4 pp | Table 3, LPP, Experiment F, Model Purity |
| ConfigF | CO2 recovery vs paper model | paper model | 97.69 % | 94.6 % | 3.09 pp abs. error | 3.8 pp | Table 3, LPP, Experiment F, Model Recovery |

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

*System:* CO2/CH4 and H2/CO2/CH4 breakthrough on zeolite 13X. *Compared against:* digitised experimental outlet compositions and thermocouples (Figs. 6, 7, 11).

**Findings.**

- CO2 front positions agree with the digitised experiments to 0.5-9 %, CH4 fronts to 0.3-15 %. Curve NRMSE is large (13-83 % for CH4) because the fronts are nearly vertical: a 10-20 s shift on a 20 s front dominates the RMSE. Thermocouple RMSE is 11-33 K, driven by the timing and height of the 30-100 K temperature fronts at the centre line.
- adsim uses extended Sips competition and the MFC set-point velocity; the paper shows IAST/RAST predict the binary cases better than extended Sips and fits the feed velocity per experiment within the MFC uncertainty. The ~10 m of sampling line to the MS is not modelled.
- A1/A3 previously compared simulated CO2 against the digitised CH4 curve (alphabetical file order) with a bound of 0.8 in mole fraction; fixed.

| Case | Quantity | Against | Simulated | Reference | Error | Tolerance | Source |
|---|---|---|---|---|---|---|---|
| ExperimentA1 | outlet y_CO2 | experiment | curve | digitised data | 44.7 % NRMSE | 53 % | Fig. 7 (digitised) |
| ExperimentA1 | outlet y_CH4 | experiment | curve | digitised data | 52.4 % NRMSE | 83 % | Fig. 7 (digitised) |
| ExperimentA1 | CO2 front time (y_CO2 = y_feed/2) | experiment | 3929 s | 3663 s | 7.28 % rel. error | 9.8 % | Fig. 7 (digitised) |
| ExperimentA1 | CH4 front time (y_CH4 = y_feed/2) | experiment | 475.7 s | 458.6 s | 3.74 % rel. error | 12 % | Fig. 7 (digitised) |
| ExperimentA1 | thermocouple T (5 positions, pooled) | experiment | curve | digitised data | 20.3 K RMSE | 25 K | Fig. 7 (digitised) |
| ExperimentA2 | outlet y_CO2 | experiment | curve | digitised data | 36 % NRMSE | 37 % | Fig. 7 (digitised) |
| ExperimentA2 | outlet y_CH4 | experiment | curve | digitised data | 84.1 % NRMSE | 96 % | Fig. 7 (digitised) |
| ExperimentA2 | CO2 front time (y_CO2 = y_feed/2) | experiment | 911.3 s | 965.7 s | 5.63 % rel. error | 5.8 % | Fig. 7 (digitised) |
| ExperimentA2 | CH4 front time (y_CH4 = y_feed/2) | experiment | 326.3 s | 365 s | 10.6 % rel. error | 13 % | Fig. 7 (digitised) |
| ExperimentA2 | thermocouple T (5 positions, pooled) | experiment | curve | digitised data | 19.9 K RMSE | 22 K | Fig. 7 (digitised) |
| ExperimentA3 | outlet y_CO2 | experiment | curve | digitised data | 17.8 % NRMSE | 27 % | Fig. 6 (digitised) |
| ExperimentA3 | outlet y_CH4 | experiment | curve | digitised data | 60.9 % NRMSE | 73 % | Fig. 6 (digitised) |
| ExperimentA3 | CO2 front time (y_CO2 = y_feed/2) | experiment | 433.3 s | 414.9 s | 4.44 % rel. error | 8 % | Fig. 6 (digitised) |
| ExperimentA3 | CH4 front time (y_CH4 = y_feed/2) | experiment | 151.3 s | 156.5 s | 3.33 % rel. error | 4.1 % | Fig. 6 (digitised) |
| ExperimentA3 | thermocouple T (5 positions, pooled) | experiment | curve | digitised data | 11.4 K RMSE | 19 K | Fig. 6 (digitised) |
| ExperimentA3 | CH4 breakthrough time (y &gt; 0.01) | experiment (text, approx.) | 149.7 s | 140 s | 9.75 s abs. error | info | Sect. 4.2.1 text (Fig. 6, reference experiment A3) |
| ExperimentA3 | CO2 breakthrough time (y &gt; 0.01) | experiment (text, approx.) | 419.6 s | 400 s | 19.6 s abs. error | info | Sect. 4.2.1 text (Fig. 6, reference experiment A3) |
| ExperimentA4 | outlet y_CO2 | experiment | curve | digitised data | 29.5 % NRMSE | 35 % | Fig. 7 (digitised) |
| ExperimentA4 | outlet y_CH4 | experiment | curve | digitised data | 79.4 % NRMSE | 91 % | Fig. 7 (digitised) |
| ExperimentA4 | CO2 front time (y_CO2 = y_feed/2) | experiment | 202.6 s | 226.2 s | 10.4 % rel. error | 10.5 % | Fig. 7 (digitised) |
| ExperimentA4 | CH4 front time (y_CH4 = y_feed/2) | experiment | 70.27 s | 85.66 s | 18 % rel. error | 19.5 % | Fig. 7 (digitised) |
| ExperimentA4 | thermocouple T (5 positions, pooled) | experiment | curve | digitised data | 16.4 K RMSE | 20 K | Fig. 7 (digitised) |
| ExperimentA5 | outlet y_CO2 | experiment | curve | digitised data | 18.5 % NRMSE | 19 % | Fig. 7 (digitised) |
| ExperimentA5 | outlet y_CH4 | experiment | curve | digitised data | 66.3 % NRMSE | 76 % | Fig. 7 (digitised) |
| ExperimentA5 | CO2 front time (y_CO2 = y_feed/2) | experiment | 177.2 s | 189.5 s | 6.48 % rel. error | 6.6 % | Fig. 7 (digitised) |
| ExperimentA5 | CH4 front time (y_CH4 = y_feed/2) | experiment | 80.82 s | 94.98 s | 14.9 % rel. error | 17 % | Fig. 7 (digitised) |
| ExperimentA5 | thermocouple T (5 positions, pooled) | experiment | curve | digitised data | 14.3 K RMSE | 21 K | Fig. 7 (digitised) |
| ExperimentA6 | outlet y_CO2 | experiment | curve | digitised data | 17 % NRMSE | 17.5 % | Fig. 7 (digitised) |
| ExperimentA6 | outlet y_CH4 | experiment | curve | digitised data | 75.2 % NRMSE | 85 % | Fig. 7 (digitised) |
| ExperimentA6 | CO2 front time (y_CO2 = y_feed/2) | experiment | 434.8 s | 447.9 s | 2.92 % rel. error | 3 % | Fig. 7 (digitised) |
| ExperimentA6 | CH4 front time (y_CH4 = y_feed/2) | experiment | 141.9 s | 160.6 s | 11.6 % rel. error | 14 % | Fig. 7 (digitised) |
| ExperimentA6 | thermocouple T (5 positions, pooled) | experiment | curve | digitised data | 11.6 K RMSE | 13 K | Fig. 7 (digitised) |
| ExperimentA7 | outlet y_CO2 | experiment | curve | digitised data | 2.44 % NRMSE | 14 % | Fig. 7 (digitised) |
| ExperimentA7 | outlet y_CH4 | experiment | curve | digitised data | 69.1 % NRMSE | 72 % | Fig. 7 (digitised) |
| ExperimentA7 | CO2 front time (y_CO2 = y_feed/2) | experiment | 427.8 s | 426.5 s | 0.31 % rel. error | 4 % | Fig. 7 (digitised) |
| ExperimentA7 | CH4 front time (y_CH4 = y_feed/2) | experiment | 131.7 s | 146.2 s | 9.94 % rel. error | 12 % | Fig. 7 (digitised) |
| ExperimentA7 | thermocouple T (5 positions, pooled) | experiment | curve | digitised data | 8.71 K RMSE | 14 K | Fig. 7 (digitised) |
| ExperimentB1 | outlet y_CO2 | experiment | curve | digitised data | 37.7 % NRMSE | 45 % | Fig. 11 (digitised) |
| ExperimentB1 | outlet y_CH4 | experiment | curve | digitised data | 12.7 % NRMSE | 15 % | Fig. 11 (digitised) |
| ExperimentB1 | CO2 front time (y_CO2 = y_feed/2) | experiment | 9019 s | 8474 s | 6.43 % rel. error | 8.1 % | Fig. 11 (digitised) |
| ExperimentB1 | CH4 front time (y_CH4 = y_feed/2) | experiment | 337.1 s | 337.7 s | 0.181 % rel. error | 0.4 % | Fig. 11 (digitised) |
| ExperimentB1 | outlet y_H2 | experiment | curve | digitised data | 37.2 % NRMSE | 42 % | Fig. 11 (digitised) |
| ExperimentB1 | thermocouple T (5 positions, pooled) | experiment | curve | digitised data | 17.8 K RMSE | 20 K | Fig. 11 (digitised) |
| ExperimentB2 | outlet y_CO2 | experiment | curve | digitised data | 45 % NRMSE | 57 % | Fig. 11 (digitised) |
| ExperimentB2 | outlet y_CH4 | experiment | curve | digitised data | 32.8 % NRMSE | 35 % | Fig. 11 (digitised) |
| ExperimentB2 | CO2 front time (y_CO2 = y_feed/2) | experiment | 2405 s | 2694 s | 10.7 % rel. error | 11 % | Fig. 11 (digitised) |
| ExperimentB2 | CH4 front time (y_CH4 = y_feed/2) | experiment | 630.8 s | 713.3 s | 11.6 % rel. error | 15 % | Fig. 11 (digitised) |
| ExperimentB2 | outlet y_H2 | experiment | curve | digitised data | 39.8 % NRMSE | 46 % | Fig. 11 (digitised) |
| ExperimentB2 | thermocouple T (5 positions, pooled) | experiment | curve | digitised data | 28.5 K RMSE | 38 K | Fig. 11 (digitised) |

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

- The decks use the paper's model forms (WP10 audit of the paper text, 2026-10-31): helium carrier, the paper's dispersion coefficient (Eq. 4 with D_m = 0.605e-5 m2/s, times the factor 4 the paper fitted), -dU_1 + RT as the heat of adsorption, the wall ends held at ambient (Table 4) and the Ergun pressure drop. The CO2 front now arrives 0.1-1.9 % late (4-6 % early before), the curve NRMSE is 3-18 % (34-38 % before), the thermocouple RMSE 1.2-4.0 K, and the temperature rise per Table 9 is 0.9-2.9 K from the measured value and inside the paper's 95 % interval at 100 and 200 ccm.
- Attribution (one change at a time, README): the fitted dispersion factor halves the temperature rise, and the wall end condition moves the front from early to slightly late; the carrier gas, the RT term and the pressure-drop constant are each worth under 1 % in the front time. The earlier explanation (a single site-1 heat of adsorption heating the bed more) was wrong: the paper's own model uses that value plus RT.
- Every deck value is traced to the paper in the README's provenance table. The deck keeps the particle voidage at 0 in the gas holdup (the paper's balances count the bed voidage only) and uses 60 cells to the paper's 30.

| Case | Quantity | Against | Simulated | Reference | Error | Tolerance | Source |
|---|---|---|---|---|---|---|---|
| Q100ccm | outlet y_CO2 | experiment | curve | digitised data | 26.5 % NRMSE | 29 % | Fig. 3 / data of ref. [5] (digitised) |
| Q100ccm | CO2 front time (y_CO2 = 0.5) | experiment | 1304 s | 1259 s | 3.56 % rel. error | 3.9 % | Fig. 3 / data of ref. [5] (digitised) |
| Q100ccm | thermocouple T at z = 0.0715 m | experiment | curve | digitised data | 1.15 K RMSE | 2.4 K | Fig. 3 / data of ref. [5] (digitised) |
| Q100ccm | breakthrough time (y_CO2,out = 0.5) | paper model (95% CI) | 1304 s | 1250 s | 54.3 s abs. error | info | Table 10, adsorption, f_in = 100 ccm, t_br 95% CI |
| Q100ccm | temperature rise at z = 0.0715 m | paper model (95% CI) | 13.65 K | 15.05 K | 1.4 K abs. error | info | Table 10, adsorption, f_in = 100 ccm, Delta T_ads 95% CI |
| Q200ccm | outlet y_CO2 | experiment | curve | digitised data | 14.1 % NRMSE | 17 % | Fig. 3 / data of ref. [5] (digitised) |
| Q200ccm | CO2 front time (y_CO2 = 0.5) | experiment | 601.3 s | 591.1 s | 1.74 % rel. error | 1.8 % | Fig. 3 / data of ref. [5] (digitised) |
| Q200ccm | thermocouple T at z = 0.0715 m | experiment | curve | digitised data | 2.56 K RMSE | 4.6 K | Fig. 3 / data of ref. [5] (digitised) |
| Q200ccm | breakthrough time (y_CO2,out = 0.5) | paper model (95% CI) | 601.3 s | 575 s | 26.3 s abs. error | info | Table 10, adsorption, f_in = 200 ccm, t_br 95% CI |
| Q200ccm | temperature rise at z = 0.0715 m | paper model (95% CI) | 23.41 K | 23.6 K | 0.195 K abs. error | info | Table 10, adsorption, f_in = 200 ccm, Delta T_ads 95% CI |
| Q50ccm | outlet y_CO2 | experiment | curve | digitised data | 35.2 % NRMSE | 39 % | Fig. 3 / data of ref. [5] (digitised) |
| Q50ccm | CO2 front time (y_CO2 = 0.5) | experiment | 2689 s | 2574 s | 4.48 % rel. error | 4.9 % | Fig. 3 / data of ref. [5] (digitised) |
| Q50ccm | thermocouple T at z = 0.0715 m | experiment | curve | digitised data | 1.24 K RMSE | 1.4 K | Fig. 3 / data of ref. [5] (digitised) |
| Q50ccm | breakthrough time (y_CO2,out = 0.5) | paper model (95% CI) | 2689 s | 2540 s | 149 s abs. error | info | Table 10, adsorption, f_in = 50 ccm, t_br 95% CI |
| Q50ccm | temperature rise at z = 0.0715 m | paper model (95% CI) | 5.851 K | 8.445 K | 2.59 K abs. error | info | Table 10, adsorption, f_in = 50 ccm, Delta T_ads 95% CI |

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
| FourStepPSA | CO2 purity | paper model | 96.11 % | 98.5 % | 2.39 pp abs. error | 2.7 pp | Table 4, Maximum separation effectiveness, Mg-MOF-74 (single-adsorbent) |
| FourStepPSA | CO2 recovery | paper model | 91.76 % | 89 % | 2.76 pp abs. error | 3.4 pp | Table 4, Maximum separation effectiveness, Mg-MOF-74 (single-adsorbent) |
| FourStepPSA | CO2 productivity | paper model | 0.8455 mol/(m3 s) | 0.828 mol/(m3 s) | 2.12 % rel. error | 2.7 % | Table 4, Maximum separation effectiveness, Mg-MOF-74 (single-adsorbent) |

??? note "Notes on these comparisons"

    - paper's optimised detailed-model result (simulation)
    - paper's detailed-model result
    - CO2 from the evacuation step per bed volume (V_bed = L A) per total cycle time (Sect. 4.2)

Input decks, assumptions and known discrepancies: [`core/test/papers/ward_pini_2024/README.md`](https://github.com/samaffleck/skarstrom/blob/main/core/test/papers/ward_pini_2024/README.md).

## Reproducing this page

```bash
cmake --preset linux-release
cmake --build --preset linux-release
ctest --preset linux-release -R "(Arora2021|Casas2012|Ferre2026|Nguyen2022|Nguyen2023|Perez2019|Streb2021|Vysyaraju2024|Ward2022|WardPini2024)" -j2
python scripts/docs/gen_validation_page.py --merge build/linux-release/validation
```

Each test writes `<build>/validation/<paper>__<case>.json`; the script merges them into `core/test/papers/validation_summary.json` and regenerates this page. The test suite also writes about 4 GB of results to `core/Results/`; delete it afterwards.
