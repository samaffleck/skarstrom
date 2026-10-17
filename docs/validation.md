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

48 cases from 10 publications; 194 asserted comparisons (0 outside their recorded tolerance) and 15 informational ones. Summary generated on 2026-10-17.

| Publication | Cases | Against | Comparisons | Median error | Worst error |
|---|---|---|---|---|---|
| [Arora & Hasan (2021)](#arora-hasan-2021) | 2 | paper model | 4 | 4.09 pp (abs. error) | 45.6 pp (abs. error) |
| [Casas et al. (2012) / Schell et al. (2012)](#casas-2012) | 4 | paper model, isotherm fit | 8 | 3.68 % (rel. error)<br/>0.595 % (max rel. error) | 18.8 % (rel. error)<br/>6.99 % (max rel. error) |
| [Ferré et al. (2026)](#ferre-2026) | 4 | paper model | 20 | 1.11 pp (abs. error)<br/>7.37 % (rel. error) | 10.7 pp (abs. error)<br/>58.6 % (rel. error) |
| [Nguyen et al. (2022)](#nguyen-2022) | 6 | experiment | 12 | 9.94 % (NRMSE)<br/>1.69 % (rel. error) | 12.4 % (NRMSE)<br/>4.98 % (rel. error) |
| [Nguyen et al. (2023)](#nguyen-2023) | 10 | experiment, paper model | 60 | 2.07 pp (abs. error)<br/>7.89 % (rel. error) | 9.85 pp (abs. error)<br/>17.5 % (rel. error) |
| [Estupiñan Perez et al. (2019)](#perez-2019) | 6 | experiment, paper model | 24 | 3.42 pp (abs. error) | 9.63 pp (abs. error) |
| [Streb & Mazzotti (2021)](#streb-2021) | 9 | experiment | 47 | 36.1 % (NRMSE)<br/>7.27 % (rel. error)<br/>16.5 K (RMSE) | 82.9 % (NRMSE)<br/>14.8 % (rel. error)<br/>32.6 K (RMSE) |
| [Vysyaraju et al. (2024)](#vysyarajua-2024) | 3 | experiment, isotherm (derived), paper model | 7 | 1.81 % (rel. error)<br/>0.253 K (abs. error)<br/>1.91 pp (abs. error) | 3.81 % (rel. error)<br/>0.253 K (abs. error)<br/>14 pp (abs. error) |
| [Ward & Pini (2022)](#ward-2022) | 3 | experiment | 9 | 37.9 % (NRMSE)<br/>4.58 % (rel. error)<br/>3.54 K (RMSE) | 38.2 % (NRMSE)<br/>6.11 % (rel. error)<br/>4.21 K (RMSE) |
| [Ward & Pini (2024)](#ward-pini-2024) | 1 | paper model | 3 | 1.79 pp (abs. error)<br/>0.718 % (rel. error) | 2.3 pp (abs. error)<br/>0.718 % (rel. error) |

## Arora & Hasan (2021) {#arora-hasan-2021}

A. Arora, M.M.F. Hasan, Flexible oxygen concentrators for medical applications, Scientific Reports 11 (2021) 14317. [doi:10.1038/s41598-021-93796-3](https://doi.org/10.1038/s41598-021-93796-3)

*System:* O2/N2 (air) on LiLSX, 4-step PSA and PVSA. *Compared against:* paper's optimised model results (Table 4).

**Findings.**

- O2 purity is within about 4 percentage points of Table 4, but recovery is far below the paper (4 % vs 30 % for PSA, 20 % vs 66 % for PVSA). This is an open disagreement, not a tuned result.
- The paper's recovery is a model value from an optimiser; the t1/t4 rows of the PSA column of Table 4 were lost in the text extraction (t1 = 3.0 s is an estimate). The isotherm form is the non-competitive dual-site Langmuir; the paper does not state whether its DSL is competitive, so it was left as is.
- The PSA deck now vents the depressurisation step counter-currently (feed end), as the paper states for every Table 4 case; it previously vented co-currently.

| Case | Quantity | Against | Simulated | Reference | Error | Tolerance | Source |
|---|---|---|---|---|---|---|---|
| LiLSX_PSA | O2 purity at cycle 50 | paper model | 88.56 % | 84.72 % | 3.84 pp abs. error | 4.5 pp | Table 4, PSA cycle with LiLSX adsorbent, column PC_O2 &gt;= 5 L/min, P_O2 &gt;= 85% |
| LiLSX_PSA | O2 recovery at cycle 50 | paper model | 4.238 % | 30.14 % | 25.9 pp abs. error | 30 pp | Table 4, PSA cycle with LiLSX adsorbent, column PC_O2 &gt;= 5 L/min, P_O2 &gt;= 85% |
| LiLSX_PVSA | O2 purity at cycle 50 | paper model | 94.06 % | 89.97 % | 4.09 pp abs. error | 4.7 pp | Table 4, PVSA cycle with LiLSX adsorbent, column PC_O2 &gt;= 5 L/min, P_O2 &gt;= 90% |
| LiLSX_PVSA | O2 recovery at cycle 50 | paper model | 19.88 % | 65.51 % | 45.6 pp abs. error | 53 pp | Table 4, PVSA cycle with LiLSX adsorbent, column PC_O2 &gt;= 5 L/min, P_O2 &gt;= 90% |

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
| Fig4_Reference_50_50_15bar_25C | H2 breakthrough time (y &gt; 0.002) | paper model | 88.87 s | 105 s | 15.4 % rel. error | 18 % | Casas et al. 2012, Table 5 (50/50, 25 C, 15 bar) |
| Fig4_Reference_50_50_15bar_25C | CO2 breakthrough time (y &gt; 0.002) | paper model | 445.8 s | 440 s | 1.32 % rel. error | 2.4 % | Casas et al. 2012, Table 5 (50/50, 25 C, 15 bar) |
| Fig4_Reference_50_50_15bar_25C | H2 breakthrough time vs experiment (approx.) | experiment | 88.87 s | 110 s | 19.2 % rel. error | info | Casas et al. 2012, Sect. 4.1 text (Fig. 4) |
| Fig6_25_75_15bar_25C | H2 breakthrough time (y &gt; 0.002) | paper model | 64.97 s | 80 s | 18.8 % rel. error | 22 % | Casas et al. 2012, Table 5 (25/75, 25 C, 15 bar) |
| Fig6_25_75_15bar_25C | CO2 breakthrough time (y &gt; 0.002) | paper model | 613.3 s | 610 s | 0.543 % rel. error | 1.6 % | Casas et al. 2012, Table 5 (25/75, 25 C, 15 bar) |
| Fig6_75_25_15bar_25C | H2 breakthrough time (y &gt; 0.002) | paper model | 134.8 s | 130 s | 3.68 % rel. error | 4.7 % | Casas et al. 2012, Table 5 (75/25, 25 C, 15 bar) |
| Fig6_75_25_15bar_25C | CO2 breakthrough time (y &gt; 0.002) | paper model | 360.8 s | 320 s | 12.7 % rel. error | 15 % | Casas et al. 2012, Table 5 (75/25, 25 C, 15 bar) |
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
| MIL160Al | CO2 purity (heavy product) | paper model | 96.83 % | 95 % | 1.83 pp abs. error | 2.2 pp | Table 5, MIL-160(Al), recovery-recovery case, first row |
| MIL160Al | CO2 recovery | paper model | 91.11 % | 94.1 % | 2.99 pp abs. error | 3.5 pp | Table 5, MIL-160(Al), recovery-recovery case, first row |
| MIL160Al | CH4 purity (light product) | paper model | 96.5 % | 97.6 % | 1.1 pp abs. error | 1.4 pp | Table 5, MIL-160(Al), recovery-recovery case, first row |
| MIL160Al | CH4 recovery | paper model | 46.55 % | 35.8 % | 10.7 pp abs. error | 13 pp | Table 5, MIL-160(Al), recovery-recovery case, first row |
| MIL160Al | CH4 productivity | paper model | 1.427 mol/(m3 s) | 0.9 mol/(m3 s) | 58.6 % rel. error | 68 % | Table 5, MIL-160(Al), recovery-recovery case, first row |
| NaY_F | CO2 purity (heavy product) | paper model | 94.36 % | 95.1 % | 0.741 pp abs. error | 1.1 pp | Table 5, NaY (F), recovery-recovery (Re_CO2-Re_CH4) case, first row |
| NaY_F | CO2 recovery | paper model | 96.84 % | 97.8 % | 0.962 pp abs. error | 1.3 pp | Table 5, NaY (F), recovery-recovery case, first row |
| NaY_F | CH4 purity (light product) | paper model | 98.51 % | 98.8 % | 0.287 pp abs. error | 0.59 pp | Table 5, NaY (F), recovery-recovery case, first row |
| NaY_F | CH4 recovery | paper model | 61.18 % | 66 % | 4.82 pp abs. error | 5.6 pp | Table 5, NaY (F), recovery-recovery case, first row |
| NaY_F | CH4 productivity | paper model | 3.087 mol/(m3 s) | 3.2 mol/(m3 s) | 3.52 % rel. error | 4.6 % | Table 5, NaY (F), recovery-recovery case, first row |
| NaY_W | CO2 purity (heavy product) | paper model | 93.5 % | 95 % | 1.5 pp abs. error | 1.8 pp | Table 5, NaY (W), recovery-recovery case, first row |
| NaY_W | CO2 recovery | paper model | 95.32 % | 95.9 % | 0.58 pp abs. error | 0.89 pp | Table 5, NaY (W), recovery-recovery case, first row |
| NaY_W | CH4 purity (light product) | paper model | 97.74 % | 97.6 % | 0.142 pp abs. error | 0.45 pp | Table 5, NaY (W), recovery-recovery case, first row |
| NaY_W | CH4 recovery | paper model | 58.88 % | 65.4 % | 6.52 pp abs. error | 7.6 pp | Table 5, NaY (W), recovery-recovery case, first row |
| NaY_W | CH4 productivity | paper model | 2.872 mol/(m3 s) | 3.1 mol/(m3 s) | 7.37 % rel. error | 8.5 % | Table 5, NaY (W), recovery-recovery case, first row |
| Zeolite13X | CO2 purity (heavy product) | paper model | 93.89 % | 95 % | 1.11 pp abs. error | 1.5 pp | Table 5, 13X, recovery-recovery case, first row |
| Zeolite13X | CO2 recovery | paper model | 94.67 % | 95.9 % | 1.23 pp abs. error | 1.6 pp | Table 5, 13X, recovery-recovery case, first row |
| Zeolite13X | CH4 purity (light product) | paper model | 97.57 % | 98 % | 0.426 pp abs. error | 0.73 pp | Table 5, 13X, recovery-recovery case, first row |
| Zeolite13X | CH4 recovery | paper model | 43.73 % | 48.6 % | 4.87 pp abs. error | 5.7 pp | Table 5, 13X, recovery-recovery case, first row |
| Zeolite13X | CH4 productivity | paper model | 1.187 mol/(m3 s) | 1.3 mol/(m3 s) | 8.7 % rel. error | 11 % | Table 5, 13X, recovery-recovery case, first row |

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
| Fig4_BinaryBreakthrough | outlet y_CO2 | experiment | curve | digitised data | 8.62 % NRMSE | 10 % | Fig. 4 (digitised) |
| Fig4_BinaryBreakthrough | CO2 front time (y_CO2 = y_feed/2) | experiment | 1718 s | 1747 s | 1.69 % rel. error | 2.7 % | Fig. 4 (digitised) |
| Fig4_Binary_100pct | outlet y_CO2 | experiment | curve | digitised data | 9.94 % NRMSE | 12 % | Fig. 4 (digitised) |
| Fig4_Binary_100pct | CO2 front time (y_CO2 = y_feed/2) | experiment | 1250 s | 1258 s | 0.688 % rel. error | 1.7 % | Fig. 4 (digitised) |
| Fig4_Binary_50pct | outlet y_CO2 | experiment | curve | digitised data | 12.2 % NRMSE | 14 % | Fig. 4 (digitised) |
| Fig4_Binary_50pct | CO2 front time (y_CO2 = y_feed/2) | experiment | 2320 s | 2441 s | 4.98 % rel. error | 6 % | Fig. 4 (digitised) |
| Fig5_TernaryBreakthrough | outlet y_CO2 | experiment | curve | digitised data | 8.42 % NRMSE | 9.7 % | Fig. 5 (digitised) |
| Fig5_TernaryBreakthrough | CO2 front time (y_CO2 = y_feed/2) | experiment | 3200 s | 3364 s | 4.88 % rel. error | 5.9 % | Fig. 5 (digitised) |
| Fig5_Ternary_30pct | outlet y_CO2 | experiment | curve | digitised data | 10.4 % NRMSE | 13 % | Fig. 5 (digitised) |
| Fig5_Ternary_30pct | CO2 front time (y_CO2 = y_feed/2) | experiment | 1757 s | 1708 s | 2.83 % rel. error | 3.9 % | Fig. 5 (digitised) |
| Fig5_Ternary_5pct | outlet y_CO2 | experiment | curve | digitised data | 12.4 % NRMSE | 15 % | Fig. 5 (digitised) |
| Fig5_Ternary_5pct | CO2 front time (y_CO2 = y_feed/2) | experiment | 7282 s | 7397 s | 1.56 % rel. error | 2.6 % | Fig. 5 (digitised) |

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
- LPP cycles (VI-X): purity within 0.1-2.7 pp of the experiments, but recovery is 6-10 pp high and exceeds 100 % for IX and X. The single-bed emulation couples the light-product pressurisation to the feed step by composition, not by moles, so CO2 enters through the LPP stream that is not drawn from the raffinate; this is a limitation of the deck/connection model, recorded rather than hidden.
- The decks now use the competitive DSL form (as the paper's ref. [10]); the paper's column length (0.267 m) is inconsistent with its adsorbent mass and particle density, and the deck keeps the mass (L = 0.304 m).

| Case | Quantity | Against | Simulated | Reference | Error | Tolerance | Source |
|---|---|---|---|---|---|---|---|
| ExperimentI | CO2 purity vs experiment | experiment | 95.7 % | 97.7 % | 2 pp abs. error | 2.4 pp | Table 2, Exp. I, Experimental Pu |
| ExperimentI | CO2 recovery vs experiment | experiment | 69.33 % | 66.9 % | 2.43 pp abs. error | 2.8 pp | Table 2, Exp. I, Experimental Re |
| ExperimentI | CO2 productivity vs experiment | experiment | 0.6074 mol CO2/(m3 s) | 0.65 mol CO2/(m3 s) | 6.55 % rel. error | 7.6 % | Table 2, Exp. I, Experimental Pr |
| ExperimentI | CO2 purity vs paper model | paper model | 95.7 % | 96.7 % | 1 pp abs. error | 1.4 pp | Table 2, Exp. I, Model Pu |
| ExperimentI | CO2 recovery vs paper model | paper model | 69.33 % | 70.1 % | 0.77 pp abs. error | 1.1 pp | Table 2, Exp. I, Model Re |
| ExperimentI | CO2 productivity vs paper model | paper model | 0.6074 mol CO2/(m3 s) | 0.63 mol CO2/(m3 s) | 3.59 % rel. error | 4.6 % | Table 2, Exp. I, Model Pr |
| ExperimentII | CO2 purity vs experiment | experiment | 89.34 % | 91.2 % | 1.86 pp abs. error | 2.2 pp | Table 2, Exp. II, Experimental Pu |
| ExperimentII | CO2 recovery vs experiment | experiment | 72.93 % | 71.7 % | 1.23 pp abs. error | 1.6 pp | Table 2, Exp. II, Experimental Re |
| ExperimentII | CO2 productivity vs experiment | experiment | 0.4711 mol CO2/(m3 s) | 0.52 mol CO2/(m3 s) | 9.39 % rel. error | 11 % | Table 2, Exp. II, Experimental Pr |
| ExperimentII | CO2 purity vs paper model | paper model | 89.34 % | 90.3 % | 0.957 pp abs. error | 1.3 pp | Table 2, Exp. II, Model Pu |
| ExperimentII | CO2 recovery vs paper model | paper model | 72.93 % | 72.3 % | 0.634 pp abs. error | 0.94 pp | Table 2, Exp. II, Model Re |
| ExperimentII | CO2 productivity vs paper model | paper model | 0.4711 mol CO2/(m3 s) | 0.54 mol CO2/(m3 s) | 12.8 % rel. error | 15 % | Table 2, Exp. II, Model Pr |
| ExperimentIII | CO2 purity vs experiment | experiment | 83.9 % | 85.6 % | 1.7 pp abs. error | 2.1 pp | Table 2, Exp. III, Experimental Pu |
| ExperimentIII | CO2 recovery vs experiment | experiment | 74.29 % | 71.8 % | 2.49 pp abs. error | 2.9 pp | Table 2, Exp. III, Experimental Re |
| ExperimentIII | CO2 productivity vs experiment | experiment | 0.4393 mol CO2/(m3 s) | 0.48 mol CO2/(m3 s) | 8.47 % rel. error | 9.8 % | Table 2, Exp. III, Experimental Pr |
| ExperimentIII | CO2 purity vs paper model | paper model | 83.9 % | 85.7 % | 1.8 pp abs. error | 2.2 pp | Table 2, Exp. III, Model Pu |
| ExperimentIII | CO2 recovery vs paper model | paper model | 74.29 % | 73.4 % | 0.892 pp abs. error | 1.2 pp | Table 2, Exp. III, Model Re |
| ExperimentIII | CO2 productivity vs paper model | paper model | 0.4393 mol CO2/(m3 s) | 0.51 mol CO2/(m3 s) | 13.9 % rel. error | 16 % | Table 2, Exp. III, Model Pr |
| ExperimentIV | CO2 purity vs experiment | experiment | 79.12 % | 82.3 % | 3.18 pp abs. error | 3.7 pp | Table 2, Exp. IV, Experimental Pu |
| ExperimentIV | CO2 recovery vs experiment | experiment | 75.16 % | 72.5 % | 2.66 pp abs. error | 3.1 pp | Table 2, Exp. IV, Experimental Re |
| ExperimentIV | CO2 productivity vs experiment | experiment | 0.4277 mol CO2/(m3 s) | 0.46 mol CO2/(m3 s) | 7.02 % rel. error | 8.1 % | Table 2, Exp. IV, Experimental Pr |
| ExperimentIV | CO2 purity vs paper model | paper model | 79.12 % | 81.3 % | 2.18 pp abs. error | 2.6 pp | Table 2, Exp. IV, Model Pu |
| ExperimentIV | CO2 recovery vs paper model | paper model | 75.16 % | 74.3 % | 0.86 pp abs. error | 1.2 pp | Table 2, Exp. IV, Model Re |
| ExperimentIV | CO2 productivity vs paper model | paper model | 0.4277 mol CO2/(m3 s) | 0.5 mol CO2/(m3 s) | 14.5 % rel. error | 17 % | Table 2, Exp. IV, Model Pr |
| ExperimentIX | CO2 purity vs experiment | experiment | 84.69 % | 85.67 % | 0.982 pp abs. error | 1.3 pp | Table 2, Exp. IX, Experimental Pu |
| ExperimentIX | CO2 recovery vs experiment | experiment | 100.3 % | 90.5 % | 9.85 pp abs. error | 12 pp | Table 2, Exp. IX, Experimental Re |
| ExperimentIX | CO2 productivity vs experiment | experiment | 0.1452 mol CO2/(m3 s) | 0.13 mol CO2/(m3 s) | 11.7 % rel. error | 14 % | Table 2, Exp. IX, Experimental Pr |
| ExperimentIX | CO2 purity vs paper model | paper model | 84.69 % | 88.3 % | 3.61 pp abs. error | 4.2 pp | Table 2, Exp. IX, Model Pu |
| ExperimentIX | CO2 recovery vs paper model | paper model | 100.3 % | 91.6 % | 8.75 pp abs. error | 11 pp | Table 2, Exp. IX, Model Re |
| ExperimentIX | CO2 productivity vs paper model | paper model | 0.1452 mol CO2/(m3 s) | 0.14 mol CO2/(m3 s) | 3.69 % rel. error | 4.7 % | Table 2, Exp. IX, Model Pr |
| ExperimentV | CO2 purity vs experiment | experiment | 73.53 % | 78.4 % | 4.87 pp abs. error | 5.7 pp | Table 2, Exp. V, Experimental Pu |
| ExperimentV | CO2 recovery vs experiment | experiment | 76.78 % | 74.2 % | 2.58 pp abs. error | 3 pp | Table 2, Exp. V, Experimental Re |
| ExperimentV | CO2 productivity vs experiment | experiment | 0.3803 mol CO2/(m3 s) | 0.4 mol CO2/(m3 s) | 4.92 % rel. error | 6 % | Table 2, Exp. V, Experimental Pr |
| ExperimentV | CO2 purity vs paper model | paper model | 73.53 % | 75.6 % | 2.07 pp abs. error | 2.4 pp | Table 2, Exp. V, Model Pu |
| ExperimentV | CO2 recovery vs paper model | paper model | 76.78 % | 75.4 % | 1.38 pp abs. error | 1.7 pp | Table 2, Exp. V, Model Re |
| ExperimentV | CO2 productivity vs paper model | paper model | 0.3803 mol CO2/(m3 s) | 0.44 mol CO2/(m3 s) | 13.6 % rel. error | 16 % | Table 2, Exp. V, Model Pr |
| ExperimentVI | CO2 purity vs experiment | experiment | 96.27 % | 97.2 % | 0.935 pp abs. error | 1.3 pp | Table 2, Exp. VI, Experimental Pu |
| ExperimentVI | CO2 recovery vs experiment | experiment | 89.85 % | 83.1 % | 6.75 pp abs. error | 7.8 pp | Table 2, Exp. VI, Experimental Re |
| ExperimentVI | CO2 productivity vs experiment | experiment | 0.2802 mol CO2/(m3 s) | 0.3 mol CO2/(m3 s) | 6.59 % rel. error | 7.6 % | Table 2, Exp. VI, Experimental Pr |
| ExperimentVI | CO2 purity vs paper model | paper model | 96.27 % | 96.7 % | 0.435 pp abs. error | 0.74 pp | Table 2, Exp. VI, Model Pu |
| ExperimentVI | CO2 recovery vs paper model | paper model | 89.85 % | 84.8 % | 5.05 pp abs. error | 5.9 pp | Table 2, Exp. VI, Model Re |
| ExperimentVI | CO2 productivity vs paper model | paper model | 0.2802 mol CO2/(m3 s) | 0.28 mol CO2/(m3 s) | 0.0773 % rel. error | 1.1 % | Table 2, Exp. VI, Model Pr |
| ExperimentVII | CO2 purity vs experiment | experiment | 93.6 % | 94.3 % | 0.701 pp abs. error | 1.1 pp | Table 2, Exp. VII.a, Experimental Pu |
| ExperimentVII | CO2 recovery vs experiment | experiment | 96.15 % | 88.4 % | 7.75 pp abs. error | 9 pp | Table 2, Exp. VII.a, Experimental Re |
| ExperimentVII | CO2 productivity vs experiment | experiment | 0.1595 mol CO2/(m3 s) | 0.14 mol CO2/(m3 s) | 13.9 % rel. error | 16 % | Table 2, Exp. VII.a, Experimental Pr |
| ExperimentVII | CO2 purity vs paper model | paper model | 93.6 % | 95.5 % | 1.9 pp abs. error | 2.3 pp | Table 2, Exp. VII, Model Pu |
| ExperimentVII | CO2 recovery vs paper model | paper model | 96.15 % | 88.1 % | 8.05 pp abs. error | 9.3 pp | Table 2, Exp. VII, Model Re |
| ExperimentVII | CO2 productivity vs paper model | paper model | 0.1595 mol CO2/(m3 s) | 0.14 mol CO2/(m3 s) | 13.9 % rel. error | 16 % | Table 2, Exp. VII, Model Pr |
| ExperimentVII | CO2 purity vs repeat experiment VII.b | experiment | 93.6 % | 93.8 % | 0.201 pp abs. error | info | Table 2, Exp. VII.b, Experimental Pu |
| ExperimentVII | CO2 recovery vs repeat experiment VII.b | experiment | 96.15 % | 87.5 % | 8.65 pp abs. error | info | Table 2, Exp. VII.b, Experimental Re |
| ExperimentVIII | CO2 purity vs experiment | experiment | 90.08 % | 92.7 % | 2.62 pp abs. error | 3.1 pp | Table 2, Exp. VIII, Experimental Pu |
| ExperimentVIII | CO2 recovery vs experiment | experiment | 97.04 % | 87.8 % | 9.24 pp abs. error | 11 pp | Table 2, Exp. VIII, Experimental Re |
| ExperimentVIII | CO2 productivity vs experiment | experiment | 0.1566 mol CO2/(m3 s) | 0.17 mol CO2/(m3 s) | 7.89 % rel. error | 9.1 % | Table 2, Exp. VIII, Experimental Pr |
| ExperimentVIII | CO2 purity vs paper model | paper model | 90.08 % | 93.6 % | 3.52 pp abs. error | 4.1 pp | Table 2, Exp. VIII, Model Pu |
| ExperimentVIII | CO2 recovery vs paper model | paper model | 97.04 % | 89.9 % | 7.14 pp abs. error | 8.3 pp | Table 2, Exp. VIII, Model Re |
| ExperimentVIII | CO2 productivity vs paper model | paper model | 0.1566 mol CO2/(m3 s) | 0.16 mol CO2/(m3 s) | 2.13 % rel. error | 3.2 % | Table 2, Exp. VIII, Model Pr |
| ExperimentX | CO2 purity vs experiment | experiment | 80.84 % | 80.6 % | 0.237 pp abs. error | 0.54 pp | Table 2, Exp. X, Experimental Pu |
| ExperimentX | CO2 recovery vs experiment | experiment | 100.8 % | 92.4 % | 8.38 pp abs. error | 9.7 pp | Table 2, Exp. X, Experimental Re |
| ExperimentX | CO2 productivity vs experiment | experiment | 0.1527 mol CO2/(m3 s) | 0.15 mol CO2/(m3 s) | 1.82 % rel. error | 2.9 % | Table 2, Exp. X, Experimental Pr |
| ExperimentX | CO2 purity vs paper model | paper model | 80.84 % | 80.6 % | 0.237 pp abs. error | 0.54 pp | Table 2, Exp. X, Model Pu |
| ExperimentX | CO2 recovery vs paper model | paper model | 100.8 % | 93.6 % | 7.18 pp abs. error | 8.3 pp | Table 2, Exp. X, Model Re |
| ExperimentX | CO2 productivity vs paper model | paper model | 0.1527 mol CO2/(m3 s) | 0.13 mol CO2/(m3 s) | 17.5 % rel. error | 21 % | Table 2, Exp. X, Model Pr |

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

- Four-step cycles: CO2 purity within 0.1-1.2 pp of the experiments, but recovery 8-10 pp high. LPP cycles: purity 3-9 pp low and recovery 3-4 pp high.
- adsim's macropore LDF coefficient omits the c/q* capacity factor used in the paper (much faster kinetics), and its Carman-Kozeny constant is 180 against the paper's 150. The decks now use the competitive DSL (Sect. 4.1) and the paper's K_z = 0.0903 W/(m K).

| Case | Quantity | Against | Simulated | Reference | Error | Tolerance | Source |
|---|---|---|---|---|---|---|---|
| ConfigA | CO2 purity vs experiment | experiment | 98.94 % | 99.7 % | 0.761 pp abs. error | 1.1 pp | Table 3, 4-step, Experiment A, Experimental Purity |
| ConfigA | CO2 recovery vs experiment | experiment | 68.33 % | 58.7 % | 9.63 pp abs. error | 12 pp | Table 3, 4-step, Experiment A, Experimental Recovery |
| ConfigA | CO2 purity vs paper model | paper model | 98.94 % | 99.3 % | 0.361 pp abs. error | 0.67 pp | Table 3, 4-step, Experiment A, Model Purity |
| ConfigA | CO2 recovery vs paper model | paper model | 68.33 % | 60 % | 8.33 pp abs. error | 9.6 pp | Table 3, 4-step, Experiment A, Model Recovery |
| ConfigB | CO2 purity vs experiment | experiment | 96 % | 96.1 % | 0.103 pp abs. error | 0.41 pp | Table 3, 4-step, Experiment B, Experimental Purity |
| ConfigB | CO2 recovery vs experiment | experiment | 70.32 % | 61.8 % | 8.52 pp abs. error | 9.9 pp | Table 3, 4-step, Experiment B, Experimental Recovery |
| ConfigB | CO2 purity vs paper model | paper model | 96 % | 95.2 % | 0.797 pp abs. error | 1.1 pp | Table 3, 4-step, Experiment B, Model Purity |
| ConfigB | CO2 recovery vs paper model | paper model | 70.32 % | 64.1 % | 6.22 pp abs. error | 7.2 pp | Table 3, 4-step, Experiment B, Model Recovery |
| ConfigC1 | CO2 purity vs experiment | experiment | 86.88 % | 88.1 % | 1.22 pp abs. error | 1.6 pp | Table 3, 4-step, Experiment C1, Experimental Purity |
| ConfigC1 | CO2 recovery vs experiment | experiment | 72.78 % | 64.6 % | 8.18 pp abs. error | 9.5 pp | Table 3, 4-step, Experiment C1, Experimental Recovery |
| ConfigC1 | CO2 purity vs paper model | paper model | 86.88 % | 83.4 % | 3.48 pp abs. error | 4.1 pp | Table 3, 4-step, Experiment C1/C2, Model Purity |
| ConfigC1 | CO2 recovery vs paper model | paper model | 72.78 % | 68.2 % | 4.58 pp abs. error | 5.3 pp | Table 3, 4-step, Experiment C1/C2, Model Recovery |
| ConfigC1 | CO2 purity vs repeat experiment C2 | experiment | 86.88 % | 86.9 % | 0.0159 pp abs. error | info | Table 3, 4-step, Experiment C2, Experimental Purity |
| ConfigC1 | CO2 recovery vs repeat experiment C2 | experiment | 72.78 % | 71.9 % | 0.88 pp abs. error | info | Table 3, 4-step, Experiment C2, Experimental Recovery |
| ConfigD1 | CO2 purity vs experiment | experiment | 93.72 % | 97.1 % | 3.38 pp abs. error | 3.9 pp | Table 3, LPP, Experiment D1, Experimental Purity |
| ConfigD1 | CO2 recovery vs experiment | experiment | 90.22 % | 86.1 % | 4.12 pp abs. error | 4.8 pp | Table 3, LPP, Experiment D1, Experimental Recovery |
| ConfigD1 | CO2 purity vs paper model | paper model | 93.72 % | 96.1 % | 2.38 pp abs. error | 2.8 pp | Table 3, LPP, Experiment D1/D2, Model Purity |
| ConfigD1 | CO2 recovery vs paper model | paper model | 90.22 % | 86.8 % | 3.42 pp abs. error | 4 pp | Table 3, LPP, Experiment D1/D2, Model Recovery |
| ConfigD1 | CO2 purity vs repeat experiment D2 | experiment | 93.72 % | 94.5 % | 0.779 pp abs. error | info | Table 3, LPP, Experiment D2, Experimental Purity |
| ConfigD1 | CO2 recovery vs repeat experiment D2 | experiment | 90.22 % | 83.4 % | 6.82 pp abs. error | info | Table 3, LPP, Experiment D2, Experimental Recovery |
| ConfigE | CO2 purity vs experiment | experiment | 86.56 % | 88.5 % | 1.94 pp abs. error | 2.3 pp | Table 3, LPP, Experiment E, Experimental Purity |
| ConfigE | CO2 recovery vs experiment | experiment | 94.87 % | 91.7 % | 3.17 pp abs. error | 3.7 pp | Table 3, LPP, Experiment E, Experimental Recovery |
| ConfigE | CO2 purity vs paper model | paper model | 86.56 % | 90.4 % | 3.84 pp abs. error | 4.5 pp | Table 3, LPP, Experiment E, Model Purity |
| ConfigE | CO2 recovery vs paper model | paper model | 94.87 % | 91.1 % | 3.77 pp abs. error | 4.4 pp | Table 3, LPP, Experiment E, Model Recovery |
| ConfigF | CO2 purity vs experiment | experiment | 75.41 % | 77.7 % | 2.29 pp abs. error | 2.7 pp | Table 3, LPP, Experiment F, Experimental Purity |
| ConfigF | CO2 recovery vs experiment | experiment | 98.33 % | 95.6 % | 2.73 pp abs. error | 3.2 pp | Table 3, LPP, Experiment F, Experimental Recovery |
| ConfigF | CO2 purity vs paper model | paper model | 75.41 % | 79.9 % | 4.49 pp abs. error | 5.2 pp | Table 3, LPP, Experiment F, Model Purity |
| ConfigF | CO2 recovery vs paper model | paper model | 98.33 % | 94.6 % | 3.73 pp abs. error | 4.3 pp | Table 3, LPP, Experiment F, Model Recovery |

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
| ExperimentA1 | outlet y_CO2 | experiment | curve | digitised data | 45.7 % NRMSE | 53 % | Fig. 7 (digitised) |
| ExperimentA1 | outlet y_CH4 | experiment | curve | digitised data | 72.1 % NRMSE | 83 % | Fig. 7 (digitised) |
| ExperimentA1 | CO2 front time (y_CO2 = y_feed/2) | experiment | 3972 s | 3663 s | 8.44 % rel. error | 9.8 % | Fig. 7 (digitised) |
| ExperimentA1 | CH4 front time (y_CH4 = y_feed/2) | experiment | 503.1 s | 458.6 s | 9.71 % rel. error | 12 % | Fig. 7 (digitised) |
| ExperimentA1 | thermocouple T (5 positions, pooled) | experiment | curve | digitised data | 21.1 K RMSE | 25 K | Fig. 7 (digitised) |
| ExperimentA2 | outlet y_CO2 | experiment | curve | digitised data | 30.2 % NRMSE | 35 % | Fig. 7 (digitised) |
| ExperimentA2 | outlet y_CH4 | experiment | curve | digitised data | 82.9 % NRMSE | 96 % | Fig. 7 (digitised) |
| ExperimentA2 | CO2 front time (y_CO2 = y_feed/2) | experiment | 927 s | 965.7 s | 4 % rel. error | 5.1 % | Fig. 7 (digitised) |
| ExperimentA2 | CH4 front time (y_CH4 = y_feed/2) | experiment | 337.1 s | 365 s | 7.64 % rel. error | 8.8 % | Fig. 7 (digitised) |
| ExperimentA2 | thermocouple T (5 positions, pooled) | experiment | curve | digitised data | 14.1 K RMSE | 17 K | Fig. 7 (digitised) |
| ExperimentA3 | outlet y_CO2 | experiment | curve | digitised data | 22.9 % NRMSE | 27 % | Fig. 6 (digitised) |
| ExperimentA3 | outlet y_CH4 | experiment | curve | digitised data | 21.2 % NRMSE | 25 % | Fig. 6 (digitised) |
| ExperimentA3 | CO2 front time (y_CO2 = y_feed/2) | experiment | 443.7 s | 414.9 s | 6.94 % rel. error | 8 % | Fig. 6 (digitised) |
| ExperimentA3 | CH4 front time (y_CH4 = y_feed/2) | experiment | 156.9 s | 156.5 s | 0.268 % rel. error | 1.3 % | Fig. 6 (digitised) |
| ExperimentA3 | thermocouple T (5 positions, pooled) | experiment | curve | digitised data | 15.8 K RMSE | 19 K | Fig. 6 (digitised) |
| ExperimentA3 | CH4 breakthrough time (y &gt; 0.01) | experiment (text, approx.) | 152.4 s | 140 s | 12.4 s abs. error | info | Sect. 4.2.1 text (Fig. 6, reference experiment A3) |
| ExperimentA3 | CO2 breakthrough time (y &gt; 0.01) | experiment (text, approx.) | 425.1 s | 400 s | 25.1 s abs. error | info | Sect. 4.2.1 text (Fig. 6, reference experiment A3) |
| ExperimentA4 | outlet y_CO2 | experiment | curve | digitised data | 25.1 % NRMSE | 29 % | Fig. 7 (digitised) |
| ExperimentA4 | outlet y_CH4 | experiment | curve | digitised data | 78.6 % NRMSE | 91 % | Fig. 7 (digitised) |
| ExperimentA4 | CO2 front time (y_CO2 = y_feed/2) | experiment | 208.9 s | 226.2 s | 7.65 % rel. error | 8.8 % | Fig. 7 (digitised) |
| ExperimentA4 | CH4 front time (y_CH4 = y_feed/2) | experiment | 72.95 s | 85.66 s | 14.8 % rel. error | 18 % | Fig. 7 (digitised) |
| ExperimentA4 | thermocouple T (5 positions, pooled) | experiment | curve | digitised data | 16.5 K RMSE | 20 K | Fig. 7 (digitised) |
| ExperimentA5 | outlet y_CO2 | experiment | curve | digitised data | 13.1 % NRMSE | 16 % | Fig. 7 (digitised) |
| ExperimentA5 | outlet y_CH4 | experiment | curve | digitised data | 65.4 % NRMSE | 76 % | Fig. 7 (digitised) |
| ExperimentA5 | CO2 front time (y_CO2 = y_feed/2) | experiment | 181.5 s | 189.5 s | 4.22 % rel. error | 5.3 % | Fig. 7 (digitised) |
| ExperimentA5 | CH4 front time (y_CH4 = y_feed/2) | experiment | 83.48 s | 94.98 s | 12.1 % rel. error | 14 % | Fig. 7 (digitised) |
| ExperimentA5 | thermocouple T (5 positions, pooled) | experiment | curve | digitised data | 17.5 K RMSE | 21 K | Fig. 7 (digitised) |
| ExperimentA6 | outlet y_CO2 | experiment | curve | digitised data | 4.77 % NRMSE | 5.8 % | Fig. 7 (digitised) |
| ExperimentA6 | outlet y_CH4 | experiment | curve | digitised data | 73.9 % NRMSE | 85 % | Fig. 7 (digitised) |
| ExperimentA6 | CO2 front time (y_CO2 = y_feed/2) | experiment | 445.8 s | 447.9 s | 0.463 % rel. error | 1.5 % | Fig. 7 (digitised) |
| ExperimentA6 | CH4 front time (y_CH4 = y_feed/2) | experiment | 146.6 s | 160.6 s | 8.68 % rel. error | 10 % | Fig. 7 (digitised) |
| ExperimentA6 | thermocouple T (5 positions, pooled) | experiment | curve | digitised data | 11.2 K RMSE | 13 K | Fig. 7 (digitised) |
| ExperimentA7 | outlet y_CO2 | experiment | curve | digitised data | 11.7 % NRMSE | 14 % | Fig. 7 (digitised) |
| ExperimentA7 | outlet y_CH4 | experiment | curve | digitised data | 62.3 % NRMSE | 72 % | Fig. 7 (digitised) |
| ExperimentA7 | CO2 front time (y_CO2 = y_feed/2) | experiment | 438.9 s | 426.5 s | 2.92 % rel. error | 4 % | Fig. 7 (digitised) |
| ExperimentA7 | CH4 front time (y_CH4 = y_feed/2) | experiment | 135.6 s | 146.2 s | 7.27 % rel. error | 8.4 % | Fig. 7 (digitised) |
| ExperimentA7 | thermocouple T (5 positions, pooled) | experiment | curve | digitised data | 11.6 K RMSE | 14 K | Fig. 7 (digitised) |
| ExperimentB1 | outlet y_CO2 | experiment | curve | digitised data | 38.8 % NRMSE | 45 % | Fig. 11 (digitised) |
| ExperimentB1 | outlet y_CH4 | experiment | curve | digitised data | 12.9 % NRMSE | 15 % | Fig. 11 (digitised) |
| ExperimentB1 | CO2 front time (y_CO2 = y_feed/2) | experiment | 9068 s | 8474 s | 7.01 % rel. error | 8.1 % | Fig. 11 (digitised) |
| ExperimentB1 | CH4 front time (y_CH4 = y_feed/2) | experiment | 335.5 s | 337.7 s | 0.657 % rel. error | 1.7 % | Fig. 11 (digitised) |
| ExperimentB1 | outlet y_H2 | experiment | curve | digitised data | 36.1 % NRMSE | 42 % | Fig. 11 (digitised) |
| ExperimentB1 | thermocouple T (5 positions, pooled) | experiment | curve | digitised data | 16.6 K RMSE | 20 K | Fig. 11 (digitised) |
| ExperimentB2 | outlet y_CO2 | experiment | curve | digitised data | 48.9 % NRMSE | 57 % | Fig. 11 (digitised) |
| ExperimentB2 | outlet y_CH4 | experiment | curve | digitised data | 29.7 % NRMSE | 35 % | Fig. 11 (digitised) |
| ExperimentB2 | CO2 front time (y_CO2 = y_feed/2) | experiment | 2448 s | 2694 s | 9.13 % rel. error | 11 % | Fig. 11 (digitised) |
| ExperimentB2 | CH4 front time (y_CH4 = y_feed/2) | experiment | 626.2 s | 713.3 s | 12.2 % rel. error | 15 % | Fig. 11 (digitised) |
| ExperimentB2 | outlet y_H2 | experiment | curve | digitised data | 39.9 % NRMSE | 46 % | Fig. 11 (digitised) |
| ExperimentB2 | thermocouple T (5 positions, pooled) | experiment | curve | digitised data | 32.6 K RMSE | 38 K | Fig. 11 (digitised) |

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
- He purity matches the paper model to 0.01 pp (1.9 pp above the experiment), but He recovery is 17 % against 31 %: the deck is a single bed with self-equalisation, the experiment and the paper's model are a 3-bed cycle.

| Case | Quantity | Against | Simulated | Reference | Error | Tolerance | Source |
|---|---|---|---|---|---|---|---|
| DCB_Ambient_N2 | N2 equilibrium loading | experiment | 0.5386 mol/kg | 0.56 mol/kg | 3.81 % rel. error | 4.9 % | Sect. 4.2 text (Fig. 4a) |
| DCB_Ambient_N2 | peak bed temperature | experiment | 304.3 K | 304 K | 0.253 K abs. error | 0.56 K | Sect. 4.2 text (Fig. 4) |
| DCB_HighPressure_N2 | N2 equilibrium loading | isotherm (derived) | 2.128 mol/kg | 2.167 mol/kg | 1.81 % rel. error | 2.9 % | Table 2 DSL parameters evaluated at 11.2 bar, 298 K (Sect. 4.2, Fig. 5b) |
| PSA_Exp1 | He purity vs experiment | experiment | 99.91 % | 98 % | 1.91 pp abs. error | 2.3 pp | Table 3, Exp. 1, Exp Pu_He |
| PSA_Exp1 | He purity vs paper model | paper model | 99.91 % | 99.9 % | 0.0124 pp abs. error | 0.32 pp | Table 3, Exp. 1, Model Pu_He |
| PSA_Exp1 | He recovery vs experiment | experiment | 17.32 % | 30.6 % | 13.3 pp abs. error | 16 pp | Table 3, Exp. 1, Exp Re_He |
| PSA_Exp1 | He recovery vs paper model | paper model | 17.32 % | 31.3 % | 14 pp abs. error | 17 pp | Table 3, Exp. 1, Model Re_He |

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

- The CO2 front (y = 0.5) arrives 4-6 % earlier than measured. The simulated front is much sharper than the measured one, so the curve NRMSE is 34-38 %; the thermocouple RMSE is 2.6-4.2 K, and the simulated temperature rise (15-33 K) is above the paper model's 95 % interval (7-27 K).
- No transcription error was found (column length, radius, voidage, bed density, flows and isotherm match Tables 1-2). Likely causes: a single constant heat of adsorption of 29.98 kJ/mol (the site-1 internal energy, while most of the loading sits on site 2 at 20.4 kJ/mol) heats the bed more and lowers its capacity; the deck uses air as carrier where the paper uses He; the sampling line is not modelled. The ambient temperature was corrected to Table 1's 291.4 K.

| Case | Quantity | Against | Simulated | Reference | Error | Tolerance | Source |
|---|---|---|---|---|---|---|---|
| Q100ccm | outlet y_CO2 | experiment | curve | digitised data | 37.9 % NRMSE | 44 % | Fig. 3 / data of ref. [5] (digitised) |
| Q100ccm | CO2 front time (y_CO2 = 0.5) | experiment | 1182 s | 1259 s | 6.11 % rel. error | 7.2 % | Fig. 3 / data of ref. [5] (digitised) |
| Q100ccm | thermocouple T at z = 0.0715 m | experiment | curve | digitised data | 3.54 K RMSE | 4.1 K | Fig. 3 / data of ref. [5] (digitised) |
| Q100ccm | breakthrough time (y_CO2,out = 0.5) | paper model (95% CI) | 1182 s | 1250 s | 67.5 s abs. error | info | Table 10, adsorption, f_in = 100 ccm, t_br 95% CI |
| Q100ccm | temperature rise at z = 0.0715 m | paper model (95% CI) | 25.26 K | 15.05 K | 10.2 K abs. error | info | Table 10, adsorption, f_in = 100 ccm, Delta T_ads 95% CI |
| Q200ccm | outlet y_CO2 | experiment | curve | digitised data | 33.6 % NRMSE | 39 % | Fig. 3 / data of ref. [5] (digitised) |
| Q200ccm | CO2 front time (y_CO2 = 0.5) | experiment | 564 s | 591.1 s | 4.58 % rel. error | 5.6 % | Fig. 3 / data of ref. [5] (digitised) |
| Q200ccm | thermocouple T at z = 0.0715 m | experiment | curve | digitised data | 4.21 K RMSE | 4.9 K | Fig. 3 / data of ref. [5] (digitised) |
| Q200ccm | breakthrough time (y_CO2,out = 0.5) | paper model (95% CI) | 564 s | 575 s | 11 s abs. error | info | Table 10, adsorption, f_in = 200 ccm, t_br 95% CI |
| Q200ccm | temperature rise at z = 0.0715 m | paper model (95% CI) | 32.53 K | 23.6 K | 8.93 K abs. error | info | Table 10, adsorption, f_in = 200 ccm, Delta T_ads 95% CI |
| Q50ccm | outlet y_CO2 | experiment | curve | digitised data | 38.2 % NRMSE | 44 % | Fig. 3 / data of ref. [5] (digitised) |
| Q50ccm | CO2 front time (y_CO2 = 0.5) | experiment | 2464 s | 2574 s | 4.27 % rel. error | 5.3 % | Fig. 3 / data of ref. [5] (digitised) |
| Q50ccm | thermocouple T at z = 0.0715 m | experiment | curve | digitised data | 2.59 K RMSE | 3.1 K | Fig. 3 / data of ref. [5] (digitised) |
| Q50ccm | breakthrough time (y_CO2,out = 0.5) | paper model (95% CI) | 2464 s | 2540 s | 76.9 s abs. error | info | Table 10, adsorption, f_in = 50 ccm, t_br 95% CI |
| Q50ccm | temperature rise at z = 0.0715 m | paper model (95% CI) | 14.9 K | 8.445 K | 6.46 K abs. error | info | Table 10, adsorption, f_in = 50 ccm, Delta T_ads 95% CI |

??? note "Notes on these comparisons"

    - RMSE over 29 digitised points; normalised by 1 -
    - first time the curve reaches 0.5 (reference time interpolated from the digitised points)
    - RMSE over 37 digitised points
    - 95% credibility interval of the paper's model prediction; reference = CI midpoint; simulated value inside the paper's 95% CI [1048, 1451]; t_br read as the time y_CO2,out reaches 0.5 (Table 9 only partly legible)
    - 95% credibility interval of the paper's model prediction; reference = CI midpoint; simulated value outside the paper's 95% CI [12.9, 17.2]; simulated = max T at the thermocouple minus initial T
    - RMSE over 28 digitised points; normalised by 1 -
    - RMSE over 39 digitised points
    - 95% credibility interval of the paper's model prediction; reference = CI midpoint; simulated value inside the paper's 95% CI [482, 668]; t_br read as the time y_CO2,out reaches 0.5 (Table 9 only partly legible)
    - 95% credibility interval of the paper's model prediction; reference = CI midpoint; simulated value outside the paper's 95% CI [20, 27.2]; simulated = max T at the thermocouple minus initial T
    - RMSE over 27 digitised points; normalised by 1 -
    - RMSE over 31 digitised points
    - 95% credibility interval of the paper's MODEL prediction (not a measurement); t_br definition (Table 9) is only partly legible in the extraction; reference = CI midpoint; simulated value inside the paper's 95% CI [2132, 2949]; t_br read as the time y_CO2,out reaches 0.5 (Table 9 only partly legible)
    - 95% credibility interval of the paper's model prediction; reference = CI midpoint; simulated value outside the paper's 95% CI [7.36, 9.53]; simulated = max T at the thermocouple minus initial T

Input decks, assumptions and known discrepancies: [`core/test/papers/ward_2022/README.md`](https://github.com/samaffleck/skarstrom/blob/main/core/test/papers/ward_2022/README.md).

## Ward & Pini (2024) {#ward-pini-2024}

A. Ward, R. Pini, Design and Performance Evaluation of Multisorbent Vacuum-Swing Adsorption Processes for Postcombustion Carbon Capture, Industrial & Engineering Chemistry Research 63 (2024) 13787-13800. [doi:10.1021/acs.iecr.4c01345](https://doi.org/10.1021/acs.iecr.4c01345)

*System:* CO2/N2 4-step VSA on Mg-MOF-74. *Compared against:* paper's optimised detailed-model results (Table 4).

**Findings.**

- CO2 purity within 2.3 pp, recovery within 1.8 pp and productivity within 0.7 % of Table 4.
- The deck now uses the extended (competitive) DSL the paper states (Sect. 2.1). Energy is not compared (no pump model in the deck).

| Case | Quantity | Against | Simulated | Reference | Error | Tolerance | Source |
|---|---|---|---|---|---|---|---|
| FourStepPSA | CO2 purity | paper model | 96.2 % | 98.5 % | 2.3 pp abs. error | 2.7 pp | Table 4, Maximum separation effectiveness, Mg-MOF-74 (single-adsorbent) |
| FourStepPSA | CO2 recovery | paper model | 90.79 % | 89 % | 1.79 pp abs. error | 2.1 pp | Table 4, Maximum separation effectiveness, Mg-MOF-74 (single-adsorbent) |
| FourStepPSA | CO2 productivity | paper model | 0.8339 mol/(m3 s) | 0.828 mol/(m3 s) | 0.718 % rel. error | 1.8 % | Table 4, Maximum separation effectiveness, Mg-MOF-74 (single-adsorbent) |

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
