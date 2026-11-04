# Cross-code regression

Single-column breakthrough and purge cases that were run in two other open-source adsorption
simulators, with their outputs frozen here as data. `test_crosscode_breakthrough.cpp` replays
each case in Skarstrom and checks the agreement. This turns a one-off code-to-code comparison
into a CI regression: if a change to Skarstrom moves a front away from the other code, the test
fails.

The cases, the runs in every tool, the grid studies and the attribution of each difference are
in the front-end repository, `benchmarks/open-source-psa/matrix/` of samaffleck/skarstrom-gui:
the harness and `REPORT.md`. This folder holds only what the test needs.

| File | What it is |
|---|---|
| `test_crosscode_breakthrough.cpp` | GoogleTest suites `CrossCodePeese` (BT01–BT10) and `CrossCodeToPSAil` (TP01, TP02) |
| `decks/<id>.json` | The Skarstrom deck of each case, written by the harness before the run (no results) |
| `data/<id>_<tool>_N<cells>.csv` | The reference tool's outlet curve: `t_s, y1_out` (CO₂ mole fraction), `T_out_K`; every output sample (2 s for PEESE, 5 s for toPSAil) where y or T moves, every 50 s elsewhere |
| `references.json` | Per case: the full case specification, the reference tool, version, licence and grid, its front metrics (and those of its finest grid, where it ran), the deck's Skarstrom options and the metrics the harness measured for Skarstrom on engine `main` |

## Reference tools and licences

| Tool | Version | Licence | What is committed |
|---|---|---|---|
| PEESEgroup/PSA column model (Yancy-Caballero et al., *Mol. Syst. Des. Eng.* 5 (2020) 1329, doi:10.1039/D0ME00060D), run through PSASimulator.py's `FuncAdsorption`, its line-for-line Python port, unchanged | PEESEgroup/PSA `3de0832`, PSASimulator.py `d77d0fa`; SciPy `solve_ivp` BDF, rtol = atol = 1e-6 | MIT | outlet curves |
| toPSAil (Georgia Tech), under GNU Octave 8.4 with syntax-only patches (`benchmarks/open-source-psa/topsail/` in skarstrom-gui) | `6f6e409`; `ode15s`, RelTol 1e-4, AbsTol 1e-6 | GPL-3.0 | numerical outputs only, no code (the outputs of a program are not covered by its licence) |

## Cases

All are one 1 m column at 1 bar: a bed in a uniform equilibrium state is fed at a constant molar
flux (interstitial 0.30 m/s at the feed conditions) from the bottom, with the top held at 1 bar.
They use the competitive dual-site Langmuir form of the PEESE code and its parameters (zeolite
13X, Yancy-Caballero et al. 2020), unless stated. Constant LDF, constant axial dispersion
D_L = 0.7 D_m + v0 r_p, axial conduction 0.09 W/(m K), Ergun.

| ID | Case | Reference (cells) | t50 reference [s] | Skarstrom − reference (t05 / t50 / t95 [s], T_ext [K]) |
|---|---|---|---|---|
| BT01 | 15 % CO₂ in N₂, isothermal 313 K | PEESE (100) | 3306.3 | −0.45 / −0.48 / −0.51, 0 |
| BT02 | 15 % CO₂, adiabatic | PEESE (100) | 1702.3 | +10.5 / +10.8 / –, −0.15 |
| BT03 | purge of a 15 % bed with N₂, isothermal | PEESE (100) | 1397.0 | −0.14 / +0.15 / –, 0 |
| BT04 | purge, adiabatic | PEESE (100) | 635.3 | −0.33 / −0.31 / –, −0.43 |
| BT05 | 50 % CO₂, single-site Langmuir, **pressure** form, isothermal | PEESE (100) | 1795.2 | −0.04 / −0.09 / −0.14, 0 |
| BT06 | 50 % CO₂, pressure-form Langmuir, adiabatic | PEESE (100) | 962.3 | +14.5 / +14.5 / +15.9, −0.87 |
| BT07 | 1 % CO₂, Henry isotherm (the 13X slope at 393 K), N₂ inert, isothermal 393 K | PEESE (100) | 1292.6 | +3.0 / +0.41 / −3.4, 0 |
| BT08 | 15 % CO₂, 0.3 mm particles (Ergun drop 25 % of P₀), isothermal | PEESE (85) | 3380.0 | +2.9 / +2.8 / +2.7, 0 |
| BT09 | 40 % CO₂ in CH₄ on LiX (Ferré et al. 2026 parameters, library entry `lix_ferre2026`), k = 46.7 1/s, isothermal | PEESE (100) | 978.6 | −0.24 / −0.37 / −0.59, 0 |
| BT10 | as BT09, adiabatic (heats chosen for the comparison: 30 and 15 kJ/mol) | PEESE (100) | 558.9 | +2.0 / +1.8 / +1.8, 0 |
| TP01 | BT01 with toPSAil's model choices | toPSAil (100) | 3304.4 | −0.51 / −0.51 / −0.55, 0 |
| TP02 | BT02 with a 5 mm steel wall, h_in = 20, h_out = 5 W/(m² K), ambient 313 K, toPSAil's choices | toPSAil (40) | 2031.4 | +7.2 / +6.4 / –, −0.25 |

"–": the front does not reach 95 % in the window (the adiabatic fronts' plateau, and the
purge's long tail).

## Why the codes differ where they do

Every difference was attributed in the harness by changing one modelling choice at a time
(REPORT.md has the tables):

- **Isothermal cases (BT01, BT03, BT05, BT07, BT09, TP01).** Agreement within 0.6 s on t50 at
  the same grid; t05 and t95 within 3 s (BT07). Both codes change their t05/t95 by more than
  that between 100 and 400 cells (e.g. BT01 t05 3274.5 → 3292.1 s in Skarstrom), so what is
  left is the difference between two discretisations at the same cell count (cell-centred
  finite volumes against the reference's nodes, two WENO variants), not a model difference.
- **Adiabatic cases (BT02, BT04, BT06, BT10).** Skarstrom takes the gas heat capacity from
  CoolProp (about 29.2 J/(mol K) for N₂ and 36 for CH₄ at 313 K) and the heat of adsorption as a
  constant −ΔU + R T₀; the reference uses 30.7 J/(mol K) and the local −(ΔU − R T). Switching
  Skarstrom to the reference's choices (the fluid switches of skarstrom#45) moves BT02 from
  +10.8 to −0.5 s, BT04's outlet minimum from −0.43 to −0.03 K, BT06 from +14.5 to +2.9 s
  (+1.0 s at 400 cells) and BT10 from +1.8 to +0.2 s. The tests assert the offsets
  measured without the switches, since those are not on `main` yet.
- **BT08 (pressure drop).** Skarstrom's viscosity polynomial gives 1.85e-5 Pa s against the
  reference's 1.72e-5. The larger pressure drop raises the mean bed pressure and the capacity:
  +2.8 s. With the reference's viscosity: −0.6 s.
- **toPSAil (TP01, TP02).** The decks take toPSAil's model: no axial dispersion or conduction
  (its CSTRs in series are first-order upwind cells), one heat per component for both
  isotherm sites and the heat release, the adsorbed phase heat capacity c_p,g q, and a
  negligible pressure drop. TP02's +6.4 s is again the CoolProp gas heat capacity (+0.1 s with
  30.7 J/(mol K)). toPSAil exchanges heat with the wall through the cross-sectional area
  π r² of each CSTR rather than its lateral area 2π r Δz (`getEnergyBalanceParams.m` lines
  121 and 154–163), so its wall heat transfer grows with the cell count: t50 is 1908, 2031 and
  2243 s at 20, 40 and 100 cells. At N cells this equals a physical coefficient h r N / (2 L),
  and the TP02 deck uses that (40 cells: ×1.0 inside, ×1.1 outside).

No Skarstrom defect was found on these cases.

## Bounds

From the harness run of each deck on engine `main` (the numbers in the table):

- where the codes agree: |Skarstrom − reference| ≤ 1.25 |measured| + floor;
- where the difference is attributed to a named property choice (BT02, BT04, BT06, BT08, BT10,
  TP02): Skarstrom − reference = measured ± (0.25 |measured| + floor);
- floors 0.25 s and 0.05 K; curve RMS ≤ 1.25 × measured + 0.002 (f = (y − y_init)/(y_feed −
  y_init) on a 1 s grid where either curve is between 0.01 and 0.99).

Guard: each t50 band is below 0.5 % of t50, so a 1 % error in the feed flux, bed capacity or
voidage (about 1 % on t50) fails the test.

When skarstrom#45 (constant gas heat capacity, constant viscosity, −ΔU + R T heat) is merged,
the decks of BT02, BT04, BT06, BT08, BT10 and TP02 can switch to the reference's property
choices and their bounds tighten to the "agree" form (measured with that branch: within 0.6 s
on t50 except BT06, +2.9 s).

Run time (Release build): about 1–3 s per case, 15 s for both suites.

```bash
ctest --preset linux-release -R "CrossCode"
```

To regenerate the data, run the harness in skarstrom-gui (`benchmarks/open-source-psa/matrix/`,
`python freeze_refs.py <engine>/core/test/crosscode`); it rewrites the decks, the curves and
`references.json`, and prints the measured differences from which the bounds are set.
