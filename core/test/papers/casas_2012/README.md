# Casas / Schell (2012) — CO₂/H₂/N₂ on activated carbon (AP3-60)

Two companion papers from the same ETH Zurich group: the equilibrium study
supplies the isotherms, the fixed-bed study supplies the breakthrough
experiments.

## Citation

> Schell, J., Casas, N., Pini, R., Mazzotti, M. *Pure and binary adsorption of
> CO₂, H₂, and N₂ on activated carbon.* Adsorption 18(1):49–65 (2012).
> DOI: [10.1007/s10450-011-9382-y](https://doi.org/10.1007/s10450-011-9382-y)

> Casas, N., Schell, J., Pini, R., Mazzotti, M. *Fixed bed adsorption of CO₂/H₂
> mixtures on activated carbon: experiments and modeling.* Adsorption
> 18(2):143–161 (2012).
> DOI: [10.1007/s10450-012-9389-z](https://doi.org/10.1007/s10450-012-9389-z)

## Files

| File | Purpose |
|------|---------|
| `casas_2012.h`   | `Casas2012Test` fixture: shared physical setup + `RunBreakthroughCase(...)` |
| `casas_2012.cpp` | GoogleTest cases (one `TEST_F` per reproduced experiment) |

## Reproduced cases

Binary CO₂/H₂ breakthrough on AP3-60 activated carbon. The bed is pre-pressurised
with inert helium, then fed a CO₂/H₂ mixture at a constant superficial flow of
10 cm³/s (paper's setpoint). All three cases pass with mass, mole and
component-mole conservation.

| `TEST_F` | Paper | Feed CO₂/H₂ | P | T |
|----------|-------|-------------|---|---|
| `Fig4_Reference_50_50_15bar_25C` | Fig. 4 (reference / fitting experiment) | 50/50 | 15 bar | 25 °C |
| `Fig6_25_75_15bar_25C` | Fig. 6 (prediction) | 25/75 | 15 bar | 25 °C |
| `Fig6_75_25_15bar_25C` | Fig. 6 (prediction) | 75/25 | 15 bar | 25 °C |

> **Status.** Each breakthrough case asserts its H₂ and CO₂ breakthrough times
> (outlet mole fraction > 0.002, the paper's definition) against Table 5 as a
> relative error. Table 5 lists the paper's *simulated* times; the
> experimental curves (Figs. 4, 6) are not digitised yet.

## Key parameters

### Geometry & packing (Casas 2012, Table 2)

| Parameter | Value |
|-----------|-------|
| Column length `L` | 1.20 m |
| Inner diameter | 0.025 m |
| Bed void fraction `ε_b` | 0.403 |
| Particle density `ρ_p` | 850 kg/m³ |
| Skeletal density (He pycnometer) | 1970 kg/m³ → `ε_p = 1 − 850/1970 ≈ 0.568` |
| Particle diameter `d_p` | 0.003 m |
| Solid heat capacity `C_s` | 1000 J/(kg·K) |
| Lumped wall heat capacity `C_w` | 4×10⁶ J/(m³·K) (as `ρ_w·Cp_w = 8000·500`) |
| Lumped outer radius `R_o` | 0.02 m → wall thickness `R_o − R_i = 0.0075 m` |
| Heat of adsorption `ΔH_CO2 / ΔH_H2` | 26.0 / 9.8 kJ/mol |
| Mass transfer `k_CO2 / k_H2` (fitted) | 0.15 / 1.0 s⁻¹ |
| Internal HTC `h_L` | Leva (1947), `η1 = 41.13, η2 = 0.32` (fitted) |
| Wall HTC `h_W` (fitted) | 5.0 W/(m²·K) |

### Isotherm — exact competitive Sips via `ArrheniusSipsIsotherm`

The equilibrium paper (Schell 2012) describes each pure component with a Sips
isotherm whose parameters vary with temperature via Eqs. (5)–(7): Arrhenius `n∞`
and `k`, and (for CO₂) an S-shaped inverse-tangent exponent `c`. This exact form
is implemented as **`adsim::ArrheniusSipsIsotherm`** (new — added for this
validation), evaluated competitively (Ruthven 1984 empirical binary Sips):

```
q_i     = n∞_i · z_i / (1 + Σ_j z_j),   z_i = (k_i·P·y_i)^{c_i}
n∞_i(T) = a_i·exp(b_i/(R·T))                          [Arrhenius]
k_i(T)  = A_i·exp(B_i/(R·T))                          [Arrhenius]
c_i(T)  = c_ref,i + amp_i·(2/π)·arctan(β_i·(T − Tref,i))
```

The per-component constants (in `casas_2012.h → SetSchellSipsParams`) were
obtained by fitting these exact forms to the paper's per-temperature **Table 3**
values across **all five** measured temperatures (25–140 °C). The fit:

- **recovers the paper's own Table 4 Arrhenius magnitudes** to 3–4 sig figs
  (e.g. CO₂ `a = 1.383` vs paper 1.38, `b = 5626` vs ≈5628 J/mol; `A = 1.68e-8`
  vs `16.8e-9` 1/Pa, `B = 9160` vs 9159 J/mol) — confirming an earlier PDF
  extraction had simply flipped the sign of `b, B`;
- **reproduces Table 3 loadings** to ≤7 % (CO₂) and <1 % (H₂) across 0.5–40 bar
  and all five temperatures — the same fidelity as the paper's own Table-3→Table-4
  smoothing (CO₂'s `k` is non-monotonic at 45/65 °C).

Table 3 (Schell 2012), n∞ [mol/kg], k [1/MPa], c [-]:

| T | CO₂ n∞ / k / c | H₂ n∞ / k / c |
|---|---|---|
| 25 °C | 13.27 / 0.666 / 0.74 | 6.66 / 0.037 / 0.96 |
| 45 °C | 12.11 / 0.488 / 0.77 | 6.66 / 0.029 / 0.96 |
| 65 °C | 9.85 / 0.498 / 0.88 | 6.66 / 0.023 / 0.96 |
| 100 °C | 8.38 / 0.340 / 0.93 | — |
| 140 °C | 7.22 / 0.225 / 0.94 | — |

**Acceptance test** — `casas_2012_isotherm.cpp` gates the isotherm before any
process validation: `MatchesClosedForm` checks the C++ evaluates its documented
equation; `ReproducesTable3AcrossTemperatures` checks the continuous model
reproduces the paper's per-temperature Sips loadings across 25–140 °C.

### Heat transfer

- **Internal `h_L`** — Casas estimated `h_L` from **Leva's (1947)** packed-bed
  wall correlation, then fitted its coefficient/exponent to the reference run.
  Implemented as **`adsim::InternalHeatTransferModel_Leva`** (new):
  `Nu = h_L·D_t/k_g = η1·Re^{η2}·exp(−6·d_p/D_t)`. `η1, η2` default to Leva's
  canonical 0.813 / 0.90 and are set to the paper's fitted 41.13 / 0.32 here.
- **Wall `h_W`** — a fitted constant (5.0 W/m²K) in the paper's final model, so
  adsim's constant external-HTC model is exact; no correlation needed.

## Remaining work

- **Digitise breakthrough figures** (Figs. 4, 6–9) into `data/casas_2012/…` to
  compare against the measured curves (currently the comparison is against the
  paper's simulated Table 5 breakthrough times).
- **N₂** is measured in the equilibrium paper but is not part of the CO₂/H₂
  breakthrough study, so it is omitted (its Table 3 fit is available if needed).
- Wall/adsorbent thermal conductivity is not tabulated; a nominal 0.15 W/(m·K)
  is used and axial wall conduction is neglected (`thermal_conductivity_ = 0`).
- The exact Leva form (the `exp(−6·d_p/D_t)` geometric factor and whether Nu uses
  `k_g` or `K_L`) is reconstructed from the standard correlation — the paper's
  equation did not survive PDF extraction. The `η1/η2` interpretation is
  confirmed by the paper's "from correlations" values (0.813, 0.90).


## Validation

The tests compare the simulation with the paper using the helpers in
`../../fixtures/validation_metrics.h`; the published numbers they use are in
[`reference_values.csv`](reference_values.csv) (each with its table / figure /
page). Every tolerance is the measured error rounded up, with the measured
value in a comment next to it. Results of the last run:
[validation page](../../../../docs/validation.md#casas-2012)
(generated from `../validation_summary.json`).

- CO2 breakthrough times agree with the paper's simulated Table 5 values to 0.5-1.3 % for the 50/50 and 25/75 feeds and 13 % for 75/25; H2 breakthrough is 15-19 % early for the two CO2-lean feeds.
- Table 5 is the paper's simulation (feed velocity fitted per experiment, piping to the mass spectrometer modelled); adsim uses the 10 cm3/s set-point and no piping. The experimental curves (Figs. 4, 6) have not been digitised.
- The competitive Arrhenius-Sips isotherm reproduces the paper's per-temperature Sips fits (Schell et al., Table 3) to 7 % (CO2) and 0.6 % (H2).
