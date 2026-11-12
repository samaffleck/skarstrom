# Data: Nguyen et al. (2022)

Breakthrough curves from Nguyen, Lin, Shimizu, Rajendran, *Separation of CO₂
and N₂ on a hydrophobic metal organic framework CALF-20*, Chemical
Engineering Journal 442:136263 (2022),
doi:[10.1016/j.cej.2022.136263](https://doi.org/10.1016/j.cej.2022.136263).
A 7.86 cm × 2.82 cm column holds 16.7 g of structured CALF-20 and starts full of
helium (Sect. 3.2).

- `Fig4/`: single-component CO₂ breakthrough, 50, 70 and 100 % CO₂ in helium at
  50 cm³/min (Fig. 4).
- `Fig5/`: CO₂/N₂ breakthrough, 5, 15 and 30 % CO₂ in nitrogen (no helium) at
  100 cm³/min (Fig. 5).

The flows, temperatures and pressures of each run are in Table A2.1 of the first
author's PhD thesis (University of Alberta); see `core/test/papers/nguyen_2022/README.md`.

## What the files are

**`Fig4/` and `Fig5/` are the measured markers**, re-digitised on 2026-11-12 from the
300 dpi raster figures in the paper's PDF by `digitise_figures.py` (this folder). The
script rerun on the same images gives these files byte for byte.

**`paper_model/Fig4/` and `paper_model/Fig5/` are the paper's model curves.** These
files were in `Fig4/` and `Fig5/` until 2026-11-12, described as the experimental
data. Overlaid on the figures, they follow the paper's solid simulation lines, not
its markers (most visibly the 50 % CO₂ front of Fig. 4, where the markers lie about
15 % earlier than the line). The paper test uses them for its "vs paper model"
(code-to-code) metrics. Their digitiser and method were not recorded (the files
predate the repository's history).

## Method (`digitise_figures.py`)

- **Axes.** Each panel is calibrated on its tick marks, found in the image: Fig. 4
  t̄ 0–300 and y 0–1 or 20–60 °C; Fig. 5 temperature t̄ 0–1000 and 20–50 °C. The
  time axes of Fig. 5 (a) and (c) are broken (0–20, then 200–1000 to the right of a
  gap), so they are calibrated as two linear segments.
- **Markers.** The fill colour of each series (red circles, blue triangles, green
  squares) is opened with a 7-px disc. This removes the model lines, which are about
  4–5 px wide, and keeps the marker fills. An isolated marker gives one point, at
  the centroid of its core.
- **Dense stretches.** Where markers overlap into a band of fill and black outlines
  (plateaus, knees, tails), the band is opened with an 11-px disc. Each band run is
  given to the series whose fill colour is in it, and its centre is sampled every
  6 px.
- **Legend boxes** are masked.
- **Baseline.** Before breakthrough the y_CO₂ markers of the three series overlap
  at y = 0, so colour cannot tell them apart. Points at y = 0 are therefore
  **placed**, not digitised: every 10 units of t̄, up to 5 units before each
  series' first resolved marker.
- **Time.** The figures plot dimensionless time t̄ = t v / L (Sect. 5.3). It is
  converted with each run's interstitial feed velocity, v = Q / (ε A), with ε 0.4
  and Q at the feed conditions: L / v = 23.56 s (50 cm³/min, Fig. 4) and 11.78 s
  (100 cm³/min, Fig. 5). The paper does not say which v it used. A standard-condition
  flow would change v by about 3 %.
- **Temperature** is converted from °C to K (+273.15).

Limits:

- Where many markers overlap, the band centre is a curve through them, not
  individual markers.
- Where the three temperature tails of Fig. 4 overlap (after t̄ ≈ 100), few points
  could be given to one series. The 100 % and 50 % thermocouple files end at
  t̄ ≈ 91 and 118.
- The plateau of the 100 % CO₂ curve at y = 1 lies on the plot frame and is not
  digitised.
- Expected error: about 1 px, i.e. about 0.5 units of t̄ (Fig. 4), 0.002 in y and
  0.1 K.

## Files

Every file is two columns, comma-separated: time in s and the quantity named in
the header. *Points* counts the data rows.

| File | Quantity (header) | Points | Time range [s] | Data range |
|---|---|---|---|---|
| `Fig4/T_CO2_100%.csv` | Time [s], T [K] | 28 | 451.9–2148 | 299.6–325.4 |
| `Fig4/T_CO2_50%.csv` | Time [s], T [K] | 39 | 606.8–2776 | 297.6–314.5 |
| `Fig4/T_CO2_70%.csv` | Time [s], T [K] | 34 | 637.5–6340 | 295.4–320.1 |
| `Fig4/Y_CO2_100%.csv` | Time [s], Y_CO2 [-] | 18 | 0–1275 | 0–0.9635 |
| `Fig4/Y_CO2_50%.csv` | Time [s], Y_CO2 [-] | 67 | 0–6815 | 0–0.5041 |
| `Fig4/Y_CO2_70%.csv` | Time [s], Y_CO2 [-] | 70 | 0–6964 | 0–0.7252 |
| `Fig5/T_CO2_15%.csv` | Time [s], T [K] | 38 | 121.9–7680 | 294.5–309.9 |
| `Fig5/T_CO2_30%.csv` | Time [s], T [K] | 35 | 224.9–3394 | 297–322.8 |
| `Fig5/T_CO2_5%.csv` | Time [s], T [K] | 70 | 150–11160 | 295.5–300.7 |
| `Fig5/Y_CO2_15%.csv` | Time [s], Y_CO2 [-] | 46 | 0–7870 | 0–0.1485 |
| `Fig5/Y_CO2_30%.csv` | Time [s], Y_CO2 [-] | 27 | 0–7229 | 0–0.2977 |
| `Fig5/Y_CO2_5%.csv` | Time [s], Y_CO2 [-] | 83 | 0–11180 | 0–0.0499 |
| `Fig5/Y_N2_15%.csv` | Time [s], Y_N2 [-] | 28 | 46.81–5467 | 0.0291–1.001 |
| `Fig5/Y_N2_30%.csv` | Time [s], Y_N2 [-] | 14 | 37.88–7241 | 0.0991–0.9589 |
| `Fig5/Y_N2_5%.csv` | Time [s], Y_N2 [-] | 48 | 39.85–11370 | 0.0738–0.9969 |

The thermocouple sits 1.82 cm from the column outlet (Sect. 3.2; the figure
captions say z = 0.75 L).

## Use

- `core/test/papers/nguyen_2022/nguyen_2022.cpp` compares each case with:
  - the measured files: outlet y_CO₂ (NRMSE over the feed fraction, and the front
    time at half the feed fraction), the thermocouple (RMSE) and, in Fig. 5, the
    outlet y_N₂ (RMSE);
  - `paper_model/`: the y_CO₂ NRMSE and the front time, code to code.
- Each of the six example decks `Nguyen2022_Fig4`, `_Fig4_50pct`, `_Fig4_100pct`,
  `Nguyen2022_Fig5`, `_Fig5_5pct` and `_Fig5_30pct` embeds every file of its
  figure's folder. It binds its own case's curves to the outlet mole fractions and
  to the temperature 1.82 cm from the outlet.
