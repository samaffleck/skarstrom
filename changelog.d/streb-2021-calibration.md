### Changed — Streb 2021 decks calibrated on experiment A3
- The `Streb2021_A1` and `Streb2021_A3` decks, and the paper test's nine cases, are calibrated on the reference experiment A3 alone:
  - the CO2 and CH4 LDF coefficients: 0.0592 and 0.137 1/s (Table 4 gives 0.06 and 0.2, also fitted to A3);
  - the wall-to-thermofluid coefficient: 49.7 W/(m² K) (the paper takes 220 from a reference);
  - the sampling line to the mass spectrometer: 0.142 L of gas at column pressure. It is an inert, nearly empty layer above the bed (bed voidage 0.99, 40 cells), so the line carries the outlet gas as plug flow. A void region is one well-mixed cell and smears the CH4 front.
- Each case's feed flow is its Table 2 set point times a flow factor, fitted to the speed of its own CO2 temperature front between the thermocouples (A1 1.051, A3 1.020). The paper also fits the feed velocity per experiment. Six of the nine factors are more than 2 % from the set point (up to −11.3 % for B2).
- The other eight cases are predictions. Their CO2 fronts are now within 1.4 % of the measurement (before: 0.3–11 %) and their CH4 fronts within 2.4 %, except B1 at 4.7 % (before: 0.2–18 %). The thermocouple RMSE is now 5.5–11.4 K (before: 8.7–28.5 K). `core/test/papers/streb_2021/README.md`, "Calibration on A3", has the fit, the flow factors, the result held to 2 % of the set point and the limits. `streb_calibration.py` in the same folder repeats the calibration.
- `CellVariableAt` (test fixture `validation_metrics.h`) takes the log index of the layer's first cell, for beds of more than one layer.
