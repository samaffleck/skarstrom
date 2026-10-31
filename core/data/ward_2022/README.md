# Digitised data: Ward & Pini (2022)

Curves digitised from Fig. 2 of Ward & Pini, Adsorption 28:161–183 (2022),
doi:10.1007/s10450-022-00361-z, which plots the measured outlet CO₂ mole
fraction, outlet volumetric flow and the thermocouple at z = L/4 = 0.0715 m
for the three adsorption experiments of their ref. [5] (Pini, Joss, Hejazi
2021). One folder per feed flow (`Qin=50`, `Qin=100`, `Qin=200` cm³/min),
each with `outlet_CO2_mole_fraction[-]_...csv`, `outlet_volumetric_flowrate_...csv`
and `temperature[K]_z=0.0715[m]_...csv` (23–42 points per curve; time in s).

Used by `core/test/papers/ward_2022/ward_2022.cpp` (outlet composition and
thermocouple curves, front time at y = 0.5) and embedded in the bundled
`Ward2022_Q*ccm` decks as experimental data. The outlet flow series is not
compared yet.

Known limits: points were read off the published figure (the original
measurements of ref. [5] were not available to the repository); a digitised
point on the steep part of a front carries a time uncertainty of a few
seconds, which is below the differences the validation reports.
