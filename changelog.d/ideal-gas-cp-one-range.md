### Changed — one ideal-gas heat capacity set per component, to 1000 K
- Each component's ideal-gas heat capacity polynomial (`ideal_gas_heat_capacity.h`) now holds one NASA TM-4513 coefficient set, the low-temperature one, valid to 1000 K. The 1000–6000 K set is dropped: adsorption processes stay well below 1000 K. Above 1000 K the value at 1000 K is used (as below the lower limit). Results do not change below 1000 K.
- C++: `NasaPolynomial` is `{t_min, t_max, a}` (a0..a6); `t_mid`, `low` and `high` are gone.
- The display name of the default gas heat capacity mode is now "Ideal gas" (was "Ideal gas Cp(T)"); the registry name is unchanged (`ideal_gas`).
