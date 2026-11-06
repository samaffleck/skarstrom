### Fixed — gas heat capacity of condensable components
- In the CoolProp gas heat capacity mode, a component that is condensed in its
  pure state at the initial pressure and temperature (water, ethanol, hexane,
  toluene, ... at 1 bar and 298 K) now takes its ideal-gas heat capacity
  instead of CoolProp's liquid value: water vapour is about 1865 J/(kg·K),
  where the liquid's 4181 J/(kg·K) was used before. Components in a gas or
  supercritical state, including the permanent gases, are unchanged, as is
  the user-defined heat capacity mode. Viscosity and thermal conductivity
  already used vapour-phase fits. No bundled example deck has a condensable
  component, so their results do not change.
