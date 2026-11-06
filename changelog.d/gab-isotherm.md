### Added — GAB and BET water isotherms
- Two isotherm kinds for water and other condensable vapours, in the relative
  pressure x = p_i / p_sat(T): `gab` (Guggenheim-Anderson-de Boer with
  C = C0 exp(dH_C/RT), K = K0 exp(dH_K/RT), q_m = q_m0 exp(dH_m/RT); BET is
  K0 = 1, dH_K = 0) and `gab_young` (the layer-energy form of Young et al.
  2021 and Low et al. 2023, with the papers' own parameters). The saturation
  pressure comes from CoolProp's saturation ancillary (water: within 0.02 % of
  IAPWS-95). The loading is exact for 0 <= x <= min(1, 0.95/K) and continued
  smoothly beyond (the Henry line for x < 0, the tangent above), so it stays
  finite at and above saturation. Not competitive: water adsorbs on its own.
- `adsim::IsostericHeat(isotherm, component, T, P, y)`: the isosteric heat
  implied by any isotherm's temperature dependence, to help choose the
  constant heat of adsorption the energy balance uses.
- The GAB kinds work under the IAST mixing rule (closed-form spreading
  pressure).

### Fixed — IAST with an isotherm kind it does not know
- An isotherm kind without an IAST form now raises an error instead of being
  treated as inert (which silently gave that component zero loading).
