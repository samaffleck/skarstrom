### Added — axial heat dispersion as a stagnant term plus a dynamic term
- Two new packed-bed heat dispersion kinds, `krupiczka` and `zehner_schlunder`: lambda = k0 + c k_g Pr Re. The kind chooses the stagnant term k0: the Krupiczka-type term of Yagi and Kunii's correlation, or Zehner–Schlünder (Zehner and Schlünder 1970, shape factor C = 1.25). The new parameter `dynamic_coefficient` gives c. Literature values: 0.75 Yagi–Kunii (default of `krupiczka`), 0.5 Wakao–Kaguei (default of `zehner_schlunder`), 1/K_ax Bauer–Schlünder. The Zehner–Schlünder term is also available to other correlations as `ZehnerSchlunderConductivityRatio`; near kappa = B it uses a series form instead of the closed form.

### Changed — `yagi_kunii` heat dispersion is legacy
- `yagi_kunii` is `krupiczka` at `dynamic_coefficient` 0.75. Decks that use it load and give the same results, bit for bit. The builder still accepts it, but the kind lists and the manifest no longer offer it (`legacy`).
