### Changed — layer compositions are checked
- **Behaviour change:** a deck whose layer mass fractions do not sum to 1 now
  fails the deck checks instead of running with scaled properties. The
  message names the layer, the sum and the fix, e.g. `Layer "Bottom": mass
  fractions sum to 0.8; they must sum to 1`. No bundled deck is affected.
- A layer whose mass fractions do not sum to 1 (within 1e-6), or that holds a
  negative fraction, is now refused by the deck checks, with a message naming
  the layer and the sum. Every property of a layer's solid scales with the
  sum, so a sum of 0.8 used to mean 80 % of the solid without notice. The
  fractions are never renormalised.
- A solid at mass fraction 0 is allowed but reported as a warning (it takes no
  part in the run).
- The rules live in one engine function, `CheckLayerComposition`
  (`porous_media.h`), which the deck checks (C++ runs, Python `sim.check()`),
  `LayerCompositionWarnings` (`deck_warnings.h`) and front ends share.
- `RemoveAdsorbentFromDeck` refuses to remove a solid that is part of a blend
  (a layer that also holds another solid) and says which layers use it, like
  removing a component that a reaction uses. A layer that holds only that
  solid is still re-seeded with another one. New checked overload
  `RemoveAdsorbentFromDeck(data, id, err)` and `AdsorbentInUseByMixtures`.
- New blend edits in `deck_edit.h` that keep a layer's fractions summing to 1
  (the other solids are rescaled in proportion): `AddSolidToLayer`,
  `SetLayerMassFraction`, `RemoveSolidFromLayer`. `SeedLayerComposition` no
  longer adds a second solid at 1.0 to a layer that already holds one.
- The solver logs the composition warnings at the start of a run, with the
  other deck warnings.

### Added — verification of mixed-layer edges
- `MixedLayers.*`: a solid at mass fraction 0 changes nothing; a layer of an
  inert solid alone holds nothing and balances; splitting a layer into two of
  the same blend gives the same solution bit for bit; the mole balance closes
  across the interface between two layers of different blends.
