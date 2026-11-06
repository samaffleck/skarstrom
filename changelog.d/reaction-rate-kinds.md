### Added — reactions on a catalyst solid (data model and API)
- A solid (adsorbent) can carry heterogeneous reactions
  (`Adsorbent::reaction_source_`, `adsim/objects/reaction.h`): stoichiometry,
  a rate law per kg of the solid, a rate basis (partial pressure in bar, or
  concentration), a heat of reaction (constant, or Kirchhoff with the gas
  heat capacities the energy balance uses) and an effectiveness factor.
  A catalyst mixed into a layer with a sorbent makes a sorption-enhanced
  reactor. The reactions do not enter the balances of a run yet.
- Rate kinds (`adsim/equations/reaction_kinds.h`): `power_law`,
  `reversible_power_law`, `lhhw` (Langmuir-Hinshelwood-Hougen-Watson with any
  number of adsorption terms; reproduces the Xu and Froment steam-reforming
  rates) and `expression` (a formula of T, P and each component's `p_<C>`,
  `c_<C>`, `y_<C>`, which expressions can now use).
- Element formulas for the built-in components and an element-balance check
  of every reaction.
- Builder: `AddReaction`, `SetReactionKind`, `SetReactionStoichiometry`,
  `SetReactionBasis`, `SetEffectiveness`, `AddAdsorptionTerm`,
  `RemoveAdsorptionTerm`, `RemoveReaction`, `ReactionRateKinds`. Registry
  paths `adsorbent[ID].reaction[ID].*` (category `reaction`).
- Deck checks (`ReactionChecks`) and warnings (`ReactionWarnings`).
- `RemoveComponentFromDeck(data, name, err)` refuses to remove a component a
  reaction uses.

### Changed — adsorbent archive version 5
- `Adsorbent` is saved at cereal version 5 with a `Reaction_Source` node.
  Decks and library files at version 4 load unchanged, with no reactions.
- The library files in `data/adsorbents/` and `examples/decks/Illustrative_TSA_13X.json`
  are re-saved by the test suite at version 5 (version stamp and an empty
  `Reaction_Source` only).
