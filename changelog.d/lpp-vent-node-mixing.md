### Fixed — light-product decks amplified a trace component through the vent
- When a step draws gas from a bed end through one stream (for example a light-product draw to
  another step) and another stream at the same end flows back in (a vent valve), the streams that
  leave now carry a mix: the bed's net outflow at the bed-end composition, plus the back-fed gas at
  the vent's own composition (its `y`, normalised). Before, the back-fed gas took the bed-end
  composition, so the draw carried about twice the gas that left the bed, all at the top composition.
  In `Nguyen2023_Exp_VI`–`X`, a component present only at round-off (for example water added to the
  deck) grew about 1.9 times per cycle and changed the cyclic steady state (#270). A vent whose `y`
  is all zero keeps the old behaviour.
- The vents of the light-product decks (`Nguyen2023_Exp_VI`–`X`, `Perez2019_D1`/`E`/`F`,
  `FiveStepProductPurge`, `Perez2019_CaseStudy_Base_LPP`) now carry N2, the feed's inert, as their
  backflow composition. Their results change slightly (see the pull request).
