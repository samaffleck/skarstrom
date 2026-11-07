### Fixed — slow-tier CSS tests on the Perez 2019 decks
- `CSSAndersonDecks` and `CSSDistanceDecks.Perez2019EReportsTheRemainingDistance`
  use stops and reference KPIs re-measured on the current `Perez2019_A` and
  `Perez2019_E` decks. The decks' constant gas heat capacity moved the plain
  stops (A: 48 → 47, E: 205 → 206) and the long-run reference KPIs, so the
  tests' earlier constants no longer matched. No tolerance changed.
