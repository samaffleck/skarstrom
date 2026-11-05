### Fixed — flaky tests reading example decks under `ctest --parallel`
- The tests that re-save `examples/decks/*.json` now write each deck to a
  temporary file in the same directory and rename it over the deck, so a test
  loading the deck at the same time always reads a complete file instead of
  truncated JSON (`rapidjson internal assertion failure: IsObject()`). On
  Windows, where the replace fails while a reader holds the file open, the
  rename is retried for up to about 5 s. The deck contents are unchanged.
