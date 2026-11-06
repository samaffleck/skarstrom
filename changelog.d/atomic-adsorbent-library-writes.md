### Fixed — library re-saves under `ctest --parallel`
- The tests that re-save `data/adsorbents/*.json` now write through the
  suite's atomic-write helper (`test::SaveAdsorbentToJsonAtomically`): the
  rename over the target is retried while a reader holds the file open
  (Windows), and a failed write throws instead of leaving the old file in
  place unnoticed. The library files are unchanged.
