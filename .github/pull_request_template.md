## What and why

<!-- What does this change, and why? Link the issue it fixes. -->

## Checklist

- [ ] Commits are signed off (`git commit -s`, see CONTRIBUTING.md)
- [ ] Tests added or updated; `ctest` (fast subset at least) and `pytest python/tests` pass
- [ ] Physics changes are checked against a derived or published reference, not recorded output
- [ ] Serialised classes changed? `cereal_class_version` bumped and old decks still load
- [ ] Adsorbent/validation data: every value cited; generated pages regenerated (`scripts/docs/`), not hand-edited
- [ ] Regenerated example decks (`examples/decks/`) committed if the tests changed them
- [ ] Docs (`docs/`) and `CHANGELOG.md` updated
- [ ] `core/Results/` not committed
