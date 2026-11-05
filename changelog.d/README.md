# Changelog fragments

Do not edit `CHANGELOG.md` in a pull request. Add one new file here instead:

```
changelog.d/<topic>.md        # e.g. changelog.d/macropore-capacity-factor.md
```

It holds your entry exactly as it would have appeared under *Unreleased*: one
or more subsections whose kind is `Added`, `Changed`, `Deprecated`, `Removed`,
`Fixed` or `Security` (the kinds choose the next version, see RELEASING.md),
optionally with a topic:

```markdown
### Fixed — breakthrough decks stall on fine grids
- ...
```

Two pull requests never add the same new file, so they cannot conflict
(GitHub ignores the union merge driver `.gitattributes` sets on
`CHANGELOG.md`). To change an entry that is still a fragment, edit that
fragment.

`scripts/release.py` treats the fragments as part of *Unreleased* (newest
first, the order git added them). The release pull request moves them into
`CHANGELOG.md` under the new version and deletes them;
`python scripts/release.py collect` does the same without a release, and
`python scripts/release.py fragments` checks them (CI runs it through
`scripts/test_release.py`).
