### Changed — changelog entries are fragments in `changelog.d/`
- A pull request adds its changelog entry as a new file
  `changelog.d/<topic>.md` instead of editing `CHANGELOG.md`, so pull
  requests in flight together no longer conflict on the top of *Unreleased*
  (GitHub ignores the union merge driver). `scripts/release.py` reads the
  fragments as part of *Unreleased*; the release pull request moves them into
  `CHANGELOG.md` and deletes them. New commands: `release.py fragments`
  (check) and `release.py collect`. CONTRIBUTING.md, RELEASING.md and the pull
  request template describe it.
