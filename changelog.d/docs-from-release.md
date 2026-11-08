### Changed — documentation site shows the latest release
- The GitHub Pages root is built from the latest `v<semver>` tag with that tag's own configuration, so it
  describes the installable wheel. `main` is previewed under `/dev/` with a "not on PyPI yet" announcement bar.
  The site is rebuilt on every push to `main` and after every release.
- Convention: documented new names carry a "New in X.Y" note (AGENTS.md); applied to the run-health page and
  the `recovered_steps` statistic.
