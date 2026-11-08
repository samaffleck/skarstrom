# Releasing Skarstrom

Releases are automatic. Once a day, `.github/workflows/daily-release.yml`
releases `main` if it has changelog fragments and its CI is green: no pull
request to merge, no tag to push. The version comes from the fragments (see
*Versions*), and the release goes to PyPI within the day.

`.github/workflows/release.yml` does the build. It builds the wheels
(cibuildwheel: manylinux x86_64/aarch64, macOS arm64/x86_64 at deployment
target 13.3, Windows x64; CPython 3.10–3.15), an sdist that contains the
submodule sources, and a full source tarball. It tags the commit
`v<version>` and attaches them to a GitHub Release whose notes are the
version's `CHANGELOG.md` section. It publishes to PyPI or TestPyPI when
enabled.

## Documentation site

The site at <https://samaffleck.github.io/skarstrom/> is the documentation of the latest release tag, built
from that tag's own tree, so it matches what `pip install skarstrom` gives. `main` is previewed under `/dev/`
behind a "not on PyPI yet" banner. `.github/workflows/docs.yml` rebuilds both on every push to `main` and,
through `workflow_run`, after every successful *Release* run (the tag is created with the workflow token, which
starts no other workflow; the daily release dispatches *Release* as a top-level run, so this still fires), so a release refreshes the root automatically. To rebuild by hand, run the *Docs*
workflow.

## Versions

SemVer, written as PEP 440. The version lives in the root `CMakeLists.txt`:

```cmake
set(ADSIM_PROJECT_VERSION 0.9.1)
set(ADSIM_VERSION_SUFFIX ".dev0")   # "" released, "rc1"/"a1"/"b2" pre-release, ".dev0" between releases
```

The two lines must stay adjacent (pyproject.toml reads them with a regex).
Never edit them by hand, in a feature pull request or elsewhere: the release
workflow sets them. Changes accumulate as changelog fragments in
`changelog.d/` (one file per pull request, see `changelog.d/README.md`). The
release commit moves the fragments into `CHANGELOG.md`.
`scripts/release.py` reads the fragments as part of *Unreleased* everywhere
below.

### What the version strings mean

**`main` is never a released version.** Right after releasing 0.9.0, the
workflow sets `main` to the next development version, **0.9.1.dev0**. That is
PEP 440's "on the way to 0.9.1", which sorts after 0.9.0 and before any
release that follows it. A development build also names its commit, as a
PEP 440 local version:

| What you built | `skarstrom --version`, `skarstrom.__version__`, `adsim::kVersionString` |
|---|---|
| a wheel or sdist from PyPI, or a build at a release tag | `0.9.0` (or `0.10.0rc1`) |
| a checkout of `main` (git available) | `0.9.1.dev0+g0a19084` |
| `main` without git (e.g. the full source tarball of a dev build) | `0.9.1.dev0` |

- The local part (`+g<commit>`) is read by CMake when it configures
  (`git rev-parse --short HEAD`). Re-run CMake to refresh it after new
  commits. Only development versions get one, so the release build at the
  tag is always clean. `-DADSIM_VERSION_LOCAL=...` overrides it.
- The Python package metadata (`pip show skarstrom`) is `0.9.1.dev0`,
  without the commit: scikit-build-core reads it from the two lines above.
  PyPI does not accept local versions, which is right, because development
  builds are never published.
- The CMake package version (`find_package(Adsim 0.9)`) is the numeric part,
  `0.9.1`. CMake has no dev marker.
- After a pre-release `0.10.0rc1`, `main` is `0.10.0rc2.dev0`.
- CI fails (`python scripts/release.py check-dev`) if `CMakeLists.txt` has a
  version that is already tagged without a dev marker, so `main` can never
  quietly claim a released version again.

### Which version is next

The next version is chosen from the fragments and the *Unreleased* subsections
(`scripts/release.py`). The development version is a floor: from 0.9.1.dev0,

| Unreleased has | Next version (from 0.9.1.dev0, after 0.9.0) |
|---|---|
| any `### Added`, `### Changed`, `### Removed` or `### Deprecated` | minor: 0.10.0 |
| only `### Fixed` or `### Security` | patch: 0.9.1 |
| nothing | no release |

A topic after the kind (`### Added — wall design`) is fine. Before 1.0 an
incompatible change is a minor release too (SemVer allows anything in 0.y). A
major release is never chosen automatically: pick it by hand (below).

## The daily release

Every day at 05:23 UTC (and whenever it is run by hand) the *Daily release*
workflow:

1. **Checks.** It stops with a notice, and releases nothing, when:
   - main's latest CI run (`ci.yml`, for the newest commit that is not one of
     the workflow's own release commits) did not succeed, or is still running;
   - there are no changelog fragments and nothing under *Unreleased*;
   - on a scheduled run, `main` is on a pre-release line (`0.10.0rc2.dev0`);
     a pre-release holds the automatic releases until it is finished by hand;
   - on a scheduled run, the repository variable `AUTO_RELEASE` is `false`.
2. **Cuts the release.** It runs `python scripts/release.py prepare`, which
   sets the version in `CMakeLists.txt` and `CITATION.cff`, moves the
   fragments under `## [X.Y.Z] — <date>` in `CHANGELOG.md` and deletes them.
   `release.py check` must then pass. It commits that as **"Release X.Y.Z"**,
   then commits the development version (`release.py develop`, "Development
   version X.Y.(Z+1).dev0"). It pushes both to `main` (a plain push, never
   forced: if `main` moved during the run the push fails, and the next run
   releases from the new `main`).
3. **Starts the build.** It dispatches *Release* on `main` with
   `publish=true` and `ref=<the release commit>`. *Release* builds and tests
   every wheel, creates the tag `vX.Y.Z` at the release commit and the GitHub
   Release, then publishes. Then the *Docs* workflow rebuilds the site.

Why a dispatch: pushes and tags made with the workflow token (`GITHUB_TOKEN`)
start no workflows, but a `workflow_dispatch` made with it does (GitHub docs,
"Triggering a workflow from a workflow"). The daily workflow has
`actions: write` for this. *Release* could not simply be called as a reusable
workflow (`workflow_call`), because PyPI trusted publishing does not work
from a reusable workflow: "Reusable workflows cannot currently be used as the
workflow in a Trusted Publisher" (PyPI docs, *Troubleshooting*). As a
top-level run, *Release* keeps its trusted publisher (`release.yml`), its
environments (the deployment ref is `main`) and the *Docs* refresh after it,
all unchanged. No personal access token is needed.

The workflow also deletes the branch `bot/release`, left over from the
retired *Prepare release* workflow, which closes its pull request.

### By hand

- **Release now, or pick the bump:** Actions → Daily release → Run workflow
  on `main`, with `bump` (`patch`, `minor`, `major`; `auto` reads the
  fragments). The same checks apply, except the pre-release hold.
- **A pre-release:** the same, with `suffix` `rc1` (or `b1`, `a1`). It
  changes only `CMakeLists.txt` (the changelog is cut at the final release),
  publishes to TestPyPI, and sets `main` to `X.Y.Zrc2.dev0`. Scheduled runs
  then wait. Run the workflow again with an empty suffix to finish it
  (`X.Y.Z`, with the changelog), or with `rc2` for another pre-release.
- **Dry run:** tick `dry_run` (works on any branch). It does everything except
  push, dispatch and delete, and writes the release commit's diff and notes to
  the run summary.
- **Stop releases:** set the repository variable `AUTO_RELEASE` to `false`
  (Settings → Secrets and variables → Actions → Variables) to pause the
  scheduled runs, or disable the workflow (Actions → Daily release → ⋯ →
  Disable workflow). To stop one release in flight, cancel its *Release* run
  before the publish job. The release commits are then on `main` without a
  tag. Retry as below, or let the next fix ship as the next version.
- **Without the workflow:** `python scripts/release.py prepare` (`--bump`,
  `--suffix`) makes the same three-file change. Commit it as "Release X.Y.Z"
  in a pull request of its own. Merging it releases (the push to `main`
  starts *Release*). Then set the development version with
  `python scripts/release.py develop` in a second pull request.
- **If the Release run fails:** fix the cause and use *Re-run failed jobs* on
  that run, or run Actions → Release → Run workflow on `main` with *publish*
  ticked and `ref` set to the "Release X.Y.Z" commit. That releases the
  version at that commit if it has no tag yet; `ref` must be on `main`. If the
  release commit itself cannot be built, leave it: the fix ships as the next
  version, and the unreleased X.Y.Z simply has no tag.
- **Pushing a tag** `v<version>` by hand still releases that commit (the tag
  must equal the version in `CMakeLists.txt` there, and it may not be a
  development version).

`release.yml` refuses to publish a final version unless `CHANGELOG.md` has its
section, nothing is left under *Unreleased* and `CITATION.cff` carries it
(`python scripts/release.py check X.Y.Z`). So a version is never published
without a changelog that lists everything in it. A development version
(`.devN`) is never released, and a push to `main` that does not change the
version, or sets a development one, only builds.

- **Final versions** (`0.4.0`) publish to **PyPI**, environment `pypi`, only
  when the repository variable `PYPI_PUBLISH` is `true`.
- **Pre-releases** (`0.4.0rc1`, `0.4.0a1`, `0.4.0b1`) publish to
  **TestPyPI**, environment `testpypi`, only when `TESTPYPI_PUBLISH` is `true`,
  and mark the GitHub Release as a pre-release.
- A manual run of Release without *publish* builds everything but publishes
  nothing.

## PyPI storage

A release is about 145 MB (31 files: a wheel per platform and CPython version,
plus the sdist). PyPI's default project limit is 10 GB, which daily releasing
fills in about 70 releases. Before that happens:

- Request a larger project size limit: open an issue with the *Limit request*
  template at <https://github.com/pypi/support>, giving the project
  (`skarstrom`), the size per release and the release cadence.
- Old releases can be deleted on PyPI (Manage project → Releases); their files
  stay on the GitHub Release.
- Fewer files per release: stable-ABI (`abi3`) wheels would replace the six
  per-Python wheels of each platform with one.

## One-time setup: trusted publishing (no API tokens)

On PyPI (<https://pypi.org/manage/account/publishing/>) and on TestPyPI
(<https://test.pypi.org/manage/account/publishing/>) add a **pending publisher**
for the project `skarstrom` with exactly these values:

| Field | PyPI | TestPyPI |
|---|---|---|
| PyPI project name | `skarstrom` | `skarstrom` |
| Owner | `samaffleck` | `samaffleck` |
| Repository name | `skarstrom` | `skarstrom` |
| Workflow name | `release.yml` | `release.yml` |
| Environment name | `pypi` | `testpypi` |

Then, in the GitHub repository settings:

1. **Environments:** create `pypi` and `testpypi`. If you restrict their
   deployment branches and tags, allow both `main` (the daily release
   dispatches *Release* on `main`) and `v*` tags (the manual route); a
   `v*`-only rule makes PyPI publishing from `main` fail. Required reviewers
   on `pypi` would put a human approval back in front of every publish.
2. **Variables** (Settings → Secrets and variables → Actions → Variables): set
   `PYPI_PUBLISH` = `true` and/or `TESTPYPI_PUBLISH` = `true` when ready.
3. **Pages:** set Settings → Pages → Source to *GitHub Actions* for the docs
   (`.github/workflows/docs.yml`).
4. **Daily release:** nothing to set up. The workflow pushes to `main` with
   the workflow token (`contents: write`) and dispatches *Release*
   (`actions: write`). If `main` is ever protected by a branch protection rule
   or ruleset, it must let GitHub Actions push: add the *GitHub Actions* app
   (or the repository's workflows) to the bypass list, or allow pushes that do
   not come through a pull request for it. (Today `main` has no protection.)
   Optional: the variable `AUTO_RELEASE` = `false` pauses scheduled releases.
   The secret `RELEASE_PR_TOKEN`, used by the retired Prepare release
   workflow, is no longer needed and can be deleted.

No API tokens or secrets are needed: the publish jobs request a short-lived
OIDC token (`id-token: write`) that PyPI exchanges for an upload credential.
