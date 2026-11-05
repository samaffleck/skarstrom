# Releasing Skarstrom

Releases are automatic once a version change reaches `main`.
`.github/workflows/release.yml` builds the wheels (cibuildwheel: manylinux
x86_64/aarch64, macOS arm64/x86_64 at deployment target 13.3, Windows x64;
CPython 3.10–3.15), an sdist that contains the submodule sources, and a full
source tarball; tags the commit `v<version>`; attaches them to a GitHub Release
whose notes are the version's `CHANGELOG.md` section; and publishes to PyPI or
TestPyPI when enabled.

## Versions

SemVer. The version lives in the root `CMakeLists.txt`:

```cmake
set(ADSIM_PROJECT_VERSION 0.1.0)
set(ADSIM_VERSION_SUFFIX "")        # "rc1", "a1", "b2" for a pre-release
```

The two lines must stay adjacent (pyproject.toml reads them with a regex).
**The version on `main` is the last released version** (or the pre-release
being tested): do not bump it for development, and do not bump it by hand in a
feature pull request. Changes accumulate as changelog fragments in
`changelog.d/` (one file per pull request, see `changelog.d/README.md`) and
under *Unreleased* in `CHANGELOG.md` until the next release; the release
commit moves the fragments into `CHANGELOG.md`. `scripts/release.py` reads the
fragments as part of *Unreleased* everywhere below.

The next version is chosen from the *Unreleased* subsections
(`scripts/release.py`):

| Unreleased has | Next version (from 0.3.0) |
|---|---|
| any `### Added`, `### Changed`, `### Removed` or `### Deprecated` | minor: 0.4.0 |
| only `### Fixed` or `### Security` | patch: 0.3.1 |
| nothing | no release |

A topic after the kind (`### Added — wall design`) is fine. Before 1.0 an
incompatible change is a minor release too (SemVer allows anything in 0.y). A
major release is never chosen automatically: pick it by hand.

## Making a release

1. **The release pull request is already open.** Whenever there are
   changelog fragments or entries under *Unreleased*, `.github/workflows/prepare-release.yml` keeps one
   pull request, "Release X.Y.Z", open on the bot branch `bot/release`. It sets
   the version in `CMakeLists.txt` and `CITATION.cff` and moves the
   fragments and *Unreleased* entries under `## [X.Y.Z] — <date>` (deleting
   the fragment files). It is rebuilt from `main` whenever the changelog
   changes there (so it carries every merged change and
   today's date), and weekly. Never push to `bot/release`.
2. **Merge it.** That is the release: the push to `main` changes the version,
   `release.yml` builds and tests every wheel, creates the tag and the GitHub
   Release, then publishes. Check the run and the PyPI page.

Release when there is something worth shipping, or on a regular cadence (for
example, merge the release pull request weekly). Nothing else needs doing
after a release: the next pull request's changelog entry reopens the release
pull request.

- **Another bump or a pre-release:** Actions → Prepare release → Run workflow,
  with `bump` (patch, minor, major) and/or `suffix` (`rc1`, `b1`, `a1`). The
  pull request keeps that choice until it is merged or closed. A pre-release
  only changes `CMakeLists.txt` (the changelog is cut at the final release) and
  publishes to TestPyPI; the next automatic refresh then proposes the final
  version.
- **By hand, without the workflow:** `python scripts/release.py prepare`
  (`--bump`, `--suffix`) makes the same three-file change; commit it as
  "Release X.Y.Z" in a pull request of its own.
- **If the release run fails** after the merge, fix the cause and use *Re-run
  failed jobs* on that run, or, once a fix is on `main`, Actions → Release →
  Run workflow on `main` with *publish* ticked (it releases the version in
  `CMakeLists.txt` if it has no tag yet).
- **Pushing a tag** `v<version>` by hand still releases that commit, as before
  (the tag must equal the version in `CMakeLists.txt` there).

`release.yml` refuses to publish a final version unless `CHANGELOG.md` has its
section, nothing is left under *Unreleased* and `CITATION.cff` carries it
(`python scripts/release.py check X.Y.Z`), so `main` is never published under
a version whose changelog does not list everything in it. A push to `main` that
does not change the version, or changes it to one already tagged, only builds.

- **Final versions** (`0.4.0`) publish to **PyPI**, environment `pypi`, only
  when the repository variable `PYPI_PUBLISH` is `true`.
- **Pre-releases** (`0.4.0rc1`, `0.4.0a1`, `0.4.0b1`) publish to
  **TestPyPI**, environment `testpypi`, only when `TESTPYPI_PUBLISH` is `true`,
  and mark the GitHub Release as a pre-release.
- A manual run of Release without *publish* builds everything but publishes
  nothing.

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

1. **Environments:** create `pypi` and `testpypi` (optionally with required
   reviewers). If you restrict their deployment branches and tags, allow both
   `main` (the normal release route) and `v*` tags (the manual one); a `v*`-only
   rule makes PyPI publishing from `main` fail.
2. **Variables** (Settings → Secrets and variables → Actions → Variables): set
   `PYPI_PUBLISH` = `true` and/or `TESTPYPI_PUBLISH` = `true` when ready.
3. **Pages:** set Settings → Pages → Source to *GitHub Actions* for the docs
   (`.github/workflows/docs.yml`).
4. **Release pull request:** enable Settings → Actions → General → *Allow
   GitHub Actions to create and approve pull requests*, so the Prepare release
   workflow can open its pull request. Pull requests opened with the workflow
   token start no CI; to have CI run on the release pull request, add a secret
   `RELEASE_PR_TOKEN` (a fine-grained personal access token or GitHub App token
   with *Contents* and *Pull requests* read and write on this repository).

No API tokens or secrets are needed: the publish jobs request a short-lived
OIDC token (`id-token: write`) that PyPI exchanges for an upload credential.
