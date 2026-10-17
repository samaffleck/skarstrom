# Releasing Skarstrom

Releases are made from tags. `.github/workflows/release.yml` builds the wheels
(cibuildwheel: manylinux x86_64/aarch64, macOS arm64/x86_64 at deployment target
13.3, Windows x64; CPython 3.10–3.13), an sdist that contains the submodule
sources, and a full source tarball; attaches them to a GitHub Release; and
publishes to PyPI or TestPyPI when enabled.

## Versions

SemVer. The version lives in the root `CMakeLists.txt`:

```cmake
set(ADSIM_PROJECT_VERSION 0.1.0)
set(ADSIM_VERSION_SUFFIX "")        # "rc1", "a1", "b2" for a pre-release
```

The two lines must stay adjacent (pyproject.toml reads them with a regex). The
tag must be `v` + version + suffix, e.g. `v0.1.0` or `v0.2.0rc1`; the workflow
fails if they differ.

## Making a release

1. Update `CHANGELOG.md` (move *Unreleased* under the new version with the date)
   and `version` in `CITATION.cff`.
2. Set the version (and suffix for a pre-release) in `CMakeLists.txt`.
3. Merge to `main` with CI green (the push to `main` runs the full test suite).
4. Tag and push: `git tag -s v0.2.0 -m "Skarstrom 0.2.0" && git push origin v0.2.0`.
5. Check the GitHub Release the workflow creates, and the PyPI/TestPyPI upload.
6. Set the suffix back to `""` (or bump the version) for development.

- **Final tags** (`v0.2.0`) publish to **PyPI**, environment `pypi`, only when
  the repository variable `PYPI_PUBLISH` is `true`.
- **Pre-release tags** (`v0.2.0rc1`, `v0.2.0a1`, `v0.2.0b1`) publish to
  **TestPyPI**, environment `testpypi`, only when `TESTPYPI_PUBLISH` is `true`,
  and mark the GitHub Release as a pre-release.
- A manual run (Actions → Release → Run workflow) builds everything but
  publishes nothing.

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
   reviewers and a `v*` tag rule).
2. **Variables** (Settings → Secrets and variables → Actions → Variables): set
   `PYPI_PUBLISH` = `true` and/or `TESTPYPI_PUBLISH` = `true` when ready.
3. **Pages:** set Settings → Pages → Source to *GitHub Actions* for the docs
   (`.github/workflows/docs.yml`).

No API tokens or secrets are needed: the publish jobs request a short-lived
OIDC token (`id-token: write`) that PyPI exchanges for an upload credential.
