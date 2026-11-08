#!/usr/bin/env bash
# Build the Pages site: the latest release tag at the root, main under dev/.
#   scripts/docs/build_site.sh [SITE_DIR [PREBUILT_DEV_DIR]]   (default site; run from a full clone with tags)
#   scripts/docs/build_site.sh --dev-version                    print the banner's version string
# With PREBUILT_DEV_DIR (a --strict build of mkdocs.dev.yml made elsewhere, e.g. by the docs-checks
# action) it is moved to SITE_DIR/dev instead of building main again.
# The release is built from its own tree (its mkdocs.yml, requirements and pages) in a worktree;
# main is built from this checkout through mkdocs.dev.yml (banner, site_url .../dev/).
# The released root gets a footer note through a throwaway config that sets theme.custom_dir; the tag is not edited.
# SKIP_INSTALL=1 skips pip install (local runs).
set -euo pipefail

root=$(git rev-parse --show-toplevel)
cd "$root"
# Version for the banner: CMake version (+ suffix) and short commit.
ver=$(sed -n 's/^set(ADSIM_PROJECT_VERSION \([^ )]*\)).*/\1/p' CMakeLists.txt)
suffix=$(sed -n 's/^set(ADSIM_VERSION_SUFFIX "\([^"]*\)").*/\1/p' CMakeLists.txt)
export DEV_VERSION="${ver}${suffix}+$(git rev-parse --short HEAD)"

if [ "${1:-}" = "--dev-version" ]; then echo "$DEV_VERSION"; exit 0; fi
site=$(realpath -m "${1:-site}")
base=https://samaffleck.github.io/skarstrom

tag=$(git tag -l 'v[0-9]*' | grep -E '^v[0-9]+\.[0-9]+\.[0-9]+$' | sort -V | tail -n 1 || true)
if [ -z "$tag" ]; then
  echo "error: no release tag v<major>.<minor>.<patch> found (fetch tags: git fetch --tags)" >&2
  exit 1
fi
echo "Latest release: $tag"

install() { [ -n "${SKIP_INSTALL:-}" ] || python -m pip install -r "$1"; }

rm -rf "$site"
src=$(mktemp -d)/release
git worktree add --detach "$src" "$tag" >/dev/null
trap 'git worktree remove --force "$src" >/dev/null 2>&1 || true' EXIT

# Released docs -> root. Footer note via main's overrides (custom_dir) and extra.dev_url, both set in a throwaway config.
install "$src/requirements-docs.txt"
cat > "$src/mkdocs.root.yml" <<YML
INHERIT: mkdocs.yml
theme:
  custom_dir: $root/overrides
extra:
  dev_url: $base/dev/
YML
(cd "$src" && mkdocs build --strict -f mkdocs.root.yml --site-dir "$site")

# main -> dev/ (strict, with the banner).
if [ -n "${2:-}" ]; then
  mv "$2" "$site/dev"
else
  install requirements-docs.txt
  mkdocs build --strict -f mkdocs.dev.yml --site-dir "$site/dev"
fi
