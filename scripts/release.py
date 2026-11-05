#!/usr/bin/env python3
"""Release helper: read, compute and set the Skarstrom version.

The version lives in the root CMakeLists.txt (ADSIM_PROJECT_VERSION plus
ADSIM_VERSION_SUFFIX); CHANGELOG.md and CITATION.cff follow it. Used by
.github/workflows/prepare-release.yml (which opens the release pull request)
and .github/workflows/release.yml (which publishes a version once it reaches
main). See RELEASING.md. Standard library only.

  python scripts/release.py current             # 0.3.0
  python scripts/release.py next [--bump B] [--suffix rc1]
  python scripts/release.py prepare [--bump B] [--suffix rc1] [--date D]
  python scripts/release.py check VERSION       # ready to publish VERSION?
  python scripts/release.py notes VERSION       # CHANGELOG section of VERSION
  python scripts/release.py fragments           # check the changelog.d/ fragments
  python scripts/release.py collect             # fold them into CHANGELOG.md

Changelog fragments: a pull request does not edit CHANGELOG.md; it adds one
new file, changelog.d/<topic>.md, holding its "### Kind — topic" entry (see
changelog.d/README.md), so parallel pull requests never conflict on the top
of Unreleased. Every command here reads the fragments as part of Unreleased
(newest first, the order git added them), and `prepare` (a final release) and
`collect` move them into CHANGELOG.md and delete the files.

--bump is auto (default), patch, minor or major. auto reads the subsections
under "## [Unreleased]" in CHANGELOG.md: Added, Changed, Removed or Deprecated
make a minor release, Fixed or Security alone a patch release. A major release
is never chosen automatically. A pre-release (current version has a suffix)
is finished without a further bump.

Exit codes: 0 success, 1 error, 3 nothing to release (Unreleased is empty).
"""

from __future__ import annotations

import argparse
import datetime as _dt
import re
import subprocess
import sys
from pathlib import Path

ROOT = Path(__file__).resolve().parent.parent

VERSION_RE = re.compile(
    r'^set\(ADSIM_PROJECT_VERSION (?P<value>[0-9]+\.[0-9]+\.[0-9]+)\)\s*\n'
    r'set\(ADSIM_VERSION_SUFFIX "(?P<suffix>(?:(?:a|b|rc)[0-9]+)?)"\)',
    re.MULTILINE,
)
SUFFIX_RE = re.compile(r"^(?:(?:a|b|rc)[0-9]+)?$")
FULL_RE = re.compile(r"^(?P<value>[0-9]+\.[0-9]+\.[0-9]+)(?P<suffix>(?:a|b|rc)[0-9]+)?$")
UNRELEASED = "## [Unreleased]"
FRAGMENTS = ROOT / "changelog.d"
FRAGMENT_NAME_RE = re.compile(r"^[a-z0-9][a-z0-9._-]*\.md$")
MINOR_KINDS = {"added", "changed", "removed", "deprecated"}
PATCH_KINDS = {"fixed", "security"}

NOTHING_TO_RELEASE = 3


class ReleaseError(Exception):
    pass


# --- reading -----------------------------------------------------------------

def read_version(cmake_text: str) -> tuple[str, str]:
    m = VERSION_RE.search(cmake_text)
    if not m:
        raise ReleaseError(
            "CMakeLists.txt: the adjacent ADSIM_PROJECT_VERSION / ADSIM_VERSION_SUFFIX "
            "lines were not found")
    return m["value"], m["suffix"]


def split_changelog(text: str) -> tuple[str, str, str]:
    """(before and including the Unreleased heading, Unreleased body, rest)."""
    start = text.find(UNRELEASED + "\n")
    if start < 0:
        raise ReleaseError(f"CHANGELOG.md: no '{UNRELEASED}' heading")
    body_start = start + len(UNRELEASED) + 1
    nxt = re.compile(r"^## \[", re.MULTILINE).search(text, body_start)
    body_end = nxt.start() if nxt else len(text)
    return text[:body_start], text[body_start:body_end], text[body_end:]


def unreleased_kinds(changelog_text: str) -> list[str]:
    """Lower-case kinds of the '### Kind[ — topic]' subsections under Unreleased."""
    _, body, _ = split_changelog(changelog_text)
    kinds = []
    for line in body.splitlines():
        m = re.match(r"^### +([A-Za-z]+)", line)
        if m:
            kinds.append(m.group(1).lower())
    return kinds


def unreleased_is_empty(changelog_text: str) -> bool:
    _, body, _ = split_changelog(changelog_text)
    return not body.strip()


def section(changelog_text: str, version: str) -> str | None:
    m = re.search(rf"^## \[{re.escape(version)}\][^\n]*\n", changelog_text, re.MULTILINE)
    if not m:
        return None
    nxt = re.compile(r"^## \[", re.MULTILINE).search(changelog_text, m.end())
    return changelog_text[m.end():nxt.start() if nxt else len(changelog_text)].strip("\n")


# --- changelog fragments -----------------------------------------------------

def fragment_files(root: Path = ROOT) -> list[Path]:
    """changelog.d/*.md (not README.md), newest first by the commit that added each."""
    d = root / "changelog.d"
    if not d.is_dir():
        return []
    files = [f for f in d.glob("*.md") if f.name != "README.md"]

    def added(f: Path) -> float:
        try:
            out = subprocess.run(
                ["git", "log", "--diff-filter=A", "--format=%ct", "--", str(f.relative_to(root))],
                cwd=root, capture_output=True, text=True, check=True).stdout.split()
            return float(min(out, key=float)) if out else float("inf")  # untracked: newest
        except (OSError, subprocess.CalledProcessError):
            return float("inf")

    return sorted(files, key=lambda f: (-added(f), f.name))


def fragment_problems(name: str, text: str) -> list[str]:
    problems = []
    where = f"changelog.d/{name}"
    if not FRAGMENT_NAME_RE.match(name):
        problems.append(f"{where}: name it <topic>.md (lower case, digits, . _ -)")
    if not text.strip():
        problems.append(f"{where}: empty")
    elif not text.startswith("### "):
        problems.append(f'{where}: must start with a "### Kind — topic" heading')
    if re.search(r"^#{1,2} ", text, re.MULTILINE):
        problems.append(f'{where}: only "###" (and deeper) headings')
    for m in re.finditer(r"^### +([A-Za-z]+)", text, re.MULTILINE):
        if m.group(1).lower() not in MINOR_KINDS | PATCH_KINDS:
            problems.append(f"{where}: unknown kind {m.group(1)!r} "
                            "(use Added, Changed, Deprecated, Removed, Fixed or Security)")
    return problems


def with_fragments(changelog_text: str, fragments: list[str]) -> str:
    """CHANGELOG.md text with the fragments put at the top of Unreleased."""
    if not fragments:
        return changelog_text
    head, body, rest = split_changelog(changelog_text)
    block = "".join(f.strip("\n") + "\n\n" for f in fragments)
    return f"{head}\n{block}{body.lstrip(chr(10))}{rest}"


# --- computing ---------------------------------------------------------------

def bump_kind(kinds: list[str]) -> str | None:
    """'minor', 'patch' or None (nothing releasable) from the Unreleased kinds."""
    unknown = sorted(set(kinds) - MINOR_KINDS - PATCH_KINDS)
    if unknown:
        raise ReleaseError(
            "CHANGELOG.md: unknown subsection kind(s) under Unreleased: "
            + ", ".join(unknown)
            + " (use Added, Changed, Deprecated, Removed, Fixed or Security)")
    if any(k in MINOR_KINDS for k in kinds):
        return "minor"
    if kinds:
        return "patch"
    return None


def next_version(value: str, suffix: str, bump: str, new_suffix: str,
                 kinds: list[str]) -> str:
    if not SUFFIX_RE.match(new_suffix):
        raise ReleaseError(f"suffix {new_suffix!r} must be empty or aN, bN, rcN")
    major, minor, patch = (int(x) for x in value.split("."))
    if suffix and bump == "auto":
        # Finishing (or continuing) a pre-release: 0.4.0rc1 -> 0.4.0 / 0.4.0rc2.
        if new_suffix == suffix:
            raise ReleaseError(f"{value}{suffix} is already the current version")
        return f"{value}{new_suffix}"
    if bump == "auto":
        kind = bump_kind(kinds)
        if kind is None:
            raise ReleaseError("nothing to release: Unreleased in CHANGELOG.md is empty")
        bump = kind
    if bump == "major":
        major, minor, patch = major + 1, 0, 0
    elif bump == "minor":
        minor, patch = minor + 1, 0
    elif bump == "patch":
        patch += 1
    else:
        raise ReleaseError(f"unknown bump {bump!r}")
    return f"{major}.{minor}.{patch}{new_suffix}"


# --- writing -----------------------------------------------------------------

def set_cmake_version(cmake_text: str, full: str) -> str:
    m = FULL_RE.match(full)
    if not m:
        raise ReleaseError(f"bad version {full!r}")
    out, n = VERSION_RE.subn(
        f'set(ADSIM_PROJECT_VERSION {m["value"]})\n'
        f'set(ADSIM_VERSION_SUFFIX "{m["suffix"] or ""}")',
        cmake_text, count=1)
    if n != 1:
        raise ReleaseError("CMakeLists.txt: version lines not found")
    return out


def cut_changelog(changelog_text: str, full: str, date: str) -> str:
    head, body, rest = split_changelog(changelog_text)
    if not body.strip():
        raise ReleaseError("nothing to release: Unreleased in CHANGELOG.md is empty")
    if section(changelog_text, full) is not None:
        raise ReleaseError(f"CHANGELOG.md already has a section for {full}")
    return f"{head}\n## [{full}] — {date}\n\n{body.strip(chr(10))}\n\n{rest}"


def set_citation(citation_text: str, full: str, date: str) -> str:
    out, n1 = re.subn(r'^version: "[^"]*"', f'version: "{full}"', citation_text,
                      count=1, flags=re.MULTILINE)
    out, n2 = re.subn(r'^date-released: "[^"]*"', f'date-released: "{date}"', out,
                      count=1, flags=re.MULTILINE)
    if n1 != 1 or n2 != 1:
        raise ReleaseError("CITATION.cff: version or date-released not found")
    return out


def check_ready(full: str, cmake_text: str, changelog_text: str,
                citation_text: str) -> list[str]:
    """Problems that stop `full` from being published from this tree."""
    problems = []
    value, suffix = read_version(cmake_text)
    if f"{value}{suffix}" != full:
        problems.append(f"CMakeLists.txt has {value}{suffix}, not {full}")
    if suffix:
        return problems  # pre-releases do not cut the changelog
    if section(changelog_text, full) is None:
        problems.append(f"CHANGELOG.md has no '## [{full}]' section")
    if not unreleased_is_empty(changelog_text):
        problems.append(
            "CHANGELOG.md has entries under Unreleased: they would ship in "
            f"{full} without being listed in it")
    m = re.search(r'^version: "([^"]*)"', citation_text, re.MULTILINE)
    if not m or m.group(1) != full:
        problems.append(f"CITATION.cff version is {m.group(1) if m else '?'}, not {full}")
    return problems


# --- command line ------------------------------------------------------------

def _read(name: str) -> str:
    return (ROOT / name).read_text(encoding="utf-8")


def _write(name: str, text: str) -> None:
    (ROOT / name).write_text(text, encoding="utf-8", newline="\n")


def main(argv: list[str] | None = None) -> int:
    ap = argparse.ArgumentParser(description=__doc__.split("\n\n")[0])
    sub = ap.add_subparsers(dest="cmd", required=True)
    sub.add_parser("current", help="print the version in CMakeLists.txt")
    for name in ("next", "prepare"):
        p = sub.add_parser(name)
        p.add_argument("--bump", choices=["auto", "patch", "minor", "major"], default="auto")
        p.add_argument("--suffix", default="", help="pre-release suffix: aN, bN or rcN")
        if name == "prepare":
            p.add_argument("--date", default=_dt.date.today().isoformat())
    for name in ("check", "notes"):
        sub.add_parser(name).add_argument("version")
    sub.add_parser("fragments", help="check the changelog.d/ fragments")
    sub.add_parser("collect", help="move the changelog.d/ fragments into CHANGELOG.md")
    args = ap.parse_args(argv)

    try:
        cmake = _read("CMakeLists.txt")
        files = fragment_files()
        texts = [f.read_text(encoding="utf-8") for f in files]
        problems = [p for f, t in zip(files, texts) for p in fragment_problems(f.name, t)]
        if problems:
            raise ReleaseError("\n".join(problems))
        if args.cmd == "fragments":
            print(f"changelog.d: {len(files)} fragment(s), all well formed")
            return 0
        changelog = with_fragments(_read("CHANGELOG.md"), texts)
        if args.cmd == "collect":
            if files:
                _write("CHANGELOG.md", changelog)
                for f in files:
                    f.unlink()
            print(f"collected {len(files)} fragment(s) into CHANGELOG.md")
            return 0
        value, suffix = read_version(cmake)
        if args.cmd == "current":
            print(f"{value}{suffix}")
            return 0
        if args.cmd in ("next", "prepare"):
            # A final release needs changelog entries; a pre-release does not.
            if not args.suffix and unreleased_is_empty(changelog):
                print("nothing to release: Unreleased in CHANGELOG.md is empty",
                      file=sys.stderr)
                return NOTHING_TO_RELEASE
            full = next_version(value, suffix, args.bump, args.suffix,
                                unreleased_kinds(changelog))
            if args.cmd == "prepare":
                _write("CMakeLists.txt", set_cmake_version(cmake, full))
                if not args.suffix:
                    _write("CHANGELOG.md", cut_changelog(changelog, full, args.date))
                    _write("CITATION.cff", set_citation(_read("CITATION.cff"), full, args.date))
                    for f in files:  # now listed under the version's heading
                        f.unlink()
            print(full)
            return 0
        if args.cmd == "check":
            problems = check_ready(args.version, cmake, changelog, _read("CITATION.cff"))
            for p in problems:
                print(f"error: {p}", file=sys.stderr)
            return 1 if problems else 0
        if args.cmd == "notes":
            body = section(changelog, args.version)
            if body is None:
                if FULL_RE.match(args.version)["suffix"]:
                    print(f"Pre-release {args.version}. Changes since the last release "
                          "are under Unreleased in CHANGELOG.md.")
                    return 0
                raise ReleaseError(f"CHANGELOG.md has no section for {args.version}")
            print(body)
            return 0
    except ReleaseError as e:
        print(f"error: {e}", file=sys.stderr)
        return 1
    return 1


if __name__ == "__main__":
    sys.exit(main())
