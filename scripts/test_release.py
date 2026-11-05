"""Tests for scripts/release.py (python -m unittest discover -s scripts -p 'test_release.py')."""

import sys
import unittest
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parent))
import release as r  # noqa: E402

CMAKE = 'x\nset(ADSIM_PROJECT_VERSION 0.3.0)\nset(ADSIM_VERSION_SUFFIX "")\ny\n'
CITATION = 'cff-version: 1.2.0\nversion: "0.3.0"\ndate-released: "2026-10-25"\n'


def changelog(unreleased: str) -> str:
    return ("# Changelog\n\nintro\n\n## [Unreleased]\n\n" + unreleased
            + "## [0.3.0] — 2026-10-25\n\n### Added\n- old\n")


class ReleaseTests(unittest.TestCase):
    def test_read_version(self):
        self.assertEqual(r.read_version(CMAKE), ("0.3.0", ""))
        self.assertEqual(r.read_version(CMAKE.replace('""', '"rc2"')), ("0.3.0", "rc2"))

    def test_bump_from_changelog(self):
        self.assertEqual(r.unreleased_kinds(changelog("### Added — x\n- a\n\n### Fixed\n- b\n\n")),
                         ["added", "fixed"])
        self.assertEqual(r.next_version("0.3.0", "", "auto", "", ["fixed", "added"]), "0.4.0")
        self.assertEqual(r.next_version("0.3.0", "", "auto", "", ["changed"]), "0.4.0")
        self.assertEqual(r.next_version("0.3.0", "", "auto", "", ["fixed", "security"]), "0.3.1")
        self.assertEqual(r.next_version("1.2.3", "", "auto", "", ["removed"]), "1.3.0")
        with self.assertRaises(r.ReleaseError):
            r.next_version("0.3.0", "", "auto", "", [])
        with self.assertRaises(r.ReleaseError):
            r.next_version("0.3.0", "", "auto", "", ["improved"])

    def test_explicit_bumps_and_prereleases(self):
        self.assertEqual(r.next_version("0.3.4", "", "patch", "", []), "0.3.5")
        self.assertEqual(r.next_version("0.3.4", "", "minor", "", []), "0.4.0")
        self.assertEqual(r.next_version("0.3.4", "", "major", "", []), "1.0.0")
        self.assertEqual(r.next_version("0.3.0", "", "auto", "rc1", ["added"]), "0.4.0rc1")
        self.assertEqual(r.next_version("0.4.0", "rc1", "auto", "rc2", []), "0.4.0rc2")
        self.assertEqual(r.next_version("0.4.0", "rc1", "auto", "", ["added"]), "0.4.0")
        with self.assertRaises(r.ReleaseError):
            r.next_version("0.4.0", "rc1", "auto", "rc1", [])
        with self.assertRaises(r.ReleaseError):
            r.next_version("0.3.0", "", "minor", "beta", [])

    def test_prepare_then_check(self):
        log = changelog("### Added — topic\n- new thing\n\n")
        self.assertFalse(r.unreleased_is_empty(log))
        self.assertTrue(r.check_ready("0.4.0", CMAKE, log, CITATION))  # not ready yet
        cmake = r.set_cmake_version(CMAKE, "0.4.0")
        log2 = r.cut_changelog(log, "0.4.0", "2026-11-05")
        cit = r.set_citation(CITATION, "0.4.0", "2026-11-05")
        self.assertEqual(r.read_version(cmake), ("0.4.0", ""))
        self.assertIn("## [Unreleased]\n\n## [0.4.0] — 2026-11-05\n\n### Added — topic\n"
                      "- new thing\n\n## [0.3.0]", log2)
        self.assertTrue(r.unreleased_is_empty(log2))
        self.assertEqual(r.section(log2, "0.4.0"), "### Added — topic\n- new thing")
        self.assertIn('version: "0.4.0"\ndate-released: "2026-11-05"', cit)
        self.assertEqual(r.check_ready("0.4.0", cmake, log2, cit), [])
        with self.assertRaises(r.ReleaseError):
            r.cut_changelog(log2, "0.4.1", "2026-11-05")  # Unreleased is empty

    def test_check_refuses_unlisted_changes(self):
        # The state main was in before this script: version 0.3.0 with
        # unreleased entries. Publishing it as 0.3.0 must be refused.
        problems = r.check_ready("0.3.0", CMAKE, changelog("### Fixed\n- x\n\n"), CITATION)
        self.assertEqual(len(problems), 1)
        self.assertIn("Unreleased", problems[0])

    def test_prerelease_check_skips_changelog(self):
        cmake = r.set_cmake_version(CMAKE, "0.4.0rc1")
        self.assertEqual(r.read_version(cmake), ("0.4.0", "rc1"))
        self.assertEqual(r.check_ready("0.4.0rc1", cmake, changelog("### Added\n- x\n\n"),
                                       CITATION), [])

    def test_fragments(self):
        log = changelog("### Fixed\n- old fix\n\n")
        merged = r.with_fragments(log, ["### Added — b\n- newest\n", "### Fixed — a\n- older\n"])
        self.assertIn("## [Unreleased]\n\n### Added — b\n- newest\n\n### Fixed — a\n- older\n\n"
                      "### Fixed\n- old fix\n\n## [0.3.0]", merged)
        self.assertEqual(r.unreleased_kinds(merged), ["added", "fixed", "fixed"])
        self.assertEqual(r.with_fragments(log, []), log)
        # Fragments alone make a release.
        only = r.with_fragments(changelog(""), ["### Fixed\n- x\n"])
        self.assertFalse(r.unreleased_is_empty(only))
        self.assertEqual(r.next_version("0.3.0", "", "auto", "", r.unreleased_kinds(only)), "0.3.1")
        self.assertEqual(r.fragment_problems("ok-name.md", "### Added — x\n- y\n"), [])
        self.assertTrue(r.fragment_problems("Bad Name.md", "### Added\n- y\n"))
        self.assertTrue(r.fragment_problems("x.md", "- no heading\n"))
        self.assertTrue(r.fragment_problems("x.md", "### Improved\n- y\n"))
        self.assertTrue(r.fragment_problems("x.md", "### Added\n## [0.5.0]\n"))

    def test_repository_files_parse(self):
        root = r.ROOT
        value, suffix = r.read_version((root / "CMakeLists.txt").read_text(encoding="utf-8"))
        text = (root / "CHANGELOG.md").read_text(encoding="utf-8")
        r.bump_kind(r.unreleased_kinds(text))  # no unknown subsection kinds
        for f in r.fragment_files():  # changelog.d/ fragments are well formed
            self.assertEqual(r.fragment_problems(f.name, f.read_text(encoding="utf-8")), [])
        if not suffix:
            self.assertIsNotNone(r.section(text, value),
                                 f"CHANGELOG.md has no section for the current version {value}")


if __name__ == "__main__":
    unittest.main()
