"""Tests for scripts/release.py (python -m unittest discover -s scripts -p 'test_release.py')."""

import re
import subprocess
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

    def test_dev_versions(self):
        dev = CMAKE.replace('0.3.0)', '0.3.1)').replace('""', '".dev0"')
        self.assertEqual(r.read_version(dev), ("0.3.1", ".dev0"))
        self.assertEqual(r.read_version(CMAKE.replace('""', '"rc2.dev0"')), ("0.3.0", "rc2.dev0"))
        self.assertEqual(r.dev_version("0.3.0"), "0.3.1.dev0")
        self.assertEqual(r.dev_version("0.4.0rc1"), "0.4.0rc2.dev0")
        self.assertEqual(r.dev_version("1.0.0b3"), "1.0.0b4.dev0")
        with self.assertRaises(r.ReleaseError):
            r.dev_version("0.3.1.dev0")  # already a development version
        self.assertEqual(r.read_version(r.set_cmake_version(CMAKE, "0.3.1.dev0")), ("0.3.1", ".dev0"))
        self.assertEqual(r.read_version(r.set_cmake_version(dev, "0.4.0")), ("0.4.0", ""))
        self.assertTrue(r.is_dev("0.3.1.dev0") and r.is_dev("0.4.0rc2.dev0"))
        self.assertFalse(r.is_dev("0.3.1") or r.is_dev("0.4.0rc2"))

    def test_bump_from_dev_floor(self):
        # 0.3.1.dev0 follows 0.3.0: the dev version is the floor, never skipped over.
        self.assertEqual(r.next_version("0.3.1", ".dev0", "auto", "", ["fixed"]), "0.3.1")
        self.assertEqual(r.next_version("0.3.1", ".dev0", "auto", "", ["added", "fixed"]), "0.4.0")
        self.assertEqual(r.next_version("0.3.1", ".dev0", "major", "", []), "1.0.0")
        self.assertEqual(r.next_version("0.3.1", ".dev0", "auto", "rc1", ["added"]), "0.4.0rc1")
        self.assertEqual(r.next_version("0.4.0", ".dev0", "minor", "", []), "0.4.0")
        self.assertEqual(r.next_version("1.0.0", ".dev0", "major", "", []), "1.0.0")
        self.assertEqual(r.next_version("0.3.1", ".dev0", "patch", "", []), "0.3.1")
        with self.assertRaises(r.ReleaseError):
            r.next_version("0.3.1", ".dev0", "auto", "", [])  # nothing to release
        # A pre-release line: 0.4.0rc2.dev0 follows 0.4.0rc1.
        self.assertEqual(r.next_version("0.4.0", "rc2.dev0", "auto", "", ["added"]), "0.4.0")
        self.assertEqual(r.next_version("0.4.0", "rc2.dev0", "auto", "rc2", []), "0.4.0rc2")
        self.assertEqual(r.next_version("0.4.0", "rc2.dev0", "auto", "rc3", []), "0.4.0rc3")
        with self.assertRaises(r.ReleaseError):
            r.next_version("0.4.0", "rc2.dev0", "auto", "rc1", [])  # already out
        with self.assertRaises(r.ReleaseError):
            r.next_version("0.4.0", "rc2.dev0", "auto", "b1", [])

    def test_check_dev(self):
        tags = {"v0.3.0", "v0.4.0rc1"}
        self.assertEqual(r.check_dev("0.3.1.dev0", tags), [])
        self.assertEqual(r.check_dev("0.4.0rc2.dev0", tags), [])
        self.assertEqual(r.check_dev("0.3.1", tags), [])  # a release being cut by hand
        problems = r.check_dev("0.3.0", tags)
        self.assertEqual(len(problems), 1)
        self.assertIn("0.3.1.dev0", problems[0])
        self.assertTrue(r.check_dev("0.4.0rc1", tags))
        # Development versions are never published.
        dev = r.set_cmake_version(CMAKE, "0.3.1.dev0")
        self.assertTrue(r.check_ready("0.3.1.dev0", dev, changelog(""), CITATION))

    def test_pyproject_regex_reads_dev_versions(self):
        # scikit-build-core reads the Python version with this regex; it must
        # accept every form release.py writes and give the same string.
        text = (r.ROOT / "pyproject.toml").read_text(encoding="utf-8")
        m = re.search(r"^regex = '(.*)'$", text, re.MULTILINE)
        self.assertIsNotNone(m, "pyproject.toml: [tool.scikit-build.metadata.version] regex")
        rx = re.compile(m.group(1), re.MULTILINE)
        for full in ("0.3.0", "0.4.0rc1", "0.3.1.dev0", "0.4.0rc2.dev0"):
            g = rx.search(r.set_cmake_version(CMAKE, full))
            self.assertIsNotNone(g, full)
            self.assertEqual(g["value"] + g["suffix"], full)

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

    def test_notes_survive_a_closed_pipe(self):
        # `release.py notes | head -n 0` must not fail with BrokenPipeError.
        script = Path(__file__).with_name("release.py")
        p = subprocess.Popen([sys.executable, str(script), "notes", "0.8.0"],
                             stdout=subprocess.PIPE, stderr=subprocess.PIPE)
        p.stdout.close()
        _, err = p.communicate()
        self.assertEqual(p.returncode, 0, err)


if __name__ == "__main__":
    unittest.main()
