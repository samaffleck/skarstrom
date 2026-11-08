#!/usr/bin/env python3
"""Run the code blocks of the documentation against the installed package.

Every fenced ```python and ```bash / ```sh / ```console block in docs/*.md and
AGENTS.md is run in a fresh temporary directory (one per page), so a call that
the engine rejects (a renamed method, a removed option) fails CI instead of
misleading a reader. Python blocks run in a subprocess of the interpreter
running this script, shell blocks in bash with that interpreter's ``skarstrom``
on PATH.

Conventions (documented in AGENTS.md, "Documentation code blocks"):

  <!-- no-run -->    on the line before the fence (blank lines between are
                     fine): the block is shown but never run. Say why in the
                     comment, e.g. ``<!-- no-run: needs a GPU -->``; use
                     ``<!-- no-run: engine bug #123 -->`` for a known bug.
  ```python no-run   the same, as an info string word.
  <!-- long -->       the block is slow: it runs only with --full.
  <!-- doc-setup     a hidden Python block (the lines up to a closing --> line)
  ...                that runs first on its page and defines what the visible
  -->                fragments assume (``sim``, ``ads``, input data); it is not
                     rendered. Visible blocks that define everything they use
                     need none.
  ```text            (and any other language) is never run.

All Python blocks of a page share one interpreter namespace and one temporary
directory, in page order after the setup blocks, as a reader would type them
into one session; every shell block is its own bash process in that directory
after the Python blocks, so a deck the Python blocks save is there for the
command line. A page never sees another page's state.

Cycle cap: unless --full is given, every run of a deck is limited to --cycles
cycles (default 3). The runner puts a ``sitecustomize.py`` on PYTHONPATH that
clamps ``Simulation.run`` (which the CLI and the Python API both go through)
and sets ``SKARSTROM_DOC_RUN=1`` and ``SKARSTROM_DOC_CYCLES=<n>``, which a
block may read. A block that must not be capped sets nothing: mark it long.

Usage:
  python scripts/docs/run_doc_blocks.py [--full] [--cycles N] [--timeout S]
                                        [--jobs N] [--only SUBSTRING] [FILE ...]

Exit status 0 if every block that ran passed, 1 otherwise.
"""
from __future__ import annotations

import argparse
import json
import os
import re
import shutil
import subprocess
import sys
import tempfile
import textwrap
import time
from concurrent.futures import ThreadPoolExecutor
from dataclasses import dataclass
from pathlib import Path

ROOT = Path(__file__).resolve().parents[2]
PY_LANGS = {"python", "py"}
SH_LANGS = {"bash", "sh", "console", "shell"}
FENCE = re.compile(r"^(\s*)(`{3,}|~{3,})\s*([^\s`]*)\s*(.*)$")
MARKER = re.compile(r"^\s*<!--\s*(no-run|long)\b(.*?)-->\s*$")

SITECUSTOMIZE = '''\
"""Cycle cap for documentation runs (written by run_doc_blocks.py)."""
import os

_cap = int(os.environ.get("SKARSTROM_DOC_CYCLES", "0") or 0)
if _cap > 0:
    try:
        import skarstrom as _sk

        _run = _sk.Simulation.run

        def _capped(self, cycles=None, *args, **kwargs):
            n = cycles if cycles is not None else int(self.get("cycle.number_of_cycles"))
            return _run(self, min(int(n), _cap), *args, **kwargs)

        _capped.__doc__ = _run.__doc__
        _sk.Simulation.run = _capped
    except Exception as exc:  # the block then fails on its own import
        import sys
        print(f"doc runner: cycle cap not installed: {exc!r}", file=sys.stderr)
'''


@dataclass
class Block:
    file: Path
    line: int  # 1-based line of the opening fence
    lang: str  # "python" or "shell"
    code: str
    skip: str = ""  # reason, when the block is not run
    long: bool = False
    setup: bool = False  # a hidden <!-- doc-setup --> block: runs first, shares the namespace

    @property
    def label(self) -> str:
        return f"{self.file.relative_to(ROOT)}:{self.line}"


def extract(path: Path) -> list[Block]:
    blocks: list[Block] = []
    lines = path.read_text(encoding="utf-8").splitlines()
    pending: list[tuple[str, str]] = []  # markers seen since the last non-blank line
    i = 0
    while i < len(lines):
        if re.match(r"^\s*<!--\s*doc-setup\s*$", lines[i]):
            start, body = i, []
            i += 1
            while i < len(lines) and lines[i].strip() != "-->":
                body.append(lines[i])
                i += 1
            i += 1
            blocks.append(Block(path, start + 1, "python", "\n".join(body) + "\n", setup=True))
            pending = []
            continue
        m = FENCE.match(lines[i])
        if not m:
            mk = MARKER.match(lines[i])
            if mk:
                pending.append((mk.group(1), mk.group(2).strip(" :")))
            elif lines[i].strip():
                pending = []
            i += 1
            continue
        indent, fence, lang, rest = m.groups()
        start = i
        body: list[str] = []
        i += 1
        while i < len(lines) and not (lines[i].strip().startswith(fence) and lines[i].strip().strip(fence[0]) == ""):
            body.append(lines[i])
            i += 1
        i += 1  # closing fence
        lang_l = lang.lower()
        words = rest.lower().replace("{", " ").replace("}", " ").split()
        kind = "python" if lang_l in PY_LANGS else "shell" if lang_l in SH_LANGS else None
        if kind:
            code = textwrap.dedent("\n".join(l[len(indent):] if l.startswith(indent) else l for l in body)) + "\n"
            b = Block(path, start + 1, kind, code)
            for name, why in pending:
                if name == "no-run":
                    b.skip = why or "marked no-run"
                elif name == "long":
                    b.long = True
            if "no-run" in words:
                b.skip = b.skip or "marked no-run"
            blocks.append(b)
        pending = []
    return blocks


DRIVER = '''\
import json as _json, sys as _sys, time as _time, traceback as _tb, signal as _signal
_blocks = _json.load(open(_sys.argv[1]))
_timeout = int(_sys.argv[2])
_results = []

def _alarm(*_):
    raise TimeoutError(f"block exceeded {_timeout} s")

_signal.signal(_signal.SIGALRM, _alarm)
# The blocks run in this module's own namespace (``__main__``), so a function
# they define can be pickled for a process pool.
for _i, _code in enumerate(_blocks):
    _t0 = _time.monotonic()
    _signal.alarm(_timeout)
    try:
        exec(compile(_code, f"<block {_i}>", "exec"), globals())
        _st, _err = "PASS", ""
    except BaseException:
        _st, _err = "FAIL", _tb.format_exc()
    _signal.alarm(0)
    _results.append([_st, _time.monotonic() - _t0, _err])
    _json.dump(_results, open(_sys.argv[3], "w"))
'''


def _env(args, shim: Path, bindir: Path) -> dict:
    env = dict(os.environ)
    env["PYTHONPATH"] = str(shim) + os.pathsep + env.get("PYTHONPATH", "")
    env["PATH"] = str(bindir) + os.pathsep + env["PATH"]
    env["MPLBACKEND"] = "Agg"
    env["PYTHONIOENCODING"] = "utf-8"
    env["SKARSTROM_DOC_RUN"] = "1"
    env["SKARSTROM_DOC_REPO"] = str(ROOT)
    env["SKARSTROM_DOC_CYCLES"] = "0" if args.full else str(args.cycles)
    return env


def run_file(blocks: list[Block], args, shim: Path, bindir: Path) -> dict[int, tuple[str, float, str]]:
    """Run the runnable blocks of one file in one fresh temporary directory.

    Python blocks share one interpreter namespace in file order, as the page
    reads (a later block uses ``sim`` from an earlier one); each shell block is
    its own bash process in the same directory. Returns {index: (status,
    seconds, output)} with status PASS, FAIL or TIMEOUT.
    """
    env = _env(args, shim, bindir)
    timeout = args.timeout * (5 if args.full else 1)
    res: dict[int, tuple[str, float, str]] = {}
    with tempfile.TemporaryDirectory(prefix="docblock-") as tmp:
        py = [(n, b) for n, b in enumerate(blocks) if b.lang == "python"]
        py.sort(key=lambda nb: not nb[1].setup)  # setup blocks first, otherwise file order
        if py:
            (Path(tmp) / ".blocks.json").write_text(json.dumps([b.code for _, b in py]))
            (Path(tmp) / ".driver.py").write_text(DRIVER)
            try:
                p = subprocess.run([sys.executable, ".driver.py", ".blocks.json", str(timeout), ".results.json"],
                                   cwd=tmp, env=env, capture_output=True, text=True,
                                   timeout=timeout * len(py) + 30)
                crash = p.stderr[-2000:]
            except subprocess.TimeoutExpired:
                crash = "the interpreter did not finish"
            try:
                done = json.loads((Path(tmp) / ".results.json").read_text())
            except (OSError, ValueError):
                done = []
            for k, (n, b) in enumerate(py):
                st, sec, err = done[k] if k < len(done) else ("FAIL", 0.0, "not run: " + crash)
                res[n] = ("TIMEOUT" if "TimeoutError" in err else st, sec, err)
        for n, b in enumerate(blocks):
            if b.lang != "shell":
                continue
            t0 = time.monotonic()
            (Path(tmp) / ".block.sh").write_text("set -eo pipefail\n" + b.code, encoding="utf-8")
            try:
                p = subprocess.run(["bash", ".block.sh"], cwd=tmp, env=env, capture_output=True, text=True,
                                   timeout=timeout)
                st, out = ("PASS" if p.returncode == 0 else "FAIL"), p.stdout + p.stderr
            except subprocess.TimeoutExpired as e:
                st, out = "TIMEOUT", ""
            res[n] = (st, time.monotonic() - t0, out if st != "PASS" else "")
    return res


def main() -> int:
    ap = argparse.ArgumentParser(description=__doc__.split("\n\n")[0])
    ap.add_argument("files", nargs="*", help="Markdown files (default: docs/*.md and AGENTS.md)")
    ap.add_argument("--full", action="store_true", help="no cycle cap, include blocks marked long (nightly)")
    ap.add_argument("--cycles", type=int, default=3, help="cycle cap per run (default 3)")
    ap.add_argument("--timeout", type=int, default=120, help="seconds per block (x5 with --full)")
    ap.add_argument("--jobs", type=int, default=2, help="files run in parallel (default 2)")
    ap.add_argument("--only", help="run only blocks whose 'file:line' contains this text")
    ap.add_argument("--list", action="store_true", help="list the blocks and exit")
    args = ap.parse_args()

    files = [Path(f).resolve() for f in args.files] or sorted((ROOT / "docs").glob("*.md")) + [ROOT / "AGENTS.md"]
    blocks = [b for f in files for b in extract(f)]
    if args.only:
        blocks = [b for b in blocks if args.only in b.label]
    if args.list:
        for b in blocks:
            print(f"{b.label}  {b.lang}  {'skip: ' + b.skip if b.skip else 'long' if b.long else 'run'}")
        return 0

    try:
        import skarstrom  # noqa: F401
    except ImportError:
        print("skarstrom is not installed in this interpreter; install the built wheel first.", file=sys.stderr)
        return 2
    bindir = Path(sys.executable).parent
    if not (bindir / "skarstrom").exists() and not shutil.which("skarstrom"):
        print("the skarstrom command is not on PATH", file=sys.stderr)
        return 2

    with tempfile.TemporaryDirectory(prefix="docshim-") as shim_dir:
        shim = Path(shim_dir)
        (shim / "sitecustomize.py").write_text(SITECUSTOMIZE, encoding="utf-8")
        todo = []
        results: dict[int, tuple[str, float, str]] = {}
        by_file: dict[Path, list[int]] = {}
        for n, b in enumerate(blocks):
            if b.skip:
                results[n] = ("SKIP", 0.0, b.skip)
            elif b.long and not args.full:
                results[n] = ("SKIP", 0.0, "long (runs with --full)")
            else:
                todo.append(n)
                by_file.setdefault(b.file, []).append(n)
        t0 = time.monotonic()

        def job(f: Path):
            idx = by_file[f]
            r = run_file([blocks[n] for n in idx], args, shim, bindir)
            return {n: r[k] for k, n in enumerate(idx)}

        with ThreadPoolExecutor(max_workers=max(1, args.jobs)) as pool:
            for r in pool.map(job, list(by_file)):
                results.update(r)
        total = time.monotonic() - t0
        for n in todo:
            st, sec, _ = results[n]
            print(f"{st:7} {sec:6.1f}s  {blocks[n].label}  [{blocks[n].lang}]")

    failed = [n for n in todo if results[n][0] != "PASS"]
    for n in failed:
        st, _, out = results[n]
        print(f"\n--- {st}: {blocks[n].label} ---")
        print(textwrap.indent(blocks[n].code.rstrip(), "  | "))
        print("  output:")
        print(textwrap.indent("\n".join(out.rstrip().splitlines()[-25:]), "  > "))
    skipped = [n for n in results if results[n][0] == "SKIP"]
    print(f"\n{len(todo) - len(failed)} passed, {len(failed)} failed, {len(skipped)} skipped "
          f"of {len(blocks)} blocks in {total:.0f} s"
          + ("" if args.full else f" (cycle cap {args.cycles}; --full for the uncapped run)"))
    for n in skipped:
        print(f"  skipped {blocks[n].label}: {results[n][2]}")
    return 1 if failed else 0


if __name__ == "__main__":
    sys.exit(main())
