"""The ``skarstrom`` command-line tool.

    skarstrom run deck.json [--cycles N] [--set PATH=VALUE ...] [--csv DIR] [--json]
    skarstrom params deck.json [--all] [--search TEXT] [--json]
    skarstrom examples [list] [--json]
    skarstrom examples run NAME [run options]
    skarstrom examples path NAME
    skarstrom examples copy NAME [DEST]
    skarstrom --version
"""

from __future__ import annotations

import argparse
import csv
import json
import math
import os
import re
import sys
from typing import List, Optional, Sequence, Tuple

from . import Simulation, __version__, examples


def _fmt(x) -> str:
    if isinstance(x, (bool, str)):
        return str(x).lower() if isinstance(x, bool) else x
    if isinstance(x, float) and math.isnan(x):
        return "n/a"
    if isinstance(x, int) or x == 0:
        return str(x)
    if 1e-3 <= abs(x) < 1e5:
        return f"{x:.6g}"
    return f"{x:.4e}"


def _resolve_deck(path: str) -> str:
    if os.path.isfile(path):
        return path
    try:
        return examples.get(path).path
    except KeyError:
        pass
    raise SystemExit(f"skarstrom: no such deck: {path} (not a file, nor a bundled example; see 'skarstrom examples')")


def _split_assignment(sim: Simulation, text: str) -> Tuple[str, str]:
    """Split PATH=VALUE at the first '=' whose left side names a parameter.

    Paths may contain '=' (display names such as "Q=100"), so the split point
    is found by asking the deck which prefix is a parameter.
    """
    positions = [i for i, c in enumerate(text) if c == "="]
    if not positions:
        raise SystemExit(f"skarstrom: --set expects PATH=VALUE, got {text!r}")
    for i in positions:
        path = text[:i].strip()
        try:
            sim.get_text(path)
        except (KeyError, ValueError, RuntimeError):
            continue
        return path, text[i + 1 :].strip()
    raise SystemExit(
        f"skarstrom: unknown parameter in {text!r}; list them with 'skarstrom params DECK'"
    )


def _apply_sets(sim: Simulation, assignments: Sequence[str]) -> None:
    for a in assignments:
        path, value = _split_assignment(sim, a)
        try:
            number = float(value)
        except ValueError:
            number = None
        try:
            if number is not None:
                sim.set(path, number)
            else:
                sim.set(path, value)
        except (ValueError, RuntimeError) as e:
            raise SystemExit(f"skarstrom: cannot set {path} = {value!r}: {e}")


def _print_table(rows: List[List[str]], header: List[str], out=None) -> None:
    out = out or sys.stdout  # resolved per call, so redirected stdout is honoured
    widths = [len(h) for h in header]
    for r in rows:
        widths = [max(w, len(c)) for w, c in zip(widths, r)]
    line = "  ".join(h.ljust(w) for h, w in zip(header, widths))
    print(line, file=out)
    print("  ".join("-" * w for w in widths), file=out)
    for r in rows:
        print("  ".join(c.ljust(w) for c, w in zip(r, widths)), file=out)


KPI_FIELDS = [
    ("purity", "purity [mol/mol]"),
    ("recovery", "recovery [mol/mol]"),
    ("productivity_mass", "productivity [kg/(kg s)]"),
    ("productivity_volume", "productivity [kg/(m3 s)]"),
    ("energy_consumption", "energy [J/kg]"),
    ("production_rate", "production [mol/s]"),
]


def cmd_run(args: argparse.Namespace) -> int:
    sim = Simulation.from_file(_resolve_deck(args.deck))
    _apply_sets(sim, args.set or [])
    if args.cycles is not None:
        sim.set("cycle.number_of_cycles", args.cycles)

    def progress(p):
        if not args.quiet:
            print(
                f"\r  cycle {p['cycle'] + 1}/{p['total_cycles']}  step {p['step'] + 1}/{p['total_steps']}",
                end="",
                file=sys.stderr,
                flush=True,
            )
        return True

    def on_message(level: str, text: str) -> None:
        print(f"[{level}] {text}", file=sys.stderr)

    try:
        result = sim.run(
            progress=progress,
            log_level=args.log or "low",
            on_message=on_message if args.log else None,
        )
    except KeyboardInterrupt:
        if not args.quiet:
            print("", file=sys.stderr)
        print("skarstrom: interrupted", file=sys.stderr)
        return 130
    if not args.quiet:
        print("", file=sys.stderr)

    kpis = sim.kpis() if sim.has_results else []
    css = sim.css_history() if sim.has_results else []
    last = sim.last_cycle

    if args.json:
        payload = {
            "deck": args.deck,
            "status": result.status,
            "message": result.message,
            "css_reached": result.css_reached,
            "css_cycle": result.css_cycle,
            "cycles_run": result.cycles_run,
            "wall_time_s": result.wall_time_s,
            "stats": result.stats,
            "kpis": [{k: (None if isinstance(v, float) and math.isnan(v) else v) for k, v in r.as_dict().items()}
                     for r in kpis],
            "css": [{"cycle": c.cycle, "converged": c.converged, "error": c.error, "tolerance": c.tolerance,
                     "method": c.method} for c in css],
        }
        print(json.dumps(payload, indent=2))
    else:
        print(f"Deck:    {args.deck}")
        print(f"Status:  {result.status} ({result.message})")
        print(f"Cycles:  {result.cycles_run} run in {result.wall_time_s:.2f} s "
              f"({result.stats['steps']} integrator steps)")
        if result.css_reached:
            print(f"CSS:     reached at cycle {result.css_cycle + 1}")
        elif css and css[-1].method != "None":
            c = css[-1]
            min_cycles = sim.get("css.min_cycles")
            if result.cycles_run < min_cycles:
                print(f"CSS:     not reached ({c.method}; the check starts at cycle {min_cycles}, "
                      f"css.min_cycles)")
            else:
                print(f"CSS:     not reached ({c.method}: error {c.error:.3e}, tolerance {c.tolerance:.3e})")
        else:
            print("CSS:     not checked (the deck's CSS method is None)")
        final = [k for k in kpis if k.cycle == last]
        if final:
            print(f"\nKPIs at cycle {last + 1}:")
            rows = []
            for k in final:
                rows.append([k.component] + ([_fmt(getattr(k, f)) for f, _ in KPI_FIELDS] if k.valid
                                             else ["(no feed/product endpoints)"] + [""] * (len(KPI_FIELDS) - 1)))
            _print_table(rows, ["component"] + [h for _, h in KPI_FIELDS])
        elif sim.has_results:
            print("\nKPIs:    none (no component is tracked in the deck's KPI configuration)")

    if args.csv:
        os.makedirs(args.csv, exist_ok=True)
        with open(os.path.join(args.csv, "kpis.csv"), "w", newline="") as f:
            w = csv.writer(f)
            w.writerow(["cycle", "component", "valid"] + [h for _, h in KPI_FIELDS])
            for k in kpis:
                w.writerow([k.cycle + 1, k.component, k.valid] + [getattr(k, fld) for fld, _ in KPI_FIELDS])
        with open(os.path.join(args.csv, "css.csv"), "w", newline="") as f:
            w = csv.writer(f)
            w.writerow(["cycle", "method", "converged", "error", "tolerance"])
            for c in css:
                w.writerow([c.cycle + 1, c.method, c.converged, c.error, c.tolerance])
        with open(os.path.join(args.csv, "run.json"), "w") as f:
            json.dump({"status": result.status, "message": result.message, "css_reached": result.css_reached,
                       "css_cycle": result.css_cycle, "cycles_run": result.cycles_run,
                       "wall_time_s": result.wall_time_s, "stats": result.stats}, f, indent=2)
        if sim.has_results:
            sim.save_logs(os.path.join(args.csv, "logs"))
        if not args.json:
            print(f"\nWrote kpis.csv, css.csv, run.json and logs/ to {args.csv}")

    return 0 if result.ok else 1


def cmd_params(args: argparse.Namespace) -> int:
    sim = Simulation.from_file(_resolve_deck(args.deck))
    params = sim.parameters(tunable_only=not args.all, search=args.search)
    if args.json:
        out = []
        for p in params:
            d = dict(p.__dict__)
            for k in ("lower", "upper", "value"):
                if isinstance(d[k], float) and math.isnan(d[k]):
                    d[k] = None
            out.append(d)
        print(json.dumps(out, indent=2))
        return 0
    rows = []
    for p in params:
        if p.kind == "string":
            value = repr(p.text)
        elif p.kind == "expression" and p.text != _fmt(p.value):
            value = f"{p.text} (= {_fmt(p.value)})"  # a formula: show both
        else:
            value = _fmt(p.value)
        b = p.bounds
        rows.append([p.path, value, p.unit, f"{_fmt(b[0])} .. {_fmt(b[1])}" if b else "", p.alias if p.alias != p.path else ""])
    _print_table(rows, ["path", "value", "unit", "bounds", "alias"])
    if not args.all:
        print(f"\n{len(rows)} tunable parameters. Use --all to include solver settings, integers and text.")
    return 0


def _example_family(name: str) -> str:
    """The paper / study an example deck belongs to: its AuthorYear prefix."""
    m = re.match(r"[A-Za-z]+\d{4}", name)
    return m.group(0) if m else name.split("_")[0]


def _example(name: str) -> "examples.Example":
    try:
        return examples.get(name)
    except KeyError as e:
        raise SystemExit(f"skarstrom: {e.args[0]}") from None


def cmd_examples(args: argparse.Namespace) -> int:
    action = getattr(args, "action", None) or "list"
    if action == "list":
        names = examples.names()
        if args.json:
            print(json.dumps([{"name": e.name, "title": e.title, "family": _example_family(e.name), "path": e.path}
                              for e in examples.ALL], indent=2))
            return 0
        for name in names:
            print(name)
        if not names:
            print("skarstrom: no example decks are installed with this package", file=sys.stderr)
            return 1
        return 0
    if action == "path":
        print(_example(args.name).path)
        return 0
    if action == "copy":
        print(_example(args.name).copy_to(args.dest))
        return 0
    if action == "run":
        args.deck = args.name
        return cmd_run(args)
    raise SystemExit(f"skarstrom: unknown examples action {action!r}")


def _add_run_arguments(r: argparse.ArgumentParser) -> None:
    r.add_argument("--cycles", type=int, help="override the number of cycles")
    r.add_argument("--set", action="append", metavar="PATH=VALUE",
                   help="set a parameter before running (repeatable); VALUE may be an expression such as '5 [bar]'")
    r.add_argument("--csv", metavar="DIR", help="write kpis.csv, css.csv, run.json and the solver logs to DIR")
    r.add_argument("--json", action="store_true", help="print the run summary as JSON")
    r.add_argument("--log", choices=["low", "medium", "high"], help="print solver messages to stderr")
    r.add_argument("-q", "--quiet", action="store_true", help="no progress line")


def build_parser() -> argparse.ArgumentParser:
    ap = argparse.ArgumentParser(prog="skarstrom", description="Skarstrom PSA/TSA cycle simulator")
    ap.add_argument("--version", action="version", version=f"skarstrom {__version__}")
    sub = ap.add_subparsers(dest="command", required=True)

    r = sub.add_parser("run", help="run a deck and print its KPIs and CSS status")
    r.add_argument("deck", help="deck file (or the name of a shipped example, see 'skarstrom examples')")
    _add_run_arguments(r)
    r.set_defaults(func=cmd_run)

    p = sub.add_parser("params", help="list a deck's settable parameters")
    p.add_argument("deck")
    p.add_argument("--all", action="store_true", help="include solver settings, integers and text fields")
    p.add_argument("--search", metavar="TEXT", help="only entries whose path, alias or name contains TEXT")
    p.add_argument("--json", action="store_true", help="print as JSON")
    p.set_defaults(func=cmd_params)

    e = sub.add_parser("examples", help="list, run or copy the example decks shipped with the package")
    e.add_argument("--json", action="store_true", help="list as JSON (name, family, path)")
    e.set_defaults(func=cmd_examples, action="list")
    esub = e.add_subparsers(dest="action")
    el = esub.add_parser("list", help="list the bundled decks (the default)")
    el.add_argument("--json", action="store_true", help="list as JSON (name, family, path)")
    er = esub.add_parser("run", help="run a bundled deck (same options as 'skarstrom run')")
    er.add_argument("name", help="deck name, with or without .json (see 'skarstrom examples')")
    _add_run_arguments(er)
    ep = esub.add_parser("path", help="print the absolute path of a bundled deck")
    ep.add_argument("name")
    ec = esub.add_parser("copy", help="copy a bundled deck to DEST (default: the current directory)")
    ec.add_argument("name")
    ec.add_argument("dest", nargs="?", default=".")
    return ap


def main(argv: Optional[Sequence[str]] = None) -> int:
    args = build_parser().parse_args(argv)
    try:
        return args.func(args)
    except BrokenPipeError:  # e.g. piped into `head`
        try:
            sys.stdout = open(os.devnull, "w")
        except OSError:
            pass
        return 0
    except (RuntimeError, ValueError, OSError) as e:
        print(f"skarstrom: {e}", file=sys.stderr)
        return 1


if __name__ == "__main__":
    sys.exit(main())
