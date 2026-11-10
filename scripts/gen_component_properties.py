#!/usr/bin/env python3
# Skarstrom - pressure and temperature swing adsorption simulator
# Copyright (c) 2024-2026 Sam Affleck and contributors
# Licensed under the terms in the LICENSE file at the repository root.
"""Generate core/include/adsim/objects/component_properties.h.

The engine reads a few fixed per-component constants: the formula, CAS number,
aliases and molar mass (the component manifest, the molar mass of the gas) and
the critical point, triple point and vapour-pressure ancillary (the saturation
pressure of the GAB isotherms). This script reads them once from the CoolProp
fluid files (dev/fluids/*.json) and writes them as a C++ table, so the engine
does not call CoolProp at run time.

The component list is Components::componentNames in
core/include/adsim/objects/fluid.h. Each name is resolved to a fluid file by
the file's INFO.NAME or one of its INFO.ALIASES, as CoolProp resolves it.

The values are the ones CoolProp reads for the same calls:
  molar mass      EOS[0].molar_mass                     (Props1SI "molar_mass")
  T_c, p_c        STATES.critical.T / .p                (T_critical, p_critical)
  T_triple        EOS[0].STATES.sat_min_liquid.T        (Ttriple)
  p_sat ancillary ANCILLARIES.pS, else ANCILLARIES.pV   (saturation_ancillary(iP, 1, iT, T))
  formula         INFO.FORMULA with the "_{n}" markup removed ("N/A" when absent)

Usage:
    python scripts/gen_component_properties.py                # write the header
    python scripts/gen_component_properties.py --check        # fail if the committed header differs
    python scripts/gen_component_properties.py --fluids DIR   # another dev/fluids directory
"""
from __future__ import annotations

import argparse
import difflib
import json
import re
import sys
from pathlib import Path

ROOT = Path(__file__).resolve().parents[1]
FLUID_H = ROOT / "core/include/adsim/objects/fluid.h"
OUT = ROOT / "core/include/adsim/objects/component_properties.h"
COOLPROP = ROOT / "external/CoolProp"
DEFAULT_FLUIDS = COOLPROP / "dev/fluids"

# Most terms of any vapour-pressure ancillary in the table (std::array size).
MAX_TERMS = 6


def component_names() -> list[str]:
    """Components::componentNames, read from fluid.h."""
    text = FLUID_H.read_text(encoding="utf-8")
    m = re.search(r"componentNames\s*=\s*\{(.*?)\};", text, re.S)
    if not m:
        sys.exit(f"cannot find Components::componentNames in {FLUID_H}")
    return re.findall(r'"([^"]+)"', m.group(1))


def coolprop_version() -> str:
    cmake = COOLPROP / "CMakeLists.txt"
    if not cmake.exists():
        return "unknown"
    text = cmake.read_text(encoding="utf-8")
    parts = [re.search(rf"set\(COOLPROP_VERSION_{k} ([^)]*)\)", text) for k in ("MAJOR", "MINOR", "PATCH")]
    if not all(parts):
        return "unknown"
    return ".".join(p.group(1).strip() for p in parts)


def plain_formula(f: str) -> str:
    """'C_{1}O_{2}' -> 'CO2' (drop the subscript markup and the ones)."""
    return re.sub(r"_\{(\d+)\}", lambda m: "" if m.group(1) == "1" else m.group(1), f).replace("{", "").replace("}", "")


def num(x: float) -> str:
    """Shortest text that round-trips to the same double."""
    r = repr(float(x))
    return r if any(c in r for c in ".einf") else r + ".0"


def cstr(s: str) -> str:
    return '"' + s.replace("\\", "\\\\").replace('"', '\\"') + '"'


def load_index(fluids: Path) -> dict[str, Path]:
    index: dict[str, Path] = {}
    for f in sorted(fluids.glob("*.json")):
        info = json.loads(f.read_text(encoding="utf-8"))["INFO"]
        for name in [info["NAME"], *info.get("ALIASES", [])]:
            index.setdefault(name, f)
    return index


def row(name: str, path: Path) -> str:
    d = json.loads(path.read_text(encoding="utf-8"))
    info, eos, anc = d["INFO"], d["EOS"][0], d["ANCILLARIES"]
    if eos.get("molar_mass_units", "kg/mol") != "kg/mol":
        sys.exit(f"{path.name}: molar mass is not in kg/mol")
    crit = d["STATES"]["critical"]
    formula = plain_formula(info["FORMULA"]) if "FORMULA" in info else "N/A"
    aliases = ",".join(info.get("ALIASES", []))
    key = "pS" if "pS" in anc else ("pV" if "pV" in anc else None)
    lines = [
        f"    // {path.name}" + (f", ANCILLARIES.{key}" if key else ", no vapour-pressure ancillary"),
        "    {" + f"{cstr(name)}, {cstr(info['NAME'])}, {cstr(formula)}, {cstr(info['CAS'])},",
        f"     {cstr(aliases)},",
        f"     {num(eos['molar_mass'])}, {num(crit['T'])}, {num(crit['p'])}, "
        f"{num(eos['STATES']['sat_min_liquid']['T'])},",
    ]
    if key:
        a = anc[key]
        if a["type"] not in ("pL", "pV"):
            sys.exit(f"{path.name}: ancillary {key} has type {a['type']}, not an exponential pL/pV form")
        n, t = a["n"], a["t"]
        if len(n) != len(t) or len(n) > MAX_TERMS:
            sys.exit(f"{path.name}: ancillary {key} has {len(n)} terms (at most {MAX_TERMS})")
        pad = [0.0] * (MAX_TERMS - len(n))
        lines.append(
            f"     {{true, {len(n)}, {{{{{', '.join(num(v) for v in n + pad)}}}}},\n"
            f"      {{{{{', '.join(num(v) for v in t + pad)}}}}},\n"
            f"      {num(a['T_r'])}, {num(a['reducing_value'])}, {'true' if a['using_tau_r'] else 'false'}, "
            f"{num(a['Tmin'])}, {num(a['Tmax'])}}}}},"
        )
    else:
        lines.append("     {}},")
    return "\n".join(lines)


def generate(fluids: Path) -> str:
    index = load_index(fluids)
    names = component_names()
    missing = [n for n in names if n not in index]
    if missing:
        sys.exit(f"no CoolProp fluid file for: {', '.join(missing)}")
    rows = "\n".join(row(n, index[n]) for n in names)
    return f"""/*
 * Skarstrom - pressure and temperature swing adsorption simulator
 * Copyright (c) 2024-2026 Sam Affleck and contributors
 *
 * Licensed under the terms in the LICENSE file at the repository root.
 */
// GENERATED by scripts/gen_component_properties.py. Do not edit.
#pragma once

#include <array>
#include <cmath>
#include <cstddef>
#include <numeric>
#include <string_view>

// ---------------------------------------------------------------------------
// Fixed per-component constants of every component in
// Components::componentNames (fluid.h).
//
// Source: the CoolProp {coolprop_version()} fluid files (external/CoolProp/dev/fluids/*.json,
// one file per row, named in the comment above the row). CoolProp is a
// development-time data source only: run the generator again to update.
//   molar_mass        EOS[0].molar_mass [kg/mol]
//   T_critical        STATES.critical.T [K]
//   p_critical        STATES.critical.p [Pa]
//   T_triple          EOS[0].STATES.sat_min_liquid.T [K] (CoolProp's Ttriple)
//   p_sat             the vapour-pressure ancillary ANCILLARIES.pS (pV when
//                     the file has separate pL and pV; Air):
//                       p = reducing_value * exp((T_r / T) * sum_i n_i theta^t_i),
//                       theta = 1 - T / T_r  (T_r / T -> 1 when !using_tau_r)
//   formula           INFO.FORMULA without the subscript markup ("CO2")
// ---------------------------------------------------------------------------

namespace adsim {{

  struct VaporPressureAncillary {{
    bool valid = false;
    std::size_t terms = 0;
    std::array<double, {MAX_TERMS}> n{{}};
    std::array<double, {MAX_TERMS}> t{{}};
    double T_r = 0.0;             // [K]
    double reducing_value = 0.0;  // [Pa]
    bool using_tau_r = false;
    double T_min = 0.0, T_max = 0.0;  // [K] fit range of the ancillary (provenance)

    // Vapour pressure [Pa] at T [K], evaluated as CoolProp evaluates it
    // (SaturationAncillaryFunction::evaluate), term by term in the same order.
    double Evaluate(double T) const {{
      const double theta = 1.0 - T / T_r;
      std::array<double, {MAX_TERMS}> s{{}};
      for (std::size_t i = 0; i < terms; ++i) s[i] = n[i] * std::pow(theta, t[i]);
      const double sum = std::accumulate(s.begin(), s.begin() + static_cast<std::ptrdiff_t>(terms), 0.0);
      const double tau_r = using_tau_r ? T_r / T : 1.0;
      return reducing_value * std::exp(tau_r * sum);
    }}
  }};

  struct ComponentProperties {{
    std::string_view name;           // engine (CoolProp) name, as decks write it
    std::string_view coolprop_name;  // INFO.NAME of the CoolProp fluid file
    std::string_view formula;        // plain formula, "CO2" ("N/A" when none)
    std::string_view cas;
    std::string_view aliases;        // CoolProp's aliases, comma-separated
    double molar_mass;               // [kg/mol]
    double T_critical;               // [K]
    double p_critical;               // [Pa]
    double T_triple;                 // [K]
    VaporPressureAncillary p_sat;
  }};

  inline constexpr const char* kComponentPropertiesSource =
      "CoolProp {coolprop_version()} fluid files (generated table)";

  inline const std::array<ComponentProperties, {len(names)}> kComponentProperties = {{{{
{rows}
  }}}};

  // The row of `name`, or null for a component the table does not have.
  inline const ComponentProperties* FindComponentProperties(std::string_view name) {{
    for (const auto& c : kComponentProperties)
      if (c.name == name) return &c;
    return nullptr;
  }}

}} // namespace adsim
"""


def main() -> int:
    ap = argparse.ArgumentParser(description=__doc__.splitlines()[0])
    ap.add_argument("--check", action="store_true", help="fail if the committed header differs from the output")
    ap.add_argument("--fluids", type=Path, default=DEFAULT_FLUIDS, help="CoolProp dev/fluids directory")
    args = ap.parse_args()
    if not args.fluids.is_dir():
        sys.exit(f"{args.fluids} is not a directory (run: git submodule update --init external/CoolProp)")
    text = generate(args.fluids)
    if args.check:
        old = OUT.read_text(encoding="utf-8") if OUT.exists() else ""
        if old != text:
            sys.stdout.writelines(difflib.unified_diff(old.splitlines(True), text.splitlines(True),
                                                       str(OUT.relative_to(ROOT)), "generated"))
            print(f"\n{OUT.relative_to(ROOT)} is out of date: run python scripts/gen_component_properties.py",
                  file=sys.stderr)
            return 1
        print(f"{OUT.relative_to(ROOT)} is up to date")
        return 0
    OUT.write_text(text, encoding="utf-8")
    print(f"wrote {OUT.relative_to(ROOT)}")
    return 0


if __name__ == "__main__":
    sys.exit(main())
