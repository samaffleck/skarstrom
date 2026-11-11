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
does not use CoolProp at all: CoolProp is a development-time data source.

The fluid files are pinned: FLUID_FILES names the file of each component (the
file whose INFO.NAME or INFO.ALIASES has the name, as CoolProp resolves it)
and its SHA-256, at the CoolProp commit COOLPROP_COMMIT. By default the script
downloads them from raw.githubusercontent.com into a cache directory; --fluids
reads them from a local dev/fluids directory instead (a CoolProp checkout at
that commit). Either way a file whose SHA-256 differs is an error, so the
output depends only on the pinned files. To move to another CoolProp version,
change COOLPROP_VERSION, COOLPROP_COMMIT and the hashes, and regenerate.

The component list is Components::componentNames in
core/include/adsim/objects/fluid.h. A component without a FLUID_FILES entry is
an error.

The values are the ones CoolProp reads for the same calls:
  molar mass      EOS[0].molar_mass                     (Props1SI "molar_mass")
  T_c, p_c        STATES.critical.T / .p                (T_critical, p_critical)
  T_triple        EOS[0].STATES.sat_min_liquid.T        (Ttriple)
  p_sat ancillary ANCILLARIES.pS, else ANCILLARIES.pV   (saturation_ancillary(iP, 1, iT, T))
  formula         INFO.FORMULA with the "_{n}" markup removed ("N/A" when absent)

Usage:
    python scripts/gen_component_properties.py                # write the header
    python scripts/gen_component_properties.py --check        # fail if the committed header differs
    python scripts/gen_component_properties.py --fluids DIR   # local dev/fluids directory (no download)
"""
from __future__ import annotations

import argparse
import difflib
import hashlib
import json
import os
import re
import sys
import tempfile
import urllib.request
from pathlib import Path

ROOT = Path(__file__).resolve().parents[1]
FLUID_H = ROOT / "core/include/adsim/objects/fluid.h"
OUT = ROOT / "core/include/adsim/objects/component_properties.h"

COOLPROP_VERSION = "6.7.0"
COOLPROP_COMMIT = "e84ec7a8ab16a23e7fa8d498bc04db937bd93ab3"
URL = "https://raw.githubusercontent.com/CoolProp/CoolProp/{commit}/dev/fluids/{file}"

# Component -> (CoolProp fluid file, SHA-256 of the file at COOLPROP_COMMIT).
FLUID_FILES = {
    "Air": ("Air.json", "c5e0d064d83fc91ed35d7354cea1fa38fea832a09d2e85faf4336e0017c33616"),
    "Ammonia": ("Ammonia.json", "f7f4ce9a4227caca6bde4195b3ec484a9f9d084e6ce9035d3a67a756aebff312"),
    "Argon": ("Argon.json", "f52c19ef34818fcbae0f29444c8bfa6bc444e7c87bbaa58de08f392617315b0d"),
    "Benzene": ("Benzene.json", "0482070857b11b7ef87b2f50343f9952ffede0a98a037acac9eb55265e8521d8"),
    "Butane": ("n-Butane.json", "4af56f624fbf3bf65783e49937f6b35257be7eddbf7611c475bb0fc6720fe400"),
    "CarbonDioxide": ("CarbonDioxide.json", "2e56e1bdaa884c6ff0c18dc246ef9adfd7f4e0fceddc93f552dc5289483d847d"),
    "CarbonMonoxide": ("CarbonMonoxide.json", "41e59afc15a68db2ff377e8476dda7bdfd257dc13aa7e82251b800cabdf57f40"),
    "Cyclohexane": ("CycloHexane.json", "8b60f402a79f52734838f554324860a5f2d7da789c471cdc3609372edd3e181b"),
    "Ethane": ("Ethane.json", "b8d1d30be352688c2c40b043f9fd54859d27e2252a283229232c8fefea1004dc"),
    "Ethanol": ("Ethanol.json", "10f7b661c63e6a84eaa4fc6eeac21b4a417c76761df9ec0929af23f8b01aa512"),
    "Ethylene": ("Ethylene.json", "35fc2fa44976e10d48608c6cd1162de0b4da55ce51c3e4302b8866013c9f4a1b"),
    "Fluorine": ("Fluorine.json", "e77966c033aa420e799d2a29d184785fd00d77b245df1b20b6ba4102a6b583df"),
    "Helium": ("Helium.json", "ae864f5eb2d359926b2d4b4c158740e16ec611e034d4cbc33d13f33f4ae8cc51"),
    "Hydrogen": ("Hydrogen.json", "c1ee7e5f1dc86f14ba77b38fd1378d017475a7e39833a559b45f16584e1db9a0"),
    "HydrogenSulfide": ("HydrogenSulfide.json", "e776b1cc7937049b783dea9812a4b7cd619efa8788a2ad7da615a39ba0ee3aa9"),
    "Isobutane": ("IsoButane.json", "b946fc427c7ca3fd6e7c46aae388e0292e04369da071289ac9200d0e469a2c7e"),
    "Isopentane": ("Isopentane.json", "9b3e601afe29f0ff76ad834e6dae1f86077f08f6c06a0906721f0b4aba11aef1"),
    "Krypton": ("Krypton.json", "33efe15f9e68ae736798943a838faa50aa0049687f3ed31b9bfbbba93d98e353"),
    "Methane": ("Methane.json", "edd228248de4a610592483cc917d03267198780b3574a3092124817a275c4f45"),
    "Methanol": ("Methanol.json", "b746c99df235794267a199c2c205e8e6ec45d7d8907ecb1e54490c1ded533f10"),
    "Neon": ("Neon.json", "50a2f035e0b03a915840541d95f1c2f766cb0a32c9084503eee35db493952990"),
    "Nitrogen": ("Nitrogen.json", "8e9ad6c9f2c009541ddb00607f537ee66e642ae785aa73729b092115d5b71695"),
    "NitrousOxide": ("NitrousOxide.json", "165947802bd67d14ca063d84caf9f018275c301a008248fb3553844ca66dd439"),
    "Octane": ("n-Octane.json", "ba2d41e834a2d425ff6f372f015a15a02706ab92fbb88019e6853a2bc613781d"),
    "Oxygen": ("Oxygen.json", "a6f0c406add67e4cfd36caefc8f3fe178636c050e9d36e854cb7be692a0f0f35"),
    "Pentane": ("n-Pentane.json", "36a4858b2d5c570da373baee0619341973c994cd6d2228e97eb50403b8a8e3ab"),
    "Propane": ("n-Propane.json", "12942833c8e39820851d7f45d80356a4ce0f8102515008630e2cc25090fb7381"),
    "Propylene": ("Propylene.json", "f1265ea365bfd389128f2c673bce2124cc49aa94c8187962d425a38e4913ad01"),
    "SulfurDioxide": ("SulfurDioxide.json", "05a26680ef172bf8a3b6b6f4f4e5529a9f1e089e5b77d77487112b8cbdb0d371"),
    "SulfurHexafluoride": ("SulfurHexafluoride.json", "c3f86354900ae1377e9f56bfa240ef63368361368fa73f4cf93cccd790cee1c1"),
    "Toluene": ("Toluene.json", "b58b8925818c835141b4eb13e3dda88b971483d65da22f50efabf8ed4f5ac261"),
    "Water": ("Water.json", "a841fdcf36f7d5bc6b5e03ec6f94a552d89ac2881da3fe03cbd32ffaa8077863"),
    "Xenon": ("Xenon.json", "7b959f5bb55b70bd95977656e069bb76aad5316d8fd0d3177a270e037ff6589f"),
}

# Most terms of any vapour-pressure ancillary in the table (std::array size).
MAX_TERMS = 6


def component_names() -> list[str]:
    """Components::componentNames, read from fluid.h."""
    text = FLUID_H.read_text(encoding="utf-8")
    m = re.search(r"componentNames\s*=\s*\{(.*?)\};", text, re.S)
    if not m:
        sys.exit(f"cannot find Components::componentNames in {FLUID_H}")
    return re.findall(r'"([^"]+)"', m.group(1))


def plain_formula(f: str) -> str:
    """'C_{1}O_{2}' -> 'CO2' (drop the subscript markup and the ones)."""
    return re.sub(r"_\{(\d+)\}", lambda m: "" if m.group(1) == "1" else m.group(1), f).replace("{", "").replace("}", "")


def num(x: float) -> str:
    """Shortest text that round-trips to the same double."""
    r = repr(float(x))
    return r if any(c in r for c in ".einf") else r + ".0"


def cstr(s: str) -> str:
    return '"' + s.replace("\\", "\\\\").replace('"', '\\"') + '"'


def fluid_files(fluids: Path | None, names: list[str]) -> dict[str, Path]:
    """The pinned fluid file of each component: from `fluids`, or downloaded
    into a cache directory. Each file's SHA-256 is checked."""
    missing = [n for n in names if n not in FLUID_FILES]
    if missing:
        sys.exit(f"no pinned CoolProp fluid file for: {', '.join(missing)} (add it to FLUID_FILES)")
    if fluids is None:
        cache = Path(os.environ.get("SKARSTROM_COOLPROP_CACHE",
                                    Path(tempfile.gettempdir()) / f"skarstrom-coolprop-{COOLPROP_COMMIT[:12]}"))
        cache.mkdir(parents=True, exist_ok=True)
    out: dict[str, Path] = {}
    for n in names:
        file, sha = FLUID_FILES[n]
        path = (fluids if fluids is not None else cache) / file
        if fluids is None and not (path.exists() and hashlib.sha256(path.read_bytes()).hexdigest() == sha):
            url = URL.format(commit=COOLPROP_COMMIT, file=file)
            try:
                with urllib.request.urlopen(url, timeout=60) as r:
                    data = r.read()
            except OSError as e:
                sys.exit(f"cannot download {url}: {e} (or pass --fluids DIR)")
            path.write_bytes(data)
        if not path.exists():
            sys.exit(f"{path} does not exist")
        got = hashlib.sha256(path.read_bytes()).hexdigest()
        if got != sha:
            sys.exit(f"{path}: SHA-256 {got} is not the pinned {sha} (CoolProp {COOLPROP_COMMIT[:12]})")
        out[n] = path
    return out


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


def generate(fluids: Path | None) -> str:
    names = component_names()
    files = fluid_files(fluids, names)
    rows = "\n".join(row(n, files[n]) for n in names)
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
// Source: the CoolProp {COOLPROP_VERSION} fluid files, dev/fluids/*.json at CoolProp
// commit {COOLPROP_COMMIT[:12]} (one file per row, named in the comment above the
// row; pinned by SHA-256 in the generator). CoolProp is a development-time
// data source only: the engine does not build or link it.
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
      "CoolProp {COOLPROP_VERSION} fluid files (generated table)";

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
    ap.add_argument("--fluids", type=Path, default=None,
                    help="local CoolProp dev/fluids directory at the pinned commit (default: download)")
    args = ap.parse_args()
    if args.fluids is not None and not args.fluids.is_dir():
        sys.exit(f"{args.fluids} is not a directory")
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
