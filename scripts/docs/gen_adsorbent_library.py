#!/usr/bin/env python3
"""Validate the curated adsorbent library and generate its docs page.

The library has ONE source of truth:
  data/adsorbents/library.json   catalogue: citation + per-parameter provenance
  data/adsorbents/<key>.json     cereal single-Adsorbent blobs, named after the key

This script VALIDATES the catalogue against the blobs (every entry cited or
marked legacy, metadata values equal to the blob values, no orphan blob files)
and the MATERIALS rule: every entry belongs to one material, every material has
exactly one recommended entry, and a material's second entry is a variant with
a name that sets it apart (the entry's label is the material name, plus the
variant in brackets). Each isotherm records its data quality
and enforces the COMPLETENESS rule: for a sourced entry every required field
(DOI, isotherm fit range, heat of adsorption, particle density, heat capacity,
particle size, example decks) is either sourced/derived or carries a ref
beginning "Not stated in source (checked: ...)"; a field that could not be
checked against the source must be listed in the entry's `unverified_fields`.
Nothing may be silently missing. Tabulated source data (`source_points`, per
isotherm) are validated and summarised.

It then GENERATES three files (never edit them by hand): docs/adsorbent-library.md,
data/adsorbents/library.merged.json (the validated catalogue with every entry's
serialised adsorbent inlined, in the order of the docs page) and
data/adsorbents/library.resolved.json (every entry with the engine's own
resolution applied: kind tags and labels, parameters in engine order with
units, source points in SI, completeness, property units and provenance, mass
transfer). Tools that embed the library read these files, so they need no
validator of their own and no copy of the engine's rules. The resolved file
is checked against data/adsorbents/library.resolved.schema.json.

The resolved file takes its labels, units and parameter order from the engine
(skarstrom.model_kinds(), the parameter registry, sim.describe()), so writing
it needs the skarstrom package (pip install skarstrom, or pip install . from a
checkout). Without the package, --check still validates the committed file
against the schema and the catalogue's entry list; python/tests/
test_library_resolved.py compares a full regeneration.

    python scripts/docs/gen_adsorbent_library.py          # validate + write the page
    python scripts/docs/gen_adsorbent_library.py --check  # validate + fail on drift

Standard library only (the skarstrom package only for the resolved file).
"""

from __future__ import annotations

import argparse
import json
import math
import re
import sys
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parent))
from _fmt import js_exp, js_prec  # noqa: E402

ROOT = Path(__file__).resolve().parents[2]
LIB_DIR = ROOT / "data" / "adsorbents"
CATALOG = LIB_DIR / "library.json"
DECK_DIR = ROOT / "examples" / "decks"
DOCS_OUT = ROOT / "docs" / "adsorbent-library.md"
MERGED_OUT = LIB_DIR / "library.merged.json"
RESOLVED_OUT = LIB_DIR / "library.resolved.json"
RESOLVED_SCHEMA = LIB_DIR / "library.resolved.schema.json"
RESOLVED_VERSION = 1

# adsim::IsothermType enum order (core/include/adsim/equations/isotherm.h) ->
# kind tag (adsim::IsothermKindString in adsim/equations/isotherm_kinds.h).
ISOTHERM_KINDS = [
    "single_site_langmuir_p", "dual_site_langmuir_p", "inert", "henry",
    "single_site_langmuir_c", "dual_site_langmuir_c",
    "extended_single_site_langmuir_p", "extended_dual_site_langmuir_p",
    "extended_single_site_langmuir_c", "extended_dual_site_langmuir_c",
    "anti_langmuir", "freundlich", "sips", "sips_2", "extended_sips_2",
    "dual_site_sips", "redlich_peterson", "toth", "toth_2", "arrhenius_sips",
    "extended_langmuir_freundlich", "gab", "gab_young", "toth_wadst", "toth_humid",
]
# Display labels (mirror isothermModelNames in isotherm.h).
ISOTHERM_LABELS = {
    "single_site_langmuir_p": "Single Site Langmuir",
    "dual_site_langmuir_p": "Dual Site Langmuir",
    "inert": "Inert",
    "henry": "Henry",
    "single_site_langmuir_c": "Single Site Langmuir (Concentration)",
    "dual_site_langmuir_c": "Dual Site Langmuir (Concentration)",
    "extended_single_site_langmuir_p": "Extended Single Site Langmuir",
    "extended_dual_site_langmuir_p": "Extended Dual Site Langmuir",
    "extended_single_site_langmuir_c": "Extended Single Site Langmuir (Concentration)",
    "extended_dual_site_langmuir_c": "Extended Dual Site Langmuir (Concentration)",
    "anti_langmuir": "Anti-Langmuir",
    "freundlich": "Freundlich",
    "sips": "Sips",
    "sips_2": "Sips 2",
    "extended_sips_2": "Extended Sips 2",
    "dual_site_sips": "Dual Site Sips",
    "redlich_peterson": "Redlich-Peterson",
    "toth": "Toth",
    "toth_2": "Toth 2",
    "arrhenius_sips": "Arrhenius Sips (Schell 2012)",
    "extended_langmuir_freundlich": "Extended Langmuir-Freundlich",
    "gab": "GAB",
    "gab_young": "GAB (Young 2021 energies)",
    "toth_wadst": "Toth WADST (CO2 with water, Young 2021)",
    "toth_humid": "Toth humid (CO2 with water, Stampi-Bombelli 2020)",
}

FAMILIES = ["zeolite", "carbon", "MOF", "silica", "alumina", "amine", "other", "support"]
# A "support" material is a non-adsorbing solid (a monolith wall, a laminate
# sheet, ceramic paper) mixed by mass with a sorbent in a structured layer. It
# has no adsorbing component, its density is the apparent density of the wall
# or sheet (particle voidage 0 by convention, so no gas sits in its pores), and
# it has no particle size.
SUPPORT_FAMILY = "support"
STATUSES = ["sourced", "legacy"]
SCHEMA_VERSION = 3
# Per-isotherm data quality, best first. measured_tabulated: measured by the
# cited authors, and the entry holds their tabulated data (source_points);
# measured_fit: measured by the cited authors, fit only; refit: the cited paper
# fitted data measured by others; secondary: the cited paper copied the
# parameters from another paper; not_checked: the source has not yet been read
# for this (needs a data_quality_note); unknown: legacy entries only.
DATA_QUALITY = ["measured_tabulated", "measured_fit", "refit", "secondary", "not_checked", "unknown"]
DATA_QUALITY_LABELS = {
    "measured_tabulated": "measured, data tabulated",
    "measured_fit": "measured, fit only",
    "refit": "fitted to other authors' data",
    "secondary": "copied from an earlier paper",
    "not_checked": "not yet checked",
    "unknown": "unknown (legacy)",
}
MATERIAL_ID_RE = re.compile(r"^[a-z0-9]+(?:-[a-z0-9]+)*$")
ASCII_RE = re.compile(r"^[ -~]+$")
PROVENANCE = ["source", "derived", "repo", "default", "convention", "not_stated", "unknown"]
NOT_STATED_RE = re.compile(r"^Not stated in source(?: as a [^(]+)? \(checked: [^)]+\)")
UNVERIFIABLE_FIELDS = ["T_range_K", "P_range_Pa", "heat_of_adsorption", "particle_density", "heat_capacity", "particle_diameter"]
P_UNITS = {"Pa": 1, "kPa": 1e3, "MPa": 1e6, "bar": 1e5}
Q_UNITS = {"mol/kg": 1, "mmol/g": 1}
T_UNITS = {"K": 0, "C": 273.15}
POINT_KINDS = ["measured", "regenerated_from_fit"]
POINT_QUANTITIES = ["absolute", "as_measured", "excess"]
# Catalogue property key -> (label, unit, blob field it must equal or None).
PROPERTIES = {
    "particle_density": ("Particle density", "kg/m3", "Density"),
    "skeletal_density": ("Skeletal density", "kg/m3", None),
    "heat_capacity": ("Heat capacity", "J/(kg K)", "Specific_Heat_Capacity"),
    "particle_voidage": ("Particle voidage", "-", "Particle_Void_Fraction"),
    "particle_diameter": ("Particle diameter", "m", None),
    "thermal_conductivity": ("Thermal conductivity", "W/(m K)", "Thermal_Conductivity"),
}
REQUIRED_PROPERTIES = ["particle_density", "heat_capacity", "particle_voidage", "particle_diameter", "thermal_conductivity"]
# A support solid has no particle size; its voidage is 0 by convention.
SUPPORT_REQUIRED_PROPERTIES = ["particle_density", "heat_capacity", "particle_voidage", "thermal_conductivity"]
SUPPORT_SOURCED_PROPERTIES = ["particle_density", "heat_capacity", "thermal_conductivity"]
DOI_RE = re.compile(r"^10\.\d{4,9}/\S+$")
KEY_RE = re.compile(r"^[a-z0-9_]+$")


def si(x: float) -> float:
    return float(f"{x:.12g}")


def rel_eq(a: float, b: float) -> bool:
    return abs(a - b) <= 1e-9 * max(1.0, abs(a), abs(b))


def is_num(x) -> bool:
    return isinstance(x, (int, float)) and not isinstance(x, bool) and math.isfinite(x)


def blob_view(blob: dict) -> dict:
    v = blob["value0"]
    src = v["Adsorption_Source"]
    iso_map = {m["key"]: m["value"] for m in src.get("Isotherm_Map", [])}
    hoa = {h["key"]: h["value"]["Value"] for h in src.get("Heat_Of_Adsorption", [])}
    params = {}
    for iso in src.get("Isotherms", []):
        data = ((iso.get("value") or {}).get("ptr_wrapper") or {}).get("data") or {}
        for p in data.get("Isotherm_Params", []):
            params[f"{iso['key']}:{p['key']}"] = [
                {"name": k[3:], "value": x["Value"]}
                for k, x in p["value"].items()
                if k.startswith("IP_") and isinstance(x, dict) and "Value" in x
            ]
    return {
        "name": v["Name"],
        "scalars": {
            "Density": v["Density"]["Value"],
            "Specific_Heat_Capacity": v["Specific_Heat_Capacity"]["Value"],
            "Thermal_Conductivity": v["Thermal_Conductivity"]["Value"],
            "Particle_Void_Fraction": v["Particle_Void_Fraction"]["Value"],
        },
        "iso_map": iso_map,
        "hoa": hoa,
        "params_for": lambda comp: params.get(f"{iso_map.get(comp)}:{comp}", []),
        "mtc_components": [m["key"] for m in src.get("Mass_Transfer_Models", [])],
    }


def check_temperature_table(err, t, value, entry_source_ids):
    """Validate a property's optional temperature_table: {"T_K": [...],
    "values": [...], "source": <entry source id>, "ref": "..."}, values in the
    property's SI unit at increasing temperatures, transcribed from a table
    (never a plot). The property's own value must be one of the tabulated
    values or lie within their range (it is the value at the reference
    temperature the engine uses). Returns the table for the docs, or None."""
    if not isinstance(t, dict):
        err("must be an object")
        return None
    T, v = t.get("T_K"), t.get("values")
    n = len(T) if isinstance(T, list) else -1
    if n < 2 or not isinstance(v, list) or len(v) != n or not all(is_num(x) for x in T + v):
        err("needs equal-length arrays T_K and values (two or more finite numbers)")
        return None
    if any(b <= a for a, b in zip(T, T[1:])) or T[0] <= 0:
        err("T_K must be positive and increasing")
    if any(x <= 0 for x in v):
        err("values must be positive")
    if t.get("source") not in entry_source_ids:
        err(f"source '{t.get('source')}' is not one of this entry's sources")
    if not t.get("ref"):
        err("needs a ref (the table it was transcribed from)")
    if is_num(value) and not (min(v) * (1 - 1e-9) <= value <= max(v) * (1 + 1e-9)):
        err(f"the property's value {value} is outside the tabulated range")
    return {"T_K": T, "values": v, "ref": t.get("ref"), "note": t.get("note") or ""}


def load_library(catalog=None):
    """Validate the catalogue against the blobs. Returns (entries, errors)."""
    catalog = catalog if catalog is not None else json.loads(CATALOG.read_text(encoding="utf-8"))
    errors: list = []

    def err(key, msg):
        errors.append(f"{key}: {msg}")

    if catalog.get("schema_version") != SCHEMA_VERSION:
        errors.append(f"schema_version must be {SCHEMA_VERSION} (got {catalog.get('schema_version')})")
    sources = catalog.get("sources") or {}
    for sid, s in sources.items():
        if not s.get("citation"):
            err(f"source {sid}", "missing citation")
        if not DOI_RE.match(s.get("doi") or ""):
            err(f"source {sid}", f"invalid DOI '{s.get('doi')}'")

    materials = catalog.get("materials")
    if not isinstance(materials, dict) or not materials:
        errors.append("catalog needs a 'materials' object (material id -> name, family, recommended entry)")
        materials = {}
    for mid, mat in materials.items():
        where = f"material {mid}"
        if not MATERIAL_ID_RE.match(mid):
            err(where, "id must be lower case words joined by hyphens ([a-z0-9]+(-[a-z0-9]+)*)")
        if not mat.get("name") or not ASCII_RE.match(mat.get("name") or ""):
            err(where, "needs a name (printable ASCII: it becomes the adsorbent name in a deck)")
        if mat.get("family") not in FAMILIES:
            err(where, f"family '{mat.get('family')}' not in {'/'.join(FAMILIES)}")
        if not mat.get("recommended"):
            err(where, "needs a recommended entry key")
        if not mat.get("recommended_reason"):
            err(where, "needs a recommended_reason (one line: why this entry)")
        if not isinstance(mat.get("aliases"), list):
            err(where, "aliases must be an array (removed entry keys that now resolve to this material)")

    seen_keys, seen_files, entries = set(), set(), []
    for e in catalog.get("entries") or []:
        k = e.get("key") or "<no key>"
        if not KEY_RE.match(k):
            err(k, "key must match [a-z0-9_]+")
        if k in seen_keys:
            err(k, "duplicate key")
        seen_keys.add(k)
        if not e.get("label"):
            err(k, "missing label")
        if not e.get("material"):
            err(k, "missing material")
        if "family" in e:
            err(k, "family belongs to the material (materials.<id>.family), not the entry")
        mat = materials.get(e.get("material_id"))
        if mat is None:
            err(k, f"material_id '{e.get('material_id')}' is not in materials")
            mat = {}
        variant = e.get("variant") or ""
        if variant and not ASCII_RE.match(variant):
            err(k, "variant must be printable ASCII (it becomes part of the adsorbent name)")
        if mat.get("name"):
            want = f"{mat['name']} ({variant})" if variant else mat["name"]
            if e.get("label") != want:
                err(k, f"label must be the material name plus the variant in brackets: '{want}'")
        if e.get("status") not in STATUSES:
            err(k, f"status '{e.get('status')}' not in {'/'.join(STATUSES)}")
        if not isinstance(e.get("applications"), list) or not e.get("applications"):
            err(k, "needs at least one application")
        if e.get("status") == "sourced":
            if not isinstance(e.get("sources"), list) or not e.get("sources"):
                err(k, "sourced entry needs at least one source")
            for r in e.get("sources") or []:
                if r.get("id") not in sources:
                    err(k, f"unknown source id '{r.get('id')}'")
                if not r.get("refs"):
                    err(k, f"source '{r.get('id')}' needs table/figure refs")
        elif not re.search(r"source unknown", e.get("legacy_note") or "", re.I):
            err(k, "legacy entry must carry a legacy_note starting 'Source unknown'")

        # The blob is named after the key: shell-safe by construction (see
        # cmake/check_data_names.cmake).
        if e.get("file") != f"{k}.json":
            err(k, f"blob file '{e.get('file')}' must be named '{k}.json' (the entry key; shell-safe file name)")
        f = LIB_DIR / (e.get("file") or "")
        if not e.get("file") or not f.is_file():
            err(k, f"blob file '{e.get('file')}' not found")
            continue
        if e["file"] in seen_files:
            err(k, f"blob file '{e['file']}' used twice")
        seen_files.add(e["file"])
        blob_text = f.read_text(encoding="utf-8")
        try:
            blob = json.loads(blob_text)
        except json.JSONDecodeError as x:
            err(k, f"blob is not JSON: {x}")
            continue
        if "value0" not in blob or len(blob) != 1:
            err(k, "blob must be a single cereal 'value0' Adsorbent")
            continue
        b = blob_view(blob)

        adsorbing = [c for c, t in b["iso_map"].items() if ISOTHERM_KINDS[t] != "inert"]
        comps = e.get("components") or []
        if sorted(comps) != sorted(adsorbing):
            err(k, f"components {comps} != blob isotherm components {adsorbing}")
        legacy = e.get("status") == "legacy"
        entry_source_ids = {r.get("id") for r in e.get("sources") or []}

        unverified = e.get("unverified_fields") or []
        if not isinstance(unverified, list):
            err(k, "unverified_fields must be an array")
            unverified = []
        for fld in unverified:
            if fld not in UNVERIFIABLE_FIELDS:
                err(k, f"unverified_fields: unknown field '{fld}' (allowed: {', '.join(UNVERIFIABLE_FIELDS)})")
        if unverified and legacy:
            err(k, "unverified_fields is for sourced entries (legacy entries are incomplete by definition)")
        if unverified and not e.get("verification_note"):
            err(k, "unverified_fields needs a verification_note saying why the fields could not be checked")
        missing: list = []
        satisfied_unverified: set = set()

        def need(field, ok, what=None):
            label = f"{field} ({what})" if what else field
            if ok:
                if field in unverified:
                    satisfied_unverified.add(field)
                return
            missing.append(label)
            if not legacy and field not in unverified:
                err(k, f"{label} is neither sourced nor marked 'Not stated in source (checked: ...)' "
                       "(list it in unverified_fields only if the source could not be checked)")

        def not_stated(text):
            return isinstance(text, str) and bool(NOT_STATED_RE.match(text))

        need("doi", not legacy and bool(e.get("sources"))
             and all(DOI_RE.match((sources.get(r.get("id")) or {}).get("doi") or "") for r in e.get("sources")))

        iso_meta = {i["component"]: i for i in e.get("isotherms") or []}
        points_summary = {}
        for comp in comps:
            im = iso_meta.get(comp)
            if not im:
                err(k, f"no isotherms[] entry for {comp}")
                continue
            if not im.get("ref"):
                err(k, f"isotherm {comp} needs a ref")
            dq = im.get("data_quality")
            if dq not in DATA_QUALITY:
                err(k, f"isotherm {comp}: data_quality '{dq}' not one of {'/'.join(DATA_QUALITY)}")
            if (dq == "unknown") != legacy:
                err(k, f"isotherm {comp}: data_quality 'unknown' is for legacy entries, and only it")
            if dq == "not_checked" and not im.get("data_quality_note"):
                err(k, f"isotherm {comp}: data_quality 'not_checked' needs a data_quality_note saying what is left to check")
            sp_kind = (im.get("source_points") or {}).get("kind")
            if (dq == "measured_tabulated") != (sp_kind == "measured"):
                err(k, f"isotherm {comp}: data_quality 'measured_tabulated' goes with, and only with, measured source_points")
            for rk, nk in (("T_range_K", "T_range_not_stated"), ("P_range_Pa", "P_range_not_stated")):
                r = im.get(rk)
                if r is None:
                    if nk in im and not not_stated(im[nk]):
                        err(k, f"isotherm {comp}: {nk} must begin 'Not stated in source (checked: ...)'")
                    need(rk, not_stated(im.get(nk)), comp)
                    continue
                if nk in im:
                    err(k, f"isotherm {comp}: {nk} given although {rk} is stated")
                if not (isinstance(r, list) and len(r) == 2 and all(is_num(x) for x in r) and r[0] <= r[1]):
                    err(k, f"isotherm {comp}: {rk} must be [min, max] or null")
                need(rk, True, comp)
            sp = im.get("source_points")
            if sp is not None:
                where = f"isotherm {comp}: source_points"
                if sp.get("kind") not in POINT_KINDS:
                    err(k, f"{where}.kind must be one of {'/'.join(POINT_KINDS)}")
                if sp.get("quantity") not in POINT_QUANTITIES:
                    err(k, f"{where}.quantity must be one of {'/'.join(POINT_QUANTITIES)}")
                u = sp.get("units") or {}
                if u.get("p") not in P_UNITS:
                    err(k, f"{where}.units.p '{u.get('p')}' not one of {'/'.join(P_UNITS)}")
                if u.get("q") not in Q_UNITS:
                    err(k, f"{where}.units.q '{u.get('q')}' not one of {'/'.join(Q_UNITS)}")
                if u.get("T") not in T_UNITS:
                    err(k, f"{where}.units.T '{u.get('T')}' not one of {'/'.join(T_UNITS)}")
                if sp.get("source") not in entry_source_ids:
                    err(k, f"{where}.source '{sp.get('source')}' is not one of this entry's sources")
                if not sp.get("ref"):
                    err(k, f"{where} needs a ref (table / SI file it was transcribed from)")
                if not sp.get("note"):
                    err(k, f"{where} needs a note (measured vs regenerated, what the data are)")
                if sp.get("kind") == "regenerated_from_fit" and not re.search(r"regenerat", sp.get("note") or "", re.I):
                    err(k, f"{where}: regenerated points must say so in the note")
                n = len(sp["p"]) if isinstance(sp.get("p"), list) else -1
                arrays_ok = n > 0 and all(isinstance(sp.get(a), list) and len(sp[a]) == n and all(is_num(x) for x in sp[a])
                                          for a in ("T", "p", "q"))
                if not arrays_ok:
                    err(k, f"{where}: T, p, q must be equal-length non-empty arrays of finite numbers")
                if arrays_ok and u.get("p") in P_UNITS and u.get("q") in Q_UNITS and u.get("T") in T_UNITS:
                    T_K = [si(t + T_UNITS[u["T"]]) for t in sp["T"]]
                    p_Pa = [si(x * P_UNITS[u["p"]]) for x in sp["p"]]
                    if any(t <= 0 for t in T_K) or any(x < 0 for x in p_Pa):
                        err(k, f"{where}: T must be > 0 K and p >= 0")
                    points_summary[comp] = {
                        "kind": sp["kind"], "quantity": sp["quantity"], "n": n,
                        "temperatures_K": sorted(set(T_K)), "ref": sp["ref"], "note": sp["note"],
                    }
        for comp in iso_meta:
            if comp not in comps:
                err(k, f"isotherms[] lists {comp} which is not a component")

        hoa_out = []
        for comp in comps:
            h = (e.get("heat_of_adsorption") or {}).get(comp)
            if not h:
                err(k, f"heat_of_adsorption missing for {comp}")
                continue
            prov = h.get("provenance")
            if prov not in PROVENANCE:
                err(k, f"heat_of_adsorption {comp}: bad provenance '{prov}'")
            if not h.get("ref"):
                err(k, f"heat_of_adsorption {comp}: needs ref")
            if prov == "unknown" and not legacy:
                err(k, f"heat_of_adsorption {comp}: provenance 'unknown' is only for legacy entries")
            if prov == "not_stated" and not not_stated(h.get("ref")):
                err(k, f"heat_of_adsorption {comp}: not_stated needs a ref beginning 'Not stated in source (checked: ...)'")
            if prov == "not_stated" and h.get("value") is not None:
                err(k, f"heat_of_adsorption {comp}: not_stated needs value null")
            if (h.get("value") is None and prov != "not_stated" and not (legacy and prov == "unknown")
                    and "heat_of_adsorption" not in unverified):
                err(k, f"heat_of_adsorption {comp}: a null value needs provenance 'not_stated'")
            bv = b["hoa"].get(comp)
            if h.get("value") is None:
                if bv is not None and bv != 0:
                    err(k, f"heat_of_adsorption {comp}: catalog null but blob has {bv}")
            elif bv is None or not rel_eq(h["value"], bv):
                err(k, f"heat_of_adsorption {comp}: catalog {h['value']} != blob {bv}")
            need("heat_of_adsorption", (prov in ("source", "derived") and h.get("value") is not None) or not_stated(h.get("ref")), comp)
            hoa_out.append({"component": comp, "value": h.get("value"), "provenance": prov, "ref": h.get("ref")})

        support = mat.get("family") == SUPPORT_FAMILY
        if support and comps:
            err(k, "a support solid adsorbs nothing: its blob must give no component a non-inert isotherm")
        props = e.get("properties") or {}
        for p in (SUPPORT_REQUIRED_PROPERTIES if support else REQUIRED_PROPERTIES):
            if p not in props:
                err(k, f"property {p} missing")
        if support and "particle_diameter" in props:
            err(k, "a support solid has no particle_diameter")
        if support and (props.get("particle_voidage") or {}).get("value") != 0:
            err(k, "a support solid has particle_voidage 0 (its density is the apparent density of the wall or sheet)")
        props_out = []
        for p, rec in props.items():
            d = PROPERTIES.get(p)
            if not d:
                err(k, f"unknown property {p}")
                continue
            prov = rec.get("provenance")
            if prov not in PROVENANCE:
                err(k, f"property {p}: bad provenance '{prov}'")
            if not rec.get("ref"):
                err(k, f"property {p}: needs ref")
            if rec.get("value") is not None and not is_num(rec.get("value")):
                err(k, f"property {p}: value must be a number or null")
            if prov == "unknown" and not legacy:
                err(k, f"property {p}: provenance 'unknown' is only for legacy entries")
            if prov == "convention" and not (support and p == "particle_voidage"):
                err(k, f"property {p}: provenance 'convention' is only for the particle_voidage (0) of a support solid")
            table_out = None
            if "temperature_table" in rec:
                table_out = check_temperature_table(lambda m: err(k, f"property {p}: temperature_table {m}"),
                                                    rec["temperature_table"], rec.get("value"), entry_source_ids)
            if prov == "not_stated" and (rec.get("value") is not None or not not_stated(rec.get("ref"))):
                err(k, f"property {p}: not_stated needs value null and a ref beginning 'Not stated in source (checked: ...)'")
            if rec.get("value") is None and prov != "not_stated" and not (legacy and prov == "unknown") and p not in unverified:
                err(k, f"property {p}: a null value needs provenance 'not_stated'")
            if d[2]:
                bv = b["scalars"][d[2]]
                if rec.get("value") is None or not rel_eq(rec["value"], bv):
                    err(k, f"property {p}: catalog {rec.get('value')} != blob {d[2]} {bv}")
            props_out.append({"name": p, "label": d[0], "unit": d[1], "value": rec.get("value"),
                              "provenance": prov, "ref": rec.get("ref"), "temperature_table": table_out})
        for p in (SUPPORT_SOURCED_PROPERTIES if support else ("particle_density", "heat_capacity", "particle_diameter")):
            rec = props.get(p)
            need(p, bool(rec) and ((rec.get("provenance") in ("source", "derived") and rec.get("value") is not None)
                                   or not_stated(rec.get("ref"))))

        decks = e.get("example_decks")
        if not isinstance(decks, list):
            err(k, "example_decks must be an array (empty when no deck uses the entry)")
            decks = []
        for dk in decks:
            if not (DECK_DIR / dk).is_file():
                err(k, f"example deck '{dk}' not found in examples/decks/")
        if decks and e.get("example_deck_note"):
            err(k, "example_deck_note is only for entries no deck uses")
        need("example_decks", bool(decks) or bool(e.get("example_deck_note")))
        for dc in e.get("deck_checks") or []:
            if dc.get("deck") not in decks:
                err(k, f"deck_check deck '{dc.get('deck')}' not listed in example_decks")
            if not dc.get("adsorbent"):
                err(k, "deck_check needs an adsorbent name")
        for fld in unverified:
            if fld in satisfied_unverified and not any(m == fld or m.startswith(f"{fld} (") for m in missing):
                err(k, f"unverified_fields lists '{fld}' but it is satisfied - remove it")

        meta = {
            "key": k,
            "label": e.get("label"),
            "material": e.get("material"),
            "material_id": e.get("material_id"),
            "material_name": mat.get("name", ""),
            "variant": variant,
            "recommended": mat.get("recommended") == k,
            "family": mat.get("family"),
            "support": support,
            "status": e.get("status"),
            "legacy_note": e.get("legacy_note") or "",
            "applications": e.get("applications") or [],
            "components": comps,
            "sources": [{
                "id": r.get("id"),
                "citation": (sources.get(r.get("id")) or {}).get("citation", ""),
                "doi": (sources.get(r.get("id")) or {}).get("doi", ""),
                "url": f"https://doi.org/{sources[r['id']]['doi']}" if (sources.get(r.get("id")) or {}).get("doi") else "",
                "refs": r.get("refs"),
            } for r in e.get("sources") or []],
            "isotherms": [{
                "component": comp,
                "kind_label": ISOTHERM_LABELS.get(ISOTHERM_KINDS[b["iso_map"][comp]], "unknown") if comp in b["iso_map"] else "unknown",
                "params": b["params_for"](comp),
                "ref": (iso_meta.get(comp) or {}).get("ref", ""),
                "data_quality": (iso_meta.get(comp) or {}).get("data_quality", ""),
                "data_quality_note": (iso_meta.get(comp) or {}).get("data_quality_note", ""),
                "T_range_K": (iso_meta.get(comp) or {}).get("T_range_K"),
                "P_range_Pa": (iso_meta.get(comp) or {}).get("P_range_Pa"),
                "T_range_not_stated": (iso_meta.get(comp) or {}).get("T_range_not_stated", ""),
                "P_range_not_stated": (iso_meta.get(comp) or {}).get("P_range_not_stated", ""),
                "range_note": (iso_meta.get(comp) or {}).get("range_note", ""),
                "source_points": points_summary.get(comp),
            } for comp in comps],
            "heat_of_adsorption": hoa_out,
            "properties": props_out,
            "example_decks": decks,
            "example_deck_note": e.get("example_deck_note") or "",
            "completeness": {"complete": not missing, "missing": missing,
                             "verification_note": e.get("verification_note") or ""},
            "notes": e.get("notes") or [],
        }
        entries.append({"meta": meta, "file": e["file"], "points": points_summary})

    for f in sorted(LIB_DIR.iterdir()):
        if f.suffix == ".json" and f.name not in ("library.json", MERGED_OUT.name, RESOLVED_OUT.name, RESOLVED_SCHEMA.name) and f.name not in seen_files:
            errors.append(f"orphan blob '{f.name}' is not listed in library.json")

    # Materials: one recommended entry each; variants only where a material
    # has more than one entry, each with its own name; ids, keys and aliases
    # never collide (a material id or an alias resolves to the recommended entry).
    by_key = {en["meta"]["key"]: en["meta"] for en in entries}
    seen_aliases: set = set()
    for mid, mat in materials.items():
        where = f"material {mid}"
        members = [m for m in by_key.values() if m["material_id"] == mid]
        if not members:
            err(where, "has no entries")
            continue
        rec = by_key.get(mat.get("recommended"))
        if rec is None or rec["material_id"] != mid:
            err(where, f"recommended '{mat.get('recommended')}' is not one of its entries")
        elif rec["status"] == "legacy" and any(m["status"] == "sourced" for m in members):
            err(where, "a legacy entry cannot be recommended while a sourced entry exists")
        if len(members) == 1 and members[0]["variant"]:
            err(where, f"its only entry {members[0]['key']} must not have a variant")
        if len(members) > 1:
            names = [m["variant"] for m in members]
            for m in members:
                if not m["variant"]:
                    err(m["key"], f"material {mid} has {len(members)} entries, so this one needs a variant naming what sets it apart")
            if len(set(names)) != len(names):
                err(where, "variants must be unique within the material")
        if mid in seen_keys:
            err(where, "id collides with an entry key")
        for a in mat.get("aliases") or []:
            if a in seen_keys or a in materials or a in seen_aliases:
                err(where, f"alias '{a}' collides with an entry key, a material id or another alias")
            seen_aliases.add(a)
    return entries, errors


# ---- rendering -------------------------------------------------------------------

def md(s) -> str:
    s = "" if s is None else str(s)
    return (s.replace("&", "&amp;").replace("<", "&lt;").replace(">", "&gt;")
             .replace("{", "&#123;").replace("}", "&#125;").replace("$", "&#36;")
             .replace("|", "\\|").replace("*", "\\*").replace("_", "\\_"))


def fmt_num(x) -> str:
    if x is None:
        return "—"
    if x == 0:
        return "0"
    a = abs(x)
    if 1e-3 <= a < 1e6:
        return js_prec(x, 6)
    return re.sub(r"\.?0+e", "e", js_exp(x, 4))


def fmt_range(r, scale=1.0, unit="") -> str:
    return f"{fmt_num(r[0] * scale)}–{fmt_num(r[1] * scale)}{unit}" if r else "not stated"


PROV = {"source": "source", "derived": "derived", "repo": "repo deck/test", "default": "engine default",
        "convention": "library convention", "not_stated": "not stated", "unknown": "unknown"}


def anchor(key: str) -> str:
    return key.replace("_", "-")


def material_groups(entries):
    """[(material_id, name, family, [entries, recommended first])] in catalogue order of materials."""
    catalog = json.loads(CATALOG.read_text(encoding="utf-8"))
    groups = []
    for mid, mat in (catalog.get("materials") or {}).items():
        members = [e for e in entries if e["meta"]["material_id"] == mid]
        members.sort(key=lambda e: not e["meta"]["recommended"])
        groups.append((mid, mat, members))
    return groups


def render_docs(entries) -> str:
    L: list = []
    out = L.append
    n_sourced = sum(1 for e in entries if e["meta"]["status"] == "sourced")
    groups = material_groups(entries)
    out("<!-- GENERATED by scripts/docs/gen_adsorbent_library.py from data/adsorbents/library.json — do not edit. -->")
    out("")
    out("# Adsorbent Library")
    out("")
    out(f"Skarstrom ships a curated library of **{len(groups)}** materials ({len(entries)} entries: {n_sourced} sourced, {len(entries) - n_sourced} legacy) in `data/adsorbents/`. Every material has one **recommended** entry, chosen for the quality of its experimental data; that is the entry you get when you add the material by its id. A material has a second entry only as a *variant*: a different product, form or data range, named in the entry's label. Every sourced entry cites the publication (DOI and table/figure) its parameters were transcribed from, and every property records its provenance:")
    out("")
    out("- **source** — transcribed from the cited table or figure;")
    out("- **derived** — computed from cited values (the formula is given);")
    out("- **repo deck/test** — the value used by the repository's validation deck or test for that paper, but not stated in the paper itself;")
    out("- **engine default** — not given by the source; the engine default is kept;")
    out("- **not stated** — the source does not state the quantity (no value).")
    out("")
    out("A quantity the source does not state is never filled in silently: its reference says *Not stated in source (checked: …)* and names what was checked. Fields that could not be checked against the source at all are listed per entry as *unverified*.")
    out("")
    out("Legacy entries predate the curated library; their origin has not been verified and they are labelled *source unknown*.")
    out("")
    out("Each isotherm also records the **quality of its data**, best first:")
    out("")
    out("- **measured, data tabulated** — measured by the cited authors, and their tabulated data are stored with the entry;")
    out("- **measured, fit only** — measured by the cited authors, who publish the fit but no table of the data;")
    out("- **fitted to other authors' data** — the cited paper fitted data measured by another group (cited there);")
    out("- **copied from an earlier paper** — the cited paper took the parameters from earlier work;")
    out("- **not yet checked** — the source has not yet been read for this; the note says what is left.")
    out("")
    out("The recommended entry is the one with the best data quality, then the best coverage of the material's usual components, then the widest range. A legacy entry is never recommended while a sourced one exists.")
    out("")
    complete = sum(1 for e in entries if e["meta"]["completeness"]["complete"])
    with_points = sum(1 for e in entries if any(i["source_points"] for i in e["meta"]["isotherms"]))
    out(f"**Completeness:** {complete} of {len(entries)} entries carry every required field (DOI, isotherm fit range, heat of adsorption, particle density, heat capacity, particle size, example decks; for a support solid DOI, density, heat capacity and thermal conductivity) as a sourced value or an explicit *not stated in source* note; the others are listed as incomplete below. **Tabulated source data:** {with_points} entries include equilibrium points tabulated in their source (tables, SI files or datasets — never digitised from plots), stored with the entry in `library.json`.")
    out("")
    out("**Support solids** (family *support*; the library ships none at present) are the non-adsorbing walls and sheets of structured contactors. A support has no adsorbing component (every component uses the inert isotherm); mix it by mass with a sorbent in a structured layer's composition. Its density is the apparent density of the wall or sheet, with particle voidage 0 by convention, and it has no particle size; its required fields are density, heat capacity and thermal conductivity. Where the source tabulates a property against temperature, the table is listed under the entry; the engine uses the constant value at the reference temperature (300 K unless stated).")
    out("")
    out("## Using the library")
    out("")
    out("Each entry is a catalogue record in `library.json` plus one serialised adsorbent, `<key>.json`. The Python package ships the library:")
    out("")
    out("```python")
    out("import skarstrom")
    out("lib = skarstrom.adsorbent_library()               # the catalogue as a dict")
    out('list(lib["materials"])                             # material ids')
    out('path = skarstrom.adsorbent_path("zeolite-13x")      # the recommended 13X entry')
    out('path = skarstrom.adsorbent_path("zeolite_13x_streb2021")  # one variant, by its key')
    out('sim = skarstrom.Simulation.new(title="Library demo")')
    out('sim.add_library_adsorbent("zeolite-13x")           # add the recommended entry to a deck')
    out("```")
    out("")
    out("In C++, `adsim/objects/adsorbent_library.h` reads the catalogue (`AdsorbentLibrary::Resolve` turns a material id, an entry key or the alias of a removed entry into an entry), `LoadAdsorbentFromJSON` loads an entry's file, and `AddLibraryAdsorbent` (`adsim/objects/deck_edit.h`) adds an entry, or a material's recommended entry, to a deck with default physics for any component the entry does not cover. The library is curated: see `CONTRIBUTING.md` for the data-integrity rules a new entry must meet (every parameter cited).")
    out("")
    out("### The resolved library")
    out("")
    out('!!! note "New in 0.9"')
    out("    `library.resolved.json` and its schema are new in 0.9.")
    out("")
    out("`data/adsorbents/library.resolved.json` is the library with the engine's own resolution applied, so a front end can embed it as it is and keep no copy of the engine's rules. It is **generated** by `scripts/docs/gen_adsorbent_library.py` (never edit it) and its shape is `data/adsorbents/library.resolved.schema.json`, which the generator checks in `--check`. `library.merged.json` (the catalogue with each serialised adsorbent) is unchanged. Writing the resolved file needs the `skarstrom` package, because the labels, units and parameter order come from the engine (`skarstrom.model_kinds()`, the parameter registry), not from a table in the script.")
    out("")
    out("The top level has `schema_version`, `catalogue_schema_version`, `materials` (id to name, family, recommended entry, reason, aliases) and `entries`, in the order of this page. Each entry holds:")
    out("")
    out("- identity: `key`, `label`, `material_id`, `variant`, `family`, `support`, `recommended`, `status`, `components`, `aliases`, and `sources` with a DOI URL;")
    out("- `isotherms`, one per component: the kind tag and label, `parameters` in the engine's order with name, label, unit and value, the data quality, the fit ranges, and `source_points` in SI (Pa, mol/kg, K) with a summary (count, temperatures, T and P ranges);")
    out("- `heat_of_adsorption` and `properties` with unit, value and provenance (`source`, `derived`, `repo`, `default`, `convention`, `not_stated`), and the registry path of each property;")
    out("- `mass_transfer`, one per component: the rate model, the coefficient kind and its parameters with units;")
    out("- `completeness`: `complete`, the `missing` fields and the `unverified_fields` with the verification note. A gap is reported here, never filled in.")
    out("")
    out("```python")
    out("import json, os")
    out("import skarstrom")
    out('path = os.path.join(skarstrom.adsorbent_library_dir(), "library.resolved.json")')
    out("with open(path, encoding=\"utf-8\") as f:")
    out("    resolved = json.load(f)")
    out('entry = next(e for e in resolved["entries"] if e["recommended"] and e["family"] == "zeolite")')
    out('iso = entry["isotherms"][0]')
    out('print(entry["label"], iso["kind_label"], [(p["name"], p["unit"]) for p in iso["parameters"]])')
    out("```")
    out("")
    out("!!! warning \"Before you run\"")
    out("    Check the heat of adsorption (some sources report only site energies, so the value is left at 0) and the heat-of-adsorption convention of your deck (isosteric vs internal energy) against the source. Fit ranges are listed per component; outside them the isotherm is an extrapolation.")
    out("")
    out("| Material | Id | Family | Recommended entry | Components | Source | Other variants |")
    out("|---|---|---|---|---|---|---|")
    for mid, mat, members in groups:
        if not members:
            continue
        r = members[0]["meta"]
        src = "source unknown (legacy)" if r["status"] == "legacy" else ", ".join(f"[{md(s['doi'])}]({s['url']})" for s in r["sources"])
        others = ", ".join(f"[{md(e['meta']['variant'])}](#{anchor(e['meta']['key'])})" for e in members[1:]) or "—"
        out(f"| [{md(mat['name'])}](#material-{mid}) | `{mid}` | {md(mat['family'])} | [{md(r['variant'] or r['label'])}](#{anchor(r['key'])}) | {', '.join(md(c) for c in r['components']) or 'none (inert)'} | {src} | {others} |")
    out("")
    for mid, mat, members in groups:
        if not members:
            continue
        out("---")
        out("")
        out(f"## {md(mat['name'])} {{#material-{mid}}}")
        out("")
        rec_line = f"**Material id:** `{mid}` · **Family:** {md(mat['family'])} · **Recommended:** [{md(members[0]['meta']['label'])}](#{anchor(members[0]['meta']['key'])}). {md(mat['recommended_reason'])}"
        if mat.get("aliases"):
            rec_line += f" Removed keys that now resolve here: {', '.join(f'`{a}`' for a in mat['aliases'])}."
        out(rec_line)
        out("")
        for e in members:
            render_entry(out, e)
    return "\n".join(L)


def render_entry(out, e) -> None:
    m, f = e["meta"], e["file"]
    tag = "Recommended" if m["recommended"] else "Variant"
    out(f"### {md(m['label'])} {{#{anchor(m['key'])}}}")
    out("")
    out(f"**{tag}** · **Material:** {md(m['material'])} · **Applications:** {', '.join(md(a) for a in m['applications'])} · **Library key:** `{m['key']}` · **File:** `{f}`")
    out("")
    if m["status"] == "legacy":
        out(f"> **{md(m['legacy_note'])}**")
        out("")
    else:
        out("**Sources**")
        out("")
        for s in m["sources"]:
            out(f"- {md(s['citation'])} [doi:{md(s['doi'])}]({s['url']}) — {md(s['refs'])}")
        out("")
    if m["support"]:
        out("*Non-adsorbing support solid:* every component uses the inert isotherm (a deck gives each component the inert default when the entry is added). Mix it by mass with a sorbent in a layer's composition.")
        out("")
    else:
        out("| Component | Isotherm | Parameters | T range of fit | P range of fit | Data quality | Ref. |")
        out("|---|---|---|---|---|---|---|")
    for i in m["isotherms"]:  # (none for a support solid)
        params = "; ".join(f"{md(p['name'])} = {fmt_num(p['value'])}" for p in i["params"])
        out(f"| {md(i['component'])} | {md(i['kind_label'])} | {params} | {fmt_range(i['T_range_K'], 1, ' K')} | {fmt_range(i['P_range_Pa'], 1e-5, ' bar')} | {DATA_QUALITY_LABELS.get(i['data_quality'], md(i['data_quality']))} | {md(i['ref'])} |")
    if not m["support"]:
        out("")
    dq_notes = [f"{md(i['component'])}: {md(i['data_quality_note'])}" for i in m["isotherms"] if i["data_quality_note"]]
    if dq_notes:
        out(f"*Data quality.* {' '.join(dq_notes)}")
        out("")
    notes = []
    for i in m["isotherms"]:
        parts = [x for x in (i["T_range_not_stated"] and f"T range: {i['T_range_not_stated']}",
                             i["P_range_not_stated"] and f"P range: {i['P_range_not_stated']}",
                             i["range_note"]) if x]
        if parts:
            notes.append(" ".join(f"{md(i['component'])}: {md(t)}" for t in parts))
    if notes:
        out(f"*Fit ranges.* {' '.join(notes)}")
        out("")
    out("| Quantity | Value | Provenance | Reference |")
    out("|---|---|---|---|")
    for h in m["heat_of_adsorption"]:
        val = "not set (0)" if h["value"] is None else f"{fmt_num(h['value'] / 1000)} kJ/mol"
        out(f"| Heat of adsorption, {md(h['component'])} | {val} | {PROV.get(h['provenance'], h['provenance'])} | {md(h['ref'])} |")
    for p in m["properties"]:
        val = "—" if p["value"] is None else f"{fmt_num(p['value'])}{'' if p['unit'] == '-' else ' ' + md(p['unit'])}"
        out(f"| {md(p['label'])} | {val} | {PROV.get(p['provenance'], p['provenance'])} | {md(p['ref'])} |")
    out("")
    tables = [p for p in m["properties"] if p["temperature_table"]]
    for p in tables:
        t = p["temperature_table"]
        out(f"*{md(p['label'])} against temperature* ({md(p['unit'])}; {md(t['ref'])}; the engine uses the constant value above)")
        out("")
        out("| T [K] | " + " | ".join(fmt_num(x) for x in t["T_K"]) + " |")
        out("|---|" + "---|" * len(t["T_K"]))
        out(f"| {md(p['label'])} | " + " | ".join(fmt_num(x) for x in t["values"]) + " |")
        out("")
        if t["note"]:
            out(f"- {md(t['note'])}")
            out("")
    sp = [i for i in m["isotherms"] if i["source_points"]]
    if m["support"]:
        pass
    elif sp:
        out("**Tabulated source data** (stored in `library.json`)")
        out("")
        out("| Component | Points | Temperatures [K] | Quantity | Data | Transcribed from |")
        out("|---|---|---|---|---|---|")
        for i in sp:
            s = i["source_points"]
            temps = ", ".join(fmt_num(t) for t in s["temperatures_K"])
            out(f"| {md(i['component'])} | {s['n']} | {temps} | {md(s['quantity'].replace('_', ' ', 1))} | {md(s['kind'].replace('_', ' '))} | {md(s['ref'])} |")
        out("")
        for i in sp:
            out(f"- {md(i['component'])}: {md(i['source_points']['note'])}")
        out("")
    else:
        out("*No tabulated source data* (the source gives its equilibrium data only as plots, which are not digitised).")
        out("")
    if m["example_decks"]:
        out(f"**Used by example decks:** {', '.join(f'`{d}`' for d in m['example_decks'])}")
    else:
        out(f"**Used by example decks:** none. {md(m['example_deck_note'])}")
    out("")
    if not m["completeness"]["complete"]:
        vn = m["completeness"]["verification_note"]
        out(f"> **Incomplete entry.** Missing or unverified: {', '.join(md(x) for x in m['completeness']['missing'])}.{' ' + md(vn) if vn else ''}")
        out("")
    if m["notes"]:
        for n in m["notes"]:
            out(f"- {md(n)}")
        out("")


def render_merged(entries) -> str:
    """The validated catalogue as one JSON document: materials, sources and every
    entry (catalogue record, provenance and the adsorbent blob) in docs-page order."""
    catalog = json.loads(CATALOG.read_text(encoding="utf-8"))
    raw = {e["key"]: e for e in catalog["entries"]}
    merged = []
    for _mid, _mat, members in material_groups(entries):
        for m in members:
            rec = dict(raw[m["meta"]["key"]])
            rec["recommended"] = m["meta"]["recommended"]
            rec["adsorbent"] = json.loads((LIB_DIR / m["file"]).read_text(encoding="utf-8"))
            merged.append(rec)
    doc = {
        "_about": "GENERATED by scripts/docs/gen_adsorbent_library.py from library.json and the <key>.json blobs; do not edit. "
                  "Validated, merged catalogue: each entry is its library.json record plus `recommended` and its serialised adsorbent (`adsorbent`), in the order of docs/adsorbent-library.md.",
        "schema_version": catalog["schema_version"],
        "materials": catalog["materials"],
        "sources": catalog["sources"],
        "entries": merged,
    }
    return json.dumps(doc, indent=2, ensure_ascii=False) + "\n"


# ---- resolved metadata -----------------------------------------------------------

# Catalogue property key -> leaf of the engine's parameter registry that holds it
# (adsorbent[ID].<leaf>, or layer[ID].<leaf> for the particle diameter). Properties
# with no engine field (the skeletal density is information only) keep the
# catalogue's label and unit.
PROPERTY_REGISTRY = {
    "particle_density": ("adsorbent", "solid_density"),
    "heat_capacity": ("adsorbent", "heat_capacity"),
    "thermal_conductivity": ("adsorbent", "thermal_conductivity"),
    "particle_voidage": ("adsorbent", "particle_void_fraction"),
    "particle_diameter": ("layer", "particle_diameter"),
}


class EngineTables:
    """The engine's own labels, units and parameter order, read through the
    skarstrom package (model_kinds(), the parameter registry, describe())."""

    def __init__(self, sk):
        self.sk = sk
        kinds = sk.model_kinds()
        self.isotherm = {k["tag"]: k for k in kinds["isotherm"]}
        self.isotherm_order = [k["tag"] for k in kinds["isotherm"]]  # enum order (the blob's stored index)
        self.rate_model = {k["tag"]: k["label"] for k in kinds["rate_model"]}
        self.rate_coefficient = {k["tag"]: k["label"] for k in kinds["rate_coefficient"]}

    def resolve_entry(self, blob_names, comps, support, library_key):
        """Add the adsorbent to a scratch deck and read what the registry resolves:
        {"properties": {catalogue key: {name, path, unit}}, "heat_unit", "mass_transfer": [...]}."""
        sk = self.sk
        sim = sk.Simulation.new(title="resolve")
        for c in blob_names:
            sim.add_component(c)
        aid = sim.add_library_adsorbent(library_key)
        name = sim.adsorbents[0].name
        layer_id = None
        if not support:
            sim.add_layer(name, composition={name: 1.0}, length=0.1, cells=5, particle_diameter=1e-3)
            layer_id = 1
        reg = {p.path: p for p in sim.parameters(tunable_only=False)}
        props = {}
        for key, (family, leaf) in PROPERTY_REGISTRY.items():
            ident = aid if family == "adsorbent" else layer_id
            p = reg.get(f"{family}[{ident}].{leaf}") if ident is not None else None
            if p is not None:
                props[key] = {"engine_name": p.name, "registry_path": f"{family}[ID].{leaf}", "unit": p.unit}
        hp = [p for path, p in reg.items() if path.startswith(f"adsorbent[{aid}].heat_of_adsorption[")]
        rate_models = sim.describe()["adsorbents"][0].get("rate_models", {})
        mass_transfer = []
        for comp in comps:
            rm = rate_models.get(comp) or {}
            prefix = f"adsorbent[{aid}].kinetics[{comp}]."
            params, resistances = [], []
            for path, p in reg.items():
                if not path.startswith(prefix):
                    continue
                leaf = path[len(prefix):]
                if leaf.startswith("resistances."):
                    if p.value:
                        resistances.append(leaf.split(".", 1)[1])
                else:
                    params.append({"name": leaf, "label": p.name, "unit": p.unit, "value": p.value})
            mass_transfer.append({
                "component": comp,
                "model": rm.get("model", ""),
                "model_label": self.rate_model.get(rm.get("model", ""), ""),
                "coefficient": rm.get("coefficient", ""),
                "coefficient_label": self.rate_coefficient.get(rm.get("coefficient", ""), ""),
                "resistances": resistances,
                "parameters": params,
            })
        return {"properties": props, "heat_unit": hp[0].unit if hp else "J/mol", "mass_transfer": mass_transfer}


def load_engine():
    """The engine tables, or None when the skarstrom package (with model_kinds) is missing."""
    try:
        import skarstrom  # noqa: PLC0415
        if not hasattr(skarstrom, "model_kinds"):
            return None
        return EngineTables(skarstrom)
    except ImportError:
        return None


def resolved_source_points(sp, source_ids):
    """A catalogue source_points record converted to SI, with its summary."""
    u = sp["units"]
    T = [si(t + T_UNITS[u["T"]]) for t in sp["T"]]
    p = [si(x * P_UNITS[u["p"]]) for x in sp["p"]]
    q = [si(x * Q_UNITS[u["q"]]) for x in sp["q"]]
    return {
        "kind": sp["kind"], "quantity": sp["quantity"], "source": sp["source"], "ref": sp["ref"], "note": sp["note"],
        "summary": {
            "n": len(p), "temperatures_K": sorted(set(T)), "T_range_K": [min(T), max(T)],
            "P_range_Pa": [min(p), max(p)], "q_range_mol_per_kg": [min(q), max(q)],
        },
        "T_K": T, "p_Pa": p, "q_mol_per_kg": q,
    }


def render_resolved(entries, engine) -> str:
    """Every entry with the engine's resolution applied (see the module docstring)."""
    catalog = json.loads(CATALOG.read_text(encoding="utf-8"))
    raw = {e["key"]: e for e in catalog["entries"]}
    sources = catalog["sources"]
    out_entries = []
    for _mid, mat, members in material_groups(entries):
        for en in members:
            m, key = en["meta"], en["meta"]["key"]
            blob = json.loads((LIB_DIR / en["file"]).read_text(encoding="utf-8"))
            view = blob_view(blob)
            r = engine.resolve_entry(list(view["iso_map"]), m["components"], m["support"], key)
            iso_raw = {i["component"]: i for i in raw[key].get("isotherms") or []}
            isotherms = []
            for i in m["isotherms"]:
                comp = i["component"]
                tag = engine.isotherm_order[view["iso_map"][comp]]
                desc = engine.isotherm[tag]
                values = {p["name"]: p["value"] for p in i["params"]}
                sp = iso_raw[comp].get("source_points")
                isotherms.append({
                    "component": comp, "kind": tag, "kind_label": desc["label"],
                    "parameters": [{"name": p["name"], "label": p["label"], "unit": p["unit"], "value": values[p["name"]]}
                                   for p in desc["parameters"]],
                    "data_quality": i["data_quality"],
                    "data_quality_label": DATA_QUALITY_LABELS.get(i["data_quality"], i["data_quality"]),
                    "data_quality_note": i["data_quality_note"],
                    "ref": i["ref"],
                    "T_range_K": i["T_range_K"], "P_range_Pa": i["P_range_Pa"],
                    "T_range_not_stated": i["T_range_not_stated"], "P_range_not_stated": i["P_range_not_stated"],
                    "range_note": i["range_note"],
                    "source_points": resolved_source_points(sp, None) if sp else None,
                })
            props = []
            for p in m["properties"]:
                eng = r["properties"].get(p["name"])
                props.append({
                    "name": p["name"], "label": p["label"],
                    "unit": eng["unit"] if eng else p["unit"],
                    "registry_path": eng["registry_path"] if eng else None,
                    "engine_name": eng["engine_name"] if eng else None,
                    "value": p["value"], "provenance": p["provenance"],
                    "provenance_label": PROV.get(p["provenance"], p["provenance"]),
                    "ref": p["ref"], "temperature_table": p["temperature_table"],
                })
            out_entries.append({
                "key": key, "label": m["label"], "file": en["file"],
                "material_id": m["material_id"], "material_name": m["material_name"], "variant": m["variant"],
                "family": m["family"], "support": m["support"], "recommended": m["recommended"],
                "status": m["status"], "legacy_note": m["legacy_note"],
                "applications": m["applications"], "components": m["components"],
                "aliases": mat.get("aliases") or [],
                "sources": m["sources"],
                "isotherms": isotherms,
                "heat_of_adsorption": [{**h, "unit": r["heat_unit"], "provenance_label": PROV.get(h["provenance"], h["provenance"])}
                                       for h in m["heat_of_adsorption"]],
                "properties": props,
                "mass_transfer": r["mass_transfer"],
                "example_decks": m["example_decks"], "example_deck_note": m["example_deck_note"],
                "completeness": {
                    "complete": m["completeness"]["complete"], "missing": m["completeness"]["missing"],
                    "unverified_fields": raw[key].get("unverified_fields") or [],
                    "verification_note": m["completeness"]["verification_note"],
                },
                "notes": m["notes"],
            })
    doc = {
        "_about": "GENERATED by scripts/docs/gen_adsorbent_library.py from library.json and the <key>.json blobs; do not edit. "
                  "Resolved catalogue: each entry with the engine's labels, units and parameter order applied, source points in SI, "
                  "completeness and provenance, in the order of docs/adsorbent-library.md. Shape: library.resolved.schema.json.",
        "schema_version": RESOLVED_VERSION,
        "catalogue_schema_version": catalog["schema_version"],
        "materials": {mid: {"name": mat["name"], "family": mat["family"], "recommended": mat["recommended"],
                            "recommended_reason": mat["recommended_reason"], "aliases": mat["aliases"]}
                      for mid, mat in catalog["materials"].items()},
        "entries": out_entries,
    }
    return json.dumps(doc, indent=2, ensure_ascii=False) + "\n"


def validate_schema(value, schema, root, path="$"):
    """A small JSON Schema checker (type, enum, required, properties,
    additionalProperties, items, $ref into $defs, minimum): the subset
    library.resolved.schema.json uses. Returns a list of error strings."""
    if "$ref" in schema:
        node = root
        for part in schema["$ref"].lstrip("#/").split("/"):
            node = node[part]
        return validate_schema(value, node, root, path)
    errs = []
    types = schema.get("type")
    if types is not None:
        types = [types] if isinstance(types, str) else types
        ok = {"object": isinstance(value, dict), "array": isinstance(value, list), "string": isinstance(value, str),
              "boolean": isinstance(value, bool), "null": value is None,
              "integer": isinstance(value, int) and not isinstance(value, bool),
              "number": isinstance(value, (int, float)) and not isinstance(value, bool)}
        if not any(ok[t] for t in types):
            return [f"{path}: expected {'/'.join(types)}, got {type(value).__name__}"]
    if "enum" in schema and value not in schema["enum"]:
        errs.append(f"{path}: {value!r} not in {schema['enum']}")
    if "minimum" in schema and isinstance(value, (int, float)) and value < schema["minimum"]:
        errs.append(f"{path}: {value} < {schema['minimum']}")
    if isinstance(value, dict):
        props = schema.get("properties", {})
        for req in schema.get("required", []):
            if req not in value:
                errs.append(f"{path}: missing '{req}'")
        extra = schema.get("additionalProperties", True)
        for k, v in value.items():
            if k in props:
                errs += validate_schema(v, props[k], root, f"{path}.{k}")
            elif isinstance(extra, dict):
                errs += validate_schema(v, extra, root, f"{path}.{k}")
            elif extra is False:
                errs.append(f"{path}: unexpected '{k}'")
    if isinstance(value, list) and "items" in schema:
        for i, v in enumerate(value):
            errs += validate_schema(v, schema["items"], root, f"{path}[{i}]")
    return errs


def check_resolved_schema(text):
    schema = json.loads(RESOLVED_SCHEMA.read_text(encoding="utf-8"))
    return validate_schema(json.loads(text), schema, schema)


def completeness_summary(entries) -> str:
    inc = [e for e in entries if not e["meta"]["completeness"]["complete"]]
    pts = sum(1 for e in entries if e["points"])
    lines = [f"[adsorbent-library] completeness: {len(entries) - len(inc)}/{len(entries)} complete; {pts} with tabulated source data"]
    for e in inc:
        lines.append(f"  - incomplete {e['meta']['key']} ({e['meta']['status']}): {', '.join(e['meta']['completeness']['missing'])}")
    return "\n".join(lines)


def main(argv=None) -> int:
    ap = argparse.ArgumentParser(description=__doc__.splitlines()[0])
    ap.add_argument("--check", action="store_true", help="validate and fail if the docs page or the merged catalogue is out of date")
    args = ap.parse_args(argv)
    entries, errors = load_library()
    if errors:
        print(f"[adsorbent-library] {len(errors)} catalog error(s):", file=sys.stderr)
        for e in errors:
            print(f"  - {e}", file=sys.stderr)
        return 1
    outputs = [(DOCS_OUT, render_docs(entries)), (MERGED_OUT, render_merged(entries))]
    engine = load_engine()
    if engine is not None:
        resolved = render_resolved(entries, engine)
        problems = check_resolved_schema(resolved)
        if problems:
            print(f"[adsorbent-library] resolved metadata violates {RESOLVED_SCHEMA.name}:", file=sys.stderr)
            for p in problems[:20]:
                print(f"  - {p}", file=sys.stderr)
            return 1
        outputs.append((RESOLVED_OUT, resolved))
    if args.check:
        stale = [p for p, text in outputs
                 if (p.read_text(encoding="utf-8") if p.is_file() else None) != text]
        for p in stale:
            print(f"[adsorbent-library] DRIFT: {p.relative_to(ROOT)} is out of date — run "
                  "python scripts/docs/gen_adsorbent_library.py", file=sys.stderr)
        if engine is None:
            # No engine package: validate the committed resolved file against the
            # schema and the catalogue's entry list (full regeneration: python/tests).
            if not RESOLVED_OUT.is_file():
                print(f"[adsorbent-library] DRIFT: {RESOLVED_OUT.relative_to(ROOT)} is missing", file=sys.stderr)
                return 1
            committed = RESOLVED_OUT.read_text(encoding="utf-8")
            problems = check_resolved_schema(committed)
            keys = [e["key"] for e in json.loads(committed)["entries"]]
            want = [en["meta"]["key"] for _m, _mat, members in material_groups(entries) for en in members]
            if keys != want:
                problems.append(f"entries {keys} do not match the catalogue {want}")
            for p in problems[:20]:
                print(f"[adsorbent-library] {RESOLVED_OUT.name}: {p}", file=sys.stderr)
            if problems:
                return 1
            print("[adsorbent-library] skarstrom package not installed: resolved file checked against its schema "
                  "and the entry list only (python/tests/test_library_resolved.py compares a regeneration).")
        if stale:
            return 1
        print(f"[adsorbent-library] {len(entries)} entries; generated files up to date.")
    else:
        if engine is None:
            print("[adsorbent-library] the resolved file needs the skarstrom package (pip install skarstrom, "
                  "or pip install . from a checkout)", file=sys.stderr)
            return 1
        for p, text in outputs:
            p.write_text(text, encoding="utf-8")
            print(f"[adsorbent-library] wrote {p.relative_to(ROOT)} ({len(entries)} entries)")
    print(completeness_summary(entries))
    return 0


if __name__ == "__main__":
    sys.exit(main())
