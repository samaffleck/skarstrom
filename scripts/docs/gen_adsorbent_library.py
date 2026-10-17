#!/usr/bin/env python3
"""Validate the curated adsorbent library and generate its docs page.

The library has ONE source of truth:
  data/adsorbents/library.json   catalogue: citation + per-parameter provenance
  data/adsorbents/<key>.json     cereal single-Adsorbent blobs, named after the key

This script VALIDATES the catalogue against the blobs (every entry cited or
marked legacy, metadata values equal to the blob values, no orphan blob files)
and enforces the COMPLETENESS rule: for a sourced entry every required field
(DOI, isotherm fit range, heat of adsorption, particle density, heat capacity,
particle size, example decks) is either sourced/derived or carries a ref
beginning "Not stated in source (checked: ...)"; a field that could not be
checked against the source must be listed in the entry's `unverified_fields`.
Nothing may be silently missing. Tabulated source data (`source_points`, per
isotherm) are validated and summarised.

It then GENERATES docs/adsorbent-library.md (never edit it by hand).

    python scripts/docs/gen_adsorbent_library.py          # validate + write the page
    python scripts/docs/gen_adsorbent_library.py --check  # validate + fail on drift

Standard library only.
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

# adsim::IsothermType enum order (core/include/adsim/equations/isotherm.h) ->
# kind tag (adsim::IsothermKindString in adsim/equations/isotherm_kinds.h).
ISOTHERM_KINDS = [
    "single_site_langmuir_p", "dual_site_langmuir_p", "inert", "henry",
    "single_site_langmuir_c", "dual_site_langmuir_c",
    "extended_single_site_langmuir_p", "extended_dual_site_langmuir_p",
    "extended_single_site_langmuir_c", "extended_dual_site_langmuir_c",
    "anti_langmuir", "freundlich", "sips", "sips_2", "extended_sips_2",
    "dual_site_sips", "redlich_peterson", "toth", "toth_2", "arrhenius_sips",
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
}

FAMILIES = ["zeolite", "carbon", "MOF", "silica", "alumina", "amine", "other"]
STATUSES = ["sourced", "legacy"]
SCHEMA_VERSION = 2
PROVENANCE = ["source", "derived", "repo", "default", "not_stated", "unknown"]
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
    "bulk_density": ("Bed (bulk) density", "kg/m3", None),
    "skeletal_density": ("Skeletal density", "kg/m3", None),
    "heat_capacity": ("Heat capacity", "J/(kg K)", "Specific_Heat_Capacity"),
    "particle_voidage": ("Particle voidage", "-", "Particle_Void_Fraction"),
    "particle_diameter": ("Particle diameter", "m", None),
    "thermal_conductivity": ("Thermal conductivity", "W/(m K)", "Thermal_Conductivity"),
}
REQUIRED_PROPERTIES = ["particle_density", "heat_capacity", "particle_voidage", "particle_diameter", "thermal_conductivity"]
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
        if e.get("family") not in FAMILIES:
            err(k, f"family '{e.get('family')}' not in {'/'.join(FAMILIES)}")
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

        props = e.get("properties") or {}
        for p in REQUIRED_PROPERTIES:
            if p not in props:
                err(k, f"property {p} missing")
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
            if prov == "not_stated" and (rec.get("value") is not None or not not_stated(rec.get("ref"))):
                err(k, f"property {p}: not_stated needs value null and a ref beginning 'Not stated in source (checked: ...)'")
            if rec.get("value") is None and prov != "not_stated" and not (legacy and prov == "unknown") and p not in unverified:
                err(k, f"property {p}: a null value needs provenance 'not_stated'")
            if d[2]:
                bv = b["scalars"][d[2]]
                if rec.get("value") is None or not rel_eq(rec["value"], bv):
                    err(k, f"property {p}: catalog {rec.get('value')} != blob {d[2]} {bv}")
            props_out.append({"name": p, "label": d[0], "unit": d[1], "value": rec.get("value"),
                              "provenance": prov, "ref": rec.get("ref")})
        for p in ("particle_density", "heat_capacity", "particle_diameter"):
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
            "family": e.get("family"),
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
        if f.suffix == ".json" and f.name != "library.json" and f.name not in seen_files:
            errors.append(f"orphan blob '{f.name}' is not listed in library.json")
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
        "not_stated": "not stated", "unknown": "unknown"}


def anchor(key: str) -> str:
    return key.replace("_", "-")


def render_docs(entries) -> str:
    L: list = []
    out = L.append
    n_sourced = sum(1 for e in entries if e["meta"]["status"] == "sourced")
    out("<!-- GENERATED by scripts/docs/gen_adsorbent_library.py from data/adsorbents/library.json — do not edit. -->")
    out("")
    out("# Adsorbent Library")
    out("")
    out(f"Skarstrom ships a curated library of **{len(entries)}** adsorbents ({n_sourced} sourced, {len(entries) - n_sourced} legacy) in `data/adsorbents/`. Every sourced entry cites the publication (DOI and table/figure) its parameters were transcribed from, and every property records its provenance:")
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
    complete = sum(1 for e in entries if e["meta"]["completeness"]["complete"])
    with_points = sum(1 for e in entries if any(i["source_points"] for i in e["meta"]["isotherms"]))
    out(f"**Completeness:** {complete} of {len(entries)} entries carry every required field (DOI, isotherm fit range, heat of adsorption, particle density, heat capacity, particle size, example decks) as a sourced value or an explicit *not stated in source* note; the others are listed as incomplete below. **Tabulated source data:** {with_points} entries include equilibrium points tabulated in their source (tables, SI files or datasets — never digitised from plots), stored with the entry in `library.json`.")
    out("")
    out("## Using the library")
    out("")
    out("Each entry is a catalogue record in `library.json` plus one serialised adsorbent, `<key>.json`. The Python package ships the library:")
    out("")
    out("```python")
    out("import skarstrom")
    out("lib = skarstrom.adsorbent_library()               # the catalogue as a dict")
    out('[e["key"] for e in lib["entries"]]')
    out('path = skarstrom.adsorbent_path("zeolite_13x_perez2019")')
    out("```")
    out("")
    out("In C++, `adsim/objects/adsorbent_library.h` reads the catalogue, `LoadAdsorbentFromJSON` loads an entry's file, and `AddLibraryAdsorbent` (`adsim/objects/deck_edit.h`) adds an entry to a deck with default physics for any component the entry does not cover. The library is curated: see `CONTRIBUTING.md` for the data-integrity rules a new entry must meet (every parameter cited).")
    out("")
    out("!!! warning \"Before you run\"")
    out("    Check the heat of adsorption (some sources report only site energies, so the value is left at 0) and the heat-of-adsorption convention of your deck (isosteric vs internal energy) against the source. Fit ranges are listed per component; outside them the isotherm is an extrapolation.")
    out("")
    out("| Material | Family | Components | Applications | Source |")
    out("|---|---|---|---|---|")
    for e in entries:
        m = e["meta"]
        src = "source unknown (legacy)" if m["status"] == "legacy" else ", ".join(f"[{md(s['doi'])}]({s['url']})" for s in m["sources"])
        out(f"| [{md(m['label'])}](#{anchor(m['key'])}) | {md(m['family'])} | {', '.join(md(c) for c in m['components'])} | {', '.join(md(a) for a in m['applications'])} | {src} |")
    out("")
    for e in entries:
        m, f = e["meta"], e["file"]
        out("---")
        out("")
        out(f"### {md(m['label'])} {{#{anchor(m['key'])}}}")
        out("")
        out(f"**Material:** {md(m['material'])} · **Family:** {md(m['family'])} · **Applications:** {', '.join(md(a) for a in m['applications'])} · **Library key:** `{m['key']}` · **File:** `{f}`")
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
        out("| Component | Isotherm | Parameters | T range of fit | P range of fit | Ref. |")
        out("|---|---|---|---|---|---|")
        for i in m["isotherms"]:
            params = "; ".join(f"{md(p['name'])} = {fmt_num(p['value'])}" for p in i["params"])
            out(f"| {md(i['component'])} | {md(i['kind_label'])} | {params} | {fmt_range(i['T_range_K'], 1, ' K')} | {fmt_range(i['P_range_Pa'], 1e-5, ' bar')} | {md(i['ref'])} |")
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
        sp = [i for i in m["isotherms"] if i["source_points"]]
        if sp:
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
    return "\n".join(L)


def completeness_summary(entries) -> str:
    inc = [e for e in entries if not e["meta"]["completeness"]["complete"]]
    pts = sum(1 for e in entries if e["points"])
    lines = [f"[adsorbent-library] completeness: {len(entries) - len(inc)}/{len(entries)} complete; {pts} with tabulated source data"]
    for e in inc:
        lines.append(f"  - incomplete {e['meta']['key']} ({e['meta']['status']}): {', '.join(e['meta']['completeness']['missing'])}")
    return "\n".join(lines)


def main(argv=None) -> int:
    ap = argparse.ArgumentParser(description=__doc__.splitlines()[0])
    ap.add_argument("--check", action="store_true", help="validate and fail if the docs page is out of date")
    args = ap.parse_args(argv)
    entries, errors = load_library()
    if errors:
        print(f"[adsorbent-library] {len(errors)} catalog error(s):", file=sys.stderr)
        for e in errors:
            print(f"  - {e}", file=sys.stderr)
        return 1
    text = render_docs(entries)
    if args.check:
        cur = DOCS_OUT.read_text(encoding="utf-8") if DOCS_OUT.is_file() else None
        if cur != text:
            print(f"[adsorbent-library] DRIFT: {DOCS_OUT.relative_to(ROOT)} is out of date — run "
                  "python scripts/docs/gen_adsorbent_library.py", file=sys.stderr)
            return 1
        print(f"[adsorbent-library] {len(entries)} entries; generated page up to date.")
    else:
        DOCS_OUT.write_text(text, encoding="utf-8")
        print(f"[adsorbent-library] wrote {DOCS_OUT.relative_to(ROOT)} ({len(entries)} entries)")
    print(completeness_summary(entries))
    return 0


if __name__ == "__main__":
    sys.exit(main())
