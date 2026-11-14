#!/usr/bin/env python3
# Skarstrom - pressure and temperature swing adsorption simulator
# Copyright (c) 2024-2026 Sam Affleck and contributors
# Licensed under the terms in the LICENSE file at the repository root.
"""Generate the pure-gas transport tables of the engine:

    core/include/adsim/objects/viscosity_coefficients.h
    core/include/adsim/objects/thermal_conductivity_coefficients.h

Each component has one degree-4 polynomial in T for its dilute-gas viscosity
[Pa s] and one for its dilute-gas thermal conductivity [W/(m K)], for the
window [T_MIN, T_MAX] of Fluid (kViscTmin, kViscTmax). Each row has one of
these sources:

  NASA      The transport database of the NASA CEA program, data/trans.inp
            (R. A. Svehla, Transport Coefficients for the NASA Lewis Chemical
            Equilibrium Program, NASA TM-4647, 1995, and the primary sources
            that each species names). The form is
                ln(eta) = A ln T + B/T + C/T^2 + D
            with eta in micropoise (viscosity) or microwatt/(cm K)
            (conductivity). The script takes the interval that starts at the
            lowest temperature. Below the lower limit T0 of that interval (200 K
            or 300 K), it uses a kinetic-theory shape E(T) (Chapman-Enskog for
            the viscosity, Eucken for the conductivity) that joins the NASA
            curve N(T) smoothly: eta = N(T0) (E(T)/E(T0))^r, with r the ratio of
            the log-slopes d ln N/d ln T and d ln E/d ln T at T0 (0.7-1.2 for
            these components). The polynomial is a least-squares fit, in
            relative error, to the result on a 1 K grid over the window.
  CoolProp  A polynomial fitted before this script existed to CoolProp's
            dilute-gas value. The script carries the coefficients unchanged
            (CARRIED_VISCOSITY, CARRIED_CONDUCTIVITY), with the BibTeX key of
            the CoolProp 6.7.0 correlation, so that the header names the
            primary source. A refit with CoolProp 6.7.0 does not reproduce
            them exactly (0.2-3 %), so they are not refitted.
  Eucken    k = mu (cp/M + 1.25 R/M), with mu the viscosity row of the same
            component and cp the NASA-7 polynomial of the engine
            (ideal_gas_heat_capacity.h, continued below its lower limit). Only
            for a component that neither NASA nor CoolProp has a conductivity
            for.

The NASA rows replace estimates (Chapman-Enskog viscosity with the Lennard-Jones
parameters of Svehla 1962, and the Eucken conductivity). A component in a
CoolProp row keeps it also when NASA has the species: those CoolProp
correlations are more recent, dedicated reference correlations.

trans.inp is pinned: CEA_COMMIT names the nasa/cea commit and TRANS_SHA256 the
SHA-256 of the file. By default the script downloads the file into a cache
directory; --trans reads a local copy instead. A file whose SHA-256 differs is
an error, so the output depends only on the pinned file and the engine headers
that the script reads (component list, Lennard-Jones table, molar masses, NASA-7
heat capacities).

Usage:
    python scripts/gen_transport_coefficients.py              # write the headers
    python scripts/gen_transport_coefficients.py --check      # fail if a committed header differs
    python scripts/gen_transport_coefficients.py --trans FILE # local trans.inp (no download)
"""
from __future__ import annotations

import argparse
import difflib
import hashlib
import math
import os
import re
import sys
import tempfile
import urllib.request
from pathlib import Path

ROOT = Path(__file__).resolve().parents[1]
OBJ = ROOT / "core/include/adsim/objects"
FLUID_H = OBJ / "fluid.h"
PROPS_H = OBJ / "component_properties.h"
CP_H = OBJ / "ideal_gas_heat_capacity.h"
OUT_VISC = OBJ / "viscosity_coefficients.h"
OUT_COND = OBJ / "thermal_conductivity_coefficients.h"

CEA_COMMIT = "2f79a647e85737742fb08a127da00d1ea1d069da"
TRANS_URL = f"https://raw.githubusercontent.com/nasa/cea/{CEA_COMMIT}/data/trans.inp"
TRANS_SHA256 = "379c514a7f7638371d8d9254a1d653cd078084a9371061aa0b45c5a5c3334d41"

T_MIN, T_MAX = 150.0, 650.0     # Fluid::kViscTmin, Fluid::kViscTmax
R = 8.314462618                  # [J/(mol K)]

# Component -> trans.inp species. NASA gives each row of these components that
# has no CoolProp entry (for HydrogenSulfide, only the conductivity).
NASA_SPECIES = {
    "CarbonMonoxide": "CO",
    "Ethylene": "C2H4",
    "Fluorine": "F2",
    "HydrogenSulfide": "H2S",
    "Krypton": "Kr",
    "Neon": "Ne",
    "NitrousOxide": "N2O",
    "SulfurDioxide": "SO2",
    "Xenon": "Xe",
}

# Largest relative error of a fitted polynomial: where its source has data
# (the NASA interval, or the full window), and in the band below the lower
# limit of a NASA interval (the kinetic-theory shape).
MAX_FIT_ERROR = 0.01
MAX_SHAPE_ERROR = 0.03

# CoolProp rows: (a0..a4 as written in the header, BibTeX key of the CoolProp
# 6.7.0 correlation; a note in brackets where it is not a reference correlation).
CARRIED_VISCOSITY = {
    "Air": ('-3.975668308e-07', '8.355550932e-08', '-8.946761945e-11', '8.105793652e-14', '-3.291047444e-17',
        "Lemmon-IJT-2004"),
    "Ammonia": ('5.662838731e-06', '-2.243235441e-08', '1.897691980e-10', '-2.536813623e-13', '1.230607119e-16',
        "Fenghour-JPCRD-1995"),
    "Argon": ('-8.761247680e-07', '9.881118067e-08', '-8.412117788e-11', '6.410174481e-14', '-2.296012007e-17',
        "Lemmon-IJT-2004"),
    "Benzene": ('1.254467318e-06', '1.638907665e-08', '2.301358987e-11', '-2.659481940e-14', '1.084363444e-17',
        "Avgeri-JPCRD-2014-Benzene"),
    "Butane": ('1.156191447e-07', '2.297894162e-08', '1.056158257e-11', '-2.219275149e-14', '1.244823019e-17',
        "Vogel-HTHP-1999"),
    "CarbonDioxide": ('1.308709777e-07', '4.820694113e-08', '1.759231875e-11', '-5.203688805e-14', '2.874773551e-17',
        "Laesecke-JPCRD-2017-CO2"),
    "Cyclohexane": ('4.473278869e-06', '-9.547866920e-09', '8.972370936e-11', '-1.123234660e-13', '4.951796980e-17',
        "Tariq-JPCRD-2014-Cyclohexane"),
    "Ethane": ('-3.428385448e-08', '3.184711463e-08', '6.571797844e-12', '-3.240126156e-14', '2.075877142e-17',
        "Friend-JPCRD-1991"),
    "Ethanol": ('-1.059293726e-06', '3.562551730e-08', '-1.129025515e-11', '1.013812599e-14', '-6.916014690e-18',
        "Kiselev-IECR-2005"),
    "Helium": ('3.258489119e-06', '7.019132557e-08', '-6.643826261e-11', '6.845952015e-14', '-3.105348082e-17',
        "Arp-NIST-1998"),
    "Hydrogen": ('1.063502253e-06', '3.469806338e-08', '-3.932579965e-11', '4.337354026e-14', '-2.025121124e-17',
        "Muzny-JCED-2013"),
    "HydrogenSulfide": ('2.570644577e-06', '7.734227358e-09', '1.359101931e-10', '-2.202066452e-13', '1.187091121e-16',
        "QuinonesCisneros-JCED-2012"),
    "Isobutane": ('-1.114114530e-07', '2.646426604e-08', '-1.489597263e-12', '-6.509845136e-15', '4.182420795e-18',
        "Vogel-IJT-2000"),
    "Isopentane": ('7.793850349e-07', '1.309729644e-08', '4.315222167e-11', '-6.834421373e-14', '3.426868483e-17',
        "Chung-IECR-1988 (Chung estimate)"),
    "Methane": ('-1.113330524e-06', '5.500995567e-08', '-6.266647772e-11', '6.466274551e-14', '-2.931415673e-17',
        "QuinonesCisneros-JPCB-2006"),
    "Methanol": ('8.255187980e-07', '2.844654012e-08', '-3.736111742e-12', '3.362927905e-14', '-3.085984234e-17',
        "Xiang-JPCRD-2006"),
    "Nitrogen": ('-2.489909434e-07', '8.016074213e-08', '-8.622093890e-11', '7.807414792e-14', '-3.158186321e-17',
        "Lemmon-IJT-2004"),
    "Octane": ('-1.079308585e-06', '3.187610322e-08', '-4.949795456e-11', '7.746855436e-14', '-4.370967753e-17',
        "Huber-FPE-2004"),
    "Oxygen": ('-6.424959437e-07', '9.219506868e-08', '-9.156634332e-11', '7.888483229e-14', '-3.102898027e-17',
        "Lemmon-IJT-2004"),
    "Pentane": ('-3.193560173e-07', '2.589662993e-08', '-1.395470524e-11', '2.753700815e-14', '-2.045304206e-17',
        "QuinonesCisneros-JPCB-2006"),
    "Propane": ('9.477457721e-07', '1.723091423e-08', '4.052645690e-11', '-6.898194921e-14', '3.552897634e-17',
        "Vogel-JPCRD-1998"),
    "Propylene": ('1.040721774e-06', '1.361304903e-08', '6.932026272e-11', '-1.227180457e-13', '6.835768150e-17',
        "Huber-IECR-2003 (extended corresponding states)"),
    "SulfurHexafluoride": ('8.115148060e-07', '4.938162982e-08', '3.676082818e-12', '-2.907000788e-14', '1.630163748e-17',
        "QuinonesCisneros-JPCRD-2012"),
    "Toluene": ('3.460682436e-07', '2.562857369e-08', '-2.334484737e-11', '5.045965298e-14', '-3.465751228e-17',
        "Avgeri-JPCRD-2015-Toluene"),
    "Water": ('6.685582872e-06', '-2.114915253e-08', '1.421511602e-10', '-1.395610837e-13', '4.941363610e-17',
        "Huber-JPCRD-2009"),
}

CARRIED_CONDUCTIVITY = {
    "Air": ('-1.103689343e-03', '1.150872000e-04', '-1.017278775e-07', '8.892823023e-11', '-3.489346002e-14',
        "Lemmon-IJT-2004"),
    "Ammonia": ('3.364189134e-02', '-1.526062879e-04', '3.784640531e-07', '2.804148391e-10', '-5.422946728e-13',
        "Tufeu-BBPC-1984"),
    "Argon": ('-7.125512209e-04', '7.752521344e-05', '-6.538144016e-08', '4.867387913e-11', '-1.701578264e-14',
        "Lemmon-IJT-2004"),
    "Benzene": ('1.298818600e-02', '-1.453267222e-04', '6.602809907e-07', '-7.751658646e-10', '3.470083372e-13',
        "Assael-JPCRD-2012-Benzene"),
    "Butane": ('2.171007845e-03', '-3.984701616e-06', '1.842877141e-07', '-3.543541947e-11', '1.860946374e-14',
        "Perkins-JCED-2002-nButane"),
    "CarbonDioxide": ('9.416515217e-04', '1.366165474e-05', '1.900929996e-07', '-2.284300458e-10', '9.819116439e-14',
        "Huber-JPCRD-2016-CO2"),
    "Ethane": ('1.284237111e-03', '9.409162784e-06', '1.599536977e-07', '1.661591530e-10', '-2.083007400e-13',
        "Friend-JPCRD-1991"),
    "Ethanol": ('1.183679899e-02', '-1.147692558e-04', '5.990535289e-07', '-6.806635415e-10', '3.033441519e-13',
        "Assael-JPCRD-2013-Ethanol"),
    "Helium": ('2.277101233e-02', '5.654886831e-04', '-5.498894004e-07', '5.592509200e-10', '-2.506192097e-13',
        "Hands-CRYO-1981"),
    "Hydrogen": ('-2.815427557e-02', '1.071286869e-03', '-1.651895725e-06', '1.779587551e-09', '-7.351566094e-13',
        "Assael-JPCRD-2011-Hydrogen"),
    "Isobutane": ('-1.916044286e-03', '2.117536815e-05', '1.457239015e-07', '-1.832318819e-11', '6.720988992e-15',
        "Perkins-JCED-2002-Isobutane"),
    "Isopentane": ('1.780281893e-03', '-2.285042130e-05', '2.226711585e-07', '9.485397529e-12', '-7.823995723e-14',
        "Vassiliou-JPCRD-2015-pentanes"),
    "Methane": ('-1.342352025e-03', '1.309864327e-04', '-1.853001437e-07', '6.135535740e-10', '-4.171245103e-13',
        "Friend-JPCRD-1989"),
    "Methanol": ('7.149426974e-03', '-4.860312275e-05', '3.068597586e-07', '-1.962327540e-10', '5.129961260e-14',
        "Sykioti-JPCRD-2013-Methanol"),
    "Nitrogen": ('-1.140235134e-03', '1.150931704e-04', '-1.074526406e-07', '9.459915977e-11', '-3.721974272e-14',
        "Lemmon-IJT-2004"),
    "Octane": ('1.514313723e-02', '-1.412034281e-04', '5.769547248e-07', '-5.727057956e-10', '2.249106792e-13',
        "Huber-FPE-2005"),
    "Oxygen": ('-1.337308378e-03', '1.110772666e-04', '-7.704117649e-08', '5.946003230e-11', '-2.125278802e-14',
        "Lemmon-IJT-2004"),
    "Pentane": ('6.083391917e-03', '-5.726343072e-05', '3.451353097e-07', '-2.059411437e-10', '3.085914190e-14',
        "Vassiliou-JPCRD-2015-pentanes"),
    "Propane": ('-6.332823724e-04', '1.426336110e-05', '1.779276446e-07', '-5.130067726e-11', '2.878203637e-14',
        "Marsh-JCED-2002"),
    "Propylene": ('5.618233006e-03', '-4.293905616e-05', '3.428177340e-07', '-2.533491288e-10', '7.086168179e-14',
        "Huber-IECR-2003 (extended corresponding states)"),
    "SulfurHexafluoride": ('-2.893111203e-03', '9.567804820e-06', '2.379813880e-07', '-3.598758976e-10', '1.800158661e-13',
        "Assael-JPCRD-2012-SF6"),
    "Toluene": ('2.437868018e-03', '-2.370139831e-05', '1.914645592e-07', '-5.045649611e-11', '-3.937756297e-14',
        "Assael-JPCRD-2012-Toluene"),
    "Water": ('1.992543704e-02', '-1.304602429e-04', '6.227676609e-07', '-8.012497117e-10', '4.136292535e-13',
        "Huber-JPCRD-2012"),
}


# ---------------------------------------------------------------------------
# Inputs: engine headers and trans.inp
# ---------------------------------------------------------------------------

def component_names() -> list[str]:
    """Components::componentNames, read from fluid.h."""
    text = FLUID_H.read_text(encoding="utf-8")
    m = re.search(r"componentNames\s*=\s*\{(.*?)\};", text, re.S)
    if not m:
        sys.exit(f"cannot find Components::componentNames in {FLUID_H}")
    return re.findall(r'"([^"]+)"', m.group(1))


def check_window() -> None:
    """T_MIN and T_MAX must be the evaluation window of Fluid."""
    text = FLUID_H.read_text(encoding="utf-8")
    for name, value in (("kViscTmin", T_MIN), ("kViscTmax", T_MAX)):
        m = re.search(rf"{name}\s*=\s*([0-9.eE+-]+)\s*;", text)
        if not m or float(m.group(1)) != value:
            sys.exit(f"Fluid::{name} in {FLUID_H} is not {value}: change T_MIN/T_MAX")


def lennard_jones() -> dict[str, tuple[float, float]]:
    """Components::componentData: name -> (eps/k [K], sigma [Angstrom])."""
    text = FLUID_H.read_text(encoding="utf-8")
    return {m.group(1): (float(m.group(2)), float(m.group(3)))
            for m in re.finditer(r'\{\s*"(\w+)",\s*\{\s*([0-9.]+),\s*([0-9.]+)\s*[,}]', text)}


def molar_masses() -> dict[str, float]:
    """kComponentProperties: name -> molar mass [kg/mol]."""
    text = PROPS_H.read_text(encoding="utf-8")
    return {m.group(1): float(m.group(2))
            for m in re.finditer(r'\{"(\w+)", "[^"]*", "[^"]*", "[^"]*",\s*"[^"]*",\s*([0-9.eE+-]+),', text)}


def nasa7() -> dict[str, tuple[float, float, list[float]]]:
    """kIdealGasHeatCapacity: name -> (t_min, t_max, a0..a4)."""
    text = CP_H.read_text(encoding="utf-8")
    out = {}
    for m in re.finditer(r'\{"(\w+)", \{([0-9.]+), ([0-9.]+),\s*\{([^}]*)\}', text):
        a = [float(x) for x in m.group(4).split(",")]
        out[m.group(1)] = (float(m.group(2)), float(m.group(3)), a[:5])
    return out


def trans_file(path: Path | None) -> Path:
    """The pinned trans.inp: `path`, or downloaded into a cache directory."""
    if path is None:
        cache = Path(os.environ.get("SKARSTROM_CEA_CACHE",
                                    Path(tempfile.gettempdir()) / f"skarstrom-cea-{CEA_COMMIT[:12]}"))
        cache.mkdir(parents=True, exist_ok=True)
        path = cache / "trans.inp"
        if not (path.exists() and hashlib.sha256(path.read_bytes()).hexdigest() == TRANS_SHA256):
            try:
                with urllib.request.urlopen(TRANS_URL, timeout=60) as r:
                    path.write_bytes(r.read())
            except OSError as e:
                sys.exit(f"cannot download {TRANS_URL}: {e} (or pass --trans FILE)")
    if not path.exists():
        sys.exit(f"{path} does not exist")
    got = hashlib.sha256(path.read_bytes()).hexdigest()
    if got != TRANS_SHA256:
        sys.exit(f"{path}: SHA-256 {got} is not the pinned {TRANS_SHA256} (nasa/cea {CEA_COMMIT[:12]})")
    return path


def fortran_float(s: str) -> float:
    """'0.61205763E 00' (a blank exponent sign) -> 0.61205763."""
    return float(s.replace("E ", "E+"))


def read_trans(path: Path) -> dict[str, dict]:
    """Pure species of trans.inp: name -> {'source', 'V': [...], 'C': [...]},
    each interval (t0, t1, (A, B, C, D)). Records for a pair of species (a
    second name in columns 17-32) are skipped."""
    species: dict[str, dict] = {}
    current = None
    for i, line in enumerate(path.read_text(encoding="ascii").splitlines()):
        if i == 0 or line.strip() in ("", "end"):
            continue
        if line[0] != " ":
            name, partner = line[:16].strip(), line[16:32].strip()
            current = None if partner else name
            if current:
                species[current] = {"source": " ".join(line[38:].split()), "V": [], "C": []}
        elif current:
            kind = line[1]
            t0, t1 = float(line[2:11]), float(line[11:20])
            coeffs = tuple(fortran_float(line[20 + 15 * j:35 + 15 * j]) for j in range(4))
            species[current][kind].append((t0, t1, coeffs))
    return species


# ---------------------------------------------------------------------------
# Property models
# ---------------------------------------------------------------------------

def nasa_value(interval, T: float) -> float:
    """The CEA form, in micropoise or microwatt/(cm K)."""
    A, B, C, D = interval[2]
    return math.exp(A * math.log(T) + B / T + C / (T * T) + D)


def nasa_log_slope(interval, T: float) -> float:
    """d ln(eta) / d ln T of the CEA form."""
    A, B, C, _ = interval[2]
    return A - B / T - 2.0 * C / (T * T)


def log_slope(f, T: float, h: float = 1e-4) -> float:
    """d ln f / d ln T, central difference."""
    return (math.log(f(T * (1.0 + h))) - math.log(f(T * (1.0 - h)))) / (math.log(1.0 + h) - math.log(1.0 - h))


def joined(interval, scale: float, shape):
    """The NASA value (times `scale`, to SI) at T >= T0, and below T0 the
    shape that meets it in value and log-slope (see the module docstring)."""
    t0 = interval[0]
    n0 = nasa_value(interval, t0) * scale
    e0 = shape(t0)
    r = nasa_log_slope(interval, t0) / log_slope(shape, t0)

    def f(T):
        if T >= t0:
            return nasa_value(interval, T) * scale
        return n0 * (shape(T) / e0) ** r
    return f, t0


def chapman_enskog_viscosity(eps_k: float, sigma: float, M: float, T: float) -> float:
    """Dilute-gas viscosity [Pa s], with the Neufeld et al. (1972) collision
    integral Omega(2,2) (Poling, Prausnitz and O'Connell, eq. 9-4.3)."""
    ts = T / eps_k
    omega = 1.16145 * ts ** -0.14874 + 0.52487 * math.exp(-0.77320 * ts) + 2.16178 * math.exp(-2.43787 * ts)
    return 26.69e-7 * math.sqrt(M * 1000.0 * T) / (sigma * sigma * omega)


def cp_molar(p: tuple[float, float, list[float]], T: float) -> float:
    """NASA-7 C_p0 [J/(mol K)]. Below t_min (200 K or 300 K) the polynomial
    continues, so that the Eucken shape has no kink at t_min. (The engine holds
    C_p0 at t_min there; this value only shapes a fit below the data.)"""
    t_min, t_max, a = p
    Tc = min(T, t_max)
    return R * (a[0] + Tc * (a[1] + Tc * (a[2] + Tc * (a[3] + Tc * a[4]))))


def eucken(mu: float, cp: float, M: float) -> float:
    """k = mu (cp/M + 1.25 R/M) [W/(m K)]."""
    return mu * (cp + 1.25 * R) / M


def poly(a, T: float) -> float:
    return a[0] + T * (a[1] + T * (a[2] + T * (a[3] + T * a[4])))


# ---------------------------------------------------------------------------
# Fit
# ---------------------------------------------------------------------------

GRID = [T_MIN + i for i in range(int(T_MAX - T_MIN) + 1)]
T_MID, T_HALF = 0.5 * (T_MIN + T_MAX), 0.5 * (T_MAX - T_MIN)


def fit_quartic(f, t0: float = T_MIN) -> tuple[list[float], float, float]:
    """Least squares in relative error of a0 + a1 T + ... + a4 T^4 to f on the
    1 K grid. Solved in x = (T - T_MID)/T_HALF, then expanded to powers of T.
    Returns the coefficients and the largest relative error at T >= t0 and at
    T < t0."""
    rows = []
    for T in GRID:
        y = f(T)
        x = (T - T_MID) / T_HALF
        rows.append(([x ** j / y for j in range(5)], 1.0))
    n = 5
    a = [[sum(r[0][i] * r[0][j] for r in rows) for j in range(n)] for i in range(n)]
    b = [sum(r[0][i] * r[1] for r in rows) for i in range(n)]
    for c in range(n):                                   # Gaussian elimination, partial pivoting
        p = max(range(c, n), key=lambda r: abs(a[r][c]))
        a[c], a[p], b[c], b[p] = a[p], a[c], b[p], b[c]
        for r in range(c + 1, n):
            m = a[r][c] / a[c][c]
            for k in range(c, n):
                a[r][k] -= m * a[c][k]
            b[r] -= m * b[c]
    s = [0.0] * n
    for r in reversed(range(n)):
        s[r] = (b[r] - sum(a[r][k] * s[k] for k in range(r + 1, n))) / a[r][r]
    # sum_j s_j ((T - T_MID)/T_HALF)^j  ->  sum_k c_k T^k
    coeffs = [0.0] * n
    for j in range(n):
        for k in range(j + 1):
            coeffs[k] += s[j] * math.comb(j, k) * (-T_MID) ** (j - k) / T_HALF ** j
    errs = [(abs(poly(coeffs, T) / f(T) - 1.0), T) for T in GRID]
    err_data = max(e for e, T in errs if T >= t0)
    err_low = max((e for e, T in errs if T < t0), default=0.0)
    return coeffs, err_data, err_low


def check_fit(name: str, what: str, t0: float, err_data: float, err_low: float, report: list[str]) -> None:
    if err_data > MAX_FIT_ERROR:
        sys.exit(f"{name}: {what} fit error {err_data:.2%} at T >= {t0:.0f} K exceeds {MAX_FIT_ERROR:.0%}")
    if err_low > MAX_SHAPE_ERROR:
        sys.exit(f"{name}: {what} fit error {err_low:.2%} below {t0:.0f} K exceeds {MAX_SHAPE_ERROR:.0%}")
    low = f", below {t0:.0f} K {err_low:.2%}" if t0 > T_MIN else ""
    report.append(f"{name:16} {what:12} fit error {err_data:.2%}{low}")


# ---------------------------------------------------------------------------
# Rows
# ---------------------------------------------------------------------------

def build(trans: Path) -> tuple[list[tuple], list[tuple], list[str]]:
    names = component_names()
    check_window()
    lj, mm, cp = lennard_jones(), molar_masses(), nasa7()
    db = read_trans(trans)
    visc, cond, report = [], [], []
    mu_fn = {}                                           # name -> viscosity [Pa s] (the row's function)

    for n in names:
        if n in CARRIED_VISCOSITY:
            *a, key = CARRIED_VISCOSITY[n]
            vals = [float(x) for x in a]
            mu_fn[n] = lambda T, v=vals: poly(v, T)
            visc.append((n, list(a), f"CoolProp: {key}"))
        elif n in NASA_SPECIES:
            sp = db[NASA_SPECIES[n]]
            iv = min(sp["V"], key=lambda r: r[0])
            eps, sig = lj[n]
            f, t0 = joined(iv, 1e-7,                    # micropoise -> Pa s
                           lambda T, eps=eps, sig=sig, M=mm[n]: chapman_enskog_viscosity(eps, sig, M, T))
            c, err, err_low = fit_quartic(f, t0)
            check_fit(n, "viscosity", t0, err, err_low, report)
            mu_fn[n] = f
            note = f"NASA {NASA_SPECIES[n]}, {sp['source']}" + (f"; Chapman-Enskog shape below {t0:.0f} K" if t0 > T_MIN else "")
            visc.append((n, [f"{x:.9e}" for x in c], note))
        else:
            sys.exit(f"{n}: no viscosity source (add it to CARRIED_VISCOSITY or NASA_SPECIES)")

    for n in names:
        if n in CARRIED_CONDUCTIVITY:
            *a, key = CARRIED_CONDUCTIVITY[n]
            cond.append((n, list(a), f"CoolProp: {key}"))
            continue
        if n not in cp:
            sys.exit(f"{n}: no NASA-7 heat capacity for the Eucken relation")
        M, mu, p = mm[n], mu_fn[n], cp[n]

        def eu(T, M=M, mu=mu, p=p):
            return eucken(mu(T), cp_molar(p, T), M)
        if n in NASA_SPECIES:
            sp = db[NASA_SPECIES[n]]
            iv = min(sp["C"], key=lambda r: r[0])
            f, t0 = joined(iv, 1e-4, eu)                 # microwatt/(cm K) -> W/(m K)
            note = f"NASA {NASA_SPECIES[n]}, {sp['source']}" + (f"; Eucken shape below {t0:.0f} K" if t0 > T_MIN else "")
        else:
            f, t0 = eu, T_MIN
            note = "Eucken (no NASA or CoolProp conductivity)"
        c, err, err_low = fit_quartic(f, t0)
        check_fit(n, "conductivity", t0, err, err_low, report)
        cond.append((n, [f"{x:.9e}" for x in c], note))
    return visc, cond, report


# ---------------------------------------------------------------------------
# Output
# ---------------------------------------------------------------------------

LICENSE = """/*
 * Skarstrom - pressure and temperature swing adsorption simulator
 * Copyright (c) 2024-2026 Sam Affleck and contributors
 *
 * Licensed under the terms in the LICENSE file at the repository root.
 */
"""

SOURCES = f"""//   NASA      The transport database of NASA CEA (data/trans.inp at nasa/cea
//             commit {CEA_COMMIT[:12]}, Apache-2.0): R. A. Svehla, Transport
//             Coefficients for the NASA Lewis Chemical Equilibrium Program, NASA
//             TM-4647 (1995), with the primary source that each species names
//             (in the row comment). Below the lower limit of the NASA interval
//             (200 K or 300 K), a kinetic-theory shape (Chapman-Enskog for mu,
//             Eucken for k) that joins the NASA curve in value and slope.
//   CoolProp  A fit to CoolProp's dilute-gas value. The row comment gives the
//             BibTeX key of the CoolProp 6.7.0 correlation (the primary source).
//   Eucken    k = mu (cp/M + 1.25 R/M), with the viscosity row of the same
//             component and the NASA-7 heat capacity of ideal_gas_heat_capacity.h.
//             Only where NASA and CoolProp have no conductivity."""


def header(kind: str, rows: list[tuple]) -> str:
    if kind == "V":
        title, var, unit, short = "VISCOSITY", "kViscosityCoefficients", "Pa-s", "mu"
        at300 = lambda a: f"mu300={poly([float(x) for x in a], 300.0) * 1e6:.2f} uPa.s"
        floor = "1.0e-5 Pa-s (Fluid::GetComponentViscosity)"
    else:
        title, var, unit, short = "THERMAL CONDUCTIVITY", "kThermalConductivityCoefficients", "W/(m.K)", "k"
        at300 = lambda a: f"k300={poly([float(x) for x in a], 300.0) * 1e3:.2f} mW/m.K"
        floor = "0.024 W/(m.K) (Fluid::GetComponentThermalConductivity)"
    body = "\n".join(f'    {{"{n}", {{{{{", ".join(a)}}}}}}}, // {at300(a)}; {note}' for n, a, note in rows)
    return f"""{LICENSE}#pragma once

// GENERATED by scripts/gen_transport_coefficients.py. Do not edit.

#include <array>
#include <string>
#include <unordered_map>

// ---------------------------------------------------------------------------
// Pure-gas dilute {title.lower()} {short}_i(T) [{unit}].
//
// Each entry maps a component to the coefficients [a0..a4] of the polynomial
//   {short}_i(T) = a0 + a1*T + a2*T^2 + a3*T^3 + a4*T^4,
// for T in [{T_MIN:.0f}, {T_MAX:.0f}] K (Fluid::kViscTmin/kViscTmax; held at the nearer limit
// outside). A polynomial costs about 1 ns per evaluation, so the solver can
// evaluate it in every cell at every residual. The sources (row comments):
//
{SOURCES}
//
// Each fitted row is within {MAX_FIT_ERROR:.0%} of its source where the source has data, and
// within {MAX_SHAPE_ERROR:.0%} of the kinetic-theory shape below a NASA interval. The value is
// temperature-only: the pressure dependence over PSA/TSA swings is small. A
// component without a row falls back to {floor}.
// ---------------------------------------------------------------------------

namespace adsim {{

  // C++17 inline variable: one shared instance across all translation units.
  inline const std::unordered_map<std::string, std::array<double, 5>> {var} = {{
{body}
  }};

}} // namespace adsim
"""


def main() -> int:
    ap = argparse.ArgumentParser(description=__doc__.split("\n\n")[0])
    ap.add_argument("--check", action="store_true", help="fail if a committed header differs")
    ap.add_argument("--trans", type=Path, help="local trans.inp (no download)")
    ap.add_argument("-v", "--verbose", action="store_true", help="print the fit error of each fitted row")
    args = ap.parse_args()
    visc, cond, report = build(trans_file(args.trans))
    if args.verbose:
        print("\n".join(report))
    stale = False
    for path, text in ((OUT_VISC, header("V", visc)), (OUT_COND, header("K", cond))):
        if args.check:
            old = path.read_text(encoding="utf-8")
            if old != text:
                stale = True
                sys.stdout.writelines(difflib.unified_diff(
                    old.splitlines(True), text.splitlines(True),
                    str(path.relative_to(ROOT)), "generated", n=1))
                print(f"\n{path.relative_to(ROOT)} is out of date: run python scripts/gen_transport_coefficients.py",
                      file=sys.stderr)
            else:
                print(f"{path.relative_to(ROOT)} is up to date")
        else:
            path.write_text(text, encoding="utf-8")
            print(f"wrote {path.relative_to(ROOT)}")
    return 1 if stale else 0


if __name__ == "__main__":
    sys.exit(main())
