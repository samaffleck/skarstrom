"""Number formatting shared by the docs generators (JavaScript-compatible output)."""

from __future__ import annotations

import math


def js_str(x: float) -> str:
    """JavaScript Number.prototype.toString for the magnitudes used here."""
    if x == 0:
        return "0"
    s = repr(float(x))
    if "e" in s:
        m, e = s.split("e")
        if m.endswith(".0"):
            m = m[:-2]
        n = int(e)
        return f"{m}e{'+' if n >= 0 else '-'}{abs(n)}"
    if s.endswith(".0"):
        s = s[:-2]
    return s


def js_exp(x: float, digits: int) -> str:
    """JavaScript Number.prototype.toExponential(digits)."""
    m, e = f"{x:.{digits}e}".split("e")
    n = int(e)
    return f"{m}e{'+' if n >= 0 else '-'}{abs(n)}"


def js_prec(x: float, digits: int) -> str:
    """String(Number(x.toPrecision(digits)))."""
    return js_str(float(f"{x:.{digits}g}"))


def fmt(v, digits: int = 3) -> str:
    """`digits` significant figures; exponent form for very small or large values."""
    if v is None or (isinstance(v, float) and math.isnan(v)):
        return "—"
    a = abs(v)
    if a != 0 and (a < 1e-3 or a >= 1e5):
        return js_exp(v, 2)
    return js_prec(v, digits)


def is_fraction_metric(m: dict) -> bool:
    """Metrics whose error (and tolerance) is a fraction, shown as a percentage."""
    return m.get("metric") in ("rel_error", "nrmse", "max_rel_error")
