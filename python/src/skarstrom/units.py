"""Units: how skarstrom reads and writes dimensional values.

A bare number is in the parameter's SI unit: seconds, kelvin, pascal, metre,
mol/s and so on (the unit is in every docstring, in the repr of every object
and in :attr:`skarstrom.Parameter.unit`). To use another unit, say so::

    sim.steps["Feed"].time = "2 min"              # a string quantity
    sim.steps["Feed"].time = "2 [min]"            # the same, in deck syntax
    sim.steps["Feed"].time = 2 * ureg.minute      # a pint Quantity
    sim.set("step[Feed].time", 2, unit="min")     # an explicit unit argument
    sim.steps["Feed"].get("time", unit="min")     # read back in minutes

Quantities are checked against the parameter's dimension, so ``"2 bar"`` on a
step time raises ``ValueError`` instead of setting 200 000 s. String units are
parsed by the engine (the same parser decks and the CLI use): SI symbols with
prefixes (``kPa``, ``mmol``, ``cm``), compound units (``kJ/mol``,
``W/(m^2*K)``, ``mmol/g``), and common named units (``bar``, ``atm``,
``degC``, ``min``, ``h``, ``L/min``, ``SCCM``, ``Nm3/h``, ``%``). Everything
after a ``/`` is in the denominator, so ``W/m2-K`` is W/(m²·K).

pint Quantities work without any setup: skarstrom converts them with pint
itself (``to_base_units``), so any unit pint knows is accepted.

To make units mandatory, create the simulation with
``Simulation(deck, require_units=True)``: a bare number for a dimensional
parameter then raises ``TypeError``.
"""

from __future__ import annotations

from dataclasses import dataclass
from typing import Any, Tuple

from . import _core

__all__ = ["Unit", "convert", "dimension_name", "is_quantity", "parse"]

# Exponents of (length, mass, time, temperature, amount of substance).
Dimension = Tuple[int, int, int, int, int]

_PINT_DIMENSIONS = {
    "[length]": 0,
    "[mass]": 1,
    "[time]": 2,
    "[temperature]": 3,
    "[substance]": 4,
}


@dataclass(frozen=True)
class Unit:
    """A parsed unit: ``SI value = value * factor + offset``."""

    text: str
    factor: float
    offset: float
    dimension: Dimension
    dimension_name: str  # "time", "pressure", ... or a base-unit formula

    @property
    def dimensionless(self) -> bool:
        return self.dimension == (0, 0, 0, 0, 0)


def parse(unit: str) -> Unit:
    """Parse a unit expression (``"kJ/mol"``, ``"W/(m^2*K)"``...); ``ValueError`` if unknown."""
    d = _core.parse_unit(unit)
    return Unit(unit, d["factor"], d["offset"], tuple(d["dimension"]), d["dimension_name"])


def convert(value: float, from_unit: str, to_unit: str) -> float:
    """Convert ``value`` between two units of the same dimension.

    ``convert(25, "degC", "K")`` is 298.15; ``convert(1, "bar", "s")`` raises
    ``ValueError``.
    """
    return _core.convert_unit(float(value), from_unit, to_unit)


def dimension_name(unit: str) -> str:
    """The physical dimension of a unit: ``dimension_name("kPa")`` is ``"pressure"``."""
    return parse(unit).dimension_name


def is_quantity(value: Any) -> bool:
    """True for a pint-style Quantity (anything with ``magnitude``, ``units`` and ``to_base_units``)."""
    return hasattr(value, "magnitude") and hasattr(value, "units") and hasattr(value, "to_base_units")


def _quantity_to(value: Any, unit: str, what: str) -> float:
    """A pint Quantity as a number in ``unit`` (dimension-checked)."""
    target = parse(unit)
    base = value.to_base_units()
    dims = [0, 0, 0, 0, 0]
    for name, power in dict(base.dimensionality).items():
        if name not in _PINT_DIMENSIONS:
            raise ValueError(f"{what} expects {_describe(target)}; {value!s} has dimension {name}")
        dims[_PINT_DIMENSIONS[name]] = int(power)
    if tuple(dims) != target.dimension:
        raise ValueError(f"{what} expects {_describe(target)}; {value!s} is {_name_of(tuple(dims))}")
    magnitude = base.magnitude
    try:
        magnitude = float(magnitude)
    except TypeError:
        raise TypeError(f"{what} takes a scalar quantity, got {value!r}") from None
    return (magnitude - target.offset) / target.factor


def _describe(unit: Unit) -> str:
    if unit.dimensionless:
        return "a dimensionless value"
    return f"a {unit.dimension_name} [{unit.text}]"


def _name_of(dimension: Dimension) -> str:
    """Dimension name for a dimension tuple, via the base-unit expression."""
    symbols = ("m", "kg", "s", "K", "mol")
    parts = [f"{s}^{p}" for s, p in zip(symbols, dimension) if p]
    if not parts:
        return "dimensionless"
    return "a " + parse("*".join(parts)).dimension_name
