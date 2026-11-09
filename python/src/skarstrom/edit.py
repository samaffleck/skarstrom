"""Structural deck edits: the engine's cascades (``adsim/objects/deck_edit.h``).

``sim.edit`` changes a deck's structure the way a front end does: adding,
removing and moving layers, steps and streams, and switching a model to
another kind. Every edit keeps the rest of the deck consistent (connections,
KPI references, buffer tanks, mole fractions), and a kind switch keeps the
parameters the old and new kinds share (the same name and unit; see
:func:`skarstrom.model_kinds`). Values are then set through the parameter
registry as usual (``sim["layer[3].length"] = 0.2``).

Objects are given by ID or by their Python object (``sim.steps["Feed"]``,
``sim.layers[0]``, an endpoint); kinds by their tags (``"valve"``,
``"macropore"``, ...). A refused edit raises ``ValueError`` with the engine's
message and changes nothing.

>>> sim.edit.add_layer()                        # "Layer 2", seeded composition
>>> sim.edit.set_isotherm_kind("13X", "CarbonDioxide", "extended_dual_site_langmuir_p")
>>> sim.edit.set_stream_kind("Feed", "bottom", 0, "valve")
"""

from __future__ import annotations

from typing import TYPE_CHECKING, Any, Dict, List

from . import _core

if TYPE_CHECKING:  # pragma: no cover
    from . import Simulation

__all__ = ["Editor", "model_kinds"]


def model_kinds() -> Dict[str, List[Dict[str, Any]]]:
    """Every model kind's descriptor, by family.

    Each kind is ``{"tag", "label", "legacy", "parameters": [{"name", "unit",
    "default", "lower", "upper", "tunable", "label"}]}`` in its family's stored
    order. A legacy kind is not offered (the former rate coefficients, now
    presets). The rate-coefficient kind ``series_resistance`` also has
    ``resistances`` (tag, label, parent, scaled_by_capacity, parameters, LaTeX
    ``fragments`` per shape, source), ``capacity_modes`` and ``presets``.
    Families: rate_model, rate_coefficient, momentum, dispersion,
    heat_dispersion, internal_heat_transfer, external_heat_transfer,
    particle_heat_transfer, boundary, stream, profile, flow, pump,
    pump_energy, pump_efficiency, pump_power and isotherm (with ``competitive``,
    ``water_coupled``, ``fittable``, ``iast`` and ``paired``). A profile
    parameter whose unit is ``"@"`` takes the unit of the profiled quantity.
    """
    return dict(_core.model_kinds())


class Editor:
    """The structural cascades of one :class:`~skarstrom.Simulation` (``sim.edit``)."""

    __slots__ = ("_sim",)

    def __init__(self, sim: "Simulation") -> None:
        self._sim = sim

    # ---- plumbing ------------------------------------------------------------------
    def _do(self, op: str, *args: Any) -> Any:
        try:
            return self._sim._sim.edit(op, list(args))
        finally:
            self._sim._invalidate()

    def _step(self, step: Any) -> int:
        if isinstance(step, int):
            return step
        if hasattr(step, "id") and not isinstance(step, str):
            return int(step.id)
        return int(self._sim.steps[step].id)

    def _layer(self, layer: Any) -> int:
        return int(self._sim._layer_id(layer))

    def _adsorbent(self, adsorbent: Any) -> int:
        return int(self._sim._adsorbent_id(adsorbent))

    @staticmethod
    def _endpoint(endpoint: Any) -> int:
        if isinstance(endpoint, int):
            return endpoint
        eid = int(getattr(endpoint, "id", -1))
        if eid < 0:
            raise ValueError("give one endpoint (a stream or a pressure boundary), not a step end with several")
        return eid

    # ---- layers ----------------------------------------------------------------------
    def add_layer(self, name: str = "") -> int:
        """Append a layer ("Layer N" when unnamed; a taken name is made unique) with a
        seeded composition; returns its layer ID."""
        return int(self._do("add_layer", name))

    def remove_layer(self, layer: Any) -> None:
        """Remove a layer and the steps' energy sources on it."""
        self._do("remove_layer", self._layer(layer))

    def move_layer(self, layer: Any, offset: int) -> None:
        """Move a layer one place towards the product end (+1) or the feed end (-1)."""
        self._do("move_layer", self._layer(layer), int(offset))

    # ---- steps -------------------------------------------------------------------------
    def add_step(self, name: str = "") -> int:
        """Append a step with the default ends (a one-stream flow inlet at the bottom, a
        constant-pressure outlet at the top; compositions sum to 1); returns its step ID."""
        return int(self._do("add_step", name))

    def remove_step(self, step: Any) -> None:
        """Remove a step; its connections, KPI references and the buffer tanks drawing from it go too."""
        self._do("remove_step", self._step(step))

    def duplicate_step(self, step: Any) -> int:
        """Copy a step (inserted after it, disconnected); returns the copy's step ID."""
        return int(self._do("duplicate_step", self._step(step)))

    def move_step(self, step: Any, offset: int) -> None:
        """Move a step one place later (+1) or earlier (-1)."""
        self._do("move_step", self._step(step), int(offset))

    # ---- streams and connections ----------------------------------------------------------
    def add_stream(self, step: Any, side: str, kind: str = "flow") -> int:
        """Add a ``"flow"``, ``"valve"`` or ``"pump"`` stream to a step's flow end; returns its endpoint ID."""
        return int(self._do("add_stream", self._step(step), side, kind))

    def remove_stream(self, step: Any, side: str, index: int) -> None:
        """Remove stream ``index`` of a step's flow end, with its connection and references."""
        self._do("remove_stream", self._step(step), side, int(index))

    def connect(self, a: Any, b: Any) -> None:
        """Connect two endpoints under a fresh pair id."""
        self._do("connect", self._endpoint(a), self._endpoint(b))

    # ---- sources and sinks -----------------------------------------------------------------
    def _reservoir(self, reservoir: Any) -> int:
        if isinstance(reservoir, int):
            return reservoir
        if isinstance(reservoir, str):
            return int(self._sim.reservoirs[reservoir].id)
        return int(reservoir.id)

    def add_reservoir(self, kind: str, name: str = "") -> int:
        """Add a ``"source"`` or ``"sink"`` ("Source N" / "Sink N" when unnamed); returns its reservoir ID."""
        return int(self._do("add_reservoir", kind, name))

    def reservoir_from_endpoint(self, endpoint: Any, name: str = "") -> int:
        """Add a source (for an inlet) or sink (for an outlet) with the endpoint's temperature,
        mole fractions and pressure, and attach the endpoint to it; returns its reservoir ID."""
        return int(self._do("reservoir_from_endpoint", self._endpoint(endpoint), name))

    def remove_reservoir(self, reservoir: Any) -> None:
        """Remove a source or sink; the endpoints attached to it keep its values as their own."""
        self._do("remove_reservoir", self._reservoir(reservoir))

    def rename_reservoir(self, reservoir: Any, name: str) -> None:
        """Rename a source or sink."""
        self._do("rename_reservoir", self._reservoir(reservoir), name)

    def attach_reservoir(self, endpoint: Any, reservoir: Any) -> None:
        """Attach an endpoint to a source or sink; ``None`` or 0 detaches it."""
        self._do("attach_reservoir", self._endpoint(endpoint), 0 if reservoir in (None, 0) else self._reservoir(reservoir))

    def detach_mismatched_reservoirs(self) -> int:
        """Detach every endpoint whose direction no longer agrees with its reservoir's kind; returns the count."""
        return int(self._do("detach_mismatched_reservoirs"))

    # ---- kind switches ---------------------------------------------------------------------
    def set_isotherm_kind(self, adsorbent: Any, component: str, kind: str) -> None:
        """The isotherm of one component; shared parameters (same name and unit) are kept."""
        self._do("set_isotherm_kind", self._adsorbent(adsorbent), self._sim._component(component), kind)

    def set_mass_transfer_kind(self, adsorbent: Any, component: str, model: str, coefficient: str) -> None:
        """The rate model (``"ldf"``, ``"ldf2"``, ``"dual_rate_ldf"``) of one component and a coefficient
        preset (``"constant"``, ``"arrhenius"``, ``"macropore"``, ``"macropore_knudsen"``, ``"haynes"``,
        ``"macropore_capacity"``) or ``"series_resistance"``, which keeps the current switches."""
        self._do("set_mass_transfer_kind", self._adsorbent(adsorbent), self._sim._component(component), model,
                 coefficient)

    def set_mass_transfer_resistance(self, adsorbent: Any, component: str, tag: str, on: bool,
                                     site: str = "fast") -> None:
        """Switch one resistance (``"film"``, ``"macropore"``, ``"knudsen"``, ``"micropore"``,
        ``"lumped"``) of a component's coefficient on or off, on the ``"fast"`` or ``"slow"`` site.
        Refused (ValueError): the last of film, macropore, micropore and lumped off, knudsen without
        macropore."""
        self._do("set_mass_transfer_resistance", self._adsorbent(adsorbent), self._sim._component(component),
                 site, tag, bool(on))

    def set_capacity_factor_mode(self, adsorbent: Any, component: str, mode: str, site: str = "fast") -> None:
        """The capacity-factor mode (``"none"``, ``"slope"``, ``"secant"``, ``"fixed"``) of a
        component's coefficient, on the ``"fast"`` or ``"slow"`` site."""
        self._do("set_capacity_factor_mode", self._adsorbent(adsorbent), self._sim._component(component),
                 site, mode)

    def set_dispersion_kind(self, layer: Any, kind: str) -> None:
        self._do("set_dispersion_kind", self._layer(layer), kind)

    def set_heat_dispersion_kind(self, layer: Any, kind: str) -> None:
        self._do("set_heat_dispersion_kind", self._layer(layer), kind)

    def set_internal_heat_transfer_kind(self, layer: Any, kind: str) -> None:
        self._do("set_internal_heat_transfer_kind", self._layer(layer), kind)

    def set_external_heat_transfer_kind(self, layer: Any, kind: str) -> None:
        self._do("set_external_heat_transfer_kind", self._layer(layer), kind)

    def set_particle_heat_transfer_kind(self, kind: str) -> None:
        self._do("set_particle_heat_transfer_kind", kind)

    def set_boundary_kind(self, step: Any, side: str, kind: str) -> None:
        """Make a step end a ``"pressure"``, ``"flow"`` or ``"closed"`` boundary."""
        self._do("set_boundary_kind", self._step(step), side, kind)

    def set_stream_kind(self, step: Any, side: str, index: int, kind: str) -> None:
        """The device of a stream (``"flow"``, ``"valve"``, ``"pump"``); it keeps its endpoint ID."""
        self._do("set_stream_kind", self._step(step), side, int(index), kind)

    def set_flow_kind(self, endpoint: Any, kind: str) -> None:
        """The flow quantity of a set-flow stream (velocity, volumetric_flow, molar_flow, mass_flow)."""
        self._do("set_flow_kind", self._endpoint(endpoint), kind)

    def set_profile_kind(self, endpoint: Any, target: str, kind: str) -> None:
        """The profile of an endpoint's ``"pressure"``, ``"flow"`` or ``"temperature"``."""
        self._do("set_profile_kind", self._endpoint(endpoint), target, kind)

    def set_pump_kind(self, endpoint: Any, kind: str) -> None:
        """``"none"``, ``"compressor"`` or ``"vacuum_pump"``; efficiency and energy models carry over."""
        self._do("set_pump_kind", self._endpoint(endpoint), kind)

    def set_pump_efficiency_kind(self, endpoint: Any, kind: str) -> None:
        self._do("set_pump_efficiency_kind", self._endpoint(endpoint), kind)

    def set_pump_power_kind(self, endpoint: Any, kind: str) -> None:
        self._do("set_pump_power_kind", self._endpoint(endpoint), kind)

    # ---- identity and repair -----------------------------------------------------------------
    def ensure_ids(self) -> None:
        """Stable IDs by first appearance, connection pair ids renumbered 1..N."""
        self._do("ensure_ids")

    def repair_connections(self) -> int:
        """Re-link connections from their pair ids; returns how many orphans were disconnected."""
        return int(self._do("repair_connections"))

    def sync_kpi_config(self) -> None:
        """One KPI entry per component; references to missing endpoints dropped."""
        self._do("sync_kpi_config")

    def fill_gaps(self) -> None:
        """Fill every entry run-time initialisation would otherwise fill (idempotent)."""
        self._do("fill_gaps")
