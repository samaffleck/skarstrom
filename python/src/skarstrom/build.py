"""Building a simulation from scratch.

Everything a deck holds can be created in code, so a script can state every
choice and every value explicitly::

    import skarstrom as sk
    from skarstrom import ureg as u

    sim = sk.Simulation.new(title="CO2/N2 VSA")
    sim.add_component("CarbonDioxide")
    sim.add_component("Nitrogen")

    ads = sim.add_adsorbent("13X", solid_density=1138 * u.kg / u.m**3, heat_capacity=900 * u.J / (u.kg * u.K),
                            thermal_conductivity=0.63 * u.W / (u.m * u.K), particle_void_fraction=0.35)
    ads.set_isotherm("extended_dual_site_langmuir_c", {
        "CO2": {"m1": 2.768, "b1": 3.934e-7, "q1": 40.98 * u.kJ / u.mol, "m2": 2.758, "b2": 1.256e-7, "q2": 34.19 * u.kJ / u.mol},
        "N2":  {"m1": 2.768, "b1": 3.134e-6, "q1": 15.71 * u.kJ / u.mol, "m2": 2.758, "b2": 3.134e-6, "q2": 15.72 * u.kJ / u.mol},
    })
    ads.set_rate_model("CO2", "ldf", "macropore", tortuosity=3, macropore_void_fraction=0.35)
    ads.heat_of_adsorption["CO2"] = 30.85 * u.kJ / u.mol

    sim.bed.set_shape("cylinder", diameter=3.2 * u.cm)
    layer = sim.add_layer("13X", composition={"13X": 1.0}, length=0.3 * u.m, cells=30,
                          bed_void_fraction=0.41, particle_diameter=1.5 * u.mm)

    feed = sim.add_step("Feed", time=99 * u.s)
    feed.bottom.set_streams("inlet", [sk.Stream("Feed inlet", flow=sk.velocity(8.5 * u.cm / u.s),
                                                y={"CO2": 0.15, "N2": 0.85}, temperature=298 * u.K)])
    feed.top.set_pressure("outlet", 1 * u.bar)
    ...
    sim.kpi.track("CO2", feed=[feed.bottom], product=[evacuation.bottom])

Structure is created in the engine (core/include/adsim/objects/deck_builder.h), each model
chosen by a kind name (``skarstrom.build_kinds()`` lists them); values are
then set through the parameter registry, with the usual unit handling: a bare
number is SI, a pint Quantity carries its unit (a string is an expression or a
deck-syntax quantity such as "2 [min]").

:func:`to_script` writes this code for any deck, e.g. one made in a front end.
"""

from __future__ import annotations

import keyword
import math
import re
from dataclasses import dataclass, field
from typing import TYPE_CHECKING, Any, Dict, Iterable, List, Optional, Union

from . import _core

if TYPE_CHECKING:
    from . import Simulation
    from .model import Endpoint

__all__ = [
    "Flow",
    "KPIConfig",
    "PolynomialPower",
    "Profile",
    "Pump",
    "PumpStream",
    "Stream",
    "ValveStream",
    "build_kinds",
    "constant",
    "exponential",
    "linear",
    "mass_flow",
    "molar_flow",
    "polynomial_power",
    "to_script",
    "velocity",
    "volumetric_flow",
]


def build_kinds() -> Dict[str, List[str]]:
    """The model kinds the builder accepts, by family (isotherm, rate_model, dispersion...)."""
    return dict(_core.build_kinds())


# ---------------------------------------------------------------------------
# Value specifications
# ---------------------------------------------------------------------------


@dataclass(frozen=True)
class Profile:
    """How a boundary value changes during a step: constant, or a linear or
    exponential ramp from ``start`` to ``value`` over ``time``.

    For a pressure the ramp starts from the bed pressure and lasts the step,
    so only its end is given; temperature and flow ramps need all three.
    """

    kind: str
    value: Any = None
    start: Any = None
    time: Any = None


def constant(value: Any) -> Profile:
    """A value held for the whole step."""
    return Profile("constant", value)


def linear(end: Any, start: Any = None, time: Any = None) -> Profile:
    """A linear ramp to ``end`` (from ``start`` over ``time`` for temperatures and flows)."""
    return Profile("linear", end, start, time)


def exponential(end: Any, start: Any = None, time: Any = None) -> Profile:
    """An exponential approach to ``end`` (from ``start`` over ``time`` for temperatures and flows)."""
    return Profile("exponential", end, start, time)


def _profile(x: Any) -> Profile:
    return x if isinstance(x, Profile) else Profile("constant", x)


@dataclass(frozen=True)
class Flow:
    """A set flow: ``kind`` is velocity, volumetric_flow, molar_flow or mass_flow."""

    kind: str
    profile: Profile


def velocity(value: Any) -> Flow:
    """A superficial velocity (a value, or a :class:`Profile`)."""
    return Flow("velocity", _profile(value))


def volumetric_flow(value: Any) -> Flow:
    """A volumetric flow at the inlet conditions."""
    return Flow("volumetric_flow", _profile(value))


def molar_flow(value: Any) -> Flow:
    """A molar flow."""
    return Flow("molar_flow", _profile(value))


def mass_flow(value: Any) -> Flow:
    """A mass flow."""
    return Flow("mass_flow", _profile(value))


class Pump:
    """A compressor or vacuum pump on a stream or pressure boundary.

    ``kind`` is ``"compressor"`` or ``"vacuum_pump"``; ``energy`` the
    compression model (``"isothermal"`` or ``"adiabatic"``);
    ``efficiency_model`` ``"constant"`` (efficiency), ``"vacuum_pump"``
    (eta_max, alpha) or ``"polynomial"`` (poly_a ... poly_e). Other keyword
    values: polytropic_index, suction_pressure (compressor),
    discharge_pressure (vacuum pump).
    """

    def __init__(self, kind: str = "compressor", energy: str = "isothermal", efficiency_model: str = "constant",
                 **values: Any) -> None:
        self.kind = kind
        self.energy = energy
        self.efficiency_model = efficiency_model
        self.values = values

    def __repr__(self) -> str:
        extra = "".join(f", {k}={v!r}" for k, v in self.values.items())
        return f"Pump({self.kind!r}, energy={self.energy!r}, efficiency_model={self.efficiency_model!r}{extra})"


@dataclass
class Stream:
    """A stream with a set flow (:func:`velocity`, :func:`molar_flow`...).

    ``y`` and ``temperature`` are the gas that enters (for an outlet stream,
    when the flow reverses); ``temperature`` may be a :class:`Profile`.
    """

    name: str
    flow: Flow
    y: Optional[Dict[str, Any]] = None
    temperature: Any = None
    pump: Optional[Pump] = None


@dataclass(frozen=True)
class PolynomialPower:
    """Pump shaft power as a polynomial of the boundary pressure P [Pa]:
    ``a + b P + c P**2 + d P**3 + e P**4 + f P**5`` [W]."""

    a: Any = 0.0
    b: Any = 0.0
    c: Any = 0.0
    d: Any = 0.0
    e: Any = 0.0
    f: Any = 0.0


def polynomial_power(a: Any = 0.0, b: Any = 0.0, c: Any = 0.0, d: Any = 0.0, e: Any = 0.0,
                     f: Any = 0.0) -> PolynomialPower:
    """Pump power that depends on the boundary pressure (see :class:`PolynomialPower`)."""
    return PolynomialPower(a, b, c, d, e, f)


@dataclass
class PumpStream:
    """A stream driven by a pump of given shaft power: a compressor on an inlet,
    a vacuum pump on an outlet. The flow follows from the power, the pump's
    efficiency and the pressures.

    ``power`` is a constant power (W, or a pint quantity) or a
    :func:`polynomial_power` of the boundary pressure. ``pump`` gives the
    compression model and efficiency (a :class:`Pump`; its ``kind`` follows
    the direction and is ignored). ``y`` and ``temperature`` are the gas a
    compressor delivers.
    """

    name: str
    power: Any
    pump: Optional[Pump] = None
    y: Optional[Dict[str, Any]] = None
    temperature: Any = None


@dataclass
class ValveStream:
    """A stream through a valve against a pressure (a value or a :class:`Profile`)."""

    name: str
    pressure: Any
    cv: Any = None
    xt: Any = None
    y: Optional[Dict[str, Any]] = None
    temperature: Any = None


# ---------------------------------------------------------------------------
# Applying specifications (used by the model objects)
# ---------------------------------------------------------------------------


def _set_endpoint(sim: "Simulation", eid: int, leaf: str, value: Any) -> None:
    if value is not None:
        sim.set(f"endpoint[{eid}].{leaf}", value)


def _apply_profile(sim: "Simulation", eid: int, head: str, prof: Profile) -> None:
    if prof.kind == "constant":
        _set_endpoint(sim, eid, head, prof.value)
        return
    _set_endpoint(sim, eid, f"{head}_{prof.kind}_end", prof.value)
    if head == "pressure":
        if prof.start is not None or prof.time is not None:
            raise ValueError("a pressure ramp starts from the bed pressure and lasts the whole step; "
                             "give only its end value")
        return
    _set_endpoint(sim, eid, f"{head}_{prof.kind}_start", prof.start)
    _set_endpoint(sim, eid, f"{head}_{prof.kind}_time", prof.time)


def _apply_inlet(sim: "Simulation", eid: int, temperature: Any, y: Optional[Dict[str, Any]]) -> None:
    if temperature is not None:
        prof = _profile(temperature)
        sim._build("set_temperature_profile", eid, prof.kind)
        _apply_profile(sim, eid, "temperature", prof)
    for component, value in (y or {}).items():
        _set_endpoint(sim, eid, f"y[{sim._component(component)}]", value)


def _apply_pump(sim: "Simulation", eid: int, pump: Optional[Pump]) -> None:
    if pump is None:
        sim._build("set_pump", eid, "none", "isothermal", "constant")
        return
    sim._build("set_pump", eid, pump.kind, pump.energy, pump.efficiency_model)
    for name, value in pump.values.items():
        _set_endpoint(sim, eid, f"pump.{name}", value)


def _set_pressure(sim: "Simulation", step_id: int, location: str, direction: str, pressure: Any, temperature: Any,
                  y: Optional[Dict[str, Any]], pump: Optional[Pump]) -> "Endpoint":
    prof = _profile(pressure)
    eid = sim._build("set_pressure_boundary", step_id, location, direction, prof.kind)
    _apply_profile(sim, eid, "pressure", prof)
    _apply_inlet(sim, eid, temperature, y)
    if pump is not None:
        _apply_pump(sim, eid, pump)
    return sim._endpoint(eid)


def _set_streams(sim: "Simulation", step_id: int, location: str, direction: str,
                 streams: Iterable[Union[Stream, ValveStream, PumpStream]]) -> List["Endpoint"]:
    sim._build("set_flow_boundary", step_id, location, direction)
    ids = []
    for st in streams:
        if isinstance(st, ValveStream):
            prof = _profile(st.pressure)
            eid = sim._build("add_valve_stream", step_id, location, st.name, prof.kind)
            _apply_profile(sim, eid, "pressure", prof)
            _set_endpoint(sim, eid, "valve.cv", st.cv)
            _set_endpoint(sim, eid, "valve.xt", st.xt)
        elif isinstance(st, Stream):
            eid = sim._build("add_stream", step_id, location, st.name, st.flow.kind, st.flow.profile.kind)
            _apply_profile(sim, eid, st.flow.kind, st.flow.profile)
            if st.pump is not None:
                _apply_pump(sim, eid, st.pump)
        elif isinstance(st, PumpStream):
            pump = st.pump or Pump()
            poly = isinstance(st.power, PolynomialPower)
            eid = sim._build("add_pump_stream", step_id, location, st.name, "polynomial" if poly else "constant",
                             pump.energy, pump.efficiency_model)
            if poly:
                for coefficient in "abcdef":
                    _set_endpoint(sim, eid, f"power_poly_{coefficient}", getattr(st.power, coefficient))
            else:
                _set_endpoint(sim, eid, "power", st.power)
            for name, value in pump.values.items():
                _set_endpoint(sim, eid, f"pump.{name}", value)
        else:
            raise TypeError("streams are skarstrom.Stream, skarstrom.ValveStream or skarstrom.PumpStream, "
                            f"not {type(st).__name__}")
        _apply_inlet(sim, eid, st.temperature, st.y)
        ids.append(eid)
    return [sim._endpoint(e) for e in ids]


def _endpoint_ids(x: Any) -> List[int]:
    """Endpoint IDs from an Endpoint, a StepEnd (all its endpoints), an int, or a list of those."""
    from .model import Endpoint, StepEnd

    if isinstance(x, StepEnd):
        return [e.id for e in x]
    if isinstance(x, Endpoint):
        return [x.id]
    if isinstance(x, int):
        return [x]
    out: List[int] = []
    for item in x:
        out += _endpoint_ids(item)
    return out


class KPIConfig:
    """Which components get KPIs, and which streams are their feed, product and recycle."""

    def __init__(self, sim: "Simulation") -> None:
        self._sim = sim

    def track(self, component: str, feed: Any, product: Any, recycle: Any = ()) -> None:
        """Compute ``component``'s purity, recovery... with ``feed``, ``product`` and ``recycle``
        given as endpoints, step ends (every stream of that end) or endpoint IDs.

        ``recycle`` lists the streams that carry product back into the bed (a purge or a
        pressurisation fed from the product). Purity is the product's; recovery, productivity,
        production rate and specific energy count the product minus the recycle.
        """
        self._sim._build("set_kpi", self._sim._component(component), _endpoint_ids(feed), _endpoint_ids(product),
                         _endpoint_ids(recycle))

    def __repr__(self) -> str:
        rows = self._sim.describe()["kpis"]
        return "<kpi: " + (", ".join(r["component"] for r in rows) or "none tracked") + ">"


# ---------------------------------------------------------------------------
# Writing a deck as a script
# ---------------------------------------------------------------------------


def _literal(sim: "Simulation", path: str) -> str:
    """The value at ``path`` as Python source: a number, a bool, or the expression text."""
    p = sim._index().by_path.get(path)
    value = sim.get(path)
    if p is not None and p.choices and isinstance(value, int) and not isinstance(value, bool) \
            and 0 <= value < len(p.choices):
        return repr(p.choices[value])  # a model or mode switch, by name
    if isinstance(value, bool) or isinstance(value, int):
        return repr(value)
    if isinstance(value, str):
        return repr(value)
    if p is not None and p.kind == "expression":
        text = sim.get_text(path)
        try:
            float(text)
        except ValueError:
            return repr(text)  # a formula or a quantity: keep what the deck says
    if isinstance(value, float) and not math.isfinite(value):
        return f"float({str(value)!r})"
    return repr(value)


def _kwargs(sim: "Simulation", node: Any, skip: Iterable[str] = ()) -> str:
    skip = set(skip)
    return ", ".join(f"{f}={_literal(sim, node.path(f))}" for f in node.fields if f not in skip)


def _var(name: str, taken: set) -> str:
    v = re.sub(r"\W+", "_", name.strip().lower()).strip("_") or "step"
    if v[0].isdigit() or keyword.iskeyword(v):
        v = "s_" + v
    base, n = v, 2
    while v in taken or v in ("sim", "sk", "ads", "layer"):
        v = f"{base}_{n}"
        n += 1
    taken.add(v)
    return v


def to_script(sim: "Simulation", name: str = "sim") -> str:
    """Python source that builds ``sim`` from scratch with the builder API, every value explicit.

    Running the script gives a deck with the same structure and the same
    values (results and experimental data are not included).
    """
    d = sim.describe()
    out: List[str] = []
    w = out.append
    w("import skarstrom as sk")
    w("")
    w(f"{name} = sk.Simulation.new()")
    for f in ("title", "description", "author"):
        w(f"{name}.info.{f} = {_literal(sim, 'info.' + f)}")
    if d["variables"]:
        w("")
        w("# ---- user variables (other values' expressions may use them) ------------------")
    for var in d["variables"]:
        extra = f", description={var['description']!r}" if var["description"] else ""
        if var["formula"]:
            w(f"{name}.add_variable({var['name']!r}, category={var['category']!r}, formula={var['formula']!r}{extra})")
        else:
            value = _literal(sim, f"variable[{var['name']}]")
            w(f"{name}.add_variable({var['name']!r}, {value}, category={var['category']!r}{extra})")
    w("")
    w("# ---- fluid --------------------------------------------------------------------")
    for c in d["components"]:
        w(f"{name}.add_component({c!r})")
    fluid = sim.fluid
    for f in fluid.fields:
        w(f"{name}.fluid.{f} = {_literal(sim, fluid.path(f))}")
    for c, path in fluid.adsorbed_heat_capacity._paths.items():
        w(f"{name}.fluid.adsorbed_heat_capacity[{c!r}] = {_literal(sim, path)}")

    w("")
    w("# ---- adsorbents -----------------------------------------------------------------")
    ads_var: Dict[int, str] = {}
    taken: set = set()
    for a in d["adsorbents"]:
        node = sim.adsorbents.by_id(a["id"])
        v = _var("ads_" + a["name"], taken)
        ads_var[a["id"]] = v
        w(f"{v} = {name}.add_adsorbent({a['name']!r}, {_kwargs(sim, node)})")
        for iso in a["isotherms"]:
            params = []
            for c in iso["components"]:
                if c not in node.isotherm:  # a model without parameters (inert)
                    continue
                inode = node.isotherm[c]
                vals = ", ".join(f"{f!r}: {_literal(sim, inode.path(f))}" for f in inode.fields)
                params.append(f"    {c!r}: {{{vals}}},")
            if params:
                w(f"{v}.set_isotherm({iso['kind']!r}, {{")
                out.extend(params)
                w(f"}}, components={iso['components']!r})")
            else:
                w(f"{v}.set_isotherm({iso['kind']!r}, components={iso['components']!r})")
        for c, r in a["rate_models"].items():
            kn = node.kinetics[c]
            kw = _kwargs(sim, kn)
            w(f"{v}.set_rate_model({c!r}, {r['model']!r}, {r['coefficient']!r}{', ' + kw if kw else ''})")
        try:
            heats = node.heat_of_adsorption._paths
        except AttributeError:  # an adsorbent no component has been set up on
            heats = {}
        for c, path in heats.items():
            w(f"{v}.heat_of_adsorption[{c!r}] = {_literal(sim, path)}")

    w("")
    w("# ---- column ----------------------------------------------------------------------")
    bed, wall = sim.bed, sim.wall
    b = d["bed"]
    shape_fields = [f for f in ("diameter", "width", "length") if f in bed.fields]
    shape_kw = ", ".join(f"{f}={_literal(sim, bed.path(f))}" for f in shape_fields)
    w(f"{name}.bed.set_shape({b['shape']!r}{', ' + shape_kw if shape_kw else ''})")
    w(f"{name}.bed.ambient_temperature = {_literal(sim, bed.path('ambient_temperature'))}")
    w(f"{name}.bed.set_thermal_equilibrium({b['thermal_equilibrium']!r})")
    pkw = f", particle_htc={_literal(sim, bed.path('particle_htc'))}" if "particle_htc" in bed.fields else ""
    w(f"{name}.bed.set_particle_heat_transfer({b['particle_heat_transfer']!r}{pkw})")
    for void in ("bottom_void", "top_void"):
        vnode = getattr(bed, void)
        for f in vnode.fields:
            w(f"{name}.bed.{void}.{f} = {_literal(sim, vnode.path(f))}")
    wall_props = ("density", "heat_capacity", "thermal_conductivity")
    w(f"{name}.wall.set_material({b['wall_material']!r}, "
      + ", ".join(f"{f}={_literal(sim, wall.path(f))}" for f in wall_props) + ")")
    w(f"{name}.wall.thickness = {_literal(sim, wall.path('thickness'))}")
    fixed = f", fixed_temperature={_literal(sim, wall.path('fixed_temperature'))}" if "fixed_temperature" in wall.fields else ""
    w(f"{name}.wall.set_boundary({b['wall_boundary']!r}{fixed})")

    w("")
    w("# ---- layers (feed end first) ------------------------------------------------------")
    layer_var: Dict[int, str] = {}
    model_fields = {
        "dispersion": ("dispersion_coefficient", "molecular_diffusion_coefficient"),
        "heat_dispersion": ("heat_dispersion_coefficient",),
        "internal_heat_transfer": ("internal_htc", "leva_eta1", "leva_eta2", "leva_exp_coeff"),
        "external_heat_transfer": ("external_htc", "holman_c", "holman_m", "holman_n", "external_velocity"),
    }
    model_owned = {f for fs in model_fields.values() for f in fs} | {"momentum_model"}
    for l in d["layers"]:
        node = sim.layers.by_id(l["id"])
        v = _var("layer_" + l["name"], taken)
        layer_var[l["id"]] = v
        names = {i: sim.adsorbents.by_id(i).name for i, _ in l["composition"]}
        comp = "{" + ", ".join(f"{names[i]!r}: {fr!r}" for i, fr in l["composition"]) + "}"
        kw = _kwargs(sim, node, skip=model_owned)
        w(f"{v} = {name}.add_layer({l['name']!r}, composition={comp}, momentum={l['momentum']!r}, {kw})")
        for model, fs in model_fields.items():
            mkw = ", ".join(f"{f}={_literal(sim, node.path(f))}" for f in fs if f in node.fields)
            w(f"{v}.set_{model}({l[model]!r}{', ' + mkw if mkw else ''})")

    w("")
    w("# ---- initial conditions --------------------------------------------------------")
    ic = sim.initial_conditions
    for f in ic.fields:
        w(f"{name}.initial_conditions.{f} = {_literal(sim, ic.path(f))}")
    for c, path in ic.y._paths.items():
        w(f"{name}.initial_conditions.y[{c!r}] = {_literal(sim, path)}")

    w("")
    w("# ---- cycle -------------------------------------------------------------------------")
    cyc = sim.cycle
    for f in cyc.fields:
        w(f"{name}.cycle.{f} = {_literal(sim, cyc.path(f))}")
    ep_expr: Dict[int, str] = {}
    step_taken: set = set(taken)
    for s in d["steps"]:
        node = sim.steps.by_id(s["id"])
        v = _var(s["name"], step_taken)
        timing = ", ".join(f"{f}={_literal(sim, node.path(f))}" for f in ("time", "fraction") if f in node.fields)
        w(f"{v} = {name}.add_step({s['name']!r}{', ' + timing if timing else ''})")
        for loc in ("bottom", "top"):
            end = s[loc]
            if end["kind"] == "closed":
                w(f"{v}.{loc}.set_closed()")
            elif end["kind"] == "pressure":
                eid = end["endpoint"]
                ep_expr[eid] = f"{v}.{loc}"
                w(f"{v}.{loc}.set_pressure({end['direction']!r}, {_profile_src(sim, eid, 'pressure', end['profile'])}"
                  f"{_inlet_src(sim, eid, end['temperature_profile'])}{_pump_src(sim, eid, end['pump'])})")
            else:
                items = []
                names = [st["name"] for st in end["streams"]]
                for i, st in enumerate(end["streams"]):
                    eid = st["endpoint"]
                    key = repr(st["name"]) if st["name"] and names.count(st["name"]) == 1 else str(i)
                    ep_expr[eid] = f"{v}.{loc}[{key}]"
                    if st["type"] == "valve":
                        items.append(f"    sk.ValveStream({st['name']!r}, "
                                     f"pressure={_profile_src(sim, eid, 'pressure', st['profile'])}, "
                                     f"cv={_literal(sim, f'endpoint[{eid}].valve.cv')}, "
                                     f"xt={_literal(sim, f'endpoint[{eid}].valve.xt')}"
                                     f"{_inlet_src(sim, eid, st['temperature_profile'])}),")
                    elif st["type"] == "stream":
                        flow = f"sk.{st['flow']}({_profile_src(sim, eid, st['flow'], st['profile'])})"
                        items.append(f"    sk.Stream({st['name']!r}, flow={flow}"
                                     f"{_inlet_src(sim, eid, st['temperature_profile'])}"
                                     f"{_pump_src(sim, eid, st['pump'])}),")
                    elif st["type"] == "pump":
                        if st["power"] == "polynomial":
                            coeffs = ", ".join(f"{c}={_literal(sim, f'endpoint[{eid}].power_poly_{c}')}" for c in "abcdef")
                            power = f"sk.polynomial_power({coeffs})"
                        else:
                            power = _literal(sim, f"endpoint[{eid}].power")
                        items.append(f"    sk.PumpStream({st['name']!r}, power={power}"
                                     f"{_inlet_src(sim, eid, st['temperature_profile'])}"
                                     f"{_pump_src(sim, eid, st['pump'])}),")
                    else:
                        items.append(f"    # stream {st['name']!r}: a type the builder does not create")
                w(f"{v}.{loc}.set_streams({end['direction']!r}, [")
                out.extend(items)
                w("])")
        for layer_id in s["energy_sources"]:
            base = f"step[{s['id']}].energy[{layer_id}]"
            w(f"{v}.add_energy_source({layer_var.get(layer_id, layer_id)}, power={_literal(sim, base + '.power')}, "
              f"enabled={_literal(sim, base + '.enabled')})")
    for a_id, b_id in d["connections"]:
        w(f"{name}.connect({ep_expr.get(a_id, a_id)}, {ep_expr.get(b_id, b_id)})")
    for s in d["steps"]:
        for loc in ("bottom", "top"):
            for e in [s[loc]] + s[loc].get("streams", []):
                if e.get("buffer_tank"):
                    eid = e["endpoint"]
                    w(f"{ep_expr.get(eid, eid)}.set_buffer_tank({ep_expr.get(e['buffer_tank'], e['buffer_tank'])}, "
                      f"percentage={_literal(sim, f'endpoint[{eid}].buffer_tank.percentage')})")
    for k in d["kpis"]:
        feed = "[" + ", ".join(ep_expr.get(e, str(e)) for e in k["feed"]) + "]"
        product = "[" + ", ".join(ep_expr.get(e, str(e)) for e in k["product"]) + "]"
        recycle = ", ".join(ep_expr.get(e, str(e)) for e in k.get("recycle", []))
        recycle = f", recycle=[{recycle}]" if recycle else ""
        w(f"{name}.kpi.track({k['component']!r}, feed={feed}, product={product}{recycle})")

    w("")
    w("# ---- solver, cyclic steady state, energy balance --------------------------------------")
    for group in ("settings", "css", "energy_balance"):
        node = getattr(sim, group)
        for f in node.fields:
            w(f"{name}.{group}.{f} = {_literal(sim, node.path(f))}")
    return "\n".join(out) + "\n"


def _profile_src(sim: "Simulation", eid: int, head: str, kind: str) -> str:
    p = f"endpoint[{eid}].{head}"
    if kind == "constant":
        return _literal(sim, p)
    end = _literal(sim, f"{p}_{kind}_end")
    if head == "pressure":
        return f"sk.{kind}({end})"
    return f"sk.{kind}({end}, start={_literal(sim, f'{p}_{kind}_start')}, time={_literal(sim, f'{p}_{kind}_time')})"


def _inlet_src(sim: "Simulation", eid: int, temperature_kind: str) -> str:
    node_paths = [p.path for p in sim.parameters() if p.path.startswith(f"endpoint[{eid}].y[")]
    t = _profile_src(sim, eid, "temperature", temperature_kind)
    y = "{" + ", ".join(f"{p[len(f'endpoint[{eid}].y['):-1]!r}: {_literal(sim, p)}" for p in node_paths) + "}"
    return f", y={y}, temperature={t}"


def _pump_src(sim: "Simulation", eid: int, pump: Optional[Dict[str, str]]) -> str:
    if not pump:
        return ""
    prefix = f"endpoint[{eid}].pump."
    values = [(p.path[len(prefix):], p.path) for p in sim.parameters() if p.path.startswith(prefix)]
    kw = "".join(f", {k}={_literal(sim, path)}" for k, path in values)
    return (f", pump=sk.Pump({pump['kind']!r}, energy={pump['energy']!r}, "
            f"efficiency_model={pump['efficiency']!r}{kw})")
