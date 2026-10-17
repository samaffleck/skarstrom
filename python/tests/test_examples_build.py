"""Every example script builds its cases from scratch; each must equal its shipped deck.

An example module (examples/<family>.py) that defines ``build(case)`` and
``DECKS`` (case -> deck name in examples/decks/) is checked case by case: the
simulation it builds, without loading any deck, must have the deck's structure
and values (title and description aside). So the examples are complete,
explicit descriptions of the validated decks, and they cannot drift from them.
"""

import importlib
import sys
from pathlib import Path

import pytest

import skarstrom
from deck_compare import diff, normalise, parameter_diff

EXAMPLES = Path(__file__).resolve().parents[2] / "examples"
pytestmark = pytest.mark.skipif(not EXAMPLES.is_dir(), reason="examples/ not found (not a repository checkout)")


def _families():
    if not EXAMPLES.is_dir():
        return []
    sys.path.insert(0, str(EXAMPLES))
    out = []
    for path in sorted(EXAMPLES.glob("*.py")):
        if path.name.startswith("_"):
            continue
        module = importlib.import_module(path.stem)
        if hasattr(module, "build") and hasattr(module, "DECKS"):
            out += [(path.stem, case, deck) for case, deck in module.DECKS.items()]
    return out


@pytest.mark.parametrize("family,case,deck", _families(), ids=lambda x: str(x))
def test_example_builds_its_deck(family, case, deck):
    module = importlib.import_module(family)
    built = module.build(case)
    shipped = skarstrom.examples.get(deck).load()
    differences = diff(normalise(built.to_json(), ignore_info=True), normalise(shipped.to_json(), ignore_info=True))
    assert differences == [], f"{family}.build({case!r}) differs from {deck}: {differences[:8]}"
    params = [d for d in parameter_diff(built, shipped) if not d[0].startswith("info.")]
    assert params == [], params[:8]


@pytest.mark.parametrize("family,deck", [("perez_2019", "Perez2019_E"), ("ward_2022", "Ward2022_Q100ccm"),
                                         ("tsa_illustrative", "Illustrative_TSA_13X")], ids=lambda x: str(x))
def test_example_runs_like_its_deck(family, deck):
    """Built and shipped simulations give the same results (2 cycles).

    Hand-written units ("100 ccm") can round differently from the deck's stored
    number in the last bit, which changes IDA's adaptive step sequence, so the
    runs are compared physically (interpolated onto the deck's time grid, to
    0.1 % of each signal's range, which covers interpolating a steep front
    between two grids) rather than bit for bit as for skarstrom.to_script.
    """
    import numpy as np

    module = importlib.import_module(family)
    if not hasattr(module, "DECKS"):
        pytest.skip(f"{family} does not build its cases yet")
    case = next(c for c, d in module.DECKS.items() if d == deck)
    shipped = skarstrom.examples.get(deck).load()
    built = module.build(case)
    cycles = min(2, shipped.cycle.number_of_cycles)
    ra, rb = shipped.run(cycles=cycles), built.run(cycles=cycles)
    assert ra.ok and rb.ok and ra.cycles_run == rb.cycles_run, (ra.message, rb.message)
    for var, comps in (("pressure", [None]), ("temperature", [None]), ("mole_fraction", shipped.components)):
        for comp in comps:
            for loc in ("top", "bottom"):
                a, b = shipped.series(var, comp, location=loc), built.series(var, comp, location=loc)
                scale = max(float(np.ptp(a.values)), 1e-4 * float(np.max(np.abs(a.values))), 1e-30)
                np.testing.assert_allclose(np.interp(a.time, b.time, b.values), a.values, rtol=0,
                                           atol=1e-3 * scale, err_msg=f"{var} {comp} {loc}")
    for ka, kb in zip(shipped.kpis(), built.kpis()):
        for key, value in ka.as_dict().items():
            if isinstance(value, float):
                np.testing.assert_allclose(kb.as_dict()[key], value, rtol=1e-4, err_msg=key)
