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
