"""The Example type behind :mod:`skarstrom.examples`."""

from __future__ import annotations

import os
import shutil
from typing import TYPE_CHECKING, Union

if TYPE_CHECKING:
    from . import Simulation

PathLike = Union[str, "os.PathLike[str]"]


def _decks_dir() -> str:
    return os.path.join(os.path.dirname(__file__), "data", "decks")


class Example:
    """One example deck shipped with the package.

    An ``Example`` works wherever a deck path does (it is ``os.PathLike``), so
    ``skarstrom.Simulation(skarstrom.examples.Perez2019_A)`` loads it, as does
    ``skarstrom.examples.Perez2019_A.load()``.
    """

    __slots__ = ("name", "title", "description")

    def __init__(self, name: str, title: str = "", description: str = "") -> None:
        self.name = name
        #: The deck's title (its System_Information title).
        self.title = title or name
        #: The deck's description: source, case and notes.
        self.description = description

    @property
    def path(self) -> str:
        """Absolute path of the deck file."""
        return os.path.join(_decks_dir(), self.name + ".json")

    def __fspath__(self) -> str:
        return self.path

    def load(self) -> "Simulation":
        """A new :class:`skarstrom.Simulation` of this deck."""
        from . import Simulation

        return Simulation(self.path)

    def copy_to(self, destination: PathLike = ".") -> str:
        """Copy the deck to ``destination`` (a directory or a file path) and return the copy's path."""
        dst = os.fspath(destination)
        if os.path.isdir(dst):
            dst = os.path.join(dst, self.name + ".json")
        shutil.copyfile(self.path, dst)
        return dst

    def __eq__(self, other: object) -> bool:
        return isinstance(other, Example) and other.name == self.name

    def __hash__(self) -> int:
        return hash(self.name)

    def __repr__(self) -> str:
        return f"<skarstrom example {self.name}>"
