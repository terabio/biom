from enum import Enum
from typing import Literal, Union

from .strand import Strand


class Direction(Enum):
    fwd = "+"
    rev = "-"
    unknown = "."

    def __repr__(self):
        return self.value

    def __str__(self):
        return self.value

    @classmethod
    def normalize(cls, direction: 'DirectionLike') -> 'Direction':
        match direction:
            case "+" | 1:
                return Direction.fwd
            case "-" | -1:
                return Direction.rev
            case "." | 0:
                return Direction.unknown
            case _ if isinstance(direction, cls):
                return direction
            case _ if isinstance(direction, Strand):
                return cls.from_strand(direction)
            case _:
                raise ValueError(f"Unknown strand: {direction}")

    @classmethod
    def from_strand(cls, strand: Strand) -> 'Direction':
        match strand:
            case Strand.fwd:
                return cls.fwd
            case Strand.rev:
                return cls.rev
            case _:
                raise ValueError(f"Unknown strand: {strand}")


DirectionLike = Union[Strand, Direction, Literal["+", "-", ".", 1, -1, 0]]
