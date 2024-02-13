from enum import Enum
from typing import Literal


class Strand(Enum):
    fwd = "+"
    rev = "-"
    unknown = "."

    def __str__(self):
        return self.value

    def __repr__(self):
        return str(self)

    @classmethod
    def normalize(cls, strand: 'str | int | Strand') -> 'Strand':
        match strand:
            case "+" | 1:
                return Strand.fwd
            case "-" | -1:
                return Strand.rev
            case "." | 0:
                return Strand.unknown
            case _ if isinstance(strand, cls):
                return strand
            case _:
                raise ValueError(f"Unknown strand: {strand}")


StrandLike = Strand | Literal["+", "-", ".", 1, -1, 0]
