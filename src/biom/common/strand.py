from enum import Enum
from typing import Literal, TypeVar, Generic

from attrs import define, field


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

T = TypeVar("T")


@define(slots=True, frozen=True, eq=True, repr=True, hash=True)
class Stranded(Generic[T]):
    fwd: T = field()
    rev: T = field()

    def with_rev(self, value: T) -> 'Stranded[T]':
        return Stranded(self.fwd, value)

    def with_fwd(self, value: T) -> 'Stranded[T]':
        return Stranded(value, self.rev)

    def __iter__(self):
        yield from (self.fwd, self.rev)
