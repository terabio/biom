from typing import Protocol, TypeVar, Iterable

from attrs import define, field

from core import Interval
from ..resolve import Counts
from ..source import Source

_T = TypeVar('_T')
_K = TypeVar('_K')


@define(slots=True, frozen=True)
class CountingStats:
    time: float
    partition: Interval
    inside: float
    outside: float

    extra: dict[str, float] = field(factory=dict)

    def to_dict(self) -> dict:
        return {
            "Time(s)": self.time,
            "Partition": f"{self.partition.contig}:{self.partition.rng.start}-{self.partition.rng.end}",
            "Inside ROIs": self.inside,
            "Outside ROIs": self.outside,
            **self.extra,
        }


class MultiReadsCounter(Protocol[_T, _K]):
    def count(self, data: Iterable[_T], intervals: Iterable[Interval]):
        ...

    def counts(self) -> dict[_K, Counts[_T]]:
        ...

    def sources(self) -> dict[_K, Source]:
        ...

    def stats(self) -> Iterable[CountingStats]:
        ...

    def reset(self):
        ...
