from typing import Protocol, TypeVar, Iterable

from attrs import define, field

from biom.primitives import Interval
from ..resolve import Counts
from ..source import Source

_T = TypeVar('_T')


@define(slots=True, frozen=True)
class CountingStats:
    time: float
    partition: Interval
    has_overlap: float
    no_overlap: float

    extra: dict[str, float] = field(factory=dict)

    def to_dict(self) -> dict:
        return {
            "Time(s)": self.time,
            "Partition": f"{self.partition.contig}:{self.partition.rng.start}-{self.partition.rng.end}",
            "Has overlap": self.has_overlap,
            "No overlap": self.no_overlap,
            **self.extra,
        }


class MultiReadsCounter(Protocol[_T]):
    def count(self, data: Iterable[_T], intervals: Iterable[Interval]):
        ...

    def counts(self) -> dict[str, Counts[_T]]:
        ...

    def sources(self) -> dict[str, Source]:
        ...

    def stats(self) -> Iterable[CountingStats]:
        ...

    def reset(self):
        ...
