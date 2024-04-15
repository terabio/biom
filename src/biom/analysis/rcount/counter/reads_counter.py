from typing import Protocol, TypeVar, Iterable

from biom.primitives import Interval
from ..resolve import Counts
from ..source import Source

_T = TypeVar('_T')


class MultiReadsCounter(Protocol[_T]):
    def count(self, data: Iterable[_T], intervals: Iterable[Interval]):
        ...

    def counts(self) -> dict[str, Counts[_T]]:
        ...

    def sources(self) -> dict[str, Source]:
        ...

    def stats(self) -> dict:
        ...

    def reset(self):
        ...
