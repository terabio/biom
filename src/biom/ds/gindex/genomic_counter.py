from typing import TypeVar, Generic, Iterable

from biom.primitives import OrientationLike, Range
from .genomic_index import GenomicIndex
from .overlap import Overlap
from .resolve import Resolution, Counts

_T = TypeVar("_T")


class GenomicCounter(Generic[_T]):
    index: GenomicIndex[_T]
    resolution: Resolution[list[Overlap[_T]], Counts[_T]]
    _counts: Counts[_T]

    def __init__(
            self,
            index: GenomicIndex[_T],
            resolution: Resolution[list[Overlap[_T]], Counts[_T]]
    ):
        self.index = index
        self.resolution = resolution
        self._counts = {}

    def count(self, contig: str, orientation: OrientationLike, blocks: Iterable[Range]):
        overlaps = [self.index.overlap(contig, orientation, rng=rng) for rng in blocks]
        for k, v in self.resolution(overlaps).items():
            self._counts[k] = self._counts.get(k, 0) + v

    def counts(self) -> Counts[_T]:
        return self._counts

    def reset(self):
        self._counts = {}
