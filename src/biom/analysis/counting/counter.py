from typing import TypeVar, Generic

from biom.ds.gindex import GenomicIndex, Overlap
from .resolve import Resolution, AnnotationWeights
from .source import Source

_T = TypeVar("_T")


class Counter(Generic[_T]):
    def __init__(
            self,
            index: GenomicIndex[_T],
            source: Source,
            resolution: Resolution[list[Overlap[_T]], AnnotationWeights[_T]]
    ):
        self.index: GenomicIndex[_T] = index
        self.source: Source = source
        self.resolution: Resolution[list[Overlap[_T]], AnnotationWeights[_T]] = resolution
        self.processed: int = 0
        self.results: AnnotationWeights[_T] = {}

    def annotate(self, contig: str, start: int, end: int) -> AnnotationWeights[_T]:
        for strand, blocks in self.source.fetch(contig, start, end):
            overlaps = [self.index.overlap(contig, strand, rng=rng) for rng in blocks]
            for k, v in self.resolution(overlaps).items():
                self.results[k] = self.results.get(k, 0) + v
            self.processed += 1
        return self.results
