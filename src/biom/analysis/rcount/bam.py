import time
from itertools import chain
from typing import TypeVar

from attrs import define, field

from biom.ds.gindex import GenomicCounter, Counts
from biom.io import bam
from biom.primitives import Range
from biom.primitives.strdeductor import StrandDeductor
from .read_counter import IndexedReadCounter

_T = TypeVar("_T")


@define(frozen=True, slots=True)
class BAMCounter(IndexedReadCounter[_T]):
    """
    A BAM-based reads counter. The counter can be either used to count all reads in a given file or, if the BAM file
    is indexed, it can be directed to process reads from a specific region. The counts can be retrieved using the
    `counts` method, and the statistics about the counting process can be retrieved using the `stats` method.
    """

    _inner: bam.Reader = field(alias="inner")
    _strander: StrandDeductor = field(alias="strander")
    _counter: GenomicCounter[_T] = field(alias="counter")
    _stats: dict = field(factory=dict, init=False)

    def process_all(self) -> Counts[_T]:
        start = time.time()
        for segment in self._inner:
            contig = segment.reference_name
            assert contig is not None, f"Segment reference name is unknown: {segment}"
            strand = self._strander(segment)
            blocks = tuple(Range(start, end) for start, end in segment.get_blocks())
            self._counter.count(contig, strand, blocks)
        end = time.time()

        self._stats["Time"] = self._stats.get("Time", 0) + (end - start)
        return self._counter.counts()

    def process_within(self, contig: str, start: int, end: int) -> Counts[_T]:
        self._inner.fetch(contig, start, end)
        return self.process_all()

    def stats(self) -> dict:
        stats = {
            "Time": self._stats.get("Time", 0),
        }
        if self._inner.statistics:
            assert self._inner.discarded is not None and self._inner.consumed is not None
            stats.update({
                "Discarded reads": self._inner.discarded.total,
                "Processed reads": self._inner.consumed.total,
            })
        return stats

    def counts(self) -> Counts[_T]:
        return self._counter.counts()

    def reset(self):
        self._counter.reset()


class PEBAMCounter(IndexedReadCounter[_T]):
    """
    A paired-end BAM-based reads counter. The counter can be either used to count all reads in a given file or, if the BAM
    file is indexed, it can be directed to process reads from a specific region. The counts can be retrieved using the
    `counts` method, and the statistics about the counting process can be retrieved using the `stats` method.

    Attributes:
        inner: Preconfigured paired-end BAM reader.
        strander: Strand deduction function.
        counter: Preconfigured genomic counter.
    """

    _inner: bam.PEReadsBundler = field(alias="inner")
    _strander: StrandDeductor = field(alias="strander")
    _counter: GenomicCounter[_T] = field(alias="counter")
    _stats: dict = field(factory=dict, init=False)

    def process_all(self) -> Counts[_T]:
        start = time.time()
        for lmate, rmate in self._inner:
            contig = lmate.reference_name
            assert contig is not None, f"Segment reference name is unknown: {lmate}, {rmate}"

            strand = self._strander(lmate)
            assert strand == self._strander(rmate)

            blocks = tuple(
                Range(start, end) for start, end in set(chain(lmate.get_blocks(), rmate.get_blocks()))
            )
            self._counter.count(contig, strand, blocks)
        end = time.time()

        self._stats["Time"] = self._stats.get("Time", 0) + (end - start)
        return self._counter.counts()

    def process_within(self, contig: str, start: int, end: int) -> Counts[_T]:
        self._inner.fetch(contig, start, end)
        return self.process_all()

    def stats(self) -> dict:
        stats = {
            "Time": self._stats.get("Time", 0),
            "Unpaired": self._inner.unpaired
        }

        if self._inner.reader.statistics:
            assert self._inner.reader.discarded is not None and self._inner.reader.consumed is not None
            stats.update({
                "Discarded reads": self._inner.reader.discarded.total,
                "Processed reads": self._inner.reader.consumed.total,
            })
        return stats

    def counts(self) -> Counts[_T]:
        return self._counter.counts()

    def reset(self):
        self._counter.reset()
