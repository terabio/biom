from itertools import chain
from typing import Self, Iterator, Protocol

from attrs import define

from biom.io import bam
from biom.primitives import Range, Strand
from biom.primitives.strdeductor import StrandDeductor


class Source(Protocol):
    """
    Abstraction that represents a source of reads. The source can be directed to fetch reads from a specific region
    and then iterated over to yield the reads in that region. Each read comes as a tuple of strand and a list of
    blocks (ranges) that denote the read's alignment to the reference genome.
    """

    def fetch(self, contig: str, start: int, end: int) -> Self:
        """
        Fetch reads from a specific region of the reference genome.
        :param contig: The contig to fetch reads from.
        :param start: The start position of the region.
        :param end: The end position of the region.
        :return: Self, for chaining.
        """
        pass

    def __iter__(self) -> Iterator[tuple[Strand, tuple[Range, ...]]]:
        """
        Iterate over the reads in the current region.
        :return: Iterator of reads, each read is a tuple of strand and a list of blocks (ranges). Each block denotes
        the read's alignment to the reference genome.
        """
        pass


@define(frozen=True, slots=True)
class SingleEndBAM(Source):
    obj: bam.Reader
    strander: StrandDeductor

    def fetch(self, contig: str, start: int, end: int) -> Self:
        self.obj.fetch(contig, start, end)
        return self

    def __iter__(self) -> Iterator[tuple[Strand, tuple[Range, ...]]]:
        for segment in self.obj:
            strand = self.strander(segment)

            blocks = tuple(Range(start, end) for start, end in segment.get_blocks())
            yield strand, blocks


@define(frozen=True, slots=True)
class PairedEndBAM(Source):
    obj: bam.PEReadsBundler
    strander: StrandDeductor

    def fetch(self, contig: str, start: int, end: int) -> Self:
        self.obj.fetch(contig, start, end)
        return self

    def __iter__(self) -> Iterator[tuple[Strand, tuple[Range, ...]]]:
        for lmate, rmate in self.obj:
            strand = self.strander(lmate)
            assert strand == self.strander(rmate)

            blocks = tuple(
                Range(start, end) for start, end in set(chain(lmate.get_blocks(), rmate.get_blocks()))
            )
            yield strand, blocks
