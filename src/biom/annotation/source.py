from itertools import chain
from typing import Self, Iterator, Protocol

from biom.common import Range, Strand
from biom.common.strdeductor import StrandDeductor
from biom.io import bam


class Source(Protocol):
    def fetch(self, contig: str, start: int, end: int) -> Self:
        pass

    def __iter__(self) -> Iterator[tuple[Strand, list[Range]]]:
        pass


def wrap(strander: StrandDeductor, bampe: bam.PEReadsBundler | None = None) -> Source:
    if bampe is not None:
        return _PEBundledBAM(strander, bampe)
    else:
        raise ValueError("No source provided")


class _PEBundledBAM(Source):
    def __init__(self, strander: StrandDeductor, obj: bam.PEReadsBundler):
        self.strander = strander
        self.obj = obj

    def fetch(self, contig: str, start: int, end: int) -> Self:
        self.obj.fetch(contig, start, end)
        return self

    def __iter__(self) -> Iterator[tuple[Strand, list[Range]]]:
        for lmate, rmate in self.obj:
            strand = self.strander(lmate)
            assert strand == self.strander(rmate)

            blocks = {Range(start, end) for start, end in chain(lmate.get_blocks(), rmate.get_blocks())}
            yield strand, list(blocks)
