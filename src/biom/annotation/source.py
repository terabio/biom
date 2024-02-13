from itertools import chain
from typing import Self, Iterator, Protocol

from biom.common import Range, StrandLike
from biom.io import bam


class Source(Protocol):
    def fetch(self, contig: str, strand: StrandLike, start: int, end: int) -> Self:
        pass

    def __iter__(self) -> Iterator[list[Range]]:
        pass


def wrap(bampe: bam.PEReadsBundler | None = None) -> Source:
    if bampe is not None:
        return _PEBundledBAM(bampe)
    else:
        raise ValueError("No source provided")


class _PEBundledBAM(Source):
    def __int__(self, obj: bam.PEReadsBundler):
        self.obj = obj

    def fetch(self, contig: str, strand: StrandLike, start: int, end: int) -> Self:
        self.obj.fetch(contig, start, end)
        return self

    def __iter__(self) -> Iterator[list[Range]]:
        for lmate, rmate in self.obj:
            blocks = {Range(start, end) for start, end in chain(lmate.get_blocks(), rmate.get_blocks())}
            yield list(blocks)
