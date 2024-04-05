from collections import defaultdict
from pathlib import Path
from typing import Callable, TypeVar, Generic

from intervaltree import IntervalTree
from pybedtools import BedTool, Interval as BedInterval

from biom.primitives import Range, Direction, DirectionLike
from .overlap import Overlap

_T = TypeVar("_T")


def bedname(it: BedInterval) -> str:
    return it.name


def merge(*indices: "GenomicIndex[_T]") -> "GenomicIndex[_T]":
    if not indices:
        raise ValueError("Can't merge empty list of genomic indexes")
    elif len(indices) == 1:
        return indices[0]

    keys: set[tuple[str, Direction]] = set()
    for ind in indices:
        for key in ind.itrees.keys():
            keys.add(key)

    itrees = {}
    for key in keys:
        merged = set()
        for ind in indices:
            tree = ind.itrees.get(key, None)
            if tree is not None:
                merged |= set(tree)

        itrees[key] = IntervalTree(merged)
    return GenomicIndex(itrees)


class GenomicIndex(Generic[_T]):
    itrees: dict[tuple[str, Direction], IntervalTree]

    def __init__(self, itrees: dict[tuple[str, Direction], IntervalTree] | None = None):
        self.itrees = itrees if itrees else {}

    def add(self, contig: str, direction: DirectionLike, index: IntervalTree):
        self.itrees[(contig, Direction.normalize(direction))] = index

    def overlap(
            self,
            contig: str,
            direction: DirectionLike,
            start: int | None = None,
            end: int | None = None,
            rng: Range | None = None,
    ) -> Overlap[_T]:
        dr = Direction.normalize(direction)

        if (start is None or end is None) and rng is None:
            raise ValueError("Either start and end or range must be provided")
        elif (start is not None or end is not None) and rng is not None:
            raise ValueError("Either start and end or range must be provided, not both")

        if rng is None:
            if start is None or end is None:
                raise ValueError("Both start and end must be provided")
            rng = Range(start, end)

        index = self.itrees.get((contig, dr), None)
        if index is None:
            return Overlap(rng, [], [])

        hits, annotation = [], []
        for it in sorted(index.overlap(rng.start, rng.end), key=lambda x: x.begin):
            ov = rng.overlap(start=it.begin, end=it.end)
            assert ov is not None
            hits.append(ov.rng)
            annotation.append(it.data)
        return Overlap(rng, hits, annotation)

    def merge(self, *indices: 'GenomicIndex[_T]', inplace: bool = False) -> 'GenomicIndex[_T]':
        merged = merge(self, *indices)
        if inplace:
            self.itrees = merged.itrees
            return self
        else:
            return merged

    @staticmethod
    def from_bed(
            bed: Path,
            datafn: Callable[[BedInterval], _T] = bedname,  # type: ignore
            filterfn: Callable[[BedInterval], bool] | None = None,
    ) -> "GenomicIndex[_T]":
        filterfn = filterfn if filterfn is not None else (lambda _: True)

        itrees: defaultdict[tuple[str, Direction], IntervalTree] = defaultdict(IntervalTree)
        for it in BedTool(bed):
            if not filterfn(it):
                continue
            key = (it.chrom, Direction.normalize(it.strand))
            itrees[key].addi(it.start, it.end, data=datafn(it))
        result: dict[tuple[str, Direction], IntervalTree] = dict(itrees)
        return GenomicIndex(result)
