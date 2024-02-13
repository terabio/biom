from pathlib import Path
from typing import Callable, Self

from intervaltree import IntervalTree
from pybedtools import BedTool, Interval as BedInterval

from biom.common import Range, Strand, StrandLike
from .overlap import Overlap


class Index[T]:
    itrees: dict[tuple[str, Strand], IntervalTree]

    def __init__(self, itrees: dict[tuple[str, Strand], IntervalTree] | None = None):
        self.itrees = itrees if itrees else {}

    def overlap(
            self,
            contig: str,
            strand: StrandLike,
            start: int | None = None,
            end: int | None = None,
            rng: Range | None = None,
    ) -> Overlap[T]:
        strand = Strand.normalize(strand)

        if (start is None or end is None) and rng is None:
            raise ValueError("Either start and end or range must be provided")
        elif (start is not None or end is not None) and rng is not None:
            raise ValueError("Either start and end or range must be provided, not both")

        if rng is None:
            rng = Range(start, end)

        index = self.itrees.get((contig, strand), None)
        if index is None:
            return Overlap(rng, [], [])

        hits, annotation = [], []
        for it in sorted(index.overlap(rng.start, rng.end), key=lambda x: x.begin):
            ov = rng.overlap(start=it.begin, end=it.end)
            assert ov is not None
            hits.append(ov.rng)
            annotation.append(it.data)
        return Overlap(rng, hits, annotation)

    def merge(self, *indices: Self, inplace: bool = False) -> Self:
        merged = merge(self, *indices)
        if inplace:
            self.itrees = merged.itrees
            return self
        else:
            return merged


def merge[T](*indices: Index[T]) -> Index[T]:
    if not indices:
        raise ValueError("Can't merge empty list of genomic indexes")
    elif len(indices) == 1:
        return indices[0]

    keys = set()
    for ind in indices:
        keys |= ind.itrees.keys()

    itrees = {}
    for key in keys:
        merged = set()
        for ind in indices:
            tree = ind.itrees.get(key, None)
            if tree is not None:
                merged |= set(tree)

        itrees[key] = IntervalTree(merged)
    return Index(itrees)


def bedname(it: BedInterval) -> str:
    return it.name


def from_bed[T](
        bed: Path,
        datafn: Callable[[BedInterval], T] = bedname,
        filterfn: Callable[[BedInterval], bool] | None = None,
) -> Index[T]:
    filterfn = filterfn if filterfn is not None else (lambda _: True)

    itrees = {}
    for it in BedTool(bed):
        if not filterfn(it):
            continue
        key = (it.chrom, Strand.normalize(it.strand))
        tree = itrees.setdefault(key, IntervalTree())
        tree.addi(it.start, it.end, data=datafn(it))
    return Index(itrees)
