from dataclasses import dataclass
from pathlib import Path
from typing import Callable, Any

from intervaltree import IntervalTree
from pybedtools import BedTool, Interval as BedInterval
from sortedcontainers import SortedList

from ..range import Range


def bedname(it: BedInterval) -> Any:
    return it.name


def from_bed(bed: Path, datafn: Callable[[BedInterval], Any] = bedname) -> 'Index':
    itrees = {}
    for it in BedTool(bed):
        key = (it.chrom, it.strand)
        tree = itrees.get(key, None)
        if tree is None:
            tree = IntervalTree()
            tree.addi(it.start, it.end, data=datafn(it))
            itrees[key] = tree
        else:
            tree.addi(it.start, it.end, data=datafn(it))
    return Index(itrees)


def merge(*indices: 'Index') -> 'Index':
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


@dataclass(frozen=True, slots=True)
class AnnotationSteps:
    boundaries: SortedList[int]
    annotation: list[set[Any]]

    def __len__(self) -> int:
        return len(self.annotation)


@dataclass(frozen=True, slots=True)
class AnnotationIntervals:
    intervals: list[Range]
    annotation: list[Any]

    def to_steps(self, rng: Range) -> AnnotationSteps:
        boundaries = SortedList()
        for it in self.intervals:
            boundaries.add(it.start)
            boundaries.add(it.end)

        if len(boundaries) == 0 or boundaries[0] != rng.start:
            boundaries.add(rng.start)

        if boundaries[-1] != rng.end:
            boundaries.add(rng.end)

        annotation = [set() for _ in range(len(boundaries) - 1)]

        for it, anno in zip(self.intervals, self.annotation):
            st, en = boundaries.bisect_left(it.start), boundaries.bisect_left(it.end)
            for stanno in annotation[st:en]:
                stanno.add(anno)

        return AnnotationSteps(boundaries, annotation)

    def __len__(self) -> int:
        return len(self.intervals)


class Index:
    itrees: dict[tuple[str, str], IntervalTree]

    def __init__(self, itrees=None):
        if itrees is None:
            itrees = {}
        self.itrees = itrees

    def overlap(self, contig: str, strand: str, start: int, end: int) -> AnnotationIntervals:
        index = self.itrees.get((contig, strand), None)
        if index is None:
            return AnnotationIntervals([], [])

        hits, annotation = [], []
        for it in index.overlap(start, end):
            hits.append(Range(
                max(start, it.begin), min(end, it.end)
            ))
            annotation.append(it.data)
        return AnnotationIntervals(hits, annotation)
