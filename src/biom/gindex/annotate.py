import math
from collections import defaultdict
from typing import Iterable, Any, Union, Literal, Callable

from .index import Index, AnnotationIntervals
from ..range import Range

CollapseFn = Callable[['Annotator', dict[Any, float]], Any]
CollapseStrategy = Union[
    Literal["nms", "norm-coverage", "coverage"],
    tuple[Literal["priority"], tuple[Any, ...]],
    CollapseFn
]

OverlapsFn = Callable[['Annotator', AnnotationIntervals, Range, dict[Any, float]], None]
OverlapsStrategy = Union[
    Literal["proportional"],
    tuple[Literal["new-category"], str],
    OverlapsFn
]


class Annotator:
    index: Index
    empty: Any
    overlapfn: OverlapsFn
    collapsefn: CollapseFn

    def __init__(self, index: Index, empty: Any = "NA",
                 overlap: OverlapsStrategy = "proportional", collapse: CollapseStrategy = "nms"):
        self.index = index
        self.empty = empty

        match overlap:
            case "proportional":
                self.overlapfn = proportional
            case ("new-category", category):
                self.overlapfn = overlap_to_new_category(category)
            case _:
                self.overlapfn = overlap

        match collapse:
            case "nms":
                self.collapsefn = nms
            case "norm-coverage":
                self.collapsefn = norm_coverage
            case "coverage":
                self.collapsefn = coverage
            case ("priority", scoring):
                self.collapsefn = priority(scoring)
            case _:
                self.collapsefn = collapse

    def annotate(self, contig: str, strand: str, blocks: Iterable[Range]) -> Any:
        overlap = defaultdict(int)
        for bl in blocks:
            intervals = self.index.overlap(contig, strand, bl.start, bl.end)
            match len(intervals):
                case 0:
                    overlap[self.empty] += len(bl)
                case 1:
                    anno, rng = intervals.annotation[0], intervals.intervals[0]
                    length = len(rng)
                    overlap[anno] += length
                    overlap[self.empty] += len(bl) - length
                case _:
                    self.overlapfn(self, intervals, bl, overlap)

        assert math.isclose(sum(overlap.values()), sum(len(bl) for bl in blocks), abs_tol=1e-6)
        return self.collapsefn(self, overlap)


def nms(_: Annotator, overlap: dict[Any, float]) -> Any:
    return max(overlap.items(), key=lambda x: x[0])[1]


def coverage(_: Annotator, overlap: dict[Any, float]) -> Any:
    return overlap


def norm_coverage(_: Annotator, overlap: dict[Any, float]) -> Any:
    total = sum(overlap.values())
    return {k: v / total for k, v in overlap.items()}


def priority(scoring: tuple[Any, ...]) -> CollapseFn:
    scoring = {k: ind for ind, k in enumerate(scoring)}

    def job(_: Annotator, overlap: dict[Any, float]) -> Any:
        return min(overlap.items(), key=lambda x: scoring[x[0]])

    return job


def proportional(self: Annotator, intervals: AnnotationIntervals, rng: Range, overlap: dict[Any, float]) -> None:
    stepwise = intervals.to_steps(rng)
    for ind, anno in enumerate(stepwise.annotation):
        length = stepwise.boundaries[ind + 1] - stepwise.boundaries[ind]
        if len(anno) == 0:
            overlap[self.empty] += length
        else:
            weight = length / len(anno)
            for a in anno:
                overlap[a] += weight


def overlap_to_new_category(category: Any) -> OverlapsFn:
    def job(_: Annotator, __: AnnotationIntervals, ___: Range, ____: dict[Any, float]) -> Any:
        return category

    return job
