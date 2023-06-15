import math
from collections import defaultdict
from typing import Iterable, Any, Union, Literal, Callable

from .index import Index, AnnotationIntervals
from ..range import Range

AnnotateFn = Callable[['Annotator', dict[Any, float]], Any]
AnnotateStrategy = Union[
    Literal["nms", "frac-overlap"],
    tuple[Literal["priority"], tuple[Any, ...]],
    AnnotateFn
]

DisambiguateFn = Callable[['Annotator', AnnotationIntervals, Range, dict[Any, float]], None]
DisambiguateStrategy = Union[
    Literal["proportional"],
        # tuple[Literal["new-category"], str],
    DisambiguateFn
]


class Annotator:
    index: Index
    empty: Any
    disambigfn: DisambiguateFn
    annotatefn: AnnotateFn

    def __init__(
            self, index: Index, empty: Any = "NA",
            disambiguation: DisambiguateStrategy = "proportional",
            annotation: AnnotateStrategy = "nms"
    ):
        self.index = index
        self.empty = empty

        match disambiguation:
            case "proportional":
                self.disambigfn = proportional
            # case ("new-category", category):
            #     self.disambigfn = overlap_to_new_category(category)
            case _:
                self.disambigfn = disambiguation

        match annotation:
            case "nms":
                self.annotatefn = nms
            case "frac-overlap":
                self.annotatefn = frac_overlap
            # case "coverage":
            #     self.collapsefn = coverage
            case ("priority", scoring):
                self.annotatefn = priority(scoring)
            case _:
                self.annotatefn = annotation

    def annotate(self, contig: str, strand: str, blocks: Iterable[Range]) -> Any:
        overlap = defaultdict(int)
        for bl in blocks:
            intervals = self.index.overlap(contig, strand, bl.start, bl.end)
            match len(intervals):
                case 0:
                    overlap[self.empty] += len(bl)
                case 1:
                    anno, length = intervals.annotation[0], len(intervals.intervals[0])
                    overlap[anno] += length

                    diff = len(bl) - length
                    if diff > 0:
                        overlap[self.empty] += diff
                case _:
                    self.disambigfn(self, intervals, bl, overlap)

        assert math.isclose(sum(overlap.values()), sum(len(bl) for bl in blocks), abs_tol=1e-6)
        return self.annotatefn(self, overlap)

    def annotatei(self, contig: str, strand: str, blocks: Iterable[tuple[int, int]]) -> Any:
        blocks = [Range(start, end) for start, end in blocks]
        return self.annotate(contig, strand, blocks)


def nms(_: Annotator, overlap: dict[Any, float]) -> Any:
    return max(overlap.items(), key=lambda x: x[1])[0]


# def coverage(_: Annotator, overlap: dict[Any, float]) -> Any:
#     return overlap


def frac_overlap(_: Annotator, overlap: dict[Any, float]) -> Any:
    total = sum(overlap.values())
    return {k: v / total for k, v in overlap.items()}


def priority(scoring: tuple[Any, ...]) -> AnnotateFn:
    scoring = {k: ind for ind, k in enumerate(scoring)}

    def job(_: Annotator, overlap: dict[Any, float]) -> Any:
        return min(overlap.items(), key=lambda x: scoring[x[0]])[0]

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

# def overlap_to_new_category(category: Any) -> DisambiguateFn:
#     def job(_: Annotator, __: AnnotationIntervals, ___: Range, ____: dict[Any, float]) -> Any:
#         return category
#
#     return job
