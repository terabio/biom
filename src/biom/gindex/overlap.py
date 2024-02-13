from dataclasses import dataclass
from typing import Iterable

from sortedcontainers import SortedList

from biom.common import Range


@dataclass(frozen=True, slots=True)
class OverlapSteps[T]:
    rng: Range
    boundaries: list[Range]
    annotations: list[set[T]]

    def __post_init__(self):
        if len(self.boundaries) != len(self.annotations):
            raise ValueError("Boundaries and annotations must have the same length")
        if len(self.boundaries) == 0:
            raise ValueError("Boundaries and annotations must have at least one element")
        if self.boundaries[0].start != self.rng.start or self.boundaries[-1].end != self.rng.end:
            raise ValueError("Boundaries must start at the beginning of the range and end at the end of the range")

    def __len__(self) -> int:
        return len(self.annotations)

    def __iter__(self) -> Iterable[tuple[Range, set[T]]]:
        return zip(self.boundaries, self.annotations)


@dataclass(frozen=True, slots=True)
class Overlap[T]:
    rng: Range
    intervals: list[Range]
    annotations: list[T]

    def to_steps(self) -> OverlapSteps[T]:
        boundaries = {self.rng.start, self.rng.end}
        for it in self.intervals:
            boundaries.add(it.start)
            boundaries.add(it.end)

        boundaries = SortedList(boundaries)

        annotation = [set() for _ in range(len(boundaries) - 1)]

        for it, anno in zip(self.intervals, self.annotations):
            st, en = boundaries.bisect_left(it.start), boundaries.bisect_left(it.end)
            for stanno in annotation[st:en]:
                stanno.add(anno)
        boundaries = [Range(boundaries[i], boundaries[i + 1]) for i in range(len(boundaries) - 1)]

        return OverlapSteps(self.rng, boundaries, annotation)

    def __iter__(self) -> Iterable[tuple[Range, T]]:
        return zip(self.intervals, self.annotations)

    def __len__(self) -> int:
        return len(self.intervals)
