from abc import ABC, abstractmethod
from typing import Callable, Any

from biom.gindex import Overlap

__all__ = ["Resolution", "Proportional", "Binary", "Priority", "Normalize", "AnnotationWeights"]

type AnnotationWeights[T] = dict[T, float]


class Resolution[I, O](ABC):
    # Resolution is required to turn a list of overlapping annotations into weights for one or more annotation categories.
    # This is necessary when a sequenced read/read pair isn't entirely within a known annotation 'type'.
    # All resolution strategies are chainable, so that multiple strategies can be applied in sequence.
    @abstractmethod
    def __call__(self, data: I) -> O:
        pass

    def apply[T](self, nxt: Callable[[O], T] | 'Resolution[O, T]') -> 'Resolution[I, T]':
        if isinstance(nxt, Resolution):
            return _Chain(self, nxt)
        elif isinstance(nxt, Callable):
            return _Chain(self, _Lambda(nxt))
        else:
            raise TypeError(f"Cannot apply {nxt} to {self}")

    def __or__[T](self, nxt: Callable[[O], T] | 'Resolution[O, T]') -> 'Resolution[I, T]':
        return self.apply(nxt)


class _Chain[I, O, T](Resolution[I, T]):
    def __init__(self, first: Resolution[I, O], second: Resolution[O, T]):
        self.first = first
        self.second = second

    def __call__(self, data: I) -> T:
        return self.second(self.first(data))


class _Lambda[I, O](Resolution[I, O]):
    def __init__(self, fn: Callable[[I], O]):
        self.fn = fn

    def __call__(self, data: I) -> O:
        return self.fn(data)


class Proportional[A](Resolution[list[Overlap[A]], AnnotationWeights[A]]):
    def __init__(self, empty: Any = None):
        self.empty = empty

    def __call__(self, overlaps: list[Overlap[A]]) -> AnnotationWeights[A]:
        """
        The proportional strategy assigns weights to annotation categories based on the overlap of the target ROI
        with each category
        """
        if len(overlaps) == 0:
            return {self.empty: 1}

        weights, covered, total = {}, 0, 0
        for roiov in overlaps:
            total += len(roiov.rng)  # original ROI length
            for rng, annotations in roiov.to_steps():  # step-wise ROI annotation
                length = len(rng)
                if not annotations:
                    covered += length
                    weights[self.empty] = weights.get(self.empty, 0) + length
                else:
                    covered += length
                    for anno in annotations:
                        weights[anno] = weights.get(anno, 0) + length

        assert covered == total, f"Covered length {covered} does not match total length {total}"
        return weights


class Binary[A](Resolution[list[Overlap[A]], AnnotationWeights[A]]):
    def __init__(self, empty: Any = None):
        self.empty = empty

    def __call__(self, overlaps: list[Overlap[A]]) -> AnnotationWeights[A]:
        """
        The binary strategy assigns weights to annotation categories based on the presence of at least a single base pair
        overlap between the read block and each category
        """
        categories = set()
        length, covered = 0, 0
        for blockov in overlaps:
            length += len(blockov.rng)
            for rng, step in blockov.to_steps():
                if step:
                    covered += len(rng)
                    categories |= step

        if not categories:
            return {self.empty: 1}

        if covered != length:
            categories.add(self.empty)

        weights = {cat: 1 for cat in categories}
        return weights


class Priority[A, D](Resolution[AnnotationWeights[A], AnnotationWeights[A]]):
    def __init__(self, priority: list[D], key: Callable[[D], A] = lambda x: x):
        self.priority = priority
        self.key = key

    def __call__(self, weights: AnnotationWeights[A]) -> AnnotationWeights[A]:
        """
        The priority strategy picks the first category in the priority list that is also present in the weights.
        If more than one category is mapped to the same priority, all such categories are returned.
        """
        keys = {}
        for cat in weights:
            keys.setdefault(self.key(cat), []).append(cat)

        categories = None
        for cat in self.priority:
            if cat in keys:
                categories = keys[cat]
                break

        if categories is None:
            raise ValueError(f"None of the priority categories {self.priority} are present in the weights {weights}")

        return {cat: weights[cat] for cat in categories}


class Normalize[A](Resolution[AnnotationWeights[A], AnnotationWeights[A]]):
    def __init__(self, total: float = 1.0):
        self.total = total

    def __call__(self, weights: AnnotationWeights[A]) -> AnnotationWeights[A]:
        """
        Normalize annotation weights to sum to a given total
        """
        if len(weights) == 0:
            return weights

        total = sum(weights.values())
        if total == 0:
            raise ValueError(f"Total weight is 0, cannot normalize {weights} to {self.total}")

        factor = self.total / total
        return {k: v * factor for k, v in weights.items()}
