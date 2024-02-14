from abc import ABC, abstractmethod
from typing import Callable, Any

from biom.gindex import Overlap

__all__ = ["Resolution", "Proportional", "Binary", "Priority", "AnnotationWeights"]


class Resolution[I, O](ABC):
    # Resolution is required to turn a list of overlapping annotations into weights for one or more annotation categories.
    # This is necessary when a sequenced read/read pair isn't entirely within a known annotation 'type'.
    # All resolution strategies are chainable, so that multiple strategies can be applied in sequence.
    @abstractmethod
    def __call__(self, data: I) -> O:
        pass

    def chain[T](self, other: 'Resolution[O, T]') -> 'Resolution[I, T]':
        return ResolutionChain(self, other)

    def apply[T](self, nxt: Callable[[O], T] | 'Resolution[O, T]') -> 'Resolution[I, T]':
        if isinstance(nxt, Resolution):
            return ResolutionChain(self, nxt)
        elif isinstance(nxt, Callable):
            return ResolutionChain(self, ResolutionLambda(nxt))
        else:
            raise TypeError(f"Cannot apply {nxt} to {self}")


class ResolutionChain[I, O, T](Resolution[I, T]):
    def __init__(self, first: Resolution[I, O], second: Resolution[O, T]):
        self.first = first
        self.second = second

    def __call__(self, data: I) -> T:
        return self.second(self.first(data))


class ResolutionLambda[I, O](Resolution[I, O]):
    def __init__(self, fn: Callable[[I], O]):
        self.fn = fn

    def __call__(self, data: I) -> O:
        return self.fn(data)


type AnnotationWeights[T] = dict[T, float]


class Proportional[A](Resolution[list[Overlap[A]], AnnotationWeights[A]]):
    def __init__(self, empty: Any = None):
        self.empty = empty

    def __call__(self, overlaps: list[Overlap[A]]) -> AnnotationWeights[A]:
        """
        The proportional strategy assigns weights to annotation categories based on the proportion of the target ROI
        that overlaps with each category.
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
        weights = {anno: weight / covered for anno, weight in weights.items()}
        return weights


class Binary[A](Resolution[list[Overlap[A]], AnnotationWeights[A]]):
    def __init__(self, empty: Any = None, normalize: bool = False):
        self.empty = empty
        self.normalize = normalize

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

        weight = 1 / len(categories) if self.normalize else 1
        weights = {cat: weight for cat in categories}
        return weights


class Priority[A](Resolution[AnnotationWeights[A], AnnotationWeights[A]]):
    def __init__(self, priority: list[A], use_original_weights: bool = False):
        self.priority = priority
        self.use_original_weights = use_original_weights

    def __call__(self, weights: AnnotationWeights) -> AnnotationWeights:
        """
        The priority strategy picks the first category in the priority list that is also present in the weights. Then
        it assigns a weight of 1.0 to this category and discards all other categories. If use_original_weights is True,
        the original weights will be used for the category instead of 1.0.
        """
        category = None
        for cat in self.priority:
            if cat in weights:
                category = cat
                break
        if category is None:
            raise ValueError(f"None of the priority categories {self.priority} are present in the weights {weights}")

        if self.use_original_weights:
            return {category: weights[category]}
        else:
            return {category: 1.0}
