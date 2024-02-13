from dataclasses import dataclass


@dataclass(slots=True)
class OverlapInfo:
    size: int
    rng: "Range"


@dataclass(slots=True)
class Range:
    start: int
    end: int

    def __post_init__(self):
        if self.start > self.end:
            raise ValueError(
                f"Start must be less than or equal to end. Got {self.start} > {self.end}."
            )

    def __add__(self, other):
        match other:
            case int():
                return Range(self.start + other, self.end + other)
            case _:
                raise TypeError()

    def __sub__(self, other):
        return self.__add__(-other)

    def __len__(self):
        return self.end - self.start

    def __eq__(self, other):
        return (
            isinstance(other, Range)
            and self.start == other.start
            and self.end == other.end
        )

    def __contains__(self, item: 'Range') -> bool:
        return self.start <= item.start < item.end <= self.end

    def overlap(
        self,
        other: 'Range' | None = None,
        start: int | None = None,
        end: int | None = None,
    ) -> OverlapInfo | None:
        if other is None:
            if start is None or end is None:
                raise ValueError("Either other or start and end must be provided")
            other = Range(start, end)

        start = max(self.start, other.start)
        end = min(self.end, other.end)
        size = end - start
        if size <= 0:
            return None
        return OverlapInfo(size, Range(start, end))
