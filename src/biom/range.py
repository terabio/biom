from dataclasses import dataclass
from typing import Optional


@dataclass(slots=True)
class Range:
    start: int
    end: int

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
        return isinstance(other, Range) and self.start == other.start and self.end == other.end

    def contains(self, other: 'Range') -> bool:
        return self.start <= other.start < other.end <= self.end

    def overlap(self, other: 'Range') -> Optional[int]:
        overlap = min(self.end, other.end) - max(self.start, other.start)
        return overlap if overlap > 0 else None
