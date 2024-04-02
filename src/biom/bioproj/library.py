from enum import Enum

from attrs import define, field


class Stranding(Enum):
    Unstranded = "unstranded"
    Forward = "forward"
    Reverse = "reverse"
    Unknown = "unknown"

    def __eq__(self, other):
        if isinstance(other, Stranding):
            return self.value == other.value
        elif isinstance(other, str):
            return self.value == other
        else:
            return False


@define(slots=True, frozen=True, eq=True, order=True, repr=True, hash=True)
class Library:
    """
    A class to describe meta information about a sequencing library made from a biological sample.

    Attributes
    ----------
    source : tuple[str, ...]
        What molecules were used to generate the library?
    selection : tuple[str, ...]
        Were there any selection/enrichment steps during library generation?
    stranding : Stranding
        What is the stranding of the library?
    """
    source: tuple[str, ...] = field(converter=lambda x: tuple(x))
    selection: tuple[str, ...] = field(converter=lambda x: tuple(x))
    stranding: Stranding = field(converter=Stranding)

    @source.validator
    def check_source(self, _, value):
        if not value:
            raise ValueError("Library source must be specified")

    @selection.validator
    def check_selection(self, _, value):
        if not value:
            raise ValueError("Library selection method must be specified")
