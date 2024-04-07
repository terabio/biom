from attrs import define

from biom.primitives import Range, Orientation


@define(slots=True, frozen=True, eq=True, order=True, repr=True, hash=True)
class Interval:
    """
    Coordinates of an interval underlying a region of interest / reference region.

    :param contig: Parental contig of the region piece.
    :param strand: Strand of the region piece.
    :param start: Start position of the region piece.
    :param end: End position of the region piece.
    """
    contig: str
    strand: Orientation
    start: int
    end: int

    def range(self) -> Range:
        return Range(self.start, self.end)


@define(slots=True, frozen=True, eq=True, order=True, repr=True, hash=True)
class Region:
    """
    A region of interest/reference region composed of multiple intervals. Note that the intervals are not required to be
    on the same contig or strand.

    :param ind: Identifier of the region. Unique within the group.
    :param intervals: A tuple of non-overlapping intervals defining the region.
    """

    ind: str
    intervals: tuple[Interval, ...]


@define(slots=True, frozen=True, eq=True, order=True, repr=True, hash=True)
class Group:
    """
    Group of regions of interest sharing the same reference region. Differential usage analysis is performed on
    regions of interest relative to the reference region within each group.

    :param ind: Identifier of the group, must be globally unique.
    :param reference: Reference region of the group.
    :param regions: Regions of interest within the group. Each region will be compared to the reference region.
    """

    ind: str
    reference: Region
    regions: tuple[Region, ...]
