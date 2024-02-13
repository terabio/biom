from typing import Callable

from pysam import AlignedSegment

StrandDeductor = Callable[[AlignedSegment, AlignedSegment], str]


def get(protocol: str):
    if protocol == "f/s":
        return _fs
    elif protocol == "s/f":
        return _sf
    else:
        raise ValueError(f"Unknown stranding protocol: {protocol}")


def _fs(lmate: AlignedSegment, rmate: AlignedSegment) -> str:
    assert lmate.is_read1 and rmate.is_read2
    # (-/+) -> +
    if lmate.is_reverse:
        assert not rmate.is_reverse
        return "+"
    # (+/-) -> -
    else:
        assert rmate.is_reverse
        return "-"


def _sf(lmate: AlignedSegment, rmate: AlignedSegment) -> str:
    assert lmate.is_read1 and rmate.is_read2
    # (-/+) -> -
    if lmate.is_reverse:
        assert not rmate.is_reverse
        return "-"
    # (+/-) -> +
    else:
        assert rmate.is_reverse
        return "+"
