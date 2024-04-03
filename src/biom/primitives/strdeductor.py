from typing import Callable, Literal

from pysam import AlignedSegment

StrandDeductor = Callable[[AlignedSegment], Literal["+", "-"]]


def get(protocol: Literal["f/s", "s/f", "f", "s"]) -> StrandDeductor:
    match protocol:
        case "f/s":
            return _fs
        case "s/f":
            return _sf
        case "f":
            return _f
        case "s":
            return _s
        case _:
            raise ValueError(f"Unknown stranding protocol: {protocol}")


def _f(read: AlignedSegment) -> Literal["+", "-"]:
    # Flip the mapped strand
    match read.is_reverse:
        case True:
            return "+"
        case False:
            return "-"
        case _:
            assert False, "Unreachable"


def _s(read: AlignedSegment) -> Literal["+", "-"]:
    # Keep the mapped strand
    match read.is_reverse:
        case True:
            return "-"
        case False:
            return "+"
        case _:
            assert False, "Unreachable"


def _fs(read: AlignedSegment) -> Literal["+", "-"]:
    match (read.is_read1, read.is_reverse):
        # First mate -> flip its strand
        case (True, True):
            return "+"
        case (True, False):
            return "-"
        # Second mate -> keep its strand
        case (False, True):
            return "-"
        case (False, False):
            return "+"
        case _:
            assert False, "Unreachable"


def _sf(read: AlignedSegment) -> Literal["+", "-"]:
    match (read.is_read1, read.is_reverse):
        # First mate -> keep its strand
        case (True, True):
            return "-"
        case (True, False):
            return "+"
        # Second mate -> flip its strand
        case (False, True):
            return "+"
        case (False, False):
            return "-"
        case _:
            assert False, "Unreachable"
