import pytest

from biom.common import strdeductor


def test_get():
    assert strdeductor.get("f/s") == strdeductor._fs
    assert strdeductor.get("s/f") == strdeductor._sf
    assert strdeductor.get("f") == strdeductor._f
    assert strdeductor.get("s") == strdeductor._s
    with pytest.raises(ValueError):
        strdeductor.get("invalid")  # type: ignore


def test_f():
    assert strdeductor._f(_AlignedSegmentStub(True)) == "+"
    assert strdeductor._f(_AlignedSegmentStub(False)) == "-"


def test_s():
    assert strdeductor._s(_AlignedSegmentStub(True)) == "-"
    assert strdeductor._s(_AlignedSegmentStub(False)) == "+"


def test_fs():
    assert strdeductor._fs(_AlignedSegmentStub(True, True)) == "+"
    assert strdeductor._fs(_AlignedSegmentStub(True, False)) == "-"
    assert strdeductor._fs(_AlignedSegmentStub(False, True)) == "-"
    assert strdeductor._fs(_AlignedSegmentStub(False, False)) == "+"


def test_sf():
    assert strdeductor._sf(_AlignedSegmentStub(True, True)) == "-"
    assert strdeductor._sf(_AlignedSegmentStub(True, False)) == "+"
    assert strdeductor._sf(_AlignedSegmentStub(False, True)) == "+"
    assert strdeductor._sf(_AlignedSegmentStub(False, False)) == "-"


class _AlignedSegmentStub:
    def __init__(self, is_reverse, is_read1=None):
        self.is_reverse = is_reverse
        self.is_read1 = is_read1
