import pytest

from biom.primitives import Direction


def test_strand_str_repr():
    assert str(Direction.fwd) == "+"
    assert repr(Direction.fwd) == "+"
    assert str(Direction.rev) == "-"
    assert repr(Direction.rev) == "-"
    assert str(Direction.unknown) == "."
    assert repr(Direction.unknown) == "."


def test_strand_normalize():
    assert Direction.normalize("+") == Direction.fwd
    assert Direction.normalize(1) == Direction.fwd
    assert Direction.normalize(Direction.fwd) == Direction.fwd
    assert Direction.normalize("-") == Direction.rev
    assert Direction.normalize(-1) == Direction.rev
    assert Direction.normalize(Direction.rev) == Direction.rev
    assert Direction.normalize(".") == Direction.unknown
    assert Direction.normalize(0) == Direction.unknown
    assert Direction.normalize(Direction.unknown) == Direction.unknown


def test_strand_normalize_invalid():
    with pytest.raises(ValueError):
        Direction.normalize("invalid")
