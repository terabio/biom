import pytest

from biom.primitives import Orientation


def test_strand_str_repr():
    assert str(Orientation.fwd) == "+"
    assert repr(Orientation.fwd) == "+"
    assert str(Orientation.rev) == "-"
    assert repr(Orientation.rev) == "-"
    assert str(Orientation.unstranded) == "."
    assert repr(Orientation.unstranded) == "."


def test_strand_normalize():
    assert Orientation.normalize("+") == Orientation.fwd
    assert Orientation.normalize(1) == Orientation.fwd
    assert Orientation.normalize(Orientation.fwd) == Orientation.fwd
    assert Orientation.normalize("-") == Orientation.rev
    assert Orientation.normalize(-1) == Orientation.rev
    assert Orientation.normalize(Orientation.rev) == Orientation.rev
    assert Orientation.normalize(".") == Orientation.unstranded
    assert Orientation.normalize(0) == Orientation.unstranded
    assert Orientation.normalize(Orientation.unstranded) == Orientation.unstranded


def test_strand_normalize_invalid():
    with pytest.raises(ValueError):
        Orientation.normalize("invalid")
