import pytest

from biom.primitives import Strand


def test_strand_str_repr():
    assert str(Strand.fwd) == "+"
    assert repr(Strand.fwd) == "+"
    assert str(Strand.rev) == "-"
    assert repr(Strand.rev) == "-"
    assert str(Strand.unknown) == "."
    assert repr(Strand.unknown) == "."


def test_strand_normalize():
    assert Strand.normalize("+") == Strand.fwd
    assert Strand.normalize(1) == Strand.fwd
    assert Strand.normalize(Strand.fwd) == Strand.fwd
    assert Strand.normalize("-") == Strand.rev
    assert Strand.normalize(-1) == Strand.rev
    assert Strand.normalize(Strand.rev) == Strand.rev
    assert Strand.normalize(".") == Strand.unknown
    assert Strand.normalize(0) == Strand.unknown
    assert Strand.normalize(Strand.unknown) == Strand.unknown


def test_strand_normalize_invalid():
    with pytest.raises(ValueError):
        Strand.normalize("invalid")
