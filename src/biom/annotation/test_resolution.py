import pytest

from biom.annotation.resolve import ResolutionChain, ResolutionLambda, Proportional, Binary, Priority
from biom.common import Range
from biom.gindex import Overlap


# Test ResolutionChain
def test_resolution_chain():
    # Define two simple resolutions
    res1 = ResolutionLambda(lambda x: x + 1)
    res2 = ResolutionLambda(lambda x: x * 2)

    # Chain them
    chain = ResolutionChain(res1, res2)

    # Test
    assert chain(2) == 6  # (2 + 1) * 2
    assert res1.chain(res2)(2) == 6


# Test ResolutionLambda
def test_resolution_lambda():
    # Define a simple resolution
    res = ResolutionLambda(lambda x: x + 1)

    # Test
    assert res(2) == 3
    assert res.apply(lambda x: x * 2)(2) == 6


# Test Proportional
def test_proportional():
    # Define a simple proportional resolution
    prop = Proportional(empty="None")

    # Test with no overlaps
    assert prop([]) == {"None": 1}
    assert prop([Overlap(Range(0, 10), [], [])]) == {"None": 1}

    # Test with a single overlap
    overlaps = [Overlap(Range(0, 10), [Range(0, 7)], ['a'])]
    assert prop(overlaps) == {"a": 0.7, "None": 0.3}

    overlaps = [Overlap(Range(0, 7), [Range(0, 7)], ['a'])]
    assert prop(overlaps) == {"a": 1}

    # Test with several overlaps
    overlaps = [
        Overlap(Range(0, 10), [Range(0, 7), Range(3, 9)], ['a', 'b']),
        Overlap(Range(11, 15), [Range(14, 15)], ['a']),
    ]
    weights = {"a": 8, "b": 6, "None": 4}
    weights = {k: v / 14 for k, v in weights.items()}

    assert prop(overlaps) == weights


def test_binary():
    def _test(normalize):
        binary = Binary(empty=None, normalize=normalize)

        # Test with no overlaps
        assert binary([]) == {None: 1}
        assert binary([Overlap(Range(0, 10), [], [])]) == {None: 1}

        # Test with a single overlap
        overlaps = [Overlap(Range(0, 10), [Range(0, 7)], ['a'])]
        assert binary(overlaps) == {"a": 1, None: 1} if not normalize else {"a": 0.5, None: 0.5}

        overlaps = [Overlap(Range(0, 7), [Range(0, 7)], ['a'])]
        assert binary(overlaps) == {"a": 1}

        # Test with several overlaps
        overlaps = [
            Overlap(Range(0, 10), [Range(0, 7), Range(3, 9)], ['a', 'b']),
            Overlap(Range(11, 15), [Range(14, 15)], ['a']),
        ]
        assert binary(overlaps) == {"a": 1, "b": 1, None: 1} if not normalize else {"a": 1 / 3, "b": 1 / 3, None: 1 / 3}

    _test(normalize=False)
    _test(normalize=True)


def test_priority():
    priority_original_weight = Priority(['a', 'b'], use_original_weights=True)
    priority_renormalized = Priority(['a', 'b'], use_original_weights=False)

    # Test unknown category
    with pytest.raises(ValueError):
        priority_original_weight({'c': 1})
        priority_renormalized({'c': 1})

    # Test with a single category
    assert priority_original_weight({'a': 1}) == {'a': 1}
    assert priority_renormalized({'a': 1}) == {'a': 1}

    # Test with several categories
    weights = {'a': 0.5, 'b': 0.5}
    assert priority_original_weight(weights) == {'a': 0.5}
    assert priority_renormalized(weights) == {'a': 1.0}
