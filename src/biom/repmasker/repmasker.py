import gzip
from functools import lru_cache
from pathlib import Path
from typing import Optional, Set, Tuple

_CLASSIFICATION = Path(__file__).parent.joinpath("hg38-repeat-masker-classification.tsv.gz")


@lru_cache(maxsize=1)
def _load_classification():
    mapping = {}
    with gzip.open(_CLASSIFICATION, 'rt') as stream:
        for line in stream:
            line = line.strip()
            if len(line) == 0:
                continue
            name, cls, family = line.split()
            assert name not in mapping, name
            mapping[name] = (name, family, cls)
    return mapping


# def family_to_cls(family: str) -> str:
#     mapping = {x[1]: x[2] for x in _load_classification(_REPEAT_MASKER_CLASSIFICATION).values()}
#     return mapping[family]


def classify(repname: str) -> Optional[Tuple[str, str, str]]:
    return _load_classification().get(repname, None)


@lru_cache(maxsize=1)
def names() -> Set[str]:
    return set(x[0] for x in _load_classification().values())


@lru_cache(maxsize=1)
def families() -> Set[str]:
    return set(x[1] for x in _load_classification().values())


@lru_cache(maxsize=1)
def classes() -> Set[str]:
    return set(x[2] for x in _load_classification().values())
