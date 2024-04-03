import os
from io import TextIOBase

import cattrs.preconf.pyyaml

from ..experiment import Experiment
from ..library import Library
from ..project import Project
from ..sample import Sample
from ..seqrun import SeqRun

_YAML_CONVERTER = cattrs.preconf.pyyaml.make_converter()


def _unstructure_hook(exp: Experiment) -> dict:
    return {
        "ind": exp.ind,
        "sample": exp.sample.ind,
        "library": _YAML_CONVERTER.unstructure(exp.library),
        "runs": _YAML_CONVERTER.unstructure(exp.runs)
    }


def _structure_hook(data: dict, ttype: type) -> Project:
    assert ttype is Project

    samples = _YAML_CONVERTER.structure(data["samples"], tuple[Sample, ...])
    if len(samples) != len(set(s.ind for s in samples)):
        raise ValueError("Sample IDs must be unique")

    samples_mapping = {s.ind: s for s in samples}
    experiments = tuple(Experiment(
        data["ind"],
        samples_mapping[data["sample"]],
        _YAML_CONVERTER.structure(data["library"], Library),
        _YAML_CONVERTER.structure(data["runs"], tuple[SeqRun, ...])
    ) for data in data["experiments"])

    return Project(data["ind"], tuple(experiments), samples)


_YAML_CONVERTER.register_unstructure_hook(Experiment, _unstructure_hook)
_YAML_CONVERTER.register_structure_hook(Project, _structure_hook)


def load(file: os.PathLike[str] | TextIOBase) -> Project:
    if isinstance(file, TextIOBase):
        return _YAML_CONVERTER.loads(file.read(), Project)
    else:
        with open(file) as f:
            return _YAML_CONVERTER.loads(f.read(), Project)


def loads(string: str) -> Project:
    return _YAML_CONVERTER.loads(string, Project)


def dump(project: Project, saveto: os.PathLike[str] | TextIOBase) -> str:
    string = dumps(project)
    if isinstance(saveto, TextIOBase):
        saveto.write(string)
    else:
        with open(saveto, "w") as f:
            f.write(string)
    return string


def dumps(project: Project) -> str:
    return _YAML_CONVERTER.dumps(project, sort_keys=False)
