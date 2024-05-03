from collections.abc import Iterable
from pathlib import Path
from typing import TypeVar, Callable

from biom import io
from ...nfcore import rnaseq
from core import strdeductor, bioproj
from .bam import SingleEndBAM, PairedEndBAM
from .source import Source

_K = TypeVar('_K')


def default_se_bam_reader(path: Path) -> biom.M.io.bam.Reader:
    return biom.M.io.bam.Reader(path, inflags=0, exflags=2564, minmapq=0, statistics=True)


def default_pe_bam_reader(path: Path) -> biom.M.io.bam.PEReadsBundler:
    reader = biom.M.io.bam.Reader(path, inflags=3, exflags=2564, minmapq=0, statistics=True)
    return biom.M.io.bam.PEReadsBundler(reader)


def default_key(prj: rnaseq.Project, experiment: rnaseq.Experiment) -> tuple[str, str]:
    return prj.ind, experiment.ind


def from_rnaseq(
        projects: Iterable[rnaseq.Project], *,
        se_bam: Callable[[Path], biom.M.io.bam.Reader] = default_se_bam_reader,
        pe_bam: Callable[[Path], biom.M.io.bam.PEReadsBundler] = default_pe_bam_reader,
        key: Callable[[rnaseq.Project, rnaseq.Experiment], _K] = default_key,  # type: ignore
) -> dict[_K, Source]:
    results: dict[_K, Source] = {}
    for prj in projects:
        for experiment in prj.experiments:
            stranding = strdeductor.from_bioproj(experiment.source)

            runs = list(experiment.runs)
            layout = runs[0].layout
            for r in runs[1:]:
                if r.layout != layout:
                    raise ValueError(
                        f"All runs must have the same layout. Found different layouts: {layout} and {r.layout} in {experiment}"
                    )

            record = key(prj, experiment)
            if record in results:
                raise ValueError(f"Duplicate key: {record}")

            match layout:
                case bioproj.SeqLayout.Single:
                    results[record] = SingleEndBAM(se_bam(experiment.bam), stranding)
                case bioproj.SeqLayout.Paired:
                    results[record] = PairedEndBAM(pe_bam(experiment.bam), stranding)
                case _:
                    raise ValueError(f"Unsupported layout: {layout}")
    return results
