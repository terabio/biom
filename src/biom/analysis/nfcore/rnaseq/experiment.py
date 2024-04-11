import os
from itertools import chain
from pathlib import Path
from typing import Iterable

from attrs import define, field

from biom.primitives import Stranded, bioproj


@define(slots=True, frozen=True, eq=True, order=True, repr=True, hash=True, init=False)
class Experiment:
    descriptor: str = field()
    source: bioproj.Experiment = field()
    # Files generated by the pipeline
    bam: Path = field(converter=lambda x: Path(x))
    bigwig: Stranded[Path] = field(converter=lambda x: Stranded(*[Path(p) for p in x]))
    salmon: Path = field(converter=lambda x: Path(x))

    def __init__(self, descriptor: str, source: bioproj.Experiment, resfolder: os.PathLike[str]):
        # Folder with the results of the pipeline
        resfolder = Path(resfolder)
        # if not resfolder.is_dir():
        #     raise ValueError(f"Results folder {resfolder} does not exist")

        # BAM file
        bam = resfolder / "star_salmon" / f"{descriptor}.markdup.sorted.bam"

        # Bigwig
        if source.library.stranding in {bioproj.Stranding.Unstranded, bioproj.Stranding.Unknown}:
            raise ValueError(
                f"Unstranded experiments or experiments with unknown strandedness are not supported: {descriptor}"
            )

        fwd = resfolder / "star_salmon" / "bigwig" / f"{descriptor}.forward.bigWig"
        rev = resfolder / "star_salmon" / "bigwig" / f"{descriptor}.reverse.bigWig"

        bigwig = Stranded(fwd, rev)

        # Salmon file
        salmon = resfolder / "star_salmon" / descriptor / "quant.sf"

        self.__attrs_init__(descriptor, source, bam, bigwig, salmon)

    @property
    def ind(self) -> str:
        return self.source.ind

    @property
    def library(self) -> bioproj.Library:
        return self.source.library

    @property
    def sample(self) -> bioproj.Sample:
        return self.source.sample

    @property
    def runs(self) -> Iterable[bioproj.SeqRun]:
        return self.source.runs

    @property
    def ensure_exists(self):
        allfiles = chain(
            [self.bam, self.bam.with_suffix(".bam.bai"), self.bigwig.fwd, self.bigwig.rev, self.salmon],
            *[r.files for r in self.runs]
        )
        missing = [f for f in allfiles if not f.exists()]
        if missing:
            report = "\n".join(f"\t{f}" for f in missing)
            raise FileNotFoundError(
                f"nf-core/rnaseq project ({self.source.ind}) is not complete, missing files: {report}"
            )
