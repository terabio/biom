from enum import Enum
from pathlib import Path

from attrs import define, field, converters


class SeqLayout(Enum):
    Paired = "paired"
    Single = "single"

    def __eq__(self, other) -> bool:
        if isinstance(other, SeqLayout):
            return self.value == other.value
        elif isinstance(other, str):
            return self.value == other
        else:
            return False

    def __repr__(self) -> str:
        return f"SeqLayout({self.value})"


@define(slots=True, frozen=True, eq=True, order=True, repr=True, hash=True)
class SeqRun:
    """
    A class to represent a sequencing run in a biological experiment.

    Attributes
    ----------
    ind : str
        Index of the sequencing run, should be unique within the project.
    machine : str
        Sequencing machine used for the run. E.g. 'Illumina NovaSeq 6000', 'Oxford Nanopore MinION', etc.
    layout : SeqLayout
        Layout of the sequencing run. E.g. SeqLayout.Paired, SeqLayout.Single or just 'paired', 'single'.
    length : int
        Average read length of the sequencing run. None if the length is unknown.
    files : tuple[Path, ...]
        Paths to the sequencing files. Should be one (single) or two files (paired) depending on the layout.
    description : str
        Description of the sequencing run, if available.
    """
    ind: str = field()
    machine: str = field()
    layout: SeqLayout = field(converter=SeqLayout)
    files: tuple[Path, ...] = field(converter=lambda x: tuple(Path(f) for f in x))
    # Optional metadata
    reads: int | None = field(default=None, converter=converters.optional(lambda x: int(x)))
    bases: int | None = field(default=None, converter=converters.optional(lambda x: int(x)))
    description: str | None = field(default=None)

    @ind.validator
    def check_ind(self, _, value):
        if not value:
            raise ValueError("Sequencing run index must be specified")

    @machine.validator
    def check_machine(self, _, value):
        if not value:
            raise ValueError("Sequencing machine must be specified")

    @files.validator
    def check_files(self, _, value):
        if not value:
            raise ValueError("Sequencing files must be specified")

    @layout.validator
    def check_layout(self, _, value):
        if not value:
            raise ValueError("Sequencing layout must be specified")
        match value:
            case SeqLayout.Paired:
                if len(self.files) != 2:
                    raise ValueError("Paired-end sequencing requires two files")
            case SeqLayout.Single:
                if len(self.files) != 1:
                    raise ValueError("Single-end sequencing requires one file")
            case _:
                raise ValueError(f"Unknown sequencing layout: {value}")

    @reads.validator
    def check_reads(self, _, value):
        if value is not None and (value <= 0 or not isinstance(value, int)):
            raise ValueError("Total number of reads must be a positive integer")

    @bases.validator
    def check_bases(self, _, value):
        if value is not None and (value <= 0 or not isinstance(value, int)):
            raise ValueError("Total number of bases must be a positive integer")

    @description.validator
    def check_description(self, _, value):
        if value is not None and not value:
            raise ValueError("If specified, description must be non-empty. Use None to indicate lack of description")
