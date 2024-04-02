from attrs import define, field


@define(slots=True, frozen=True, eq=True, order=True, repr=True, hash=True)
class Sample:
    """
    A class to represent a sample in a biological experiment.

    Attributes
    ----------
    ind : str
        Sample ID, must be unique within the project.
    organism : tuple[str, ...]
        Organism(s) from which the sample was derived. E.g. ('Homo sapiens',), ('Mus musculus', 'HSV-1'), etc.
    source: tuple[str, ...]
        Cells/tissues used to generate the sample. E.g. ('HeLa',), ('brain', 'neurons'), etc.
    replicate: str
        Replicate number or identifier.
    genotype: tuple[str, ...]
        Genotype of the sample, optional. E.g. ['WT'], ['MDA5-KO', 'RIGI-KO'], etc.
    treatment: tuple[str, ...]
        Treatments applied to the sample, optional. E.g. ['Untreated'], ['DMSO', 'IFNb'], etc.
    tags: tuple[str, ...]
        Additional descriptive tags for the sample, optional. E.g. ['Confluence 75%'], etc.
    description: str
        A description of the sample, if available.
    """
    ind: str = field()
    organism: tuple[str, ...] = field(converter=lambda x: tuple(x))
    # Optional meta information
    source: tuple[str, ...] = field(converter=lambda x: tuple(x), factory=tuple)
    genotype: tuple[str, ...] = field(converter=lambda x: tuple(x), factory=tuple)
    treatment: tuple[str, ...] = field(converter=lambda x: tuple(x), factory=tuple)
    replicate: str | None = field(default=None)
    tags: tuple[str, ...] = field(converter=lambda x: tuple(x), factory=tuple)
    description: str | None = field(default=None)

    @ind.validator
    def check_ind(self, _, value):
        if not value:
            raise ValueError("Sample ID must be specified")

    @organism.validator
    def check_organism(self, _, value):
        if not value:
            raise ValueError("Organism must be specified")
