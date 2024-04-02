from biom import bioproj

__all__ = ["from_bioexp", "to_expind"]


def from_bioexp(experiment: bioproj.Experiment) -> str:
    """
    Converts a bioproj experiment into a human-readable descriptor that can be used for
    nf-core/rnaseq pipeline design files.

    :param experiment: The experiment object to be converted.
    :return: The descriptor string.
    """
    tags: list[str] = []
    for tag in [
        experiment.sample.source, experiment.sample.genotype, experiment.sample.treatment, experiment.sample.tags
    ]:
        tags.extend(tag)
    alltags = "_".join(tags)

    descriptor = f"{experiment.ind}${alltags}"
    if experiment.sample.replicate is not None:
        descriptor += f"${experiment.sample.replicate}"

    return descriptor


def to_expind(descriptor: str) -> str:
    """
    Extracts the experiment index from a descriptor string.

    :param descriptor: The descriptor string.
    :return: The experiment index.
    """
    return descriptor.split("$")[0]
