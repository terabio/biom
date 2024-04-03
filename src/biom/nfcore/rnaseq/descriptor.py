from biom import bioproj

__all__ = ["from_bioexp", "to_expind"]


def from_bioexp(experiment: bioproj.Experiment, separator: str = "+") -> str:
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

    # If no tags are present, use the sample description as a descriptor
    if len(tags) == 0:
        if experiment.sample.description is None:
            raise ValueError("No tags or description available for the sample.")
        tags.append(experiment.sample.description)

    if experiment.sample.replicate is not None:
        tags.append(str(experiment.sample.replicate))

    alltags = "_".join(tags)

    descriptor = f"{experiment.ind}{separator}{alltags}"

    return descriptor


def to_expind(descriptor: str, separator: str = "+") -> str:
    """
    Extracts the experiment index from a descriptor string.

    :param descriptor: The descriptor string.
    :return: The experiment index.
    """
    return descriptor.split(separator)[0]
