from biom.primitives import bioproj

__all__ = ["from_bioexp", "to_expind"]


def from_bioexp(experiment: bioproj.Experiment, *, sort_attributes: bool = True, separator: str = "+") -> str:
    """
    Converts a bioproj experiment into a human-readable descriptor that can be used for
    nf-core/rnaseq pipeline design files.

    :param experiment: The experiment object to be converted.
    :param sort_attributes: Whether to sort the attributes alphabetically.
    :param separator: The separator to use between the experiment index and the attributes.
    :return: The descriptor string.
    """
    attributes: list[str] = []
    items = sorted(experiment.sample.attributes.items()) if sort_attributes else experiment.sample.attributes.items()
    for atrribute, value in items:
        attributes.append(f"{atrribute}={value}")

    # If no tags are present, use the sample description as a descriptor
    if len(attributes) == 0:
        if experiment.sample.description is None:
            raise ValueError("No tags or description available for the sample.")
        attributes.append(experiment.sample.description)

    all_attributes = "_".join(attributes)

    descriptor = f"{experiment.ind}{separator}{all_attributes}"

    return descriptor


def to_expind(descriptor: str, separator: str = "+") -> str:
    """
    Extracts the experiment index from a descriptor string.

    :param descriptor: The descriptor string.
    :return: The experiment index.
    """
    return descriptor.split(separator)[0]
