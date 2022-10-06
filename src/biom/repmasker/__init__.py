from pathlib import Path

from .repmasker import RepmaskerClassification

GRCh38 = RepmaskerClassification(Path(__file__).parent.joinpath("GRCh38.tsv.gz"))
GRCm39 = RepmaskerClassification(Path(__file__).parent.joinpath("GRCm39.tsv.gz"))
