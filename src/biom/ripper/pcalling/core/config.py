from dataclasses import dataclass
from pathlib import Path
from typing import Dict, List, Optional, Tuple

import numpy as np


@dataclass(frozen=True)
class Scaling:
    treatment: np.float32
    control: np.float32


@dataclass(frozen=True)
class PeakCallingConfig:
    @dataclass(frozen=True)
    class Saveto:
        title: str
        pileup: Optional[Path]
        enrichment: Optional[Path]
        pvpeaks: Optional[Path]
        fdrpeaks: Optional[Path]

    @dataclass(frozen=True)
    class ProcessingParams:
        # Stranding (only f/s is supported)
        stranding: str
        # Scaling options
        scaling: Scaling
        # Extension size for each contig (both transcriptomic & genomic coordinates)
        extsize: Dict[str, List[int]]
        # Number of reads to parallelize the pileup computation
        threads: int
        backend: str
        # Bam flags
        inflags: int
        exflags: int
        minmapq: int

    @dataclass()
    class PeakCallingParams:
        qvcutoff: Optional[float] = 0.01
        pvcutoff: Optional[float] = None
        fecutoff: Optional[float] = 2
        minsize: int = 50
        maxgap: int = 25
        mintrtfrag: int = 10

        def __post_init__(self):
            assert self.qvcutoff or self.pvcutoff or self.fecutoff

    # Input bam files
    treatment: List[Path]
    control: List[Path]
    # Target contigs
    contigs: Optional[Tuple[str]]
    # Genome size
    geffsize: int

    process: ProcessingParams
    callp: PeakCallingParams
    saveto: Saveto

    def __post_init__(self):
        assert all(x.is_file() for x in self.treatment + self.control)
