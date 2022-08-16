import json
import os
import subprocess
import tempfile
from pathlib import Path
from typing import List

import pysam


def fragments(bam: Path, threads: int):
    flagstat = pysam.flagstat('-@', str(threads), '-O', 'json', bam.as_posix())
    return json.loads(flagstat)['QC-passed reads']['properly paired'] / 2


def merge(bams: List[Path], threads: int) -> Path:
    if len(bams) == 1:
        return bams[0]
    # Merge bam files in a TMP file
    fd, merged = tempfile.mkstemp()
    os.close(fd)

    # Calculate coverage
    bams = [x.as_posix() for x in bams]
    pysam.merge('-O', 'BAM', '-@', str(threads), '-f', merged, *bams)
    subprocess.run(['samtools', 'index', "-@", str(threads), merged], check=True)
    return Path(merged)


def coverage(bams: List[Path], title: str, factor: float, saveto: Path, threads: int):
    allbam = None

    for strand in 'forward', 'reverse':
        o = saveto.joinpath(f"{title}.{strand}.bigWig")
        # Cache results
        if o.exists():
            continue

        # Merge only if needed
        if allbam is None:
            allbam = merge(bams, threads)
        assert allbam is not None

        command = [
            'bamCoverage', '-p', str(threads), '--scaleFactor', str(factor), '--filterRNAstrand', strand,
            '--skipNonCoveredRegions', '--minMappingQuality', '1', '-o', o.as_posix(), "--samFlagInclude", "3",
            "--samFlagExclude", "2816", '-b', allbam.as_posix(), '--binSize', '5'
        ]
        subprocess.run(command, check=True)

    # Remove merged BAM if needed
    if allbam is not None and len(bams) > 1:
        os.unlink(allbam.as_posix())
        os.unlink(allbam.as_posix() + ".bai")


def cpm(bams: List[Path], title: str, saveto: Path, threads: int):
    factor, allbam = None, None

    for strand in 'forward', 'reverse':
        o = saveto.joinpath(f"{title}.{strand}.bigWig")
        # Cache results
        if o.exists():
            continue

        # Calculate factor & merge only if needed
        if factor is None:
            allbam = merge(bams, threads)
            total = fragments(allbam, threads)
            factor = 1_000_000 / total
        assert allbam is not None and factor is not None

        command = [
            'bamCoverage', '-p', str(threads), '--scaleFactor', str(factor), '--filterRNAstrand', strand,
            '--skipNonCoveredRegions', '--minMappingQuality', '1', '-o', o.as_posix(), "--samFlagInclude", "3",
            "--samFlagExclude", "2816", '-b', allbam.as_posix(), '--binSize', '5'
        ]
        subprocess.run(command, check=True)

    # Remove merged BAM if needed
    if allbam is not None and len(bams) > 1:
        os.unlink(allbam.as_posix())
        os.unlink(allbam.as_posix() + ".bai")
