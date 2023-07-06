from biom.ripper.core import fragments
from pathlib import Path

deductor = fragments.strdeductors.get("f/s")
folder = Path("/home/alnfedorov/tmp")
for bam in folder.glob("*.bam"):
    print(bam)
    ab = fragments.loadfrom(
        [bam],
        deductor, contig='16', inflags=3, exflags=2564, minmapq=1
    )
