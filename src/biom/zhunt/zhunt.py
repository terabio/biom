import os
import subprocess
import tempfile
from pathlib import Path
from subprocess import DEVNULL, PIPE

import numpy as np
import pandas as pd

ROOT = Path(__file__).parent
EXECUTABLE = ROOT.joinpath("zhunt3-alan").as_posix()
NEXECUTABLE = ROOT.joinpath("zhunt3-optimized").as_posix()

# def fetch(contig: str, start: int, end: int) -> np.ndarray:
#     bw = pyBigWig.open(INDEX)
#     return bw.values(ensembl.contigs.to_UCSC(contig), start, end, numpy=True)

def _run(query: str, windowsize, minsize, maxsize) -> Path:
    fd, temp = tempfile.mkstemp()
    os.close(fd)
    with open(temp, 'w') as stream:
        stream.write(query)

    subprocess.run([EXECUTABLE, str(windowsize), str(minsize), str(maxsize), temp],
                   check=True, stdout=subprocess.PIPE, stderr=subprocess.DEVNULL,
                   input=query, encoding='ascii')
    file = Path(temp + ".Z-SCORE")
    assert file.exists()
    os.unlink(temp)
    return file


def predict(query: str, windowsize: int = 6, minsize: int = 3, maxsize: int = 6):
    file = _run(query, windowsize, minsize, maxsize)
    with open(file, 'r') as stream:
        df = pd.read_csv(stream,
                         names=['Start', 'End', 'nu-1', 'nu-2', 'nu-3', 'ZH-Score', 'Sequence', 'Conformation'],
                         usecols=['Start', 'End', 'ZH-Score', 'Sequence', 'Conformation'],
                         dtype={"Start": np.int32, "End": np.int32, "ZH-Score": np.float32,
                                "Sequence": str, "Conformation": str},
                         skiprows=1, sep='\s+')
    file.unlink()
    return df


def iterpredict(query: str, windowsize: int = 6, minsize: int = 3, maxsize: int = 6):
    cmd = [NEXECUTABLE, str(windowsize), str(minsize), str(maxsize)]
    popen = subprocess.Popen(cmd, stdout=PIPE, stdin=PIPE, stderr=DEVNULL)
    popen.stdin.write(query.encode("ascii"))
    popen.stdin.close()

    iterator = iter(popen.stdout.readline, "")
    length = int(next(iterator).decode().split()[0])
    assert length == len(query)

    yielded = 0
    for line in iterator:
        try:
            line = float(line)
        except ValueError:
            break
        yield line
        yielded += 1
    popen.stdout.close()

    assert yielded == length == len(query), \
        f"For query {query} \n Z-score has length {yielded} vs query length {len(query)}"
    return_code = popen.wait()
    if return_code:
        raise subprocess.CalledProcessError(return_code, cmd)
