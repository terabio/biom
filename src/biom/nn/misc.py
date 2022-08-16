from typing import List

from pybedtools import Interval
from math import ceil


def makewindows(interval: Interval, wsize: int, overlapfrac: float) -> List[Interval]:
    assert interval.length >= wsize

    step = int(wsize * (1 - 2 * overlapfrac))
    total = int(ceil(interval.length / step))
    result = []
    for i in range(total):
        start = step * i
        end = start + wsize

        assert end - start == wsize
        start, end = start + interval.start, end + interval.start

        # Make sure that we are inside the interval
        if end > interval.end:
            offset = end - interval.end
            start -= offset
            end -= offset
        assert end - start == wsize and start >= interval.start and end <= interval.end

        result.append(Interval(interval.chrom, start, end))
    assert result[0].start == interval.start and result[-1].end == interval.end
    return result
