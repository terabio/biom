import tempfile
from pathlib import Path

from biom import gindex
from biom.gindex import Annotator


def make_annotator(disambigfn, annotatefn) -> Annotator:
    intervals = {
        ('1', '+'): [(1, 3), (5, 8)],
        ('2', '-'): [(0, 2), (3, 6), (9, 10), (5, 10)],
        ('2', '+'): [
            (1, 2), (3, 4), (5, 6), (7, 8),
            (1, 4), (6, 10),
            (2, 9),
            (0, 10)
        ],
    }

    trees = []
    # Make bed files (with dups intentionally)
    for keys in [('2', '-'), ('1', '+')], [('2', '+'), ('1', '+')]:
        with tempfile.NamedTemporaryFile(delete=False, mode='w') as bed:
            for contig, strand in keys:
                for ind, (start, end) in enumerate(intervals[(contig, strand)]):
                    bed.write(
                        f"{contig}\t{start}\t{end}\t{ind}\t.\t{strand}\n"
                    )
            bed.close()
            trees.append(
                gindex.from_bed(Path(bed.name))
            )
    index = gindex.merge(*trees)
    return Annotator(index, disambiguation=disambigfn, annotation=annotatefn)


# Annotation blocks:
# * Single element
# * Multiple elements

# Relevant to an interval
# * Inside interval
# * Envelop interval
# * Equal to interval

# Number of intervals
# * Single
# * Multiple

# NMS
nms = make_annotator(disambigfn="proportional", annotatefn="nms")
for contig in '1', '2':
    for strand in '+', '-':
        assert nms.annotatei(contig, strand, [(0, 100)]) == "NA"
        assert nms.annotatei(contig, strand, [(0, 3), (5, 10), (50, 100)]) == "NA"

assert nms.annotatei('1', '+', [(1, 3)]) == '0'
assert nms.annotatei('1', '+', [(0, 3)]) == '0'
assert nms.annotatei('1', '+', [(1, 4)]) == '0'
assert nms.annotatei('1', '+', [(0, 5)]) == 'NA'

assert nms.annotatei('2', '+', [(0, 10)]) == '7'
assert nms.annotatei('2', '+', [(1, 5), (5, 10)]) == '7'
assert nms.annotatei('2', '+', [(1, 2), (3, 4), (5, 11)]) == '7'

assert nms.annotatei('2', '-', [(0, 2), (3, 6)]) == '1'
assert nms.annotatei('2', '-', [(0, 6)]) == '1'
assert nms.annotatei('2', '-', [(0, 10)]) == '3'
assert nms.annotatei('2', '-', [(0, 15)]) == 'NA'

# Priority
priority = make_annotator(
    disambigfn="proportional",
    annotatefn=("priority", ['1', '3', '5', '7', '0', '2', '4', '6', 'NA'])
)
for contig, strand in ('1', '+'), ('2', '+'), ('2', '-'):
    assert priority.annotatei(contig, strand, [(0, 100)]) == "1"
    assert priority.annotatei(contig, strand, [(0, 4), (5, 10), (50, 100)]) == "1"
    assert priority.annotatei(contig, strand, [(10, 50), (55, 100)]) == "NA"

priority.annotatei('1', '+', [(0, 3)]) == '0'
priority.annotatei('1', '+', [(0, 3), (5, 8)]) == '0'
priority.annotatei('2', '+', [(0, 1), (1, 2), (2, 3)]) == '7'

# Fractional
fractional = make_annotator(
    disambigfn="proportional",
    annotatefn="frac-overlap"
)

assert fractional.annotatei('2', '-', [(0, 2)]) == {'0': 1.0}
assert fractional.annotatei('2', '-', [(0, 2), (3, 5)]) == {'0': 0.5, '1': 0.5}
assert fractional.annotatei('2', '-', [(0, 2), (3, 5), (6, 9)]) == {'0': 2 / 7, '1': 2 / 7, '3': 3 / 7}
