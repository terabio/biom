# from pathlib import Path
# from typing import List
#
# from .pileup import Pileup
#
#
# def to_bedgraph(saveto: Path, pileups: List[Pileup], nozeros: bool, ndigits: int = 5):
#     pileups = sorted(pileups, key=lambda x: x.id)
#     sensitivity = 10 ** -ndigits
#     with open(saveto, 'w') as stream:
#         stream.write("track type=bedGraph name=placeholder description=placeholder\n")
#
#         for p in pileups:
#             assert p.interend.size == p.values.size
#             start, end, val = 0, p.interend[0], p.values[0]
#             for newend, newval in zip(p.interend[1:], p.values[1:]):
#                 if abs(newval - val) > sensitivity:
#                     roundedval = round(val, ndigits)
#                     if nozeros and roundedval == 0:
#                         roundedval = sensitivity
#                         assert roundedval > 0
#                     stream.write(f"{p.id}\t{start}\t{end}\t{roundedval}\n")
#                     start = end
#                     end = newend
#                     val = newval
#                 else:
#                     end = newend
#
#             # Final region
#             roundedval = round(val, ndigits)
#             if nozeros and roundedval == 0:
#                 roundedval = sensitivity
#                 assert roundedval > 0
#             stream.write(f"{p.id}\t{start}\t{end}\t{roundedval}\n")
#             assert end == p.interend[-1]
