from . import eisland, misc, roi, site

# for file in path.EDITS_PER_SITE.iterdir():
#     df = pd.read_csv(file)
#     df['contig'] = df['contig'].astype(str)
#
#     mask = [key in included_esites for key in zip(df['contig'], df['pos'], df['trstrand'])]
#     df = df[mask]
#     assert len(df) == len(included_esites), set(included_esites) - set(zip(df['contig'], df['pos'], df['trstrand']))
#     genotype, _ = misc.parse_name(file.name)
#
#     edits = np.where(df['trstrand'] == "+", df['G'], df['C'])
#     # coverage = df['A'] + df['C'] + df['G'] + df['T']
#     coverage = np.where(df['trstrand'] == "+", df['A'] + df['G'], df['T'] + df['C'])
#     for contig, pos, trstrand, cov, ed in zip(df['contig'], df['pos'], df['trstrand'], coverage, edits):
#         included_esites[(contig, pos, trstrand)].append((genotype, SiteStat(cov, ed)))
#
# result = []
# for (contig, pos, trstrand), edits in included_esites.items():
#     p150KO = [stat for (genotype, stat) in edits if genotype == "p150-KO"]
#     WT = [stat for (genotype, stat) in edits if genotype == "WT"]
#     assert len(WT) == len(p150KO) == 3
#     if np.mean([x.edits / x.coverage for x in WT]) > mean_freq_threshold or \
#             np.mean([x.edits / x.coverage for x in p150KO]) > mean_freq_threshold:
#         result.append(EditingSite(contig, pos, trstrand, tuple(WT), tuple(p150KO)))
# return result
