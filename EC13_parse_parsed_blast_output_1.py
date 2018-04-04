from sys import argv
import pandas as pd
import itertools

script, input_parsed_blast, output_file = argv

blastdata = pd.read_table(input_parsed_blast, sep='\t')
blast_contigs = list(blastdata.iloc[:,1])

print blastdata

lengths = []
coverages = []
for i in blast_contigs:
	contig = i.strip().split('_')
	lengths.append(contig[3])
	coverages.append(contig[5])

lenseries = pd.Series(lengths)
covseries = pd.Series(coverages)

blastdata.insert(4, column="Coverage", value=coverages)
blastdata.insert(4, column="Length", value=lengths)

blastdata.to_csv(path_or_buf=output_file, sep='\t')