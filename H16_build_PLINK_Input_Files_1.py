import pandas as pd
import numpy as np
from sys import argv


# Input is a binary_biallelic_csv file from first running "filter_snp_alignment_2.py" then "biallelic_fasta_to_binary_matrix_1.py":
script, input_csv, output_ped, output_map, output_ped2, output_map2 = argv

print "Generating PED file..."

snps = pd.read_csv(input_csv, sep='\t', index_col=0)
snps2 = pd.read_csv(input_csv, sep='\t', index_col=0)

sex = ["1"]*155

phenotype =  ["-9"]*155

rownames = list(snps.index)

colnames = list(snps.dtypes.index)

newcolnames1 = []
for i in colnames:
	newname = float(i) + 0.1
	newcolnames1.append(newname)

newcolnames2 = []
for i in colnames:
	newname = float(i) + 0.2
	newcolnames2.append(newname)

snps.columns = newcolnames1
snps2.columns = newcolnames2

frames = [snps, snps2]
snps3 = pd.concat(frames, axis=1)

snps4 = snps3.reindex_axis(sorted(snps3.columns), axis=1)
snps4.replace([0, 1], ["A", "G"], inplace=True)

snps5 = snps3.reindex_axis(sorted(snps3.columns), axis=1)
snps5.replace([0, 1], ["A", "G"], inplace=True)

snps4.insert(0, "phenotype", value=phenotype)
snps4.insert(0, "sex", value=sex)
snps4.insert(0, "isolate", value=rownames)

print "Generating MAP file..."

mapfile = pd.DataFrame()
mapfile.insert(0, "bp_position", value=colnames)
mapfile.insert(0, "var_id", value=colnames)
mapfile.insert(0, "chrom", value=["1"]*len(colnames))
# print mapfile

# snps4.to_csv(path_or_buf=output_ped, index=False, header=False, sep='\t')
# mapfile.to_csv(path_or_buf=output_map, index=False, header=False, sep='\t')

# Generate ped/map files with all fields (dummy):
print "Generating complete (dummy) PED file..."

snps5.insert(0, "phenotype", value=phenotype)
snps5.insert(0, "sex", value=sex)
pid = ["0"]*155
mid = ["0"]*155
snps5.insert(0, "pid", value=pid)
snps5.insert(0, "mid", value=mid)
snps5.insert(0, "isolate_id", value=rownames)
snps5.insert(0, "family_id", value=rownames)

print "Generating complete (dummy) MAP file..."

mapfile2 = pd.DataFrame()
mapfile2.insert(0, "bp_position", value=colnames)
pos = ["0"]*len(colnames)
mapfile2.insert(0, "pos", value=pos)
mapfile2.insert(0, "var_id", value=colnames)
mapfile2.insert(0, "chrom", value=["1"]*len(colnames))

snps5.to_csv(path_or_buf=output_ped2, index=False, header=False, sep='\t')
mapfile2.to_csv(path_or_buf=output_map2, index=False, header=False, sep='\t')

print "Done."