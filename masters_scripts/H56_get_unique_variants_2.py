from sys import argv
import pandas as pd
import argparse

parser = argparse.ArgumentParser()
parser.add_argument("allele_matrix", help="Output matrix from H56_fasta_to_allele_matrix_1.py that will be used as input.")
parser.add_argument("output_prefix", help="Prefix to append to all output files. These include the ped/map files, a binary matrix suitable \
	for treeWAS (etc.), and a variants key that serves as a guide as to which variants were collapsed.")
args = parser.parse_args()

print "Reading input matrix..."
inputmatrix_original = pd.read_csv(filepath_or_buffer=args.allele_matrix, sep='\t', index_col=0, dtype=str)

# Convert df to binary format:
def convert_to_binary(column_vector):
    key_set = list(set(column_vector))
    allele_a = key_set[0]
    allele_b = key_set[1]
    column_vector1 = column_vector.replace(allele_a, "0")
    column_vector2 = column_vector1.replace(allele_b, "1")
    return column_vector2

inputmatrix_binary = inputmatrix_original.apply(convert_to_binary, axis=0)

inputcol = list(inputmatrix_binary.columns)

print "Original number of columns: ",len(inputcol)

colstring_dict = {}
colstring_set = set()
for i in inputcol:
	col = list(inputmatrix_binary[i])
	colstring = ''.join(col)
	if colstring not in colstring_set:
		colstring_set.add(colstring)
		colstring_dict[colstring] = [i]
	else:
		colstring_dict[colstring].append(i)

print "Number of uniquely segregating columns: ", len(colstring_dict)
to_drop = []
for key in colstring_dict:
	key_value = colstring_dict[key]
	if len(key_value) > 1:
		key_value_minus_first = key_value[1:]
		for i in key_value_minus_first:
			to_drop.append(i)
	else:
		pass
print "Number of columns to remove: ", len(to_drop)

outmatrix1 = inputmatrix_original.drop(labels=to_drop, axis=1, inplace=False)
outmatrix_binary = inputmatrix_binary.drop(labels=to_drop, axis=1, inplace=False)

# Write new matrix file:
print "Generating binary matrix..."
output_matrix = args.output_prefix + ".binary_matrix.txt"
outmatrix_binary.to_csv(path_or_buf=output_matrix, sep='\t')

# Create ped/map files:

colnames_new = list(outmatrix1.columns)

colnames_new1 = []
colnames_new2 = []
for i in colnames_new:
	colnames_new1.append(float(i) + 0.1)
	colnames_new2.append(float(i) + 0.2)

outmatrix1.columns = colnames_new1

outmatrix2 = outmatrix1.copy()
outmatrix2.columns = colnames_new2

frames = [outmatrix1, outmatrix2]
outmatrix3 = pd.concat(frames, axis=1)

outmatrix3 = outmatrix3.reindex_axis(sorted(outmatrix3.columns), axis=1)

## Generate ped/map files with all fields (dummy):
print "Generating PED file..."

rownames = list(outmatrix3.index)
sex = ["0"]*155
phenotype =  ["-9"]*155
pid = ["0"]*155
mid = ["0"]*155

outmatrix3.insert(0, "phenotype", value=phenotype)
outmatrix3.insert(0, "sex", value=sex)
outmatrix3.insert(0, "pid", value=pid)
outmatrix3.insert(0, "mid", value=mid)
outmatrix3.insert(0, "isolate_id", value=rownames)
outmatrix3.insert(0, "family_id", value=rownames)

print "Generating MAP file..."

colnames3 = list(outmatrix3.dtypes.index)
colnames3.remove("phenotype")
colnames3.remove("sex")
colnames3.remove("pid")
colnames3.remove("mid")
colnames3.remove("isolate_id")
colnames3.remove("family_id")

newcolnames3 = []
for i in colnames3:
	if ".2" in str(i):
		pass
	else:
		newcolnames3.append(str(i)[0:-2])

mapfile2 = pd.DataFrame()
mapfile2.insert(0, "bp_position", value=newcolnames3)
pos = ["0"]*len(newcolnames3)
mapfile2.insert(0, "pos", value=pos)
mapfile2.insert(0, "var_id", value=newcolnames3)
mapfile2.insert(0, "chrom", value=["1"]*len(newcolnames3))

# Write ped/map to file:
output_ped = args.output_prefix + ".ped"
output_map = args.output_prefix + ".map"
outmatrix3.to_csv(path_or_buf=output_ped, index=False, header=False, sep='\t')
mapfile2.to_csv(path_or_buf=output_map, index=False, header=False, sep='\t')

# Write collapsed variants to file:
print "Writing variants key..."

outvariantskey = args.output_prefix + ".variants_key.txt"
with open(outvariantskey, 'w') as outfile1:
	for i in colstring_dict:
		out = '\t'.join(colstring_dict[i])
		outfile1.write(colstring_dict[i][0] + '\t' + out + '\n')

print "Done."