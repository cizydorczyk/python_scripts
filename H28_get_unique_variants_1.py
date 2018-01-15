from sys import argv
import pandas as pd

script, input_matrix, output_matrix, output_ped, output_map, output_variants_key = argv

inputmatrix = pd.read_csv(filepath_or_buffer=input_matrix, sep='\t', index_col=0)

inputcol = list(inputmatrix.columns)
print "Original number of columns: ",len(inputcol)

colstring_dict = {}
colstring_set = set()
for i in inputcol:
	col = list(inputmatrix[i])
	col = map(str, col)
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

outmatrix1 = inputmatrix.drop(labels=to_drop, axis=1, inplace=False)

# Write new matrix file:
outmatrix1.to_csv(path_or_buf=output_matrix, sep='\t')

# Create ped/map files:
outmatrix1 = outmatrix1.replace([0, 1], ["A", "G"])

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

outmatrix3 = outmatrix3.reindex(sorted(outmatrix3.columns), axis=1)

## Generate ped/map files with all fields (dummy):
print "Generating complete (dummy) PED file..."

rownames = list(outmatrix3.index)
sex = ["1"]*155
phenotype =  ["-9"]*155
pid = ["0"]*155
mid = ["0"]*155

outmatrix3.insert(0, "phenotype", value=phenotype)
outmatrix3.insert(0, "sex", value=sex)
outmatrix3.insert(0, "pid", value=pid)
outmatrix3.insert(0, "mid", value=mid)
outmatrix3.insert(0, "isolate_id", value=rownames)
outmatrix3.insert(0, "family_id", value=rownames)

print "Generating complete (dummy) MAP file..."

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
outmatrix3.to_csv(path_or_buf=output_ped, index=False, header=False, sep='\t')
mapfile2.to_csv(path_or_buf=output_map, index=False, header=False, sep='\t')

# Write collapsed variants to file:
print "Writing variants key..."

with open(output_variants_key, 'w') as outfile1:
	for i in colstring_dict:
		out = '\t'.join(colstring_dict[i])
		outfile1.write(colstring_dict[i][0] + '\t' + out + '\n')