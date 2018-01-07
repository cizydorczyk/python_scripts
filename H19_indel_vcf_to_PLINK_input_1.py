from sys import argv
import pandas as pd

script, input_vcf, output_ped2, output_map2 = argv

input_indels = pd.read_csv(filepath_or_buffer=input_vcf, sep='\t', header=0, skiprows=35, index_col=1)
input_indels_2 = pd.read_csv(filepath_or_buffer=input_vcf, sep='\t', header=0, skiprows=35, index_col=1)

del input_indels['#CHROM'] #, 'ID', 'REF', 'ALT', 'QUAL', 'FILTER', 'INFO', 'FORMAT']
del input_indels['ID']
del input_indels['REF']
del input_indels['ALT']
del input_indels['QUAL']
del input_indels['FILTER']
del input_indels['INFO']
del input_indels['FORMAT']

del input_indels_2['#CHROM'] #, 'ID', 'REF', 'ALT', 'QUAL', 'FILTER', 'INFO', 'FORMAT']
del input_indels_2['ID']
del input_indels_2['REF']
del input_indels_2['ALT']
del input_indels_2['QUAL']
del input_indels_2['FILTER']
del input_indels_2['INFO']
del input_indels_2['FORMAT']

indels3 = input_indels.transpose()
indels4 = input_indels_2.transpose()

rownames = list(indels3.index)

colnames = list(indels3.dtypes.index)

newcolnames1 = []
for i in colnames:
	newname = float(i) + 0.1
	newcolnames1.append(newname)

newcolnames2 = []
for i in colnames:
	newname = float(i) + 0.2
	newcolnames2.append(newname)

indels3.columns = newcolnames1
indels4.columns = newcolnames2

frames = [indels3, indels4]
indels5 = pd.concat(frames, axis=1)

indels6 = indels5.reindex_axis(sorted(indels5.columns), axis=1)

indels6.replace(["."], ["A"], inplace=True)

indels_dict = indels6.to_dict(orient="list")

out_dict = {}
for key, value in indels_dict.iteritems():
	out_list = []
	for item in value:
		if item == "A":
			out_list.append("A")
		else:
			out_list.append("G")
	out_dict[key] = out_list

indels7 = pd.DataFrame.from_dict(out_dict)

# Get and drop identical columns:
def get_identical_columns(column_vector):
    key_set = set(column_vector)
    if len(key_set) == 1:
        outcome = 'drop'
    else:
        outcome = 'keep'
    return outcome

drop_or_not = indels7.apply(get_identical_columns, axis=0)

to_drop_identical_cols = []
for i, j in drop_or_not.iteritems():
    if j == 'drop':
        to_drop_identical_cols.append(i)

indels7.drop(to_drop_identical_cols, axis=1, inplace=True)
print "Number of identical columns to remove: ", len(to_drop_identical_cols)
print "Number of identically distributed variants to remove: ", len(to_drop_identical_cols)/2
print "Genes df shape: ", indels7.shape

# Get and drop multiallelic columns:
def get_multiallelic_columns(column_vector):
    key_set = set(column_vector)
    if len(key_set) > 2:
        outcome = 'drop'
    else:
        outcome = 'keep'
    return outcome

out3 = indels7.apply(get_multiallelic_columns, axis=0)

to_drop_multiallelics = []
for i,j in out3.iteritems():
    if j == 'drop':
        to_drop_multiallelics.append(i)

indels7.drop(to_drop_multiallelics, axis=1, inplace=True)
print "Number of multi-allelic columns to remove: ", len(to_drop_multiallelics)
print "Number of multi-allelic variants to remove: ", len(to_drop_multiallelics)/2
print "Genes df shape: ", indels7.shape

# Get and drop singleton columns:
def get_singleton_columns(column_vector):
    key_set = set(column_vector)
    if len(key_set) == 2:
        counts = list(column_vector.value_counts())
        if 1 in counts:
            outcome = "drop"
        else:
            outcome = "keep"
    else:
        print "Something isn't right"
    return outcome

out4 = indels7.apply(get_singleton_columns, axis=0)

to_drop_singletons = []
for i, j in out4.iteritems():
    if j == 'drop':
        to_drop_singletons.append(i)

indels7.drop(to_drop_singletons, axis=1, inplace=True)
print "Number of singleton columns to remove: ", len(to_drop_singletons)
print "Number of singleton variants to remove: ", len(to_drop_singletons)/2
print "Genes df shape: ", indels7.shape

# print "Length of dict before removing indels identical among all isolates: ", len(out_dict)

# to_drop = []
# for key, value in out_dict.iteritems():
# 	key_set = set(value)
# 	if len(key_set) == 1:
# 		to_drop.append(key)
# 	elif len(key_set) > 2:
# 		print "something isn't right..."

# print "Number of columns to drop (double the number of indels): ", len(to_drop)

# for i in to_drop:
# 	del out_dict[i]

# print "Length of dict after removing indels identical among all isolates: ", len(out_dict)

# indels7 = pd.DataFrame.from_dict(out_dict)

# indels7.insert(0, "isolates", value=rownames)
# indels7.set_index("isolates", inplace=True)

# Generate ped/map files with all fields (dummy):
print "Generating complete (dummy) PED file..."

sex = ["1"]*155
phenotype =  ["-9"]*155

indels7.insert(0, "phenotype", value=phenotype)
indels7.insert(0, "sex", value=sex)
pid = ["0"]*155
mid = ["0"]*155
indels7.insert(0, "pid", value=pid)
indels7.insert(0, "mid", value=mid)
indels7.insert(0, "isolate_id", value=rownames)
indels7.insert(0, "family_id", value=rownames)

print indels7.shape

print "Generating complete (dummy) MAP file..."

colnames3 = list(indels7.dtypes.index)
colnames3.remove("phenotype")
colnames3.remove("sex")
colnames3.remove("pid")
colnames3.remove("mid")
colnames3.remove("isolate_id")
colnames3.remove("family_id")

newcolnames = []
for i in colnames3:
	if ".2" in str(i):
		pass
	else:
		newcolnames.append(str(i)[0:-2])

mapfile2 = pd.DataFrame()
mapfile2.insert(0, "bp_position", value=newcolnames)
pos = ["0"]*len(newcolnames)
mapfile2.insert(0, "pos", value=pos)
mapfile2.insert(0, "var_id", value=newcolnames)
mapfile2.insert(0, "chrom", value=["1"]*len(newcolnames))

indels7.to_csv(path_or_buf=output_ped2, index=False, header=False, sep='\t')
mapfile2.to_csv(path_or_buf=output_map2, index=False, header=False, sep='\t')

print "Done."