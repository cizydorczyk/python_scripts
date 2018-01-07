from sys import argv
import pandas as pd
import itertools

script, input_gene_presence_absence_rtab, output_ped2, output_map2, output_genes_key = argv

input_genes = pd.read_csv(filepath_or_buffer=input_gene_presence_absence_rtab, sep='\t', header=0, index_col=0)
input_genes_2 = pd.read_csv(filepath_or_buffer=input_gene_presence_absence_rtab, sep='\t', header=0, index_col=0)

genes1 = input_genes.transpose()
genes2 = input_genes.transpose()

rownames = list(genes1.index)

colnames = list(genes1.dtypes.index)

colnames_dict = {}
num_colnames = []
for i, j in itertools.izip(colnames, range(1, len(colnames)+1)):
	colnames_dict[j] = i
	num_colnames.append(j)

newcolnames1 = []
for i in num_colnames:
	newname = float(i) + 0.1
	newcolnames1.append(newname)

newcolnames2 = []
for i in num_colnames:
	newname = float(i) + 0.2
	newcolnames2.append(newname)

genes1.columns = newcolnames1
genes2.columns = newcolnames2

frames = [genes1, genes2]
genes3 = pd.concat(frames, axis=1)

genes4 = genes3.reindex_axis(sorted(genes3.columns), axis=1)

genes4.replace([0, 1], ["A", "G"], inplace=True)

colnames2 = list(genes4.dtypes.index)

# Get and drop identical columns:
def get_identical_columns(column_vector):
    key_set = set(column_vector)
    if len(key_set) == 1:
        outcome = 'drop'
    else:
        outcome = 'keep'
    return outcome

genes4_drop_or_not = genes4.apply(get_identical_columns, axis=0)

to_drop_identical_cols = []
for i, j in genes4_drop_or_not.iteritems():
    if j == 'drop':
        to_drop_identical_cols.append(i)

genes4.drop(to_drop_identical_cols, axis=1, inplace=True)
print "Number of identical columns to remove: ", len(to_drop_identical_cols)
print "Number of identically distributed variants to remove: ", len(to_drop_identical_cols)/2
print "Genes df shape: ", genes4.shape

# Get and drop multiallelic columns:
def get_multiallelic_columns(column_vector):
    key_set = set(column_vector)
    if len(key_set) > 2:
        outcome = 'drop'
    else:
        outcome = 'keep'
    return outcome

out3 = genes4.apply(get_multiallelic_columns, axis=0)

to_drop_multiallelics = []
for i,j in out3.iteritems():
    if j == 'drop':
        to_drop_multiallelics.append(i)

genes4.drop(to_drop_multiallelics, axis=1, inplace=True)
print "Number of multi-allelic columns to remove: ", len(to_drop_multiallelics)
print "Number of multi-allelic variants to remove: ", len(to_drop_multiallelics)/2
print "Genes df shape: ", genes4.shape

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

out4 = genes4.apply(get_singleton_columns, axis=0)

to_drop_singletons = []
for i, j in out4.iteritems():
    if j == 'drop':
        to_drop_singletons.append(i)

genes4.drop(to_drop_singletons, axis=1, inplace=True)
print "Number of singleton columns to remove: ", len(to_drop_singletons)
print "Number of singleton variants to remove: ", len(to_drop_singletons)/2
print "Genes df shape: ", genes4.shape

# Generate ped/map files with all fields (dummy):
print "Generating complete (dummy) PED file..."

sex = ["1"]*155
phenotype =  ["-9"]*155

genes4.insert(0, "phenotype", value=phenotype)
genes4.insert(0, "sex", value=sex)
pid = ["0"]*155
mid = ["0"]*155
genes4.insert(0, "pid", value=pid)
genes4.insert(0, "mid", value=mid)
genes4.insert(0, "isolate_id", value=rownames)
genes4.insert(0, "family_id", value=rownames)

print genes4.shape

print "Generating complete (dummy) MAP file..."

colnames3 = list(genes4.dtypes.index)
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

genes4.to_csv(path_or_buf=output_ped2, index=False, header=False, sep='\t')
mapfile2.to_csv(path_or_buf=output_map2, index=False, header=False, sep='\t')

print "Writing genes key..."

with open(output_genes_key, 'w') as outfile1:
	for i in colnames_dict:
		outfile1.write(str(i) + '\t' + str(colnames_dict[i]) + '\n')

print "Done."