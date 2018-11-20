from sys import argv
import pandas as pd
import itertools

script, input_gene_presence_absence_rtab, output_ped, output_map, output_matrix, output_genes_key = argv

## Load pangenome:
input_genes = pd.read_csv(filepath_or_buffer=input_gene_presence_absence_rtab, sep='\t', header=0, index_col=0)

genes1 = input_genes.transpose()

########################### Define required functions ################################################

# Get and drop columns with N:
def get_n_columns(column_vector):
    key_set = set(column_vector)
    if 'N' in key_set:
        outcome = "drop"
    else:
        outcome = "keep"
    return outcome

# Get and drop identical columns:
def get_identical_columns(column_vector):
    key_set = set(column_vector)
    if len(key_set) == 1:
        outcome = 'drop'
    else:
        outcome = 'keep'
    return outcome

# Get and drop multiallelic columns:
def get_multiallelic_columns(column_vector):
    key_set = set(column_vector)
    if len(key_set) > 2:
        outcome = 'drop'
    else:
        outcome = 'keep'
    return outcome

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

# Function to drop columns after running above functions:
def drop_columns(df, possible_col_to_drop):
	to_drop = []
	for i, j in possible_col_to_drop.iteritems():
		if j == 'drop':
			to_drop.append(i)
	df.drop(to_drop, axis=1, inplace=True)
	return df

# Convert df to binary format with abscence = "A" for one allele (ref base) and presence = "G" for the other (alt base):
def convert_to_binary(column_vector):
    key_set = list(set(column_vector))
    allele_a = key_set[0]
    allele_b = key_set[1]
    column_vector1 = column_vector.replace(allele_a, "A")
    column_vector2 = column_vector1.replace(allele_b, "G")
    return column_vector2

def indels_to_binary(column_vector):
	key_set = set(column_vector)
	for key in key_set:
		if key == "A":
			# column_vector.replace([key], ["A"], inplace=True)
			pass
		else:
			column_vector.replace([key], ["G"], inplace=True)
	return column_vector

############################ Recode variants with unique numbers ####################################

## Code accessory genes with a unique number:

# Pangenome:
pangenome_guide_dict = {}
pangenome_guide_list = []
pangenome_colnames = list(genes1.dtypes.index)
pangenome_numbers = range(1, len(pangenome_colnames)+1)

for i, j in itertools.izip(pangenome_colnames, pangenome_numbers):
	pangenome_guide_dict[i] = j
	pangenome_guide_list.append(j)

# Write gene key to file:
print "Writing variants key..."

with open(output_genes_key, 'w') as outfile1:
	for j in pangenome_guide_dict:
		outfile1.write(str(j) + '\t' + str(pangenome_guide_dict[j]) + '\n')

# Rename variants:
genes1.rename(columns=pangenome_guide_dict, inplace=True)

########################### Convert variant dataframes to binary format and filter to remove ############################
########################### identical columns, singletons, multiallelic columns, columns with N, etc. ############################

## Pangenome:

print "Number of genes before removing identical, multiallelic, and singleton columns: ", len(genes1.columns)

# Make binary:
genes1.replace([0, 1], ["A", "G"], inplace=True)

# Remove identical columns (there should be some for the pangenome):
possible_identical_col = genes1.apply(get_identical_columns, axis=0)
genes1 = drop_columns(genes1, possible_identical_col)
print "Number of genes after removing identical columns: ", len(genes1.columns)

# Remove multiallelic columns (shouldn't be any for pangenome):
possible_multiallelics = genes1.apply(get_multiallelic_columns, axis=0)
genes1 = drop_columns(genes1, possible_multiallelics)
print "Number of genes after removing multiallelic columns: ", len(genes1.columns)

# Remove singleton columns (some may be present):
possible_singletons = genes1.apply(get_singleton_columns, axis=0)
genes1 = drop_columns(genes1, possible_singletons)
print "Number of genes after removing singleton columns: ", len(genes1.columns)

print "Final number of genes: ", len(genes1.columns)

############ Format and write variants dataframe to file (for treeWAS/other software requiring basic matrix) ##########
print "Writing variants matrix to file..."
variants_matrix = genes1.replace(["A", "G"], [0, 1])

variants_matrix.to_csv(path_or_buf=output_matrix, sep='\t')

#################### Generate ped/map files #########################
genes2 = genes1.copy()

# Rename columns to allow for easy sorting:
gene_colnames = genes1.columns

genes1_new_colnames = []
genes2_new_colnames = []

for colname in gene_colnames:
    genes1_new_colnames.append(float(colname) + 0.1)
    genes2_new_colnames.append(float(colname) + 0.2)

genes1.columns = genes1_new_colnames
genes2.columns = genes2_new_colnames

# Combine and sort the two genes dfs:
frames = [genes1, genes2]
genes_final = pd.concat(frames, axis=1)

genes_final = genes_final.reindex(columns=sorted(genes_final.columns))

# Generate ped file:
print "Generating PED file..."

# Create fields for PED file:
isolate_family_ids = list(genes_final.index)
sex = ["0"]*155
phenotype =  ["-9"]*155
pid = ["0"]*155
mid = ["0"]*155

genes_final.insert(0, "phenotype", value=phenotype)
genes_final.insert(0, "sex", value=sex)
genes_final.insert(0, "pid", value=pid)
genes_final.insert(0, "mid", value=mid)
genes_final.insert(0, "isolate_id", value=isolate_family_ids)
genes_final.insert(0, "family_id", value=isolate_family_ids)

print "Generating MAP file..."

colnames = list(genes_final.dtypes.index)
colnames.remove("phenotype")
colnames.remove("sex")
colnames.remove("pid")
colnames.remove("mid")
colnames.remove("isolate_id")
colnames.remove("family_id")

newcolnames = []
for i in colnames:
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

######################## Write ped/map and variants key to file #####################################

genes_final.to_csv(path_or_buf=output_ped, index=False, header=False, sep='\t')
mapfile2.to_csv(path_or_buf=output_map, index=False, header=False, sep='\t')

print "Done."