from sys import argv
import pandas as pd

script, input_vcf, output_ped, output_map, output_matrix = argv

## Load indels:
input_indels = pd.read_csv(filepath_or_buffer=input_vcf, sep='\t', header=0, skiprows=35, index_col=1)

del input_indels['#CHROM'] #, 'ID', 'REF', 'ALT', 'QUAL', 'FILTER', 'INFO', 'FORMAT']
del input_indels['ID']
del input_indels['REF']
del input_indels['ALT']
del input_indels['QUAL']
del input_indels['FILTER']
del input_indels['INFO']
del input_indels['FORMAT']

indels1 = input_indels.transpose()
indel_rownames = list(indels1.index)

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
    column_vector1 = column_vector.replace(allele_a, "0")
    column_vector2 = column_vector1.replace(allele_b, "1")
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

########################### Convert variant dataframes to binary format and filter to remove ############################
########################### identical columns, singletons, multiallelic columns, columns with N, etc. ############################

## Indels:

print "Number of indels before removing identical, multiallelic, and singleton columns: ", len(indels1.columns)

# Make binary:

indels1.replace(["."], ["A"], inplace=True)
indels1 = indels1.apply(indels_to_binary, axis=0)

# Remove identical columns (there may be some for indels):
possible_identical_col = indels1.apply(get_identical_columns, axis=0)
indels1 = drop_columns(indels1, possible_identical_col)
print "Number of indels after removing identical columns: ", len(indels1.columns)

# Remove multiallelic columns (there shouldn't be any for indels):
possible_multiallelics = indels1.apply(get_multiallelic_columns, axis=0)
indels1 = drop_columns(indels1, possible_multiallelics)
print "Number of indels after removing multiallelic columns: ", len(indels1.columns)


# Remove singleton columns (there may be some for indels):
possible_singletons = indels1.apply(get_singleton_columns, axis=0)
indels1 = drop_columns(indels1, possible_singletons)
print "Number of indels after removing singleton columns: ", len(indels1.columns)

print "Final number of indels: ", len(indels1.columns)

########## Format and write variants dataframe to file (for treeWAS/other software requiring basic matrix) ##########
print "Writing variants matrix to file..."
variants_matrix = indels1.replace(["A", "G"], [0, 1])

variants_matrix.to_csv(path_or_buf=output_matrix, sep='\t')

#################### Generate ped/map files #########################
indels2 = indels1.copy()

# Rename columns to allow for easy sorting:
indel_colnames = indels1.columns

indels1_new_colnames = []
indels2_new_colnames = []

for colname in indel_colnames:
    indels1_new_colnames.append(float(colname) + 0.1)
    indels2_new_colnames.append(float(colname) + 0.2)

indels1.columns = indels1_new_colnames
indels2.columns = indels2_new_colnames

# Combine and sort the two indels dfs:
frames = [indels1, indels2]
indels_final = pd.concat(frames, axis=1)

indels_final = indels_final.reindex(columns=sorted(indels_final.columns))

print indels_final

# Generate ped file:
print "Generating PED file..."

# Create fields for PED file:
isolate_family_ids = list(indels_final.index)
sex = ["0"]*155
phenotype =  ["-9"]*155
pid = ["0"]*155
mid = ["0"]*155

indels_final.insert(0, "phenotype", value=phenotype)
indels_final.insert(0, "sex", value=sex)
indels_final.insert(0, "pid", value=pid)
indels_final.insert(0, "mid", value=mid)
indels_final.insert(0, "isolate_id", value=isolate_family_ids)
indels_final.insert(0, "family_id", value=isolate_family_ids)

print "Generating MAP file..."

colnames = list(indels_final.dtypes.index)
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

indels_final.to_csv(path_or_buf=output_ped, index=False, header=False, sep='\t')
mapfile2.to_csv(path_or_buf=output_map, index=False, header=False, sep='\t')

print "Done."