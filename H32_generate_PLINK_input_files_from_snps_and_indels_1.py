from sys import argv
import pandas as pd
import itertools

script, input_vcf, fasta, positions_file, output_matrix, output_ped, output_map = argv

# Dicts for renaming isolates so that all variants have the same isolate name formats:
snp_isolates_dict = {"298C":"298", "404C":"404","405C":"405", "406C":"406", "410C":"410", "411C":"411", "417C":"417", "418C":"418", "419C":"419", "420N2":"420", "459C":"459", "541-2":"541", "553C":"553", "559-2":"559", "572-1":"572", "573-2":"573", "573-3":"574"}

############################### Load data #################################################

## Load indels:
input_indels = pd.read_csv(filepath_or_buffer=input_vcf, sep='\t', header=0, skiprows=35, index_col=1)
# input_indels_2 = pd.read_csv(filepath_or_buffer=input_vcf, sep='\t', header=0, skiprows=35, index_col=1)

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

## Load snps:
# Create dictionary of fasta sequences, with header as key and sequence as value:
infile_dict = {}
with open(fasta, 'r') as infile1:
    for line in infile1:
        if line.startswith(">"):
            infile_dict[line.strip()[1:]] = list(next(infile1).strip())

# Create list of high quality positions:
positions_list = []

with open(positions_file, 'r') as infile2:
    for line in infile2:
        positions_list.append(int(line.strip().split('-')[1]))

outnames = pd.Series(positions_list)

# Delete reference sequence from infile_dict:
del infile_dict["REF"]
del infile_dict["528"]
del infile_dict["529"]
del infile_dict["530"]
del infile_dict["531"]
del infile_dict["532"]
del infile_dict["548"]
del infile_dict["549"]
del infile_dict["550"]
del infile_dict["538"]
del infile_dict["542"]

# Turn sequence dictionary into dataframe, with positions as the column names:
snps1 = pd.DataFrame.from_dict(infile_dict, orient='index')
snps1.columns = positions_list
snps1.rename(index=snp_isolates_dict, inplace=True)

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

##################################################################################################################
########## Get total number of variants (indels + snps) and check to make sure there are none in common ##########

print "Original number of indels: ", len(indels1.columns)
print "Original number of snps: ", len(snps1.columns)
print "Total number of variants: ", len(indels1.columns) + len(snps1.columns)

snpcolnames_set = set(list(snps1.columns))
indelcolnames_set = set(list(indels1.columns))

unique_variants = snpcolnames_set|indelcolnames_set
print "Number of variants with unique positions: ", len(unique_variants)

common_variants = snpcolnames_set&indelcolnames_set
print "Number of variants with identical positions (should be zero, 0): ", len(common_variants)

if len(common_variants) != 0:
    print "NOTE: Number of variants with identical positions should be zero (0)..."
    print "NOTE: Indels sharing a position with a snp will be recoded with a position # + 0.5."

# Renumber indels that overlap with snps:
new_variants = {}
for i in common_variants: 
    new_variants[i] = i + 0.5

indels1.rename(columns=new_variants, inplace=True)

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

## Snps:

print "Number of SNPs before removing identical, multiallelic, singleton, and N-containing columns: ", len(snps1.columns)

# Remove columns with 'N' (there may be some):
possible_ncol = snps1.apply(get_n_columns, axis=0)
snps1 = drop_columns(snps1, possible_ncol)
print "Number of snps after removing columns with N: ", len(snps1.columns)

# Remove identical columns (some possible):
possible_identical = snps1.apply(get_identical_columns, axis=0)
snps1 = drop_columns(snps1, possible_identical)
print "Number of snps after removing identical columns: ", len(snps1.columns)

# Remove multiallelic columns (some possible):
possible_multiallelics = snps1.apply(get_multiallelic_columns, axis=0)
snps1 = drop_columns(snps1, possible_multiallelics)
print "Number of snps after removing multiallelic columns: ", len(snps1.columns)

# Remove singletons (some possible):
possible_singletons = snps1.apply(get_singleton_columns, axis=0)
snps1 = drop_columns(snps1, possible_singletons)
print "Number of snps after removing singleton columns: ", len(snps1.columns)

# Make binary:
snps1 = snps1.apply(convert_to_binary, axis=0)

print "Final number of snps: ", len(snps1.columns)

## Create dataframe from combined snps, indels, and genes:
frames = [indels1, snps1]
variants1 = pd.concat(frames, axis=1, join='outer')
variants2 = variants1.copy()

############ Format and write variants dataframe to file (for treeWAS/other software requiring basic matrix) ##########

variants_matrix = variants1.replace(["A", "G"], [0, 1])

variants_matrix.to_csv(path_or_buf=output_matrix, sep='\t')

#################### Generate ped/map files #########################

colnames1 = list(variants1.dtypes.index)

newcolnames1 = []
newcolnames2 = []

for i in colnames1:
	newcolnames1.append(float(i) + 0.1)
	newcolnames2.append(float(i) + 0.2)

variants1.columns = newcolnames1
variants2.columns = newcolnames2

frames2 = [variants1, variants2]
variants3 = pd.concat(frames2, axis=1)

variants3 = variants3.reindex(sorted(variants3.columns), axis=1)

## Generate ped/map files with all fields (dummy):
print "Generating complete PED file..."

rownames = list(variants3.index)
sex = ["1"]*155
phenotype =  ["-9"]*155
pid = ["0"]*155
mid = ["0"]*155

variants3.insert(0, "phenotype", value=phenotype)
variants3.insert(0, "sex", value=sex)
variants3.insert(0, "pid", value=pid)
variants3.insert(0, "mid", value=mid)
variants3.insert(0, "isolate_id", value=rownames)
variants3.insert(0, "family_id", value=rownames)

print "Generating complete MAP file..."

colnames3 = list(variants3.dtypes.index)
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
    elif ".7" in str(i):
        pass
    elif ".6" in str(i):
        newcolnames3.append(str(i)[0:-2] + str(0.5))
    else:
        newcolnames3.append(str(i)[0:-2])

mapfile2 = pd.DataFrame()
mapfile2.insert(0, "bp_position", value=newcolnames3)
pos = ["0"]*len(newcolnames3)
mapfile2.insert(0, "pos", value=pos)
mapfile2.insert(0,  "var_id", value=newcolnames3)
mapfile2.insert(0, "chrom", value=["1"]*len(newcolnames3))

######################## Write ped/map and variants key to file #####################################

variants3.to_csv(path_or_buf=output_ped, index=False, header=False, sep='\t')
mapfile2.to_csv(path_or_buf=output_map, index=False, header=False, sep='\t')

print "Done."
