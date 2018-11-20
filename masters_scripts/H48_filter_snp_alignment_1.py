from sys import argv
import pandas

script, fasta, positions_file, output_multiall_fasta, output_multiall_positions, output_biall_fasta, output_biall_positions, output_biall_matrix = argv

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
        positions_list.append(line.strip().split('-')[1])

outnames = pandas.Series(positions_list)

# reference = infile_dict["REF"]

# # Convert reference list into a series and label positions with high quality positions list:
# ref_series = pandas.Series(reference, index=positions_list)

# # Delete reference sequence from infile_dict:
# del infile_dict["REF"]
# del infile_dict["528"]
# del infile_dict["529"]
# del infile_dict["530"]
# del infile_dict["531"]
# del infile_dict["532"]
# del infile_dict["548"]
# del infile_dict["549"]
# del infile_dict["550"]
# del infile_dict["538"]
# del infile_dict["542"]

# Turn sequence dictionary into dataframe, with positions as the column names:
df1 = pandas.DataFrame.from_dict(infile_dict, orient='index')
df1.columns = positions_list

print "Original number of SNPs: ", len(df1.columns)

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
    # If multiallelic, don't remove, as there may be a count of 1 but it's multiallelic, not a singleton:        
    else:
        outcome = "keep"
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

#########################################################################################################################

# Remove identical columns:
possible_identical_col = df1.apply(get_identical_columns, axis=0)
df1 = drop_columns(df1, possible_identical_col)
print "Number of SNPs after removing identical sites: ", len(df1.columns)

# Remove columns with N:
possible_n_col = df1.apply(get_n_columns, axis=0)
df1 = drop_columns(df1, possible_n_col)
print "Number of SNPs after removing sites with 'N': ", len(df1.columns)

# # Remove singletons:
# possible_singletons = df1.apply(get_singleton_columns, axis=0)
# df1 = drop_columns(df1, possible_singletons)
# print "Number of SNPs after removing singeltons: ", len(df1.columns)

# Get number of multiallelic (tri- & tetra-) sites:
# possible_multiallelic_col = df1.apply(get_multiallelic_columns, axis=0)
# df2 = df1.copy()
# df2 = drop_columns(df2, possible_multiallelic_col)
# df3 = df2.apply(convert_to_binary, axis=0)

# multiallelic_count = 0
# for i in possible_multiallelic_col:
# 	if i == 'drop':
# 		multiallelic_count += 1
# print "Number of multiallelic SNP sites: ", multiallelic_count
# print "Number of SNPs after removing multiallelic (tri- and tetra-) sites: ", len(df2.columns)

##########################################################################################################################

##### Write multiallelic SNPs to fasta/positions files:

output_multiallelic_positions = list(df1.columns)

# Write multiallelic positions to file:
with open(output_multiall_positions, 'w') as outfile2:
    for position in output_multiallelic_positions:
        outfile2.write(str(position) + '\n')

# Write multiallelic isolate sequences to file:
df1_dict = df1.transpose().to_dict(orient='list')

with open(output_multiall_fasta, 'w') as outfile:
    for key in df1_dict:
        outfile.write('>' + key + '\n' + ''.join(df1_dict[key]) + '\n')

##### Write biallelic SNPs to fasta/positions/matrix files:

# # Write biallelic positions to file:
# output_biallelic_positions = list(df2.columns)

# with open(output_biall_positions, 'w') as outfile3:
# 	for position in output_biallelic_positions:
# 		outfile3.write(str(position) + '\n')

# # Write biallelic isolate sequences to file:
# df2_dict = df2.transpose().to_dict(orient='list')

# with open(output_biall_fasta, 'w') as outfile4:
# 	for key in df2_dict:
# 		outfile4.write('>' + key + '\n' + ''.join(df2_dict[key]) + '\n')

# # Write biallelic SNPs to matrix:
# df3.to_csv(path_or_buf=output_biall_matrix, sep='\t')