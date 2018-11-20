from sys import argv
import pandas

script, fasta, positions_file, output_fasta, output_positions = argv

# Create dictionary of fasta sequences, with header as key and sequence as value:
infile_dict = {}
with open(fasta, 'r') as infile1:
    for line in infile1:
        if line.startswith(">"):
            infile_dict[line.strip()] = list(next(infile1).strip())

# Create list of high quality positions:
positions_list = []

with open(positions_file, 'r') as infile2:
    for line in infile2:
        positions_list.append(line.strip().split('-')[1])

outnames = pandas.Series(positions_list)

reference = infile_dict[">REF"]

# Convert reference list into a series and label positions with high quality positions list:
ref_series = pandas.Series(reference, index=positions_list)

# Delete reference sequence from infile_dict:
del infile_dict[">REF"]
del infile_dict[">528"]
del infile_dict[">529"]
del infile_dict[">530"]
del infile_dict[">531"]
del infile_dict[">532"]
del infile_dict[">548"]
del infile_dict[">549"]
del infile_dict[">550"]
del infile_dict[">538"]
del infile_dict[">542"]

# Turn sequence dictionary into dataframe, with positions as the column names:
df1 = pandas.DataFrame.from_dict(infile_dict, orient='index')
df1.columns = positions_list

print "Original number of df columns: ", len(df1.columns)

# Get and drop columns with N:
def get_n_columns(column_vector):
    key_set = set(column_vector)
    if 'N' in key_set:
        outcome = "drop"
    else:
        outcome = "keep"
    return outcome

out1 = df1.apply(get_n_columns, axis=0)
out1.rename(index=outnames, inplace=True)

to_drop_N_colnames = []
for i, j in out1.iteritems():
    if j == 'drop':
        to_drop_N_colnames.append(i)

df1.drop(to_drop_N_colnames, axis=1, inplace=True)
print "New df length after removing columns with N: ", len(df1.columns)

positions_list2 = [i for i in positions_list if i not in to_drop_N_colnames]
outnames2 = pandas.Series(positions_list2)


# Get and drop identical columns:
def get_identical_columns(column_vector):
    key_set = set(column_vector)
    if len(key_set) == 1:
        outcome = 'drop'
    else:
        outcome = 'keep'
    return outcome

out2 = df1.apply(get_identical_columns, axis=0)
out2.rename(index=outnames2, inplace=True)

to_drop_identical_cols = []
for i, j in out2.iteritems():
    if j == 'drop':
        to_drop_identical_cols.append(i)

df1.drop(to_drop_identical_cols, axis=1, inplace=True)
print "New df length after removing identical columns: ", len(df1.columns)

positions_list3 = [i for i in positions_list2 if i not in to_drop_identical_cols]
outnames3 = pandas.Series(positions_list3)


# Get and drop multiallelic columns:
def get_multiallelic_columns(column_vector):
    key_set = set(column_vector)
    if len(key_set) > 2:
        outcome = 'drop'
    else:
        outcome = 'keep'
    return outcome

out3 = df1.apply(get_multiallelic_columns, axis=0)
out3.rename(index=outnames3, inplace=True)

to_drop_multiallelics = []
for i,j in out3.iteritems():
    if j == 'drop':
        to_drop_multiallelics.append(i)

df1.drop(to_drop_multiallelics, axis=1, inplace=True)
print "New df length after removing multiallelic columns: ", len(df1.columns)

positions_list4 = [i for i in positions_list3 if i not in to_drop_multiallelics]
outnames4 = pandas.Series(positions_list4)

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

out4 = df1.apply(get_singleton_columns, axis=0)
out4.rename(index=outnames4, inplace=True)

to_drop_singletons = []
for i, j in out4.iteritems():
    if j == 'drop':
        to_drop_singletons.append(i)

df1.drop(to_drop_singletons, axis=1, inplace=True)
print "New df length after removing singletons: ", len(df1.columns)

df1.to_csv(path_or_buf="/home/conrad/Data/h_notebook_files/H8/biallelic_snps_matrix.txt", sep='\t')

positions_list5 = [i for i in positions_list4 if i not in to_drop_singletons]

# Convert dataframe and reference sequence (a series) to a dictionary and list, respectively (for easier writing to file):
df1_dict = df1.transpose().to_dict(orient='list')

# Write isolate sequences to file:
with open(output_fasta, 'w') as outfile:
    for key in df1_dict:
        outfile.write(key + '\n' + ''.join(df1_dict[key]) + '\n')

# Write positions to new list:
with open(output_positions, 'w') as outfile2:
    for position in positions_list5:
        outfile2.write(str(position) + '\n')