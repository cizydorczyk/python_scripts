from sys import argv
import pandas

script, fasta, positions_file, output_ped2, output_map2 = argv

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

###############################From here is a section specific to my fasta file; change as required############

reference = infile_dict["REF"]

# Convert reference list into a series and label positions with high quality positions list:
ref_series = pandas.Series(reference, index=positions_list)

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

############################Case-specific section ends##########################################################

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

to_drop_N_colnames = []
for i, j in out1.iteritems():
    if j == 'drop':
        to_drop_N_colnames.append(i)

df1.drop(to_drop_N_colnames, axis=1, inplace=True)
print "New df length after removing columns with N: ", len(df1.columns)

# Get and drop identical columns:
def get_identical_columns(column_vector):
    key_set = set(column_vector)
    if len(key_set) == 1:
        outcome = 'drop'
    else:
        outcome = 'keep'
    return outcome

out2 = df1.apply(get_identical_columns, axis=0)

to_drop_identical_cols = []
for i, j in out2.iteritems():
    if j == 'drop':
        to_drop_identical_cols.append(i)

df1.drop(to_drop_identical_cols, axis=1, inplace=True)
print "New df length after removing identical columns: ", len(df1.columns)

# Get and drop multiallelic columns:
def get_multiallelic_columns(column_vector):
    key_set = set(column_vector)
    if len(key_set) > 2:
        outcome = 'drop'
    else:
        outcome = 'keep'
    return outcome

out3 = df1.apply(get_multiallelic_columns, axis=0)

to_drop_multiallelics = []
for i,j in out3.iteritems():
    if j == 'drop':
        to_drop_multiallelics.append(i)

df1.drop(to_drop_multiallelics, axis=1, inplace=True)
print "New df length after removing multiallelic columns: ", len(df1.columns)

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

to_drop_singletons = []
for i, j in out4.iteritems():
    if j == 'drop':
        to_drop_singletons.append(i)

df1.drop(to_drop_singletons, axis=1, inplace=True)
print "New df length after removing singletons: ", len(df1.columns)

# Convert df1 to binary format with abscence = "A" for one allele (ref base) and presence = "G" for the other (alt base):
def convert_to_binary(column_vector):
    key_set = list(set(column_vector))
    allele_a = key_set[0]
    allele_b = key_set[1]
    column_vector1 = column_vector.replace(allele_a, "A")
    column_vector2 = column_vector1.replace(allele_b, "G")
    return column_vector2

# Create 2 columns per variant (first create 2 dfs & rename columns appropriately):
out1 = df1.apply(convert_to_binary, axis=0)
out2 = df1.apply(convert_to_binary, axis=0)

colnames = list(out1.dtypes.index)

newcolnames1 = []
for i in colnames:
    newname = float(i) + 0.1
    newcolnames1.append(newname)

newcolnames2 = []
for i in colnames:
    newname = float(i) + 0.2
    newcolnames2.append(newname)

out1.columns = newcolnames1
out2.columns = newcolnames2

frames = [out1, out2]
snps1 = pandas.concat(frames, axis=1)

snps2 = snps1.reindex_axis(sorted(snps1.columns), axis=1)

# Generate ped/map files with all fields (dummy):
print "Generating complete (dummy) PED file..."

rownames = list(snps2.index)

sex = ["1"]*155
phenotype =  ["-9"]*155

snps2.insert(0, "phenotype", value=phenotype)
snps2.insert(0, "sex", value=sex)
pid = ["0"]*155
mid = ["0"]*155
snps2.insert(0, "pid", value=pid)
snps2.insert(0, "mid", value=mid)
snps2.insert(0, "isolate_id", value=rownames)
snps2.insert(0, "family_id", value=rownames)

print snps2.shape

print "Generating complete (dummy) MAP file..."

colnames3 = list(snps2.dtypes.index)
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

mapfile2 = pandas.DataFrame()
mapfile2.insert(0, "bp_position", value=newcolnames)
pos = ["0"]*len(newcolnames)
mapfile2.insert(0, "pos", value=pos)
mapfile2.insert(0, "var_id", value=newcolnames)
mapfile2.insert(0, "chrom", value=["1"]*len(newcolnames))

snps2.to_csv(path_or_buf=output_ped2, index=False, header=False, sep='\t')
mapfile2.to_csv(path_or_buf=output_map2, index=False, header=False, sep='\t')

print "Done."