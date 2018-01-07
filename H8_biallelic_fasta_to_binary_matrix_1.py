from sys import argv
import pandas

script, fasta, positions_file = argv

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
        positions_list.append(line.strip().split('-')[0])

outnames = pandas.Series(positions_list)

# Turn sequence dictionary into dataframe, with positions as the column names:
df1 = pandas.DataFrame.from_dict(infile_dict, orient='index')
df1.columns = positions_list

def convert_to_binary(column_vector):
    key_set = list(set(column_vector))
    allele_a = key_set[0]
    allele_b = key_set[1]
    column_vector1 = column_vector.replace(allele_a, 0)
    column_vector2 = column_vector1.replace(allele_b, 1)
    return column_vector2

out = df1.apply(convert_to_binary, axis=0)
print out
# out.to_csv(path_or_buf="/home/conrad/grad_school_scripts/h_notebook_files/test_binary_output.txt", sep='\t')