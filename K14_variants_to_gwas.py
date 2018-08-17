from sys import argv
import pandas as pd
import itertools
import argparse

parser = argparse.ArgumentParser()

# Required arguments:
parser.add_argument("input_file", help="See --file_type. Fasta for SNPs, VCF for indels or SNPs, Roary Rtab for Roary pangenome output.")
parser.add_argument("output_prefix", help="Prefix to use for all output files.")
parser.add_argument("file_type", default="fasta", choices=["vcf", "roary", "fasta"], help="Input file type: 'fasta' (default), \
	'vcf', or 'roary' (roary Rtab).", nargs=1)
parser.add_argument("variant_type", default="snp", choices=['snp', 'indel', 'roary'], help="Input variant type: 'snp' (default), 'indel', \
	# or 'roary'.", nargs=1)
parser.add_argument("number_isolates", help="Number of isolates from which variants were collected.")

# File-specific options (may be required with certain input file types):
parser.add_argument("--positions_file", nargs=1, default=None, help="Required if --file_type is 'fasta'. '-' separated file, one entry per line, where the second field must be the SNP position, \
	corresponding (in order) to the SNP positions in the SNP fasta alignment. The first field can be anything. Eg. LESB58-254, where \
	LESB58 can be replaced by anything, and 254 is the SNP position.")

# Filtering options:
parser.add_argument("--unique_seg", help="Collapse identically-segregating variants.", action="store_true")
parser.add_argument("--maf", nargs=1, default=None, help="Minimum allele frequency (as decimal).", action="store")
parser.add_argument("--identical", help="Remove variants identical among all isolates.", action="store_true")
parser.add_argument("--multiallelic", help="Remove multiallelic variants (tri- and tetra-).", action="store_true")
parser.add_argument("--singletons", help="Remove singleton variants.", action="store_true")
parser.add_argument("--ncol", help="Remove variant columns with 'N'.", action="store_true")
args = parser.parse_args()

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

# Get columns with maf > specified maf:
def get_maf(column_vector):
	threshold = float(''.join(args.maf))
	key_set = set(column_vector)
	value_counts = list(column_vector.value_counts())
	sum_counts = float(sum(value_counts))
	value_freq = [float(i)/sum_counts for i in value_counts]

	if all(i > threshold for i in value_freq) == True:
		outcome = "keep"

	elif all(i > threshold for i in value_freq) == False:
		outcome = "drop"

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

def count_comments(filename):
    """Count comment lines (those that start with "#") in an optionally
    gzipped file.
    :param filename:  An optionally gzipped file.
    """
    comments = 0
    fn_open = gzip.open if filename.endswith('.gz') else open
    with fn_open(filename) as fh:
        for line in fh:
            if line.startswith('#'):
                comments += 1
            else:
                break
	return comments

def get_unique_seg(variants_df):
	inputcol = list(variants_df.columns)

	colstring_dict = {}
	colstring_set = set()
	for i in inputcol:
		# print type(i)
		col = list(variants_df[i])
		colstring = ''.join(col)
		if colstring not in colstring_set:
			colstring_set.add(colstring)
			colstring_dict[colstring] = [i]
		else:
			colstring_dict[colstring].append(i)

	print "Number of uniquely segregating variants: ", len(colstring_dict)
	to_drop = []
	for key in colstring_dict:
		key_value = colstring_dict[key]
		if len(key_value) > 1:
			key_value_minus_first = key_value[1:]
			for i in key_value_minus_first:
				to_drop.append(i)
		else:
			pass
	print "Number of identically-segregating variants to remove: ", len(to_drop)

	# output_df = variants_df.drop(labels=to_drop, axis=1, inplace=False)

	return [to_drop, colstring_dict]

#########################################################################################################################

if args.positions_file == None:
	if args.file_type == ["fasta"]:
		parser.error("Fasta input file (--file_type fasta) requires a positions file (--positions_file).")

if args.file_type == ['fasta'] and args.variant_type != ['snp']:
	parser.error("Variant type for file type 'fasta' must be 'snp'.")

if args.file_type == ['vcf'] and args.variant_type != ['indel']:
	parser.error("Variant type for file type 'vcf' must be 'indel'. SNPs are not currently handled in VCF format.")

if args.file_type == ['roary'] and args.variant_type != ['roary']:
	parser.error("Variant type for file type 'roary' must be 'roary'.")

# script, fasta, positions_file, output_matrix = argv

if args.file_type == ['fasta']:
	
	print 'Reading fasta and positions files...'

	infile_dict = {}
	
	with open(args.input_file, 'r') as infile1:
		for line in infile1:
			if line.startswith('>'):
				infile_dict[line.strip()[1:]] = list(next(infile1).strip())

	positions_list = []
	
	with open(''.join(args.positions_file), 'r') as infile2:
	    for line in infile2:
	        positions_list.append(line.strip().split('-')[1])

	outnames = pd.Series(positions_list)

	# Turn sequence dictionary into dataframe, with positions as the column names:
	variants1 = pd.DataFrame.from_dict(infile_dict, orient='index')
	variants1.columns = positions_list

if args.file_type == ['vcf']:

	print 'Reading vcf file...'

	comments = count_comments(''.join(args.input_file))

	variants_transposed = pd.read_csv(filepath_or_buffer=''.join(args.input_file), sep='\t', header=0, skiprows=(comments-1), index_col=1)

	del variants_transposed['#CHROM']
	del variants_transposed['ID']
	del variants_transposed['REF']
	del variants_transposed['ALT']
	del variants_transposed['QUAL']
	del variants_transposed['FILTER']
	del variants_transposed['INFO']
	del variants_transposed['FORMAT']

	variants1 = variants_transposed.transpose()

	variants1.replace(["."], ["A"], inplace=True)
	variants1 = variants1.apply(indels_to_binary, axis=0)

if args.file_type == ['roary']:

	print "Reading Roary Rtab file..."

	input_genes = pd.read_csv(filepath_or_buffer=''.join(args.input_file), sep='\t', header=0, index_col=0)
	variants1 = input_genes.transpose()

	## Code accessory genes with a unique number:
	# Pangenome:
	pangenome_guide_dict = {}
	pangenome_guide_list = []
	pangenome_colnames = list(variants1.dtypes.index)
	pangenome_numbers = range(1, len(pangenome_colnames)+1)

	for i, j in itertools.izip(pangenome_colnames, pangenome_numbers):
		pangenome_guide_dict[i] = j
		pangenome_guide_list.append(j)

	# Write gene key to file:
	print "Writing pangenome variants key to file..."

	output_genes_key = args.output_prefix + '_accessory_genes_key.txt'
	with open(output_genes_key, 'w') as outfile1:
		for j in pangenome_guide_dict:
			outfile1.write(str(j) + '\t' + str(pangenome_guide_dict[j]) + '\n')

	# Rename variants:
	variants1.rename(columns=pangenome_guide_dict, inplace=True)
	variants1 = variants1.applymap(str)

########################## Filtering part of the script ######################################

print "Number of variants before filtering: ", len(variants1.columns)

if args.identical == True:
	# Remove identical columns:
	possible_identical = variants1.apply(get_identical_columns, axis=0)
	variants1 = drop_columns(variants1, possible_identical)
	print "Number of variants after removing identical columns: ", len(variants1.columns)

if args.ncol == True:
	# Remove columns with 'N':
	possible_ncol = variants1.apply(get_n_columns, axis=0)
	variants1 = drop_columns(variants1, possible_ncol)
	print "Number of variants after removing columns with N: ", len(variants1.columns)

if args.multiallelic == True:
	# Remove multiallelic columns:
	possible_multiallelics = variants1.apply(get_multiallelic_columns, axis=0)
	variants1 = drop_columns(variants1, possible_multiallelics)
	print "Number of variants after removing multiallelic (tri- and tetra-) sites: ", len(variants1.columns)

if args.singletons == True:
	# Remove singletons:
	possible_singletons = variants1.apply(get_singleton_columns, axis=0)
	variants1 = drop_columns(variants1, possible_singletons)
	print "Number of variants after removing singleton sites: ", len(variants1.columns)
	
if args.maf != None:
	# Apply maf filter:
	possible_maf = variants1.apply(get_maf, axis=0)
	variants1 = drop_columns(variants1, possible_maf)
	print "Number of variants after applying maf filter: ", len(variants1.columns)

if args.unique_seg == True:

	if args.variant_type == ['snp']:

		# Variants1 is not binary prior to here because we want to have the real alleles at each SNP in the final PLINK files
		variants_binary = variants1.copy()
		variants_binary = variants_binary.apply(convert_to_binary, axis=0)

		unique_seg_output = get_unique_seg(variants_binary)

		to_drop = unique_seg_output[0]
		coldict = unique_seg_output[1]

		variants1 = variants1.drop(labels=to_drop, axis=1, inplace=False)

	if args.variant_type == ['indel'] or args.variant_type == ['roary']:

		# For indels and roary, variants1 is already binary at this stage
		unique_seg_output = get_unique_seg(variants1)

		to_drop = unique_seg_output[0]
		coldict = unique_seg_output[1]

		variants1 = variants1.drop(labels=to_drop, axis=1, inplace=False)

	# Write collapsed variants to file:
	print "Writing collapsed variants key to file..."

	output_collapsed_var_key = args.output_prefix + '_collapsed_variants_key.txt'
	with open(output_collapsed_var_key, 'w') as outfile1:
		outfile1.write("snp_kept" + '\t' + 'collapsed_snps' + '\n')
		for key in coldict:

			key_list = map(str, coldict[key])
			out_string = '\t'.join(key_list)

			
			outfile1.write(key_list[0] + '\t' + out_string + '\n')

print "Final number of variants: ", len(variants1.columns)

if args.variant_type == ['snp']:
	print "Writing binary matrix to file..."
	output_matrix = args.output_prefix + "_binary_matrix.txt"
	variants1_binary = variants1.copy()
	variants1_binary = variants1_binary.apply(convert_to_binary, axis=0)
	variants1_binary.to_csv(path_or_buf=output_matrix, sep='\t')

if args.variant_type == ['indel']:
	print "Writing binary matrix to file..."
	indels_binary = variants1.replace(["A", "G"], ["0", "1"])
	output_matrix = args.output_prefix + "_binary_matrix.txt"
	indels_binary.to_csv(path_or_buf=output_matrix, sep='\t')

if args.variant_type == ['roary']:
	print "Writing binary matrix to file..."
	output_matrix = args.output_prefix + "_binary_matrix.txt"
	variants1.to_csv(path_or_buf=output_matrix, sep='\t')

#################### Generate ped/map files #########################
if args.variant_type == ['roary']:
	variants1.replace(["0", "1"], ["A", "G"], inplace=True)

variants2 = variants1.copy()

# Rename columns to allow for easy sorting:
variants1_colnames = variants1.columns

variants1_new_colnames = []
variants2_new_colnames = []

for colname in variants1_colnames:
    variants1_new_colnames.append(float(colname) + 0.1)
    variants2_new_colnames.append(float(colname) + 0.2)

variants1.columns = variants1_new_colnames
variants2.columns = variants2_new_colnames

# Combine and sort the two variants dfs:
frames = [variants1, variants2]
variants_final = pd.concat(frames, axis=1)

variants_final = variants_final.reindex(columns=sorted(variants_final.columns))

# Generate ped file:
print "Generating PED file..."

# Create fields for PED file:
isolate_family_ids = list(variants_final.index)
sex = ["0"]*int(args.number_isolates)
phenotype =  ["-9"]*int(args.number_isolates)
pid = ["0"]*int(args.number_isolates)
mid = ["0"]*int(args.number_isolates)

variants_final.insert(0, "phenotype", value=phenotype)
variants_final.insert(0, "sex", value=sex)
variants_final.insert(0, "pid", value=pid)
variants_final.insert(0, "mid", value=mid)
variants_final.insert(0, "isolate_id", value=isolate_family_ids)
variants_final.insert(0, "family_id", value=isolate_family_ids)

print "Generating MAP file..."

colnames = list(variants_final.dtypes.index)
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

############################# Write ped/map files #################################

output_ped = args.output_prefix + ".ped"
output_map = args.output_prefix + ".map"

variants_final.to_csv(path_or_buf=output_ped, index=False, header=False, sep='\t')
mapfile2.to_csv(path_or_buf=output_map, index=False, header=False, sep='\t')

print "Done."