from sys import argv
import pandas

script, pphe_file_list = argv

pphe_files = []
with open(pphe_file_list, 'r') as infile1:
	for line in infile1:
		pphe_files.append(line.strip())


pheno_set = set()
for i in pphe_files:
	print "Reading file: ", i.strip().split("/")[-1]
	dat = pandas.read_table(i, sep='\t', header=None, index_col=0)
	del dat[1]
	colnames = map(str, range(1,1001))
	dat.columns = colnames

	for col in colnames:
		collist = map(str, list(dat[col]))
		colstr = ''.join(collist)
		pheno_set.add(colstr)

print "Number of unique phenotype permutations: ", len(pheno_set)



# dat = pandas.read_table(pphe, sep='\t', header=None, index_col=0)
# del dat[1]
# colnames = map(str, range(1,1001))
# dat.columns = colnames

# def get_unique_seg(variants_df):
# 	inputcol = list(variants_df.columns)


# 	colstring_dict = {}
# 	colstring_set = set()
# 	for i in inputcol:
# 		# print type(i)
# 		col = list(variants_df[i])
# 		col = map(str, col)
# 		colstring = ''.join(col)
# 		if colstring not in colstring_set:
# 			colstring_set.add(colstring)
# 			colstring_dict[colstring] = [i]
# 		else:
# 			colstring_dict[colstring].append(i)

# 	print "Number of uniquely segregating variants: ", len(colstring_dict)
# 	to_drop = []
# 	for key in colstring_dict:
# 		key_value = colstring_dict[key]
# 		if len(key_value) > 1:
# 			key_value_minus_first = key_value[1:]
# 			for i in key_value_minus_first:
# 				to_drop.append(i)
# 		else:
# 			pass
# 	print "Number of identically-segregating variants to remove: ", len(to_drop)

# 	# output_df = variants_df.drop(labels=to_drop, axis=1, inplace=False)

# 	return [to_drop, colstring_dict]

# unique_pheno = get_unique_seg(dat)