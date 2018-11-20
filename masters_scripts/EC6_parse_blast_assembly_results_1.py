from sys import argv
import pandas
import collections

script, inputcontigstats, inputblastfile, outputfile = argv
print "Working on file: " + inputblastfile.split('/')[-1]

blast_raw = pandas.read_table(inputblastfile, comment='#', header=None)

contig_stats = pandas.read_table(inputcontigstats, sep='\t')

contig_set = set(list(blast_raw[0]))

output_colors = {}
for i in sorted(contig_set):
	df1 = blast_raw.loc[blast_raw[0] == i]
	if "Escherichia coli" in set(list(df1[5])):
		output_colors[i+"_mapping"] = [i+"_mapping", 'deepskyblue', 'blue']
		# print i, 'deepskyblue'
	else:
		output_colors[i+"_mapping"] = [i+"_mapping", 'red', 'red4']
		# print i, 'yellow'

colors_df = pandas.DataFrame.from_dict(output_colors, orient='index')
colors_df.columns = ['Name', 'color1', 'color2']

output_df = pandas.merge(contig_stats, colors_df, on='Name')
# print output_df

output_df.to_csv(path_or_buf=outputfile, sep='\t', header=True)



# contig_list = []
# no_result_contig_list = []

# with open(inputblastfile, 'r') as infile:
#     contents = list(infile)

# print "Working on file: " + str(inputblastfile.split('/')[-1])
# data = pandas.read_table(inputcontigstats, sep="\t")
# length_list = list(data["Consensus length"])
# coverage_list = list(data["Average coverage"])

# query_list = []
# fill_color = []
# out_color = []

# for line in contents:
#     if "# Query" in line:
#         # if (lnum+4) < len(contents):
#         contig_num = line.strip().split(' ')[-1]
#         query_list.append(contig_num)

# query_dict = {}
# query_set = set(query_list)
# contents_set = set(contents)

# print contents_set
# for i in query_set:
# 	for j in contents_set:
# 		if i not in query_dict:
# 			if i in j:
# 				query_dict[i] = [j]
# 		else:
# 			if i in j:
# 				query_dict[i].append(j)

