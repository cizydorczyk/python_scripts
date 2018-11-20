from sys import argv
import pandas
import collections

script, inputblastfile, outputfile = argv

blast_raw = pandas.read_table(inputblastfile, comment='#', header=None)

# contig_stats = pandas.read_table(inputcontigstats, sep='\t')

contig_set = set(list(blast_raw[0]))

output_colors = {}
for i in sorted(contig_set):
	df1 = blast_raw.loc[blast_raw[0] == i]
	if "Escherichia coli" in set(list(df1[5])):
		output_colors[i+"_mapping"] = [i+"_mapping", 'deepskyblue', 'blue']
		print i, 'deepskyblue'
	else:
		output_colors[i+"_mapping"] = [i+"_mapping", 'red', 'red4']
		print i, 'red'

colors_df = pandas.DataFrame.from_dict(output_colors, orient='index')
colors_df.columns = ['Name', 'color1', 'color2']

colors_df.to_csv(path_or_buf=outputfile, sep='\t', header=True)

# output_df = pandas.merge(contig_stats, colors_df, on='Name')

# output_df.to_csv(path_or_buf=outputfile, sep='\t', header=True)
