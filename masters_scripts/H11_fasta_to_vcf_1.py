from sys import argv
import pandas
import numpy
import itertools

script, fasta, positions_file, outputdir = argv

# Create dictionary of fasta sequences, with header as key and sequence as value:
infile_list = []
infile_dict = {}
with open(fasta, 'r') as infile1:
    for line in infile1:
        if line.startswith(">"):
            infile_dict[line.strip().split('>')[1]] = list(next(infile1).strip())
            infile_list.append(line.strip().split('>')[1])

# Create list of high quality positions:
positions_list = []

with open(positions_file, 'r') as infile2:
    for line in infile2:
        positions_list.append(line.strip().split('-')[1])

outnames = pandas.Series(positions_list)

reference = infile_dict["REF"]

# Convert reference list into a series and label positions with high quality positions list:
ref_series = pandas.Series(reference, index=positions_list)

# Delete reference sequence from infile_dict:
del infile_dict["REF"]
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

chrom = ["ENA|FM209186|FM209186.1"]*len(ref_series)

vcfid = ["."]*len(ref_series)

qual = ["."]*len(ref_series)

vcffilter = ["."]*len(ref_series)

vcfinfo = ["."]*len(ref_series)

vcfformat = ["."]*len(ref_series)

output_file_dict = {}
for i in infile_dict:
	output_handle = outputdir + str(i) + ".vcf"
	output_file_dict[i] = output_handle

for isolatenumber in infile_list:
	if isolatenumber != "REF":	
		print isolatenumber
		table = []
		table.append(chrom)
		table.append(positions_list)
		table.append(vcfid)
		table.append(reference)
		table.append(infile_dict[isolatenumber])
		table.append(qual)
		table.append(vcffilter)
		table.append(vcfinfo)
		table.append(vcfformat)

		table = numpy.array(table)

		out_df = pandas.DataFrame(table)
		out_df = out_df.transpose()

		with open(output_file_dict[isolatenumber], 'a+') as outfile:
			outfile.write("##fileformat=VCFv4.1" + '\n' + "#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT" + '\n')
			out_df.to_csv(outfile, header=False, index=False, sep='\t')

	