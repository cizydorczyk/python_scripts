import argparse
from Bio.SeqIO.FastaIO import SimpleFastaParser
import os.path

parser = argparse.ArgumentParser()

parser.add_argument("--input_mlst", help="input file from mlst program")
parser.add_argument("--mlst_genes", help='comma separated list (lower case) of mlst genes in scheme; alphabetical order')
parser.add_argument("--mlst_alleles_directory", help='directory with fasta files for each gene with all alleles, named'
                                                   '"gene.fasta" where gene is the lowercase gene name same as in '
                                                   'mlst_genes')
parser.add_argument("--output_mlsa_alignment", help="main output file")
parser.add_argument("--nontyped_isolates", help="file with isolate records from SRST2 that could not be typed with "
                                              "certainty")


args = parser.parse_args()

# Create allele sequences dicts for each gene:
genes = args.mlst_genes.split(',')

dict1 = {}
dict2 = {}
dict3 = {}
dict4 = {}
dict5 = {}
dict6 = {}
dict7 = {}

with open(os.path.join(args.mlst_alleles_directory, genes[0] + ".fasta"), 'r') as infile1:
    for title, seq in SimpleFastaParser(infile1):
        title = title.lower()
        dict1[title] = seq
with open(os.path.join(args.mlst_alleles_directory, genes[1] + ".fasta"), 'r') as infile1:
    for title, seq in SimpleFastaParser(infile1):
        title = title.lower()
        dict2[title] = seq
with open(os.path.join(args.mlst_alleles_directory, genes[2] + ".fasta"), 'r') as infile1:
    for title, seq in SimpleFastaParser(infile1):
        title = title.lower()
        dict3[title] = seq
with open(os.path.join(args.mlst_alleles_directory, genes[3] + ".fasta"), 'r') as infile1:
    for title, seq in SimpleFastaParser(infile1):
        title = title.lower()
        dict4[title] = seq
with open(os.path.join(args.mlst_alleles_directory, genes[4] + ".fasta"), 'r') as infile1:
    for title, seq in SimpleFastaParser(infile1):
        title = title.lower()
        dict5[title] = seq
with open(os.path.join(args.mlst_alleles_directory, genes[5] + ".fasta"), 'r') as infile1:
    for title, seq in SimpleFastaParser(infile1):
        title = title.lower()
        dict6[title] = seq
with open(os.path.join(args.mlst_alleles_directory, genes[6] + ".fasta"), 'r') as infile1:
    for title, seq in SimpleFastaParser(infile1):
        title = title.lower()
        dict7[title] = seq

# Parse SRST2 output file:

output_dict = {}
nontyped_isolates = []

with open(args.input_mlst, 'r') as infile1:
    raw_mlst_data = list(infile1)
    header_line = raw_mlst_data[0]
    raw_mlst_data = raw_mlst_data[1:]
#print(dict1)
for line in raw_mlst_data:
    line_ = line.strip().split('\t')
    alleles = ''.join(line_[2:9])
    if '*' not in alleles and '?' not in alleles:
        isolate = ">" + line_[0] + "_ST" + line_[1]
        gene1_seq = dict1[genes[0] + line_[2]]
        gene2_seq = dict2[genes[1] + line_[3]]
        gene3_seq = dict3[genes[2] + line_[4]]
        gene4_seq = dict4[genes[3] + line_[5]]
        gene5_seq = dict5[genes[4] + line_[6]]
        gene6_seq = dict6[genes[5] + line_[7]]
        gene7_seq = dict7[genes[6] + line_[8]]
        complete_seq = gene1_seq + gene2_seq + gene3_seq + gene4_seq + gene5_seq + gene6_seq + gene7_seq
        output_dict[isolate] = complete_seq
    else:
        nontyped_isolates.append(line.strip())

with open(args.output_mlsa_alignment, 'w') as outfile1:
    to_write = []
    for key in output_dict:
        to_write.append(key + '\n' + output_dict[key])
    outfile1.write('\n'.join(to_write))

with open(args.nontyped_isolates, 'w') as outfile2:
    outfile2.write('\n'.join(nontyped_isolates))
