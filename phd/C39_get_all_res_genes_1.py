import argparse
from Bio import SeqIO

parser = argparse.ArgumentParser()

parser.add_argument("--input_file", help="input file listing all card output files, one per line")
parser.add_argument("--output_file", help="output file with set of genes listed")

args = parser.parse_args()

# Get list of output files:
file_list =[]
with open(args.input_file, "r") as infile1:
    for line in infile1:
        file_list.append(line.strip())

# Get de-duplicated list of resistance genes:
res_genes_list = []

for i in file_list:
    with open(i, 'r') as infile2:
        for line in infile2:
            if not line.startswith("ORF"):
                line_elements = line.strip().split("\t")
                res_genes_list.append(line_elements[8])

res_genes_dedup = sorted(list(set(res_genes_list)))

# Write results to file:
to_write = '\n'.join(res_genes_dedup)
with open(args.output_file, 'w') as outfile1:
    outfile1.write(to_write)
