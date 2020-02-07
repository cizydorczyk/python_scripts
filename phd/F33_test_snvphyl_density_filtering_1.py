import argparse
import os

parser = argparse.ArgumentParser()

parser.add_argument("--input_no_density_snp_tsv", help="snvphyl snvTable.tsv file")
parser.add_argument("--density", help="# snps in window to filter")
parser.add_argument("--window", help="Window size in which to filter 'density' number of SNPs")
parser.add_argument("--ref", help="reference genome")

args = parser.parse_args()

vcf_dict = {}
vcf_pos_list = []
with open(args.input_no_density_snp_tsv, 'r') as infile1:
    for line in infile1:
        if not line.startswith("#"):
            if "filtered-coverage" not in line:
                if "filtered-mpileup" not in line:
                    if "filtered-invalid" not in line:
                        line_elements = line.strip().split("\t")
                        vcf_dict[line_elements[1]] = line_elements
                        vcf_pos_list.append(int(line_elements[1]))

from Bio import SeqIO

ref = list(SeqIO.parse(args.ref, "fasta"))[0].seq

# Function to get circuclar string...
def access_char(string, i):
    return string[i % len(string)]


to_drop = [] # final positions to remove
for i in range(1,len(ref)+1):
    end = i + int(args.window)-1

    pos_in_window = [] # temporary list of positions falling in window; if the # of these exceeds the density, add these positions to the final to drop list

    for position in vcf_pos_list:
        if i <= position <= end:
            pos_in_window.append(position)
        else:
            continue

    if len(pos_in_window) >= int(args.density):
        to_drop = to_drop + pos_in_window

print(len(vcf_dict))
print(len(to_drop))

for key in set(to_drop):
    vcf_dict.pop(str(key))

print(len(vcf_dict))

# for pos in vcf_dict:
#     vcf_dict.pop(pos)
#
# print(len(vcf_dict))
