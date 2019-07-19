import argparse
import itertools
import subprocess
import os
import numpy as np
from Bio import SeqIO
from collections import Counter

# import os.path

parser = argparse.ArgumentParser()

parser.add_argument("--wd", help="working directory (where snippy will output files)")
parser.add_argument("--out", help="output file w/ distances")
parser.add_argument("--cleanup", default=False, action="store_true", help="cleanup all files from pipeline except final output? default False (without this flag))")
parser.add_argument("--distance_xfactor", default=1, help="factor to multiply raw normalized distance by (to simply make bigger); default = 1")

args = parser.parse_args()

# set wd:
os.chdir(args.wd)

# Get file list:
snpdists_files = sorted([f for f in os.listdir(args.wd) if f.endswith('.matrix')])
for file in snpdists_files:
    print(file)

# Get pairwise distances:

output_distances = []

for file in snpdists_files:

    ## Uncomment below to get partially corrected distances for length of comparisons:
    # Read recombination-masked fasta alignment for obtaining corrected distances:
    # fasta_fname = ".".join(file.split(".")[0:-1])
    # seq_list = list(SeqIO.parse(fasta_fname, "fasta"))
    #
    # seq1 = str(seq_list[0].seq)
    # counter1 = Counter(seq1)
    # seq1_count = counter1["A"] + counter1["C"] + counter1["G"] + counter1["T"] # to get total bases, add the following: + counter1["N"] + counter1["X"] + counter1["-"]
    #
    # seq2 = str(seq_list[1].seq)
    # counter2 = Counter(seq2)
    # seq2_count = counter2["A"] + counter2["C"] + counter2["G"] + counter2["T"] # to get total bases, add the following: + counter1["N"] + counter1["X"] + counter1["-"]
    #
    # # Calculate normalization distance:
    # normalization_value = min(seq1_count, seq2_count)


    # Get SNP distances:
    with open(file, 'r') as infile1:
        for line in infile1:
            if line.startswith("snp"):
                snp_distance = next(infile1).strip().split('\t')[2]
                ## Uncomment line below to get comparison length corrected distances:
                # normalized_distance = int(snp_distance) / normalization_value * int(args.distance_xfactor)
                isolate_1 = line.strip().split('\t')[1]
                isolate_2 = line.strip().split('\t')[2]

                dist_str = isolate_2 + '\t' + isolate_1 + '\t' + snp_distance ## add the following to get distances partially corrected for length of comparisons: + '\t' + str(normalization_value) + '\t' + str(normalized_distance)
                output_distances.append(dist_str)

# Write output distances to file:
with open(args.out, 'w') as outfile1:
    outfile1.write('\n'.join(output_distances))
