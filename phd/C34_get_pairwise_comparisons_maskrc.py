import argparse
import os
from Bio import SeqIO
from Bio.Seq import Seq
from Bio.Seq import MutableSeq
from Bio.SeqRecord import SeqRecord
from Bio.Seq import MutableSeq
import csv

parser = argparse.ArgumentParser()

parser.add_argument("--wd", help="working directory (where output files will live if full paths are not provided)")
parser.add_argument("--rc_regions", help="bed file with regions to mask recombination before running SNP dists")
parser.add_argument("--fasta", help="input fasta file")
parser.add_argument("--maskchar", default="X", help="masking character, default = 'X'")
parser.add_argument("--out", help="output fasta alignment with recombinant regions masked")

args = parser.parse_args()

# set wd:
os.chdir(args.wd)

# Import CFML files and identify recombinant regions
seq_list = []
with open(args.rc_regions) as infile1:
	RCseqs = csv.reader(infile1, delimiter='\t')
	for row in RCseqs:
		seq = row[0]
		RCstart = row[1]
		RCstop = row[2]
		seq_list.append([seq, RCstart, RCstop])

# Read fasta file:
fasta_seqs = list(SeqIO.parse(args.fasta, "fasta"))

# Mask recombinant regions:
masked_seq = []

for i in fasta_seqs:
    seq = MutableSeq(str(i.seq))
    for j in seq_list:
        if i.id == j[0]: # j[0] is the sequence id in recombinant regions list
            start_mask = int(j[1]) - 1 # 1 based positions are 1 less in 0-based indexing
            end_mask = int(j[2]) # last index in range is not included in python
            len_mask = end_mask - start_mask
            seq[start_mask:end_mask] = args.maskchar * len_mask

    masked_seq.append(SeqRecord(Seq(str(seq)), i.id, description=""))

for i in masked_seq:
    print("Number of characters masked in sequence " + i.id + ": " + str(str(i.seq).count(args.maskchar)))

# Write masked sequences to file:
SeqIO.write(masked_seq, args.out, "fasta")
