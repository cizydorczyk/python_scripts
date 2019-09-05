import argparse
from Bio import SeqIO

parser = argparse.ArgumentParser()

parser.add_argument("--mincontiglen", help="minimum contig length to keep")
parser.add_argument("--input_contigs", help="input contigs file to filter")
parser.add_argument("--out", help="output filtered contigs file")

args = parser.parse_args()

contigs = list(SeqIO.parse(args.input_contigs, "fasta"))
contigs_to_keep = []

# Identify contigs shorter than specified length and remove:
for contig in contigs:
    length = len(contig.seq)
    if length >= int(args.mincontiglen):
        contigs_to_keep.append(contig)

# Write good contigs to file:
SeqIO.write(contigs_to_keep, args.out, "fasta")
