#!/home/conrad/anaconda3/bin/python

import argparse
from Bio import SeqIO

parser = argparse.ArgumentParser()

parser.add_argument("--fasta", help="alignment to remove reference from")
parser.add_argument("--out", help="output fasta file with no reference sequence")

args = parser.parse_args()

records = list(SeqIO.parse(args.fasta, "fasta"))
records2 = [i for i in records if i.id.lower() != "reference"]
for i in records2:
	print(i.id)

with open(args.out, "w") as outfile:
	SeqIO.write(records2, outfile, "fasta")
