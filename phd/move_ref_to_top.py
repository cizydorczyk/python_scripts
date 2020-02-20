#!/home/conrad/anaconda3/bin/python

import argparse
from Bio import SeqIO

parser = argparse.ArgumentParser()

parser.add_argument("--fasta", help="fasta sequence from which to get target sequence")
parser.add_argument("--out", help="fasta sequence header to get from fasta file")

args = parser.parse_args()

print("Moving reference to top of fasta file...")

fasta_seqs = list(SeqIO.parse(args.fasta, "fasta"))

ref = []
nonref_records = []
for record in fasta_seqs:
    if record.id.lower() == "reference":
        ref.append(record)
    else:
        nonref_records.append(record)

output_records = ref + nonref_records

with open(args.out, "w") as outfile:
	SeqIO.write(output_records, outfile, "fasta")

print("Done.")
