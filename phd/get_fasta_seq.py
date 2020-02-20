#!/home/conrad/anaconda3/bin/python

import argparse
from Bio import SeqIO

parser = argparse.ArgumentParser()

parser.add_argument("--fasta", help="fasta sequence from which to get target sequence")
parser.add_argument("--seqid", help="fasta sequence header to get from fasta file")

args = parser.parse_args()

fasta_seqs = list(SeqIO.parse(args.fasta, "fasta"))

for record in fasta_seqs:
    if record.id == args.seqid:
        print(">" + record.id + "\n" + str(len(record.seq)))
