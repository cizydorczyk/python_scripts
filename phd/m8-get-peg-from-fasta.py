import argparse
from Bio import SeqIO

parser = argparse.ArgumentParser()

parser.add_argument("--input_pegs")
parser.add_argument("--input_fasta")
parser.add_argument("--output_fasta")

args = parser.parse_args()

fasta_seqs = list(SeqIO.parse(args.input_fasta, "fasta"))

with open(args.input_pegs, "r") as infile1:
    pegs = set([line.strip() for line in infile1])

seq_to_write = []
for peg in pegs:
    for seq_record in fasta_seqs:
        if peg == seq_record.id.split("|")[1]:
            seq_to_write.append(seq_record)
            print(peg)
        else:
            continue

with open(args.output_fasta, "w") as outfile1:
    SeqIO.write(seq_to_write, outfile1, "fasta")
