import argparse
from Bio import SeqIO

parser = argparse.ArgumentParser()

parser.add_argument("--input_fasta", help="fasta from which to extract taret seq")
parser.add_argument("--output_fasta", help="fasta to write selected sequence to")
parser.add_argument("--target_seq", help="id of sequence to extract; no spaces, case sensitive")

args = parser.parse_args()

for record in SeqIO.parse(args.input_fasta, "fasta"):
    if record.id == args.target_seq:
        with open(args.output_fasta, "w") as outfile:
            SeqIO.write(record, outfile, "fasta")
