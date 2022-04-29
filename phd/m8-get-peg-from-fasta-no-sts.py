import argparse
from Bio import SeqIO

parser = argparse.ArgumentParser()

parser.add_argument("--input_pegs")
parser.add_argument("--input_fasta_list", help="comma separated list of dna feature fasta files from rast, one for each ST")
parser.add_argument("--output_fasta")

args = parser.parse_args()

with open(args.input_pegs, "r") as infile1:
    pegs = set([line.strip() for line in infile1])

fasta_input_list = args.input_fasta_list.split(",")
seq_to_write = []
for fasta in fasta_input_list:
    fasta_seqs = list(SeqIO.parse(fasta, "fasta"))
    for peg in pegs:
        for seq_record in fasta_seqs:
            if peg == seq_record.id.split("|")[1]:
                seq_to_write.append(seq_record)
            else:
                continue

print("Number of pegs to find: ", len(pegs))
print("Number of pegs found: ", len(seq_to_write))

# Write pegs to file:
with open(args.output_fasta, "w") as outfile1:
    SeqIO.write(seq_to_write, outfile1, "fasta")
