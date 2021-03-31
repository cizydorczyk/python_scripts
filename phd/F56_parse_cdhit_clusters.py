import argparse
from Bio import SeqIO
import os.path

parser = argparse.ArgumentParser()

parser.add_argument("--input_clusters", help="cdhit .clstr file")
parser.add_argument("--input_fasta", help="fasta used as input to cdhit")
parser.add_argument("--output_dir", help="directory to write indiv representative\
                    fasta files to")

args = parser.parse_args()

# read fasta file
seqs = list(SeqIO.parse(args.input_fasta, "fasta"))

# read clusters file
representative_seqs = []
with open(args.input_clusters, "r") as infile1:
    for line in infile1:
        if "*" in line:
            representative_seqs.append(line.strip().split(">")[1][:-5])

for seq in seqs:
    if seq.id in representative_seqs:
        output_handle = os.path.join(args.output_dir, seq.id)
        SeqIO.write(seq, output_handle, "fasta")
