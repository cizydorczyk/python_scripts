import argparse
import os.path
from Bio import SeqIO

### NOTE
# The script ignores complex variants currently. It also reports the GC content for the reference genomes and NOT the isolate genomes.


parser = argparse.ArgumentParser()

parser.add_argument("--input_fasta", help="input vcf file")
parser.add_argument("--ref", help="reference sequence fasta file")
parser.add_argument("--output_folder", help="output folder")

args = parser.parse_args()

# Parse input fasta
sequences = list(SeqIO.parse(args.input_fasta, "fasta"))

# Parse reference fasta
reference_sequence = list(SeqIO.parse(args.ref, "fasta"))[0]

for sequence in sequences:
    output_file = os.path.join(args.output_folder, sequence.id + ".fasta")
    with open(output_file, "a") as outfile:
        SeqIO.write(reference_sequence, outfile, "fasta")
        SeqIO.write(sequence, outfile, "fasta")
