import argparse
from Bio import SeqIO

parser = argparse.ArgumentParser()

parser.add_argument("--input_gff", help="input gff file from RAST")
parser.add_argument("--fasta", help="fasta file of contigs used in RAST annotation")
parser.add_argument("--output_gff", help="Formatted GFF file ready to be used with PIRATE")

args = parser.parse_args()

contigs = list(SeqIO.parse(args.fasta, "fasta"))

sequence_regions = []
fasta_seqs = []
for i in contigs:
    sequence_region = "##sequence-region " + str(i.id) + "\t1\t" + str(len(i.seq))# Create sequence-region strings
    sequence_regions.append(sequence_region)

    fasta_seq = ">" + str(i.id) + "\n" + str(i.seq).upper() # Create fasta sequences to add to end of file; uppercase all letters
    fasta_seqs.append(fasta_seq)

sequence_region_to_write = "\n".join(sequence_regions)
fasta_seqs_to_write = "\n".join(fasta_seqs)
fasta_line = "##FASTA"

# Read original gff file:
orig_gff = "\n".join([line.rstrip("\n") for line in open(args.input_gff)][1:])
gff_header = "##gff-version 3"

# Organize output file:
to_write = gff_header + "\n" + sequence_region_to_write + "\n" + orig_gff + "\n" + fasta_line + "\n" + fasta_seqs_to_write

# Write output file:
with open(args.output_gff, 'w') as outfile:
    outfile.write(to_write)
