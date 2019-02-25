import argparse
from Bio.SeqIO.FastaIO import SimpleFastaParser

parser = argparse.ArgumentParser()

parser.add_argument("input_fasta", help="input fasta file with headers to shorten")
parser.add_argument("output_fasta", help="output fasta with shortened headers")

args = parser.parse_args()

corrected_fasta_dict = {}
with open(args.input_fasta, 'r') as infile1:
    for title, seq in SimpleFastaParser(infile1):
        title = title.strip().split('\t')[0]
        corrected_fasta_dict[title] = seq

with open(args.output_fasta, 'w') as outfile1:
    to_write = []
    for key in sorted(corrected_fasta_dict):
        to_write.append(">" + key + "\n" + corrected_fasta_dict[key] + "\n")
    outfile1.write(''.join(to_write))

