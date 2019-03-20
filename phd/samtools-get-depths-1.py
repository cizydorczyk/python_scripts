import argparse
import os.path
from Bio import SeqIO

parser = argparse.ArgumentParser()

parser.add_argument("--depths", help="samtools depths file")
parser.add_argument("--out", help="output file; can be appended to")
parser.add_argument("--reference", help="reference sequence (to get ref length)")

args = parser.parse_args()

# Get length of reference sequence:
record = list(SeqIO.parse(args.reference, "fasta"))
ref_length = len(record[0].seq)

# Parse depths to get counts >0 and >10:
g0_counts = 0
g10_counts = 0
with open(args.depths, 'r') as infile1:
    for line in infile1:
        depth = int(line.strip().split('\t')[2])
        if depth > 0:
            g0_counts += 1
        if depth > 10:
            g10_counts += 1

g0_percent = round(float(g0_counts)/float(ref_length)*100.00, 2)
g10_percent = round(float(g10_counts)/float(ref_length)*100.00, 2)

filename = args.depths.strip().split("/")[-1]

to_write = "\n" + filename + "\t" + str(ref_length) + "\t" + str(g0_percent) + "\t" + str(g10_percent)

if os.path.isfile(args.out):
    with open(args.out, 'a') as outfile:
        outfile.write(to_write)
else:
    with open(args.out, 'w') as outfile:
        outfile.write("ref_length\tpercent_positions_>0\tpercent_positions_>10" + to_write)
