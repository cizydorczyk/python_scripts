from Bio.SeqIO.QualityIO import FastqGeneralIterator
from Bio import SeqIO
import argparse
import os.path
import gzip

parser = argparse.ArgumentParser()

parser.add_argument("--r1", help="r1 file for isolate")
parser.add_argument("--r2", help="r2 file for isolate")
parser.add_argument("--est_length", help='estimated length to calc coverage against', type=int)
parser.add_argument("--out", help='output file for all isolates (will be appended to)')

args = parser.parse_args()

# get sample name:
sample = args.r1.strip().split('/')[-1].split('.')[0].split('_')[0]

# R1 reads and bases counts:
bases_count_r1 = 0
reads_count_r1 = 0
with gzip.open(args.r1, 'rt') as infile1:
    for title, seq, qual in FastqGeneralIterator(infile1):
        bases_count_r1 += len(seq)
        reads_count_r1 += 1

# R2 reads and bases counts:
bases_count_r2 = 0
reads_count_r2 = 0
with gzip.open(args.r2, 'rt') as infile2:
    for title, seq, qual in FastqGeneralIterator(infile2):
        bases_count_r2 += len(seq)
        reads_count_r2 += 1

# total reads and bases counts:
total_bases_count = bases_count_r1 + bases_count_r2
total_reads_count = reads_count_r1 + reads_count_r2

# reference length:


# coverage:
coverage = round(total_bases_count/args.est_length, 2)

# write results to file:
if os.path.isfile(args.out):
    with open(args.out, 'a') as outfile:
        outfile.write('\n' + sample + '\t' + str(reads_count_r1) + '\t' + str(reads_count_r2) + '\t' + str(total_reads_count) + '\t' + str(bases_count_r1) + '\t' + str(bases_count_r2) + '\t' + str(total_bases_count) + '\t' + str(coverage))

else:
    header = "sample\tr1_reads\tr2_reads\ttotal_reads\tr1_bases\tr2_bases\ttotal_bases\tcoverage"
    with open(args.out, 'w') as outfile:
        outfile.write(header + '\n' + sample + '\t' + str(reads_count_r1) + '\t' + str(reads_count_r2) + '\t' + str(total_reads_count) + '\t' + str(bases_count_r1) + '\t' + str(bases_count_r2) + '\t' + str(total_bases_count) + '\t' + str(coverage))
