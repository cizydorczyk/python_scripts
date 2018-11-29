import argparse
from Bio.SeqIO.QualityIO import FastqGeneralIterator
import os.path

parser = argparse.ArgumentParser()

parser.add_argument("input_R1_fastq", help="input R1 fastq file")
parser.add_argument("input_R2_fastq", help="input R2 fastq file")
parser.add_argument("genome_length", help="genome length (in bp) against which to calculate coverage")
parser.add_argument("outputtsv", help="output tsv file with base counts")

args = parser.parse_args()
isolate = args.input_R1_fastq.split("/")[-1].split(".")[0]

print("Parsing isolate ", isolate)

R1_count = 0
R1_total_len = 0
with open(args.input_R1_fastq, 'r') as infile1:
    for title, seq, qual in FastqGeneralIterator(infile1):
        R1_count += 1
        R1_total_len += len(seq)

R2_count = 0
R2_total_len = 0
with open(args.input_R2_fastq, 'r') as infile2:
    for title, seq, qual in FastqGeneralIterator(infile2):
        R2_count += 1
        R2_total_len += len(seq)

total_bases = R1_total_len + R2_total_len
est_coverage = total_bases/int(args.genome_length)

print("\tNumber of R1 reads: " + str(R1_count) + "\n\tTotal number of R1 bases: " + str(R1_total_len) +
      "\n\tNumber of R2 reads: " + str(R2_count) + "\n\tTotal number of R2 bases: " + str(R2_total_len) +
      "\n\tTotal # bases: " + str(total_bases) + "\n\tEstimated coverage: " + str(est_coverage))

if os.path.isfile(args.outputtsv):
    with open(args.outputtsv, 'a') as outfile:
        outfile.write(isolate + "\t" + str(est_coverage) + "\n")
elif not os.path.isfile(args.outputtsv):
    with open(args.outputtsv, 'w') as outfile:
        outfile.write("Isolate\tCoverage\n" + isolate + "\t" + str(est_coverage) + "\n")
