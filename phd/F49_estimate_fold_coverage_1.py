import argparse
import subprocess
from Bio.SeqIO.QualityIO import FastqGeneralIterator

################### NOTE ###################
# The script currently does NOT subsample reads...but it can! Just need to
# uncomment lines 55-60 (inclusive) and 64-72 (inclusive). In its current
# state, it just prints estimated genome sizes, depths, and subsampling factors
# to the terminal.

parser = argparse.ArgumentParser()

parser.add_argument("--R1", help="input R1 fastq file")
parser.add_argument("--R2", help="input R2 fastq file")
parser.add_argument("--isolate", help='isolate number/id')
parser.add_argument("--output_file", help="output file to write cov to")

args = parser.parse_args()

#print("\nWorking on isolate " + args.isolate + "...")

R1_file = args.R1.strip().split("/")[-1]
R2_file = args.R2.strip().split("/")[-1]

def GetBasesFromFq(inputfq):
    base_count = 0
    with open(inputfq, 'r') as infile1:
        for title, seq, qual in FastqGeneralIterator(infile1):
            base_count += len(seq)
    return base_count

#print("\nGetting total bases from fastq files...")
R1_bases = GetBasesFromFq(args.R1)
R2_bases = GetBasesFromFq(args.R2)
total_bases = R1_bases + R2_bases

#print("\nEstimating genome size with 'mash'...")

mash_cmd = ["mash", "sketch", "-o", "temp.R1.mash.sketch", "-k", "32", "-m", "3", "-s", "10000", "-r", args.R1]
mash_cmd_raw_output = subprocess.run(mash_cmd, stderr=subprocess.PIPE)
subprocess.run(["rm", "temp.R1.mash.sketch.msh"])

mash_cmd_output = mash_cmd_raw_output.stderr.decode('utf-8')
mash_gsize_est = float(mash_cmd_output.strip().split('\n')[0].split(' ')[-1])
print(args.isolate)
print("Estimated genome size is ~" + str(mash_gsize_est))
# mash_cov_est = mash_cmd_output.strip().split('\n')[1].split(' ')[-1]

#print('\nEstimating sequencing depth...')
orig_depth = total_bases/mash_gsize_est

print("Original sequencing depth is ~" + str(round(orig_depth,3)))

with open(args.output_file, 'a') as outfile:
    outfile.write(args.isolate + "\t" + str(total_bases) + "\t" + str(mash_gsize_est) + "\t" + str(round(orig_depth, 3)) + "\n")
