import argparse
import os.path
from Bio import SeqIO

parser = argparse.ArgumentParser()

parser.add_argument("--input_log", help="input log file from vcftools tstv module")
parser.add_argument("--output_file", help="output file to write results to. will append if already exists.")
parser.add_argument("--isolate", help="isolate id")

args = parser.parse_args()

print(args.isolate)

with open(args.input_log, "r") as infile1:
    for line in infile1:
        if line.startswith("Ts/Tv ratio:"):
            tstv_ratio = line.strip().split(" ")[-1]

if not os.path.isfile(args.output_file):
    with open(args.output_file, "w") as outfile1:
        to_write = args.isolate + "\t" + tstv_ratio + "\n"
        outfile1.write(to_write)
elif os.path.isfile(args.output_file):
    with open(args.output_file, "a") as outfile1:
        to_write = args.isolate + "\t" + tstv_ratio + "\n"
        outfile1.write(to_write)
