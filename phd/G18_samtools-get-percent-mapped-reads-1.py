import argparse
import os.path

parser = argparse.ArgumentParser()

parser.add_argument("--flagstat", help="samtools flagstat file")
parser.add_argument("--out", help="output file; can be appended to")

args = parser.parse_args()

with open(args.flagstat, 'r') as infile1:
    for line in infile1:
        if "mapped (" in line:
            line_ = line.strip().split(" ")[-3].strip("(").strip("%")

filename = args.flagstat.strip().split("/")[-1]

to_write = '\n' + filename + '\t' + line_

if os.path.isfile(args.out):
    with open(args.out, 'a') as outfile1:
        outfile1.write(to_write)
else:
    with open(args.out, 'w') as outfile:
        outfile.write("file\tpercent_reads_mapped" + to_write)
