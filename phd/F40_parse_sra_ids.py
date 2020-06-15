import argparse
import os.path
import random
random.seed()

parser = argparse.ArgumentParser()

parser.add_argument("--input_file", help="input file")
parser.add_argument("--output_file", help="output_file")

args = parser.parse_args()

sra_list = []
with open(args.input_file, 'r') as infile1:
    for line in infile1:
        if "SRA:" in line:
            sra_list.append(line.strip().split(" ")[-1])

with open(args.output_file, 'w') as outfile1:
    outfile1.write("\n".join(sra_list))
