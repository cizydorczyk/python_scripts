import argparse

parser = argparse.ArgumentParser()

# General script options:
parser.add_argument("--fasta_file")

args = parser.parse_args()

files = []
with open(args.fasta_file, "r") as infile1:
    for line in infile1:
        files.append(line.strip().split(".")[0])

with open("empty_profiles_file.txt", "w") as outfile1:
    to_write = "ST" + "\t" + "\t".join(files) + "\n" + "1\t" + "\t".join(["1"]*17596)
    outfile1.write(to_write)
