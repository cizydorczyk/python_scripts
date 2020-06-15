import argparse

parser = argparse.ArgumentParser()

parser.add_argument("--fasta_file")
parser.add_argument("--submol_yaml")
parser.add_argument("--output_yaml")

args = parser.parse_args()

to_write = "report_usage: false\nfasta:\n\tclass: File\n\tlocation: " +  args.fasta_file + "\nsubmol:\n\tclass: File\n\tlocation: " + args.submol_yaml

with open(args.output_yaml, 'w') as outfile1:
    outfile1.write(to_write)
