import argparse
from Bio import SeqIO

parser = argparse.ArgumentParser()

parser.add_argument("--input_fasta", help="fasta from which to extract taret seq")
parser.add_argument("--output_fasta", help="fasta to write selected sequence to")
parser.add_argument("--remove_seq", help="file with seq ids to remove, one per line")
parser.add_argument("--keep_seq", help="seq to extract; cannot be used with remove_seq")

args = parser.parse_args()

if args.remove_seq != "" and args.keep_seq != "":
    print("Cannot use both keep and remove seq options")
    exit()

if args.remove_seq != "" and args.keep_seq == "":
    seq_to_remove = []
    with open(args.remove_seq, 'r') as infile1:
        for line in infile1:
            seq_to_remove.append(line.strip())

    for record in SeqIO.parse(args.input_fasta, "fasta"):
        if record.id not in seq_to_remove:
            with open(args.output_fasta, "a") as outfile:
                SeqIO.write(record, outfile, "fasta")

elif args.remove_seq == "" and args.keep_seq != "":
    seq_to_keep = []
    with open(args.keep_seq, 'r') as infile1:
        for line in infile1:
            seq_to_keep.append(line.strip())

    for record in SeqIO.parse(args.input_fasta, "fasta"):
        if record.id in seq_to_keep:
            with open(args.output_fasta, "a") as outfile:
                SeqIO.write(record, outfile, "fasta")
