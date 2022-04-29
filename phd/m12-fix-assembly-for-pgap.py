import argparse
from Bio import SeqIO

parser = argparse.ArgumentParser()

parser.add_argument("--input_assembly", help="assembly to fix for pgap")
parser.add_argument("--output_assembly", help="assembly fixed and ready for pgap")

args = parser.parse_args()

with open(args.input_assembly, "r") as infile1:
    fasta_records = list(SeqIO.parse(infile1, "fasta"))

print("Number of original fasta records: ", len(fasta_records))

records_to_write = []
for record in fasta_records:
    if len(record.seq) >= 200:
        new_record_description = record.description.replace("=", "", 2)
        record.description = new_record_description
        records_to_write.append(record)
    elif len(record.seq) < 200:
        continue

print("New number of fasta records: ", len(records_to_write))

with open(args.output_assembly, "w") as outfile1:
    SeqIO.write(records_to_write, outfile1, "fasta")